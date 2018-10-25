'''astro_reduce -- A simple CCD images reducer for the Paris Observatory.'''

from sys import exit
from glob import glob
from json import loads
from collections import defaultdic

from astropy.io import fits
from scipy.signal import fftconvolve
import numpy as np
import click

# Paths and extensions.
TMP = "tmp"
OBJ = "objects"
DARK = "darks"
FLAT = "flats"
TMP_PNG = "tmp_png"
AUX = "aux"

# Read data from list of files and return aligned and meaned version.
def align_and_mean(in_files, nX, nY):
    '''Return align and meaned image from list of file names to read data.'''
    if len(in_files) == 1:
        return fits.getdata(infiles[0])
    
    # Collect arrays and crosscorrelate all with the first (except the first).
    images = [fits.getdata(fname) for fname in infiles]
    correlations = [fftconvolve(datas[0], image[::-1,::-1], mode='same')\
            for image in images[1:]]
    shift_indices = [np.unravel_index(np.argmax(corr_array, axis=None))\
            for corr_array in correlations]
    deltas = [(ind[0] + int(nX / 2), ind[1] + int(nY / 2))\
            for ind in shift_indices]

    # Roll the images and return their mean.
    data_to_mean = [np.roll(image, deltas[i], axis=(0,1))\
            for (i, image) in enumerate(images)]
    data_to_mean.append(images[0])
    return np.array(data_to_mean).mean(axis=0)

# Return the filter and exposure (strings) from an object or flat file name.
# "NGC1000_1_V_1000_0.fits" gives ("1", "V", "1000")
def fname_bits(fname):
    '''Return the series, filter and exposure from a file name.'''
    pieces = fname.split('_')
    return (pieces[1], pieces[-3], pieces[-2])

@click.command()
@click.argument('conf_file', type=click.File('r'))
@click.option('--png', '-p', is_flag=True,
        help='Write png format of intermediary images in png folder.')
@click.option('--interpolate', '-i', is_flag=True,
        help='Interpolate existing dark field images if some are missing.')
@click.option('--verbose', '-v', is_flag=True,
        help='Enables verbose mode.')
def cli(conf_file, png, interpolate, verbose):
    '''Reduce CCD images from objects with flat and dark field images.'''

    # Parse configuration file to obtain configuration dictionary.
    conf_dic = loads("".join(conf_file.read().split()))
    dn, fn = conf_dic['dark_name'], conf_dic['flat_name']

    # Obtain list of all object, dark, flat field files.
    object_files = dict([(obj,
        glob(f"{OBJ}/{obj}*")) for obj in conf_dic['objects']])
    dark_files = dict([(exp, 
        glob(f"{DARK}/{dn}_{exp}*")) for exp in conf_dic['exposures']])
    flat_files = dict([(filt,
        glob(f"{FLAT}/{fn}_{filt}*")) for filt in conf_dic['filters']])

    # Check if files exist.
    for key in object_files:
        if not object_files[key]:
            click.echo(f"Did not find files for {key} object. Exiting.")
            click.echo(f"Did you put them in the {OBJ} directory?")
            exit(1)
    for key in dark_files:
        if not dark_files[key]:
            click.echo(f"Did not find files for {key}ms exposure darks. Exit.")
            click.echo(f"They should be in the {DARK} directory.")
            exit(1)
    for key in flat_files:
        if not flat_files[key]:
            click.echo(f"Did not find files for the {key} filter flats.")
            click.echo(f"They should be in the {FLAT} directory.")
            exit(1)

    if verbose:
        click.echo("Files found:")
        click.echo("Objects:")
        for obj in conf_dic['objects']:
            click.echo(f"{obj}: {object_files[obj]}")
        click.echo("Dark fields:")
        for exp in conf_dic['exposures']:
            click.echo(f"{exp}ms: {dark_files[exp]}")
        click.echo("Flat fields:")
        for filt in conf_dic['filters']:
            click.echo(f"{filt}: {flat_files[filt]}")
   
   # STEP 1: Write the master dark files (medians of darks) for each exposure.
   for exp in dark_files:
       mdark_data = np.median(np.array([fits.getdata(fitsfile) \
               for fitsfile in dark_files[exp]]), axis=0)
       mdark_header = fits.getheader(dark_files[exp][0])
       fits.writeto(f"{TMP}/mdark_{exp}.fits", mdark_data, mdark_header,
               overwrite=True)

    # STEP 2: Write master flat files (medians of flats) for each filter.
    for filt in flat_files:
        mflat_data = np.median(np.array([fits.getdata(fitsfile) \
                for fitsfile in flat_files[filt]]), axis=0)
        mflat_header = fits.getheader(flat_files[filt][0])
        fits.writeto(f"{TMP}/mflat_{filt}.fits", mflat_data, mflat_header,
                overwrite=True)

    # STEP 3: Write master transmissions (master flat - master dark of same exp)
    # normalized.
    for filt in flat_files:
        fits_name = f"{TMP}/mflat_{filt}.fits" # Name of mflat file.
        _, _, exp = fname_bits(fits_name) # Strings.
        
        # Data of corresponding master dark (same exposure) and mflat (filter).
        mflat_data = fits.getdata(fits_name)
        mflat_header = fits.getheader(fits_name)
        mdark_data = fits.getdata(f"{TMP}/mdark_{exp}.fits")
        
        # (Flat - Dark) normalized
        mtrans_data = (mflat_data - mdark_data)\
                    / (mflat_data - mdark_data).mean(axis = 0)
        fits.writeto(f"{TMP}/mtrans_{filt}.fits", mtrans_data, mflat_header,
                overwrite=True)

    # STEP 4: Reduce all the object images with corresponding filter mflat
    # and exposure mdark. Do this regardless of series and number in series.
    for obj in object_files:
        for fname in object_files[obj]:
            filt, exp = fname_bits(fname)

            # Corresponding darks and flats
            mdark_data = fits.getdata(f"{TMP}/mdark_{exp}.fits")
            mtrans_data = fits.getdata(f"{TMP}/mtrans_{filt}.fits")

            # (Raw - Dark) / Trans
            aux_data = (fits.getdata(fname) - mdark_data) / mtrans_data
            aux_header = fits.getheader(fname)
            fits.writeto(f"{TMP}/{fname.split('.fit')[0]}_{AUX}.fits")

    # STEP 5: Within series, realign the aux images and write the mean images
    # of those. You are left with one image per object per filter per exposure
    # per series.
    for obj in object_files:
        # Group all the object files by *tag*, i.e. by series, filter, exposure
        name_tag_hash = [(fname, f"{fname_bits(fname)}")\
                for fname in object_files[obj]]
        names_per_tag = defaultdic(list)
        for name, tag in name_tag_hash:
            names_per_tag[tag].append(name)

        # Now you align images which have the same tag.
        for tag in names_per_tag:
            # Rebuild series, filter and exposure from tag (they are those of 
            # eg the first name in the list.
            s, f, e = fname_bits(names_per_tag[tag][0])

            # Calculate aligned and meaned image from all images with same tag.
            final_data \
                = align_and_mean(glob(f"{TMP}/{obj}_{s}_{f}_{e}_*_{AUX}.fits"),
                        conf_dic['nX'], conf_dic['nY'])
            final_header = fits.getheader(names_per_tag[tag][0])
            fits.writeto(f"{RED}/{obj}_{s}_{f}_{e}.fits", final_data,
                    final_header, overwrite=True)


            
