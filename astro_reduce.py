'''astro_reduce -- A simple CCD images reducer for the Paris Observatory.'''

from sys import exit
from os.path import basename, exists
from os import mkdir
from glob import glob
from json import loads
from collections import defaultdict
import re

from astropy.io import fits
from scipy.signal import fftconvolve
import numpy as np
import click

# Paths and extensions.
TMP = 'tmp'
OBJ = 'objects'
DARK = 'darks'
FLAT = 'flats'
AUX = 'aux'
RED = 'reduced'

# Read data from list of files and return aligned and meaned version.
def align_and_median(infiles):
    '''Return the madian of the aligned images from input list of file name.'''
    if len(infiles) == 1:
        return fits.getdata(infiles[0])
    
    # Collect arrays and crosscorrelate all with the first (except the first).
    images = [fits.getdata(fname) for fname in infiles]
    nX, nY = images[0].shape
    correlations = [fftconvolve(images[0], image[::-1,::-1], mode='same')\
            for image in images[1:]]
    shift_indices = [np.unravel_index(np.argmax(corr_array, axis=None),
            corr_array.shape) for corr_array in correlations]
    deltas = [(ind[0] + int(nX / 2), ind[1] + int(nY / 2))\
            for ind in shift_indices]

    # Roll the images and return their mean.
    data_to_mean = [np.roll(image, deltas[i], axis=(0,1))\
            for (i, image) in enumerate(images[1:])]
    data_to_mean.append(images[0])
    return np.median(data_to_mean, axis=0)

# Return the filter and exposure (strings) from an object file name.
# 'NGC1000_1_V_1000_0.fits' gives ('1', 'V', '1000')
def fname_bits(fname):
    '''Return the series, filter and exposure from a file name.'''
    pieces = fname.split('_')
    return (pieces[1], pieces[-3], pieces[-2])

# Write png from fits version of image, in same directory.
def write_png(filename, plt):
    '''Write PNG version from fits file.'''
    with fits.open(filename) as fits_file:
        plt.figure(len(filename))
        plt.imshow(fits_file[0].data, aspect='auto', origin='lower',
            cmap='magma')
        plt.colorbar()
        plt.savefig(f'{filename.split(".fit")[0]}.png', bbox_inches='tight')
        plt.close(len(filename))

@click.command()
@click.argument('conf_file', type=click.File('r'))
@click.option('--verbose', '-v', is_flag=True,
        help='Enables verbose mode (recommended).')
@click.option('--tmppng', '-t', is_flag=True,
        help='Write PNG format of intermediary images after reduction.')
@click.option('--redpng', '-r', is_flag=True,
        help='Write PNG format of reduced images after reduction.')
@click.option('--interpolate', '-i', is_flag=True,
        help='Interpolate existing dark field images if some are missing.')
def cli(conf_file, verbose, tmppng, redpng, interpolate):
    '''Reduce CCD images from objects with flat and dark field images.'''

    # Parse configuration file to obtain configuration dictionary.
    if verbose:
        click.echo(f'Parsing configuration file {conf_file.name}.')
    conf_dic = loads(''.join(conf_file.read().split()))
    dn, fn = conf_dic['dark_name'], conf_dic['flat_name']

    # Obtain list of all object, dark, flat field files.
    object_files = dict([(obj,
        glob(f'{OBJ}/{obj}*')) for obj in conf_dic['objects']])
    dark_files = dict([(exp, 
        glob(f'{DARK}/{dn}_{exp}*')) for exp in conf_dic['exposures']])
    flat_files = dict([(filt,
        glob(f'{FLAT}/{fn}_{filt}*')) for filt in conf_dic['filters']])

    # Check if files exist.
    for key in object_files:
        if not object_files[key]:
            click.echo(f'Did not find files for {key} object.')
            click.echo(f'Did you put them in the `{OBJ}` directory?')
            click.echo('Exiting.')
            exit(1)
    for key in dark_files:
        if not dark_files[key] and not interpolate:
            # If the interpolate option is off and there are some darks 
            # missing, exit.
            click.echo(f'Did not find files for {key}ms exposure darks.')
            click.echo(f'They should be in the `{DARK}` directory.')
            click.echo('If you want to interplate the missing dark fields '
                       'from the existing ones, use the `--interpolate`'
                       'option.')
            click.echo('Exiting.')
            exit(1)
    for key in flat_files:
        if not flat_files[key]:
            click.echo(f'Did not find files for the {key} filter flats.')
            click.echo(f'They should be in the `{FLAT}` directory.')
            exit(1)

    # Report found files
    reg = re.compile('_\d*\.')
    if verbose:
        click.echo('Files found:')
        click.echo(f'- Objects (`{OBJ}`):')
        for obj in conf_dic['objects']:
            uniq_names = set([reg.sub('_*.', basename(name)) \
                    for name in object_files[obj]]) 
            click.echo(f'-- {obj}: {uniq_names}')

        click.echo(f'- Dark fields (`{DARK}`):')
        for exp in conf_dic['exposures']:
            uniq_names = set([reg.sub('_*.', basename(name)) \
                    for name in dark_files[exp]]) 
            click.echo(f'-- {exp}ms: {uniq_names if uniq_names else None}')

        click.echo(f'- Flat fields (`{FLAT}`):')
        for filt in conf_dic['filters']:
            uniq_names = set([reg.sub('_*.', basename(name)) \
                    for name in flat_files[filt]]) 
            click.echo(f'-- {filt}: {uniq_names}')

    # STEP 0: Create directory for tmp and reduced images if not existent
    if verbose:
        click.echo('Creating folders for reduced and intermadiary images.')
    if not exists(RED):
        mkdir(RED)
    if not exists(TMP):
        mkdir(TMP)
   
    # STEP 1: Write the master dark files (medians of darks) 
    # for each available exposure.
    all_exposures = conf_dic['exposures']
    available_exposures = [exp for exp in dark_files if dark_files[exp]]
    for exp in available_exposures:
        mdark_data = np.median(np.array([fits.getdata(fitsfile) \
                for fitsfile in dark_files[exp]]), axis=0)
        mdark_header = fits.getheader(dark_files[exp][0])
        if verbose:
            click.echo(f'Writing {exp}ms master dark image.')
        fits.writeto(f'{TMP}/mdark_{exp}.fits', mdark_data, mdark_header,
                overwrite=True)

    # STEP 1.5: If there are some missing darks and the interpolate option
    # is on, then interpolate the master darks.
    # We use least squares linear interpolation, ie we calculate `a` and `b`
    # such that (missing_dark) = a * (exposure_time) + b
    # Exit if there are no darks at all
    if not available_exposures:
        click.echo("There are no dark files at all! Cannot interpolate...")
        click.echo("Exiting.")
        exit(1)

    if len(available_exposures) == 1:
        # If there's only one, a = (only_dark) / (its exposure), b = 0
        only_exp = available_exposures[0]
        only_mdark = fits.getdata(f'{TMP}/mdark_{only_exp}.fits')
        a = only_mdark / float(only_exp)
        b = np.zeros_like(only_mdark)
    else:
        # If not, if you want to fit y = a * x + b,
        # then the LS solution is:
        # a = (<xy> - <x><y>) / (<x ** 2> - <x> ** 2)
        # b = <y> - a * <x>
        mxy = np.mean([float(exp) * fits.getdata(f'{TMP}/mdark_{exp}.fits')\
                for exp in available_exposures], axis=0)
        mx = np.mean([float(exp) for exp in available_exposures])
        my = np.mean([fits.getdata(f'{TMP}/mdark_{exp}.fits')\
                for exp in available_exposures], axis=0)
        mx2 = np.mean([float(exp) ** 2 for exp in available_exposures])

        # a and b
        a = (mxy - mx * my) / (mx2 - mx ** 2)
        b = my - mx * a
        
    # Write all the missing master darks!
    for exp in list(set(all_exposures) - set(available_exposures)):
        new_mdark_data = float(exp) * a + b
        if verbose:
            click.echo(f'Interpolating {exp}ms master dark image.')
        fits.writeto(f'{TMP}/mdark_{exp}.fits', new_mdark_data,
                    overwrite=True)

    # STEP 2: Write master transmission files for each filter:
    # (median of flats - dark of same exposure) normalized.
    for filt in flat_files:
        # Example file of series to get exposure and header.
        fits_name = flat_files[filt][0]
        exp = fits_name.split('.fit')[0].split('_')[-2] # String exposure.
        mflat_data = np.median([fits.getdata(fitsfile) \
                for fitsfile in flat_files[filt]], axis=0)
        mflat_header = fits.getheader(fits_name)

        # Data of corresponding master dark (same exposure) and mflat (filter).
        mdark_data = fits.getdata(f'{TMP}/mdark_{exp}.fits')
        
        # (Flat - Dark) normalized
        mtrans_data = (mflat_data - mdark_data)\
                    / (mflat_data - mdark_data).mean(axis = 0)
        if verbose:
            click.echo(f'Writing {filt} filter master transmission image.')
        fits.writeto(f'{TMP}/mtrans_{filt}.fits', mtrans_data, mflat_header,
                overwrite=True)

    # STEP 3: Reduce all the object images with corresponding filter mflat
    # and exposure mdark. Do this regardless of series and number in series.
    for obj in object_files:
        if verbose:
            click.echo(f'Writing auxiliary image for {obj} object.')
        for fname in object_files[obj]:
            bfname = basename(fname)
            _, filt, exp = fname_bits(bfname)

            # Corresponding darks and flats
            mdark_data = fits.getdata(f'{TMP}/mdark_{exp}.fits')
            mtrans_data = fits.getdata(f'{TMP}/mtrans_{filt}.fits')

            # (Raw - Dark) / Trans
            aux_data = (fits.getdata(fname) - mdark_data) / mtrans_data
            aux_header = fits.getheader(fname)
            fits.writeto(f'{TMP}/{bfname.split(".fit")[0]}_{AUX}.fits',
                    aux_data, aux_header, overwrite=True)

    # STEP 4: Within series, realign the aux images and write the median images
    # of those. You are left with one image per object per filter per exposure
    # per series.
    for obj in object_files:
        # Group all the object files by *tag*, i.e. by series, filter, exposure
        name_tag_hash = [(basename(fname), f'{fname_bits(basename(fname))}')\
                for fname in object_files[obj]]
        names_per_tag = defaultdict(list)
        for name, tag in name_tag_hash:
            names_per_tag[tag].append(name)

        # Now you align images which have the same tag.
        for tag in names_per_tag:
            # Rebuild series, filter and exposure from tag (they are those of 
            # eg the first name in the list.)
            s, f, e = fname_bits(names_per_tag[tag][0])

            # Calculate aligned and medianed image
            # from all images with same tag.
            to_median = glob(f'{TMP}/{obj}_{s}_{f}_{e}_*_{AUX}.fits')
            final_data = align_and_median(to_median)
            final_header = fits.getheader(f'{OBJ}/{names_per_tag[tag][0]}')
            if verbose:
                click.echo(f'Writing realigned images for series {s}/{f}/{e} '
                           f'of object {obj}.')
            fits.writeto(f'{RED}/{obj}_{s}_{f}_{e}.fits', final_data,
                    final_header, overwrite=True)

    # STEP 5: If specified in the options red-ong and tmppng, write
    # PNG versions of all the tmp and reduced images.
    if redpng or tmppng:
        import matplotlib.pyplot as plt

    if redpng:
        if verbose:
            click.echo('Writing PNG versions of reduced images.')
        for ffile in glob(f'{RED}/*.fits'):
            write_png(ffile, plt)

    if tmppng:
        if verbose:
            click.echo('Writing PNG versions of intermediate images.')
        for ffile in glob(f'{TMP}/*.fits'):
            write_png(ffile, plt)
    
    if verbose:
        click.echo('All done.')
    # ALL DONE.
