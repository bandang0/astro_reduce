'''astro_reduce -- A Simple CCD Images Reducer for the Paris Observatory.'''

from sys import exit
from os.path import basename, exists
from os import mkdir, getcwd
from shutil import copy
from glob import glob
from json import loads, dump
from collections import defaultdict
from re import sub

import click
import numpy as np
from astropy.io import fits
from scipy.signal import fftconvolve


# Paths and extensions.
di = 'dark'
fi = 'flat'
BDARK = 'DARK'
BFLAT = 'FLAT'
BOBJ = 'ORIGINAL'
OBJ = 'ar_objects'
DARK = 'ar_darks'
FLAT = 'ar_flats'
TMP = 'tmp'
RED = 'reduced'
AUX = 'aux'

# Read data from list of files and return aligned and meaned version.
def align_and_median(infiles):
    '''Return the median of the re-aligned images from list of file names.'''
    if len(infiles) == 1:
        return fits.getdata(infiles[0])

    # Collect arrays and crosscorrelate all with the first (except the first).
    images = [fits.getdata(fname) for fname in infiles]
    nX, nY = images[0].shape
    correlations = [fftconvolve(images[0], image[::-1, ::-1], mode='same')
                    for image in images[1:]]

    # For each image determine the coordinate of maximum cross-correlation.
    shift_indices = [
        np.unravel_index(
            np.argmax(
                corr_array,
                axis=None),
            corr_array.shape) for corr_array in correlations]
    deltas = [(ind[0] - int(nX / 2), ind[1] - int(nY / 2))
              for ind in shift_indices]

    # Roll the images and return their median.
    realigned_images = [np.roll(image, deltas[i], axis=(0, 1))
                      for (i, image) in enumerate(images[1:])]
    realigned_images.append(images[0])
    return np.median(realigned_images, axis=0)

# Return exposure and custom file name for a dark field image
def d_read(fname):
    '''Return exposure and custom file name for a dark field image.'''
    head = fits.getheader(fname)
    if 'EXPTIME' in head.keys():
        exp = int(1000 * head['EXPTIME'])
    elif 'EXPOSURE' in head.keys():
        exp = int(1000 * head['EXPOSURE'])
    else:
        raise IOError(f'No exposure keyword in header of `{fname}`.')
    return exp, f'{d}_{exp}_{abs(hash(fname))}.fits'

# Return filter, exposure and custom file name for a flat field image
def f_read(fname):
    '''Return filter, exposure and custom file name for a flat field image.'''
    head = fits.getheader(fname)
    # Filter.
    if 'FILTER' in head.keys():
        fil = sub('[- _]', '', head['FILTER'])
    else:
        raise IOError(f'No filter keyword in header of `{fname}`.')

    # Exposure.
    if 'EXPTIME' in head.keys():
        exp = int(1000 * head['EXPTIME'])
    elif 'EXPOSURE' in head.keys():
        exp = int(1000 * head['EXPOSURE'])
    else:
        raise IOError(f'No exposure keyword in header of `{fname}`.')

    return fil, exp, f'{f}_{fil}_{exp}_{abs(hash(fname))}.fits'

# Return object, filter, exposure and custom file name for an object image.
def o_read(fname):
    '''Return object, filter, exposure and custom file name for object image.'''
    head = fits.getheader(fname)
    # Object.
    if 'OBJECT' in head.keys():
        obj = sub('[ _]', '-', head['OBJECT'])
    else:
        raise IOError(f'No object keyword in header of `{fname}`.')

    # Filter.
    if 'FILTER' in head.keys():
        fil = sub('[- _]', '', head['FILTER'])
    else:
        raise IOError(f'No filter keyword in header of `{fname}`.')

    # Exposure.
    if 'EXPTIME' in head.keys():
        exp = int(1000 * head['EXPTIME'])
    elif 'EXPOSURE' in head.keys():
        exp = int(1000 * head['EXPOSURE'])
    else:
        raise IOError(f'No exposure keyword in header of `{fname}`.')

    return obj, fil, exp, f'{obj}_{fil}_{exp}_{abs(hash(fname))}.fits'

# Write the configuration file for images in current directory.
def write_conf_file(objects, exposures, filters, conf_file_name):
    '''Write the configuration file from list of objects, exposures, filters.'''
    print(objects)
    print(type(objects))
    conf_dic = {'objects': objects,
                'exposures': exposures,
                'filters': filters}
    with open(conf_file_name, 'w') as cdfile:
        dump(conf_dic, cdfile, indent=2)

# Return the filter and exposure (strings) from an object file name.
# 'NGC1000_V_1000_0.fits' gives ('V', '1000')
def fname_bits(fname):
    '''Return the filter and exposure from an object file name.'''
    pieces = fname.split('.fit')[0].split('_')
    return (pieces[-3], pieces[-2])

# Write png from fits version of image, in same directory.
def write_png(fname, plt):
    '''Write PNG version of image from fits file.'''
    plt.figure(1)
    plt.imshow(fits.getdata(fname), aspect='auto', origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig(f'{fname.split(".fit")[0]}.png', bbox_inches='tight')
    plt.close(1)

@click.command()
@click.option('--setup', '-s', is_flag=True,
    help='Sets up the directory for reduction. Use this option only the first '
         'time astro_reduce is run in the directory.')
@click.option('--interpolate', '-i', is_flag=True,
    help='Interpolate existing dark field images if some are missing.')
@click.option('--verbose', '-v', is_flag=True,
    help='Enables verbose mode (recommended).')
@click.option('--tmppng', '-t', is_flag=True,
    help='Write PNG format of intermediary images after reduction.')
@click.option('--redpng', '-r', is_flag=True,
    help='Write PNG format of reduced images after reduction.')
def cli(setup, interpolate, verbose, tmppng, redpng):
    '''Reduce CCD images from objects with flat and dark field images.'''
    # Initialize globals
    conf_file_name = f'{getcwd().split("/")[-1]}.json'

    # If setup option is on, set up the directory for reduction
    if setup:
        click.echo('Setting up for reduction.')
        # Initialize objects, exposure, filters lists, and working directories
        objects = list()
        exposures = list()
        filters = list()
        if not exists(OBJ):
            mkdir(OBJ)
        if not exists(DARK):
            mkdir(DARK)
        if not exists(FLAT):
            mkdir(FLAT)

        # Open all images, retrieve exposures, filters, etc. and copy files to
        # astro_reduce working directories. This way the images are backed-up
        # at the same time
        if verbose:
            click.echo('Copying dark field images... ', nl=False)
        for file in glob(f'{BDARK}/*.fit*'):
            _, fn = d_read(file)
            copy(file, f'{DARK}/{fn}')
        if verbose:
            click.echo('Done.')

        if verbose:
            click.echo('Copying flat field images... ', nl=False)
        for file in glob(f'{BFLAT}/*.fit*'):
            _, _, fn = f_read(file)
            copy(file, f'{FLAT}/{fn}')
        if verbose:
            click.echo('Done.')

        if verbose:
            click.echo('Copying object images... ', nl=False)
        for file in glob(f'{BOBJ}/*.fit*'):
            obj, fil, exp, fn = o_read(file)
            objects.append(obj)
            filters.append(fil)
            exposures.append(exp)
            copy(file, f'{OBJ}/{fn}')
        if verbose:
            click.echo('Done.')

        # End up the setup by writing the configuration file.
        click.echo(f'Writing configuration file to {conf_file_name}.')
        write_conf_file(list(set(objects)), list(set(exposures)),
                        list(set(filters)), conf_file_name)

    # Parse configuration file to obtain configuration dictionary.
    click.echo(f'Parsing configuration file {conf_file_name}.')
    with open(conf_file_name, 'r') as cfile:
        conf_dic = loads(cfile.read())

    # Obtain list of all object, dark, flat field file names.
    object_files = dict(
        [(obj, glob(f'{OBJ}/{obj}_*.fit*'))
            for obj in conf_dic['objects']])
    dark_files = dict(
        [(exp, glob(f'{DARK}/{di}_{exp}_*.fit*'))
            for exp in conf_dic['exposures']])
    flat_files = dict(
        [(filt, glob(f'{FLAT}/{fi}_{filt}_*.fit*'))
            for filt in conf_dic['filters']])

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
            click.echo('If you want to interpolate the missing dark fields '
                       'from the existing ones, use the `--interpolate` '
                       'option.')
            click.echo('Exiting.')
            exit(1)
    for key in flat_files:
        if not flat_files[key]:
            click.echo(f'Did not find files for the {key} filter flats.')
            click.echo(f'They should be in the `{FLAT}` directory.')
            click.echo('Exiting.')
            exit(1)

    # Report all files found.
    reg = re.compile(r'_\d*\.')
    tstring = '****** {:25} ******'
    sstring = '    {:8}: {}'
    if verbose:
        click.echo('Files found:')
        click.echo(tstring.format(f' Objects (`{OBJ}`) '))
        for obj in conf_dic['objects']:
            uniq_names = set([reg.sub('_*.', basename(name))
                              for name in object_files[obj]])
            click.echo(sstring.format(obj, uniq_names))

        click.echo(tstring.format(f' Dark fields (`{DARK}`) '))
        for exp in conf_dic['exposures']:
            uniq_names = set([reg.sub('_*.', basename(name))
                              for name in dark_files[exp]]) or None
            click.echo(sstring.format(f'{exp}ms', uniq_names))

        click.echo(tstring.format(f' Flat fields (`{FLAT}`) '))
        for filt in conf_dic['filters']:
            uniq_names = set([reg.sub('_*.', basename(name))
                              for name in flat_files[filt]])
            click.echo(sstring.format(filt, uniq_names))

    # STEP 0: Create directory for tmp and reduced images if not existent.
    if verbose:
        click.echo('Creating folders for reduced and intermediary images.')
    if not exists(RED):
        mkdir(RED)
    if not exists(TMP):
        mkdir(TMP)

    # STEP 1: Write the master dark files (medians of darks)
    # for each available exposure.
    click.echo('Writing master dark images.')
    all_exposures = conf_dic['exposures']
    available_exposures = [exp for exp in dark_files if dark_files[exp]]
    for exp in available_exposures:
        if verbose:
            click.echo(f'    {exp}ms... ', nl=False)
        mdark_data = np.median(np.array([fits.getdata(fitsfile)
                                         for fitsfile in dark_files[exp]]),
                               axis=0)
        mdark_header = fits.getheader(dark_files[exp][0])
        fits.writeto(f'{TMP}/mdark_{exp}.fits', mdark_data, mdark_header,
                     overwrite=True)
        if verbose:
            click.echo('Done.')

    # STEP 1.5: If there are some missing darks and the interpolate option
    # is on, then interpolate the master darks.
    # We use least squares linear interpolation, i.e., we calculate `a` and `b`
    # such that (missing_dark) = a * (exposure_time) + b.
    # Exit if there are no darks at all.
    if not available_exposures:
        click.echo('There are no dark files at all! Cannot interpolate...')
        click.echo('Exiting.')
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
        mxy = np.mean([float(exp) * fits.getdata(f'{TMP}/mdark_{exp}.fits')
                       for exp in available_exposures], axis=0)
        mx = np.mean([float(exp) for exp in available_exposures])
        my = np.mean([fits.getdata(f'{TMP}/mdark_{exp}.fits')
                      for exp in available_exposures], axis=0)
        mx2 = np.mean([float(exp) ** 2 for exp in available_exposures])

        # a and b
        a = (mxy - mx * my) / (mx2 - mx ** 2)
        b = my - mx * a

    # Write all the missing master darks!
    click.echo('Interpolating missing master dark images.')
    for exp in list(set(all_exposures) - set(available_exposures)):
        if verbose:
            click.echo(f'    {exp}ms... ', nl=False)
        new_mdark_data = float(exp) * a + b
        fits.writeto(f'{TMP}/mdark_{exp}.fits', new_mdark_data,
                     overwrite=True)
        if verbose:
            click.echo('Done.')

    # STEP 2: Write master transmission files for each filter:
    # (median of flats - dark of same exposure) normalized.
    # Handy function to extract exposure from flat file name.
    click.echo('Calculating master transmission images.')
    fexp = lambda fname: fname.split('.fit')[0].split('_')[-2]
    for filt in flat_files:
        if verbose:
            click.echo(f'    {filt}... ', nl=False)
        # Use first file of series to get header, where we set exposure to -1
        # because the master transmission may be calculated using various
        # exposure times.
        mflat_header = fits.getheader(flat_files[filt][0])
        mflat_header['EXPTIME'] = '-1'
        mflat_header['EXPOSURE'] = '-1'

        # Calculate normalized flats.
        normalized_flats = list()
        for fitsfile in flat_files[filt]:
            tmp = fits.getdata(fitsfile) \
                - fits.getdata(f'{TMP}/mdark_{fexp(fitsfile)}.fits')
            normalized_flats.append(tmp / tmp.mean(axis=0))

        mtrans_data = np.median(normalized_flats, axis=0)
        fits.writeto(f'{TMP}/mtrans_{filt}.fits', mtrans_data, mflat_header,
                     overwrite=True)
        if verbose:
            click.echo('Done.')

    # STEP 3: Reduce all the object images with corresponding filter mflat
    # and exposure mdark.
    click.echo('Writing auxiliary object images.')
    for obj in object_files:
        if verbose:
            click.echo(f'    {obj}... ', nl=False)
        for fname in object_files[obj]:
            bfname = basename(fname)
            filt, exp = fname_bits(bfname)

            # Corresponding darks and flats
            mdark_data = fits.getdata(f'{TMP}/mdark_{exp}.fits')
            mtrans_data = fits.getdata(f'{TMP}/mtrans_{filt}.fits')

            # (Raw - Dark) / Trans
            aux_data = (fits.getdata(fname) - mdark_data) / mtrans_data
            aux_header = fits.getheader(fname)
            fits.writeto(f'{TMP}/{bfname.split(".fit")[0]}_{AUX}.fits',
                         aux_data, aux_header, overwrite=True)
        if verbose:
            click.echo('Done.')

    # STEP 4: Within series, realign the aux images and write the median images
    # of those. You are left with one image per object per filter per exposure.
    click.echo('Realigning object images.')
    for obj in object_files:
        if verbose:
            click.echo(f'    {obj}:')
        # Group all the object files by *tag*, i.e. by filter, exposure.
        name_tag_hash = [(basename(fname), f'{fname_bits(basename(fname))}')
                         for fname in object_files[obj]]
        names_per_tag = defaultdict(list)
        for name, tag in name_tag_hash:
            names_per_tag[tag].append(name)

        # Now you align images which have the same tag.
        for tag in names_per_tag:
            # Rebuild series, filter and exposure from tag (they are those of,
            # e.g., the first name in the list.)
            f, e = fname_bits(names_per_tag[tag][0])
            if verbose:
                click.echo(f'      {f}/{e}... ', nl=False)
            # Calculate aligned and medianed image
            # from all images with same tag.
            aux_files = glob(f'{TMP}/{obj}_{f}_{e}_*_{AUX}.fits')
            reduced_data = align_and_median(aux_files)
            reduced_header = fits.getheader(f'{OBJ}/{names_per_tag[tag][0]}')

            fits.writeto(f'{RED}/{obj}_{f}_{e}.fits', reduced_data,
                         reduced_header, overwrite=True)
            if verbose:
                click.echo('Done.')

    # STEP 5: If options redpng or tmppng are on, write
    # PNG versions of all the tmp and reduced images.
    if redpng or tmppng:
        import matplotlib.pyplot as plt

    if redpng:
        click.echo('Writing PNG versions of reduced images... ', nl=False)
        for ffile in glob(f'{RED}/*.fits'):
            write_png(ffile, plt)
        click.echo('Done.')

    if tmppng:
        click.echo('Writing PNG versions of intermediate images... ', nl=False)
        for ffile in glob(f'{TMP}/*.fits'):
            write_png(ffile, plt)
        click.echo('Done.')

    click.echo('All done.')
    # ALL DONE.
