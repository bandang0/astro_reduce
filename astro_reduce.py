'''astro_reduce -- A Simple CCD Images Reducer for the Paris Observatory.'''

from sys import exit
from os.path import basename, exists, getsize
from os import mkdir, getcwd, listdir
from shutil import copy, rmtree
from glob import glob
from json import loads, dump, decoder
from collections import defaultdict
from re import compile, sub
from hashlib import md5
from packaging.version import parse

import click
import numpy as np
from astropy.io import fits
import scipy
from scipy.signal import fftconvolve


# Comment for header keywords
hc = 'Exposure time in seconds'

# Paths and extensions.
# User image directories.
UDARK = 'DARK'
UFLAT = 'FLAT'
UOBJ = 'ORIGINAL'

# Astro-reduce working directories.
OBJ = 'ar_objects'
DARK = 'ar_darks'
FLAT = 'ar_flats'
TMP = 'tmp'

# File names and extensions.
di = 'dark'
fi = 'flat'
RED = 'reduced'
AUX = 'aux'

# Simple hashing function for file names.
hsh = lambda x: md5(x.encode('utf-8')).hexdigest()

# Read data from list of files and return aligned and meaned version.
def align_and_median(infiles):
    '''Return the median of the re-aligned images from list of file names.'''
    if len(infiles) == 1:
        return fits.getdata(infiles[0])

    # Collect arrays and crosscorrelate all with the first (except the first).
    images = [fits.getdata(_) for _ in infiles]
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
        raise IOError('No exposure keyword in header of `{}`.'.format(fname))
    return exp, '{}_{}_{}.fits'.format(di, exp, hsh(fname))

# Return filter, exposure and custom file name for a flat field image
def f_read(fname):
    '''Return filter, exposure and custom file name for a flat field image.'''
    head = fits.getheader(fname)
    # Filter.
    if 'FILTER' in head.keys():
        fil = sub('[- _]', '', head['FILTER'])
    else:
        raise IOError('No filter keyword in header of `{}`.'.format(fname))

    # Exposure.
    if 'EXPTIME' in head.keys():
        exp = int(1000 * head['EXPTIME'])
    elif 'EXPOSURE' in head.keys():
        exp = int(1000 * head['EXPOSURE'])
    else:
        raise IOError('No exposure keyword in header of `{}`.'.format(fname))

    return fil, exp, '{}_{}_{}_{}.fits'.format(fi, fil, exp, hsh(fname))

# Return object, filter, exposure and custom file name for an object image.
def o_read(fname):
    '''Return object, filter, exposure and custom file name for object image.'''
    head = fits.getheader(fname)
    # Object.
    if 'OBJECT' in head.keys():
        obj = sub('[ _]', '-', head['OBJECT'])
    else:
        raise IOError('No object keyword in header of `{}`.'.format(fname))

    # Filter.
    if 'FILTER' in head.keys():
        fil = sub('[- _]', '', head['FILTER'])
    else:
        raise IOError('No filter keyword in header of `{}`.'.format(fname))

    # Exposure.
    if 'EXPTIME' in head.keys():
        exp = int(1000 * head['EXPTIME'])
    elif 'EXPOSURE' in head.keys():
        exp = int(1000 * head['EXPOSURE'])
    else:
        raise IOError('No exposure keyword in header of `{}`.'.format(fname))

    return obj, fil, exp, '{}_{}_{}_{}.fits'.format(obj, fil, exp, hsh(fname))

# Write the configuration file for images in current directory.
def write_conf_file(objects, exposures, filters, conf_file_name):
    '''Write the configuration file from list of objects, exposures, filters.'''
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
    plt.savefig('{}.png'.format(fname.split(".fit")[0]), bbox_inches='tight')
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
    # Before all things, check version of scipy.
    if not (parse(scipy.__version__) < parse("1.4.1")):
        click.echo('E: scipy version {} detected.\n'
                   'E: This could cause crashes during image reduction.\n'
                   'E: Please install scipy <= 1.3.3 for astro_reduce to run\n'
                   'E: correctly.'.format(scipy.__version__))
        exit(1)

    # Initialize globals
    conf_file_name = '{}.json'.format(getcwd().split("/")[-1])

    # If setup option is on, set up the directory for reduction
    if setup:
        click.echo('Setting up for reduction.')
        # Make sure the user image folders are there:

        if not (exists(UOBJ) and exists(UFLAT) and exists(UDARK)):
            click.echo('E: Could not find the raw images folders `DARK`,\n'
                       'E: `FLAT` or `ORIGINAL`. They should contain the\n'
                       'E: images to be reduced. Refer to the documentation\n'
                       'E: for details.')
            exit(1)

        # Initialize objects, exposure, filters lists, and working directories
        objects = list()
        exposures = list()
        filters = list()
        if exists(OBJ):
            rmtree(OBJ, ignore_errors=True)
        mkdir(OBJ)
        if exists(DARK):
            rmtree(DARK, ignore_errors=True)
        mkdir(DARK)
        if exists(FLAT):
            rmtree(FLAT, ignore_errors=True)
        mkdir(FLAT)

        # Open all images, retrieve exposures, filters, etc. and copy files to
        # astro_reduce working directories. This way the images are backed-up
        # at the same time
        if verbose:
            click.echo('Copying dark field images... ', nl=False)
        for file in glob('{}/*.fit*'.format(UDARK)):
            exp, fn = d_read(file)
            exposures.append(exp)
            copy(file, '{}/{}'.format(DARK, fn))
        if verbose:
            click.echo('Done.')

        if verbose:
            click.echo('Copying flat field images... ', nl=False)
        for file in glob('{}/*.fit*'.format(UFLAT)):
            fil, exp, fn = f_read(file)
            exposures.append(exp)
            filters.append(fil)
            copy(file, '{}/{}'.format(FLAT, fn))
        if verbose:
            click.echo('Done.')

        if verbose:
            click.echo('Copying object images... ', nl=False)
        for file in glob('{}/*.fit*'.format(UOBJ)):
            obj, fil, exp, fn = o_read(file)
            objects.append(obj)
            filters.append(fil)
            exposures.append(exp)
            copy(file, '{}/{}'.format(OBJ, fn))
        if verbose:
            click.echo('Done.')

        # End up the setup by writing the configuration file.
        click.echo('Writing configuration data to `{}`.'.format(conf_file_name))
        write_conf_file(list(set(objects)), list(set(exposures)),
                        list(set(filters)), conf_file_name)

    # Parse configuration file to obtain configuration dictionary.
    try:
        click.echo('Parsing configuration file `{}`.'.format(conf_file_name))
        with open(conf_file_name, 'r') as cfile:
            conf_dic = loads(cfile.read())
    except FileNotFoundError:
        click.echo('E: Configuration file `{}` not found.\n'
                   'E: If this is the first time you run astro_reduce in this\n'
                   'E: directory, use the `--setup` option to setup the\n'
                   'E: reduction and generate a configuration file.'
                   ''.format(conf_file_name))
        exit(1)
    except decoder.JSONDecodeError:
        click.echo('E: Unable to parse configuration file `{}`.\n'
                   'E: Fix by rerunning astro_reduce with the --setup option.'
                   ''.format(conf_file_name))
        exit(1)


    # Obtain list of all object, dark, flat field file names.
    object_files = dict(
        [(obj, glob('{}/{}_*.fit*'.format(OBJ, obj)))
            for obj in conf_dic['objects']])
    dark_files = dict(
        [(exp, glob('{}/{}_{}_*.fit*'.format(DARK, di, exp)))
            for exp in conf_dic['exposures']])
    flat_files = dict(
        [(filt, glob('{}/{}_{}_*.fit*'.format(FLAT, fi, filt)))
            for filt in conf_dic['filters']])

    # Check working directories are still there.
    if not (exists(OBJ) and exists(FLAT) and exists(DARK)):
        click.echo('E: Seems like astro_reduce\'s working folders\n'
                   'E: (those starting with `ar_`) were removed.\n'
                   'E: Please rerun astro_reduce with `--setup` option.')
        exit(1)

    # Check all images are same size (if not we'll have a problem).
    if len(set(map(getsize,
        glob('{}/*'.format(DARK))
        + glob('{}/*'.format(FLAT))
        + glob('{}/*'.format(OBJ))))) != 1:
        click.echo('E: Seems like all image files don\'t have the same size.\n'
                   'E: Please remove relevant files and all `ar_` folders,\n'
                   'E: and rerun astro_reduce with `--setup` option.')
        exit(1)

    # Check if files exist.
    for key in object_files:
        if not object_files[key]:
            click.echo('E: Did not find files for {} object.'.format(key))
            exit(1)
    for key in dark_files:
        if not dark_files[key] and not interpolate:
            # If the interpolate option is off and there are some darks
            # missing, exit.
            click.echo('E: Did not find dark field images for {} exposure.\n'
                       'E: If you want to interpolate the missing dark fields\n'
                       'E: from the other ones, rerun using the\n'
                       'E: `--interpolate` option.'.format(key))
            exit(1)
    for key in flat_files:
        if not flat_files[key]:
            click.echo('E: No flat field images for {} filter.'.format(key))
            exit(1)

    # Report all files found.
    reg = compile(r'_[a-z0-9]*\.')
    tstring = '****** {:25} ******'
    sstring = '    {:8}: {}'
    if verbose:
        handy = lambda x: reg.sub('_*.', basename(x))
        click.echo('Files found:')
        click.echo(tstring.format(' Objects (`{}`) '.format(OBJ)))
        for obj in conf_dic['objects']:
            uniq_names = set(map(handy, object_files[obj]))
            click.echo(sstring.format(obj, uniq_names))

        click.echo(tstring.format(' Dark fields (`{}`) '.format(DARK)))
        for exp in conf_dic['exposures']:
            uniq_names = set(map(handy, dark_files[exp])) or None
            click.echo(sstring.format('{}'.format(exp), uniq_names))

        click.echo(tstring.format(' Flat fields (`{}`) '.format(FLAT)))
        for filt in conf_dic['filters']:
            uniq_names = set(map(handy, flat_files[filt]))
            click.echo(sstring.format(filt, uniq_names))

    # STEP 0: Create directory for tmp and reduced images if not existent.
    if verbose:
        click.echo('Creating folders for reduced and intermediary images.')
    if exists(RED):
        rmtree(RED, ignore_errors=True)
    mkdir(RED)
    if exists(TMP):
        rmtree(TMP, ignore_errors=True)
    mkdir(TMP)

    # STEP 1: Write the master dark files (medians of darks)
    # for each available exposure.
    click.echo('Writing master dark images.')
    all_exposures = conf_dic['exposures']
    available_exposures = [exp for exp in dark_files if dark_files[exp]]
    for exp in available_exposures:
        if verbose:
            click.echo('    {}... '.format(exp), nl=False)
        mdark_data = np.median([fits.getdata(_) for _ in dark_files[exp]],
                               axis=0)
        mdark_header = fits.getheader(dark_files[exp][0])

        # Write fits file and header.
        nname = '{}/mdark_{}.fits'.format(TMP, exp)
        fits.writeto(nname, mdark_data, mdark_header, overwrite=True)
        fits.setval(nname, 'FILTER', value='        ')
        fits.setval(nname, 'IMAGETYP', value='Dark    ')
        fits.setval(nname, 'EXPTIME', value=float(exp / 1000.), comment=hc)
        fits.setval(nname, 'EXPOSURE', value=float(exp / 1000.), comment=hc)
        fits.setval(nname, 'OBJECT', value='DARK    ')

        if verbose:
            click.echo('Done ({} images).'.format(len(dark_files[exp])))

    # STEP 1.5: If there are some missing darks and the interpolate option
    # is on, then interpolate the master darks.
    # We use least squares linear interpolation, i.e., we calculate `a` and `b`
    # such that (missing_dark) = a * (exposure_time) + b.
    # Exit if there are no darks at all.
    if not available_exposures:
        click.echo('E: There are no dark files at all! Cannot interpolate...')
        exit(1)

    if len(available_exposures) == 1:
        # If there's only one, a = (only_dark) / (its exposure), b = 0
        only_exp = available_exposures[0]
        only_mdark = fits.getdata('{}/mdark_{}.fits'.format(TMP, only_exp))
        a = only_mdark / float(only_exp)
        b = np.zeros_like(only_mdark)
    else:
        # If not, if you want to fit y = a * x + b,
        # then the LS solution is:
        # a = (<xy> - <x><y>) / (<x ** 2> - <x> ** 2)
        # b = <y> - a * <x>
        mxy = np.mean([float(exp)
                       * fits.getdata('{}/mdark_{}.fits'.format(TMP, exp))
                       for exp in available_exposures], axis=0)
        mx = np.mean([float(exp) for exp in available_exposures])
        my = np.mean([fits.getdata('{}/mdark_{}.fits'.format(TMP, exp))
                      for exp in available_exposures], axis=0)
        mx2 = np.mean([float(exp) ** 2 for exp in available_exposures])

        # a and b
        a = (mxy - mx * my) / (mx2 - mx ** 2)
        b = my - mx * a

    # Write all the missing master darks!
    click.echo('Interpolating missing master dark images.')
    for exp in list(set(all_exposures) - set(available_exposures)):
        if verbose:
            click.echo('    {}... '.format(exp), nl=False)
        new_mdark_data = float(exp) * a + b

        # Write fits file and header.
        nname = '{}/mdark_{}.fits'.format(TMP, exp)
        fits.writeto(nname, new_mdark_data, overwrite=True)
        fits.setval(nname, 'FILTER', value='        ')
        fits.setval(nname, 'IMAGETYP', value='Interpolated dark')
        fits.setval(nname, 'EXPTIME', value=float(exp / 1000.), comment=hc)
        fits.setval(nname, 'EXPOSURE', value=float(exp / 1000.), comment=hc)
        fits.setval(nname, 'OBJECT', value='DARK    ')
        if verbose:
            click.echo('Done.')

    # STEP 2: Write master transmission files for each filter:
    # (median of flats - dark of same exposure) normalized.
    # Handy function to extract exposure from flat file name.
    click.echo('Calculating master transmission images.')
    fexp = lambda fname: fname.split('.fit')[0].split('_')[-2]
    for filt in flat_files:
        if verbose:
            click.echo('    {}... '.format(filt), nl=False)

        # Calculate normalized flats.
        normalized_flats = list()
        for fitsfile in flat_files[filt]:
            tmp = fits.getdata(fitsfile) \
                - fits.getdata('{}/mdark_{}.fits'.format(TMP, fexp(fitsfile)))
            normalized_flats.append(tmp / tmp.mean(axis=0))

        mtrans_data = np.median(normalized_flats, axis=0)
        mflat_header = fits.getheader(flat_files[filt][0])

        # Write fits file and header.
        nname = '{}/mtrans_{}.fits'.format(TMP, filt)
        fits.writeto(nname, mtrans_data, mflat_header, overwrite=True)
        fits.setval(nname, 'FILTER', value=filt)
        fits.setval(nname, 'IMAGETYP', value='Light Frame')
        fits.setval(nname, 'EXPTIME', value='-1', comment=hc)
        fits.setval(nname, 'EXPOSURE', value='-1', comment=hc)
        fits.setval(nname, 'OBJECT', value='FLAT    ')

        if verbose:
            click.echo('Done ({} images).'.format(len(flat_files[filt])))

    # STEP 3: Reduce all the object images with corresponding filter mflat
    # and exposure mdark.
    click.echo('Writing auxiliary object images.')
    for obj in object_files:
        if verbose:
            click.echo('    {}... '.format(obj), nl=False)
        for fname in object_files[obj]:
            bfname = basename(fname)
            filt, exp = fname_bits(bfname)

            # Corresponding darks and flats
            mdark_data = fits.getdata('{}/mdark_{}.fits'.format(TMP, exp))
            mtrans_data = fits.getdata('{}/mtrans_{}.fits'.format(TMP, filt))

            # (Raw - Dark) / Trans
            aux_data = (fits.getdata(fname) - mdark_data) / mtrans_data
            aux_header = fits.getheader(fname)

            # Write fits file and header.
            nname = '{}/{}_{}.fits'.format(TMP, bfname.split(".fit")[0], AUX)
            fits.writeto(nname, aux_data, aux_header, overwrite=True)
            fits.setval(nname, 'FILTER', value=filt)
            fits.setval(nname, 'IMAGETYP', value='Light Frame')
            fits.setval(nname, 'EXPTIME', value=float(exp) / 1000., comment=hc)
            fits.setval(nname, 'EXPOSURE', value=float(exp) / 1000., comment=hc)
            fits.setval(nname, 'OBJECT', value=obj)
        if verbose:
            click.echo('Done.')

    # STEP 4: For all objects realign the aux images and write the median images
    # of those. You are left with one image per object per filter per exposure.
    click.echo('Realigning object images.')
    for obj in object_files:
        if verbose:
            click.echo('    {}:'.format(obj))
        # Group all the object files by *tag*, i.e. by filter, exposure.
        name_tag_hash = [(basename(fname),
                         '{}'.format(fname_bits(basename(fname))))
                         for fname in object_files[obj]]
        names_per_tag = defaultdict(list)
        for name, tag in name_tag_hash:
            names_per_tag[tag].append(name)

        # Now you align images which have the same tag.
        for tag in names_per_tag:
            # Rebuild filter and exposure from tag (they are those of,
            # e.g., the first name in the list.)
            f, e = fname_bits(names_per_tag[tag][0])
            if verbose:
                click.echo('      {}/{}... '.format(f, e), nl=False)
            # Calculate aligned and medianed image
            # from all images with same tag.
            aux_files = glob('{}/{}_{}_{}_*_{}.fits'.format(TMP, obj, f, e,AUX))
            reduced_data = align_and_median(aux_files)
            reduced_header = fits.getheader(
                                     '{}/{}'.format(OBJ, names_per_tag[tag][0]))

            # Write fits file and header.
            nname = '{}/{}_{}_{}.fits'.format(RED, obj, f, e)
            fits.writeto(nname, reduced_data, reduced_header, overwrite=True)
            fits.setval(nname, 'FILTER', value=f)
            fits.setval(nname, 'IMAGETYP', value='Light Frame')
            fits.setval(nname, 'EXPTIME', value=float(e) / 1000., comment=hc)
            fits.setval(nname, 'EXPOSURE', value=float(e) / 1000., comment=hc)
            fits.setval(nname, 'OBJECT', value=obj)
            if verbose:
                click.echo('Done ({} images).'.format(len(aux_files)))

    # STEP 5: If options redpng or tmppng are on, write
    # PNG versions of all the tmp and reduced images.
    if redpng or tmppng:
        import matplotlib.pyplot as plt

    if redpng:
        click.echo('Writing PNG versions of reduced images... ', nl=False)
        for ffile in glob('{}/*.fits'.format(RED)):
            write_png(ffile, plt)
        click.echo('Done.')

    if tmppng:
        click.echo('Writing PNG versions of intermediate images... ', nl=False)
        for ffile in glob('{}/*.fits'.format(TMP)):
            write_png(ffile, plt)
        click.echo('Done.')

    click.echo('All done.')
    # ALL DONE.
