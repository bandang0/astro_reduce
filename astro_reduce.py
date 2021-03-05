'''astro_reduce -- A Simple CCD Images Reducer for the Paris Observatory.'''


from collections import defaultdict
from datetime import timedelta
from glob import glob
from json import loads, decoder
from os.path import basename, exists, getsize
from os import mkdir, getcwd, remove
from re import compile, sub
from shutil import copy, rmtree
from sys import exit
from time import time

import click
import numpy as np
from matplotlib import pyplot as plt

from cosmetic import align_and_median
from helpers import *


@click.command()
@click.version_option()
@click.option('--setup', '-s', is_flag=True,
              help='Set up the directory for reduction. Use this option the '
              'first time astro_reduce is run in the directory or after the '
              '--clear option was used.')
@click.option('--clear', '-c', is_flag=True,
              help='Remove all astro_reduce-related files and folders in '
              'current directory and exit.')
@click.option('--interpolate', '-i', is_flag=True,
              help='Interpolate existing dark fields if some are missing.')
@click.option('--verbose', '-v', is_flag=True,
              help='Enables verbose mode (recommended).')
@click.option('--tmppng', '-t', is_flag=True,
              help='Write PNG format of intermediate images after reduction.')
@click.option('--redpng', '-r', is_flag=True,
              help='Write PNG format of reduced images after reduction.')
def cli(setup, clear, interpolate, verbose, tmppng, redpng):
    '''Main run of astro_reduce:

    If --setup option is on: copy all user fits data to files with standard
       names in the working directories, and exit.

    If --clear option is on: clear directory of all astro_reduce working data
    and exit.

    Else, reduction procedure:
    i) Reduce dark field images to master darks for all exposures.
    ii) If necessary, interpolate darks for exposures lacking dark fields.
    iii) Reduce all flat field images to master transmission files.
    iv) Reduce all object images with the master reduction files.
    v) Realign the object images of same series.

    If --{tmp,red}png options are on: generate PNG versions of temporary
    and reduced images.

    '''
    # Current working directory.
    cwd = basename(getcwd())

    # Welcome.
    click.secho('    Welcome to astro_reduce!\n'
                '    Software is copyright 2018-2021 R. Duque.\n\n',
                fg='cyan', bold=True)
    click.secho('Currently working in directory `{}`.'.format(cwd), fg='green')

    # Initialize configuration file name and timer.
    conf_file_name = '{}.json'.format(cwd)
    t0 = time()

    # If clear option is on, remove all files and folders and exit.
    if clear:
        click.secho('Clearing directory of astro_reduce files and folders...',
                    fg='green', nl=False)
        if exists(conf_file_name):
            remove(conf_file_name)
        for folder in [OBJ, FLAT, DARK, RED, TMP]:
            if exists(folder):
                rmtree(folder, ignore_errors=True)
        click.secho(' Done.', fg='green')
        exit(0)

    # If setup option is on, set up the directory for reduction.
    if setup:
        click.secho('Setting up for reduction:', fg='green')
        # Make sure the user image folders are there:
        if not (exists(UOBJ) and exists(UFLAT) and exists(UDARK)):
            click.secho('E: Did not find folder `DARK`, `FLAT` or `ORIGINAL`\n'
                        'E: containing the raw images to be reduced.\n'
                        'E: Refer to the documentation for details.', fg='red')
            exit(1)

        # Initialize objects, exposure, filters lists, and working directories.
        objects = list()
        exposures = list()
        filters = list()
        for folder in [OBJ, DARK, FLAT]:
            if exists(folder):
                rmtree(folder, ignore_errors=True)
            mkdir(folder)

        # Open all images, retrieve exposures, filters, etc. and copy files to
        # astro_reduce working directories. This way the images are backed-up
        # at the same time.
        if verbose:
            click.secho('  Copying dark field images...', nl=False)
        for file in glob('{}/*.fit*'.format(UDARK)):
            exp, fn = dark_read_header(file)
            exposures.append(exp)
            copy(file, '{}/{}'.format(DARK, fn))
        if verbose:
            click.secho('     Done.')

        if verbose:
            click.echo('  Copying flat field images...', nl=False)
        for file in glob('{}/*.fit*'.format(UFLAT)):
            fil, exp, fn = flat_read_header(file)
            exposures.append(exp)
            filters.append(fil)
            copy(file, '{}/{}'.format(FLAT, fn))
        if verbose:
            click.echo('     Done.')

        if verbose:
            click.echo('  Copying object images...', nl=False)
        for file in glob('{}/*.fit*'.format(UOBJ)):
            obj, fil, exp, fn = obj_read_header(file)
            objects.append(obj)
            filters.append(fil)
            exposures.append(exp)
            copy(file, '{}/{}'.format(OBJ, fn))
        if verbose:
            click.echo('         Done.')

        # End up the setup by writing the configuration file.
        if verbose:
            click.echo('  Writing configuration file `{}`.'.format(
                                                                conf_file_name))
        write_conf_file(list(set(objects)), list(set(exposures)),
                        list(set(filters)), conf_file_name)

        click.secho('Done.', fg='green')
        exit(0)

    # Parse configuration file to obtain configuration dictionary.
    try:
        click.secho('Parsing configuration file `{}`.'.format(conf_file_name),
                    fg='green')
        with open(conf_file_name, 'r') as cfile:
            conf_dic = loads(cfile.read())
    except FileNotFoundError:
        click.secho('E: Configuration file `{}` not found.\n'
                    'E: If it is the first time you run astro_reduce in this\n'
                    'E: directory, use the --setup option to setup the\n'
                    'E: reduction and generate a configuration file.'
                    ''.format(conf_file_name), fg='red')
        exit(1)
    except decoder.JSONDecodeError:
        click.secho('E: Unable to parse configuration file `{}`.\n'
                    'E: Fix by rerunning astro_reduce with the --setup option.'
                    ''.format(conf_file_name), fg='red')
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
        click.secho('E: Seems like astro_reduce\'s working folders\n'
                    'E: (those starting with `ar_`) were removed.\n'
                    'E: Please rerun astro_reduce with the --setup option.',
                    fg='red')
        exit(1)

    # Check all images are same size (if not we'll have a problem).
    if len(set(map(getsize,
                   glob('{}/*'.format(DARK))
                   + glob('{}/*'.format(FLAT))
                   + glob('{}/*'.format(OBJ))))) != 1:
        click.secho('E: Seems like all image files don\'t have the same size.\n'
                    'E: Please remove offending files and rerun astro_reduce\n'
                    'E: with the --setup option.', fg='red')
        exit(1)

    # Check if files exist.
    for key in object_files:
        if not object_files[key]:
            click.secho('E: Did not find files for {} object.'.format(key),
                        fg='red')
            exit(1)
    for key in dark_files:
        if not dark_files[key] and not interpolate:
            # If the interpolate option is off and there are some darks
            # missing, exit.
            click.secho('E: Did not find dark field images for {}ms exposure.\n'
                        'E: In order to interpolate the missing dark fields\n'
                        'E: from the ones available, rerun using the\n'
                        'E: --interpolate option.'.format(key), fg='red')
            exit(1)
    for key in flat_files:
        if not flat_files[key]:
            click.secho('E: No flat field images for {} filter.'.format(key),
                        fg='red')
            exit(1)

    # Report all files found.
    reg = compile(r'_[a-z0-9]*\.')
    if verbose:
        handy = lambda x: reg.sub('_*.', basename(x))
        click.secho('Files found:', fg='green')
        click.secho('  Objects (`{}`):'.format(OBJ), fg='blue')
        for obj in conf_dic['objects']:
            uniq_names = set(map(handy, object_files[obj]))
            click.echo('    {:10}: {}'.format(obj, uniq_names))

        click.secho('  Dark fields (`{}`):'.format(DARK), fg='blue')
        for exp in conf_dic['exposures']:
            uniq_names = set(map(handy, dark_files[exp])) or None
            click.echo('    {:10}: {}'.format('{}ms'.format(exp), uniq_names))

        click.secho('  Flat fields (`{}`):'.format(FLAT), fg='blue')
        for filt in conf_dic['filters']:
            uniq_names = set(map(handy, flat_files[filt]))
            click.echo('    {:10}: {}'.format(filt, uniq_names))

    # STEP 0: Create directory for tmp and reduced images if not existent.
    if verbose:
        click.secho('Creating folders to hold reduced and intermediate images.',
                    fg='green')
    if exists(RED):
        rmtree(RED, ignore_errors=True)
    mkdir(RED)
    if exists(TMP):
        rmtree(TMP, ignore_errors=True)
    mkdir(TMP)

    # STEP 1: Write the master dark files (medians of darks)
    # for each available exposure.
    click.secho('Writing master dark images:', fg='green')
    all_exposures = conf_dic['exposures']
    available_exposures = [exp for exp in dark_files if dark_files[exp]]
    for exp in available_exposures:
        if verbose:
            click.echo('    {:14} '.format('{}ms...'.format(exp)), nl=False)
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
        click.secho('E: There are no dark files at all! Cannot interpolate...',
                    fg='red')
        exit(1)

    if len(available_exposures) == 1:
        # If there's only one available exosure time, consider that the darks
        # are dominated by the bias, which is likely. In this case:
        # a = 0, b = only_dark
        only_exp = available_exposures[0]
        only_mdark = fits.getdata('{}/mdark_{}.fits'.format(TMP, only_exp))
        a = np.zeros_like(only_mdark)
        b = only_mdark
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
    click.secho('Interpolating missing master dark images:', fg='green')
    for exp in list(set(all_exposures) - set(available_exposures)):
        if verbose:
            click.echo('    {:14} '.format('{}ms...'.format(exp)), nl=False)
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
    # mtrans = median(normalized(flat - dark of same exposure)).
    click.secho('Calculating master transmission images:', fg='green')
    # Handy function to extract exposure from flat file name.
    fexp = lambda fname: fname.split('.fit')[0].split('_')[-2]
    for filt in flat_files:
        if verbose:
            click.echo('    {:12}   '.format('{}...'.format(filt)), nl=False)

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
        fits.setval(nname, 'EXPTIME', value=-1., comment=hc)
        fits.setval(nname, 'EXPOSURE', value=-1., comment=hc)
        fits.setval(nname, 'OBJECT', value='FLAT    ')

        if verbose:
            click.echo('Done ({} images).'.format(len(flat_files[filt])))

    # STEP 3: Reduce all the object images with corresponding filter mtrans
    # and exposure mdark.
    click.secho('Writing auxiliary object images:', fg='green')
    for obj in object_files:
        if verbose:
            click.echo('    {:12}   '.format('{}...'.format(obj)), nl=False)
        for fname in object_files[obj]:
            bfname = basename(fname)
            filt, exp = fname_bits(bfname)

            # Corresponding darks and flats.
            mdark_data = fits.getdata('{}/mdark_{}.fits'.format(TMP, exp))
            mtrans_data = fits.getdata('{}/mtrans_{}.fits'.format(TMP, filt))

            # (Raw - Dark) / Trans.
            aux_data = (fits.getdata(fname) - mdark_data) / mtrans_data
            aux_header = fits.getheader(fname)

            # Write fits file and header.
            nname = '{}/{}_{}.fits'.format(TMP, bfname.split('.fit')[0], AUX)
            fits.writeto(nname, aux_data, aux_header, overwrite=True)
            fits.setval(nname, 'FILTER', value=filt)
            fits.setval(nname, 'IMAGETYP', value='Light Frame')
            fits.setval(nname, 'EXPTIME', value=float(exp) / 1000, comment=hc)
            fits.setval(nname, 'EXPOSURE', value=float(exp) / 1000, comment=hc)
            fits.setval(nname, 'OBJECT', value=obj)
        if verbose:
            click.echo('Done.')

    # STEP 4: For all objects realign and median the aux images.
    # You are left with one image per object per filter per exposure.
    click.secho('Realigning object images:', fg='green')
    for obj in object_files:
        if verbose:
            click.secho('  {}:'.format(obj), fg='blue')
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
                click.echo('    {:23}'.format('{}:{}ms...'.format(f, e)),
                           nl=False)
            # Calculate aligned and medianed image
            # from all images with same tag.
            aux_files = glob('{}/{}_{}_{}_*_{}.fits'.format(TMP, obj, f,
                                                            e, AUX))
            reduced_data = align_and_median(aux_files)
            reduced_header = fits.getheader('{}/{}'
                                            ''.format(OBJ,
                                                      names_per_tag[tag][0]))

            # Write fits file and header.
            nname = '{}/{}_{}_{}.fits'.format(RED, obj, f, e)
            fits.writeto(nname, reduced_data, reduced_header, overwrite=True)
            fits.setval(nname, 'FILTER', value=f)
            fits.setval(nname, 'IMAGETYP', value='Light Frame')
            fits.setval(nname, 'EXPTIME', value=float(e) / 1000., comment=hc)
            fits.setval(nname, 'EXPOSURE', value=float(e) / 1000., comment=hc)
            fits.setval(nname, 'OBJECT', value=obj)
            if verbose:
                click.echo('       Done ({} images).'.format(len(aux_files)))

    # STEP 5: If options redpng or tmppng are on, write
    # PNG versions of all the tmp and reduced images.
    if redpng:
        click.secho('Writing PNG versions of reduced images... ',
                    fg='green', nl=False)
        for ffile in glob('{}/*.fits'.format(RED)):
            write_png(ffile, plt)
        click.secho('Done.', fg='green')

    if tmppng:
        click.secho('Writing PNG versions of intermediate images... ',
                    fg='green', nl=False)
        for ffile in glob('{}/*.fits'.format(TMP)):
            write_png(ffile, plt)
        click.secho('Done.', fg='green')

    # Report execution time.
    t1 = time()
    click.secho('\nAll done. ({})'.format(timedelta(seconds=int(t1 - t0))),
                fg='green')
    # ALL DONE.
