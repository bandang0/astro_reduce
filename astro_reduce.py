'''astro_reduce -- A Simple CCD Images Reducer for the Paris Observatory.'''


from collections import defaultdict
from datetime import timedelta
from glob import glob
from json import loads, decoder
from os.path import basename, exists, getsize, splitext
from os import mkdir, getcwd, remove, system
from pkg_resources import resource_filename
from re import compile #, sub
from shutil import copy, rmtree
from sys import exit
from time import time

import click
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.io.votable import parse
from astropy.table import Table

from cosmetic import align_and_median
#from helpers import *
from helpers import write_png, fname_bits,write_conf_file,obj_read_header
from helpers import obj_correct_header
from helpers import flat_read_header,dark_read_header
from helpers import OPT_LIST,HC,UDARK,UFLAT,UOBJ,URED,DARK,FLAT,OBJ,TMP,MASTER
from helpers import DI,FI,RED,STK,AUX,SEX_RES,PSFEX_RES,SCAMP_RES,AR,DATA
from helpers import T120_SEX,T120_PARAM,T120_PSFEX,T120_PARAMPSFEX,DEFAULT_CONV
from helpers import T120_SCAMP,T120_AHEAD
from helpers import SEX_TMP,PSFEX_TMP,SEXAGAIN_OPT_TMP,SCAMP_TMP



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
@click.option('--sex', is_flag=True,
              help='Run the `sex` astromatic command on all reduced images '
              'after reduction.')
@click.option('--psfex', is_flag=True,
              help='Run the `psfex` astromatic command with the '
              '`sex`-determined sources in all reduced images.')
@click.option('--sexagain', is_flag=True,
              help='Run the `sex` astromatic command a second time, using the '
              '`psfex`-determined PSF data.')
#@click.option('--scamp', is_flag=True,
#              help='Run the `scamp` astromatic command on all reduced images '
#              'after reduction.')
#@click.option('--nostack', is_flag=True,
#              help='If True, skip the staking process')
def cli(setup, clear, interpolate, verbose, tmppng, redpng,
        sex, psfex, sexagain):
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

    If --{sex,psfex,sexagain} options are on: run the astromatic suite on the
    reduced images.

    '''
    # Current working directory.
    cwd = basename(getcwd())

    # Welcome.
    click.secho('    Welcome to astro_reduce!\n'
                '    Software is copyright 2018-2021 R. Duque.\n\n',
                fg='cyan', bold=True)
    click.secho('Currently working in directory `{}`.'.format(cwd), fg='green')
    click.secho('Options: {}.'.format(', '.join(filter(eval, OPT_LIST))),
                fg='green')

    # Initialize configuration file name and timer.
    conf_file_name = '{}.json'.format(cwd)
    t0 = time()
    # correspondance between original and temporary file
    u2ar = {}

    # If clear option is on, remove all files and folders and exit.
    if clear:
        click.secho('Clearing directory of astro_reduce files and folders...',
                    fg='green', nl=False)
        if exists(conf_file_name):
            remove(conf_file_name)
        for folder in [OBJ, FLAT, DARK, MASTER, RED, STK, TMP, SEX_RES, PSFEX_RES, SCAMP_RES]:
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
        for folder in [OBJ, DARK, FLAT, MASTER]:
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
        [(exp, glob('{}/{}_{}_*.fit*'.format(DARK, DI, exp)))
            for exp in conf_dic['exposures']])
    flat_files = dict(
        [(filt, glob('{}/{}_{}_*.fit*'.format(FLAT, FI, filt)))
            for filt in conf_dic['filters']])
    
    # work out correspondance between original and processed file
    for file in glob('{}/*.fit*'.format(UOBJ)):
        obj, fil, exp, fn = obj_read_header(file)
        nname = '{}'.format(fn.split('.fit')[0])
#        if verbose:
#            click.secho(' code: '+nname+' ori: '+basename(file))
        u2ar[nname] = basename(file)    
    
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
    if exists(STK):
        rmtree(STK, ignore_errors=True)
    mkdir(STK)
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
        nname = '{}/mdark_{}.fits'.format(MASTER, exp)
        fits.writeto(nname, mdark_data, mdark_header, overwrite=True)
        fits.setval(nname, 'FILTER', value='        ')
        fits.setval(nname, 'IMAGETYP', value='Dark    ')
        fits.setval(nname, 'EXPTIME', value=float(exp / 1000.), comment=HC)
        fits.setval(nname, 'EXPOSURE', value=float(exp / 1000.), comment=HC)
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
        only_mdark = fits.getdata('{}/mdark_{}.fits'.format(MASTER, only_exp))
        a = np.zeros_like(only_mdark)
        b = only_mdark
    else:
        # If not, if you want to fit y = a * x + b,
        # then the LS solution is:
        # a = (<xy> - <x><y>) / (<x ** 2> - <x> ** 2)
        # b = <y> - a * <x>
        mxy = np.mean([float(exp)
                       * fits.getdata('{}/mdark_{}.fits'.format(MASTER, exp))
                       for exp in available_exposures], axis=0)
        mx = np.mean([float(exp) for exp in available_exposures])
        my = np.mean([fits.getdata('{}/mdark_{}.fits'.format(MASTER, exp))
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
        nname = '{}/mdark_{}.fits'.format(MASTER, exp)
        fits.writeto(nname, new_mdark_data, overwrite=True)
        fits.setval(nname, 'FILTER', value='        ')
        fits.setval(nname, 'IMAGETYP', value='Interpolated dark')
        fits.setval(nname, 'EXPTIME', value=float(exp / 1000.), comment=HC)
        fits.setval(nname, 'EXPOSURE', value=float(exp / 1000.), comment=HC)
        fits.setval(nname, 'OBJECT', value='DARK    ')
        
        if verbose:
            click.echo('Done.')

    # STEP 2: Write master transmission files for each filter:
    # mtrans = median(normalized(flat - dark of same exposure)).
    click.secho('Calculating master transmission (Flat) images:', fg='green')
    # Handy function to extract exposure from flat file name.
    fexp = lambda fname: fname.split('.fit')[0].split('_')[-2]
    for filt in flat_files:
        if verbose:
            click.echo('    {:12}   '.format('{}...'.format(filt)), nl=False)

        # Calculate normalized flats.
        normalized_flats = list()
        for fitsfile in flat_files[filt]:
            tmp = fits.getdata(fitsfile) \
                - fits.getdata('{}/mdark_{}.fits'.format(MASTER, fexp(fitsfile)))
            normalized_flats.append(tmp / tmp.mean(axis=0))

        mtrans_data = np.median(normalized_flats, axis=0)
        mflat_header = fits.getheader(flat_files[filt][0])

        # Write fits file and header.
        nname = '{}/mtrans_{}.fits'.format(MASTER, filt)
        fits.writeto(nname, mtrans_data, mflat_header, overwrite=True)
        fits.setval(nname, 'FILTER', value=filt)
        fits.setval(nname, 'IMAGETYP', value='Light Frame')
        fits.setval(nname, 'EXPTIME', value=-1., comment=HC)
        fits.setval(nname, 'EXPOSURE', value=-1., comment=HC)
        fits.setval(nname, 'OBJECT', value='FLAT    ')
        
        if verbose:
            click.echo('Done ({} images).'.format(len(flat_files[filt])))

    # STEP 3: Reduce all the object images with corresponding filter mtrans
    # and exposure mdark.
    noastro = []
    click.secho('Writing auxiliary object images:', fg='green')
    for obj in object_files:
        if verbose:
            click.echo('    {:12}   '.format('{}...'.format(obj)), nl=False)
        for fname in object_files[obj]:
            bfname = basename(fname)
            filt, exp = fname_bits(bfname)

            # Corresponding darks and flats.
            mdark_data = fits.getdata('{}/mdark_{}.fits'.format(MASTER, exp))
            mtrans_data = fits.getdata('{}/mtrans_{}.fits'.format(MASTER, filt))

            # (Raw - Dark) / Trans.
            aux_data = (fits.getdata(fname) - mdark_data) / mtrans_data
            aux_header = fits.getheader(fname)
            try:
                aux_header = obj_correct_header(aux_header)
            except:
                click.secho('W: Not enough astrometric info in file '+fname, fg='magenta')
                noastro.append(fname)

            # Write fits file and header.
            nname = '{}/{}_{}.fits'.format(TMP, bfname.split('.fit')[0], AUX)
            fits.writeto(nname, aux_data, aux_header, overwrite=True)
            fits.setval(nname, 'FILTER', value=filt)
            fits.setval(nname, 'IMAGETYP', value='Light Frame')
            fits.setval(nname, 'EXPTIME', value=float(exp) / 1000, comment=HC)
            fits.setval(nname, 'EXPOSURE', value=float(exp) / 1000, comment=HC)
            fits.setval(nname, 'OBJECT', value=obj)
            
        if verbose:
            click.echo('Done.')

    # STEP 4: For all objects realign and median the aux images.
    # You are left with one image per object per filter per exposure.
    # TODO: remove the following line and replace by program argument
    nostack = True
    stack = not nostack
    if stack:
        click.secho('Realigning/Stacking object images:', fg='green')
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
                stacked_data = align_and_median(aux_files)
                stacked_header = fits.getheader('{}/{}'
                                                ''.format(OBJ,
                                                          names_per_tag[tag][0]))
                try:
                    stacked_header = obj_correct_header(stacked_header)
                except:
                    click.secho('W: Not enough astrometric info in file '+aux_files[0], fg='magenta')            
                
                # Write fits file and header.
                nname = '{}/{}_{}_{}.fits'.format(STK, obj, f, e)
                fits.writeto(nname, stacked_data, stacked_header, overwrite=True)
                fits.setval(nname, 'FILTER', value=f)
                fits.setval(nname, 'IMAGETYP', value='Light Frame')
                fits.setval(nname, 'EXPTIME', value=float(e) / 1000., comment=HC)
                fits.setval(nname, 'EXPOSURE', value=float(e) / 1000., comment=HC)
                fits.setval(nname, 'OBJECT', value=obj)
                if verbose:
                    click.echo('       Done ({} images).'.format(len(aux_files)))

    # STEP 5: If options redpng or tmppng are on, write
    # PNG versions of all the tmp and reduced/stacked images.
    if redpng:
        click.secho('Writing PNG versions of reduced/stacked images... ',
                    fg='green', nl=False)
        for ffile in glob('{}/*.fits'.format(RED)):
            write_png(ffile, plt)
        for ffile in glob('{}/*.fits'.format(STK)):
            write_png(ffile, plt)
        click.secho('Done.', fg='green')

    if tmppng:
        click.secho('Writing PNG versions of intermediate images... ',
                    fg='green', nl=False)
        for ffile in glob('{}/*.fits'.format(TMP)):
            write_png(ffile, plt)
        click.secho('Done.', fg='green')

    # STEP 6: If options sex, psfex or sexagain are on,
    # run the astromatic suite.
    # Sexagain implies psfex and psfex implies sex.
    psfex = psfex or sexagain
    sex = sex or psfex
    scamp = True # TODO: REMOVE THIS LINE AND PUT IT AS ARGUMENT
    if sex or psfex or sexagain or scamp:
        # Setup for the astrometry: initialize empty result folders
        for folder in [SEX_RES, PSFEX_RES, SCAMP_RES]:
            if exists(folder):
                rmtree(folder, ignore_errors=True)
            mkdir(folder)

        # Setup for the astrometry: find configuration files in the file system
        t120_sex = resource_filename(AR, '{}/{}'.format(DATA, T120_SEX))
        t120_param = resource_filename(AR, '{}/{}'.format(DATA, T120_PARAM))
        t120_psfex = resource_filename(AR, '{}/{}'.format(DATA, T120_PSFEX))
        t120_parampsfex = resource_filename(AR, '{}/{}'.format(DATA,
                                                               T120_PARAMPSFEX))
        default_conv = resource_filename(AR, '{}/{}'.format(DATA, DEFAULT_CONV))
        t120_scamp = resource_filename(AR, '{}/{}'.format(DATA, T120_SCAMP))
        t120_ahead = resource_filename(AR, '{}/{}'.format(DATA, T120_AHEAD))
        
        # set the directory of images to process astrometry
        SXD = TMP
    
    
    # Run sextractor.
    if sex:
        for ffile in glob('{}/*.fits'.format(SXD)):
            if ffile in noastro:
                click.secho('W: No enough astrometry info in file '+ffile,fg='magenta')
                continue
            stem = basename(ffile.split('.fit')[0])
            sex_cmd = SEX_TMP.format(ffile, t120_sex, t120_param, default_conv,
                                  SEX_RES, stem + '.ldac',
                                  SEX_RES, stem + '-bckg.fits',
                                  SEX_RES, stem + '-obj.fits',
                                  SEX_RES, stem + '.xml')
#            if verbose:
#                click.secho(' submitting SExtractor cmd: '+sex_cmd, nl=True)
            system(sex_cmd)

    # Run PSFEx with sextractor-determined sources.
    if psfex:
        for ffile in glob('{}/*.fits'.format(SXD)):
            stem = basename(ffile.split('.fit')[0])
            psfex_cmd = PSFEX_TMP.format(t120_psfex, stem + '.xml')
#            if verbose:
#                click.secho(' submitting PSFex cmd: '+psfex_cmd, nl=True)
            system(psfex_cmd)

    # Run sextractor with PSFEx-determined PSF.
    if sexagain:
        for ffile in glob('{}/*.fits'.format(SXD)):
            if ffile in noastro:
                click.secho('W: No enough astrometry info in file '+ffile,fg='magenta')
                continue
            stem = basename(ffile.split('.fit')[0])
            sexagain_cmd = SEX_TMP.format(ffile, t120_sex, t120_parampsfex,
                                  default_conv,
                                  SEX_RES, stem + '.ldac',
                                  SEX_RES, stem + '-bckg.fits',
                                  SEX_RES, stem + '-obj.fits',
                                  SEX_RES, stem + '.xml')\
                         + SEXAGAIN_OPT_TMP.format(stem + '.psf')
#            if verbose:
#                click.secho(' submitting again SExtractor cmd: '+sexagain_cmd, nl=True)
            system(sexagain_cmd)
   
    # Run SCAMP
    if scamp:
        listldac = glob('{}/*.ldac'.format(SEX_RES))
        scamp_cmd = SCAMP_TMP.format(' '.join(listldac), t120_scamp, t120_ahead)
#        if verbose:
#            click.secho(' submitting SCAMP cmd: '+scamp_cmd, nl=True)
        system(scamp_cmd)
        
        #  read scamp result
        scampxml = 'scamp.xml'
        votable = parse(scampxml)
        table = votable.get_first_table()
        catcont = Table([(table.array['Catalog_Name']).data,
                         (table.array['XY_Contrast']).data],names=['name','contrast'])

        min_contrast = 10.0
        # put astrometry info into image fits header
        for fhed in glob('{}/*.head'.format(SEX_RES)):
            # get/set original/reduced image file name
            key = splitext(basename(fhed))[0].replace('_aux','')
            fori = '{}/{}'.format(UOBJ, u2ar[key])
            fred = '{}/{}'.format(URED,basename(fori))
            ffts = '{}/{}'.format(SXD,basename(fhed.replace('.head','.fits')))
            # check if SCAMP contrast is enough
            mask = catcont['name']==basename(fhed).replace('.head','.ldac')
            contrast = catcont[mask]['contrast'][0]
            click.secho('{} = contrast: {}'.format(fori,contrast))
            if (contrast < min_contrast):
                click.secho('W: Contrast too low ({}) for {}'.format(contrast,fori), fg='magenta')
                remove(fhed)
                remove(fhed.replace('.head','.ldac'))
                continue
            
            # update header
            red_data = fits.getdata(ffts)
            red_header = fits.Header.fromtextfile(fhed)
            fts_header = fits.getheader(ffts)
            # write first version of fits file
            fits.writeto(fred, red_data, fts_header, overwrite=True)
            # now update header
            for hdr_key in red_header:
                if hdr_key=='HISTORY' or hdr_key=='COMMENT' or 'FLXSCALE':
                    continue
                fits.setval(fred,hdr_key , value=red_header[hdr_key])
            
            if verbose:
                click.secho('Reduced data saved in '+fred)
    
    
    
    # Report execution time.
    t1 = time()
    click.secho('\nAll done. ({})'.format(timedelta(seconds=int(t1 - t0))),
                fg='green')
    # ALL DONE.
