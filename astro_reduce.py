'''astro_reduce -- A Simple CCD Images Reducer for the Paris Observatory.'''

from sys import exit
from os.path import basename, exists
from os import mkdir
from glob import glob
from json import loads
from collections import defaultdict
import re

import click
from astropy.io import fits
from scipy.signal import fftconvolve
import numpy as np

# Paths and extensions.
TMP = 'tmp'
OBJ = 'objects'
DARK = 'darks'
FLAT = 'flats'
AUX = 'aux'
RED = 'reduced'

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

# Return the filter and exposure (strings) from an object file name.
# 'NGC1000_1_V_1000_0.fits' gives ('1', 'V', '1000')
def fname_bits(fname):
    '''Return the series, filter and exposure from an object file name.'''
    pieces = fname.split('.fit')[0].split('_')
    return (pieces[1], pieces[2], pieces[3])

# Write png from fits version of image, in same directory.
def write_png(fname, plt):
    '''Write PNG version of image from fits file.'''
    plt.figure(1)
    plt.imshow(fits.getdata(fname), aspect='auto', origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig(f'{fname.split(".fit")[0]}.png', bbox_inches='tight')
    plt.close(1)


@click.command()
@click.argument('conf_file', type=click.File('r'))
@click.option('--verbose', '-v', is_flag=True,
              help='Enables verbose mode (recommended).')
@click.option('--tmppng', '-t', is_flag=True,
              help='Write PNG format of intermediary images after reduction.')
@click.option('--redpng', '-r', is_flag=True,
              help='Write PNG format of reduced images after reduction.')
@click.option(
    '--interpolate',
    '-i',
    is_flag=True,
    help='Interpolate existing dark field images if some are missing.')
@click.option(
    '--cross',
    '-c',
    is_flag=True,
    help='Realign across series of a same object, filter and exposure.')
def cli(conf_file, verbose, tmppng, redpng, interpolate, cross):
    '''Reduce CCD images from objects with flat and dark field images.'''
    # Parse configuration file to obtain configuration dictionary.
    if verbose:
        click.echo(f'Parsing configuration file {conf_file.name}.')
    conf_dic = loads(conf_file.read())
    dn, fn = conf_dic['dark_name'], conf_dic['flat_name']

    # Obtain list of all object, dark, flat field file names.
    object_files = dict(
        [(obj, glob(f'{OBJ}/{obj}_*.fit*'))
            for obj in conf_dic['objects']])
    dark_files = dict(
        [(exp, glob(f'{DARK}/{dn}_{exp}_*.fit*'))
            for exp in conf_dic['exposures']])
    flat_files = dict(
        [(filt, glob(f'{FLAT}/{fn}_{filt}_*.fit*'))
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
    tstring = "****** {:25} ******"
    sstring = "    {:8}: {}"
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
            click.echo("Done.")

    # STEP 3: Reduce all the object images with corresponding filter mflat
    # and exposure mdark. Do this regardless of series and number in series.
    click.echo('Writing auxiliary object images.')
    for obj in object_files:
        if verbose:
            click.echo(f'    {obj}... ', nl=False)
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
        if verbose:
            click.echo('Done.')

    # STEP 4: Within series, realign the aux images and write the median images
    # of those. You are left with one image per object per filter per exposure
    # per series.
    click.echo('Realigning object images.')
    for obj in object_files:
        if verbose:
            click.echo(f'    {obj}:')
        # Group all the object files by *tag*, i.e. by series, filter, exposure.
        name_tag_hash = [(basename(fname), f'{fname_bits(basename(fname))}')
                         for fname in object_files[obj]]
        names_per_tag = defaultdict(list)
        for name, tag in name_tag_hash:
            names_per_tag[tag].append(name)

        # Now you align images which have the same tag.
        for tag in names_per_tag:
            # Rebuild series, filter and exposure from tag (they are those of,
            # e.g., the first name in the list.)
            s, f, e = fname_bits(names_per_tag[tag][0])
            if verbose:
                click.echo(f'      {s}/{f}/{e}... ', nl=False)
            # Calculate aligned and medianed image
            # from all images with same tag.
            aux_files = glob(f'{TMP}/{obj}_{s}_{f}_{e}_*_{AUX}.fits')
            reduced_data = align_and_median(aux_files)
            reduced_header = fits.getheader(f'{OBJ}/{names_per_tag[tag][0]}')

            fits.writeto(f'{RED}/{obj}_{s}_{f}_{e}.fits', reduced_data,
                         reduced_header, overwrite=True)
            if verbose:
                click.echo('Done.')

    # STEP 4.5: If the cross-series option is on, realign and median
    # reduced images across the series (but within same filter and exposure).
    if cross:
        click.echo('Realigning object images across series.')
        for obj in object_files:
            if verbose:
                click.echo(f'    {obj}... ', nl=False)
            # Group all reduced files of object across series, i.e. by filter
            # and exposure ("fe").
            name_fe_hash = [(basename(fname),
                             f'{fname_bits(basename(fname))[1:]}')
                            for fname in glob(f'{RED}/{obj}_*_*_*.fits')]
            names_per_fe = defaultdict(list)
            for name, fe_hash in name_fe_hash:
                names_per_fe[fe_hash].append(name)

            # Now align images with same filter and exposure.
            for fe in names_per_fe:
                # Get corresponding filter and exposure.
                example_file = names_per_fe[fe][0]
                _, f, e = fname_bits(example_file)

                # Calculate realigned image for all the images of object
                # with filter "f" and exposure "e".
                red_files = glob(f'{RED}/{obj}_*_{f}_{e}.fits')
                if len(red_files) < 2:
                    # Only one series, no realignment to do.
                    continue
                aligned_data = align_and_median(red_files)
                aligned_header = fits.getheader(f'{RED}/{example_file}')
                fits.writeto(f'{RED}/{obj}_{f}_{e}.fits', aligned_data,
                             aligned_header, overwrite=True)
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
