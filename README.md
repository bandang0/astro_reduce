# astro_reduce: A simple CCD-image reducer for the _Observatoire de Paris-Meudon_

[![image](http://img.shields.io/pypi/v/astro-reduce.svg)](https://pypi.python.org/pypi/astro-reduce/)

## What does astro_reduce do?
`astro_reduce` carries out the cosmetic and astrometric reduction of astronomical CCD images read in FITS format. The cosmetic reduction is done in a standard fashion with dark and flat field images. The astrometric reduction is done via calling the [Astromatic Software Suite](https://www.astromatic.net/) on the images obtained from the cosmetic reduction.

Optionally, it can interpolate dark fields for exposure times missing in the available data and can also convert the intermediate and reduced images to PNG format for easy inspection. The astrometric reduction can be carried out anywhere from simply extracting sources to catalogs up to updating the image headers with astrometric output such as the linear projection matrix coefficients. For this, the relevant programs--either `sex`, `psfex` or `scamp` accordingly with the specified options--must be installed on your system. The configuration files for these commands are handled internally by `astro_reduce`, you do not need to provide them.

To operate, `astro_reduce` is to be launched from a directory containing three folders:
- A `DARK` folder, containing all the dark field images,
- A `FLAT` folder, containing all flat field images,
- A `ORIGINAL` folder, containing all the images of objects.

All files must be in FITS format, with correct filter, exposure time and object header keywords wherever relevant.

The flowchart below illustrates how `astro_reduce` runs. The boxes describe the content of folders found after the run, with directory names and name formats for files therein. The lines describe the flow of data during the run: solid lines indicate the standard run of the program, while dashed lines can be realized optionally.

When launched for the first time in a directory, use the `--setup` option. This will copy all data files to working folders (thereby backing up the original data). Then, run `astro_reduce` again with relevant options to run the cosmetic and astrometric reductions; see details below on usage.

During cosmetic reduction, intermediate images such as master dark fields or non-stacked object images are safe-kept in the `ar_masters/` and `ar_tmp/` folders. The final stacked images can be found in the `stacked/` folder. The cosmetic reduction process and the nomenclature for the naming of the files in the `ar_masters/`, `ar_tmp/` and `stacked/` folders are given in the _Cosmetic reduction method_ paragraph below. If options to perform the astrometric reduction were used, the results (such as source catalogs) can be found after running in the `SEXRES`, `PSFRES` and `SCAMPRES` folders. The images with updated astrometric header data can be found in the `reduced/` folder, see _Astrometric reduction method_ below.

__Note:__ It is highly recommended to inspect the intermediate images produced by `astro_reduce` before considering the final reduced images. Also, we highly recommend you inspect your raw images and eliminate any bad ones before launching the software.

![flowchart](info/flowchart.jpeg)

## Installing
The `astro_reduce` program can be installed from the PyPI under the name `astro-reduce`. This package provides the command with all options.

Alternatively, use the `setup.py` file in the project directory with `python3 setup.py install --user` to install the program locally. You should then include `~/.local/bin` (or `~/Library/Python/3.x/bin` for Mac users) to your path in order to invoke `astro_reduce` from the console.

### Dependencies
`astro_reduce` is written in Python 3.

`astro_reduce` depends on the `click`, `astropy`, `numpy`, `scipy` and `matplotlib` Python packages, which are all available through the PyPI.

## Usage
When used for the first time in a directory, or when the data has been modified since the last use of `astro_reduce`, the program should be invoked with the `--setup` option:

`astro_reduce --setup [--verbose]`

This will setup the reduction by copying all the raw data to `astro_reduce`'s working folders whilst changing the file names to standardized formats, as seen in the flowchart. Also, during this phase, a configuration file in the JSON format will be written. It summarizes the objects, filters and exposures times present in the original data and is used by `astro_reduce` in the reduction process. This file should be left in the directory for further reference or for later rerunning of the program.

__Note:__ This first setup step is effectively a back-up of your data, as the data left in the `DARK`, `FLAT` and `ORIGINAL` folders are no longer touched during the subsequent reduction process.

If the original data has changed or your are not sure when the `--setup` option was used last, you can remove all `astro_reduce` working data with the `--clear` option.

If the setup has already been done in the directory by prior use of the `--setup` option, `astro_reduce` should be used without this option:

`astro_reduce [OPTIONS]`

`astro_reduce` has been tested on Linux and Mac platforms, and has yet to be tested on Windows.

### Options

- `--version          Show the version and exit.`
- `-s, --setup        Set up the directory for reduction. Use this option the first time astro_reduce is run in the directory or after the '--clear' option was used.`
- `-c, --clear        Remove all astro_reduce-related files and folders in current directory and exit.`
- `-i, --interpolate  Interpolate existing dark fields if some are missing.`
- `-v, --verbose      Enables verbose mode (recommended).`
- `-t, --tmppng       Write PNG format of auxiliary and master images after cosmetic reduction.`
- `-r, --stkpng       Write PNG format of stacked images after cosmetic reduction.`
- `--sex              Run the 'sex' Astromatic command on all auxiliary images after the cosmetic reduction.`
- `--psfex            Run the 'psfex' Astromatic command with the 'sex'-determined sources in all auxiliary images.`
- `--sexagain         Run the 'sex' Astromatic command a second time, using the 'psfex'-determined PSF data.`
- `--scamp            Run the 'scamp' Astromatic command on all auxiliary images after cosmetic reduction.`
- `--nomaster         If set, do not calculate the master darks and flats (and assume they are already there!).`
- `--nostack          If set, skip the stacking process.`
- `--help             Show this message and exit.`

## Cosmetic reduction method
### Master dark images
For a given exposure time, the _master dark_ is calculated as the pixel-wise __median__ of all the dark fields of that exposure. This allows to eliminate cosmic ray traces.

#### Dark field interpolation
In the case where dark fields for some exposure times are missing, `astro_reduce` can interpolate from the available master darks to obtain master darks for all the exposures necessary to reduce the object or flat images. This is done by specifying the `--interpolate` option.

The interpolation is _least-square linear_, i.e., two images A and B are determined from the available master dark images such as to __minimize the square error__ on the linear interpolation (master dark) = (exposure time) x A + B.

Using these A and B, the missing master darks are calculated according to this equation.

>The FITS files for all the master dark images (deduced from dark fields or interpolated) can be found after run in the `ar_masters/` directory under the names `mdark_[exposure].fits`.

### Master transmission images
The _master transmission_ image for a given filter is an image which encompasses the relative transmission of each pixel in the optical setup (telescope optics through filters to CCD matrix). It is calculated for every filter as the __median__ over all flat field images, after __subtraction of corresponding master dark images__ and __normalization__.

>The FITS files for all the master transmission images can be found after reduction in the `ar_masters/` directory under the names `mtrans_[filter].fits`.

If you have already run `astro_reduce` or have your own master dark and flat images, you can skip these two last steps with the `--nomaster` option. Beware that, in this case, you must place your custom master files in the `ar_masters/` folder before running.

### Individual image reduction
An object image of given exposure and filter is cosmetically reduced by __subtracting__ the corresponding exposure master dark image, and __dividing__ by the corresponding filter master transmission image. We refer to these images as the _auxiliary_ images, as in the flowchart.

>The FITS files for object _obj_, with filter _filt_ and exposure time _exp_ are found in the `ar_tmp/` folder after reduction, under the name `[obj]_[filt]_[exp]_*_aux.fits` (for _auxiliary_), which contains an `astro_reduce`-internal reference number.

### Realignment and stacking
Finally, for each series of same exposure and filter for each object, the auxiliary images stacked: they are realigned through optimization of their __mutual cross-correlations__, and then their pixel-wise __median__ image is calculated. Using the median rather than the mean allows to efficiently remove hot pixels. This will be all the more effective as dithering has been used in acquiring the images of a series.

During the realignment, images are rolled to superimpose themselves. Therefore, if the images were too misaligned to begin with, objects that roll beyond the edge can end up in odd places, potentially producing ghost images of objects. During `astro_reduce`'s run, a warning is issued to the user if any image is rolled by more than 15% of the field size during the realignment. In this case, it may be useful to remove the offending image from the dataset. If there are many such images and they are aligned amongst themselves, it is advised to change their header `OBJECT`s in order to process them together as a separate batch.

>These stacked images are the final cosmetically reduced images and can be found in the `stacked/` folder after reduction, with the same names as in the preceding step, save the `_aux` suffix and the internal image reference number.

This stacking step can be altogether skipped using the `--nostack`, as appropriate for, e.g., fast moving objects.

## Some things to know
- As mentioned earlier, the original files in the `DARK`, `FLAT` and `ORIGINAL` folders are not used during the reduction process. Only their copies in the `ar_*` working folders are used.
- In copying the original files during the setup process, `astro_reduce` in effect backs up your data.
- `astro_reduce` will deduce the object captured in a given field from the `OBJECT` keyword entry in the file's header. In the reduction, all files of a same object, filter and exposure time are cosmetically reduced and stacked together. Therefore, in order to keep different series of images of a same object separate, it is recommended to append a tag to the object name in acquiring the images of the series. For example NGC4993-1, NGC4993-2, etc.
- During the cosmetic reduction, the headers of synthetic images are filled in with the correct `OBJECT`, `FILTER`, `EXPOSURE`, `EXPTIME` and `IMAGETYP` keyword values. The other keyword values (such as the time of acquisition) are inherited from the parent images (raw darks and flats for master dark and master transmission images; raw object images for final reduced object fields).
- The cosmetically reduced images found in `ar_tmp/` and then `stacked/` are more voluminous (typically 4 times) then the original images, because they are encoded as floats following the dark-flat manipulations.
- `fit` and other derivatives can be used as extensions for the FITS files.
- Remaining hot pixels in the reduced images can be the result of insufficient dithering.

## Astrometric reduction method

The astrometric reduction is launched with the `--sex`, `--psfex`, `--sexagain` and `--scamp` options, depending on how far you would like to take the reduction:
- Using `--sex` runs the `SExtractor` program on all auxiliary files and places the results (source catalogs, etc.) in the `SEXRES/` folder, as in the flowchart.
- Using `--psfex` (which implies `--sex`) runs the `PSFEx` program on the auxiliary images with `SExtractor`'s results to produce PSF-related data in the `PSFRES/` folder.
- Using `--sexagain` (which implies `--psfex`) runs `SExtractor` a second time on the auxiliary images, but using the PSF data from the `PSFEx` run. This updates the results in `SEXRES/`.
- Finally, using `--scamp` (which implies `--sex`) runs `SCAMP` to perform the astrometric reduction on the auxiliary images, using the source catalogs from `SEXRES/` and data on remote servers.

>If the `--scamp` option was used, the final astrometrically reduced files with updated headers can be found in `reduced/` after run.

__Note:__ `astro_reduce` ships with all configuration files for the astrometric reduction, you need not provide them. No astrometry is performed on the stacked images, only on the auxiliary images.

## Enhancing `astro_reduce`

### How can I add an option to `astro_reduce`?

The project uses the [Click](https://palletsprojects.com/p/click/) Python package to handle the CLI options. To add an option, follow these steps:
1. Add a `@click.option` decorator to the `cli` function defined in `astro_reduce.py`. This allows to define the option's command-line flag, its short version and its help sentence. You can take inspiration from the existing options;
2. The previous step will define a argument to the `cli` function, with exactly the spelling of the CLI's flag; You must add this argument to the signature of the `cli` function;
3. So that the listing of current options in the welcome message continues to work properly, you must add the option to the `OPT_LIST` in `helpers.py`;
4. If the option extends the astrometric reduction capabilities of `astro_reduce`, please also add the option to the `ASTROMATIC_LIST` in `helpers.py`.

For more information on defining options, see the [Click documentation](https://click.palletsprojects.com/en/8.0.x/#documentation).

### Making releases to the PyPi

We use the PyPi to distribute `astro_reduce`. To make a release, follow these steps:
1. Define the version number of your release in the `x.y.z` form in `setup.py`. This version must not already exist in the PyPi, or else you will run into trouble down the line. We use [semantic versioning](https://semver.org/):
- Change `x` when you make backwards-incompatible changes, such as changing an existing option or the output format;
- Change `y` when you make backwards-compatible changes, such as introducing new functionality through a new option or new documentation;
- Change `z` when you fix minor bugs without noteworthy functionality improvement.
2. Make sure you have updated the documentation (i.e., this `README`) such that it fits the new version;
3. Make a commit to the git repository where you indicate the number of the new version;
4. Run `python3 setup.py bdist_wheel` and then `twine upload dist/*`.

This last step will prompt you for your PyPi credentials. Naturally, it will only work if you have maintainer or owner status on the project's PyPi repository.

### Adding or changing a data file

We also use the PyPi to distribute data files for `astro_reduce`. To include a new file to the distribution, do the following:
1. Add the file to the `data/` folder in the project repository;
2. To make sure it is distributed with the rest of the project, list it in the `data_files` option in `setup.py`;
3. To retrieve the path to the file in the `astro_reduce.py` script in order to use the file, we use
the `resource_filename` function of the `pkg_resources` package: If you added `my_new_file.txt` to `data/` in the repo, the path to file `my_new_file.txt` in the user's system will be the result of the call to `pkg_resources.resource_filename` with arguments `astro_reduce` and `data/my_new_file.txt`. For example, see how we retrieve the paths to the Astromatic configuration files around line 570 in `astro_reduce.py`.

### Potential improvements

- Allow the user to specify their own configuration file for the Astromatic programs (thus overriding the ones distributed through `pip`).
- Allow the user to run the Astromatic suite on the stacked images as well, not only on the auxiliary images.
