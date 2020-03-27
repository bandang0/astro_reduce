# astro_reduce -- A Simple CCD Images Reducer for the Paris Observatory.

## What does astro_reduce do?

`astro_reduce` reduces raw astronomical CCD images read in FITS format with dark and flat field images. Optionally, it can interpolate missing dark fields from the ones available, and can also convert the intermediate and reduced images to PNG format for easy inspection.

To operate, `astro_reduce` is to be launched from a directory containing three folders:
- A `DARK` folder, containing all the dark field images,
- A `FLAT` folder, containing all flat field images,
- A `ORIGINAL` folder, containing all the images of objects.

All files must be in FITS format, with correct filter, exposure time and object header keywords wherever relevant.

When launched, `astro_reduce` will copy all data files to working folders (thereby backing up the original data), and apply the standard astronomical CCD image reduction process. On the way, intermediate images such as master dark fields or non-realigned object fields are safe-kept in the `tmp/` folder.

After running, the final reduced images can be found in the `reduced/` folder. The reduction process and the nomenclature for the naming of the files in the `tmp/` and `reduced/` folders are given in the _Reduction method_ paragraph below.

__Note:__ It is highly recommended to inspect the intermediate images produced by `astro_reduce` before considering the final reduced images. Also, we highly recommend you inspect your raw images and eliminate any bad ones before launching the reduction.

## Installing
The `astro_reduce` program can be installed from the PyPI under the name `astro-reduce`. This package provides the command with all options.

Alternatively, use the `setup.py` file in the project directory with `python3 setup.py install --user` to install the program locally. You should then include `~/.local/bin` (or `~/Library/Python/3.7/bin` for Mac users) to your path in order to invoke `astro_reduce` from the console.

### Dependencies
`astro_reduce` is written in Python 3.5.

`astro_reduce` depends on the `click`, `astropy`, `numpy`, `scipy` and `matplotlib` Python packages, which are all available through the PyPI.

## Usage
When used for the first time in a directory, or when the data has been modified since the last use of `astro_reduce`, the program should be invoked with the `--setup` option:

`astro_reduce --setup [OPTIONS]`

This will setup the reduction by copying all the raw data to `astro_reduce`'s working folders whilst changing the file names to standardized formats. Also, during this phase, a configuration file in the JSON format will be written. It summarizes the objects, filters and exposures times present in the original data and is used by `astro_reduce` in the reduction process. This file should be left in the directory for further reference or for later rerunning of the program.

__Note:__ This first setup step is effectively a back-up of your data, as the data left in the `DARK`, `FLAT` and `ORIGINAL` folders are no longer touched during the subsequent reduction process.

If the setup has already been done in the directory by prior use of the `--setup` option, `astro_reduce` should be used without this option:

`astro_reduce [OPTIONS]`

`astro_reduce` has been tested on Linux and Mac platforms, and has yet to be tested on Windows.

### Options
- `-s, --setup        Sets up the directory for reduction. Use this option only the first time astro_reduce is run in the directory.`
- `-i, --interpolate  Interpolate existing dark fields if some are missing.`
- `-v, --verbose      Enables verbose mode (recommended).`
- `-t, --tmppng       Write PNG format of intermediary images after reduction.`
- `-r, --redpng       Write PNG format of reduced images after reduction.`
- `--help             Show this message and exit.`


## Reduction method
### Master dark images
For a given exposure time, the _master dark_ is calculated as the pixel-wise __median__ of all the dark fields of that exposure. This allows to eliminate cosmic ray traces.

#### Dark field interpolation
In the case where some dark fields are missing, `astro_reduce` can interpolate from the available master darks to obtain master darks for all the exposures necessary to reduce the object or flat images. This is done by specifying the `--interpolate` option.

The interpolation is _least-square linear_, i.e., two images A and B are determined from the available master dark images such as to __minimize the square error__ on the linear interpolation (master dark) = (exposure time) x A + B.

Using these A and B, the missing master darks are calculated according to this equation.

>The FITS files for all the master dark images (deduced from dark fields or interpolated) can be found after reduction in the `tmp` directory under the names `mdark_[exposure].fits`.

### Master transmission images
The _master transmission_ image for a given filter is an image which encompasses the relative transmission of each pixel in the optical setup (telescope optics through filter to CCD matrix). It is calculated for every filter as the __median__ over all flat field images, after __subtraction of corresponding master dark images__ and __normalization__.

>The FITS files for all the master transmission images can be found after reduction in the `tmp/` directory under the names `mtrans_[filter].fits`.

### Individual image reduction
An object image of given exposure and filter is reduced by __subtracting__ the corresponding exposure master dark image, and __dividing__ by the corresponding filter master transmission image.

>The FITS files for object _obj_, with filter _filt_ and exposure time _exp_ are found in the `tmp/` folder after reduction, under the name `[obj]_[filt]_[exp]_*_aux.fits` (for _auxiliary_).

### Realignment and reduction
Finally, for each series of same exposure and filter for each object, the auxiliary files are realigned through optimization of their __mutual cross-correlations__, and then their pixel-wise __median__ image is calculated. Using the median rather than the mean allows to efficiently remove hot pixels. This will be all the more effective as dithering has been used in acquiring the images of a series.

During the realignment, images are rolled to superpose themselves. Therefore, if the images were too misaligned to begin with, objects that roll beyond the edge can end up in odd places, potentially producing ghost images of objects. During `astro_reduce`'s run, a warning is issued to the user if any image is rolled by more than a quarter-field during the realignment.

>These realigned images are the final reduced images and can be found in the `reduced/` folder after reduction, with the same names as in the preceding step.

## Some things to know
- As mentioned earlier, the original files in the `DARK`, `FLAT` and `ORIGINAL` folders are not used during the reduction process. Only their copies in the `ar_*` working folders are used.
- In copying the original files during the setup process, `astro_reduce` in effect backs up your data.
- `astro_reduce` will deduce the object captured in a given field from the `OBJECT` keyword entry in the file's header. In the reduction, all files of a same object, filter and exposure time are reduced and realigned together. Therefore, in order to keep different series of images of a same object separate, it is recommended to append a tag to the object name in acquiring the images of the series. For example NGC4993-1, NGC4993-2, etc.
- During the reduction, the headers of synthetic images are filled in with the correct `OBJECT`, `FILTER`, `EXPOSURE`, `EXPTIME` and `IMAGETYP` keyword values. The other keyword values (such as the time of acquisition) are inherited from the parent images (raw darks and flats for master dark and master transmission images; raw object images for final reduced object fields).
- `fit` and other derivatives can be used as extensions for the FITS files.
- Remaining hot pixels in the reduced images can be the result of insufficient dithering.
