# astro_reduce -- A Simple CCD Images Reducer for the Paris Observatory.

## What does astro_reduce do?

`astro_reduce` reduces astronomical CCD images in FITS format by treating raw images with dark and flat field images. It keeps the intermediate images produced along the reduction process.

Optionally, it can interpolate missing dark field images from the available ones, and also convert the intermediate and reduced images to PNG format for easy inspection.

To operate, `astro_reduce` needs to be launched from a directory organized in a standardized fashion (see _Nomenclature and configuration file_ paragraph below).

__Note:__ It is highly recommended to back up your raw data before using `astro_reduce`, and to inspect the intermediate images produced by `astro_reduce` before considering the final reduced images.

## Installing
Use the `setup.py` file with `python3 setup.py install --user` to install the program locally. You should then include `~/.local/bin` (or `~/Library/Python/3.7/bin` for Mac users) to your path in order to invoke `astro_reduce` from the console.

### Dependencies
`astro_reduce` depends on the `click`, `astropy`, `numpy`, `scipy` and `matplotlib` (only for the `--tmppng` and `--redpng` options) Python packages, which are all available through the PyPI.

## Usage
In the directory containing the `.json` configuration file (see the _Configuration file_ paragraph below), launch the reduction with:

`astro_reduce [OPTIONS] CONF_FILE`

You will then find the intermediate reduction files in the `tmp/` directory, and your reduced images in the `reduced/` directory.

`astro_reduce` has been tested on Linux and Mac platforms, and has yet to be tested on Windows.

### Options
- `  -v, --verbose : Enables verbose mode (recommended).`
- `  -t, --tmppng : Write PNG format of intermediary images after reduction (useful for inspection after reduction).`
- `  -r, --redpng : Write PNG format of reduced images after reduction (idem --tmppng).`
- `  -i, --interpolate : Interpolate existing dark field images if some are missing.`
- `  -c, --cross : Realign across series of same object, filter and exposure`
- `  --help : Show help message.`

## Directory structure and configuration file
### Initial folder structure
`astro_reduce` will operate in a folder which must contain:

- A `darks/` folder, containing the dark field FITS files, with names according to `[DARK_NAME]_[EXPOSURE_TIME_IN_MS]_[NUMBER].fits`. Here `[DARK_NAME]` is the dark field file identifier (read in the configuration file, see next paragraph), and `[NUMBER]` is the number of the image in a series of dark fields of same exposure.
- A `flats/` folder, containing the flat field FITS files, with names according to `[FLAT_NAME]_[FILTER]_[EXPOSURE_TIME_IN_MS]_[NUMBER].fits`. Here `[FLAT_NAME]` is the flat field file identifier (also read in the configuration file), `[FILTER]` is one of the filters specified in the configuration file, and `[NUMBER]` is the number of the image in a series of flats of same exposure and filter.
- An `objects/` folder, containing all the object images, with names according to `[OBJECT_NAME]_[SERIES_NUMBER]_[FILTER]_[EXPOSURE_TIME_IN_MS]_[NUMBER].fits`. Here, `[OBJECT_NAME]` is one of those read in the configuration file, and `[SERIES_NUMBER]` is the number of the series (images taken with same object, exposure and filter).

__Note:__
- All of the bits in brackets are mandatory for each file name,
- Seeing this nomenclature, you should refrain from using `_` in object names,
- `fit` and other derivatives can be used as extensions for the FITS files.

### Configuration file
For a series of images to reduce, the configuration file must specify:
- all the exposure times (in ms) of all the images taken (including those of dark and flat fields),
- filters,
- object names,
- flat and dark image identifiers.

All of this information is contained in a `.json` file and organized as follows:

```yaml
{
    "objects": ["BRACKETED", "LIST", "OF", "OBJECT", "NAMES"],
    "filters": ["BRACKETED", "LIST", "OF", "FILTERS"],
    "exposures": [BRACKETED, LIST, OF, EXPOSURE, TIMES, IN, MS],
    "flat_name": "THE STRING USED IN FLATS FITS FILE NAMES",
    "dark_name": "THE STRING USED IN DARKS FITS FILE NAMES"
}
```

__Note:__
- Of course, all objects or flats need not be imaged with every filter or every exposure time,
- It suffices to specify each exposure time and filter used once,
- Also, numbers of images in series need not start with 0 or 1 (see the example directory below),
- For the time being, the configuration file is not checked for consistency with the data before reduction. Thus reducing images with unspecified filters or filters with no corresponding flat field images may lead to crashes.

## Example directory structure and configuration file
Say you observed the following objects on August 17th 2017:
- the NGC4993 galaxy (at two different moments of the night, hence two series),
- the B612 asteroid,
- Regulus to serve as a reference star for photometry of NGC4993.

You will thus have two series of NGC4993 images (filter V throughout and 30s exposures), and one series for both Regulus (filters V and Clear, exposure time 1s) and the asteroid (both filters, exposure time 5s).

You also have dark field images (whose file names start with `dark`) for each exposure time (0.5s, 1s, 5s and 30s), and flat fields (whose file names start with `FLAT`) for each filter, taken with exposures of 0.5s.

Your B612 series in the V filter had 3 images starting from 0, but images 0 and 1 turned out defocalized, so you eliminated them.

Here is the directory structure you end up with:

```
├── 170817A.json
├── darks
│   ├── dark_1000_0.fits
│   ├── dark_1000_1.fits
│   ├── dark_1000_2.fits
│   ├── dark_1000_3.fits
│   ├── dark_30000_0.fits
│   ├── dark_30000_1.fits
│   ├── dark_500_0.fits
│   ├── dark_500_1.fits
│   ├── dark_5000_0.fits
│   ├── dark_5000_1.fits
│   └── dark_5000_2.fits
├── flats
│   ├── FLAT_Clear_500_0.fits
│   ├── FLAT_Clear_500_1.fits
│   ├── FLAT_V_500_0.fits
│   ├── FLAT_V_500_1.fits
│   └── FLAT_V_500_2.fits
└── objects
    ├── B612_1_Clear_5000_0.fits
    ├── B612_1_V_5000_2.fits
    ├── NGC4993_1_V_30000_0.fits
    ├── NGC4993_1_V_30000_1.fits
    ├── NGC4993_1_V_30000_2.fits
    ├── NGC4993_1_V_30000_3.fits
    ├── NGC4993_1_V_30000_4.fits
    ├── NGC4993_1_V_30000_5.fits
    ├── NGC4993_2_V_30000_0.fits
    ├── NGC4993_2_V_30000_1.fits
    ├── NGC4993_2_V_30000_2.fits
    ├── regulus_1_Clear_1000_0.fits
    ├── regulus_1_Clear_1000_1.fits
    ├── regulus_1_Clear_1000_2.fits
    ├── regulus_1_V_1000_1.fits
    └── regulus_1_V_1000_2.fits
```
Your configuration file 170817A.json will accordingly be:

```json
{
    "objects": ["NGC4993", "regulus", "B612"],
    "filters": ["V", "Clear"],
    "exposures": [1000, 5000, 30000],
    "flat_name": "FLAT",
    "dark_name": "dark"
}
```
For another example configuration file, please see the `example_config.json` file in the project directory.

## Reduction method
### Master dark images
The _master dark_ for a given exposure time is calculated as the pixel-wise __median__ of all the dark fields of that exposure. This allows to eliminate cosmic ray traces.

#### Dark field interpolation
In the case where some dark fields are missing, `astro_reduce` can interpolate the available master darks to obtain master darks for all the exposures necessary to reduce the object images. This is done by specifying the `--interpolate` option.

The interpolation is _least-square linear_, i.e., two images A and B are determined from the available master dark images such as to __minimize the square error__ on the linear interpolation (master dark) = (exposure time) x A + B.

Using these A and B, the missing master darks are calculated according to this equation.

>The FITS files for all the master dark images (deduced from dark fields or interpolated) can be found after reduction in the `tmp` directory under the names `mdark_[exposure].fits`.

### Master transmission images
The _master transmission_ image for a given filter is an image which encompasses the relative transmission of each pixel in the optical setup (telescope optics through filter to CCD matrix). It is calculated for every filter as the __median__ over all flat field images, after __subtraction of corresponding master dark images__ and __normalization__.

>The FITS files for all the master transmission images (deduced from dark fields or interpolated) can be found after reduction in the `tmp/` directory under the names `mtrans_[filter].fits`.

### Individual image reduction
An object image of given exposure and filter is reduced by __subtracting__ the corresponding exposure master dark image, and __dividing__ by the corresponding filter master transmission image.

>The FITS files for these images are found in the `tmp/` folder after reduction, with the extension `_aux.fits` (for _auxiliary_).

### Realignment and reduction
Finally, for each series of same exposure and filter for each object, the auxiliary files are realigned through optimization of their __mutual cross-correlations__, and then their pixel-wise __median__ image is calculated. Using the median rather than the mean allows to efficiently remove hot pixels. This will be all the more effective as dithering has been used in acquiring the images of a series.

>These are the final reduced images and can be found in the `reduced/` folder after reduction.

>Remaining hot pixels in the reduced images can be the result of insufficient dithering.

### Cross-series realignment
Optionally, `astro_reduce` may realign the reduced images of a same object across different series, provided the filters and exposures match. This is done with the `-c, -cross` option. The corresponding images can be found in the `reduced/` folder under the same name as the reduced files, except the series bit is no longer present.

>Beware of the interpretation of such cross-series images, as the series might well have been taken under different sky conditions.
