### Astroreduce -- A simple CCD images reducer for the Paris Observatory.

#### What does this do?

Given CCD flat field, dark field and astronomical object images, this program allows to produce images reduced in the standard fashion from the object images and using the dark and flat field. 
All it needs is a standardized organization of the folder in which lie the corresponding fits fils and a configuration file reporting the exposure times, filters, image sizes, and object names.
Along the reduction, Astroreduce will keep the intermediate images in a specific folder for later inspection and can optionally write png versions of these for easy checking of the overall reduction process. It is recommended to inspect theintermediate  images before considering the final reduced images.

#### Installing

Use the `setup.py` file with `python3 setup.py install --user` to install the program locally. You should then include `~/.local/bin` to your path in order to invoke `astro_reduce` from the console.

#### Usage

In the directory which contains the `.json` configuration file and is organized according the standard nomenclature (see the corresponding paragraph below), launch the reduction with:

`astro_reduce [OPTIONS] CONF_FILE`

You will then find the intermediate reduction files in the `tmp` directory, and your reduced images in the `reduced` directory.

##### Options

`--png`: Write png versions of the intermediate files. If on, these images will end up in a new `png` folder.
`--verbose`: Print additional reduction information to the console during run.

#### Nomenclature and configuration file

##### Initial folder structure

Astroreduce will operate in a folder which must have a `darks` folder, a `flats` folder, and a `objects` folder. An example of a valid folder structure and corresponding configuration file can be found in the Example paragraph below.

- The `darks` folder will contain only the dark field fits files, with names according to `[DARK\_NAME]\_[EXPOSURE\_TIME\_IN\_MS]\_[NUMBER].fits`. Here [DARK\_NAME] is the way you refer to darks (read in the configuration file, see next paragraph), and [NUMBER] is the number of the image in a series, if any.
- The `flats` folder will contain only the flat field fits files, with names according to `[FLAT\_NAME]\_[FILTER]\_[NUMBER].fits`. Here [FLAT\_NAME] is the way you refer to flats (also read in the configuration file), [FILTER] the filter used in acquisition (one of those read in the configuration file), and [NUMBER] is the number of the image in a series, if any.
- The `objects` folder will contain all the object images, with names according to `[OBJECT\_NAME]\_[SERIES\_NUMBER]\_[FILTER]\_[EXPOSURE\_TIME\_IN\_MS]\_[NUMBER].fits`. Here, [OBJECT\_NAME] is one of those read in the configuration file, and [SERIES\_NUMBER] is the number of the series, if any.

##### Configuration file

For a series of images to reduce, the configuration file must contain all the exposure times (in ms), filters, image sizes, object names and flat and dark image names as they appear in the different images file names.

All of this information is contained in a `json` file and organized as follows:

```json
{
    "objects": BRACKETED LIST OF OBJECT NAMES,
    "filters": BRACKETED LIST OF FILTERS,
    "exposures": BRACKETED LIST OF EXPOSURE TIMES IN MS,
    "flat_name": THE STRING USED IN IMAGE FILE NAMES TO REFER TO FLATS,
    "dark_name": THE STRING USED IN IMAGE FILE NAMES TO REFER TO DARKS,
    "nX": PIXEL SIZE OF IMAGES IN X AXIS,
    "nY: PIXEL SIZE OF IMAGES IN Y AXIS
}
```

Note: of course, all objects need not be imaged with every filter or every exposure time. It suffices to specify each expoure time and filter used once. Also, numbers of images in series need not start with 0 or 1 (see the example directory below).

##### Folder structure after reduction

After reduction, you will find two more folders in the directory: `tmp` and `reduced`.

#### Example

Say you observed on August 17th 2017 the NGC4993 galaxy (at two different moments of the night, hence two series), the B612 asteroid and Regulus to serve as a reference star for photometry of NGC4993. You will thus have two series of NGC4993 images (filter V throughout and 30s exposure time), and one series for both Regulus (filters V and Clear, exposure time 1s) and the astroid (both filters, exposure time 5s). 

You also have darks (which you call "dark") for each exposure time (1s, 5s and 30s), and flats (which you call "FLAT") for each filter. 

Your B612 series in the V filter had 3 images starting from 0, but images 0 and 1 turned out unfocalized, so you eliminated them. 

Here is the directory structure you end up with:

```
├── 170817A.json
├── darks
│   ├── dark_1000_0.fits
│   ├── dark_1000_1.fits
│   ├── dark_1000_2.fits
│   ├── dark_1000_3.fits
│   ├── dark_30000_0.fits
│   ├── dark_30000_1.fits
│   ├── dark_5000_0.fits
│   ├── dark_5000_1.fits
│   └── dark_5000_2.fits
├── flats
│   ├── FLAT_Clear_0.fits
│   ├── FLAT_Clear_1.fits
│   ├── FLAT_V_0.fits
│   ├── FLAT_V_1.fits
│   └── FLAT_V_2.fits
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
    "dark_name": "dark",
    "nX": 1024,
    "nY": 512,
}
```

#### Method
