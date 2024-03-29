'''Setup script for astro_reduce.'''

from setuptools import setup


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    # General info
    name="astro_reduce",
    version="4.1.0",
    py_modules=["astro_reduce", "helpers", "cosmetic"],

    # Dependencies
    install_requires=["click",
                      "astropy",
                      "scipy",
                      "numpy",
                      "matplotlib"],

    # Python version
    python_requires=">=3.0",

    # Main function trigger
    entry_points='''
        [console_scripts]
        astro_reduce=astro_reduce:cli
    ''',

    # We have data files packaged, for the astrometry features
    data_files=[('data', ['data/t120.sex',
                          'data/t120.param',
                          'data/t120.psfex',
                          'data/t120.scamp',
                          'data/t120-andor1024.ahead',
                          'data/t120.parampsfex',
                          'data/default.conv'])],

    # Metadata to display on PyPI
    author="Raphaël Duque",
    author_email="duque@physik.uni-frankfurt.de",
    url="https://github.com/bandang0/astro_reduce",
    description="A simple CCD-image reducer for the Observatoire de Paris-Meudon.",
    long_description_content_type="text/markdown",
    long_description=long_description,
    license="MIT",
    keywords="astronomy reduction images automatic astrometry",
    classifiers=[
        "Environment :: Console",
        "Programming Language :: Python :: 3.0",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Astronomy"
    ],
)
