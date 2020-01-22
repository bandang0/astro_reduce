'''Setup script for astro_reduce.'''
from setuptools import setup

setup(
    # General info
    name="astro_reduce",
    version="0.1",
    py_modules=["astro_reduce"],

    # Dependencies
    install_requires=["click",
        "astropy",
        "scipy",
        "numpy",
        "matplotlib"],

    # Main function trigger
    entry_points='''
        [console_scripts]
        astro_reduce=astro_reduce:cli
    ''',

    # Metadata to display on PyPI
    author="Raphaël Duque",
    author_email="duque@iap.fr",
    description="A Simple CCD Images Reducer for the Paris Observatory.",
    license="MIT",
    keywords="astronomy reduction images automatic",
    url="http://github.com/bandang0/astro_reduce",
)
