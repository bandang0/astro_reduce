'''Setup script for astro_reduce.'''
from setuptools import setup

setup(
    # General info
    name="astro_reduce",
    version="0.1",
    scripts=['astro_reduce'],

    # Dependencies
    install_requires=["astropy", "scipy", "pandas", "matplotlib"],

    # metadata to display on PyPI
    author="Raphaël Duque",
    author_email="duque@iap.fr",
    description="A simple CCD images reducer for the Paris Observatory.",
    license="MIT",
    keywords="astronomy reduction images automatic",
    url="http://github.com/bandang0/astro_reduce", # project home page, if any
)
