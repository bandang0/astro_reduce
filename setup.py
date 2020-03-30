'''Setup script for astro_reduce.'''
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    # General info
    name="astro_reduce",
    version="2.2.1",
    py_modules=["astro_reduce"],

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

    # Metadata to display on PyPI
    author="Raphaël Duque",
    author_email="duque@iap.fr",
    url="https://github.com/bandang0/astro_reduce",
    description="A Simple CCD Images Reducer for the Paris Observatory.",
    long_description_content_type="text/markdown",
    long_description=long_description,
    license="MIT",
    keywords="astronomy reduction images automatic",
    classifiers=[
        "Environment :: Console",
        "Programming Language :: Python :: 3.0",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Astronomy"
    ],
)
