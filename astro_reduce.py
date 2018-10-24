'''astro_reduce -- A Simple CCD images reducer for the Paris Observatory.'''

import click

@click.command()
@click.argument('conf_file', type=click.File('r'))
@click.option('--png', '-p', is_flag=True,
        help='Write png format of intermediary images in png folder.')
@click.option('--verbose', '-v', is_flag=True,
        help='Enables verbose mode.')
def cli(conf_file, png, verbose):
    '''Reduce CCD images from objects with flat and dark field images.'''
    click.echo(conf_file)
    click.echo(png)
    click.echo(verbose)
