# -*- coding: utf-8 -*-

from __future__ import print_function
import click

from . import segnet
from . import utils
from . import __logo_text__, __version__

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__, message=__logo_text__)
def main():
    """
    segNET

    A suite of tools for identifying gene modules from complex genome-scale networks

    """


@main.command('diffuse', options_metavar='[OPTIONS]', short_help='Run diffuse kernel')
@click.argument('netfile', metavar='netfile', type=click.File())
@click.argument('seedfile', metavar='seedfile', type=click.File())
@click.option('-m', '--multiseeds', is_flag=True)
@click.option('-c', '--gidfile', type=click.File(mode='r'), help="Gene ID file")
@click.option('-o', '--outdir', type=click.Path(exists=True, writable=True), help="Output directory")
@click.option('-v', '--verbose', count=True, help='The more times listed, the more output')
def diffuse(netfile, seedfile, multiseeds, gidfile, outdir, verbose):
    """
    Runs the diffusion kernel from each of the specified source (or positive seed) genes

    :param netfile:
    :param seedfile:
    :param multiseeds:
    :param gidfile:
    :param outdir:
    :param verbose:
    :return:
    """
    utils.configure_logging(verbose)
    if gidfile:
        idmap = segnet.get_idmap(gidfile, 0, 1)
    else:
        idmap = None
    pos_seeds, neg_seeds = segnet.prepare_seeds(seedfile)
    network = segnet.readin_network(netfile, idmap)
    if multiseeds:
        output = segnet.diffuse_multi_se(network, pos_seeds, neg_seeds, outdir)
    else:
        output = segnet.diffuse_single_seed(network, pos_seeds, neg_seeds, outdir)


if __name__ == "__main__":
    main()
