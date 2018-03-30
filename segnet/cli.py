
# -*- coding: utf-8 -*-

from __future__ import print_function
import click

from . import segnet
from . import utils
from collections import defaultdict, OrderedDict
from . import __logo_text__, __version__


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__, message=__logo_text__)
def main():
    """
    segNET

    A suite of tools for identifying gene modules from complex genome-scale networks


    Example usage: segnet diffuse -o outdir network_file source_genes

    segnet diffuse -v -o ../../data/segnet/p2test ../../data/segnet/brain_top_p2 ../../data/segnet/seedsChr2_2.txt

    """


@main.command('diffuse', options_metavar='[OPTIONS]', short_help='Run diffuse kernel')
@click.argument('netfile', metavar='netfile')
@click.argument('seedfile', metavar='seedfile')
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
    network=segnet.refine_network(netfile,idmap)
    if multiseeds:
        output = segnet.diffuse_multi_seeds(network, pos_seeds, neg_seeds, outdir)
        #segnet.draw_modules(network, network.nodes(), output)
        #print("after diffuse_multi_seeds, printing output")
        #print(output)
        #segnet.plot_histogram(output, outfile="test_histo.png")
        #print("after plot_histogram, now drawing network")
        segnet.draw_modules(network, network.nodes(), network.edges(),outdir)
        print("after draw_modules")
        segnet.plot_histogram_fdist(output, outfile="test_histo.png")
        segnet.plot_histogram_per_seed(output, outfile="test_histo_perseed.png")
    else:
        print("before segnet.diffuse_single_seed in cli.py")

        
        print("after diffuse_single_seed, printing output")
        #print(output)     
        segnet.plot_histogram_fdist(output, outfile="test_histo.png")
        #segnet.plot_histogram_per_seed(output, outfile="test_histo_perseed.png")
        print("after plot_histogram, now drawing network")
        segnet.draw_modules(network, network.nodes(), network.edges(), outdir)
        print("after draw_modules")

        #f = open(os.path.join(outdir, 'results.pickled', 'wb')
        #cPickle.dump(all_results, f)
        #f.close()
        #segnet.plot_histogram(output)
    


if __name__ == "__main__":
    main()
