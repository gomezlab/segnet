# -*- coding: utf-8 -*-

from . import utils
from collections import defaultdict, OrderedDict
import os
import math
import networkx as nx
import numpy as np
import scipy
import scipy.sparse
from scipy.sparse.linalg import minres
import matplotlib.pyplot as pp
import csv
#import gzip

import pandas as pd
import scipy.sparse as sps

LOG = utils.get_logger()


#
# Functions for loading the whole genome network
#

def get_idmap(idfile, key_column, value_column, header=False):
    """
    get_dict_from_file
    ------------------
    Use when you want to get the mappings between id's, e.g.,from MGI to Entrez Gene ID

    Args:
        idfile: Gene id file (tab-separated values)
        key_column: column index of keys
        value_column: column index of values
        header: True if datafile has a header row
    Returns:
        gene_id_map (OrderedDict)

    """
    with open(idfile) as fh:
        if header is True:
            fh.next()
        mydict = OrderedDict()
        for curline in fh:
            curline = curline.rstrip()
            items = curline.split("\t")
            key = items[key_column]
            value = items[value_column]
            if mydict.has_key(key):
                if not mydict[key] == value:
                    LOG.debug("Conflicting items: key=%s, registered value=%s <=> new value=%s" %
                              (key, mydict[key], value))
            else:
                mydict[key] = value
    return mydict


def convert_ids(id_list, idmap):
    """
    `convert_ids` convert id's of input genes

    :param id_list:
    :param idmap: id conversion table (OrderedDict)
    :return: a list that contains the converted id's

    """
    new_list = []
    for gid in id_list:
        if idmap.has_key(gid):
            new_list.append(idmap[gid])
        else:
            print "Error: " + gid + " does not exist."
    return new_list


def remove_nonexistent_ids(id_list, G):
    """
    `remove_nonexistent_ids` removes genes that do not exist in the network from a set of genes

    :param id_list:
    :param G:
    :return:
    """
    new_list = []
    for gid in id_list:
        if G.has_node(gid):
            new_list.append(gid)
    return new_list


def prepare_seeds(seedfile):
    with open(seedfile) as fh:
        pos_seeds = OrderedDict()
        neg_seeds = OrderedDict()
        for curline in fh:
            item = curline.rstrip().split('\t')
            seed_gid = item[0]
            seed_val = float(item[1])
            if seed_val > 0:
                pos_seeds[seed_gid] = seed_val
            else:
                neg_seeds[seed_gid] = seed_val
    return pos_seeds, neg_seeds


def readin_network(netfile, idmap=None, beta=None, header=False):
    """
    `readin_network` constructs networkx graph based upon input network file

    :param netfile: the name of whole genome network file (tab-delimited string)
    :param idmap: a conversion table (dictionary) of id's (dict)
    :param beta: a parameter used in gaussian kernel (float)
    :param header: whether the network file contains header or not (bool)
    :return: a whole genome network (networkx graph)

    :usage:
        gfile = '../data/mouseNET/mouseNET_out_normal_10'
        G = get_network(gfile, idmap=idmap, beta=1.0)

    """
    LOG.info("Loading the network from %s" % netfile)

    G = nx.DiGraph()
    #with gzip.open(netfile) as fh:
    with open(netfile) as fh:
        if header is True:
            fh.next()
        for curline in fh:
            items = curline.rstrip().split("\t")
            node1 = items[0]
            node2 = items[1]
            if node1 != node2:  # Not allowing self edges
                # Convert node id whenever it's desired
                if idmap is not None:
                    if idmap.has_key(node1):
                        node1 = idmap[node1]
                    if idmap.has_key(node2):
                        node2 = idmap[node2]
                if beta is not None:
                    edge_weight = math.exp(-beta * math.pow((1.0 - float(items[2])), 2)) + 1e-6
                else:
                    edge_weight = float(items[2])
                G.add_edge(node1, node2, {'weight': edge_weight})
    return G


#
# Functions for obtaining diffusion seeds
#

def sort_genes_by_coordinate(gene_list, gene_locs):
    """
    sort_genes_by_coordinate
    ------------------------
    Args:
    :param gene_list: a list of gene_name
    :param gene_locs: a dict of {gene_name: [chromosome, gene_start, gene_end]}

    :return: ordered_list of [ [gene_start, gene_end, gene_name] ]


    """
    ordered_list_final = []
    chrom = [None] * 18
    for gname in gene_list:
        chro, start, end = gene_locs[gname]
        if chrom[chro] is None:
            chrom[chro] = []
        if start <= end:
            chrom[chro].append([start, gname])
        else:
            chrom[chro].append([end, gname])
    for i in xrange(1, len(chrom)):
        if chrom[i] is not None:
            chrom[i].sort()
            ordered_loc, ordered_list = zip(*chrom[i])
            ordered_list_final += list(ordered_list)
    return ordered_list_final


def get_associated_markers(markers_and_pvals, threshold=None):
    """
    get_associated_markers

    :param gname: eQTL target gene name
    :param threshold: a threshold p-value or LOD score

    :return: a list of marker id's that satisfy the p-value threshold

    """
    if threshold is not None:
        associated_markers = []
        for k, v in markers_and_pvals.iteritems():
            if v[1] <= threshold:
                associated_markers.append(k)
    else:
        associated_markers = markers_and_pvals.keys()
    return associated_markers


def get_genes_near_markers(chromosomes, markers, gene_locs, mids_wanted, brange=10000):
    """
    `get_genes_near_markers` returns set of systematic names given a list of markers

    :param chromosomes: a list of [ [gene_start, gene_end, gene_name] ]
    :param markers: a dict of {marker_id: marker_position (chromosome, pos)}
    :param gene_locs: a dict of {gene_name: [chromosome, gene_start, gene_end]}
    :param mids_wanted: a list of marker ids
    :param brange:

    :return: a list of genes that are located near the given markers

    Note:
        Use when getting diffusion origins (regulated genes) and target genes.

    """
    gene_set = set()
    for mid in mids_wanted:
        chro, loc = markers[mid]
        range_start = loc - brange
        range_end = loc + brange
        all_gene_locs = np.array(chromosomes[chro])[:, 0:2]
        all_gene_locs = np.array(all_gene_locs, dtype='int')
        nearby_genes = np.nonzero((all_gene_locs[:, 0] <= range_end) & (all_gene_locs[:, 1] >= range_start))[0]
        new_genes = np.array(chromosomes[chro])[nearby_genes]
        gene_set = gene_set.union(set(new_genes[:, 2]))
    return list(gene_set)


def get_complement_nodeset(G, diffusion_field):
    """
    get_complement_nodeset

    :param G:
    :param diffusion_field:
    :return:
    """
    all_nodes = set(G.nodes())
    return list(all_nodes.difference(set(diffusion_field)))


def check_neighbors(nonfield_dict, neighbors, check_type='all_pos'):
    """
    check_neighbors
    ---------------

    Count neighboring positive seed nodes (or negative seed nodes)

    :param nonfield_dict:
    :param neighbors:
    :param check_type:
    :return:

    """
    neg_cnt = 0
    pos_cnt = 0
    for k, v in neighbors.iteritems():
        if nonfield_dict.has_key(k):
            neg_cnt += 1
        else:
            pos_cnt += 1
    if check_type == 'all_pos':
        return (neg_cnt == 0)
    elif check_type == 'all_neg':
        return (pos_cnt == 0)


def set_network_nodevals(nodefile, key_column, value_column, G, header=False):
    """
    set_network_nodevals

    :param nodefile:
    :param key_column:
    :param value_column:
    :param G:
    :param header:
    :return:

    :usage: Use to read in node_id and its value
        bxd_datafile = '../data/rusyn/BXD_all_nodevals'
        all_nodevals = set_network_nodevals(bxd_datafile, 0, 1)

    """
    with open(nodefile) as fh:
        if header is True:
            fh.next()
        dict_nodevals = {}
        for curline in fh:
            items = curline.rstrip().split("\t")
            nodeid = items[key_column]
            value = items[value_column]
            if G.has_node(nodeid):  # Not allowing nodes that do not exist in the network
                if not dict_nodevals.has_key(nodeid):  # Not allowing duplicate values
                    dict_nodevals[nodeid] = float(value)
    return dict_nodevals


def get_network_nodevals(gidfile, nodevals, gid_column):
    """
    get_network_nodevals

    :param gidfile:
    :param nodevals:
    :param gid_column:
    :return:

    :usage: Use to read out node values
        gidfile = '../data/rusyn/MDP_chr12_transband_genes'
        tb_nodevals = get_network_nodevals(gidfile, all_nodevals, 1)

    """
    extracted_nodevals = {}
    with open(gidfile) as fh:
        for curline in fh:
            items = curline.rstrip().split("\t")
            gid = items[gid_column]
            if nodevals.has_key(gid):
                extracted_nodevals[gid] = nodevals[gid]
    return extracted_nodevals


def assign_indices_to_seeds(posseeds_nodevals, negseeds_nodevals):
    """
    `assign_indices_to_seeds`

    :param posseeds_nodevals:
    :param negseeds_nodevals:
    :return:

    """
    indmap_seeds = {}
    nodeval_seeds = []
    curindex = 0

    # Register all positive & negative seed nodes
    for key, nval in posseeds_nodevals.iteritems():
        if not indmap_seeds.has_key(key):  # Not allowing duplicate values
            indmap_seeds[key] = curindex
            nodeval_seeds.append(nval)
            curindex += 1
    for key, nval in negseeds_nodevals.iteritems():
        if not indmap_seeds.has_key(key):  # Not allowing duplicate values
            indmap_seeds[key] = curindex
            nodeval_seeds.append(nval)
            curindex += 1
    return indmap_seeds, nodeval_seeds


def assign_indices_to_nonseeds(indmap_seeds, G):
    """
    `get_nonseeds` assigns indices to non-seed nodes

    :param indmap_seeds: index_map_of_all_seeds (dict)
    :param G: whole genome network
    :return: index_map_of_nonseed_nodes (OrderedDict)

    """
    indmap_nonseeds = OrderedDict()
    curindex = 0
    for curnode in G.nodes_iter():
        if not indmap_seeds.has_key(curnode) and \
           not indmap_nonseeds.has_key(curnode) and \
           G.degree(curnode) > 0:
            indmap_nonseeds[curnode] = curindex
            curindex += 1
    return indmap_nonseeds


#
# Functions for preparing or running a diffusion kernel
#

def construct_laplacian(G, indmap_seeds, nodeval_seeds, indmap_nonseeds):
    """
    `construct_laplacian` create weighted Laplacian (L_U and minusB_T separately)

    :param G:
    :param indmap_seeds:
    :param nodeval_seeds:
    :param indmap_nonseeds:
    :return:

    TODO:
        Think about how to use different classes (Maybe for dealiing with many TB's)

    """
    numseeds = len(indmap_seeds)
    numnonseeds = len(indmap_nonseeds)
    numclasses = 1

    X_M = scipy.zeros([numseeds, numclasses])
    for i in xrange(numseeds):
        X_M[i, 0] = nodeval_seeds[i]

    L_U = scipy.sparse.lil_matrix((numnonseeds, numnonseeds))
    minusB_T = scipy.zeros([numnonseeds, numseeds])
    for node1, node2 in G.edges_iter():
        if indmap_seeds.has_key(node1):
            node1_is_seed = True
            node1ind = indmap_seeds[node1]
        else:
            node1_is_seed = False
            node1ind = indmap_nonseeds[node1]
        if indmap_seeds.has_key(node2):
            node2_is_seed = True
            node2ind = indmap_seeds[node2]
        else:
            node2_is_seed = False
            node2ind = indmap_nonseeds[node2]

        edge_weight = G[node1][node2]['weight'] #should be the W matrix
        
        # print edge_weight, node1, node1ind, '<=>', node2, node2ind, float(edge_weight)

        if not node1_is_seed and not node2_is_seed:
            L_U[node1ind, node2ind] = -edge_weight
            L_U[node2ind, node1ind] = -edge_weight
            L_U[node1ind, node1ind] += edge_weight
            L_U[node2ind, node2ind] += edge_weight
        if node1_is_seed and not node2_is_seed:
            L_U[node2ind, node2ind] += edge_weight
            minusB_T[node2ind, node1ind] = edge_weight
        if not node1_is_seed and node2_is_seed:
            L_U[node1ind, node1ind] += edge_weight
            minusB_T[node1ind, node2ind] = edge_weight

    L_U.tocsr()
    #print("printing G the similarity graph")
    #print(list(G.edges_iter()))
    
    return X_M, L_U, minusB_T
#X_M: just the seed values 1*1, per seed
#L_U: keys are node-node, values are +/- edge weights
#minusB_T: the original seed-node values column vector


def solve_diffusion_eqn(X_M, L_U, minusB_T, tol=0.00001, maxiter=500000): #change tol to change number of steps
    #tol is just a range of telling if the diffusion needs to go further
    """
    `solve_diffusion_eqn`

    :param X_M:
    :param L_U:
    :param minusB_T:
    :param tol:
    :param maxiter:
    :return:

    """
    a = OrderedDict()
    b = scipy.dot(minusB_T, X_M).flatten()
    X_U, info = scipy.sparse.linalg.minres(L_U, b, tol=tol, maxiter=maxiter)

    L = L_U.toarray()
    n, m = L.shape
    diags = L.sum(axis=1)
    D = sps.spdiags(diags.flatten(), [0], m, n, format='csr')
    S = D - L
    P = D * S * D
    eigval, eigvec = np.linalg.eig(P)

    #print("check X_U")
    #print(X_U)
    #print(X_U.shape)
    
    #print("similarity matrix S")
    #print(S)
    #print("laplacian matrix L")
    #print(L)
    #print("doubly stochastic symmetric matrix P")
    #print(P)
    #X_M_array = X_M.toarray()
    #L_U_array = L_U.toarray()
    #print("L_U........................")
    #print(L_U_array[np.nonzero(L_U_array)])
    #print("minusB_T,,,,,,,,,,,,,,,,,,,")
    #print(minusB_T)
    #print("shapes:::::::::::::::::::::")
    #print(X_M.shape)
    #print(L_U_array.shape)
    #print(minusB_T.shape)
    #print(X_U)
    return X_U, P, eigval # X_U: row vectors - diffusion scores, per seed

    #L_U is actually the similarity matrix? def not the graph laplacian bc it is not postive semi-def

def get_diffusion_profile(indmap_allseeds, nodeval_allseeds, nonseeds, X_U):
    """
    `report_diffusion_profile` organizes diffusion profile

    :param indmap_allseeds:
    :param nodeval_allseeds:
    :param nonseeds:
    :param X_U:
    :return: a dict of {gene name: [diffusion value, if_the_gene_was_seed]}  # 0/1 for seeds/nonseeds

    Note:
        Use the following line to get the names of top scoring genes
        sorted_keys = sorted(result, key=result.__getitem__, reverse=True)

    """
    result = OrderedDict()
    for node_name, node_index in indmap_allseeds.items():
        result[node_name] = [nodeval_allseeds[node_index], 1]

    for i in xrange(len(X_U)):
        result[nonseeds[i]] = [X_U[i], 0]
    return result #ordered dictionary: diffusion profile
                  #OrderedDict([('A', [1.0, 1]), ('1', [0.9921955980991151, 0])])


def diffuse_multi_seeds(G, pos_seeds, neg_seeds, outdir):
    """
    `diffuse_multiseeds` runs diffusion kernel on the whole genome network

    :param G:
    :param pos_seeds:
    :param neg_seeds:
    :param outdir:
    :return:

    """
    indmap_seeds, nodeval_seeds = assign_indices_to_seeds(pos_seeds, neg_seeds)  # both pos & neg
    indmap_nonseeds = assign_indices_to_nonseeds(indmap_seeds, G)
    nonseeds = indmap_nonseeds.keys()

    X_M, L_U, minusB_T = construct_laplacian(G, indmap_seeds, nodeval_seeds, indmap_nonseeds)
    
    X_U, P, eigval = solve_diffusion_eqn(X_M, L_U, minusB_T)
    #print("=====+++++=~~~~~~===")

    
    #print(eigval)
    #
    # Report the result
    #
    result = get_diffusion_profile(indmap_seeds, nodeval_seeds, nonseeds, X_U)
    if outdir is not None:
        with open(os.path.join(outdir, 'diffusion_profile_%s.tsv'), 'w') as fh:
            for node_name, node_info in result.iteritems():
                fh.write("%s\t%.4f\t%d\n" % (node_name, node_info[0], node_info[1]))
        return result, eigval #added this line to have diffusion_profile saved as well as result returned
    else:
        return result, eigval


def diffuse_single_seed(G, pos_seeds, neg_seeds, outdir):
    """
    `diffuse_single_seed`

    :param G:
    :param pos_seeds:
    :param neg_seeds:
    :param outdir:
    :return:
    """
    all_results = dict()
    for src_node_name, src_node_val in pos_seeds.iteritems():
        result, eigval = diffuse_multi_seeds(G, {src_node_name: src_node_val}, neg_seeds, outdir)
        all_results[src_node_name] = result
    if outdir is not None:
        for src_node_name, result in all_results.iteritems():
            with open(os.path.join(outdir, 'diffusion_profile_%s.tsv' % src_node_name), 'w') as fh: #'diffusion_profile_%s.tsv'
                for node_name, node_info in result.iteritems():
                    fh.write("%s\t%.4f\t%d\n" % (node_name, node_info[0], node_info[1]))
        return all_results #added this line to have diffusion_profile saved as well as all_results returned
        # for src_node_name, result in all_results.iteritems():
        #     with open(os.path.join(outdir, 'diffusion_profile.%s.tsv' % src_node_name), 'w') as fh:
        #         for node_name, node_info in result.iteritems():
        #             fh.write("%s\t%.4f\t%d" % (node_name, node_info[0], node_info[1]))

    else:
        return all_results


#
# Functions for postprocessing diffusion results
#

def get_modules(nlist, G, cutoff=None, beta=None):
    """
    get_subnetwork
    --------------

    """
    H = nx.DiGraph() #changed from Graph() to DiGraph()
    H.add_nodes_from(nlist)
    for nid1 in nlist:
        for nid2 in nlist:
            if nid1 is not nid2 and not H.has_edge(nid1, nid2) and G.has_edge(nid1, nid2):
                w = G[nid1][nid2]['weight']
                if beta is not None:
                    w = math.exp(-beta * math.pow((1.0 - w), 2)) + 1e-6
                if cutoff is None or (cutoff is not None and w > cutoff):
                    H.add_edge(nid1, nid2, {'weight': w})
    return H


def draw_modules(H, nodelist, edgelist, outdir, nodelabels=None):
    """
    draw_subnetwork

    :param H:
    :param nodelist: all the nodes that we want to draw
    :param notables: nodes that we want to draw with two concentric circles
    :param nodelabels: a dictionary to convert node id's to other desired names
    :return:

    """
    # pos = nx.graphviz_layout(H)
    # nx.draw_networkx_edges(H, pos, alpha=0.5, edge_color='y')
    # if notables is not None:  # if there exist seed genes
    #     nx.draw_networkx_nodes(H, pos, nodelist=notables, node_color='white', node_size=300, alpha=0.7)
    # nx.draw_networkx_nodes(H, pos, nodelist=nodelist, node_color='lightblue', node_size=100, alpha=0.5)
    # if nodelabels is not None:
    #     nx.draw_networkx_labels(H, pos, labels=nodelabels, font_size = 6, font_weight='bold')
    # else:
    #     nx.draw_networkx_labels(H, pos, font_size = 6, font_weight='bold')

    # Test version
    # problem in the following is if nodes don't exist in a particular version
    # of the network, then this all breaks as it doesn't have a position
    
    # spectral_pos = nx.drawing.layout.spectral_layout(H, weight='weight', dim=2) #spectral layout didnt work; too many nodes (>500)
    pp.figure()

    # may want to check here if all nodes that we want to visualize are in network
    notables=nodelist
    spring_pos = nx.spring_layout(H, dim=2, k=0.2, pos=None, fixed=None, iterations=50, weight='weight', scale=1.0, center=None)
    
    # if notables is not None:  # if there exist seed genes
    nx.draw_networkx_nodes(H, spring_pos, nodelist=nodelist, node_color='red', node_size=20, alpha=0.3) #the outer circle of the nodes
#    nx.draw_networkx_nodes(H, spring_pos, nodelist=nodelist, node_color='lightblue', node_size=300, alpha=0.5) #the inner circle of the nodes
    nx.draw_networkx_edges(H, spring_pos, edgelist=edgelist, alpha=1.0, edge_color='y')
    if nodelabels is not None:
        nx.draw_networkx_labels(H, spring_pos, labels=nodelabels, font_size = 5)
   # else:
#        nx.draw_networkx_labels(H, spring_pos, labels=None, font_size = 5)
    pp.savefig('test_testnetwork.png')


def plot_histogram_fdist(result, outfile='histogram.node_values.png'):
    #
    # Plot the frequency distribution of node values / diffusion scores
    # Decompose the nested "result" to v_dict = orderedDict[seed-node combination: diffusion score]
    fig = pp.figure()
    v_dict=OrderedDict() 
    print('printing ----------------RESULT')
    for k,v in result.iteritems():
        for k_,v_ in v.iteritems():
            key=k+str(k_)
            v_dict[key]=v_[0]
    print('printing ----------------HISTO')
    pp.hist(x=v_dict.values(),bins=100,normed=False, rwidth=0.8, facecolor='red', alpha=0.7)
    #print(pp.hist(v_dict.values(),bins=100,normed=False, rwidth=0.8, facecolor='red', alpha=0.7))
    pp.savefig(outfile, dpi=300)
    #print (result)

    
def get_diffusion_profile_matrix(origins, targets, G, pickle_loc='.'):
    """
    `get_diffusion_profile_matrix`

    :param origins:
    :param targets:
    :param G:
    :param pickle_loc:
    :return:
    """
    mat = np.zeros((len(origins), len(targets)))

    for i in xrange(len(origins)):
        ogname = origins[i]
        if G.has_node(ogname):
            reader = csv.reader(open(os.path.join(pickle_loc, 'tests/diffusion_profile_%s.tsv' % ogname)))
            diffusion_profile = dict()
            for row in reader:
                r = row[0].split('\t')
                diffusion_profile [r[0]] = [r[1],r[2]]
                #print(diffusion_profile)
                for j in xrange(len(targets)):
                    tgname = targets[j]
                    if diffusion_profile.has_key(tgname):
                        #print(tgname)
                        mat[i, j] = diffusion_profile[tgname][0]
                    #else:
                        #print "%s does not exist in the diffusion field." % tgname #not efficient output
        else:
            print "%s does not exist in the network." % ogname
    return mat

def origins_and_targets(result):
    origins = list()
    targets = list()
    for k,v in result.iteritems():
        origins.append(k)
    for k,v in result.iteritems():
        for k_,v_ in v.iteritems():
            if k_ not in targets and k_ not in origins: #targets exclude seeds
                targets.append(k_)
    return origins, targets


def count_votes(mat, origins, targets, effective_vote_val=0.0):
    """
    `count_votes`

    :param mat:
    :param origins:
    :param targets:
    :param effective_vote_val:
    :return:

    """
    # Counting the top voter only
    vote_cnt = np.zeros((1, len(targets))).flatten()
    vote_data = defaultdict(list)
    vote_score = defaultdict(list) #prob can be used in visualization (stacked weighted distribution graph)
                                   #after knowing the gene IDs / the properties
    vote = mat.argmax(axis=1) #index of the max element (comparisons within one diffusion profile)
    vote_val = mat.max(axis=1) # axis = 1 --> max diffusion score per target
    for i in xrange(len(vote)):
        if effective_vote_val < vote_val[i]:    # Assuming only the origin has the value of 1 (already excluded)
                                                # some non-seeds have score over 1 as well
            # Record the name of the voters
            score = mat[i, vote[i]]
            vote_score[targets[vote[i]]].append([origins[i],score])
            vote_data[targets[vote[i]]].append(origins[i])
            vote_cnt[vote[i]] += 1
    vote_data = dict(vote_data)
    return vote, vote_cnt, vote_val, vote_data
    # vote_data =  {'70110': ['57438', '13555', '228005', '107435', '229279', '11909', '20823', '16403', '13669', '13641', '13642'],
        #'19084': ['13627'], '53895': ['13629'], '93747': ['13653', '228033'], '13200': ['70231'], '54396': ['71586']}]
    #i.e. seeds vote to the nodes



def vote_visual (G, vote_data, outfile = 'vote_network.png', beta = None, cutoff = None, nodelabels = None):
    nodelist = list()
    for k,v in vote_data.iteritems():
        nodelist.append(k)
        for v_ in v:
            nodelist.append(v_)
    H = nx.DiGraph() #DiGraph() instead of Graph()
    #print(nodelist)
    H.add_nodes_from(nodelist)
    for votee, voters in vote_data.iteritems():
        for voter in voters:
            if votee is not voter and not H.has_edge(votee, voter):
                H.add_edge(votee, voter, {'weight':1}) # = unweighted edges. Not sure what should the weighted edges be?


    #print(H.edges())                  
    pos = nx.drawing.layout.spectral_layout(H)
    pp.figure()
    nodelist = H.nodes()
    notables=nodelist
    spring_pos = nx.spring_layout(H, dim=2, k=None, pos=None, fixed=None, iterations=50, weight='weight', scale=1.0, center=None)
    
    if notables is not None:  # if there exist seed genes
        nx.draw_networkx_nodes(H, spring_pos, nodelist=nodelist, node_color='red', node_size=400, alpha=0.5) #the outer circle of the nodes
        nx.draw_networkx_nodes(H, spring_pos, nodelist=nodelist, node_color='lightblue', node_size=300, alpha=0.5) #the inner circle of the nodes
        nx.draw_networkx_edges(H, spring_pos, H.edges(), alpha=0.75, edge_color='y')
    if nodelabels is not None:
        nx.draw_networkx_labels(H, spring_pos, labels=nodelabels, font_size = 6)
    else:
        nx.draw_networkx_labels(H, spring_pos, labels=None, font_size = 6)
    pp.savefig('vote_network.png')


# uses result from the diffusion output after single_seed, which is an ordered dictionary seed: node, {score, seed boolean}
#def spectral_clustering(G,result,indmap_seeds, nodeval_seeds, indmap_nonseeds, origins, targets):
def spectral_clustering(G, result, origins, targets, pickle_loc='.'):

    """
    `spectral_clustering` create weighted Laplacian (L_U and minusB_T separately)

    :param G:
    :param result:
    :param indmap_seeds:
    :param nodeval_seeds:
    :param indmap_nonseeds:
    :return:

    TODO:
        L = D - W
        use result to create W matrix and D matrix
        get L
        confirm L is positive semi-definite (symm. obv.)
        comput eigenvalues and eigenvectors of L
        cluster eigenvalues using k-means
    """

    nodes = origins + targets
    W = scipy.zeros([len(nodes), len(nodes)])
    O = scipy.ones([len(nodes), len(nodes)])
    for i in xrange(len(nodes)):
        if G.has_node(nodes[i]):
            
            file = os.path.join(pickle_loc, 'tests/diffusion_profile_%s.tsv' % nodes[i])
            if os.path.isfile(file):
                reader = csv.reader(open(file))
                for row in reader:
                    r = row[0].split('\t')
                    j = nodes.index(r[0])
                    W [i,j] = r[1]
    for i in xrange(len(origins)):
        for j in xrange(len(origins)):
            W[i,j] = (W[i,j] + W[j,i])/2
            W[j,i] = 0

#    for i in xrange(len(nodes)):
#        for j in xrange(len(nodes)):
#            W[j,i] = 0


 #   D = W * O
#    P = np.linalg.inv(D) * W
    print("goddamn it please let this be the W matrix")
    print(W)
    print(D)
    print(P)
