# -*- coding: utf-8 -*-
import sys
from collections import defaultdict, OrderedDict
import os
import math
import networkx as nx
import numpy as np
#import scipy
#import scipy.sparse
#from scipy.sparse.linalg import minres
import matplotlib.pyplot as pp
import csv
#import pickle
#import gzip
import json
import pandas as pd
#import scipy.sparse as sps

edgelist = set()
nodelist = set()
pathlist = set()
H = nx.Graph()

diff_file = sys.argv[1]
vote_file = sys.argv[2]
region1 = str(sys.argv[3])
region2 = str(sys.argv[4])
fh = open(diff_file, 'r')
fh2 = open(vote_file, 'r')

for line2 in fh2:
    items = line2.split('\t')
    source = str(items[0].strip())
    target = str(items[1].strip())
    weight = float(items[2].strip())
    pathlist.add((source,target,weight))
    
#print(pd.read_json(fh.read()))
mat_dict = dict()
mat = json.load(fh)
for key,results in mat.iteritems():
 #   print(key)
    nodelist.add(str(key))
    for r in results.keys():
        if float(results[r][0]) > 0.0:
            edgelist.add((str(key),str(r),1/float(results[r][0])))
            nodelist.add(str(r))
H.add_nodes_from(nodelist)
H.add_weighted_edges_from(edgelist)


shortest_paths_output = list()

for tup in pathlist:
    shortest_path = nx.bidirectional_dijkstra(H,tup[0],tup[1], weight = 'weight')
    
    shortest_paths_output.append(shortest_path)
    #for path in list(shortest_paths):
#        shortest_paths_output.add(path)
df = pd.DataFrame(shortest_paths_output)
df.to_csv('inv_bidijkstra_%s_%s.txt'% (region1,region2), sep='\t',header=False,index=False)
