#!/usr/bin/env python

'''
Author: Marek Cmero
Takes a feature matrix, runs the incomplete phylogeny algorithm on it,
and if a perfect phylogeny is found, outputs a postscript file plotting the tree.
(Graphviz must be installed for plotting.)
'''

import argparse
import numpy as np
from perfect_phylogeny import phylogeny as phyl

parser = argparse.ArgumentParser(prefix_chars='--')
parser.add_argument('m_file',help='Matrix file')
parser.add_argument('outname',nargs='?',default='tree',help='Output name')
parser.add_argument('--plot',dest='plot',action='store_true',help='Plot output (if valid phylogeny). Requires Graphviz.')
args = parser.parse_args()

m_file  = args.m_file
outname = args.outname
plot    = args.plot

if __name__ == '__main__':

    m_txt = np.genfromtxt(m_file,delimiter='\t',dtype=None,invalid_raise=False)
    
    s = np.array(m_txt[1:,0],dtype='S100') # row names/samples
    c = np.array(m_txt[0,1:],dtype='S100') # column names/features
    m = np.array(m_txt[1:,1:],dtype=int)   # matrix values

    m_new, c_prime, k1, tree = phyl.solve_incomplete_phylogeny(m,s,c)
    k1 = phyl.get_k1_matrix(m_new,c_prime)

    if phyl.perfect_phylogeny_exists(k1,c_prime) and plot:
        print('Plotting tree...')
        phyl.write_dot(k1,s,outname)
