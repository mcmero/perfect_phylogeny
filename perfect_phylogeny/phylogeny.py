'''
Author: Marek Cmero
Functions for perfect phylogeny and incomplete phylogeny algorithms
'''

import numpy as np
import copy
import string
import subprocess
from mgraph import MGraph

def get_duplicates(items):
    ''' 
    returns the indices of all duplicates in a list
    @param items a 1D list
    '''
    locs = []

    for i in range(len(items)):
        loc_tmp = [i]
        item = items[i]
        start_at = i+1
        while True:
            try:
                loc = items.index(item,start_at+1)
            except ValueError:
                break
            else:
                loc_tmp.append(loc)
                start_at = loc
        
        if len(loc_tmp)>1:
            locs.append(loc_tmp)

    return(locs)

def remove_duplicates(m_prime, c_prime):
    '''
    remove any duplicate columns and merge associated features
    @param m_prime M' matrix
    @param c_prime features list corresponding to columns of M'
    '''
    m_prime_tmp = map(lambda x: '.'.join(map(str,x)),np.rot90(m_prime)[::-1])
    dups = get_duplicates(m_prime_tmp)

    n_dup = len(dups)
    to_del = []
    if n_dup > 0:
        for idx,dup in enumerate(dups):
            to_del.extend(dup[1:])
            c_prime[dup[0]] = '_'.join(c_prime[dup])

    m_prime = np.delete(m_prime,to_del,axis=1)
    c_prime = np.delete(c_prime,to_del)

    return(m_prime, c_prime)

def remove_zero_entry_cols(m3,c3):
    '''
    remove any columns that have no 0 entries
    @param m3 the M matrix
    @param c3 column features corresponding to the M matrix
    '''
    mi = np.empty(0,dtype='int')
    ncol = len(m3[0])
    idxs_to_delete = [i for i in range(ncol) if not np.any(m3[:,i]==0)]
    m3 = np.delete(m3,idxs_to_delete,axis=1)
    c3 = np.delete(c3,idxs_to_delete)
    return(m3,c3)

def get_m_prime(m,c):
    '''
    Construct M' matrix
    Obtain binary encoding for each feature, sort by binary value, 
    then rotate and reverse the resulting matrix
    @param nshared nxm matrix of samples (cols) and 
    @param svrot is true, the matrix is of form rows = svs, cols = samples
    '''
    m_prime, c_prime = remove_zero_entry_cols(m,c)
    m_prime, c_prime = remove_duplicates(m_prime, c_prime)
    m_prime = np.rot90(m_prime)
    
# count binary score of columns
    binary_strings = []
    for col in m_prime:
        col = np.array([ci if ci>0 else 0  for ci in col])
        col_string = '0b'+''.join(map(str,col))
        binary_strings.append(int(col_string,2))
        
# sort by binary score
    order = np.argsort(binary_strings)[::-1]
    m_prime = m_prime[order] 
    m_prime = np.rot90(m_prime)[::-1] #rotate again
    c_order = (len(c_prime) - 1) - order #translate order of rotated matrix to order of columns
    c_prime = c_prime[c_order]
    
    return(m_prime,c_prime)

def get_k1_matrix(m_prime,features):
    '''
    Generate k1 from m' matrix
    Allows for checking of perfect phylogeny
    @param mp the m prime matrix
    @param nodes corresponding to the matrix
    '''
    ncol = len(m_prime[0])
    k = np.empty( [0,ncol], dtype='|S15' )

    for m in m_prime:
        row_feats = features[m!=0] #features in the row
        mrow = np.zeros(ncol,dtype='|S15')
        mrow.fill('0')

        for idx,feature in enumerate(row_feats):
            mrow[idx] = feature

        n_feat = len(row_feats)    
        if n_feat < ncol: 
            mrow[n_feat]='#'

        k = np.append(k,[mrow],axis=0)

    return(k)

def perfect_phylogeny_exists(k1,features):
    '''
    Determine whether perfect phylogeny exists from a k1 matrix    
    @param k1 the k1 matrix (output from get_k1_matrix)
    @param features the column features
    '''
    locations = []
    for feature in features:
        present_at = set([])
        for k_i in k1:
            [ present_at.add(loc_list) for loc_list in list(np.where(k_i==feature)[0]) ]
        locations.append(present_at)

    loc_test = np.array([len(loc_list)>1 for loc_list in locations])
    if np.any(loc_test):
        print('No phylogeny found!')
        return(False)
    else:    
        print('Success! Found phylogeny!\nK1 matrix:')
        print(k1)
        return(True)

def get_k(k,q,m_graph):
    '''
    return the first K vector with E[K] >= 1
    @param k holds the connection vector, initialise with set() 
    @param q chain of connected vertices
    @param m_graph graph object containing matrix connections
    '''
    if not q:
        return(k)
    else:
        q1 = q.pop()
        k.add(q1) 
        pairs = m_graph.get_pairs_containing(q1)
        pairs = set([node for pair in pairs for node in pair]) #flatten set of elements
        for node in pairs:
            if node not in k:
                q.add(node)
        return(get_k(k,q,m_graph))

def solve_incomplete_phylogeny(m,s,c):
    '''
    implement Pe'er incomplete phylogeny algorithm
    @param m M matrix 
    @param s samples/species (rows)
    @param c features (columns)
    '''
    m3, c3 = get_m_prime(m,c)
    
    m_graph = MGraph()
    m_graph.build_graph(m3,s,c3)
    m_pairs = m_graph.get_edge_pairs()
    
    q = set(m_pairs[0]) #pick the first elements as k
    k = get_k(set(),q,m_graph)
    t = [set(s)] #initialise a tree

    while len(m_pairs) > 1:
        while len(k) < 3: #
            for n in k:
                m_graph.delNode(n)
            m_pairs = m_graph.get_edge_pairs()
            if len(m_pairs) > 1:
                q = set(m_pairs[0]) 
                k = get_k(set(),q,m_graph)
            else: 
                break
        
        s_prime   = set(s).intersection(k)                
        s_indexes = np.array([np.where(s_i==s)[0][0] for s_i in s_prime])
        m_tmp     = m3.copy()[s_indexes]

        print('k:')
        print(k)
        print("S':")
        print(s_prime)
        
        u = []
        c_in_k  = set(c).intersection(k)
        for c_i in c_in_k:
            c_index = np.where(c_i==c3)[0]
            m_col = m_tmp[:,c_index:c_index+1]
            if np.all(m_col!=0):
                cm = np.where(c_i==c3)[0]
                u.append(c_i)
        
        if not u:
            break
        else:
            print('u:')
            print(u)
            t.append(s_prime)
            for n in u:
                m_graph.delNode(n)
            m_pairs = m_graph.get_edge_pairs()
            k = get_k(set(),set(m_pairs[0]),m_graph)
            print('Tree:')
            print(t)

    m_new = m3.copy()

    c_indexes = [np.where(c3==x)[0][0] for x in u]
    for s_set in t[1:]:
        s_prime = set(s).intersection(s_set)
        s_indexes = np.array([np.where(s_i==s)[0][0] for s_i in s_prime])
        m_tmp = m3.copy()[s_indexes]

        for c_i in c3:
            c_idx = np.where(c_i==c3)[0]
            m_col = m_tmp[:,c_index:c_index+1]
            if np.all(m_col!=0) and np.any(m_col==-1):
                cm = np.where(c_i==c3)[0]
                m_col[m_col==-1] = 1
                m_new[s_indexes,c_idx:c_idx+1] = m_col

    for i in range(len(c3)):
        tcol = m_new[:,i:i+1]
        tcol[tcol==-1] = 0
        m_new[:,i:i+1] = tcol
    
    print("New M':")
    m_new, c_prime = get_m_prime(m_new,c3)
    k1 = get_k1_matrix(m_new,c_prime)
    print(m_new)
    
    return(m_new, c_prime, k1, t)

def write_dot(k1,s,outname):
    '''    
    takes a k1 matrix and writes a dot source file of the node 
    connections and edges, then converts the dot to a postscript file. 
    @param k1 matrix output from perfect phylogeny checking
    @param s samples (rows)
    '''
    dot_lines = []

    for r_idx,k_row in enumerate(k1):
        for c_idx,k_i in enumerate(k_row):
            if k_i=='#': break
            
            k_next = k_row[c_idx+1]
            if c_idx==0:
                dot_lines.append('\troot [label=""];\n')
                dot_lines.append('\troot -> node_%s [label="%s"];\n' % (k_i,k_i))
            if k_next=='#':
                dot_lines.append('\tnode_%s [label=""];\n' % k_i)
                dot_lines.append('\tnode_%s -> %s;\n' % (k_i,s[r_idx]))
                break                
            dot_lines.append('\tnode_%s [label=""];\n\tnode_%s [label=""];\n' % (k_i,k_next))
            dot_lines.append('\tnode_%s -> node_%s [label="%s"];\n' % (k_i,k_next,k_next))
   
    dot_lines = np.unique(np.array(dot_lines))
    with open('%s.dot'%outname,'w') as fout:
        fout.write('digraph {\n')
        fout.write('\tgraph[size="7.75,10.25"]\n')
        for line in dot_lines:
            fout.write(line)
        fout.write('}\n')

    subprocess.call(['dot','-Tps','%s.dot'%outname],stdout=open('%s.ps'%outname,'w')) 
