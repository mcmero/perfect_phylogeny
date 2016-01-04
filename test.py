'''
Author: Marek Cmero
Phylogeny algorithm unit tests
'''
import unittest
import numpy as np
from perfect_phylogeny import phylogeny as phyl
from unittest import TestCase

# test data 1 - perfect phyl with incomplete fields
m_file1 = 'test_matrix.txt'
m_txt  = np.genfromtxt(m_file1,delimiter='\t',dtype=None,invalid_raise=False)

s1 = np.array(m_txt[1:,0],dtype='S100') # row names/samples
c1 = np.array(m_txt[0,1:],dtype='S100') # column names/features
m1 = np.array(m_txt[1:,1:],dtype=int)   # matrix values

# test data 2 - no perfect phylogeny, no incomplete fields
m_file2 = 'test_matrix_no_phyl.txt'
m_txt   = np.genfromtxt(m_file2,delimiter='\t',dtype=None,invalid_raise=False)

s2 = np.array(m_txt[1:,0],dtype='S100') # row names/samples
c2 = np.array(m_txt[0,1:],dtype='S100') # column names/features
m2 = np.array(m_txt[1:,1:],dtype=int)   # matrix values

# expected output
tree_exp = [set(['s3', 's2', 's1', 's5', 's4']), set(['s3', 's1', 's5']), set(['s3', 's1']), set(['s2', 's4'])]

m_exp  = np.array([[1, 1, 0, 0],
                   [0, 0, 1, 0],
                   [1, 1, 0, 0],
                   [0, 0, 1, 1],
                   [1, 0, 0, 0]],dtype='int')

k1_exp = np.array([['c2', 'c1',  '#', '0'],
                   ['c3',  '#',  '0', '0'],
                   ['c2', 'c1',  '#', '0'],
                   ['c3', 'c4',  '#', '0'],
                   ['c2',  '#',  '0', '0']],dtype='|S15') 

k1_nophyl = np.array([['c2', 'c1',  '#',  '0', '0'],
                      ['c5', 'c3',  '#',  '0', '0'],
                      ['c2', 'c1',  'c5', '#', '0'],
                      ['c3', 'c4',  '#',  '0', '0'],
                      ['c2', 'c5',  '#',  '0', '0']],dtype='|S15') 

class test_build_phyl(unittest.TestCase):
    
    def test_solve_inc(self):
        m_eval, c_eval, k1_eval, tree_eval = phyl.solve_incomplete_phylogeny(m1,s1,c1)
        k1 = phyl.get_k1_matrix(m_eval,c_eval)
        self.assertTrue(np.all(m_exp==m_eval))
        self.assertTrue(np.all(k1_exp==k1_eval))
        self.assertTrue(np.all(tree_exp==tree_eval))
        self.assertTrue(phyl.perfect_phylogeny_exists(k1,c_eval))

    def test_no_phyl(self):
        m_eval, c_eval = phyl.get_m_prime(m2,c2)
        k1 = phyl.get_k1_matrix(m_eval,c_eval)
        self.assertTrue(np.all(k1_nophyl==k1))
        self.assertFalse(phyl.perfect_phylogeny_exists(k1,c_eval))

if __name__ == '__main__':
    unittest.main()
