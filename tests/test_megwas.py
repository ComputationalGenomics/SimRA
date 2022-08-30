from unittest import TestCase
import pandas as pd

from SimRA.statistical_genetics import _genotypes


HM_inf = pd.read_csv("tests/CEUASWMEX_fst_frq.txt",sep=u' ')

class TestModelBn(TestCase):

    def test_model_bn(self):
        result_pop_matrix, result_pop_idx, result_genetic_matrix =  _genotypes.model_bn(3, 100000, 1000, HM_inf[u'FRQ'].values, HM_inf[u'FST'].values, 3) # Uses default values for populations
        #self.assertEqual(result_pop_matrix, ) # Placeholder for real test implementation - need to input expected value
        #self.assertEqual(result_pop_idx, ) # Placeholder for real test implementation - need to input expected value
        #self.assertEqual(result_genetic_matrix,  ) # Placeholder for real test implementation - need to input expected value
        print(result_pop_matrix)
        print(result_pop_idx)
        print(result_genetic_matrix)