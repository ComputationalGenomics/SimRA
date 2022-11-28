from unittest import TestCase
import pandas as pd

from SimRA.statistical_genetics import _genotypes


HM_inf = pd.read_csv("tests/CEUASWMEX_fst_frq.txt",sep=u' ')
HGDP_PCs = pd.read_csv(u'tests/test_data/pruned_HGDP_topPops_singVecs.txt',sep=u' ',header=None)
HGDP_subpops = pd.read_csv(u'tests/test_data/pruned_HGDP_topPops_singVecs.txt',sep=u' ',header=None)
TGP_PCs = pd.read_csv(u'subpops_pruned_HGDP.txt',sep=u' ',header=None)
TGP_subpops = pd.read_csv(u'tests/test_data/subpops_pruned_TGP.txt',sep=u' ',header=None)

class TestModels(TestCase):

    def test_model_bn(self):
        result_pop_matrix, result_popidx, result_genetic_matrix =  _genotypes.model_bn(3, 100000, 1000, HM_inf[u'FRQ'].values, HM_inf[u'FST'].values, 3) # Uses default values for populations
        #self.assertEqual(result_pop_matrix, ) # Placeholder for real test implementation - need to input expected value
        #self.assertEqual(result_popidx, ) # Placeholder for real test implementation - need to input expected value
        #self.assertEqual(result_genetic_matrix,  ) # Placeholder for real test implementation - need to input expected value
        print(result_pop_matrix)
        print(result_popidx)
        print(result_genetic_matrix)

    def test_model_psd(self):
        result_pop_matrix, result_popidx, result_genetic_matrix = _genotypes.model_psd(3, 100000, 1000, HM_inf[u'FRQ'].values, HM_inf[u'FST'].values, 3)
        #self.assertEqual(result_pop_matrix, ) # Placeholder for real test implementation - need to input expected value
        #self.assertEqual(result_popidx, ) # Placeholder for real test implementation - need to input expected value
        #self.assertEqual(result_genetic_matrix,  ) # Placeholder for real test implementation - need to input expected value
        print(result_pop_matrix)
        print(result_popidx)
        print(result_genetic_matrix)

    def test_model_hgdp(self): 
        result_pop_matrix, result_popidx, result_genetic_matrix = _genotypes.model_hgdp(0, 10, 10000, 305, HGDP_PCs.values, HGDP_subpops.values, 3)
        #self.assertEqual(result_pop_matrix, ) # Placeholder for real test implementation - need to input expected value
        #self.assertEqual(result_popidx, ) # Placeholder for real test implementation - need to input expected value
        #self.assertEqual(result_genetic_matrix,  ) # Placeholder for real test implementation - need to input expected value
        print(result_pop_matrix)
        print(result_popidx)
        print(result_genetic_matrix)

    def test_model_tgp(self): 
        result_pop_matrix, result_popidx, result_genetic_matrix = _genotypes.model_tgp(0, 10, 10000, 1056, TGP_PCs.values, TGP_subpops.values, 3)
        #self.assertEqual(result_pop_matrix, ) # Placeholder for real test implementation - need to input expected value
        #self.assertEqual(result_popidx, ) # Placeholder for real test implementation - need to input expected value
        #self.assertEqual(result_genetic_matrix,  ) # Placeholder for real test implementation - need to input expected value
        print(result_pop_matrix)
        print(result_popidx)
        print(result_genetic_matrix)