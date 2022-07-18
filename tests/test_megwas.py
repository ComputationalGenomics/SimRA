from unittest import TestCase
from ..SimRA.statistical_genetics import _genotypes

print(_genotypes.model_bn)


class Test_model_bn(TestCase):
    
    def test_model_bn(self):
        result_pop_matrix, result_pop_idx = _genotypes.model_bn(3, 1000000, 1000, frq, fst, 3) # frq and fst are placeholders
        self.assertEqual(result_pop_matrix, )
        self.assertEqual(result_pop_idx, )
