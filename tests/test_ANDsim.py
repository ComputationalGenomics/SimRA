import unittest
from SimRA.statistical_genetics import ANDsim

class TestSim(unittest.TestCase):

    def test_simulate(self):
        Pn = 5
        Vn = 12
        N = 500
        v = 0.05
        k = [2, 3]
        seed = 0
        prop1 = 0.5
        phenodict, varmat = ANDsim.simulate(v, k, Pn, Vn, N, prop1, seed)
        self.assertEqual(any(ANDsim.simulate(v, k, Pn, Vn, N, prop1, seed)), any((phenodict, varmat)))


if __name__ == '__main__':
    unittest.main()
