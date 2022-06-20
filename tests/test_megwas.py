from unittest import TestCase
from SimRA.statistical_genetics import _megwas

def test_is_equal():
    assert (_megwas.is_equal(2,2))
    assert (not _megwas.is_equal(2,1))
