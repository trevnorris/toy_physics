"""Exact rational identity fallback; synthetic standard-library values only."""
from fractions import Fraction
from types import SimpleNamespace
import unittest
import S11c_d_numerical_radiating_uniform_refinement_v2 as worker
from S11c_d_numerical_radiating_uniform_refinement_test import Tests as StorageTests

class MatrixStandIn:
    shape=(1,2)
    def __init__(self,values):self.values=values
    def __iter__(self):return iter(self.values)
    def __eq__(self,other):return False

class Tests(unittest.TestCase):
    def setUp(self):worker.sp=SimpleNamespace(cancel=lambda x:x)
    def test_exact_zero_residual_accepts_distinct_representations(self):
        a=MatrixStandIn([Fraction(1,2)+Fraction(1,3),Fraction(2,7)])
        b=MatrixStandIn([Fraction(5,6),Fraction(4,14)])
        result=worker.reference_identity({'matrix':a,'physicalScale':Fraction(20,2)},b,Fraction(10))
        self.assertFalse(result['rawMatrixEqual']);self.assertEqual(result['entryResiduals'],[0,0]);self.assertEqual(result['scaleResidual'],0)
        self.assertIs(result['savedMatrix'],a);self.assertIs(result['sourceMatrix'],b)
    def test_nonzero_residual_is_not_tolerated(self):
        a=MatrixStandIn([Fraction(1,2),Fraction(2,7)])
        b=MatrixStandIn([Fraction(1,2)+Fraction(1,10**30),Fraction(2,7)])
        result=worker.reference_identity({'matrix':a,'physicalScale':10},b,11)
        self.assertNotEqual(result['entryResiduals'][0],0);self.assertEqual(result['scaleResidual'],-1)
        self.assertFalse(all(x==0 for x in result['entryResiduals']) and result['scaleResidual']==0)

if __name__=='__main__':unittest.main()
