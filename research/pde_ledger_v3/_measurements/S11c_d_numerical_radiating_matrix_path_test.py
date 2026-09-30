"""Stdlib stand-in test of actual inherited matrix/row traversal dispatch.

No scientific imports, quadrature evaluation or saved-object restoration.
This verifies source wiring and order selection, not numerical action accuracy.
"""
import ast
from pathlib import Path
from types import SimpleNamespace
import unittest

ROOT=Path(__file__).resolve().parents[3]
ENGINE=ROOT/'research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py'
INTEGRATION=Path(__file__).with_name('S11c_d_numerical_radiating_integration_v3.py')
FINITE=Path(__file__).with_name('S11c_d_finite_scattering.py')
def method(path,cls,name):
    tree=ast.parse(path.read_text())
    owner=next(n for n in ast.walk(tree) if isinstance(n,ast.ClassDef) and n.name==cls)
    return next(n for n in owner.body if isinstance(n,ast.FunctionDef) and n.name==name)
class FakeArray(list):
    observed_weights=[]
    def __matmul__(self,other):
        self.observed_weights.extend(self);return 0
class Accumulator:
    def __iadd__(self,x):return self
class FakeNP:
    @staticmethod
    def asarray(x):return FakeArray(x)
    @staticmethod
    def zeros(*a,**k):return Accumulator()
    @staticmethod
    def isfinite(x):return SimpleNamespace(all=lambda:True)
def require(value,message):
    if not value:raise AssertionError(message)
def extract(path,cls,name):
    env={'np':FakeNP,'require':require}
    exec(compile(ast.Module(body=[method(path,cls,name)],type_ignores=[]),str(path),'exec'),env)
    return env[name]
class Width:
    def subs(self,*args):return .1
class Parent:
    batches=extract(ENGINE,'ThreeMomentum','batches')
class Probe(Parent):
    row_integral=extract(INTEGRATION,'Radiating','row_integral')
    def __init__(self,setting,pairs):
        self.setting=setting;self.binding={'pairs':pairs,'abel':{'width':Width()}}
        self.r=SimpleNamespace(regulator='reg');self.batch_nodes=3;self.calls=[];self.visits=[]
    def rule(self,lower,upper,order,centers=(),width=None):
        self.calls.append((lower,upper,order,tuple(centers),width))
        shift=sum(centers)/16
        return (-.5+shift,.5+shift),(.25,.75),(lower,upper)
    def row_batch(self,row,variables,points,positions):
        self.visits.extend(tuple(x) for x in points);return 0
class Tests(unittest.TestCase):
    def test_inheritance_is_actual_selected_parent(self):
        owner=next(n for n in ast.parse(FINITE.read_text()).body if isinstance(n,ast.ClassDef) and n.name=='BasisMomentum')
        self.assertEqual(ast.unparse(owner.bases[0]),'engine.BoundedSourceFourierQuadrature.ThreeMomentum')
        self.assertFalse(any(isinstance(n,ast.FunctionDef) and n.name=='batches' for n in owner.body))
    def test_same_rule_calls_points_weights_for_native_dimensions(self):
        for outer,inner in ((48,(8,8)),(64,(12,12)),(24,(8,8)),(32,(12,12))):
            setting={'kind':'split','outerOrder':outer,'innerOrders':inner,'momentumBound':4.,'regulator':.1}
            for variables in (('a',),('a','b'),('a','b','c')):
                pairs=tuple(zip(variables,variables[1:]));p=Probe(setting,pairs)
                points=[];weights=[]
                for x,w in p.batches(variables,setting,pairs,.1):points.extend(tuple(v) for v in x);weights.extend(w)
                calls=list(p.calls);p.calls.clear();FakeArray.observed_weights=[]
                result=p.row_integral({'limits':[(v,None,None) for v in variables]},[0])
                self.assertEqual(p.calls,calls);self.assertEqual(p.visits,points);self.assertEqual(FakeArray.observed_weights,weights);self.assertEqual(result['nodes'],len(points))
                self.assertAlmostEqual(sum(weights),1.)
                for _,_,order,centers,_ in calls:self.assertEqual(order,inner[0] if centers else outer)
    def test_matrix_path_dispatches_through_worker(self):
        path=Path(__file__).with_name('S11c_d_numerical_radiating_finite.py');tree=ast.parse(path.read_text())
        fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='matrix_group')
        self.assertTrue(any(isinstance(n,ast.Call) and ast.unparse(n.func)=='worker.batches' for n in ast.walk(fn)))

if __name__=='__main__':unittest.main()
