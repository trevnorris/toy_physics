"""Stdlib-only extracted-helper dependency and unchanged-body checks."""
import ast
import operator
from pathlib import Path
from types import SimpleNamespace
import unittest

HERE=Path(__file__).parent
OLD=HERE/'S11c_d_numerical_radiating_balance.py'
NEW=HERE/'S11c_d_numerical_radiating_balance_v2.py'
def node(path,name):return next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
def helper():
    env={'operator':operator};exec(compile(ast.Module(body=[node(NEW,'complete_helper_namespace')],type_ignores=[]),str(NEW),'exec'),env);return env['complete_helper_namespace']
class ArrayStageReached(Exception):pass
class Tests(unittest.TestCase):
    def original(self):
        def array_stage(*a,**k):raise ArrayStageReached
        fake=SimpleNamespace(isfinite=lambda x:True,eye=array_stage,polynomial=SimpleNamespace(chebyshev=SimpleNamespace(chebder=array_stage)))
        env={'np':fake,'require':lambda condition,message:self.assertTrue(condition)}
        path=HERE/'S11c_d_finite_scattering.py'
        exec(compile(ast.Module(body=[node(path,'polynomial_basis')],type_ignores=[]),str(path),'exec'),env)
        return SimpleNamespace(polynomial_basis=env['polynomial_basis'])
    def test_original_missing_import_and_repaired_index_path(self):
        h=self.original()
        with self.assertRaises(NameError):h.polynomial_basis([],64,41,0)
        helper()(h)
        self.assertIs(h.polynomial_basis.__globals__['operator'],operator)
        with self.assertRaises(ArrayStageReached):h.polynomial_basis([],64,41,0)
    def test_original_noninteger_refusal_retained(self):
        h=self.original();helper()(h)
        with self.assertRaises(TypeError):h.polynomial_basis([],64,41,1.5)
    def test_numerical_run_unchanged_after_namespace_completion(self):
        old=node(OLD,'run');new=node(NEW,'run')
        self.assertIn('preservedBalanceStartupFiles',ast.unparse(new.body[0]));new.body.pop(0)
        new.body=[n for n in new.body if not(isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='complete_helper_namespace')]
        self.assertEqual(ast.dump(old),ast.dump(new))
    def test_only_standard_library_dependency_is_injected(self):
        h=self.original();before=dict(h.polynomial_basis.__globals__);helper()(h);after=h.polynomial_basis.__globals__
        self.assertEqual(set(after)-set(before),{'operator'})
        self.assertTrue(all(after[k] is v for k,v in before.items()))

if __name__=='__main__':unittest.main()
