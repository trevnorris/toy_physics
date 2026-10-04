"""Tooling tests on synthetic units and constructor shapes, not native science."""
import ast,importlib.util,json,sys,tempfile,unittest
from pathlib import Path
from fractions import Fraction as F
M=Path(__file__).resolve().parent

def load(name):
 spec=importlib.util.spec_from_file_location(name,M/(name+'.py'));v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v);return v
U=load('S11c_d_defect_packet_native_units_lib');W=load('S11c_d_defect_packet_native_units')
S="Symbol('a')";T="Symbol('b')"
class UnitsTest(unittest.TestCase):
 def dim(self,text,registry=None):return U.Units(registry or {'a':[1,0,0],'b':[0,1,0]}).dimension(U.parse(text))
 def test_product(self):self.assertEqual(self.dim('Mul('+S+','+T+')'),(1,1,0))
 def test_inverse(self):self.assertEqual(self.dim('Pow('+S+',Integer(-1))'),(-1,0,0))
 def test_rational_power(self):self.assertEqual(self.dim('Pow('+S+',Rational(1,2))'),(F(1,2),0,0))
 def test_homogeneous_add(self):self.assertEqual(self.dim('Add('+S+',Mul(Integer(2),'+S+'))'),(1,0,0))
 def test_heterogeneous_add(self):
  with self.assertRaisesRegex(ValueError,'inhomogeneous'):self.dim('Add('+S+','+T+')')
 def test_exp(self):self.assertEqual(self.dim('exp(Mul('+S+',Pow('+S+',Integer(-1))))'),(0,0,0))
 def test_bad_exp(self):
  with self.assertRaisesRegex(ValueError,'dimensionful'):self.dim('exp('+S+')')
 def test_delta(self):self.assertEqual(self.dim('DiracDelta('+S+')'),(-1,0,0))
 def test_zero_unspecified(self):self.assertIsNone(self.dim('Integer(0)'))
 def test_zero_add(self):self.assertEqual(self.dim('Add(Integer(0),'+S+')'),(1,0,0))
 def test_zero_does_not_hide_unknown(self):
  with self.assertRaisesRegex(ValueError,'unknown'):self.dim("Mul(Integer(0),Symbol('unknown'))")
 def test_zero_pole(self):
  with self.assertRaisesRegex(ValueError,'zero denominator'):self.dim('Pow(Integer(0),Integer(-1))')
 def test_unknown(self):
  with self.assertRaisesRegex(ValueError,'unknown'):self.dim("Symbol('missing')")
 def test_no_float(self):
  for s in ['Float(0.5)','Rational(1.0,2)','Integer(True)']:
   with self.assertRaises(ValueError):self.dim(s)
 def test_no_executable_constructor(self):
  with self.assertRaises(ValueError):self.dim("__import__('os').system('false')")
 def test_no_symbolic_exponent(self):
  with self.assertRaises(ValueError):self.dim('Pow('+S+','+T+')')
 def test_symbol_assumptions(self):self.assertEqual(self.dim("Symbol('a',real=True)"),(1,0,0))
 def test_nonliteral_assumption(self):
  with self.assertRaises(ValueError):self.dim("Symbol('a',real=bool(1))")
 def test_unit_registry_rationals(self):self.assertEqual(U.unit_tuple([{'srepr':'Integer(-2)'},{'srepr':'Rational(1,2)'},{'srepr':'Integer(0)'}]),(-2,F(1,2),0))
 def test_nan_unresolved(self):
  with self.assertRaises(ValueError):U.unit_tuple([{'srepr':'nan'},0,0])
 def test_tag_selector(self):self.assertTrue(U.same(U.tagged(U.parse("Tuple(Tuple(Str('X'),Integer(2)))"),'X'),U.parse('Integer(2)')))
 def test_duplicate_tag(self):
  with self.assertRaisesRegex(ValueError,'unique'):U.tagged(U.parse("Tuple(Tuple(Str('X'),Integer(2)),Tuple(Str('X'),Integer(2)))"),'X')
 def test_restore_line(self):self.assertEqual(U.quoted_restore_line('  \'value\': _restore("Tuple(Integer(2))"),'),'Tuple(Integer(2))')
 def test_invalid_restore_line(self):
  for t in ["'wrong': _restore('Integer(2)'),","'value': dangerous('Integer(2)'),","'value': _restore(eval('x')),"]:
   with self.assertRaises(ValueError):U.quoted_restore_line(t)
 def test_persist_failure_trace(self):
  with tempfile.TemporaryDirectory() as d:
   j=W.Journal(Path(d))
   with self.assertRaisesRegex(ValueError,'inhomogeneous'):W.certify(j,U,'synthetic','Add('+S+','+T+')',{'a':[1,0,0],'b':[0,1,0]},[1,0,0])
   self.assertTrue((Path(d)/'synthetic-input.json').exists());self.assertTrue((Path(d)/'synthetic-complete-unit-walk.json').exists());self.assertEqual(j.active,'synthetic')
 def test_persist_success(self):
  with tempfile.TemporaryDirectory() as d:
   j=W.Journal(Path(d));W.certify(j,U,'synthetic',S,{'a':[1,0,0]},[1,0,0]);self.assertEqual(j.completed,['synthetic'])
 def test_no_scientific_imports(self):
  t=ast.parse((M/'S11c_d_defect_packet_native_units.py').read_text());names={a.name for n in ast.walk(t) if isinstance(n,(ast.Import,ast.ImportFrom)) for a in n.names};self.assertFalse(names&{'sympy','mpmath','scipy','numpy'})
 def test_guard_before_work(self):
  s=(M/'S11c_d_defect_packet_native_units.py').read_text();self.assertLess(s.index('g=verify_gate('),s.index("enforced=ns['containment']()"));self.assertLess(s.index("enforced=ns['containment']()"),s.index('result=run(m,J,U)'))
if __name__=='__main__':unittest.main()
