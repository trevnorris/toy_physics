"""Tooling tests on synthetic units and constructor shapes, not native science."""
import ast,hashlib,importlib.util,json,sys,tempfile,unittest
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
   self.assertTrue((Path(d)/'synthetic-input.json').exists());self.assertTrue((Path(d)/'synthetic-unit-walk.json').exists());self.assertEqual(j.active,'synthetic')
   decision=json.loads((Path(d)/'synthetic-decision-operands.json').read_text());self.assertEqual(decision['status'],'UNIT_WALK_REFUSED');self.assertTrue(decision['partialWalk'])
   self.assertFalse((Path(d)/'synthetic-return.json').exists())
 def test_persist_success(self):
  with tempfile.TemporaryDirectory() as d:
   j=W.Journal(Path(d));W.certify(j,U,'synthetic',S,{'a':[1,0,0]},[1,0,0]);self.assertEqual(j.completed,['synthetic'])
 def test_no_scientific_imports(self):
  t=ast.parse((M/'S11c_d_defect_packet_native_units.py').read_text());names={a.name for n in ast.walk(t) if isinstance(n,(ast.Import,ast.ImportFrom)) for a in n.names};self.assertFalse(names&{'sympy','mpmath','scipy','numpy'})
 def test_guard_before_work(self):
  s=(M/'S11c_d_defect_packet_native_units.py').read_text();self.assertLess(s.index('g=verify_gate('),s.index("enforced=ns['containment']()"));self.assertLess(s.index("enforced=ns['containment']()"),s.index('result=run(m,J,U)'))

class SourceCaseTest(unittest.TestCase):
 def table(self):return "Tuple(Tuple(Tuple(Str('HELD'),Integer(1)),Tuple(Tuple(Str('VALUE'),Integer(7)))),Tuple(Tuple(Str('HELD'),Integer(-1)),Tuple(Tuple(Str('VALUE'),Integer(7)))))"
 def saved(self):
  payload="Tuple(Tuple(Str('VALUE'),Integer(7)))"
  return {'case':['HELD',1],'caseIndex':0,'caseConstructorText':payload,'caseConstructorSha256':hashlib.sha256(payload.encode()).hexdigest(),'valueConstructorText':'Integer(7)'}
 def test_original_labels_and_index(self):
  i,p,c=U.select_case(U.parse(self.table()),['HELD',-1]);self.assertEqual(i,1);self.assertEqual(c,[['HELD',1],['HELD',-1]]);self.assertTrue(U.same(U.tagged(p,'VALUE'),U.parse('Integer(7)')))
 def test_missing_labels_even_with_same_payload(self):
  with self.assertRaisesRegex(ValueError,'unique original labeled case'):U.select_case(U.parse(self.table()),['ADVECTED',1])
 def test_duplicate_labels_refused(self):
  with self.assertRaisesRegex(ValueError,'unique original labeled case'):U.select_case(U.parse(self.table().replace('Integer(-1)','Integer(1)')),['HELD',1])
 def test_boolean_label_refused(self):
  with self.assertRaisesRegex(ValueError,'literal expected'):U.select_case(U.parse(self.table()),['HELD',True])
 def test_noninteger_label_refused(self):
  with self.assertRaisesRegex(ValueError,'integer native'):U.select_case(U.parse(self.table().replace('Integer(-1)','Rational(1,2)')),['HELD',1])
 def test_executable_label_refused(self):
  with self.assertRaises(ValueError):U.select_case(U.parse(self.table().replace("Str('HELD')","Str(str('HELD'))")),['HELD',1])
 def test_nested_table(self):
  node=U.tagged(U.parse("Tuple(Tuple(Str('CASES'),"+self.table()+"))"),'CASES');self.assertEqual(U.select_case(node,['HELD',1])[0],0)
 def test_join_full_case_evidence(self):
  with tempfile.TemporaryDirectory() as d:
   j=W.Journal(Path(d));v=W.join_native_case(j,U,'join',self.table(),self.saved());self.assertTrue(U.same(v,U.parse('Integer(7)')));self.assertEqual(j.completed,['join'])
   result=json.loads((Path(d)/'join-return.json').read_text());self.assertEqual(result['originalLabels'],[['HELD',1],['HELD',-1]])
 def test_wrong_saved_index_saved_before_refusal(self):
  with tempfile.TemporaryDirectory() as d:
   j=W.Journal(Path(d));case=self.saved();case['caseIndex']=1
   with self.assertRaisesRegex(ValueError,'actual native case'):W.join_native_case(j,U,'join',self.table(),case)
   self.assertTrue((Path(d)/'join-decision-operands.json').exists());self.assertFalse((Path(d)/'join-return.json').exists())
 def test_wrong_saved_payload_refused(self):
  with tempfile.TemporaryDirectory() as d:
   j=W.Journal(Path(d));case=self.saved();case['caseConstructorText']=case['caseConstructorText'].replace('7','8');case['caseConstructorSha256']=hashlib.sha256(case['caseConstructorText'].encode()).hexdigest()
   with self.assertRaisesRegex(ValueError,'actual native case'):W.join_native_case(j,U,'join',self.table(),case)
 def test_refused_selector_saves_input_and_decision(self):
  with tempfile.TemporaryDirectory() as d:
   j=W.Journal(Path(d));case=self.saved();case['case']=['OTHER',1]
   with self.assertRaisesRegex(ValueError,'unique original'):W.join_native_case(j,U,'join',self.table(),case)
   self.assertEqual(json.loads((Path(d)/'join-decision-operands.json').read_text())['status'],'SOURCE_CASE_REFUSED');self.assertEqual(j.active,'join')
 def test_keyed_dtn_literal(self):
  value,node=U.export_restore_literal(ast.parse("_LEDGER={'dtn_kernel':{'value':_restore('Integer(7)')}}"),'dtn_kernel');self.assertEqual(value,'Integer(7)');self.assertEqual(U.kind(node),'_restore')
 def test_wrong_key_same_literal_refused(self):
  with self.assertRaisesRegex(ValueError,'unique original dictionary key'):U.export_restore_literal(ast.parse("_LEDGER={'other':{'value':_restore('Integer(7)')}}"),'dtn_kernel')
 def test_duplicate_export_key_refused(self):
  with self.assertRaisesRegex(ValueError,'unique original dictionary key'):U.export_restore_literal(ast.parse("_LEDGER={'dtn_kernel':{'value':_restore('Integer(7)')},'dtn_kernel':{'value':_restore('Integer(7)')}}"),'dtn_kernel')
 def test_unquoted_restore_refused(self):
  with self.assertRaisesRegex(ValueError,'original keyed'):U.export_restore_literal(ast.parse("_LEDGER={'dtn_kernel':{'value':_restore(payload)}}"),'dtn_kernel')
if __name__=='__main__':unittest.main()
