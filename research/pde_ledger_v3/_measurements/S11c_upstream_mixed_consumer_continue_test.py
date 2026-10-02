#!/usr/bin/env python3
"""Stdlib restoration/persistence/source-equivalence tests; no scientific import."""
import ast
import hashlib
import json
from pathlib import Path
import tempfile
import textwrap
import unittest
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
SOURCE=(HERE/'S11c_upstream_mixed_consumer_continue.py').read_text()
TREE=ast.parse(SOURCE)
BASE=(HERE/'S11c_upstream_mixed_consumer.py').read_text()

def definition(name):return next(n for n in TREE.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name)
def load(names,ns=None):
 ns={} if ns is None else ns
 for name in names:exec(compile(ast.get_source_segment(SOURCE,definition(name)),'<tooling>','exec'),ns)
 return ns

class Continuation(unittest.TestCase):
 def test_tail_compiles_and_original_science_assignments_are_unchanged(self):
  tail=load(['resumed_tail'],{'ast':ast})['resumed_tail'](BASE)
  new=ast.parse('def resume():\n'+tail).body[0]
  old=next(n for n in ast.parse(BASE).body if isinstance(n,ast.FunctionDef) and n.name=='science')
  first=next(i for i,n in enumerate(old.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='downstream' for t in n.targets))
  old.body=old.body[first:]
  def vals(node):
   result=[]
   for n in ast.walk(node):
    if isinstance(n,ast.Assign):result.append((';'.join(ast.unparse(t) for t in n.targets),ast.dump(n.value,include_attributes=False)))
   return result
  a,b=vals(old),vals(new);self.assertEqual(len(a),len(b))
  differences=[(x[0],y[0]) for x,y in zip(a,b) if x!=y]
  self.assertEqual(differences,[('downstream','downstream'),('consumer_controls','consumer_controls')])
  self.assertIn("source-consumer-control-%03d-input",tail)
  self.assertIn("pressure-consumer-control-%03d-input",tail)
  compile('def resume():\n'+tail,'<unfinished>','exec')
 def test_restore_uses_constructor_not_display(self):
  calls=[]
  ns=load(['load_saved_json'],{'Path':Path,'json':json,'Str':object(),
   'sp':SimpleNamespace(sympify=lambda expr,locals: calls.append(expr) or ('decoded',expr))})
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'saved.json';p.write_text(json.dumps({'item':{'text':'WRONG DISPLAY','srepr':'Integer(7)'},'other':['unchanged']}))
   raw,value=ns['load_saved_json'](p)
   self.assertEqual(calls,['Integer(7)']);self.assertEqual(value['item'],('decoded','Integer(7)'))
   self.assertEqual(raw['item']['text'],'WRONG DISPLAY')
 def test_persist_each_record_before_append(self):
  events=[]
  j=SimpleNamespace(emit=lambda name,value:events.append((name,value)))
  values=load(['persist_control_list'])['persist_control_list'](j,'control')
  values.append({'first':1});values.append({'second':2})
  self.assertEqual(events,[('control-000-return',{'first':1}),('control-001-return',{'second':2})])
  def fail(name,value):raise OSError('disk write refused')
  blocked=load(['persist_control_list'])['persist_control_list'](SimpleNamespace(emit=fail),'control')
  with self.assertRaises(OSError):blocked.append({'not_saved':3})
  self.assertEqual(blocked,[])
 def test_residual_evidence_precedes_fatal_guard(self):
  code=ast.get_source_segment(SOURCE,definition('resumed_operator_form'))
  self.assertLess(code.index("J.emit(name+'-input'"),code.index('jets=sorted'))
  self.assertLess(code.index("J.emit(name+'-raw-reconstruction'"),code.index('sp.cancel(sp.together(raw_remainder))'))
  self.assertLess(code.index("J.emit(name+'-exact-residual'"),code.index("require(remainder==0"))
  self.assertIn("require(not any(dependent_coefficients.values())",code)
  self.assertNotIn('evalf',code);self.assertNotIn('tolerance=',code)
 def test_unmodified_nested_helpers_and_no_producer_call(self):
  ns=load(['nested_fragment'],{'ast':ast,'textwrap':textwrap})
  for name in ['bind','wave_atoms','restrict']:
   fragment=ns['nested_fragment'](BASE,'science',name)
   compile(fragment,'<native>','exec')
  code=ast.get_source_segment(SOURCE,definition('science'))
  self.assertNotIn("['science'](",code)
  self.assertNotIn('build_records()',code)
  self.assertIn("originalScienceCalled=False",code)
 def test_same_original_selected_increment_and_resource_helpers(self):
  old={n.name:n for n in ast.parse(BASE).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
  for name in ['containment','selected_increment','exact_nonzero','Journal','save']:
   self.assertEqual(ast.dump(old[name],include_attributes=False),ast.dump(definition(name),include_attributes=False))
 def test_no_science_import_before_containment_or_deadlines(self):
  imports=[a.name for n in TREE.body if isinstance(n,(ast.Import,ast.ImportFrom)) for a in n.names]
  self.assertNotIn('sympy',imports)
  main=ast.get_source_segment(SOURCE,definition('main'))
  self.assertLess(main.index('containment()'),main.index('import sympy'))
  self.assertNotIn('timeout=',SOURCE)
  self.assertNotIn('signal.alarm',SOURCE)
 def test_selected_prior_files_exist_and_failure_is_preserved(self):
  p=HERE.parents[2]/'_scratch/s11c/s11c-mixed-consumer-20261001/diagnostic-01/complete'
  self.assertEqual(len(list(p.iterdir())),130)
  self.assertEqual(json.loads((p/'operation-index.json').read_text()),[])
  self.assertEqual(json.loads((p/'checks.json').read_text())['executionStatus'],'FAILED_PRESERVED')
  for n in ['end-to-end-controls-input.json','binding-context.json','source-controls.json','weak-directions.json']:
   self.assertIsInstance(json.loads((p/n).read_text()),dict)

if __name__=='__main__':unittest.main()
