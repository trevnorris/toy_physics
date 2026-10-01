#!/usr/bin/env python3
"""Stdlib-only parser, schema, persistence and routing regression tests. No CAS."""
import ast
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
SOURCE=(HERE/'S11c_upstream_mixed_consumer.py').read_text()
TREE=ast.parse(SOURCE)
SN={ '__file__':str(HERE/'S11c_upstream_mixed_consumer_source.py'),'__name__':'metadata_test'}
exec(compile((HERE/'S11c_upstream_mixed_consumer_source.py').read_text(),'<metadata>','exec'),SN)
DATA=json.loads((HERE/'S11c_upstream_mixed_consumer_native.json').read_text())

def definition(name):
 return next(n for n in TREE.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name)

class Tooling(unittest.TestCase):
 def test_tuple_navigation_preserves_nested_strings(self):
  s="Tuple(Str('x,y'), Tuple(Integer(-1), Symbol('a', real=True)))"
  self.assertEqual(SN['tuple_arguments'](s),["Str('x,y')","Tuple(Integer(-1), Symbol('a', real=True))"])
  self.assertEqual(SN['literal_key']("Tuple(Str('LAB_HELD'), Integer(-1))"),('LAB_HELD',-1))
  with self.assertRaises(ValueError):SN['literal_key']("Symbol('not_a_label')")
 def test_pressure_census_rejects_unsupported_derivative(self):
  row="Tuple(Tuple(Str('EXPANDED'), Add(Symbol('delta_p_plus'), Symbol('unrelated'))))"
  r=SN['row_census']('THETA_BALANCE',row)[0]
  self.assertEqual(r['totalChildren'],2)
  self.assertEqual(len(r['selectedPressureJetChildren']),1)
  with self.assertRaises(ValueError):SN['row_census']('THETA_BALANCE',row.replace('delta_p_plus','d_1_delta_p_plus'))
 def test_actual_schema_routes_and_symbol_names(self):
  responses=DATA['faceResponseSources']['cases']
  self.assertEqual([r['case'][1] for r in responses],[1,-1])
  for r in responses:
   names={ast.literal_eval(n.args[0]) for n in ast.walk(ast.parse(r['valueConstructorText'],mode='eval'))
     if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='Symbol'}
   label='plus' if r['case'][1]==1 else 'minus'
   self.assertIn('s11cc1_mu_theta_lab_held_'+label,names)
   self.assertIn('s11cc1_V_lab_held_'+label,names)
  names=[]
  for r in DATA['slabConsumers']['rows']:
   a=r['address'];names.append('U'+str(a[-1]) if a[-2]=='EXPANDED' else a[-2])
  self.assertEqual(names,['U0','U1','U2','THETA_BALANCE','E_W_BALANCE'])
 def test_selected_child_hashes_and_full_index(self):
  for r in DATA['slabConsumers']['rows']:
   indices=r['completeTopLevelChildIndex']
   self.assertEqual([x['childIndex'] for x in indices],list(range(r['totalChildren'])))
   selected=r['selectedPressureJetChildren']
   self.assertEqual([x['childIndex'] for x in indices if x['pressureSymbolOccurrences']],
                    [x['childIndex'] for x in selected])
   for x in selected:self.assertEqual(hashlib.sha256(x['constructorText'].encode()).hexdigest(),x['constructorSha256'])
 def test_constructor_strings_all_parse(self):
  def walk(x):
   if isinstance(x,dict):
    for k,v in x.items():
     if k.endswith('ConstructorText') or k=='constructorText':ast.parse(v,mode='eval')
     else:walk(v)
   elif isinstance(x,list):
    for v in x:walk(v)
  walk(DATA)
 def test_native_helper_namespace(self):
  c2=(HERE.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
  nodes={n.name:n for n in ast.parse(c2).body if isinstance(n,ast.FunctionDef)}
  for name in ['field','wave_jet']:
   compile(ast.get_source_segment(c2,nodes[name]),'<native>','exec')
  science=ast.get_source_segment(SOURCE,definition('science'))
  self.assertIn('ns = dict(sp=sp,re=re,X=X,TIME=time_symbol,NEW_DIMENSIONS={})',science)
  self.assertIn("for name in ('field','wave_jet')",science)
 def test_no_science_import_before_containment_or_timer(self):
  names=[]
  for n in TREE.body:
   if isinstance(n,(ast.Import,ast.ImportFrom)):
    names.extend(a.name for a in n.names)
  self.assertNotIn('sympy',names)
  main=ast.get_source_segment(SOURCE,definition('main'))
  self.assertLess(main.index('containment()'),main.index('import sympy'))
  calls=[ast.unparse(n.func) for n in ast.walk(TREE) if isinstance(n,ast.Call)]
  self.assertFalse(any(x.endswith(('.alarm','.setitimer')) for x in calls))
  self.assertNotIn('timeout=',SOURCE)
 def test_exclusive_json_persistence(self):
  ns={'Path':Path,'json':json,'os':__import__('os')}
  exec(compile(ast.get_source_segment(SOURCE,definition('save')),'<save>','exec'),ns)
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'result.json';ns['save'](p,{'first':1})
   with self.assertRaises(FileExistsError):ns['save'](p,{'replacement':2})
   self.assertEqual(json.loads(p.read_text()),{'first':1})
 def test_failed_stage_keeps_input_and_active_name(self):
  from types import SimpleNamespace
  ns={'Path':Path,'json':json,'os':__import__('os'),'hashlib':hashlib,
      'sp':SimpleNamespace(Basic=type('Basic',(),{}),MatrixBase=type('MatrixBase',(),{}))}
  for name in ['require','sha','save','Journal']:
   exec(compile(ast.get_source_segment(SOURCE,definition(name)),'<tooling>','exec'),ns)
  with tempfile.TemporaryDirectory() as d:
   j=ns['Journal'](Path(d))
   def fail():raise RuntimeError('deliberate test')
   with self.assertRaises(RuntimeError):j.stage('stage',{'operand':7},fail)
   self.assertEqual(json.loads((Path(d)/'stage-input.json').read_text()),{'operand':7})
   self.assertEqual(j.active,'stage-input');self.assertEqual(j.completed,[])
 def test_scope_keeps_selected_offshell_boundary(self):
  scope=(HERE/'S11c_upstream_mixed_consumer_scope.md').read_text()
  for phrase in ['not an','on-shell','no second middle integration','No full on-shell']:
   self.assertIn(phrase,scope)
  self.assertIn("'UNRESOLVED: nonzero/unknown source",SOURCE)

if __name__=='__main__':unittest.main()
