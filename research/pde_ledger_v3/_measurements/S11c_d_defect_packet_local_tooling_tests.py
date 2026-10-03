#!/usr/bin/env python3
"""Stdlib metadata/source/synthetic checks only; no scientific imports/restoration."""
import ast,hashlib,json,runpy,sys,unittest
from fractions import Fraction as F
from pathlib import Path
M=Path(__file__).resolve().parent;P='S11c_d_defect_packet_local';R=M.parents[2]
load=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
class Tooling(unittest.TestCase):
 @classmethod
 def setUpClass(c):
  c.m=load(M/(P+'_inputs.json'));c.ns=runpy.run_path(str(M/(P+'_lib.py')),run_name='tooling_only');c.worker=ast.parse((M/(P+'.py')).read_text());c.lib=ast.parse((M/(P+'_lib.py')).read_text())
 def test_01_pins(self):
  for p,h in self.m['sourcePins'].items():self.assertEqual(sha(p),h)
 def test_02_complete_files(self):
  for r in self.m['savedInputs'].values():self.assertEqual((sha(r['path']),Path(r['path']).stat().st_size),(r['sha256'],r['bytes']))
 def test_03_all_selection(self):
  d=lambda n:load(self.m['savedInputs'][n]['path'])
  self.assertEqual(d('selected/local-cells.json')['selected'],[c for c in d('local/all-local-cells.json') if c['row']=='THETA_BALANCE' and c['field']=='e_W'])
 def test_04_batch_coverage(self):
  d=lambda n:load(self.m['savedInputs'][n]['path']);part=d('native/THETA-partition.json');children=[]
  for name in self.m['localBatches']:
   entries=d(name);inp=d(name.replace('-return','-input'))['children'];self.assertEqual([e['childIndex'] for e in entries],[e['childIndex'] for e in inp]);children+=entries
  self.assertEqual([e['childIndex'] for e in children],part['localChildIndices'])
 def test_05_cell_actual_summands(self):
  d=lambda n:load(self.m['savedInputs'][n]['path']);children=[e for n in self.m['localBatches'] for e in d(n)]
  for cell in d('selected/local-cells.json')['selected']:
   ev=[e for e in children if e['fieldColumn']==4 and e['xOrder']==cell['xOrder']]
   self.assertEqual(cell['sourceChildren'],[e['childIndex'] for e in ev]);self.assertEqual(cell['summands'],[next(v['value'] for v in e['mappedGrades'] if v['grade']==cell['grade']) for e in ev])
 def test_06_source_identity(self):
  d=lambda n:load(self.m['savedInputs'][n]['path']);row=d('native/THETA-row.json');u=d('units/source-text-joins.json')
  self.assertEqual(hashlib.sha256(row['fullConstructor'].encode()).hexdigest(),u['rowHashes']['THETA_BALANCE']);self.assertEqual(row['source'],u['provenance'])
 def test_07_no_rule_constructor_calls(self):
  called={n.func.attr for n in ast.walk(self.lib) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)}
  self.assertFalse(called & {'kronrod15','gauss_rule','H','middle','profile','point_plan'})
 def test_08_no_science_import_at_module(self):
  for tree in [self.worker,self.lib]:
   names=[n.names[0].name for n in tree.body if isinstance(n,ast.Import)];self.assertFalse(set(names)&{'sympy','mpmath','numpy'})
 def test_09_synthetic_profile_derivative_constant(self):self.assertEqual(self.ns['profile_derivative']([['7/3','-2']]),[['0','0'],['0','0']])
 def test_10_synthetic_profile_derivative_linear(self):self.assertEqual(self.ns['profile_derivative']([['0','0'],['10','20']]),[['1','2'],['0','0'],['-1','-2']])
 def test_11_synthetic_profile_derivative_quadratic(self):self.assertEqual(self.ns['profile_derivative']([['0','0'],['0','0'],['5','0']]),[['0','0'],['1','0'],['0','0'],['-1','0']])
 def test_12_nonrational_input_rejected(self):
  with self.assertRaises(ValueError):self.ns['profile_derivative']([['0','0'],['bad','0']])
 def test_13_zero_truth_refused(self):
  for v in [False,None,1,0,[],object()]:
   with self.assertRaises(ValueError):self.ns['require'](v,'synthetic')
 def test_14_controls_actual_request(self):
  n=next(n for n in ast.walk(self.lib) if isinstance(n,ast.FunctionDef) and n.name=='run');s=ast.unparse(n)
  self.assertIn("self.request(specs[cell], 'kappa', mutation)",s);self.assertIn("changed['value'] - baseline['value']",s);self.assertIn("10 * (baseline['envelope'] + changed['envelope'])",s)
 def test_15_no_coefficient_derivative_baseline(self):
  n=next(n for n in ast.walk(self.lib) if isinstance(n,ast.FunctionDef) and n.name=='integrand');s=ast.unparse(n)
  self.assertIn("if mutant == 'Leibniz' else None",s);self.assertIn('coef * Q + correction',s)
 def test_16_final_encoder_and_gate(self):
  s=(M/(P+'.py')).read_text();self.assertIn('J.encode(result)',s);self.assertIn('verify_invocation(args,g,sys.argv)',s);self.assertLess(s.index('containment.json'),s.index('import sympy as sp'))
 def test_17_fixed_local_scope(self):
  m=self.m;self.assertEqual(m['scope']['localActions'],2);self.assertTrue(m['scope']['noPressureAction']);self.assertIsNone(m['resources']['durationLimits']);self.assertEqual(m['resources']['memoryBytes'],4*1024**3)
 def test_18_no_runtime_gate_yet(self):self.assertFalse((M/(P+'_gate.json')).exists())
 def test_19_runtime_scripts_parse(self):
  for suffix in ['.py','_lib.py','_launch.py']:ast.parse((M/(P+suffix)).read_text())
 def test_20_review_path_and_stage(self):
  self.assertEqual(self.m['reviewRecordWillBe'],str(M/(P+'_build_review_record.json')));s=(M/(P+'_launch.py')).read_text();self.assertIn("'defect_packet_local'",s);self.assertIn('01a0e01b-ef84-7192-817f-584cda5d339b',s)
if __name__=='__main__':unittest.main()
