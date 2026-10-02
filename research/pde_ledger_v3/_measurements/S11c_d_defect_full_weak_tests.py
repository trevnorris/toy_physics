#!/usr/bin/env python3
"""Stdlib tooling, syntax and source-metadata tests; no scientific restoration."""
import ast
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements';PREFIX='S11c_d_defect_full_weak'
WORKER=M/(PREFIX+'.py');SOURCE=WORKER.read_text();TREE=ast.parse(SOURCE)
NAMES={'require','sha','validate_partition','numeric_extension','chunks','checkpoint_batch','EvidenceLog','source_rules','verify_helper_paths','verify_invocation','pressure_child_coverage','inherited_slot_join','factor_address_join'}
NS={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'ROOT':ROOT,'__file__':str(WORKER),'SLOTS':('delta_p_plus','delta_p_minus','d_w_delta_p_plus','d_w_delta_p_minus')}
exec(compile(ast.Module([n for n in TREE.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in NAMES],type_ignores=[]),'inert-stdlib-definitions','exec'),NS)
read=lambda p:json.loads(Path(p).read_text())
MAN=read(M/(PREFIX+'_inputs.json'))
def saved(alias):return read(MAN['savedInputs'][alias]['path'])
def toy_partition():
 text="Add(Mul(Symbol('epsilon_shape'), Symbol('u_1')), Mul(Symbol('delta_p_plus'), Integer(2)))";args=ast.parse(text,mode='eval').body.args
 children=[]
 for i,node in enumerate(args):
  s=ast.get_source_segment(text,node);children.append({'childIndex':i,'constructorText':s,'sha256':hashlib.sha256(s.encode()).hexdigest(),'hits':[] if i==0 else [{'constructor':'Symbol','name':'delta_p_plus'}],'completeNameCoverage':True})
 return {'fullConstructor':text,'children':children,'localChildIndices':[0],'pressureChildIndices':[1],'unknownPressureAtoms':[]}
class Tests(unittest.TestCase):
 def test_partition_complete(self):self.assertEqual(NS['validate_partition'](toy_partition()),{'children':2,'local':1,'pressure':1,'exactSourcePartition':True})
 def test_partition_missing(self):
  p=toy_partition();p['localChildIndices']=[]
  with self.assertRaises(ValueError):NS['validate_partition'](p)
 def test_partition_overlap(self):
  p=toy_partition();p['pressureChildIndices']=[0,1]
  with self.assertRaises(ValueError):NS['validate_partition'](p)
 def test_partition_corrupt_source(self):
  p=toy_partition();p['children'][0]['constructorText']='Integer(0)'
  with self.assertRaises(ValueError):NS['validate_partition'](p)
 def test_partition_unknown_pressure(self):
  p=toy_partition();p['unknownPressureAtoms']=['delta_p_unknown']
  with self.assertRaises(ValueError):NS['validate_partition'](p)
 def test_actual_native_source_metadata(self):
  counts=[NS['validate_partition'](saved('inventory/'+r+'-full-native-partition.json')) for r in ['U0','U1','U2','THETA_BALANCE','E_W_BALANCE']]
  self.assertEqual(sum(c['local'] for c in counts),2908);self.assertEqual(sum(c['pressure'] for c in counts),12)
 def test_native_line_byte_convention(self):
  src=saved('inventory/U0-full-native-partition.json')['source']
  with Path(src['source']).open('rb') as f:line=next(v for i,v in enumerate(f,1) if i==src['valueLine'])
  self.assertEqual(len(line),src['sourceLineBytes']);self.assertEqual(hashlib.sha256(line).hexdigest(),src['sourceLineSha256'])
 def test_numeric_extension(self):
  n=NS['numeric_extension']({'omega':'1','eta_bg':'1/100','c_s0':'10','rho_br':'1','gamma':'3/101'},{},('eta_bg','sigma_W','epsilon_shape'))
  self.assertEqual(n,{'omega':'3','rho_br':'1','gamma':'3/101'})
 def test_actual_native_rule_ast(self):
  n=NS['source_rules'](Path(MAN['c2Source']).read_text(),Path(MAN['compositionSource']).read_text())
  self.assertIn('u_1',n['waveNames']);self.assertIn('u_1',n['dimensions'])
 def test_batch_preserves_completed_prefix(self):
  records={}
  class J:
   def emit(self,k,v):records[k]=copy.deepcopy(v)
   def stage(self,k,v,f):self.emit(k+'-input',v);return f()
  calls=[]
  def derive(v):
   calls.append(v)
   if v==3:raise ValueError('intended failure')
   return v*2
  with self.assertRaises(ValueError):NS['checkpoint_batch'](J(),'batch',[1,2,3,4],derive)
  self.assertEqual(calls,[1,2,3]);self.assertEqual(records['batch-partial-returns'],[2,4])
 def test_evidence_durable_chain(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'log';log=NS['EvidenceLog'](p,lambda v:v);log.append('input',{'x':3});log.append('return',{'result':9})
   data=[json.loads(v) for v in p.read_text().splitlines()]
   self.assertEqual(data[1]['payload']['previousSha256'],data[0]['sha256'])
   for i,r in enumerate(data):
    self.assertEqual(r['payload']['sequence'],i);self.assertEqual(r['sha256'],hashlib.sha256(json.dumps(r['payload'],sort_keys=True,allow_nan=False).encode()).hexdigest())
   with self.assertRaises(FileExistsError):NS['EvidenceLog'](p,lambda v:v)
 def test_evidence_refuses_nonfinite(self):
  with tempfile.TemporaryDirectory() as d:
   log=NS['EvidenceLog'](Path(d)/'log',lambda v:v)
   with self.assertRaises(ValueError):log.append('bad',{'x':float('nan')})
   self.assertEqual(log.count,0)
 def test_gate_routes(self):
  g={'sharedGuard':str(ROOT/'scripts/s11c_guarded_run.py'),'supervisor':str(M/'S11c_d_end_normalization_run.py')};NS['verify_helper_paths'](g);g['sharedGuard']='/tmp/other'
  with self.assertRaises(ValueError):NS['verify_helper_paths'](g)
 def test_worker_invocation(self):
  args=SimpleNamespace(out=Path('/tmp/full-weak/out'),inputs=Path('/tmp/full-weak/input'),gate=Path('/tmp/full-weak/gate'))
  argv=[str(WORKER),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
  gate={'outputDirectory':str(args.out),'command':['supervisor']+argv};NS['verify_invocation'](args,gate,argv)
  with self.assertRaises(ValueError):NS['verify_invocation'](args,gate,argv+['--bad'])
 def test_actual_pressure_ancestry_metadata(self):
  for row in ['U0','U1','U2','THETA_BALANCE','E_W_BALANCE']:
   s=saved('inventory/'+row+'-pressure-source-join-input.json');b=saved('inventory/'+row+'-pressure-bound-join-input.json');a=saved('inventory/'+row+'-affine-pressure-reconstruction-input.json');old=saved('consumer/'+row+'-consumer-input.json')
   self.assertEqual(s['right'],old['raw']);self.assertEqual(b['left'],old['bound']);self.assertEqual(b['right'],a['left'])
 def test_cancelled_identity_is_not_structural_equality(self):
  prefix='inventory/THETA_BALANCE-delta_p_plus-full-coefficient';i=saved(prefix+'-input.json');r=saved(prefix+'-return.json')
  self.assertNotEqual(i['left'],i['right'])
  out=NS['inherited_slot_join'](i['left'],i,r);self.assertFalse(out['functionCalled']);self.assertEqual(out['right'],i['right'])
 def test_slot_join_rejects_wrong_input(self):
  with self.assertRaises(ValueError):NS['inherited_slot_join']('wrong',{'left':'actual','right':'different'},{'cancelled':{'text':'0','srepr':'Integer(0)'}})
 def test_slot_join_rejects_missing_zero(self):
  with self.assertRaises(ValueError):NS['inherited_slot_join']('actual',{'left':'actual','right':'different'},{'cancelled':{'text':'1','srepr':'Integer(1)'}})
 def test_all_actual_slot_identity_arguments(self):
  for row in ['U0','U1','U2','THETA_BALANCE','E_W_BALANCE']:
   old=saved('consumer/'+row+'-consumer-input.json')
   for slot in NS['SLOTS']:
    pre='inventory/'+row+'-'+slot+'-full-coefficient';NS['inherited_slot_join'](old['slotCoefficients'][slot],saved(pre+'-input.json'),saved(pre+'-return.json'))
 def test_actual_child_hash_coverage(self):
  covered=[]
  for row in ['U0','U1','U2','THETA_BALANCE','E_W_BALANCE']:
   v=NS['pressure_child_coverage'](saved('inventory/'+row+'-full-native-partition.json'),saved('inventory/'+row+'-ordered-addresses.json'),row);covered+=v['coveredHashes']
  self.assertEqual(len(covered),12)
 def test_wrong_pressure_child_refused(self):
  p=toy_partition();a={'row':'THETA_BALANCE','slot':'pressure','face':'plus','nativeRowChildHashes':['wrong']}
  with self.assertRaises(ValueError):NS['pressure_child_coverage'](p,[a],'THETA_BALANCE')
 def test_missing_pressure_coverage_refused(self):
  with self.assertRaises(ValueError):NS['pressure_child_coverage'](toy_partition(),[],'THETA_BALANCE')
 def test_zero_row_cannot_claim_pressure_child(self):
  p={'children':[],'pressureChildIndices':[]};a={'row':'U0','slot':'pressure','face':'plus','nativeRowChildHashes':['wrong']}
  with self.assertRaises(ValueError):NS['pressure_child_coverage'](p,[a],'U0')
 def test_factor_join_and_corruption(self):
  a=saved('inventory/THETA_BALANCE-ordered-addresses.json')[0];o=saved('inventory/factors/'+a['fullFactorProof']['proof']+'-operands.json');origin=saved('inventory/U0-ordered-addresses.json')[o['addressId']];NS['factor_address_join'](a,o,origin)
  a=copy.deepcopy(a);a['normalOriginal']={'text':'wrong','srepr':"Symbol('wrong')"}
  with self.assertRaises(ValueError):NS['factor_address_join'](a,o,origin)
 def test_all_factor_actual_arguments(self):
  all_a={a['addressId']:a for row in ['U0','U1','U2','THETA_BALANCE','E_W_BALANCE'] for a in saved('inventory/'+row+'-ordered-addresses.json')}
  for a in all_a.values():
   proof=a['fullFactorProof']['proof'];v=MAN['savedInputs']['inventory/factors/'+proof+'-operands.json'];o=read(v['path']);NS['factor_address_join'](a,o,all_a[o['addressId']])
   self.assertEqual(a['epsilonCount'],0 if a['status'].startswith('EXACT_ZERO_') else 1)
 def test_factor_map_label_not_silently_ignored(self):
  a=saved('inventory/U0-ordered-addresses.json')[0];o=saved('inventory/factors/'+a['fullFactorProof']['proof']+'-operands.json');a=copy.deepcopy(a);a['responseMap']['id']='wrong-face-or-grade'
  with self.assertRaises(ValueError):NS['factor_address_join'](a,o,a)
 def test_cheap_pressure_checks_precede_local_derivation(self):
  self.assertLess(SOURCE.index("J.emit('inherited-pressure-law'"),SOURCE.index('        def derive(entry):'))
 def test_exact_whole_input_ancestry(self):
  old=saved('weak/whole-definition-inputs.json');self.assertEqual(old['tags'],saved('inventory/whole-tag-definitions.json'));self.assertEqual(old['typedDirect'],saved('inventory/typed-direct-objects.json'));self.assertEqual(old['savedReference'],saved('saved/reference/retained-response-census.json'))
 def test_actual_address_and_epsilon_metadata(self):
  routes={a['addressId']:a for a in saved('weak/weak-address-coverage.json')};seen=[];epsilon=saved('inventory/actual-binding-context.json')['restored']['epsilon']
  for row in ['U0','U1','U2','THETA_BALANCE','E_W_BALANCE']:
   splits={slot:saved('inventory/'+row+'-'+slot+'-split.json')['retained'] for slot in NS['SLOTS']}
   for a in saved('inventory/'+row+'-ordered-addresses.json'):
    wr=routes[a['addressId']];slot=('delta_p_' if a['slot']=='pressure' else 'd_w_delta_p_')+a['face'];self.assertEqual(a['consumerOriginal'],splits[slot][str(tuple(a['consumerGrade']))]);self.assertEqual(a['epsilon'],epsilon)
    self.assertEqual(wr['sourceJet'],a['jet']);self.assertEqual(wr['gradeTriple'],[a[k] for k in ['consumerGrade','responseGrade','sourceGrade']])
    n={'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL':1,'NATIVE_MIXED_ITERATION':2}.get(a['component'],0)
    self.assertEqual([t['tagCount'] for t in wr['wholeTagCheck']],[n,n]);seen.append(a['addressId'])
  self.assertEqual(sorted(seen),list(range(13260)))
 def test_saved_pressure_field_joins(self):
  certs=saved('weak/all-coefficient-certificates.json')
  for fid,f in saved('inventory/fields.json').items():
   old=saved('weak/field-'+fid+'-operands.json');p=saved('weak/field-'+fid+'-derivative-class.json');r=saved('weak/field-'+fid+'-reconstruction-input.json')
   self.assertEqual(old['original'],f['field']);self.assertEqual(old['original'],r['left']);self.assertEqual(certs[fid]['polynomial'],p['P0'])
 def test_all_runtime_pins(self):
  for p,h in MAN['sourcePins'].items():self.assertEqual(NS['sha'](p),h,p)
  for a,v in MAN['savedInputs'].items():self.assertEqual(NS['sha'](v['path']),v['sha256'],a)
 def test_no_top_level_science_or_old_run_calls(self):
  for n in TREE.body:
   if isinstance(n,ast.Import):self.assertTrue(all(x.name not in ('sympy','numpy') for x in n.names))
   if isinstance(n,ast.ImportFrom):self.assertNotIn(n.module,('sympy','numpy'))
  self.assertNotIn('runpy.run_path',SOURCE);self.assertNotIn('signal.alarm',SOURCE)
 def test_launcher_scope_and_hook(self):
  t=ast.parse(Path(MAN['launcher']).read_text());vals={x.targets[0].id:ast.literal_eval(x.value) for x in t.body if isinstance(x,ast.Assign) and len(x.targets)==1 and isinstance(x.targets[0],ast.Name) and isinstance(x.value,ast.Constant)}
  self.assertEqual(vals['THREAD'],'01a0e01b-ef84-7192-817f-584cda5d339b');self.assertEqual(vals['PREFIX'],PREFIX)
  self.assertIn('2908local/12pressure',vals['MESSAGE'])

if __name__=='__main__':unittest.main()
