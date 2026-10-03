#!/usr/bin/env python3
"""Standard-library source/JSON tests only; no scientific restoration or execution."""
import ast
import copy
import hashlib
import itertools
import json
from pathlib import Path
import unittest
M=Path(__file__).resolve().parent
WORKER=M/'S11c_d_defect_weak_ends.py'
MAN=json.loads((M/'S11c_d_defect_weak_ends_inputs.json').read_text())
NAMES=('require','inherit_zero','cell_key','address_metadata_join','literal_field_joins','source_statement_join','verify_build_assessment','verify_helper_paths','verify_invocation','native_phase_exponent','control_eligible','first_shape_grade_parts')
ns={'ast':ast,'hashlib':hashlib,'Path':Path,'ROOT':Path('/var/projects/toy_physics'),'__file__':str(WORKER),'G':((0,0),(1,0),(0,1),(1,1)),
    'ROWS':('U0','U1','U2','THETA_BALANCE','E_W_BALANCE'),'FIELDS':('u_1','u_2','u_3','theta','e_W'),'ZERO':{'text':'0','srepr':'Integer(0)'}}
tree=ast.parse(WORKER.read_text());nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in NAMES]
assert len(nodes)==len(NAMES)
exec(compile(ast.Module(body=nodes,type_ignores=[]),'metadata-only-helpers','exec'),ns)
def read(alias):return json.loads(Path(MAN['savedInputs'][alias]['path']).read_text())
FIELDS=read('inventory/fields.json');COVER={r['addressId']:r for r in read('weak/weak-address-coverage.json')};FACTORS=read('full/new-pressure-factor-arguments.json');WAVES=read('full/new-pressure-wave-arguments.json')
ADDRESSES=sum([read('inventory/'+r+'-ordered-addresses.json') for r in ns['ROWS']],[])
def join(a):return ns['address_metadata_join'](a,COVER[a['addressId']],FIELDS,FACTORS,WAVES)
class SyntheticPolynomial:
 """Integer-only two-variable stand-in; never interprets saved science."""
 def __init__(self,terms):self.terms={g:c for g,c in terms.items() if c}
 def coeff(self,var,power):
  axis=0 if var.terms=={(1,0):1} else 1
  if var.terms not in ({(1,0):1},{(0,1):1}):raise ValueError('synthetic variable')
  return SyntheticPolynomial({tuple(0 if i==axis else n for i,n in enumerate(g)):c
                              for g,c in self.terms.items() if g[axis]==power})
 def __mul__(self,other):
  result={}
  for g,c in self.terms.items():
   for h,d in other.terms.items():
    key=tuple(a+b for a,b in zip(g,h));result[key]=result.get(key,0)+c*d
  return SyntheticPolynomial(result)
 def __add__(self,other):
  result=self.terms.copy()
  for g,c in other.terms.items():result[g]=result.get(g,0)+c
  return SyntheticPolynomial(result)
 def __eq__(self,other):return isinstance(other,SyntheticPolynomial) and self.terms==other.terms

class TestMetadata(unittest.TestCase):
 def test_independent_first_shape_selection(self):
  eta=SyntheticPolynomial({(1,0):1});sigma=SyntheticPolynomial({(0,1):1})
  first=SyntheticPolynomial({(1,0):7,(0,1):11})
  parts=ns['first_shape_grade_parts'](first,eta,sigma,lambda x:x)
  self.assertEqual(parts['height'],SyntheticPolynomial({(1,0):7}))
  self.assertEqual(parts['slope'],SyntheticPolynomial({(0,1):11}))
  self.assertEqual(parts['height']+parts['slope'],first)
 def test_unselected_first_shape_grades_fail_reconstruction(self):
  eta=SyntheticPolynomial({(1,0):1});sigma=SyntheticPolynomial({(0,1):1})
  for g in ((0,0),(1,1),(2,0),(0,2),(2,1)):
   first=SyntheticPolynomial({(1,0):7,(0,1):11,g:13})
   parts=ns['first_shape_grade_parts'](first,eta,sigma,lambda x:x)
   self.assertNotEqual(parts['height']+parts['slope'],first)
 def test_common_depth_height_zero_does_not_erase_slope(self):
  eta=SyntheticPolynomial({(1,0):1});sigma=SyntheticPolynomial({(0,1):1})
  first=SyntheticPolynomial({(0,1):11})
  parts=ns['first_shape_grade_parts'](first,eta,sigma,lambda x:x)
  self.assertEqual(parts['height'],SyntheticPolynomial({}))
  self.assertEqual(parts['slope'],first)
 def test_saved_first_shape_has_two_independent_addends(self):
  # Syntax inspection only: retain the exact saved constructor tree, do not evaluate it.
  expr=ast.parse(read('reference/native-first-shape-input.json')['right']['srepr'],mode='eval').body
  self.assertIsInstance(expr,ast.Call);self.assertEqual(expr.func.id,'Add');self.assertEqual(len(expr.args),2)
  names=[]
  for term in expr.args:
   labels=[n.args[0].value for n in ast.walk(term) if isinstance(n,ast.Call)
           and isinstance(n.func,ast.Name) and n.func.id=='Symbol']
   names.append((labels.count('eta_bg'),labels.count('sigma_W')))
  self.assertEqual(sorted(names),[(0,1),(1,0)])
 def test_first_height_scope_and_persistence(self):
  run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
  source=ast.get_source_segment(WORKER.read_text(),run)
  self.assertNotIn('sigma.name:sp.S.Zero',source)
  self.assertIn('sigma.name:sigma',source)
  self.assertLess(source.index("audit.append('native-first-shape-grade-operands'"),
                  source.index("zero('native-first-shape-independent-grade-reconstruction'"))
  self.assertIn("'slopePointwiseZeroAsserted':False",source)
  self.assertIn('dtn=height_diag.replace',source)
  self.assertNotIn('constant-height-physical-first-shape-',source)
 def test_all_actual_addresses(self):
  for a in ADDRESSES:join(a)
  self.assertEqual(len(ADDRESSES),13260)
 def test_actual_field_proofs(self):
  cert=read('weak/all-coefficient-certificates.json')
  for fid in FIELDS:
   v=[read('weak/field-'+fid+'-'+tail+'.json') for tail in ['operands','polynomial','derivative-class','reconstruction-input','reconstruction-return']]
   ns['literal_field_joins'](fid,FIELDS,cert,*v)
 def test_actual_local_cell_arguments(self):
  cells=read('full/all-local-cells.json');self.assertEqual(len({ns['cell_key'](c) for c in cells}),400)
  for c in cells:
   self.assertEqual(c['identities'][0]['left'],c['coefficient']);self.assertEqual(c['coefficient'],c['polynomial']['original'])
   for z in c['identities']:ns['inherit_zero'](z)
 def test_native_source_assignments(self):
  src=Path(MAN['c2Source']).read_text();du=read('weak/new-weak-duality-and-order.json')['inheritedConvention']
  ns['source_statement_join'](src,'phase1 = '+du['sourcePhase']);ns['source_statement_join'](src,'phase = '+du['profilePhase'])
  for face in ('plus','minus'):
   for st in read('reference/'+face+'-new-native-trace.json')['nativeAssignmentSources'].values():ns['source_statement_join'](src,st)
   ns['source_statement_join'](src,read('reference/'+face+'-final-native-slot-routing.json')['source'])
 def test_all_pins(self):
  for p,h in MAN['sourcePins'].items():self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(),h,p)
 def test_speed_raw_schema(self):
  census=read('full/native-raw-symbol-census.json');self.assertTrue(census['speedAbsent']);self.assertNotIn('c_s0',census['names'])
 def test_full_published_receipts(self):
  receipts=read('full/artifact-index.json')
  for alias,v in MAN['savedInputs'].items():
   if alias.startswith('full/') and alias!='full/artifact-index.json':
    name=alias.split('/',1)[1];r=receipts[name]
    self.assertEqual((r['path'],r['sha256'],r['bytes']),(name,v['sha256'],v['bytes']))
 def test_exact_scope_and_authority(self):
  auth=json.loads(Path(MAN['executionAuthority']).read_text())
  self.assertEqual(auth['scope'],MAN['scope']);self.assertEqual(auth['scienceExecutionsAuthorized'],1)
  self.assertTrue(auth['boundedInstrumentAuthorized']);self.assertTrue(auth['noDeadline']);self.assertFalse(auth['automaticScientificRetry'])
 def test_genuine_saved_bound_hypotheses(self):
  for face in ('plus','minus'):
   self.assertTrue(read('weak/'+face+'-global-height-PV-certificate.json')['normalCoefficientBoundedDirectly'])
  env=read('weak/global-whole-kernel-envelopes.json');self.assertEqual(env['normalWeight'],'(1+|k|+|l|)^3')
  self.assertTrue(env['noComputedResponseIntegral']);self.assertFalse(env['machineMeasureTheoryProof'])
 def test_invocation_identity(self):
  class Args:pass
  a=Args();a.out=Path('/tmp/weak-end-test-out');a.inputs=Path('/tmp/input');a.gate=Path('/tmp/gate')
  argv=[str(WORKER),'--out',str(a.out),'--inputs',str(a.inputs),'--gate',str(a.gate)]
  g={'outputDirectory':str(a.out),'command':['supervisor']+argv};ns['verify_invocation'](a,g,argv)
  with self.assertRaises(ValueError):ns['verify_invocation'](a,g,argv+['--unexpected'])
  g['command'][-1]='/tmp/wrong-gate'
  with self.assertRaises(ValueError):ns['verify_invocation'](a,g,argv)
 def test_negative_route(self):
  for key,value in [('face','bogus'),('targetGrade',[1,1]),('epsilonCount',5),('normalMultiplier',ns['ZERO']),('jet',{'name':'invalid','channel':'invalid'})]:
   a=copy.deepcopy(ADDRESSES[0]);a[key]=value
   with self.assertRaises((ValueError,KeyError)):join(a)
 def test_source_field_substitution_refused(self):
  a=copy.deepcopy(ADDRESSES[0]);a['sourceField']=ns['ZERO']
  with self.assertRaises(ValueError):join(a)
 def test_old_return_must_be_zero(self):
  with self.assertRaises(ValueError):ns['inherit_zero']({'cancelled':{'text':'1','srepr':'Integer(1)'}})
 def test_stale_assignment_refused(self):
  with self.assertRaises(ValueError):ns['source_statement_join']('x=1','x=2')
 def test_field_argument_refused(self):
  fid=next(iter(FIELDS));cert=read('weak/all-coefficient-certificates.json')
  v=[read('weak/field-'+fid+'-'+tail+'.json') for tail in ['operands','polynomial','derivative-class','reconstruction-input','reconstruction-return']]
  v[0]['original']=ns['ZERO']
  with self.assertRaises(ValueError):ns['literal_field_joins'](fid,FIELDS,cert,*v)
 def test_grade_triples(self):
  triples={tuple(tuple(a[k]) for k in ('consumerGrade','responseGrade','sourceGrade')) for a in ADDRESSES};self.assertEqual(len(triples),16)
 def test_all_control_candidate_metadata(self):
  counts={}
  for kind in ('omit-height-contact','reverse-translation-phase','omit-lower-normal-sign'):
   candidates=[a for a in ADDRESSES if ns['control_eligible'](a,kind)];self.assertTrue(candidates,kind);counts[kind]=len(candidates)
   for a in candidates:self.assertEqual(a['status'],'FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED')
  lower=next(a for a in ADDRESSES if a['addressId']==10062)
  self.assertTrue(ns['control_eligible'](lower,'omit-lower-normal-sign'));self.assertEqual(lower['responseGrade'],[0,0])
  height=next(a for a in ADDRESSES if a['addressId']==8034)
  self.assertTrue(ns['control_eligible'](height,'reverse-translation-phase'))
 def test_old_empty_selector_reproduced_without_algebra(self):
  old=[a for a in ADDRESSES if a['face']=='minus' and a['slot']=='normal' and a['responseGrade']==[1,0]
       and a['status']=='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED']
  self.assertEqual(old,[])
  present=[a for a in ADDRESSES if a['face']=='minus' and a['slot']=='normal' and a['responseGrade']==[1,0]]
  self.assertTrue(present);self.assertTrue(all(a['status']=='EXACT_ZERO_CONSUMER' for a in present))
 def test_control_selector_rejects_wrong_slots_and_unknown(self):
  a=copy.deepcopy(next(a for a in ADDRESSES if a['addressId']==10062));a['slot']='pressure'
  self.assertFalse(ns['control_eligible'](a,'omit-lower-normal-sign'))
  with self.assertRaises(ValueError):ns['control_eligible'](a,'invented-control')
 def test_native_phase_binding_on_synthetic_coordinates(self):
  conv=read('weak/new-weak-duality-and-order.json')['inheritedConvention']
  b={'kout':(5,7,11),'kin':(2,7,11),'ko':(5,7,11),'ki':(2,7,11),'X':(3,4,5),'Y':(1,4,5)}
  bind=ns['native_phase_exponent'];self.assertEqual(bind(conv['sourcePhase'],b,1),13)
  self.assertEqual(bind(conv['profilePhase'],b,1),-3)
  moved={**b,'X':(7,4,5),'Y':(5,4,5)}
  self.assertEqual(bind(conv['sourcePhase'],moved,1)-bind(conv['sourcePhase'],b,1),12)
  self.assertEqual(bind(conv['profilePhase'],b,-1),3)
 def test_phase_grammar_rejects_new_functions(self):
  for text in ['sp.exp(__import__("os"))','sp.exp(sp.sin(1))','sp.exp(1/2)','sp.exp(unknown)','sp.exp(1,2)']:
   with self.assertRaises((ValueError,KeyError)):ns['native_phase_exponent'](text,{},1)
 def test_native_zip_refuses_silent_truncation(self):
  source=read('weak/new-weak-duality-and-order.json')['inheritedConvention']['sourcePhase']
  b={'kout':(1,2),'kin':(1,2,3),'X':(1,2,3),'Y':(1,2,3)}
  with self.assertRaises(ValueError):ns['native_phase_exponent'](source,b,1)
 def test_reversed_phase_code_uses_both_derived_limits(self):
  run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_science');src=ast.get_source_segment(WORKER.read_text(),run)
  self.assertIn('reversed_exponent=-translated_exponent',src)
  self.assertIn("reversed_limits['minus'],limits['plus']",src);self.assertIn("reversed_limits['plus'],limits['minus']",src)
  self.assertIn("base*rr['normal']*reversed_limits[side]*Bh",src)
  self.assertIn("('minus','plus') if kind=='reverse-translation-phase'",src)
 def test_proof_cache_groups_do_not_mix_grades(self):
  groups={}
  for a in ADDRESSES:groups.setdefault((a['fullFactorProof']['proof'],a['face'],a['slot']),set()).add(tuple(a['responseGrade']))
  self.assertTrue(all(len(v)==1 for v in groups.values()))
 def test_no_science_at_module_level(self):
  for n in tree.body:
   if isinstance(n,ast.Import):self.assertTrue(all(x.name not in ('sympy','numpy','scipy') for x in n.names))
   if isinstance(n,ast.ImportFrom):self.assertNotIn(n.module,('sympy','numpy','scipy'))
  self.assertNotIn('eval(',WORKER.read_text())
 def test_main_containment_precedes_import(self):
  main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main');source=ast.get_source_segment(WORKER.read_text(),main)
  self.assertLess(source.index("ns['containment']()"),source.index('import sympy'))
  self.assertLess(source.index('verify_gate('),source.index('mkdir('))
 def test_no_old_scientific_function_import(self):
  run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
  self.assertFalse(any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in ('exec','eval','compile','__import__') for n in ast.walk(run)))
 def test_strict_paired_review(self):
  g={k:k for k in ('workerSha256','manifestSha256','launcherSha256','sharedGuardSha256','supervisorSha256')};g['independentBuildClearance']=True
  r={**g,'methodAssessed':True,'allChecksPassed':True,'reports':{e:{'literalVerdict':'CLEAR FOR THIS TRANSLATED WEAK-END BUILD'} for e in ('claude','grok')}}
  ns['verify_build_assessment'](g,r)
  for engine in ('claude','grok'):
   bad=copy.deepcopy(r);bad['reports'][engine]['literalVerdict']='NEEDS REVISION'
   with self.assertRaises(ValueError):ns['verify_build_assessment'](g,bad)
  for key in ('workerSha256','manifestSha256','launcherSha256','sharedGuardSha256','supervisorSha256'):
   bad=copy.deepcopy(r);bad[key]='bad'
   with self.assertRaises(ValueError):ns['verify_build_assessment'](g,bad)
 def test_no_total_deadline_or_uncontained_launcher(self):
  launch=Path(MAN['launcher']).read_text();self.assertIn('s11c_guarded_run.py',launch);self.assertIn('codex_job_watch.py',launch)
  self.assertIsNone(MAN['resources']['durationLimits']);self.assertNotIn('signal.alarm',WORKER.read_text());self.assertNotIn('setitimer',WORKER.read_text())
 def test_no_gate_or_science_run(self):
  self.assertFalse((M/'S11c_d_defect_weak_ends_gate.json').exists())
  self.assertFalse(Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-weak-ends-20261002/diagnostic-01').exists())
if __name__=='__main__':unittest.main()
