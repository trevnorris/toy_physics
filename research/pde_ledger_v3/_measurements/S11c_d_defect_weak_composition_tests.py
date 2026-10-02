"""Standard-library gate, storage, metadata and source-boundary tests; no science."""
import ast,copy,hashlib,json,runpy,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
R=Path('/var/projects/toy_physics');M=R/'research/pde_ledger_v3/_measurements';W=M/'S11c_d_defect_weak_composition.py'
ns=runpy.run_path(str(W),run_name='stdlib_tests_only');sha=ns['sha']
class Tests(unittest.TestCase):
 def test_inert_imports(self):
  tree=ast.parse(W.read_text());allowed={'argparse','ast','hashlib','itertools','json','os','pathlib','re','resource','shutil','sys','time','traceback','types'}
  for n in tree.body:
   if isinstance(n,ast.Import):self.assertTrue(all(a.name in allowed for a in n.names))
   if isinstance(n,ast.ImportFrom):self.assertIn(n.module,allowed)
 def test_science_only_after_containment(self):
  t=W.read_text();self.assertLess(t.index("save(args.out/'containment.json'"),t.index('        import sympy as sp'))
  self.assertNotIn('import sympy',t[:t.index('def main():')])
 def test_no_old_producers_or_integrals(self):
  calls={n.func.attr if isinstance(n.func,ast.Attribute) else n.func.id for n in ast.walk(ast.parse(W.read_text())) if isinstance(n,ast.Call) and isinstance(n.func,(ast.Name,ast.Attribute))}
  self.assertFalse(calls & {'build_case','build_face','kernel_bridge','reference_pressure_kernels','shape_coefficients','integrate','quad','lambdify','solve','alarm','setitimer'})
 def test_atomic_no_overwrite(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'receipt.json';ns['save'](p,{'original':1})
   with self.assertRaises(FileExistsError):ns['save'](p,{'original':2})
   self.assertEqual(json.loads(p.read_text()),{'original':1})
 def test_nonboolean_refusal(self):
  for v in [False,None,1,'true']:
   with self.subTest(v=v),self.assertRaises(ValueError):ns['require'](v,'literal true only')
 def test_selector_and_refusal(self):
  base={'addressId':4,'component':'NATIVE_FLAT','face':'plus','slot':'pressure','sourceGrade':[1,0],'consumerGrade':[0,0],'status':'FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED','usable':True}
  other={**base,'addressId':2};self.assertEqual(ns['select_control_address']([base,other],'NATIVE_FLAT','plus','pressure',(1,0),(0,0),lambda a:a['usable'])['addressId'],2)
  for key,value in [('component','NATIVE_HEIGHT'),('face','minus'),('slot','normal'),('sourceGrade',[0,0]),('consumerGrade',[1,0]),('status','EXACT_ZERO_SOURCE_JET'),('usable',False)]:
   with self.subTest(key=key),self.assertRaises(ValueError):ns['select_control_address']([{**base,key:value}],'NATIVE_FLAT','plus','pressure',(1,0),(0,0),lambda a:a['usable'])
 def test_invocation(self):
  a=SimpleNamespace(out=Path('/tmp/weak-test-out'),inputs=Path('/tmp/input'),gate=Path('/tmp/gate'));argv=[str(W),'--out',str(a.out),'--inputs',str(a.inputs),'--gate',str(a.gate)];g={'outputDirectory':str(a.out),'command':['python',*argv]};ns['verify_invocation'](a,g,argv)
  for argv2 in [argv+['--extra'],argv[:-1],['wrong',*argv[1:]]]:
   with self.assertRaises(ValueError):ns['verify_invocation'](a,g,argv2)
  with self.assertRaises(ValueError):ns['verify_invocation'](a,{**g,'outputDirectory':'/tmp/wrong'},argv)
 def test_gate_fixture_and_mutations(self):
  with tempfile.TemporaryDirectory() as td:
   td=Path(td);manifest={'sourcePins':{},'launcher':str(M/'S11c_d_defect_weak_composition_launch.py'),'methodPath':str(M/'S11c_d_defect_weak_composition_method.md'),'methodRecord':str(M/'S11c_d_defect_weak_composition_review_record.json'),'reviewRecordWillBe':str(td/'synthetic-review.json'),'scope':{'syntheticTestOnly':True}}
   authority=td/'authority.json';authority.write_text(json.dumps({'boundedInstrumentAuthorized':True,'automaticScientificRetry':False}))
   mp=td/'manifest.json';mp.write_text(json.dumps(manifest));g={'status':'READY_FOR_ONE_WEAK_COMPOSITION_INSTRUMENT','workerSha256':sha(W),'manifestSha256':sha(mp),'sourcePins':{},'sharedGuard':str(R/'scripts/s11c_guarded_run.py'),'supervisor':str(M/'S11c_d_end_normalization_run.py'),'launcher':manifest['launcher'],'buildReviewRecord':manifest['reviewRecordWillBe'],'authority':str(authority),'durationLimits':None,'scientificRunsAuthorized':1,'scope':manifest['scope']}
   for k in ['sharedGuard','supervisor','launcher','authority']:g[k+'Sha256']=sha(g[k])
   review={'syntheticTestFixtureNotClearance':True,'independentBuildClearance':True,'methodAssessed':True,'allChecksPassed':True,'reports':{e:{'literalVerdict':'CLEAR FOR THIS GLOBAL WEAK-COMPOSITION BUILD'} for e in ['claude','grok']},'methodSha256':sha(manifest['methodPath'])}
   for k in ['workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256','launcherSha256']:review[k]=g[k]
   rp=Path(g['buildReviewRecord']);rp.write_text(json.dumps(review));g['buildReviewRecordSha256']=sha(rp);gp=td/'gate.json';gp.write_text(json.dumps(g));ns['verify_gate'](gp,mp,manifest)
   for k,v in [('status','WRONG'),('workerSha256','0'*64),('manifestSha256','0'*64),('sharedGuard','/tmp/fake-guard'),('supervisor','/tmp/fake-supervisor'),('durationLimits',5),('scientificRunsAuthorized',2),('scope',{})]:
    bad={**g,k:v};gp.write_text(json.dumps(bad))
    with self.subTest(key=k),self.assertRaises((ValueError,FileNotFoundError)):ns['verify_gate'](gp,mp,manifest)
   for key in ['independentBuildClearance','methodAssessed','allChecksPassed','workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256','launcherSha256','methodSha256']:
    rr=copy.deepcopy(review);rr[key]=False if isinstance(rr[key],bool) else '0'*64;rp.write_text(json.dumps(rr));gg={**g,'buildReviewRecordSha256':sha(rp)};gp.write_text(json.dumps(gg))
    with self.subTest(reviewKey=key),self.assertRaises(ValueError):ns['verify_gate'](gp,mp,manifest)
 def test_actual_metadata_contracts_no_restoration(self):
  man=json.loads((M/'S11c_d_defect_weak_composition_inputs.json').read_text());raw={n:json.loads(Path(v['path']).read_text()) for n,v in man['savedInputs'].items()};keys=('sha256','original','mapped','map','symbolAssumptions','flatSupport','frequency','positiveRegulatorContinuation');count=0
  for n,v in raw.items():
   if not n.endswith('ordered-addresses.json'):continue
   for a in v:
    p=raw['inventory/factors/'+a['fullFactorProof']['proof']+'-operands.json'];self.assertEqual(a['normalOriginal'],p['addressNormalOriginal']);self.assertEqual(a['fullFactorProof']['completeNormalMap'],p['requiredMap'])
    for key in keys:self.assertEqual(a['responseMap'][key],p['actualResponseMap'][key])
    for kind in ['source','consumer']:self.assertEqual(a[kind+'Field'],raw['inventory/fields.json'][a[kind+'Transform']['coefficientId']]['field'])
    count+=1
  self.assertEqual(count,13260)
  for p,h in man['sourcePins'].items():self.assertEqual(sha(p),h)
 def test_no_readiness_or_runtime_yet(self):
  self.assertFalse((M/'S11c_d_defect_weak_composition_gate.json').exists());self.assertFalse((R/'_scratch/s11c/s11c-defect-weak-composition-20261002/diagnostic-01').exists())
if __name__=='__main__':unittest.main()
