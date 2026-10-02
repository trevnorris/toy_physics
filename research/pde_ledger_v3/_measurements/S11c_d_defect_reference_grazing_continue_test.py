#!/usr/bin/env python3
"""Stdlib representation, source-tail and gate tests; no scientific imports."""
import ast
import copy
import hashlib
import json
from pathlib import Path
import runpy
import tempfile
import unittest
M=Path(__file__).resolve().parent
W=M/'S11c_d_defect_reference_grazing_continue.py'
BASE=M/'S11c_d_defect_reference_grazing.py'
NS=runpy.run_path(str(W),run_name='stdlib_test_only')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

class ExactConstants(unittest.TestCase):
 def test_observed_representation_shape(self):
  r=NS['exact_rational_constant']('Add(Integer(-9), Integer(9))')
  self.assertEqual((r['numerator'],r['denominator'],r['exactZero']),(0,1,True));self.assertEqual(len(r['steps']),3)
 def test_rational_nested(self):
  r=NS['exact_rational_constant']('Add(Rational(1, 2), Mul(Integer(-1), Pow(Integer(2), Integer(-1))))')
  self.assertTrue(r['exactZero'])
 def test_nonzero_not_promoted(self):
  r=NS['exact_rational_constant']('Add(Integer(-9), Integer(10))');self.assertFalse(r['exactZero']);self.assertEqual(r['numerator'],1)
 def test_reject_unsupported(self):
  for x in ["Symbol('x')","Add(Symbol('x'),Mul(Integer(-1),Symbol('x')))","Float('0.0')",'nan','oo','Integer(True)','Integer(0.0)','Rational(1,0)','Pow(Integer(0),Integer(-1))','Pow(Integer(2),Rational(1,2))','Pow(Integer(2),Integer(17))','Add()','Integer(0,1)','Integer(value=0)',"__import__('os')",'1+1']:
   with self.subTest(x=x):
    with self.assertRaises((ValueError,ZeroDivisionError)):NS['exact_rational_constant'](x)
 def test_no_eval(self):
  node=next(n for n in ast.parse(W.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='exact_rational_constant')
  calls={ast.unparse(n.func) for n in ast.walk(node) if isinstance(n,ast.Call)}
  self.assertFalse(calls&{'eval','exec','sympify','sp.sympify','float'})

class SourceReuse(unittest.TestCase):
 def test_exact_original_tail(self):
  tree=ast.parse(BASE.read_text());run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
  tail=NS['tail_ast'](BASE.read_text());self.assertEqual(len(tail),70)
  self.assertEqual(ast.dump(tail[0]),ast.dump(next(n for n in run.body if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and len(n.value.args)>1 and isinstance(n.value.args[1],ast.Constant) and n.value.args[1].value=='new-tail-polynomial')))
  self.assertEqual([ast.dump(n) for n in tail],[ast.dump(n) for n in run.body[-70:]])
 def test_tail_does_not_execute_prefix(self):
  source="def run_science(manifest,J,helpers):\n    forbidden()\n    nonnegative_polynomial(J,'new-tail-polynomial',1,[])\n    return 'tail'\n"
  calls=[];env={'J':object(),'nonnegative_polynomial':lambda *a:calls.append(a[1]),'forbidden':lambda:(_ for _ in ()).throw(AssertionError('prefix replay'))}
  self.assertEqual(NS['tail_callable'](source,env)(),'tail');self.assertEqual(calls,['new-tail-polynomial'])
 def test_marker_missing_or_duplicate_refused(self):
  for source in ["def run_science(m,J,h):\n return 0\n", "def run_science(m,J,h):\n nonnegative_polynomial(J,'new-tail-polynomial',0,[])\n nonnegative_polynomial(J,'new-tail-polynomial',0,[])\n"]:
   with self.assertRaises(ValueError):NS['tail_ast'](source)
 def test_no_top_level_science(self):
  self.assertNotIn('sp',NS)
  imports=[ast.unparse(n) for n in ast.parse(W.read_text()).body if isinstance(n,(ast.Import,ast.ImportFrom))]
  self.assertFalse(any('sympy' in n or 'numpy' in n for n in imports))
  main=next(n for n in ast.parse(W.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='main');s=ast.unparse(main)
  self.assertLess(s.index("helpers['containment']()"),s.index('import sympy'))
 def test_no_completed_body_call(self):
  tree=ast.parse(W.read_text());calls={ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n,ast.Call)}
  self.assertFalse(calls&{'run_science','sp.solve','sp.integrate','sp.Integral','sp.limit','sp.Poly','sp.nsolve'})
  self.assertFalse(any('alarm' in x or 'setitimer' in x for x in calls))
 def test_all_tail_global_context_names_supplied(self):
  tail=ast.Module(body=NS['tail_ast'](BASE.read_text()),type_ignores=[])
  loads={n.id for n in ast.walk(tail) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Load)}
  stores={n.id for n in ast.walk(tail) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Store)}
  args={n.arg for n in ast.walk(tail) if isinstance(n,ast.arg)}
  needed=loads-stores-args-{'list','str','range'}
  run=next(n for n in ast.parse(W.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='run_continuation')
  supplied=set(NS['BASE_DEFINITIONS'])
  for n in ast.walk(run):
   if isinstance(n,ast.Call) and ast.unparse(n.func)=='context.update':supplied.update(k.arg for k in n.keywords)
   if isinstance(n,ast.Assign):
    for target in n.targets:
     if isinstance(target,ast.Subscript) and ast.unparse(target.value)=='context' and isinstance(target.slice,ast.Constant):supplied.add(target.slice.value)
  symbols=next(n for n in ast.walk(run) if isinstance(n,ast.For) and ast.unparse(n.target)=='(local, saved)')
  supplied.update(ast.literal_eval(pair.elts[0]) for pair in symbols.iter.elts);supplied.update(['A','R'])
  self.assertFalse(needed-supplied,needed-supplied)
 def test_launcher_parses(self):ast.parse((M/'S11c_d_defect_reference_grazing_continue_launch.py').read_text())

class Gates(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.p=Path(self.tmp.name)
  def put(n,v):
   p=self.p/n;p.write_text(json.dumps(v));return str(p)
  self.put=put
  bw=put('baseworker.json',{'base':True});bm=put('baseManifest.json',{'source':True});launcher=put('launcher.json',{'launcher':True});guard=put('guard.json',{'guard':True});sup=put('supervisor.json',{'supervisor':True})
  pins={str(W):sha(W),bw:sha(bw),bm:sha(bm),launcher:sha(launcher),guard:sha(guard),sup:sha(sup)}
  self.man={'scope':'bounded','sourcePins':pins,'baseWorker':bw,'priorTree':{'one.json':{'sha256':'saved','bytes':10}}}
  self.mp=put('manifest.json',self.man)
  build={'independentBuildClearance':True,'workerSha256':sha(bw),'manifestSha256':sha(bm),'reports':{e:{'literalVerdict':'CLEAR FOR THIS BOUNDED REFERENCE-GRAZING BUILD'} for e in ('claude','grok')}}
  paths={'sharedGuard':guard,'supervisor':sup,'launcher':launcher,'baseGate':put('basegate.json',{'workerSha256':sha(bw),'manifestSha256':sha(bm)}),'baseManifest':bm,'baseBuildRecord':put('build.json',build),'authority':put('authority.json',{'savedEvidenceContinuationAuthorized':True,'automaticScientificRetry':False})}
  self.g={'status':'READY_FOR_ONE_SAVED_REFERENCE_CONTINUATION','workerSha256':sha(W),'manifestSha256':sha(self.mp),'sourcePins':pins,'continuationIndependentBuildClearance':False,'localToolingRepairAuthorized':True,'scientificRunsAuthorized':1,'durationLimits':None,'scope':'bounded','priorTree':self.man['priorTree'],**paths}
  self.g.update({k+'Sha256':sha(p) for k,p in paths.items()});self.gp=self.p/'gate.json'
 def verify(self,g=None):
  self.gp.write_text(json.dumps(self.g if g is None else g));return NS['verify_gate'](self.gp,self.mp,self.man)
 def test_good(self):self.assertEqual(self.verify()['scope'],'bounded')
 def test_bad_gate_and_resource_authority(self):
  for key in ['status','workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256','launcherSha256','baseGateSha256','baseManifestSha256','baseBuildRecordSha256','authoritySha256','scope','scientificRunsAuthorized','durationLimits','priorTree','sourcePins']:
   with self.subTest(key=key):
    g=copy.deepcopy(self.g);g[key]='wrong'
    with self.assertRaises((ValueError,TypeError)):self.verify(g)
 def test_no_fresh_clear_claim(self):
  self.g['continuationIndependentBuildClearance']=True
  with self.assertRaises(ValueError):self.verify()
 def test_no_local_authority(self):
  self.g['localToolingRepairAuthorized']=False
  with self.assertRaises(ValueError):self.verify()
 def test_prior_literal_required(self):
  p=Path(self.g['baseBuildRecord']);v=json.loads(p.read_text());v['reports']['grok']['literalVerdict']='NEEDS REVISION';p.write_text(json.dumps(v));self.g['baseBuildRecordSha256']=sha(p)
  with self.assertRaises(ValueError):self.verify()
 def test_prior_worker_join(self):
  p=Path(self.g['baseBuildRecord']);v=json.loads(p.read_text());v['workerSha256']='other';p.write_text(json.dumps(v));self.g['baseBuildRecordSha256']=sha(p)
  with self.assertRaises(ValueError):self.verify()
 def test_changed_actual_launcher(self):
  Path(self.g['launcher']).write_text('changed')
  with self.assertRaises(ValueError):self.verify()
 def test_no_retry_authority(self):
  p=Path(self.g['authority']);v=json.loads(p.read_text());v['automaticScientificRetry']=True;p.write_text(json.dumps(v));self.g['authoritySha256']=sha(p)
  with self.assertRaises(ValueError):self.verify()

if __name__=='__main__':unittest.main()
