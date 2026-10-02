#!/usr/bin/env python3
"""Standard-library gate/source tests only; no scientific imports or operands."""
import ast,copy,hashlib,json,runpy,tempfile,unittest
from pathlib import Path
M=Path(__file__).resolve().parent
W=M/'S11c_d_defect_closed_grazing.py'
NS=runpy.run_path(str(W),run_name='tooling_test_only')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
class ToolingTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.p=Path(self.tmp.name)
  def write(name,obj):
   f=self.p/name;f.write_text(json.dumps(obj));return str(f)
  self.man={'scope':'test-scope','sourcePins':{str(W):sha(W)},'methodRecord':write('method.json',{'jointIndependentMethodClearance':True,'allChecksPassed':True})}
  self.mp=write('manifest.json',self.man)
  self.g={'status':'READY_FOR_ONE_CLOSED_GRAZING_INSTRUMENT','workerSha256':sha(W),'manifestSha256':sha(self.mp),'sourcePins':self.man['sourcePins'],'scope':'test-scope','scientificRunsAuthorized':1,'durationLimits':None,'methodRecord':self.man['methodRecord']}
  for k in ['sharedGuard','supervisor']:self.g[k]=write(k+'.json',{'inert':k})
  self.g['authority']=write('authority.json',{'boundedInstrumentAuthorized':True,'automaticScientificRetry':False})
  for k in ['sharedGuard','supervisor','methodRecord','authority']:self.g[k+'Sha256']=sha(self.g[k])
  b={'independentBuildClearance':True,'allChecksPassed':True,'reports':{e:{'literalVerdict':'CLEAR FOR THIS BOUNDED CLOSED-GRAZING BUILD'} for e in ['claude','grok']},**{k:self.g[k] for k in ['workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256']}}
  self.g['buildReviewRecord']=write('build.json',b);self.g['buildReviewRecordSha256']=sha(self.g['buildReviewRecord']);self.gp=self.p/'gate.json'
 def verify(self,g=None):
  self.gp.write_text(json.dumps(self.g if g is None else g));return NS['verify_gate'](self.gp,self.mp,self.man)
 def test_good(self):self.assertEqual(self.verify()['scope'],'test-scope')
 def test_ten_identity_refusals(self):
  keys=['status','workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256','methodRecordSha256','buildReviewRecordSha256','authoritySha256','scope','scientificRunsAuthorized']
  for k in keys:
   with self.subTest(k=k):
    g=copy.deepcopy(self.g);g[k]='wrong'
    with self.assertRaises((ValueError,KeyError)):self.verify(g)
 def test_no_deadline(self):
  self.g['durationLimits']=1
  with self.assertRaises(ValueError):self.verify()
 def test_literal_review_required(self):
  p=Path(self.g['buildReviewRecord']);b=json.loads(p.read_text());b['reports']['grok']['literalVerdict']='NEEDS REVISION';p.write_text(json.dumps(b));self.g['buildReviewRecordSha256']=sha(p)
  with self.assertRaises(ValueError):self.verify()
 def test_helper_join_required(self):
  p=Path(self.g['buildReviewRecord']);b=json.loads(p.read_text());b['sharedGuardSha256']='other';p.write_text(json.dumps(b));self.g['buildReviewRecordSha256']=sha(p)
  with self.assertRaises(ValueError):self.verify()
 def test_no_author_clear(self):
  p=Path(self.g['buildReviewRecord']);b=json.loads(p.read_text());b['independentBuildClearance']=False;p.write_text(json.dumps(b));self.g['buildReviewRecordSha256']=sha(p)
  with self.assertRaises(ValueError):self.verify()
 def test_authority_refusal(self):
  p=Path(self.g['authority']);a=json.loads(p.read_text());a['boundedInstrumentAuthorized']=False;p.write_text(json.dumps(a));self.g['authoritySha256']=sha(p)
  with self.assertRaises(ValueError):self.verify()
 def test_source_pin_refusal(self):
  self.g['sourcePins']={str(W):'wrong'}
  with self.assertRaises(ValueError):self.verify()
 def test_helpers_exact(self):
  source=(M/'S11c_d_defect_raw_increment.py').read_text();selected=NS['helper_ast'](source)
  old={n.name:n for n in ast.parse(source).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
  self.assertEqual({n.name for n in selected.body},set(NS['HELPERS']))
  for n in selected.body:self.assertEqual(ast.dump(n),ast.dump(old[n.name]))
 def test_no_top_level_science(self):
  tree=ast.parse(W.read_text());self.assertNotIn('sp',NS)
  imports=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom))]
  self.assertFalse(any('sympy' in ast.unparse(n) or 'numpy' in ast.unparse(n) for n in imports))
  main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main');body=ast.unparse(main)
  self.assertLess(body.index("ns['containment']()"),body.index('import sympy'))
 def test_no_integral_solver_or_producer_calls(self):
  tree=ast.parse(W.read_text());calls={ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n,ast.Call)}
  forbidden={'sp.integrate','sp.Integral','sp.limit','sp.solve','sp.nsolve','sp.Poly','sp.series','sp.lambdify'}
  self.assertFalse(calls&forbidden)
  self.assertFalse(any('alarm' in x or 'setitimer' in x for x in calls))
 def test_launch_syntax(self):ast.parse((M/'S11c_d_defect_closed_grazing_launch.py').read_text())
if __name__=='__main__':unittest.main()
