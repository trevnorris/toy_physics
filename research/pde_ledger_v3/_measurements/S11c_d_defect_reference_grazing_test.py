#!/usr/bin/env python3
"""Standard-library gate/source tests only; no scientific imports or operands."""
import ast,copy,hashlib,json,runpy,tempfile,unittest
from pathlib import Path
M=Path(__file__).resolve().parent
W=M/'S11c_d_defect_reference_grazing.py'
NS=runpy.run_path(str(W),run_name='tooling_test_only')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
class ToolingTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.p=Path(self.tmp.name)
  def write(name,obj):
   f=self.p/name;f.write_text(json.dumps(obj));return str(f)
  self.man={'scope':'test-scope','sourcePins':{str(W):sha(W)},'methodPath':write('method-source.json',{'proposed':'bounded'}),'methodRecord':None,'launcher':write('launcher.json',{'inert':'launcher'})}
  self.man['methodRecord']=write('method.json',{'status':'PAIRED_INDEPENDENT_METHOD_CLEARANCE_SOURCE_ONLY','allChecksPassed':True,'methodSha256':sha(self.man['methodPath']),'reviewers':{e:{'literalVerdict':'CLEAR FOR THIS BOUNDED REFERENCE-GRAZING METHOD'} for e in ['claude','grok']}})
  self.mp=write('manifest.json',self.man)
  self.g={'status':'READY_FOR_ONE_REFERENCE_GRAZING_INSTRUMENT','workerSha256':sha(W),'manifestSha256':sha(self.mp),'sourcePins':self.man['sourcePins'],'scope':'test-scope','scientificRunsAuthorized':1,'durationLimits':None,'methodRecord':self.man['methodRecord'],'launcher':self.man['launcher']}
  for k in ['sharedGuard','supervisor']:self.g[k]=write(k+'.json',{'inert':k})
  self.g['authority']=write('authority.json',{'boundedInstrumentAuthorized':True,'automaticScientificRetry':False})
  for k in ['sharedGuard','supervisor','launcher','methodRecord','authority']:self.g[k+'Sha256']=sha(self.g[k])
  b={'independentBuildClearance':True,'allChecksPassed':True,'reports':{e:{'literalVerdict':'CLEAR FOR THIS BOUNDED REFERENCE-GRAZING BUILD'} for e in ['claude','grok']},**{k:self.g[k] for k in ['workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256','launcherSha256']}}
  self.g['buildReviewRecord']=write('build.json',b);self.g['buildReviewRecordSha256']=sha(self.g['buildReviewRecord']);self.gp=self.p/'gate.json'
 def verify(self,g=None):
  self.gp.write_text(json.dumps(self.g if g is None else g));return NS['verify_gate'](self.gp,self.mp,self.man)
 def test_good(self):self.assertEqual(self.verify()['scope'],'test-scope')
 def test_eleven_identity_refusals(self):
  keys=['status','workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256','launcherSha256','methodRecordSha256','buildReviewRecordSha256','authoritySha256','scope','scientificRunsAuthorized']
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
 def test_launcher_review_join_required(self):
  p=Path(self.g['buildReviewRecord']);original=json.loads(p.read_text())
  for change in ['missing','mismatch']:
   with self.subTest(change=change):
    b=copy.deepcopy(original)
    if change=='missing':b.pop('launcherSha256')
    else:b['launcherSha256']='different-reviewed-launcher'
    p.write_text(json.dumps(b));self.g['buildReviewRecordSha256']=sha(p)
    with self.assertRaises((KeyError,ValueError)):self.verify()
 def test_changed_launcher_refused(self):
  Path(self.g['launcher']).write_text('changed')
  with self.assertRaises(ValueError):self.verify()
 def test_launcher_route_refused(self):
  p=self.p/'other-launcher.json';p.write_text(Path(self.g['launcher']).read_text());self.g['launcher']=str(p)
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
 def test_launch_syntax(self):ast.parse((M/'S11c_d_defect_reference_grazing_launch.py').read_text())
 def test_grade_selection_tracks_actual_order(self):
  select=NS['mode_value']
  self.assertEqual(select([[1,1,'mixed'],[0,0,'flat']],['answer','other'],(1,1)),'answer')
  for modes,values in [([[0,0]],['missing']),([[1,1],[1,1]],['a','b']),([[1,1]],[])]:
   with self.subTest(modes=modes,values=values):
    with self.assertRaises(ValueError):select(modes,values,(1,1))
 def test_new_control_evidence_uses_serializable_keys(self):
  tree=ast.parse(W.read_text());run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
  # Actual literal-key dictionaries are safe; symbol maps are saved as pair lists.
  for n in ast.walk(run):
   if isinstance(n,ast.Call) and ast.unparse(n.func)=='J.emit' and len(n.args)>1 and isinstance(n.args[1],ast.Dict):
    self.assertTrue(all(isinstance(k,ast.Constant) and isinstance(k.value,str) for k in n.args[1].keys))
 def test_actual_saved_middle_symbol_census(self):
  from fractions import Fraction
  # Parse metadata syntax, never evaluate its symbolic constructors.
  manifest=json.loads((M/'S11c_d_defect_reference_grazing_inputs.json').read_text())
  saved=json.loads(Path(manifest['savedFiles']['physical-sheet.json']['path']).read_text())
  expression=ast.parse(saved['heightRoute']['srepr'],mode='eval')
  names={n.args[0].value for n in ast.walk(expression) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='Symbol'}
  self.assertEqual(names,{'k','increment_transfer','increment_effective_bulk_speed'})
  class S:
   def __init__(self,name):self.name=name
  symbols={S(n) for n in names}
  run=next(n for n in ast.parse(W.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
  node=next(n for n in ast.walk(run) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='rootbindings' for t in n.targets))
  env={'k':Fraction(3,2),'m':Fraction(30,13),'cs':Fraction(7,5)}
  exec(compile(ast.Module(body=[node],type_ignores=[]),'actual-rootbindings','exec'),env)
  mapped=NS['exact_named_map'](symbols,env['rootbindings'])
  byname={s.name:v for s,v in mapped.items()}
  self.assertEqual(byname['k']+byname['increment_transfer'],env['m'])
  old=dict(env['rootbindings']);old['H']=old.pop('increment_transfer')
  with self.assertRaises(ValueError):NS['exact_named_map'](symbols,old)
  with self.assertRaises(ValueError):NS['exact_named_map']([S('k'),S('k')],{'k':1})
  with self.assertRaises(ValueError):NS['exact_named_map']([S('unknown')],{'k':1})
 def test_trace_control_uses_actual_corrupt_coefficient(self):
  from fractions import Fraction
  from types import SimpleNamespace
  from unittest.mock import patch
  fn=NS['trace_subtraction_coefficient']
  with patch.dict(fn.__globals__,{'sp':SimpleNamespace(cancel=lambda x:x)}):
   native=fn(Fraction(6),Fraction(5),Fraction(3));corrupt=fn(Fraction(-6),Fraction(5),Fraction(3))
  self.assertEqual(native,Fraction(-10));self.assertEqual(corrupt,Fraction(10))
  tree=ast.parse(W.read_text());run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
  assignments={t.id:n.value for n in ast.walk(run) if isinstance(n,ast.Assign) for t in n.targets if isinstance(t,ast.Name)}
  call=assignments['corruptTraceCoefficient'];self.assertEqual(ast.unparse(call.func),'trace_subtraction_coefficient')
  self.assertEqual([ast.unparse(x) for x in call.args],['corruptLowerT','actualSlope','eta'])
  self.assertIn('corruptTraceCoefficient',ast.unparse(assignments['corruptTraceScalar']))
  self.assertIn('corruptTraceScalar',ast.unparse(assignments['wrongLower']))
 def test_native_height_constant_is_derived(self):
  run=next(n for n in ast.parse(W.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
  assignment=next(n for n in ast.walk(run) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='height_constant' for t in n.targets))
  self.assertEqual(ast.unparse(assignment.value),"tr['height'].subs(saved_profile, 0)")
  namespaces=[n for n in ast.walk(run) if isinstance(n,ast.Dict) and any(isinstance(k,ast.Constant) and k.value=='height_constant' for k in n.keys)]
  self.assertTrue(namespaces)
  for n in namespaces:
   value=next(v for k,v in zip(n.keys,n.values) if isinstance(k,ast.Constant) and k.value=='height_constant')
   self.assertEqual(ast.unparse(value),'height_constant')
 def test_method_literal_is_required(self):
  p=Path(self.g['methodRecord']);v=json.loads(p.read_text());v['reviewers']['grok']['literalVerdict']='NEEDS REVISION';p.write_text(json.dumps(v));self.g['methodRecordSha256']=sha(p)
  with self.assertRaises(ValueError):self.verify()
 def test_method_source_is_pinned(self):
  Path(self.man['methodPath']).write_text('changed')
  with self.assertRaises(ValueError):self.verify()
 def test_actual_native_assignment_namespaces(self):
  from types import SimpleNamespace
  helper=(M/'S11c_d_defect_raw_increment.py').read_text()
  tree=ast.parse(helper);nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('function_source','assignment_source','require')]
  env={'ast':ast};exec(compile(ast.Module(body=nodes,type_ignores=[]),'native_fragments','exec'),env)
  class E:
   def __init__(self,label):self.label=label
   def __add__(self,other):return E('sum')
   __radd__=__add__
   def __mul__(self,other):return E('product')
   __rmul__=__mul__
   def xreplace(self,mapping):return E('replaced')
   def subs(self,mapping,simultaneous=False):
    if not simultaneous:raise AssertionError('native simultaneous route')
    return E('final')
  source=(M.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
  ns={k:E(k) for k in ['qo','qi','MIDDLE_Q','normal_output','value_coefficient','height_constant','height_hat','height_kernel']}
  ns.update(sp=SimpleNamespace(Matrix=lambda rows:rows),left={},right={})
  for name in ['normal_input','normal_middle','trace_two','trace_three']:
   exec(env['assignment_source'](source,'reference_pressure_kernels',name),ns)
  self.assertEqual([len(row) for row in ns['trace_three']],[3,3,3])
  ns={'trace_map':{'REFERENCE_VALUE_SOLVE':E('solve'),'PHYSICAL_PRESSURE_TARGET':E('target')},'pressure':E('P'),'jet_slot':E('J'),'normal_jet':E('N')}
  exec(env['assignment_source'](source,'build_face','reference_pressure'),ns)
  self.assertEqual(ns['reference_pressure'].label,'final')
if __name__=='__main__':unittest.main()
