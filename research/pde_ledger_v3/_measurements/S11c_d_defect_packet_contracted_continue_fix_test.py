"""Source, closure and metadata tests only; no scientific imports/restoration."""
import ast,copy,hashlib,importlib.util,json,sys,tempfile,types,unittest
from pathlib import Path
M=Path(__file__).resolve().parent
OLD='S11c_d_defect_packet_contracted_continue';NEW=OLD+'_fix'
sys.path.insert(0,str(M))
def source(suffix='',fixed=False):return (M/((NEW if fixed else OLD)+suffix+'.py')).read_text()
def function(text,name):return next(n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name==name)
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
class Sink:
 def __init__(self):self.rows=[]
 def emit(self,k,v):self.rows.append((k,v))
class ToolingRepair(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.prep=load('fixed_prepare_metadata',M/(NEW+'_prepare.py'));cls.worker=load('fixed_worker_metadata',M/(NEW+'.py'))
 def probe(self,fixed):
  f=function(source('_prepare',fixed),'restore_prefix')
  selected=[copy.deepcopy(n) for n in f.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in {'prior','root','get','newtree','profile','outgoing_node'} for t in n.targets)]
  probe=ast.parse('def probe(J,numeric):\n pass').body[0];probe.body=selected+[ast.Return(value=ast.Name(id='get',ctx=ast.Load()))]
  tree=ast.fix_missing_locations(ast.Module(body=[probe],type_ignores=[]));ns={'ast':ast,'read':lambda p:p};exec(compile(tree,'actual-source-closure','exec'),ns)
  return ns['probe'](types.SimpleNamespace(out=Path('/tmp/source-closure-test')), 'def profile():\n pass\ndef q():\n pass\n')
 def test_original_actual_closure_reproduces_failure(self):
  with self.assertRaisesRegex(TypeError,'FunctionDef'):self.probe(False)('operand')
 def test_fixed_actual_closure_keeps_path(self):
  self.assertEqual(self.probe(True)('operand'),Path('/tmp/source-closure-test/prior/complete/operand.json'))
 def test_fixed_getter_reusable_for_multiple_literal_paths(self):
  get=self.probe(True)
  for n in ['first','nested/second','last']:self.assertEqual(get(n),Path('/tmp/source-closure-test/prior/complete')/(n+'.json'))
 def test_restore_function_only_alpha_renamed(self):
  class Names(ast.NodeTransformer):
   def visit_Name(self,n):
    if n.id=='outgoing_node':n.id='root'
    return n
  self.assertEqual(ast.dump(function(source('_prepare'),'restore_prefix')),ast.dump(Names().visit(function(source('_prepare',True),'restore_prefix'))))
 def test_scientific_prepare_body_unchanged(self):
  a=function(source('_prepare'),'prepare');b=function(source('_prepare',True),'prepare');b.body=b.body[1:]
  self.assertEqual(ast.dump(a),ast.dump(b))
 def test_numeric_run_unchanged(self):self.assertEqual(ast.dump(function(source(),'run')),ast.dump(function(source(fixed=True),'run')))
 def test_earliest_original_numeric_run_unchanged(self):
  p=M/'S11c_d_defect_packet_contracted_numeric.py';self.assertEqual(ast.dump(function(p.read_text(),'run')),ast.dump(function(source(fixed=True),'run')))
 def test_helper_and_tail_stay_original(self):
  t=ast.parse(source('_prepare',True));mods=[n.module for n in t.body if isinstance(n,ast.ImportFrom)]
  self.assertIn(OLD+'_resume',mods);self.assertIn(OLD+'_tail',mods)
 def binding_manifest(self):return {'resumeLibrary':str(M/(OLD+'_resume.py')),'tailLibrary':str(M/(OLD+'_tail.py'))}
 def binding_raw(self,p):return {'numeric-source-contracts.json':{'fragments':{'tail-contributions':{'source':str(p)}}}}
 def test_actual_base_and_module_bindings(self):
  s=Sink();self.prep.tooling_bindings(self.binding_raw(M/'S11c_d_defect_packet_preflight.py'),self.binding_manifest(),s);self.assertTrue(s.rows[0][1]['baseAssignmentsSame'])
 def test_module_mismatch_refuses_after_operand_persistence(self):
  s=Sink();m=self.binding_manifest();m['tailLibrary']='/tmp/wrong.py'
  with self.assertRaisesRegex(ValueError,'module files'):self.prep.tooling_bindings(self.binding_raw(M/'S11c_d_defect_packet_preflight.py'),m,s)
  self.assertEqual(len(s.rows),1)
 def test_base_XY_swap_refuses_after_operand_persistence(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'source.py';p.write_text("def contributions():\n for a in live:\n  z=bounds[a['addressId']]\n  cx,cy,dy=z['Y'],z['X'],z['Yprime']\n")
   s=Sink()
   with self.assertRaisesRegex(ValueError,'X/Y'):self.prep.tooling_bindings(self.binding_raw(p),self.binding_manifest(),s)
   self.assertFalse(s.rows[0][1]['baseAssignmentsSame'])
 def failure_fixture(self,d):
  root=Path(d)/'failed';root.mkdir();p=root/'operand.json';p.write_text('{"saved":true}');h=hashlib.sha256(p.read_bytes()).hexdigest();rec=Path(d)/'record.json';rec.write_text(json.dumps({'root':str(root),'counts':{'files':1,'bytes':p.stat().st_size},'records':{'operand.json':{'sha256':h,'bytes':p.stat().st_size}}}));return root,{'toolingFailureRecord':str(rec)}
 def test_preserved_tree_intact(self):
  with tempfile.TemporaryDirectory() as d:
   _,m=self.failure_fixture(d);self.assertTrue(self.worker.tooling_failure_state(m)['intact'])
 def test_preserved_tree_changed_bytes_refuse(self):
  with tempfile.TemporaryDirectory() as d:
   r,m=self.failure_fixture(d);(r/'operand.json').write_text('{"saved":null}');self.assertFalse(self.worker.tooling_failure_state(m)['intact'])
 def test_preserved_tree_missing_refuses(self):
  with tempfile.TemporaryDirectory() as d:
   r,m=self.failure_fixture(d);(r/'operand.json').unlink();self.assertFalse(self.worker.tooling_failure_state(m)['intact'])
 def test_preserved_tree_extra_refuses(self):
  with tempfile.TemporaryDirectory() as d:
   r,m=self.failure_fixture(d);(r/'extra').write_text('extra');self.assertFalse(self.worker.tooling_failure_state(m)['intact'])
 def test_scientific_imports_absent(self):self.assertNotIn('sympy',sys.modules);self.assertNotIn('mpmath',sys.modules)
 def test_all_new_sources_parse(self):
  for suffix in ('','_prepare','_launch','_test'):ast.parse(source(suffix,True))
 def test_old_pins_still_exact(self):
  m=json.loads((M/(OLD+'_inputs.json')).read_text())
  for p,h in m['sourcePins'].items():self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(),h,p)
 def test_new_gate_honest_flags(self):
  s=source(fixed=True);self.assertIn("g['independentBuildClearance'] is False",s);self.assertIn("g['localToolingExecutionAuthority'] is True",s)
if __name__=='__main__':unittest.main(verbosity=2)
