"""Task-local lossless encoder/prefix checks, stdlib only; no scientific restoration."""
import ast,importlib.util,json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
M=Path(__file__).resolve().parent;P='S11c_d_first_order_transverse_flux_continue';OLD='S11c_d_first_order_transverse_flux';spec=importlib.util.spec_from_file_location('flux_continuation',M/(P+'.py'));w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
def function(path,name):return next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
def tail(path):
 f=function(path,'scientific_work');i=next(i for i,n in enumerate(f.body) if isinstance(n,ast.For) and ast.unparse(n.iter)=='native.items()');return ast.dump(ast.Module(body=f.body[i:],type_ignores=[]),include_attributes=False)
class FunctionStub:
 def __str__(self):return 'native_f'
class Base:
 def encode(self,v):
  if isinstance(v,dict):return {k:self.encode(x) for k,x in v.items()}
  if isinstance(v,(tuple,list)):return [self.encode(x) for x in v]
  return v
class Tooling(unittest.TestCase):
 def test_no_scientific_import(self):self.assertNotIn('sympy',sys.modules);self.assertNotIn('numpy',sys.modules)
 def test_compiles(self):
  for suffix in ['.py','_launch.py']:compile((M/(P+suffix)).read_text(),P+suffix,'exec')
 def test_function_dispatch(self):
  j=w.make_journal(Base,FunctionStub,SimpleNamespace(srepr=lambda v:"Function('native_f')"))();self.assertEqual(j.encode(FunctionStub()),{'text':'native_f','srepr':"Function('native_f')"})
 def test_nested_function_registry(self):
  j=w.make_journal(Base,FunctionStub,SimpleNamespace(srepr=lambda v:"Function('native_f')"))();v=j.encode({'knownDimensions':[(FunctionStub(),(1,0,0))]});self.assertEqual(json.loads(json.dumps(v))['knownDimensions'][0],[{'text':'native_f','srepr':"Function('native_f')"},[1,0,0]])
 def test_existing_value_preservation(self):
  j=w.make_journal(Base,FunctionStub,SimpleNamespace(srepr=lambda v:None))();v={'old':[1,True,'abc',None,{'srepr':'Integer(0)','text':'0'}]};self.assertEqual(j.encode(v),v)
 def test_unknown_still_refuses(self):
  j=w.make_journal(Base,FunctionStub,SimpleNamespace(srepr=lambda v:None))()
  with self.assertRaises(TypeError):json.dumps(j.encode(object()))
 def test_exact_scientific_tail(self):self.assertEqual(tail(M/(OLD+'.py')),tail(M/(P+'.py')))
 def test_actual_binding_function_unchanged(self):
  def get(p):return next(n for n in function(p,'scientific_work').body if isinstance(n,ast.FunctionDef) and n.name=='bind')
  self.assertEqual(ast.dump(get(M/(OLD+'.py')),include_attributes=False),ast.dump(get(M/(P+'.py')),include_attributes=False))
 def test_no_original_equality_replay(self):
  calls={ast.unparse(n.func) for n in ast.walk(function(M/(P+'.py'),'scientific_work')) if isinstance(n,ast.Call)};self.assertNotIn("structural['exact_structure']",calls);self.assertNotIn('sp.integrate',calls)
 def test_gate_containment_order(self):
  f=ast.unparse(function(M/(P+'.py'),'main'));self.assertLess(f.index('verify_gate('),f.index('args.out.mkdir'));self.assertLess(f.index("ns['containment']()"),f.index('import sympy'))
 def test_prior_complete_tree_intact(self):
  m=json.loads((M/(P+'_inputs.json')).read_text());root=Path(m['priorRoot']);self.assertEqual({str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()},set(m['priorFiles']));self.assertEqual(len(m['priorFiles']),426)
  for name,r in m['priorFiles'].items():self.assertEqual(w.sha(root/name),r['sha256']);self.assertEqual((root/name).stat().st_size,r['bytes'])
 def test_partial_registry_is_preserved(self):
  m=json.loads((M/(P+'_inputs.json')).read_text());p=Path(m['priorRoot'])/'complete/LEFT-native-current-source.json'
  with self.assertRaises(json.JSONDecodeError):json.loads(p.read_text())
  self.assertIn('complete/LEFT-native-current-source.json',m['priorFiles'])
 def test_complete_prefix_count(self):
  m=json.loads((M/(P+'_inputs.json')).read_text());root=Path(m['priorRoot'])/'complete';a=json.loads((root/'artifact-index.json').read_text());self.assertEqual(len(a),6);self.assertEqual(json.loads((root/'right-native-restoration-return.json').read_text()),{'actualExpected':True,'originalExpected':True});self.assertFalse((root/'operation-index.json').exists())
 def test_all_source_and_saved_pins(self):
  m=json.loads((M/(P+'_inputs.json')).read_text());self.assertEqual(len(m['savedFiles']),344)
  for p,h in m['sourcePins'].items():self.assertEqual(w.sha(p),h)
  for r in m['savedFiles'].values():self.assertEqual(w.sha(r['path']),r['sha256'])
 def test_no_gate_or_extra_science(self):
  self.assertFalse((M/(P+'_gate.json')).exists());m=json.loads((M/(P+'_inputs.json')).read_text());self.assertIsNone(m['resources']['durationLimits']);self.assertEqual(m['resources']['memoryBytes'],4*1024**3);self.assertFalse(m['scope']['numericalIntegrals'])
if __name__=='__main__':unittest.main()
