"""Stdlib packaging, isolation and persistence checks. No native science execution."""
import ast,importlib.util,json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
M=Path(__file__).resolve().parent;P='S11c_d_first_order_transverse_flux';spec=importlib.util.spec_from_file_location('flux',M/(P+'.py'));w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
class Tooling(unittest.TestCase):
 def test_no_science_on_import(self):self.assertNotIn('sympy',sys.modules);self.assertNotIn('numpy',sys.modules)
 def test_compiles(self):
  for suffix in ['.py','_launch.py']:compile((M/(P+suffix)).read_text(),P+suffix,'exec')
 def test_strict_boolean(self):
  w.require(True,'yes')
  for v in [False,0,1,None,'true']:
   with self.assertRaises(ValueError):w.require(v,'no')
 def test_exclusive_evidence(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'record';w.save(p,{'original':1})
   with self.assertRaises(FileExistsError):w.save(p,{'overwrite':True})
   self.assertEqual(json.loads(p.read_text()),{'original':1})
 def test_gate_before_output_science_after_containment(self):
  t=ast.parse((M/(P+'.py')).read_text());f=ast.unparse(next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=='main'))
  self.assertLess(f.index('verify_gate('),f.index('args.out.mkdir'));self.assertLess(f.index("ns['containment']()"),f.index('import sympy'))
 def test_exact_argv(self):
  a=SimpleNamespace(out=Path('/tmp/flux-o'),inputs=Path('/tmp/flux-i'),gate=Path('/tmp/flux-g'));v=[str(Path(w.__file__).resolve()),'--out',str(a.out),'--inputs',str(a.inputs),'--gate',str(a.gate)];g={'command':['guard']+v,'outputDirectory':str(a.out)};w.verify_invocation(a,g,v)
  with self.assertRaises(ValueError):w.verify_invocation(a,g,v+['--retry'])
 def test_all_original_hashes_and_census(self):
  m=json.loads((M/(P+'_inputs.json')).read_text());self.assertEqual(len(m['savedFiles']),344);self.assertEqual(len(m['opaqueAliases']),7)
  for r in m['savedFiles'].values():self.assertEqual(w.sha(r['path']),r['sha256']);self.assertEqual(Path(r['path']).stat().st_size,r['bytes'])
  for p,h in m['sourcePins'].items():self.assertEqual(w.sha(p),h)
 def test_inherited_argument_context(self):
  m=json.loads((M/(P+'_inputs.json')).read_text());read=lambda n:json.loads(Path(m['savedFiles'][n]['path']).read_text());a=read('matching/matched-transverse-amplitudes.json');e=read('receiving/end-matching-prerequisite.json')
  for key in ['B','DB','deltaP']:self.assertEqual(a[key],e[key])
  self.assertEqual(a['Bprime'],e['BprimeAtIncident']);self.assertFalse(a['physicalFluxNormalizationApplied'])
 def test_full_endpoint_returns_available(self):
  m=json.loads((M/(P+'_inputs.json')).read_text())
  for folder,prefix in [('receiving','full-five-row-matched-end'),('matching','right-field-finite-matching')]:
   for i in range(5):
    for j in range(2):
     a=f'{folder}/{prefix}-{i}-{j}-return.json';v=json.loads(Path(m['savedFiles'][a]['path']).read_text());self.assertEqual(v['cancelled'],{'text':'0','srepr':'Integer(0)'})
 def test_opaque_restoration_uses_only_existing_reader(self):
  s=(M/(P+'.py')).read_text();self.assertIn("inert(manifest['codecSource'],('SavedCodec','decode'),codec)",s);self.assertIn("('exact_structure',),structural",s);self.assertNotIn('pickle.loads(',s)
 def test_no_producer_solver_or_integral(self):
  t=ast.parse((M/(P+'.py')).read_text());calls={ast.unparse(n.func) for n in ast.walk(t) if isinstance(n,ast.Call)}
  self.assertFalse(calls&{'sp.integrate','sp.limit','sp.solve','sp.nsolve','signal.alarm','signal.setitimer','subprocess.run','current.construct','balance.construct','acoustic.construct','sp.inverse'})
  self.assertFalse((M/(P+'_gate.json')).exists())
 def test_no_numeric_request_or_scope_expansion(self):
  m=json.loads((M/(P+'_inputs.json')).read_text());self.assertFalse(m['scope']['receivingField']);self.assertFalse(m['scope']['numericalIntegrals']);self.assertTrue(m['scope']['recoveryParked']);self.assertIsNone(m['resources']['durationLimits']);self.assertEqual(m['resources']['memoryBytes'],4*1024**3)
if __name__=='__main__':unittest.main()
