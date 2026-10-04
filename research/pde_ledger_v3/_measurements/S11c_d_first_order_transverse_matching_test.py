"""Stdlib gate/source/input checks only; no scientific import or restoration."""
import ast,hashlib,importlib.util,json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
M=Path(__file__).resolve().parent;P='S11c_d_first_order_transverse_matching';spec=importlib.util.spec_from_file_location('matching',M/(P+'.py'));w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
class Tooling(unittest.TestCase):
 def test_no_science_on_import(self):self.assertNotIn('sympy',sys.modules)
 def test_compiles(self):
  for s in ['.py','_launch.py']:compile((M/(P+s)).read_text(),P+s,'exec')
 def test_strict_boolean(self):
  w.require(True,'yes')
  for v in [False,0,1,None,'true']:
   with self.assertRaises(ValueError):w.require(v,'no')
 def test_exclusive_evidence(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'record';w.save(p,{'original':1})
   with self.assertRaises(FileExistsError):w.save(p,{'overwrite':True})
   self.assertEqual(json.loads(p.read_text()),{'original':1})
 def test_gate_and_containment_order(self):
  t=ast.parse((M/(P+'.py')).read_text());f=ast.unparse(next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=='main'))
  self.assertLess(f.index('verify_gate('),f.index('args.out.mkdir'));self.assertLess(f.index("ns['containment']()"),f.index('import sympy'))
 def test_exact_argv(self):
  a=SimpleNamespace(out=Path('/tmp/match-o'),inputs=Path('/tmp/match-i'),gate=Path('/tmp/match-g'));v=[str(Path(w.__file__).resolve()),'--out',str(a.out),'--inputs',str(a.inputs),'--gate',str(a.gate)];g={'command':['guard']+v,'outputDirectory':str(a.out)};w.verify_invocation(a,g,v)
  with self.assertRaises(ValueError):w.verify_invocation(a,g,v+['--retry'])
 def test_all_saved_hashes(self):
  m=json.loads((M/(P+'_inputs.json')).read_text());self.assertEqual(len(m['savedFiles']),732)
  for r in m['savedFiles'].values():self.assertEqual(w.sha(r['path']),r['sha256']);self.assertEqual(Path(r['path']).stat().st_size,r['bytes'])
 def test_original_halfline_source_context(self):
  m=json.loads((M/(P+'_inputs.json')).read_text());read=lambda n:json.loads(Path(m['savedFiles'][n]['path']).read_text());end=read('receiving/end-matching-prerequisite.json')
  for k in range(2):self.assertEqual(end['savedForces'][k],read('source/full-local-force-column-'+str(k)+'.json'))
  self.assertEqual(read('source/extended-binding-context.json')['numeric']['L_W']['srepr'],'Integer(10)')
 def test_restored_scalar_count(self):
  m=json.loads((M/(P+'_inputs.json')).read_text());n=0
  for prefix,rows,cols in [('FULL-offwave-source-correspondence',5,5),('transverse-block',2,2),('transverse-to-scalars-coupling',3,2),('scalars-to-transverse-coupling',2,3),('full-five-row-matched-end',5,2),('end-RHS-minus-operator-sign',5,2)]:
   for i in range(rows):
    for j in range(cols):
     key='receiving/'+prefix+'-'+str(i)+'-'+str(j)+'-return.json';d=json.loads(Path(m['savedFiles'][key]['path']).read_text());self.assertEqual(d['cancelled'],{'text':'0','srepr':'Integer(0)'});n+=1
  self.assertEqual(n,61)
 def test_scope_no_hidden_current_or_numeric_solve(self):
  s=(M/(P+'.py')).read_text();t=ast.parse(s);calls={ast.unparse(n.func) for n in ast.walk(t) if isinstance(n,ast.Call)}
  self.assertFalse(calls&{'sp.integrate','sp.limit','sp.solve','signal.alarm','signal.setitimer','subprocess.run'})
  self.assertIn("'G1Computed':False",s);self.assertIn("'reflectionIntegralEvaluated':False",s)
  self.assertFalse((M/(P+'_gate.json')).exists())
if __name__=='__main__':unittest.main()
