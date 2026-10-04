"""Stdlib source/metadata and refusal tests; no native scientific execution."""
import ast,importlib.util,json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
HERE=Path(__file__).resolve().parent;PREFIX='S11c_d_first_order_receiving_blocks'
spec=importlib.util.spec_from_file_location('receiving_blocks',HERE/(PREFIX+'.py'));w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
class Tooling(unittest.TestCase):
    def test_no_science_on_import(self):self.assertNotIn('sympy',sys.modules)
    def test_compile_sources(self):
        for suffix in ['.py','_launch.py']:compile((HERE/(PREFIX+suffix)).read_text(),PREFIX+suffix,'exec')
    def test_strict_predicate(self):
        w.require(True,'true')
        for v in [False,None,1,0,'yes',[True]]:
            with self.assertRaises(ValueError):w.require(v,'refuse')
    def test_source_jet(self):
        self.assertEqual(w.jet_spec('u_2_t_d1d3'),{'channel':'u_2','timeOrder':1,'spatialOrders':[1,0,1]})
        self.assertEqual(w.jet_spec('grad_theta_2'),w.jet_spec('theta_d2'))
        for n in ['u_4_d1','theta_d4','pressure']:
            with self.assertRaises(ValueError):w.jet_spec(n)
    def test_exact_invocation(self):
        args=SimpleNamespace(out=Path('/tmp/receiving-out'),inputs=Path('/tmp/receiving-in'),gate=Path('/tmp/receiving-gate'))
        argv=[str(Path(w.__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)];g={'command':['guard']+argv,'outputDirectory':str(args.out)}
        w.verify_invocation(args,g,argv)
        with self.assertRaises(ValueError):w.verify_invocation(args,g,argv+['--retry'])
        with self.assertRaises(ValueError):w.verify_invocation(args,{**g,'outputDirectory':'/tmp/other'},argv)
    def test_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'operand';w.save(p,{'a':1})
            with self.assertRaises(FileExistsError):w.save(p,{'a':2})
            self.assertEqual(json.loads(p.read_text()),{'a':1})
    def test_guard_order(self):
        t=ast.parse((HERE/(PREFIX+'.py')).read_text());main=next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=='main');s=ast.unparse(main)
        self.assertLess(s.index('verify_gate('),s.index('args.out.mkdir'))
        self.assertLess(s.index("ns['containment']()"),s.index('import sympy'))
        self.assertNotIn('pickle.loads',s)
    def test_no_scientific_deadlines(self):
        t=ast.parse((HERE/(PREFIX+'.py')).read_text());calls=[ast.unparse(n.func) for n in ast.walk(t) if isinstance(n,ast.Call)]
        self.assertFalse(set(calls)&{'signal.alarm','signal.setitimer','subprocess.run','sp.integrate','sp.limit','sp.solve'})
    def test_manifest_all_original_bytes(self):
        m=json.loads((HERE/(PREFIX+'_inputs.json')).read_text())
        for r in m['savedFiles'].values():
            self.assertEqual(w.sha(r['path']),r['sha256']);self.assertEqual(Path(r['path']).stat().st_size,r['bytes'])
        self.assertEqual(sum(n.startswith('source-result/') for n in m['savedFiles']),783)
        self.assertEqual(sum(n.startswith('source-input/') for n in m['savedFiles']),127)
    def test_actual_numerical_leaf_inventory(self):
        m=json.loads((HERE/(PREFIX+'_inputs.json')).read_text());read=lambda n:json.loads(Path(m['savedFiles'][n]['path']).read_text())
        for side,alias in [('LEFT','source-input/incident/LEFT-raw-source-binding.json'),('RIGHT','receiving/RIGHT-raw-source-binding.json')]:
            v=read(alias);tree=ast.parse(v['nativeSource']['srepr'],mode='eval');names={ast.literal_eval(n.args[0]) for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='Symbol'}
            keys={a['text'] for a,b in v['map']['mappingPairs']};self.assertLessEqual(names,keys)
        v=read('source-input/consumer/chemical-amplitude-grade-split.json');tree=ast.parse(v['zeroGrade']['srepr'],mode='eval')
        for n in ast.walk(tree):
            if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='Symbol':w.jet_spec(ast.literal_eval(n.args[0]))
    def test_saved_pressure_route_completeness(self):
        m=json.loads((HERE/(PREFIX+'_inputs.json')).read_text());a=json.loads(Path(m['savedFiles']['source-result/first-order-pressure-assembly.json']['path']).read_text())
        self.assertEqual(len(a),10)
        for n in w.ROWS:
            rows=[v for v in a if v['row']==n and v['column']==0];self.assertEqual(len(rows),1)
            self.assertEqual({(v['face'],v['slot']) for v in rows[0]['pieces']},{(f,s) for f in ['plus','minus'] for s in ['pressure','normal']})
    def test_native_faces_and_domains_are_supplied(self):
        m=json.loads((HERE/(PREFIX+'_inputs.json')).read_text());raw=json.loads(Path(m['savedFiles']['receiving/LEFT-native-1-original.json']['path']).read_text());self.assertTrue(raw['passed']);a=raw['actual']['acoustic']
        self.assertEqual({f['ORIENTATION']['text'] for f in a['FACE_RECORDS']},{'-1','1'});self.assertEqual(set(a['MEMORY_KERNELS']),{'A','V','X'})
        self.assertIn('CLOSURE_RESIDUAL',a['FACE_RECORDS'][0]);self.assertIn('MECHANICAL_LOAD',a['FACE_RECORDS'][0])
    def test_method_and_build_are_separate(self):
        m=json.loads((HERE/(PREFIX+'_inputs.json')).read_text());r=json.loads(Path(m['methodRecord']).read_text())
        self.assertTrue(r['methodAssessed']);self.assertIn('MATCHED RECEIVING AND FACE METHOD',r['literalVerdict'])
        self.assertTrue(m['scope']['T1R1Pending']);self.assertTrue(m['scope']['G1K1Pending'])
        self.assertFalse((HERE/(PREFIX+'_gate.json')).exists())
    def test_fail_closed_scope(self):
        s=(HERE/(PREFIX+'.py')).read_text();self.assertIn('STOP_NO_AUTOMATIC_SCHUR_INVERSE',s)
        self.assertIn('SINGULAR_OR_UNKNOWN_RECEIVING_THRESHOLD',s)
        self.assertIn('FULL-offwave-source-correspondence',s)
        self.assertNotIn('rawends[\'RIGHT\'].inv',s)
if __name__=='__main__':unittest.main()
