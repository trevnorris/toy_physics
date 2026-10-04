"""Stdlib admission/selector tests. No scientific imports or saved-object decoding."""
import ast,importlib.util,json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('first_source',HERE/'S11c_d_first_order_source.py');w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
class Admission(unittest.TestCase):
    def test_no_science_on_import(self):self.assertNotIn('sympy',sys.modules)
    def test_native_mixed_jet(self):self.assertEqual(w.jet_spec('u_2_t_d1d3'),{'channel':'u_2','timeOrder':1,'spatialOrders':[1,0,1]})
    def test_native_theta_alias(self):self.assertEqual(w.jet_spec('grad_theta_2'),w.jet_spec('theta_d2'))
    def test_bad_jet_refuses(self):
        for name in ['u_4_d1','u_1_ddd','pressure','theta_d4']:
            with self.assertRaises(ValueError):w.jet_spec(name)
    def test_strict_predicate(self):
        for value in [1,'yes',[True],None,False]:
            with self.assertRaises(ValueError):w.require(value,'test')
    def test_invocation(self):
        args=SimpleNamespace(out=Path('/tmp/finite-source-test'),inputs=Path('/tmp/in.json'),gate=Path('/tmp/gate.json'))
        argv=[str(Path(w.__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
        gate={'command':['guard']+argv,'outputDirectory':str(args.out)}
        w.verify_invocation(args,gate,argv)
        with self.assertRaises(ValueError):w.verify_invocation(args,gate,argv+['--retry'])
        gate['outputDirectory']='/tmp/other'
        with self.assertRaises(ValueError):w.verify_invocation(args,gate,argv)
    def test_exclusive_persistence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'operand.json';w.save(p,{'value':1})
            with self.assertRaises(FileExistsError):w.save(p,{'value':2})
            self.assertEqual(json.loads(p.read_text()),{'value':1})
    def test_sources_compile(self):
        for name in ['S11c_d_first_order_source.py','S11c_d_first_order_source_launch.py']:
            compile((HERE/name).read_text(),name,'exec')
    def test_guard_before_scientific_import(self):
        source=(HERE/'S11c_d_first_order_source.py').read_text();main=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='main');text=ast.get_source_segment(source,main)
        self.assertLess(text.index("ns['containment']()"),text.index('import sympy'))
        self.assertLess(text.index('verify_gate('),text.index('args.out.mkdir('))
    def test_exact_saved_hashes_and_route_selectors(self):
        m=json.loads((HERE/'S11c_d_first_order_source_inputs.json').read_text())
        for r in m['savedFiles'].values():self.assertEqual(w.sha(r['path']),r['sha256'])
        r=json.loads(Path(m['savedFiles']['pressure/route-inspection.json']['path']).read_text())
        for route in r['routes']:
            source=json.loads(Path(route['sourcePath']).read_text());self.assertEqual(source[int(route['pointer'][1:])],route['completeAddress'])
    def test_profile_coverage(self):
        m=json.loads((HERE/'S11c_d_first_order_source_inputs.json').read_text());read=lambda n:json.loads(Path(m['savedFiles'][n]['path']).read_text())
        profiles=read('local/profile-jet-certificates.json');names={r['name'] for r in profiles}
        for face in ['plus','minus']:
            for grade in w.GRADES:
                source=read('sources/'+face+'-source-jets-'+grade+'.json')['source']['srepr']
                tree=ast.parse(source,mode='eval')
                for n in ast.walk(tree):
                    if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='Symbol' and n.args:
                        name=ast.literal_eval(n.args[0])
                        if name.startswith(('w1_profile','m1_profile')):self.assertIn(name,names)
    def test_no_ready_gate(self):self.assertFalse((HERE/'S11c_d_first_order_source_gate.json').exists())
if __name__=='__main__':unittest.main()
