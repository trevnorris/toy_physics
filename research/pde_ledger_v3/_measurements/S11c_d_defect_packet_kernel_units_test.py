"""Synthetic unit-interpreter and static metadata tests; no native unit evaluation."""
import ast,hashlib,importlib.util,json,sys,tempfile,unittest
from pathlib import Path
from fractions import Fraction as F
M=Path(__file__).parent
spec=importlib.util.spec_from_file_location('new_kernel_unit_helpers',M/'S11c_d_defect_packet_kernel_units_lib.py');K=importlib.util.module_from_spec(spec);spec.loader.exec_module(K)
class SyntheticTests(unittest.TestCase):
    def check(self,text,expected,env=None,functions=None,params=None):self.assertEqual(K.Walk(env or {'z':[2,1,0],'w':[-1,0,0]},functions,params).dim(text),expected)
    def test_add(self):self.check('z+z',(2,1,0))
    def test_add_zero(self):self.check('z+0',(2,1,0))
    def test_zero(self):self.check('0',None)
    def test_product(self):self.check('z*w',(1,1,0))
    def test_quotient(self):self.check('z/w',(3,1,0))
    def test_integer_power(self):self.check('z**2',(4,2,0))
    def test_rational_power(self):self.check('z**(1/2)',(1,F(1,2),0))
    def test_rational_constructor_power(self):self.check('z**sp.Rational(1,2)',(1,F(1,2),0))
    def test_negative_power(self):self.check('w**(-2)',(2,0,0))
    def test_negation(self):self.check('-z',(2,1,0))
    def test_cancel_wrapper(self):self.check('sp.cancel(z)',(2,1,0))
    def test_i_pi(self):self.check('sp.I*sp.pi',K.ZERO)
    def test_function_body(self):self.check('f(w)',(-2,0,0),functions={'f':{'arguments':['u'],'body':'u*u'}})
    def test_function_arity(self):
        with self.assertRaises(ValueError):K.Walk({'w':[1,0,0]},{'f':{'arguments':['u'],'body':'u'}}).dim('f(w,w)')
    def test_function_signature(self):self.check('f(w)',(0,1,0),functions={'f':{'arguments':['u'],'inputUnits':[[-1,0,0]],'outputUnit':[0,1,0]}})
    def test_bad_function_argument(self):
        with self.assertRaises(ValueError):K.Walk({'w':[1,0,0]},{'f':{'arguments':['u'],'inputUnits':[[-1,0,0]],'outputUnit':[0,1,0]}}).dim('f(w)')
    def test_parameters(self):self.check("params['h']*w",(0,0,0),params={'h':[1,0,0]})
    def test_transcendental(self):self.check('sp.sinh(w/w)',K.ZERO)
    def test_literal_synthetic_kernel(self):
        fn=ast.parse('def f(z,w):\n a=z*w\n b=a/w\n return [a,b]\n').body[0]
        out,w=K.statements(fn,{'z':[2,1,0],'w':[-1,0,0]});self.assertEqual(out,[(1,1,0),(2,1,0)]);self.assertTrue(w.events)
    def test_partial_walk_preserved(self):
        fn=ast.parse('def f(z,w):\n a=z*w\n b=a+z\n return [b]\n').body[0];w=K.Walk({'z':[2,1,0],'w':[-1,0,0]})
        with self.assertRaises(ValueError):K.statements(fn,w.env,walk=w)
        self.assertTrue(any(v.get('refused') for v in w.events));self.assertTrue(any('assignment.a' in v['path'] for v in w.events))
    def test_unsupported_statement(self):
        fn=ast.parse('def f(z):\n z+=1\n return [z]\n').body[0]
        with self.assertRaises(ValueError):K.statements(fn,{'z':K.ZERO})
    def test_total_measure(self):
        r=K.total_unit([4,0,0],[-2,0,0],[0,1,0],[2,-1,1],[-2,0,0]);self.assertEqual(r['total'],(4,0,1))
    def test_old_zero_keeps_opposite_operands(self):
        args={'left':{'srepr':'Add(Integer(2),Integer(-2))'},'right':{'srepr':'Integer(0)'}}
        r=K.inherit_zero(args,{'cancelled':{'text':'0','srepr':'Integer(0)'}});self.assertFalse(r['functionCalled']);self.assertNotEqual(args['left'],args['right'])
    def test_refuse_false_old_zero(self):
        with self.assertRaises(ValueError):K.inherit_zero({'left':{'srepr':'Integer(1)'},'right':{'srepr':'Integer(0)'}},{'cancelled':{'text':'1','srepr':'Integer(1)'}})
    def test_changed_saved_left_refuses(self):
        with self.assertRaises(ValueError):K.inherit_zero({'left':{'srepr':'Integer(2)'},'right':{'srepr':'Integer(2)'}},{'cancelled':{'text':'0','srepr':'Integer(0)'},'left':{'srepr':'Integer(3)'}})
    def test_contract_mutation(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'sample.py';p.write_text('x=a*b\n');text='x=a*b';r={'source':str(p),'name':'x','line':1,'column':0,'kind':'Assign','text':text,'sha256':hashlib.sha256(text.encode()).hexdigest()}
            self.assertIsInstance(K.contract(r),ast.Assign);p.write_text('x=a/b\n')
            with self.assertRaises(ValueError):K.contract(r)
    def test_gate_source_only(self):
        t=ast.parse((M/'S11c_d_defect_packet_kernel_units.py').read_text());self.assertTrue(any(isinstance(n,ast.FunctionDef) and n.name=='verify_gate' for n in t.body))
    def test_no_scientific_import(self):self.assertNotIn('sympy',sys.modules);self.assertNotIn('mpmath',sys.modules)

for index,text in enumerate(['z+w','sp.exp(w)','1.5*z','True','unknown','0**0','1/0','z**w','z%w','evil(z)','sp.foo(z)',"params['missing']",'z.real','[z,w]']):
    def test(self,text=text):
        with self.assertRaises((ValueError,ZeroDivisionError)):K.Walk({'z':[2,1,0],'w':[-1,0,0]}).dim(text)
    setattr(SyntheticTests,'test_refusal_'+str(index),test)

class MetadataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=json.loads((M/'S11c_d_defect_packet_kernel_units_inputs.json').read_text())
    def test_all_pins(self):
        for p,h in self.m['sourcePins'].items():self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(),h,p)
    def test_all_saved_file_receipts(self):
        for a,r in self.m['savedInputs'].items():
            p=Path(r['path']);self.assertEqual(p.stat().st_size,r['bytes'],a);self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),r['sha256'],a)
    def test_literal_input_aliases(self):
        tree=ast.parse(Path(self.m['worker']).read_text())
        for n in ast.walk(tree):
            if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='get' and n.args and isinstance(n.args[0],ast.Constant):self.assertIn(n.args[0].value,self.m['savedInputs'])
    def test_all_source_contract_locations(self):
        c=json.loads(Path(self.m['savedInputs']['source-contracts.json']['path']).read_text())
        for r in c['fragments'].values():K.contract(r) # Source text only, never unit interpretation.
    def test_all_address_pairs_present(self):
        aliases=self.m['savedInputs'];src={a.removeprefix('source/address-').removesuffix('-input.json') for a in aliases if a.startswith('source/address-') and a.endswith('-input.json')}
        self.assertEqual(len(src),544)
        for a in src:
            self.assertIn('source/address-'+a+'-return.json',aliases);self.assertIn('preflight/address-'+a+'-wave-jet-input.json',aliases);self.assertIn('preflight/address-'+a+'-wave-jet-return.json',aliases)
    def test_full_numerical_payloads_not_used(self):self.assertTrue(all(Path(r['path']).suffix=='.json' for r in self.m['savedInputs'].values()))
    def test_no_ready_gate(self):self.assertFalse((M/'S11c_d_defect_packet_kernel_units_gate.json').exists())
    def test_no_opaque_or_scientific_execution(self):
        for name in ('S11c_d_defect_packet_kernel_units.py','S11c_d_defect_packet_kernel_units_lib.py'):
            tree=ast.parse((M/name).read_text())
            for n in ast.walk(tree):
                if isinstance(n,ast.Import):self.assertFalse(any(a.name in ('sympy','mpmath','pickle') for a in n.names))
    def test_resources(self):self.assertIsNone(self.m['resources']['durationLimits']);self.assertEqual(self.m['resources']['memoryBytes'],4*1024**3)
    def test_scope(self):self.assertTrue(self.m['scope']['kernelWaveMeasureTransportOnly']);self.assertFalse(self.m['completedFunctionsReplayed'])

if __name__=='__main__':unittest.main()
