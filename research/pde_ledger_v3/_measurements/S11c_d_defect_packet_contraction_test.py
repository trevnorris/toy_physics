"""Stdlib-only tooling/metadata tests; no native scientific arithmetic/restoration."""
import ast,hashlib,importlib.util,json,tempfile,unittest
from fractions import Fraction as F
from pathlib import Path
P=Path(__file__).with_name('S11c_d_defect_packet_contraction_lib.py')
spec=importlib.util.spec_from_file_location('new_contraction_helpers',P);C=importlib.util.module_from_spec(spec);spec.loader.exec_module(C)
M=P.parent;PREFIX='S11c_d_defect_packet_contraction'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
class Synthetic(unittest.TestCase):
    def test_strict_true(self):
        C.require(True,'ok')
        for v in (1,False,None,[],{},'yes'):
            with self.assertRaises(ValueError):C.require(v,'refuse')
    def test_units_add(self):self.assertEqual(C.add((1,2,3),('-1','1/2',0)),(F(0),F(5,2),F(3)))
    def test_units_scale(self):self.assertEqual(C.scale((1,2,3),F(1,2)),(F(1,2),F(1),F(3,2)))
    def test_unit_type_rejection(self):
        for v in ([True,0,0],[0.1,0,0],[1,2]):
            with self.assertRaises(ValueError):C.unit(v)
    def test_arithmetic_synthetic(self):
        n=ast.parse('def f(x,y):\n z=(x+y)*x/y\n return [z,x**2,-y]\n').body[0]
        a=C.Arithmetic({'x':F(2),'y':F(3)},F);self.assertEqual(a.statements(n),[F(10,3),F(4),F(-3)]);self.assertEqual(len(a.events),2)
    def test_arithmetic_missing_argument(self):
        with self.assertRaises(ValueError):C.Arithmetic({'x':F(1)},F).statements(ast.parse('def f(x,y):\n return [x]\n').body[0])
    def test_arithmetic_calls_refused(self):
        for text in ('f(x)','x.real','x[0]','x**(-1)','x**y','1.0'):
            with self.assertRaises(ValueError):C.Arithmetic({'x':F(1),'y':F(2)},F).expression(ast.parse(text,mode='eval').body)
    def test_arithmetic_statements_refused(self):
        for text in ('def f(x):\n import os\n return [x]','def f(x):\n x+=1\n return [x]','def f(x=1):\n return [x]','def f(x):\n return x'):
            with self.assertRaises(ValueError):C.Arithmetic({'x':F(1)},F).statements(ast.parse(text).body[0])
    def test_exact_source_fragment(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source.py';text='def toy(x):\n    return [x]\n';p.write_text(text);node=ast.parse(text).body[0];fragment=ast.get_source_segment(text,node)
            r={'source':str(p),'sourceSha256':sha(p),'line':1,'endLine':2,'text':fragment,'fragmentSha256':hashlib.sha256(fragment.encode()).hexdigest()}
            self.assertIsInstance(C.fragment(r),ast.FunctionDef);p.write_text(text+'# changed\n')
            with self.assertRaises(ValueError):C.fragment(r)
    def test_synthetic_scalar_constructor(self):
        class Stub:
            I='imag'
            Integer=staticmethod(lambda x:('integer',x))
            Rational=staticmethod(lambda p,q:('rational',p,q))
            Symbol=staticmethod(lambda n,**kw:('symbol',n,kw))
            Add=staticmethod(lambda *a:('add',a))
            Mul=staticmethod(lambda *a:('mul',a))
            Pow=staticmethod(lambda *a:('pow',a))
        self.assertEqual(C.restore_scalar(Stub,{'srepr':'Rational(-3, 5)'}),('rational',-3,5))
        for text in ('__import__("os")','Function("f")(Integer(1))','Float(1.0)','Rational(1,0)','Symbol("x", unknown=True)','Integer(True)','Integer(1).real'):
            with self.assertRaises((ValueError,TypeError)):C.restore_scalar(Stub,{'srepr':text})
    def test_toy_interval_central(self):self.assertEqual(C.finite_window(2,5,0)['lower'],F(-2))
    def test_toy_interval_wings(self):
        self.assertEqual(C.finite_window(2,5,6)['lower'],F(1));self.assertEqual(C.finite_window(2,5,-6)['upper'],F(-1))
    def test_toy_interval_empty(self):self.assertTrue(C.finite_window(2,5,8)['empty'])
    def test_toy_interval_boundary(self):self.assertTrue(C.finite_window(2,5,7)['degenerate'])
    def test_toy_orientation_membership(self):
        for k in range(-3,4):
            for l in range(-3,4):
                for t in range(-6,7):
                    before=C.rectangle_contains(2,5,k,l,t)
                    self.assertEqual(before,C.mapped_contains(2,5,k,l,k+t,'k'))
                    self.assertEqual(before,C.mapped_contains(2,5,k,l,l-t,'l'))
    def test_toy_wrong_variable_refuses(self):
        self.assertTrue(C.mapped_contains(2,5,-2,2,F(13,2),'l'));self.assertFalse(C.mapped_contains(2,5,-2,2,F(13,2),'k'))
    def test_toy_bad_domain_types(self):
        for args in ((2,1,0),(2,5,.5)):
            with self.assertRaises(ValueError):C.finite_window(*args)
class Metadata(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=read(M/(PREFIX+'_inputs.json'));cls.raw={k:read(v['path']) for k,v in cls.m['savedInputs'].items()}
    def test_all_saved_bytes(self):
        for r in self.m['savedInputs'].values():self.assertEqual((sha(r['path']),Path(r['path']).stat().st_size),(r['sha256'],r['bytes']))
    def test_source_pins(self):
        for p,h in self.m['sourcePins'].items():self.assertEqual(sha(p),h)
    def test_no_gate(self):self.assertFalse((M/(PREFIX+'_gate.json')).exists())
    def test_all_native_source_ancestry(self):
        for r in self.raw['source-contracts.json']['fragments'].values():
            C.fragment(r);old=self.raw['accepted-units/posthashes.json']['sources'][r['source']];self.assertEqual(old['expected'],r['sourceSha256']);self.assertTrue(old['intact'])
    def test_literal_saved_zero_return_schema(self):
        for a,d in self.raw.items():
            if a.endswith('-return.json') and a.startswith(('saved/field/','saved/factors/','saved/preflight/numeric-factor-','saved/inner/runtime-','saved/inner/new-runtime-')):
                self.assertEqual(d['cancelled'],{'text':'0','srepr':'Integer(0)'})
    def test_complete_row_selection(self):
        row=self.raw['saved/inventory/THETA_BALANCE-ordered-addresses.json'];sel=self.raw['saved/selected/pressure-addresses.json']['selected']
        self.assertEqual(len(row),2652);self.assertEqual([v for v in row if v['row']=='THETA_BALANCE' and v['jet']['channel']=='e_W'],sel);self.assertEqual(len(sel),544)
    def test_all_field_proof_operands(self):
        for k,f in self.raw['saved/pressure/fields.json'].items():self.assertEqual(f['field'],self.raw['saved/field/'+k+'-reconstruction-input.json']['left'])
    def test_all_template_proof_operands(self):
        for label,a in self.raw['saved/preflight/numeric-factor-adapters.json']['definitions'].items():
            u=self.raw['accepted-units/complete-template-'+label+'.json'];arg=self.raw['saved/preflight/numeric-factor-'+label+'-arguments.json'];proof=self.raw['saved/preflight/numeric-factor-'+label+'-input.json']
            self.assertEqual(u['actualAdapter'],a);self.assertEqual(u['actualArguments'],arg);self.assertEqual(proof['left'],a['mapped']);self.assertEqual(proof['right'],a['template']);self.assertIs(u['wholeDAdditionalResolvents'],False)
    def test_all_summand_actual_arguments(self):
        dims={r['addressId']:r for r in self.raw['accepted-units/complete-summand-dimensions.json']};adapters=self.raw['saved/preflight/numeric-factor-adapters.json']
        for i,a in enumerate(self.raw['saved/selected/pressure-addresses.json']['selected']):
            n=a['addressId'];inp=self.raw['accepted-units/summand-'+str(n)+'-input.json'];ret=self.raw['accepted-units/summand-'+str(n)+'-return.json'];label='-'.join((a['face'],a['slot'],a['component']))
            self.assertEqual(inp['address'],a);self.assertEqual(inp['sourceTransportInput']['address'],a);self.assertEqual(ret,dims[n]);self.assertEqual(inp['adapter'],adapters['definitions'][label]);self.assertEqual(inp['waveProof']['left'],a['waveMultiplier']);self.assertEqual(inp['normal'],a['normalMultiplier']);self.assertEqual(adapters['addressJoins'][i],{'addressId':n,'adapter':label})
    def test_original_physical_plan(self):
        p=self.raw['preflight/physical-plan.json'];self.assertEqual(p,self.raw['saved/preflight/physical-plan.json']);self.assertEqual(p,self.raw['accepted-units/parameter-quantity-origins.json']['physicalPlan']);self.assertEqual(p['centers'],[{'text':'-5/2','srepr':'Rational(-5, 2)'},{'text':'5/2','srepr':'Rational(5, 2)'}])
    def test_native_parameters_not_replayed(self):
        text=(M/(PREFIX+'.py')).read_text();tree=ast.parse(text)
        calls={ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n,ast.Call)}
        for forbidden in ('sp.integrate','sp.Integral','sp.limit','sp.lambdify','sp.N','mp.quad','kernel_components','profile','run_old'):self.assertNotIn(forbidden,calls)
    def test_no_science_at_import(self):
        for filename in (PREFIX+'.py',PREFIX+'_lib.py',PREFIX+'_launch.py'):
            for node in ast.parse((M/filename).read_text()).body:
                if isinstance(node,ast.Import):self.assertFalse(any(a.name.startswith(('sympy','mpmath','numpy','scipy')) for a in node.names))
                if isinstance(node,ast.ImportFrom):self.assertFalse((node.module or '').startswith(('sympy','mpmath','numpy','scipy')))
    def test_science_import_after_containment(self):
        text=(M/(PREFIX+'.py')).read_text();self.assertLess(text.index("enforced=ns['containment']()"),text.index('import sympy as sp'))
    def test_source_compile_only(self):
        for n in (PREFIX+'.py',PREFIX+'_lib.py',PREFIX+'_launch.py'):compile((M/n).read_text(),n,'exec')
    def test_declared_no_deadlines(self):
        self.assertIsNone(self.m['resources']['durationLimits']);self.assertEqual(self.m['resources']['memoryBytes'],4*1024**3);self.assertEqual(self.m['resources']['poolMemoryBytes'],16*1024**3);self.assertFalse(self.m['completedFunctionsReplayed'])
    def test_full_geometry_saved_no_reconstruction(self):self.assertEqual(len(self.m['geometryPlans']),8)
    def test_guard_and_hook_unchanged(self):
        self.assertIn('/var/projects/toy_physics/scripts/s11c_guarded_run.py',self.m['sourcePins']);text=(M/(PREFIX+'_launch.py')).read_text();self.assertIn('defect_packet_contraction',text);self.assertIn('01a0e01b-ef84-7192-817f-584cda5d339b',text)
if __name__=='__main__':unittest.main(verbosity=2)
