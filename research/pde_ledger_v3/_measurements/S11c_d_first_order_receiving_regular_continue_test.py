"""Stdlib adapter and exact continuation-boundary checks; no symbolic restoration."""
import ast,hashlib,importlib.util,json,sys,unittest
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
M=Path(__file__).resolve().parent;PRE='S11c_d_first_order_receiving_regular_continue';OLD='S11c_d_first_order_receiving_regular'
spec=importlib.util.spec_from_file_location('receiving_continuation',M/(PRE+'.py'));w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
def read(p):return json.loads(Path(p).read_text())
def manifest():return read(M/(PRE+'_inputs.json'))
def fn(path,name):return next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
def nested(f,name):return next(n for n in f.body if isinstance(n,ast.FunctionDef) and n.name==name)
def dump(nodes):return ast.dump(ast.Module(body=nodes,type_ignores=[]),include_attributes=False)
class PolynomialStandIn:
    def __init__(self,coefficients):
        c=[Fraction(v) for v in coefficients]
        while c and c[0]==0:c.pop(0)
        self.coefficients=tuple(c);self.is_zero=not c
    def LC(self):return self.coefficients[0] if self.coefficients else Fraction(0)
    def monic(self):return PolynomialStandIn([v/self.LC() for v in self.coefficients])
    def as_expr(self):return self.coefficients
class Tests(unittest.TestCase):
    def engine(self):
        calls=[]
        def gcdex(*args):calls.append(args);return 'u','v','g'
        return SimpleNamespace(S=SimpleNamespace(One=Fraction(1),Zero=Fraction(0)),gcdex=gcdex),calls
    def test_actual_constant_zero_case_does_not_call_library(self):
        e,c=self.engine();self.assertEqual(w.exact_gcd_triplet(PolynomialStandIn([30]),PolynomialStandIn([]),'r',e),(Fraction(1,30),0,(Fraction(1),)));self.assertEqual(c,[])
    def test_both_zero_refuses(self):
        e,c=self.engine()
        with self.assertRaises(ValueError):w.exact_gcd_triplet(PolynomialStandIn([0]),PolynomialStandIn([]),'r',e)
        self.assertEqual(c,[])
    def test_one_zero_monic_bezout_exact_fraction_coefficients(self):
        for coefficients in [[-7],[Fraction(3,7),-2,5],[2,0,0],[-3,6]]:
            for zero_first in [False,True]:
                with self.subTest(coefficients=coefficients,zero_first=zero_first):
                    e,c=self.engine();f=PolynomialStandIn(coefficients);z=PolynomialStandIn([]);a,b=(z,f) if zero_first else (f,z)
                    u,v,g=w.exact_gcd_triplet(a,b,'r',e)
                    live=v if zero_first else u;dead=u if zero_first else v
                    self.assertEqual(dead,0);self.assertEqual(tuple(live*x for x in f.coefficients),g);self.assertEqual(g[0],1);self.assertEqual(c,[])
    def test_nonzero_pair_delegates_unchanged_arguments(self):
        e,c=self.engine();a,b=PolynomialStandIn([1,2]),PolynomialStandIn([3,0,4]);self.assertEqual(w.exact_gcd_triplet(a,b,'r',e),('u','v','g'));self.assertEqual(c,[(a.as_expr(),b.as_expr(),'r')])
    def test_no_scientific_import_or_original_functions(self):
        self.assertNotIn('sympy',sys.modules)
        for suffix in ['.py','_launch.py']:compile((M/(PRE+suffix)).read_text(),PRE+suffix,'exec')
        body=ast.unparse(fn(M/(PRE+'.py'),'scientific_work'))
        for expression in ['D5 * Feta','JL.subs','Feta =','sp.integrate','sp.solve','sp.limit','sp.nroots']:self.assertNotIn(expression,body)
    def test_original_source_and_all_files_are_immutable(self):
        m=manifest();r=read(m['priorRecord']);self.assertEqual(w.sha(m['baseWorker']),r['workerSha256'])
        root=Path(m['priorRoot']);self.assertEqual({str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()},set(m['priorFiles']))
        for name,v in m['priorFiles'].items():self.assertEqual(w.sha(root/name),v['sha256'],name);self.assertEqual((root/name).stat().st_size,v['bytes'])
        for path,h in m['sourcePins'].items():self.assertEqual(w.sha(path),h)
        for path,h in m['librarySourcePins'].items():self.assertEqual(w.sha(path),h)
    def test_exact_saved_failure_boundary_and_zero_prefix(self):
        m=manifest();c=Path(m['priorRoot'])/'complete';a=read(c/'artifact-index.json');self.assertEqual(len(a),2466)
        name='entry-0-joined-denominator-numerator-propagating-real-imag-input.json';self.assertEqual(list(a)[-1],name)
        failed=read(c/name);self.assertEqual(failed['real']['srepr'],'Integer(30)');self.assertEqual(failed['imaginary']['srepr'],'Integer(0)')
        zeros=[read(c/n) for n in a if n.endswith('-return.json') and 'cancelled' in read(c/n)];self.assertEqual(len(zeros),639)
        self.assertTrue(all(v['cancelled']['srepr']=='Integer(0)' for v in zeros));self.assertFalse((c/'operation-index.json').exists())
    def test_remaining_post_domain_tail_is_identical(self):
        old=fn(M/(OLD+'.py'),'scientific_work');new=fn(M/(PRE+'.py'),'scientific_work')
        def tail(f):
            i=next(i for i,n in enumerate(f.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value).startswith("domain_certificate('chart-K2'"));return f.body[i:-1]
        self.assertEqual(dump(tail(old)),dump(tail(new)))
    def test_post_gcd_checks_and_sturm_are_identical(self):
        old=nested(fn(M/(OLD+'.py'),'scientific_work'),'exclude_ray');new=nested(fn(M/(PRE+'.py'),'scientific_work'),'exclude_ray')
        def remainder(f):
            i=next(i for i,n in enumerate(f.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value).startswith("J.emit(name + '-gcd-return'"));return f.body[i:]
        self.assertEqual(dump(remainder(old)),dump(remainder(new)))
    def test_nonfirst_ray_expression_construction_unchanged(self):
        old=nested(fn(M/(OLD+'.py'),'scientific_work'),'exclude_ray');new=nested(fn(M/(PRE+'.py'),'scientific_work'),'exclude_ray')
        first_require=next(i for i,n in enumerate(old.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value).startswith('require('))
        self.assertIsInstance(new.body[0],ast.If);self.assertEqual(dump(old.body[:first_require]),dump(new.body[0].orelse))
    def test_domain_group_and_entry_remaining_calls_unchanged(self):
        old=fn(M/(OLD+'.py'),'scientific_work');new=fn(M/(PRE+'.py'),'scientific_work')
        a,b=nested(old,'domain_certificate'),nested(new,'domain_certificate')
        start=next(i for i,n in enumerate(a.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='group')
        self.assertEqual(dump(a.body[start:]),dump(b.body[1:]))
        a=next(n for n in old.body if isinstance(n,ast.For) and ast.unparse(n.iter)=='enumerate(C)')
        b=next(n for n in new.body if isinstance(n,ast.For) and ast.unparse(n.iter)=='enumerate(C)')
        boundary=next(i for i,n in enumerate(a.body) if isinstance(n,ast.For) and ast.unparse(n.iter)=='enumerate(bases)')
        self.assertEqual(dump(a.body[:boundary]),dump(b.body[0].orelse))
        self.assertEqual(dump(a.body[boundary:]),dump(b.body[1:]))
    def test_gcd_incomplete_frame_uses_saved_operands(self):
        source=(M/(PRE+'.py')).read_text()
        for s in ["saved['originalPolynomial']==polynomial", "saved['interval']==[0,upper]", "a.all_coeffs()==saved['realCoefficients']", "original['expression']==expr==saved['original']", "original['q']==q"]:self.assertIn(s,source)
        self.assertIn("'acceptedSubstitutionExecutionDependency'",(M/(OLD+'.py')).read_text())
    def test_actual_completed_native_inputs_unchanged(self):
        old=read(M/(OLD+'_inputs.json'));new=manifest()
        for key in ['savedFiles','scope','resources','methodRecord','methodSha256','helperSource']:self.assertEqual(old[key],new[key])
        self.assertEqual(len(new['savedFiles']),2025)
    def test_containment_precedes_import_and_creation(self):
        main=ast.unparse(fn(M/(PRE+'.py'),'main'));self.assertLess(main.index('verify_gate'),main.index('args.out.mkdir'));self.assertLess(main.index("ns['containment']()"),main.index('import sympy'))
        source=(M/(PRE+'.py')).read_text();self.assertNotIn('signal.alarm',source);self.assertNotIn('signal.setitimer',source)
    def test_exact_launcher_worker_and_resource_command(self):
        s=(M/(PRE+'_launch.py')).read_text();self.assertIn("PREFIX='"+PRE+"'",s);self.assertIn("'first_order_receiving_regular_continue'",s)
        self.assertIn('first-order-receiving-regular-continuation-01',s);self.assertIn("'READY_FOR_ONE_FIRST_ORDER_RECEIVING_REGULAR_CONTINUATION'",s)
        self.assertIn("'--memory-gib','4','--tasks-max','32'",s)
    def test_all_prior_originals_and_copies_are_posthashed(self):
        s=ast.unparse(fn(M/(PRE+'.py'),'main'));self.assertIn("manifest['priorFiles'].items()",s);self.assertIn("'prior-copy-index.json'",s);self.assertIn("manifest['librarySourcePins']",s)
if __name__=='__main__':unittest.main()
