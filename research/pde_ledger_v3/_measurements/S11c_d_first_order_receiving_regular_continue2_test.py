"""Task-local stdlib interface/provenance tests; no native scientific execution."""
import ast,hashlib,json,runpy,unittest
from pathlib import Path
M=Path(__file__).resolve().parent;PRE='S11c_d_first_order_receiving_regular_continue2';W=M/(PRE+'.py');OLD=M/'S11c_d_first_order_receiving_regular_continue.py'
ns=runpy.run_path(str(W),run_name='tooling_test_only');read=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();manifest=read(M/(PRE+'_inputs.json'));prior=Path(manifest['priorRoot']);C=prior/'complete';A=read(C/'artifact-index.json')
def dump(n):return ast.dump(n,include_attributes=False)
def fun(src,name):return next(n for n in ast.walk(ast.parse(src)) if isinstance(n,ast.FunctionDef) and n.name==name)
def after(body,fragment):
 return body[next(i for i,n in enumerate(body) if fragment in ast.unparse(n)):]
class FakePoly:
 def __init__(self,e):self.e=e
 def as_expr(self):return self.e.reconstruction
 def all_coeffs(self):return ['exact-complex-coefficient']
class FakeEngine:
 I='exact-i'
 def __init__(self,mutation=False,unsupported=False):self.reconstruction='grouped-mutant' if mutation else 'grouped-original';self.unsupported=unsupported;self.calls=[]
 def Poly(self,p,v,extension):
  self.calls.append((p,v,extension))
  if self.unsupported:raise ValueError('unsupported coefficient')
  return FakePoly(self)
 def expand(self,v):return {'original':'canonical','grouped-original':'canonical','grouped-mutant':'different'}.get(v,v)
class Journal:
 def __init__(self):self.records=[]
 def emit(self,n,v):self.records.append((n,v))
 def zero(self,n,l,r):
  self.records.append((n,{'left':l,'right':r}))
  if l!=r:raise ValueError('nonzero exact reconstruction')
class Tests(unittest.TestCase):
 def test_grouping_is_not_rejection(self):
  j=Journal();e=FakeEngine();ns['exact_polynomial_reconstruction']('x','original','q',j,e)
  self.assertEqual(e.calls,[('original','q','exact-i')]);self.assertFalse(j.records[1][1]['structurallyEqual']);self.assertEqual(j.records[-1][1],{'left':'canonical','right':'canonical'})
 def test_actual_mismatch_refuses_after_operands(self):
  j=Journal()
  with self.assertRaises(ValueError):ns['exact_polynomial_reconstruction']('x','original','q',j,FakeEngine(mutation=True))
  self.assertEqual(len(j.records),3);self.assertEqual(j.records[1][1]['reconstructed'],'grouped-mutant')
 def test_unsupported_coefficients_refuse_without_fallback(self):
  j=Journal()
  with self.assertRaises(ValueError):ns['exact_polynomial_reconstruction']('x','original','q',j,FakeEngine(unsupported=True))
  self.assertEqual(len(j.records),1)
 def test_full_original_tail_identical(self):
  a=after(fun(OLD.read_text(),'scientific_work').body,"domain_certificate('chart-K2'");b=after(fun(W.read_text(),'scientific_work').body,"domain_certificate('chart-K2'")
  self.assertEqual([dump(n) for n in a[:-1]],[dump(n) for n in b[:-1]])
 def test_gcd_and_sturm_logic_identical(self):
  a=after(fun(OLD.read_text(),'exclude_ray').body,'require(not (a.is_zero');b=after(fun(W.read_text(),'exclude_ray').body,'require(not (a.is_zero')
  self.assertEqual([dump(n) for n in a],[dump(n) for n in b]);self.assertEqual(dump(fun(OLD.read_text(),'exact_gcd_triplet')),dump(fun(W.read_text(),'exact_gcd_triplet')))
 def test_remaining_domain_calls_unchanged(self):
  a=fun(OLD.read_text(),'domain_certificate');b=fun(W.read_text(),'domain_certificate')
  self.assertEqual(dump(a.body[0].orelse[0]),dump(b.body[0].orelse[0]));self.assertEqual([dump(n) for n in a.body[0].orelse],[dump(n) for n in b.body[0].orelse])
  self.assertEqual(dump(a.body[-2].body[-1]),dump(b.body[-2].body[-1]))
 def test_complete_failed_tree_pins(self):
  self.assertEqual(set(manifest['priorFiles']),{str(p.relative_to(prior)) for p in prior.rglob('*') if p.is_file()})
  for n,v in manifest['priorFiles'].items():self.assertEqual((sha(prior/n),(prior/n).stat().st_size),(v['sha256'],v['bytes']))
 def test_twelve_complete_returns_and_boundary(self):
  returns=[n for n in A if n.endswith('-exact-Bezout-return.json')];self.assertEqual(len(returns),12);self.assertEqual(list(A)[-1],'entry-3-base-0-domain-rational-pair.json')
  for n in returns:self.assertEqual(read(C/n)['cancelled'],{'text':'0','srepr':'Integer(0)'})
 def test_all_original_saved_input_paths(self):
  self.assertEqual(len(manifest['savedFiles']),2025)
  for alias,v in manifest['savedFiles'].items():self.assertEqual(sha(C/'prior/complete/saved'/alias),v['sha256']);self.assertEqual(sha(v['path']),v['sha256'])
 def test_completed_domains_not_called(self):
  body=fun(W.read_text(),'scientific_work').body;restore=next(n for n in body if isinstance(n,ast.For) and ast.unparse(n.iter)=='range(3)');calls={ast.unparse(n.func) for n in ast.walk(restore) if isinstance(n,ast.Call)}
  self.assertFalse(calls & {'domain_certificate','exclude_ray','sp.gcdex','zero','J.zero','sp.together','sp.Poly'})
  resume=next(n for n in body if isinstance(n,ast.For) and ast.unparse(n.iter)=='range(3, len(C))');self.assertIn('completed',ast.unparse(resume.body[1]))
 def test_published_context_not_reconstructed(self):
  body=fun(W.read_text(),'scientific_work').body;assignment=next(n for n in body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='Cq' for t in n.targets))
  self.assertEqual(ast.unparse(assignment.value),"context['reconstitutedCq']");self.assertNotIn("C_s.subs",W.read_text())
 def test_original_library_sources_and_launcher(self):
  for p,h in {**manifest['sourcePins'],**manifest['librarySourcePins']}.items():self.assertEqual(sha(p),h)
  launcher=(M/(PRE+'_launch.py')).read_text();self.assertIn("'first_order_receiving_regular_continue2'",launcher);self.assertIn("'--memory-gib','4','--tasks-max','32'",launcher);self.assertNotIn('timeout=',launcher)
 def test_containment_precedes_scientific_import(self):
  main=ast.unparse(fun(W.read_text(),'main'));self.assertLess(main.index("ns['containment']()"),main.index('import sympy as sp'))
 def test_original_first_incomplete_operands(self):
  pair=read(C/'entry-3-base-0-domain-rational-pair.json');frame=read(C/'entry-3-full-domain-operands.json')
  self.assertEqual(pair['original'],frame['originalNegativePowers'][0]['transported']);self.assertEqual(pair['denominator'],{'text':'1','srepr':'Integer(1)'})
if __name__=='__main__':unittest.main()
