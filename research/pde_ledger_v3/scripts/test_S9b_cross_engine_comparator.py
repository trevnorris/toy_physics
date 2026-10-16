#!/usr/bin/env python3
"""Synthetic native-stream controls; never read measured streams or engines.

Each mutation checks a computed product while holding the opposite operand
fixed. The companion disabling check makes that same assertion fail when its
comparison channel is disabled. No mutation assertion specifies a value.
"""
from dataclasses import replace
from contextlib import nullcontext
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import S9b_cross_engine_comparator as c


def symbol(x): return f'Symbol({x!r})'
def function(x,arg): return f'Function({x!r})({arg})'
def domain(value): return f'Tuple(Str("ATTACHED_DOMAIN"),{value})'
def py_gate(value,condition='true',outside='NOT_ESTABLISHED',outside_domain=None):
    outside_domain=outside_domain or f'Not({condition})'
    return f'Tuple(Tuple({domain(condition)},{value}),Tuple({domain(outside_domain)},Str({outside!r})))'
def wl_gate(value,condition='True',outside='NOT_ESTABLISHED',outside_domain=None):
    outside_domain=outside_domain or f'Not[{condition}]'
    return f'<|"OnDomain"->ConditionalExpression[{value},{condition}],"OutsideDomain"->ConditionalExpression["{outside}",{outside_domain}]|>'

def py_parts(value,positive='StrictGreaterThan(Symbol("rho_bulk"),Integer(0))',negative='LessThan(Symbol("rho_bulk"),Integer(0))'):
    return f'Tuple(Str("LOCAL_BULK_DENSITY_PARTS"),Tuple(Str("RHO_POSITIVE"),{positive}),Tuple(Str("RHO_NONPOSITIVE"),{negative}),Tuple(Str("OBJECT_ON_BOTH_PARTS"),{value}))'
def wl_parts(value,negative_value=None,labels=('rho > 0','rho <= 0')):
    return '<|"DensityParts"->{'+f'<|"Label"->"{labels[0]}","LabelCondition"->rho>0,"Object"->{value}|>,<|"Label"->"{labels[1]}","LabelCondition"->rho<=0,"Object"->{negative_value or value}|>'+'}|>'


def products(value):
    """Remove operands, sample values, diagnostics and raw serialization.

    Keep only mathematical products, comparison outcomes and child topology.
    No echoed operand (including a digest) can satisfy a mutation control.
    """
    if isinstance(value,list): return [products(x) for x in value]
    if not isinstance(value,dict): return value
    if value.get('samples'):
        value={k:v for k,v in value.items() if k not in ('residual','boolean_residual','zero_class')}
        # Counts/configuration text cannot rescue a disabled computation.
        value['samples']=[json.loads(x) for x in sorted({json.dumps(products(row),sort_keys=True) for row in value['samples']})]
    keys={'outcome','comparison','children','unpaired_metadata','residual','boolean_residual',
          'zero_class','samples','boolean_difference','same_exact_representation','coverage'}
    if 'outcome' in value or 'configuration' in value:
        return {k:products(v) for k,v in value.items() if k in keys}
    return {k:products(v) for k,v in value.items() if k not in ('py','wl','point','profile','py_text','wl_text','representation','unresolved','undefined','reason','precision_difference','method','reliability')}


class Controls(unittest.TestCase):
    def fixture(self,py,wl,*,names=c.NAME_TABLE,joins=('SYNTHETIC',),suffix_py='',suffix_wl=''):
        c._numeric_value.cache_clear()
        c.numerical_derivative.cache_clear()
        with tempfile.TemporaryDirectory() as tmp:
            p,w=Path(tmp)/'py.out',Path(tmp)/'wl.out'
            p.write_text(f'PY_S9B_SYNTHETIC: {py}\n'+suffix_py)
            w.write_text(f'WL_S9B_SYNTHETIC: {wl}\n'+suffix_wl)
            out,accounting=io.StringIO(),io.StringIO()
            c.run(p,w,out,accounting,names=names,joins=joins,progress=False)
        records=[json.loads(line) for line in out.getvalue().splitlines()]
        comparisons=[r['comparison'] for r in records if r['kind']=='comparison']
        return comparisons,[json.loads(line) for line in accounting.getvalue().splitlines()],records

    def mutation(self,py,wl,changed,*,channel='numeric_value'):
        before=self.fixture(py,wl); after=self.fixture(py,changed)
        self.assertEqual([r['py'] for r in before[2] if r['kind']=='operands'],
                         [r['py'] for r in after[2] if r['kind']=='operands'])
        self.assertNotEqual(products(before[0]),products(after[0]),
                            'mutation must move a computed comparison product')
        # An operand echo or parse change cannot rescue a disabled comparison.
        disabled={'numeric_value':{'unresolved':'disabled'},'symbolic_difference':c.sp.S.Zero,'text_comparison':None,'evaluate_profile_derivative':c.sp.S.Zero}
        original_translate=c.math_node
        def drop_arguments(node,bindings,engine,scope=None):
            target=node.args[0] if node.head=='Apply' else node
            if target.head in ('Symbol','FunctionName'):
                name=bindings.name(target.value,engine)
                if bindings.kinds.get(name)=='profile': return c.sp.Symbol(name)
            return original_translate(node,bindings,engine,scope)
        with patch.object(c,channel,**({'side_effect':drop_arguments} if channel=='math_node' else {'return_value':disabled[channel]})):
            with self.assertRaises(AssertionError):
                self.assertNotEqual(products(self.fixture(py,wl)[0]),products(self.fixture(py,changed)[0]))
        return before,after

    def test_one_sided_corruption(self):
        self.mutation('Add(Symbol("GM"),Integer(7))','GM+7','GM+11',channel='symbolic_difference')

    def test_form_change(self):
        self.mutation('Mul(Symbol("GM"),Symbol("b"))','GM*b','GM+b^2',channel='symbolic_difference')

    def test_sign_coefficient(self):
        self.mutation('Mul(Integer(3),Symbol("GM"))','3*GM','-5*GM',channel='symbolic_difference')

    def test_every_name_repoint(self):
        for row in c.NAME_TABLE:
            with self.subTest(binding=row.py):
                other=next(r for r in c.NAME_TABLE if r!=row and r.kind==row.kind)
                altered=tuple(replace(r,wl=other.wl) if r==row else replace(r,wl=row.wl) if r==other else r for r in c.NAME_TABLE)
                py=function(row.py,'Integer(7)') if row.kind=='profile' else f'Add(Integer(7),{symbol(row.py)})'
                wl=f'{row.wl}[7]' if row.kind=='profile' else f'7+{row.wl}'
                channel='numeric_value' if row.kind=='profile' else 'symbolic_difference'
                disabled={'unresolved':'disabled'} if row.kind=='profile' else c.sp.S.Zero
                if row.kind=='label':
                    py,wl=py_parts('Integer(7)'),wl_parts('7')
                    channel,disabled='text_comparison',None
                before=self.fixture(py,wl)[0]; after=self.fixture(py,wl,names=altered)[0]
                self.assertNotEqual(products(before),products(after),'repoint must move computed residual or label comparison')
                with patch.object(c,channel,return_value=disabled):
                    with self.assertRaises(AssertionError):
                        self.assertNotEqual(products(self.fixture(py,wl)[0]),products(self.fixture(py,wl,names=altered)[0]))

    def test_argument_stripped(self):
        self.mutation(function('s9b_V',symbol('s9b_r')),'velocityProfile[r]','velocityProfile',channel='math_node')

    def test_profile_constant(self):
        self.mutation(function('s9b_delta',symbol('s9b_r')),'deltaProfile[r]','3')

    def test_attached_domain(self):
        self.mutation(py_gate('Integer(7)','StrictGreaterThan(Symbol("GM"),Integer(0))'),
                      wl_gate('7','GM>0'),wl_gate('7','GM>=0'))

    def test_one_density_part_domain(self):
        py=py_parts(py_gate('Integer(7)','StrictGreaterThan(Symbol("GM"),Integer(0))'))
        original=wl_parts(wl_gate('7','GM>0'))
        changed=wl_parts(wl_gate('7','GM>0'),wl_gate('7','GM>=0'))
        self.mutation(py,original,changed)

    def test_density_labels_swapped(self):
        self.mutation(py_parts('Integer(7)'),wl_parts('7'),wl_parts('7',labels=('rho <= 0','rho > 0')),channel='text_comparison')

    def test_density_label_bindings_keep_both_emitted_strings(self):
        def labels(node):
            if isinstance(node,dict):
                if node.get('outcome')=='label': yield node
                else:
                    for value in node.values(): yield from labels(value)
            elif isinstance(node,list):
                for value in node: yield from labels(value)
        rows=[r for r in c.NAME_TABLE if r.kind=='label']
        self.assertEqual(len(rows),2)
        for swapped in (False,True):
            names=tuple(replace(r,wl=next(s.wl for s in rows if s!=r)) if swapped and r.kind=='label' else r for r in c.NAME_TABLE)
            result,_,records=self.fixture(py_parts('Integer(7)'),wl_parts('7'),names=names)
            self.assertEqual({(r['py_text'],r['wl_text']) for r in labels(result)},
                             {('RHO_POSITIVE','rho > 0'),('RHO_NONPOSITIVE','rho <= 0')})
            printed=next(r['name_table'] for r in records if r['kind']=='configuration')
            self.assertEqual([r for r in printed if r['kind']=='label'],[c.asdict(r) for r in names if r.kind=='label'])

    def test_brane_density_derivative_corruption_moves_formed_residual(self):
        # The reported corruption is invisible when rho=1+a/r. Use a native
        # derivative leaf, not an echoed profile formula or a production file.
        py='Derivative(Function("rho_br_profile")(Symbol("s9b_r")),Tuple(Symbol("s9b_r"),Integer(1)))'
        wl='Derivative[1][rhoBr][r]'
        changed='(1-rhoBr[r])/r'
        def check():
            before=self.fixture(py,wl); after=self.fixture(py,changed)
            self.assertEqual([r['py'] for r in before[2] if r['kind']=='operands'],
                             [r['py'] for r in after[2] if r['kind']=='operands'])
            select=lambda result:{r['configuration']:r for r in result[0][0]['samples']
                                  if r['profile'] in ('inverse','shifted','exponential')}
            a,b=select(before),select(after)
            self.assertTrue(a and a.keys()==b.keys())
            self.assertTrue(all(r['comparison']['coverage']=='residual formed' for rows in (a,b) for r in rows.values()))
            self.assertTrue(any(a[k]['residual']['computed_30']!=b[k]['residual']['computed_30'] for k in a),
                            'rho derivative corruption must move a formed residual under a generic admissible profile')
        check()
        original=c.profile_formulas
        def fixed_density(family,point):
            x,forms=original(family,point)
            return x,dict(forms,rho_br_profile=1+1/(13*x))
        with patch.object(c,'profile_formulas',side_effect=fixed_density):
            with self.assertRaises(AssertionError): check()
        with patch.object(c,'numeric_value',return_value={'unresolved':'disabled'}):
            with self.assertRaises(AssertionError): check()

    def test_every_live_profile_has_distinct_functional_forms(self):
        # Ratios at two distinct arguments remove pure rescaling. A fixed
        # opposite operand turns the sampled shapes into computed residuals.
        original=c.profile_formulas
        for binding in (r for r in c.NAME_TABLE if r.kind=='profile'):
            with self.subTest(profile=binding.py):
                py=f'Mul({function(binding.py,symbol("s9b_r"))},Pow({function(binding.py,symbol("b"))},Integer(-1)))'
                def check():
                    result=self.fixture(py,'7')[0][0]
                    selected=[r for r in result['samples'] if r['configuration'] in
                              ('p0/inverse','p0/shifted','p0/exponential')]
                    self.assertEqual(len(selected),3)
                    self.assertTrue(all(r['comparison']['coverage']=='residual formed' for r in selected))
                    self.assertEqual(len({r['residual']['computed_30'] for r in selected}),3,
                                     'each live profile must have different sampled shapes, not just amplitudes')
                check()
                with patch.object(c,'profile_formulas',side_effect=lambda family,point:original('inverse',point)):
                    with self.assertRaises(AssertionError): check()
                with patch.object(c,'numeric_value',return_value={'unresolved':'disabled'}):
                    with self.assertRaises(AssertionError): check()

    def test_marker_changed(self):
        self.mutation(py_gate('Integer(7)'),wl_gate('7'),wl_gate('7',outside='UNAVAILABLE'),channel='text_comparison')

    def test_marker_domain_changed(self):
        self.mutation(py_gate('Integer(7)'),wl_gate('7'),wl_gate('7',outside_domain='GM==0'))

    def test_boolean_flipped(self):
        self.mutation(py_gate('false'),wl_gate('False'),wl_gate('True'))

    def test_zero_parameter_stratum(self):
        before,after=self.mutation('Piecewise(ExprCondPair(Integer(7),Equality(Symbol("GM"),Integer(0))),ExprCondPair(Integer(3),true))',
                          'Piecewise[{{7,GM==0}},3]','Piecewise[{{11,GM==0}},3]')
        # Check actual zero-stratum probing, without specifying a residual.
        rows=after[0][0]['samples']
        zero=[r['residual'] for r in rows if r['point'].get('GM')=='0']
        other=[r['residual'] for r in rows if r['point'].get('GM')!='0']
        self.assertTrue(zero and other)
        self.assertNotEqual(zero[0],other[0])

    def test_zero_density_stratum(self):
        py='Piecewise(ExprCondPair(Integer(7),Equality(Add(Function("s9b_f")(Symbol("s9b_r")),Integer(1)),Integer(0))),ExprCondPair(Integer(3),true))'
        self.mutation(py,'Piecewise[{{7,1+fProfile[r]==0}},3]','Piecewise[{{11,1+fProfile[r]==0}},3]')

    def test_parsed_uncompared_leaf(self):
        py='Tuple(Integer(7),Mystery(Symbol("GM")))'
        before=self.fixture(py,'{7,Mystery[GM]}')[1]
        after=self.fixture(py,'{7,Mystery[GM,b]}')[1]
        self.assertNotEqual(before[0]['wl']['uncompared_leaves'],after[0]['wl']['uncompared_leaves'])
        with patch.object(c,'census',return_value=1000):
            with self.assertRaises(AssertionError):
                self.assertNotEqual(self.fixture(py,'{7,Mystery[GM]}')[1][0]['wl']['uncompared_leaves'],
                                    self.fixture(py,'{7,Mystery[GM,b]}')[1][0]['wl']['uncompared_leaves'])

    def test_derivative_order(self):
        py='Derivative(Function("s9b_delta")(Symbol("s9b_r")),Tuple(Symbol("s9b_r"),Integer(1)))'
        self.mutation(py,'Derivative[1][deltaProfile][r]','Derivative[2][deltaProfile][r]',channel='evaluate_profile_derivative')

    def test_integral_limits(self):
        py='Integral(Mul(Symbol("s9b_r"),Integer(2)),Tuple(Symbol("s9b_r"),Integer(1),Integer(3)))'
        self.mutation(py,'Inactive[Integrate][2*r,{r,1,3}]','Inactive[Integrate][2*r,{r,1,5}]')

    def test_binder_structure(self):
        py='Integral(Symbol("s9b_r"),Tuple(Symbol("s9b_r"),Integer(1),Integer(3)))'
        self.mutation(py,'Inactive[Integrate][r,{r,1,3}]','Inactive[Integrate][r,{b,1,3}]')

    def test_nested_sibling_removed(self):
        py='Tuple(Tuple(Str("left"),Integer(2)),Tuple(Str("right"),Tuple(Integer(3),Integer(5))))'
        wl='<|"left"->2,"right"->{3,5}|>'; changed='<|"left"->2,"right"->{3}|>'
        def missing(value):
            if isinstance(value,dict):
                return int(value.get('reason')=='nested sibling has no bound counterpart')+sum(missing(v) for v in value.values())
            return sum(missing(v) for v in value) if isinstance(value,list) else 0
        def check():
            self.assertNotEqual(missing(self.fixture(py,wl)[0]),missing(self.fixture(py,changed)[0]))
        check()
        with patch.object(c.Comparator,'absent',return_value={'outcome':'not_formed','reason':'disabled missing-sibling comparison'}):
            with self.assertRaises(AssertionError): check()

    def test_duplicate_join_and_names(self):
        with self.assertRaises(c.InputError): self.fixture('Integer(7)','7',joins=('SYNTHETIC','SYNTHETIC'))
        row=c.NAME_TABLE[0]
        for extra in (replace(row,py='another'),replace(row,wl='another')):
            with self.assertRaises(c.InputError): self.fixture('Integer(7)','7',names=(row,extra))

    def test_duplicate_stream_rejected_from_join(self):
        before=self.fixture('Integer(7)','7')[1]
        after=self.fixture('Integer(7)','7',suffix_wl='WL_S9B_SYNTHETIC: 11\n')[1]
        self.assertNotEqual(before,after)
        self.assertTrue(all(r['accounting']!='joined' for r in after))

    def test_unjoined_and_local_accounting(self):
        _,account,_=self.fixture('Integer(7)','7',suffix_py='PY_LOCAL_S9B_DETAIL: Integer(1)\nPY_S9B_UNKNOWN: Integer(2)\n',joins=('SYNTHETIC','MISSING'))
        self.assertEqual(len(account),4)
        self.assertEqual({r['accounting'] for r in account},{'joined','local','unjoined'})

    def test_native_negative_infinity(self):
        result=self.fixture('Integral(exp(Mul(Integer(-1),Pow(Symbol("s9b_r"),Integer(2)))),Tuple(Symbol("s9b_r"),-oo,oo))','Inactive[Integrate][Exp[-r^2],{r,-Infinity,Infinity}]')
        self.assertTrue(result[0])

    def test_duplicate_tuple_keys_preserved(self):
        n=c.py_parse('Tuple(Tuple(Str("x"),Integer(1)),Tuple(Str("x"),Integer(2)))')
        self.assertEqual(n.head,'Sequence')

    def test_compiled_quadrature_synthetic_reference(self):
        # Synthetic integrals only: finite, infinite, endpoint-singular,
        # complex, and cancellation. This checks the changed numerical path
        # against an independent symbolic integral, not a measured-stream target.
        x=c.sp.Symbol('x',positive=True)
        cases=[c.sp.Integral(x**2,(x,1,3)),
               c.sp.Integral(c.sp.exp(-x),(x,0,c.sp.oo)),
               c.sp.Integral(1/c.sp.sqrt(x),(x,0,1)),
               c.sp.Integral(c.sp.I*x,(x,0,1)),
               c.sp.Integral(1/((1+1/x)*x**2),(x,1,c.sp.oo))]
        for integral in cases:
            reference=complex(integral.doit().evalf(30))
            for digits in (18,30):
                self.assertLess(abs(complex(c.numerical_expression(integral,digits))-reference),1e-14)
        y=c.sp.Symbol('y',positive=True)
        nested=c.sp.Integral(x*c.sp.Integral(y,(y,0,x)),(x,0,1))
        self.assertLess(abs(complex(c.numerical_expression(nested,30))-complex(nested.doit().evalf(30))),1e-14)
        difference=cases[0]-c.sp.Integral(x**2+c.sp.Rational(1,10**12),(x,1,3))
        self.assertLess(abs(complex(c.numerical_expression(difference,30))-complex(difference.doit())),1e-25)

    def test_numerical_cache_keys_and_residual_reuse(self):
        x=c.sp.Symbol('x')
        a=c.sp.Integral(c.sp.exp(-x),(x,0,1))
        b=c.sp.Integral(c.sp.exp(-x),(x,0,2))
        c.numerical_integral.cache_clear()
        av=c.numerical_expression(a,30); bv=c.numerical_expression(b,30)
        before=c.numerical_integral.cache_info()
        c.numerical_expression(a-b,30)
        after=c.numerical_integral.cache_info()
        self.assertEqual(before.misses,after.misses)
        self.assertEqual(after.hits-before.hits,2)
        self.assertNotEqual(av,bv)
        c.numerical_expression(a,18)
        self.assertEqual(c.numerical_integral.cache_info().misses,after.misses+1)

    def test_leaf_cache_distinguishes_both_operands(self):
        # Both rows share one Comparator. The first primes the leaf cache;
        # the later row changes only one operand, in each direction in turn.
        for role in ('value','domain'):
            if role=='value':
                py,wl='Add(Symbol("GM"),Integer(7))','GM+7'
                changed_py,changed_wl='Add(Symbol("GM"),Integer(11))','GM+11'
            else:
                py=py_gate('Integer(7)','StrictGreaterThan(Symbol("GM"),Integer(0))')
                wl=wl_gate('7','GM>0')
                changed_py=py_gate('Integer(7)','GreaterThan(Symbol("GM"),Integer(0))')
                changed_wl=wl_gate('7','GM>=0')
            for side in ('py','wl'):
                with self.subTest(role=role,changed_operand=side):
                    def check():
                        def paired(p,w):
                            return self.fixture(py,wl,joins=('SYNTHETIC','ZZ_SECOND'),
                                suffix_py=f'PY_S9B_ZZ_SECOND: {p}\n',
                                suffix_wl=f'WL_S9B_ZZ_SECOND: {w}\n')
                        before=paired(py,wl)
                        after=paired(changed_py if side=='py' else py,changed_wl if side=='wl' else wl)
                        opposite='wl' if side=='py' else 'py'
                        self.assertEqual([r[opposite] for r in before[2] if r['kind']=='operands'],
                                         [r[opposite] for r in after[2] if r['kind']=='operands'])
                        self.assertEqual(products(before[0][0]),products(after[0][0]))
                        self.assertNotEqual(products(before[0][1]),products(after[0][1]),
                                            'cached result must not substitute for a different operand pair')
                    check()
                    narrowed=(lambda left,right,kind:(right,kind)) if side=='py' else (lambda left,right,kind:(left,kind))
                    with patch.object(c,'leaf_cache_key',side_effect=narrowed):
                        with self.assertRaises(AssertionError): check()
                    with patch.object(c,'symbolic_difference' if role=='value' else 'numeric_value',
                                      return_value=c.sp.S.Zero if role=='value' else {'unresolved':'disabled'}):
                        with self.assertRaises(AssertionError): check()

    def test_zero_profile_derivatives_at_distinct_arguments(self):
        py='Mul(Subs(Derivative(Function("s9b_xi")(Symbol("u")),Tuple(Symbol("u"),Integer(1))),Tuple(Symbol("u")),Tuple(Symbol("s9b_r"))),Subs(Derivative(Function("s9b_V")(Symbol("v")),Tuple(Symbol("v"),Integer(1))),Tuple(Symbol("v")),Tuple(Symbol("b"))))'
        wl='Derivative[1][xiProfile][r]*Derivative[1][velocityProfile][b]'
        comparison=self.fixture(py,wl)[0][0]
        rows=[r for r in comparison['samples'] if r['profile']=='optical_zero']
        self.assertEqual(len(rows),len(c.POINTS))
        self.assertTrue(all('exact' in r['py'] and 'exact' in r['wl'] for r in rows))

    def test_formatting_never_evaluates_numeric_integrals(self):
        x=c.sp.Symbol('x')
        integral=c.sp.Integral(c.sp.exp(-x),(x,0,1))
        held=c.sp.Add(c.sp.zoo,c.sp.Mul(2,integral,evaluate=False),evaluate=False)
        original_complex=c.sp.Expr.__complex__
        def reject_integral_complex(value):
            if value.has(c.sp.Integral):
                raise AssertionError('formatting evaluated an integral')
            return original_complex(value)
        with patch.object(c.sp.Expr,'__complex__',new=reject_integral_complex):
            printed=c.math_text(held)
            self.assertIn('Integral',printed)
            self.assertIn('zoo',printed)
            # This fixture detects the actual defect in the old printer.
            with self.assertRaises(AssertionError): str(held)
        py='Add(Mul(Integer(2),Integral(exp(Mul(Integer(-1),Symbol("s9b_r"))),Tuple(Symbol("s9b_r"),Integer(0),Integer(1)))),zoo)'
        wl='2*Inactive[Integrate][Exp[-r],{r,0,1}]+ComplexInfinity'
        # Native fixture passes through readers, extraction and every sample.
        # The patched printing-only method cannot suppress real arithmetic.
        with patch.object(c.sp.Expr,'__complex__',new=reject_integral_complex):
            rows=self.fixture(py,wl)[0][0]['samples']
        self.assertEqual(len(rows),len(c.POINTS)*len(c.PROFILE_NAMES))
        self.assertTrue(all('undefined' in r['py'] and 'undefined' in r['wl'] for r in rows))

    def test_full_configuration_set_preserved(self):
        py='Integral(Function("s9b_delta")(Symbol("s9b_r")),Tuple(Symbol("s9b_r"),Integer(2),Integer(3)))'
        wl='Inactive[Integrate][deltaProfile[r],{r,2,3}]'
        comparison=self.fixture(py,wl)[0][0]
        rows=comparison['samples']
        expected={f'p{i}/{profile}' for i in range(3) for profile in ('inverse','shifted','exponential','optical_zero','delta_zero','velocity_zero','tilt_zero','bulk_zero_at_point','bulk_negative_at_point','brane_zero_at_point','stiffness_zero_at_point')}
        self.assertEqual({r['configuration'] for r in rows},expected)
        self.assertTrue(all('py' in r and 'wl' in r and 'residual' in r for r in rows))
        # Disabling numerical evaluation must be observable in the computed
        # products, despite retaining every printed operand and configuration.
        with patch.object(c,'numeric_value',return_value={'unresolved':'disabled'}):
            self.assertNotEqual(products(comparison),products(self.fixture(py,wl)[0][0]))

    def test_negative_density_probe_computes(self):
        py='Piecewise(ExprCondPair(Integer(7),StrictLessThan(Add(Function("s9b_f")(Symbol("s9b_r")),Integer(1)),Integer(0))),ExprCondPair(Integer(3),true))'
        self.mutation(py,'Piecewise[{{7,1+fProfile[r]<0}},3]','Piecewise[{{11,1+fProfile[r]<0}},3]')

    def test_held_integral_derivative_is_not_artificially_undefined(self):
        py='Subs(Derivative(Integral(Pow(Add(Symbol("s9b_r"),Mul(Integer(-1),Symbol("b"))),Rational(-1,2)),Tuple(Symbol("s9b_r"),Symbol("b"),Add(Symbol("b"),Integer(1)))),Tuple(Symbol("b"),Integer(2))),Tuple(Symbol("b")),Tuple(Integer(2)))'
        # This derivative exists (the integral is invariant under translation).
        # Exact independence is recognized without splitting boundary terms.
        result=self.fixture(py,'0')[0][0]
        self.assertTrue(result['samples'])
        self.assertTrue(all('undefined' not in r['py'] for r in result['samples']))
        self.assertTrue(all(r['comparison'].get('coverage')!='numeric operands printed' for r in result['samples']))
        self.assertTrue(all(r['comparison']['coverage']=='residual formed' for r in result['samples']))
        with patch.object(c,'evaluate_substitution',side_effect=lambda v:v.doit(deep=False)):
            with self.assertRaises(AssertionError):
                self.assertTrue(all('undefined' not in r['py'] for r in self.fixture(py,'0')[0][0]['samples']))

    def test_interior_pole_never_prints_finite_quadrature(self):
        py='Integral(Pow(Function("rho_br_profile")(Symbol("s9b_r")),Integer(-1)),Tuple(Symbol("s9b_r"),Integer(2),Integer(5)))'
        wl='Inactive[Integrate][1/rhoBr[r],{r,2,5}]'
        def check():
            rows=self.fixture(py,wl)[0][0]['samples']
            singular=[r for r in rows if r['profile']=='brane_zero_at_point']
            self.assertEqual(len(singular),3)
            self.assertTrue(all('unresolved' in r['py'] or 'undefined' in r['py'] for r in singular))
            self.assertTrue(all(r['comparison']['coverage']!='residual formed' for r in singular))
        check()
        # Disable both safety channels to recreate the finite-value defect.
        with patch.object(c,'integral_screen'),patch.object(c,'numerical_expression',return_value=c.mp.mpf(1)):
            c.numerical_integral.cache_clear()
            with self.assertRaises(AssertionError): check()

    def test_unreliable_quadrature_is_unresolved(self):
        py='Integral(Function("s9b_delta")(Symbol("s9b_r")),Tuple(Symbol("s9b_r"),Integer(2),Integer(3)))'
        wl='Inactive[Integrate][deltaProfile[r],{r,2,3}]'
        with patch.object(c,'numerical_expression',side_effect=lambda expr,digits:c.mp.mpf(digits)):
            rows=self.fixture(py,wl)[0][0]['samples']
        self.assertTrue(all('unresolved' in r['py'] for r in rows if r['profile']!='optical_zero' and r['profile']!='delta_zero'))

    def test_quadrature_error_estimate_is_enforced(self):
        py='Integral(Function("s9b_delta")(Symbol("s9b_r")),Tuple(Symbol("s9b_r"),Integer(2),Integer(3)))'
        wl='Inactive[Integrate][deltaProfile[r],{r,2,3}]'
        def check():
            c.numerical_integral.cache_clear()
            rows=self.fixture(py,wl)[0][0]['samples']
            self.assertTrue(all('unresolved' in r['py'] for r in rows if r['profile'] not in ('optical_zero','delta_zero')))
        with patch.object(c.mp,'quad',return_value=(c.mp.mpf(1),c.mp.mpf(1))):
            check()
            with patch.object(c,'quadrature_consistent',return_value=True):
                with self.assertRaises(AssertionError): check()
        c.numerical_integral.cache_clear()

    def test_tiny_smooth_integral_forms_at_printed_precision(self):
        py='Integral(exp(Mul(Integer(-1),Symbol("s9b_r"))),Tuple(Symbol("s9b_r"),Integer(500),Integer(501)))'
        wl='Exp[-500]-Exp[-501]'
        changed=py.replace('Integer(501)','Integer(502)')
        def check():
            c.numerical_integral.cache_clear()
            before=self.fixture(py,wl); after=self.fixture(changed,wl)
            self.assertEqual([r['wl'] for r in before[2] if r['kind']=='operands'],
                             [r['wl'] for r in after[2] if r['kind']=='operands'])
            for result in (before,after):
                rows=result[0][0]['samples']
                self.assertTrue(rows)
                self.assertTrue(all(r['comparison']['coverage']=='residual formed' for r in rows))
                self.assertTrue(all(r['py'].get('reliable_digits')==12 for r in rows))
                self.assertTrue(all(r['residual']['reliability_established'] for r in rows))
            self.assertNotEqual(products(before[0]),products(after[0]))
            # An independent native antiderivative checks the actual printed
            # precision at this small operand scale, not an absolute floor.
            with c.mp.workdps(70):
                for row in before[0][0]['samples']:
                    scale=c.mp.mpf(row['residual']['operand_scale'])
                    printed=c.mp.mpc(row['py']['approx_12'])
                    reference=c.mp.mpc(*row['wl']['evaluations'][-1])
                    self.assertLess(abs(printed-reference),c.mp.mpf('1e-11')*scale)
        check()
        with patch.object(c,'quadrature_scale',return_value=c.mp.mpf(1)):
            with self.assertRaises(AssertionError): check()
        with patch.object(c,'numeric_value',return_value={'unresolved':'disabled'}):
            with self.assertRaises(AssertionError): check()
        c.numerical_integral.cache_clear()

    def test_record_gate_marker_and_domain(self):
        py='Tuple(Tuple(Tuple(Str("ATTACHED_DOMAIN"),true),Tuple(Str("condition"),true)),Tuple(Tuple(Str("ATTACHED_DOMAIN"),false),Str("NOT_ESTABLISHED")))'
        wl=wl_gate('<|"condition"->True|>')
        self.mutation(py,wl,wl_gate('<|"condition"->True|>',outside='UNAVAILABLE'),channel='text_comparison')
        self.mutation(py,wl,wl_gate('<|"condition"->True|>',outside_domain='GM==0'))
        result=self.fixture(py,wl)[0][0]['children']
        self.assertEqual(result['outside']['children']['value']['outcome'],'marker')
        self.assertEqual(result['outside']['children']['domain']['outcome'],'predicate')

    def test_empty_piece_distinct_from_unresolved(self):
        result=self.fixture('false','False')[0][0]
        self.assertTrue(any(s['coverage']=='empty' for s in result['strata']))
        self.assertTrue(any(s['coverage']=='sampled' for s in result['strata']))
        with patch.object(c,'stratum_coverage',return_value='unresolved coverage'):
            with self.assertRaises(AssertionError):
                self.assertTrue(any(s['coverage']=='empty' for s in self.fixture('false','False')[0][0]['strata']))

    def test_set_and_quantified_domains_are_sampled(self):
        condition='Equality(ConditionSet(Symbol("q"),StrictLessThan(Symbol("q"),Symbol("GM")),Interval(Integer(0),Integer(1))),EmptySet)'
        py=py_gate('Integer(7)',condition)
        original=wl_gate('7','Inactive[ForAll][r,0<=r && r<=1,r>=GM]')
        changed=wl_gate('7','Inactive[ForAll][r,0<=r && r<=1,r>GM]')
        self.mutation(py,original,changed)
        result=self.fixture(py,original)[0][0]['children']['on']['children']['domain']
        self.assertEqual(result['outcome'],'predicate')
        self.assertTrue(any('truth' in r['py'] and 'truth' in r['wl'] for r in result['samples']))

    def test_exists_domain_and_unresolved_quantifier(self):
        py=py_gate('Integer(7)','Unequality(ConditionSet(Symbol("q"),Equality(Pow(Symbol("q"),Integer(2)),Symbol("GM")),Reals),EmptySet)')
        self.mutation(py,wl_gate('7','Exists[r,r^2==GM]'),wl_gate('7','Exists[r,r^2==-GM]'))
        result=self.fixture(py,wl_gate('7','ForAll[r,Mystery[r]>GM]'))[0][0]
        domain=result['children']['on']['children']['domain']
        self.assertTrue(domain.get('samples'))
        self.assertTrue(all('unresolved' in r['wl'] for r in domain['samples']))

    def test_specific_evaluation_ablations_are_detected(self):
        with patch.object(c,'numeric_value',return_value={'unresolved':'disabled'}):
            for control in (self.test_profile_constant,self.test_integral_limits,self.test_zero_density_stratum,self.test_negative_density_probe_computes):
                with self.subTest(control=control.__name__),self.assertRaises(AssertionError): control()
        with patch.object(c,'evaluate_profile_derivative',return_value=c.sp.S.Zero):
            with self.assertRaises(AssertionError): self.test_derivative_order()
        with patch.object(c,'PROFILE_NAMES',tuple(p for p in c.PROFILE_NAMES if p not in ('bulk_zero_at_point','bulk_negative_at_point'))):
            for control in (self.test_zero_density_stratum,self.test_negative_density_probe_computes,self.test_full_configuration_set_preserved):
                with self.subTest(control=control.__name__),self.assertRaises(AssertionError): control()

    def test_unpaired_attached_domain_still_sampled(self):
        py='Tuple(Tuple(Str("value"),Integer(7)),Tuple(Str("Domain"),StrictGreaterThan(Symbol("GM"),Integer(0))))'
        result=self.fixture(py,'<|"value"->7|>')[0][0]['children']['Domain']
        rows=result['domain_comparison']['samples']
        self.assertTrue(rows)
        self.assertTrue(all('truth' in r['py'] and 'unresolved' in r['wl'] for r in rows))

    def test_completion_counts_actual_sample_outcomes(self):
        py='Integral(Function("rho_br_profile")(Symbol("s9b_r")),Tuple(Symbol("s9b_r"),Integer(2),Integer(5)))'
        records=self.fixture(py,'Inactive[Integrate][rhoBr[r],{r,2,5}]')[2]
        completion=next(r for r in records if r['kind']=='completion')
        samples=next(r for r in records if r['kind']=='comparison')['comparison']['samples']
        self.assertEqual(sum(completion['sample_outcomes'].values()),len(samples))
        self.assertEqual(set(completion['sample_outcomes']),{r['comparison']['coverage'] for r in samples})

    def test_undefined_predicate_boundary_is_unresolved(self):
        condition='StrictGreaterThan(Mul(Function("s9b_delta")(Symbol("s9b_r")),Pow(Function("s9b_V")(Symbol("s9b_r")),Integer(-1))),Integer(0))'
        result=self.fixture(py_gate('Integer(7)',condition),wl_gate('7','deltaProfile[r]/velocityProfile[r]>0'))[0][0]
        samples=result['children']['on']['children']['domain']['samples']
        boundary=[r for r in samples if r['profile']=='optical_zero']
        self.assertEqual(len(boundary),3)
        self.assertTrue(all('unresolved' in r['py'] and 'unresolved' in r['wl'] for r in boundary))
        self.assertTrue(all(r['comparison']['coverage']=='unresolved coverage' for r in boundary))

    def test_nonrational_quantifier_coefficients_stay_unresolved(self):
        condition='Equality(ConditionSet(Symbol("q"),StrictLessThan(Add(Mul(pi,Pow(Symbol("q"),Integer(6))),Mul(Integer(-1),Symbol("q")),Integer(-1)),Integer(0)),Interval(Integer(0),Integer(1))),EmptySet)'
        py=py_gate('Integer(7)',condition)
        wl=wl_gate('7','ForAll[r,0<=r && r<=1,Pi*r^6-r-1>=0]')
        with patch.object(c.Boolean,'as_set',side_effect=AssertionError('unsupported coefficient solver entered')):
            samples=self.fixture(py,wl)[0][0]['children']['on']['children']['domain']['samples']
        self.assertTrue(samples)
        self.assertTrue(all('unresolved' in r['py'] and 'unresolved' in r['wl'] for r in samples))
        self.assertTrue(all('non-rational coefficients' in r['py']['unresolved'] for r in samples))

    def test_equal_functionals_form_residual_at_operand_scale(self):
        # Integration by parts supplies two different forms of the same
        # synthetic integral. No measured-stream target is used.
        py='Integral(Mul(Symbol("s9b_r"),exp(Mul(Integer(-1),Symbol("s9b_r")))),Tuple(Symbol("s9b_r"),Integer(0),Integer(1)))'
        wl='Inactive[Integrate][Exp[-r],{r,0,1}]-Exp[-1]'
        def check():
            rows=self.fixture(py,wl)[0][0]['samples']
            self.assertTrue(rows)
            self.assertTrue(all(r['comparison'].get('coverage')=='residual formed' for r in rows))
            self.assertTrue(all(r['residual']['reliability_established'] for r in rows))
            self.assertTrue(all('computed_18' in r['residual'] and 'computed_30' in r['residual'] for r in rows))
        check()
        with patch.object(c,'numerical_residual',return_value={'unresolved':'disabled'}):
            with self.assertRaises(AssertionError): check()
        with patch.object(c,'residual_scale',side_effect=lambda a,b:abs(a[-1]-b[-1])):
            with self.assertRaises(AssertionError): check()
        # A reliability failure retains the computed values and formed outcome.
        with patch.object(c,'residual_scale',return_value=c.mp.mpf(0)):
            rows=self.fixture(py,wl)[0][0]['samples']
            self.assertTrue(all('computed_30' in r['residual'] and r['comparison']['coverage']=='residual formed' for r in rows))

    def test_symbol_exchange_moves_numeric_comparison(self):
        py='Integral(Mul(Symbol("b"),Symbol("s9b_r")),Tuple(Symbol("s9b_r"),Integer(1),Integer(2)))'
        def check():
            a=self.fixture(py,'Inactive[Integrate][b*r,{r,1,2}]')[0][0]['samples']
            b=self.fixture(py,'Inactive[Integrate][c0*r,{r,1,2}]')[0][0]['samples']
            self.assertEqual(len(a),33)
            self.assertTrue(all(products(x['comparison'])!=products(y['comparison']) for x,y in zip(a,b)))
        check()
        with patch.object(c,'POINTS',tuple(dict(p,c_0=p['b']) for p in c.POINTS)):
            with self.assertRaises(AssertionError): check()
        # Includes parameter strata, where a probe could collide with another
        # symbol even when the three generic points are distinct.
        expressions=(sum(c.sp.Symbol(n) for n in ('GM','Phi','K','n','s','rho_bulk','s9b_fractional_bulk_change')),)
        for _,point,_ in c.configurations(expressions):
            self.assertEqual(len(set(point.values())),len(point))

    def test_attempted_unformed_leaves_are_not_compared(self):
        py='Integral(Symbol("s9b_r"),Tuple(Symbol("s9b_r"),Integer(0),Symbol("unassigned")))'
        wl='Inactive[Integrate][r,{r,0,unassigned}]'
        def check():
            account=self.fixture(py,wl)[1][0]
            self.assertFalse(account['comparison_formed'])
            for engine in ('py','wl'):
                counts=account[engine]
                self.assertEqual(counts['compared_leaves'],0)
                self.assertGreater(counts['attempted_unformed_leaves'],0)
                self.assertEqual(counts['parsed_leaves'],counts['compared_leaves']+counts['attempted_unformed_leaves']+counts['unattempted_leaves'])
        check()
        with patch.object(c,'sample_formed',return_value=True):
            with self.assertRaises(AssertionError): check()
        account=self.fixture('Tuple(Integer(7),'+py+')','{7,'+wl+'}')[1][0]
        for engine in ('py','wl'):
            self.assertGreater(account[engine]['compared_leaves'],0)
            self.assertGreater(account[engine]['attempted_unformed_leaves'],0)

    def test_integration_constants_remain_engine_qualified(self):
        py='Symbol("s9b_velocity_integration_constant")'
        wl='velocityIntegrationConstant'
        result=self.fixture(py,wl)[0][0]
        self.assertIn('py::s9b_velocity_integration_constant',result['residual'])
        self.assertIn('wl::velocityIntegrationConstant',result['residual'])
        altered=c.NAME_TABLE+(c.NameRow('s9b_velocity_integration_constant','velocityIntegrationConstant','unproved binding',470,524,'parameter'),)
        self.assertNotEqual(products(result),products(self.fixture(py,wl,names=altered)[0][0]))

    def test_comparison_kinds_do_not_count_markers_or_labels_as_values(self):
        fixtures=[
            ('Str("NOT_ESTABLISHED")','"NOT_ESTABLISHED"',{'marker'}),
            (py_parts('Mystery(Integer(7))'),wl_parts('Mystery[7]'),{'label','domain'}),
            ('Tuple(Tuple(Str("Domain"),true))','<|"Domain"->True|>',{'domain'}),
        ]
        def check():
            for py,wl,kinds in fixtures:
                _,account,records=self.fixture(py,wl)
                row=account[0]
                self.assertEqual({k for k,v in row['comparison_kinds_formed'].items() if v},kinds)
                self.assertFalse(row['value_comparison_formed'])
                self.assertTrue(row['comparison_formed'])
                for engine in ('py','wl'):
                    self.assertEqual(row[engine]['compared_leaves_by_kind']['value'],0)
                    self.assertGreater(row[engine]['compared_leaves'],0)
                completion=next(r for r in records if r['kind']=='completion')
                self.assertEqual(completion['rows_without_formed_value_comparison'],1)
                self.assertEqual(completion['rows_without_formed_value_comparison_names'],['SYNTHETIC'])
                self.assertFalse(completion.get('rows_with_formed_value_comparison',0))
                self.assertFalse(completion['sample_outcomes_by_kind']['value'])
        check()
        with patch.object(c,'comparison_kind',return_value='value'):
            with self.assertRaises(AssertionError): check()
        # One-sided mutation changes the computed kind/accounting, not just
        # the marker text; the unchanged native operand is checked explicitly.
        py,wl,_=fixtures[0]
        before=self.fixture(py,wl); after=self.fixture(py,'7')
        self.assertEqual([r['py'] for r in before[2] if r['kind']=='operands'],
                         [r['py'] for r in after[2] if r['kind']=='operands'])
        self.assertNotEqual(before[1][0]['comparison_kinds_formed'],after[1][0]['comparison_kinds_formed'])
        with patch.object(c,'text_comparison',return_value=None):
            with self.assertRaises(AssertionError):
                self.assertNotEqual(self.fixture(py,wl)[1][0]['comparison_kinds_formed'],
                                    self.fixture(py,'7')[1][0]['comparison_kinds_formed'])

    def test_scalar_and_boolean_values_have_separate_value_accounting(self):
        for py,wl in [('Integer(7)','7'),('true','True'),
                      (py_gate('Integer(7)'),wl_gate('7'))]:
            _,account,records=self.fixture(py,wl)
            self.assertTrue(account[0]['value_comparison_formed'])
            self.assertTrue(account[0]['comparison_kinds_formed']['value'])
            self.assertTrue(all(account[0][e]['compared_leaves_by_kind']['value']>0 for e in ('py','wl')))
            completion=next(r for r in records if r['kind']=='completion')
            self.assertEqual(completion['rows_without_formed_value_comparison'],0)

    def test_integrable_singularity_forms_with_nonrational_profile(self):
        py='Integral(Mul(Function("s9b_delta")(Symbol("s9b_r")),Pow(Add(Symbol("s9b_r"),Integer(-1)),Rational(-1,2))),Tuple(Symbol("s9b_r"),Integer(1),Integer(3)))'
        wl='Inactive[Integrate][deltaProfile[r]/Sqrt[r-1],{r,1,3}]'
        def check():
            c.numerical_integral.cache_clear()
            rows=self.fixture(py,wl)[0][0]['samples']
            for family in ('inverse','shifted','exponential'):
                selected=[r for r in rows if r['profile']==family]
                self.assertEqual(len(selected),len(c.POINTS))
                self.assertTrue(all(r['comparison']['coverage']=='residual formed' for r in selected))
                self.assertTrue(all(r[e].get('reliable_digits')==12 for r in selected for e in ('py','wl')))
        check()
        self.mutation(py,wl,'Inactive[Integrate][deltaProfile[r]^2/Sqrt[r-1],{r,1,3}]')
        with patch.object(c,'endpoint_convergence',side_effect=c.UnresolvedEvaluation('disabled endpoint certificate')):
            with self.assertRaises(AssertionError): check()
        c.numerical_integral.cache_clear()

    def test_nonintegrable_exponential_endpoint_is_not_printed_as_numeric(self):
        py='Integral(Mul(exp(Mul(Integer(-1),Symbol("s9b_r"))),Pow(Add(Symbol("s9b_r"),Integer(-1)),Integer(-1))),Tuple(Symbol("s9b_r"),Integer(1),Integer(3)))'
        wl='Inactive[Integrate][Exp[-r]/(r-1),{r,1,3}]'
        def check():
            c.numerical_integral.cache_clear()
            rows=self.fixture(py,wl)[0][0]['samples']
            self.assertTrue(all(r['comparison']['coverage']=='unresolved coverage' for r in rows))
            self.assertTrue(all('endpoint convergence' in r[e]['unresolved'] for r in rows for e in ('py','wl')))
        with patch.object(c.mp,'quad',return_value=(c.mp.mpf(1),c.mp.mpf(0))):
            check()
            with patch.object(c,'endpoint_convergence'):
                with self.assertRaises(AssertionError): check()
        c.numerical_integral.cache_clear()

    def test_complete_variable_limit_derivative_forms(self):
        # Its second derivative exists although differentiating the lower
        # endpoint and the integrand separately produces singular terms.
        py='Subs(Derivative(Integral(Pow(Add(Symbol("s9b_r"),Mul(Integer(-1),Symbol("q"))),Rational(1,2)),Tuple(Symbol("s9b_r"),Symbol("q"),Symbol("r_E"))),Tuple(Symbol("q"),Integer(2))),Tuple(Symbol("q")),Tuple(Symbol("b")))'
        wl='1/(2*Sqrt[rE-b])'
        def check():
            rows=self.fixture(py,wl)[0][0]['samples']
            self.assertTrue(rows)
            self.assertTrue(all(r['comparison']['coverage']=='residual formed' for r in rows))
            self.assertTrue(all(r['residual']['reliability_established'] for r in rows))
            # Accuracy is checked against the independent synthetic native
            # formula, not a measured-stream target or an operand echo.
            for row in rows:
                residual=row['residual']
                with c.mp.workdps(60):
                    difference=c.mp.mpc(*row['py']['evaluations'][-1])-c.mp.mpc(*row['wl']['evaluations'][-1])
                    self.assertLessEqual(abs(difference),c.mp.mpf('1e-12')*c.mp.mpf(residual['operand_scale']))
        check()
        self.mutation(py,wl,'1/(3*Sqrt[rE-b])')
        with patch.object(c,'numerical_derivative',side_effect=c.UnresolvedEvaluation('disabled complete derivative')):
            with self.assertRaises(AssertionError): check()
        with patch.object(c,'derivative_consistent',return_value=False):
            with self.assertRaises(AssertionError): check()

    def test_variable_limit_derivative_mutation_moves_comparison(self):
        py='Subs(Derivative(Integral(exp(Mul(Integer(-1),Symbol("s9b_r"))),Tuple(Symbol("s9b_r"),Integer(1),Symbol("q"))),Tuple(Symbol("q"),Integer(1))),Tuple(Symbol("q")),Tuple(Symbol("b")))'
        changed=py.replace('Integer(1),Symbol("q")))','Integer(1),Mul(Integer(2),Symbol("q"))))')
        def check():
            a=self.fixture(py,'Exp[-b]'); b=self.fixture(changed,'Exp[-b]')
            self.assertEqual([r['wl'] for r in a[2] if r['kind']=='operands'],
                             [r['wl'] for r in b[2] if r['kind']=='operands'])
            self.assertNotEqual(products(a[0]),products(b[0]))
            self.assertTrue(all(r['comparison']['coverage']=='residual formed' for result in (a,b) for r in result[0][0]['samples']))
        check()
        with patch.object(c,'numerical_derivative',side_effect=c.UnresolvedEvaluation('disabled complete derivative')):
            with self.assertRaises(AssertionError): check()

    def test_pole_free_trigonometric_integral_forms(self):
        py='Integral(Mul(tan(Symbol("s9b_r")),Pow(Add(Integer(1),sec(Symbol("s9b_r"))),Integer(-1))),Tuple(Symbol("s9b_r"),Integer(0),Rational(1,2)))'
        wl='Inactive[Integrate][Tan[r]/(1+Sec[r]),{r,0,1/2}]'
        def check():
            c.numerical_integral.cache_clear()
            rows=self.fixture(py,wl)[0][0]['samples']
            self.assertTrue(rows)
            self.assertTrue(all(r['comparison']['coverage']=='residual formed' for r in rows))
            self.assertTrue(all(r['residual']['reliability_established'] for r in rows))
        check()
        self.mutation(py,wl,'Inactive[Integrate][Tan[r]/(1+Sec[r]),{r,0,2/3}]')
        with patch.object(c,'trig_poles',side_effect=c.UnresolvedEvaluation('disabled pole locations')):
            with self.assertRaises(AssertionError): check()
        c.numerical_integral.cache_clear()

    def test_trigonometric_pole_is_not_admitted_by_range_screen(self):
        py='Integral(sec(Symbol("s9b_r")),Tuple(Symbol("s9b_r"),Integer(0),Integer(2)))'
        wl='Inactive[Integrate][Sec[r],{r,0,2}]'
        def check():
            c.numerical_integral.cache_clear()
            rows=self.fixture(py,wl)[0][0]['samples']
            self.assertTrue(all(r['comparison']['coverage']=='unresolved coverage' for r in rows))
        # Isolate the safety check from quadrature's incidental ability to
        # notice a pole. No successful quadrature nodes certify convergence.
        with patch.object(c.mp,'quad',return_value=(c.mp.mpf(1),c.mp.mpf(0))):
            check()
            with patch.object(c,'trig_poles',return_value=set()):
                with self.assertRaises(AssertionError): check()
        c.numerical_integral.cache_clear()

    def test_every_function_binding_repoint(self):
        # Independent native fixtures, not synthesized from table targets.
        # Baseline comparisons detect an already-corrupted mapping (including
        # Sqrt -> cube root); mutations test whether the dispatch is actually used.
        cases={
            'Add':('Integer(8)','Plus[3,5]'), 'Mul':('Integer(15)','Times[3,5]'),
            'Pow':('Integer(8)','Power[2,3]'),
            'Equality':('true','Equal[3,3]'), 'Unequality':('false','Unequal[3,3]'),
            'StrictGreaterThan':('false','Greater[3,3]'), 'StrictLessThan':('false','Less[3,3]'),
            'GreaterThan':('true','GreaterEqual[3,3]'), 'LessThan':('true','LessEqual[3,3]'),
            'Sequence':('true','Element[List[2,3],Reals]'),
            'Integral':('Integral(Symbol("s9b_r"),Tuple(Symbol("s9b_r"),Integer(1),Integer(3)))','Integrate[r,{r,1,3}]'),
            'Limit':('Limit(Symbol("s9b_r"),Symbol("s9b_r"),Integer(2))','Inactive[Limit][r,r->2]'),
            'Contains':('true','Element[2,Reals]'),
            'ForAll':('true','ForAll[r,r^2>=0]'), 'Exists':('true','Exists[r,r==2]'),
            'And':('false','And[True,False]'), 'Or':('true','Or[True,False]'), 'Not':('false','Not[True]'),
            'Piecewise':('Piecewise(ExprCondPair(Integer(3),true))','Piecewise[{{3,True}},5]'),
            'Rational':('Rational(3,7)','Rational[3,7]'),
            'Min':('Integer(3)','Min[3,5]'), 'Max':('Integer(5)','Max[3,5]'),
            'Abs':('Integer(3)','Abs[-3]'), 'sin':('Rational(1,2)','Sin[Pi/6]'),
            'cos':('Rational(1,2)','Cos[Pi/3]'), 'tan':('Integer(1)','Tan[Pi/4]'),
            'sec':('Integer(2)','Sec[Pi/3]'), 'acos':('Mul(pi,Rational(1,3))','ArcCos[1/2]'),
            'log':('Integer(1)','Log[E]'), 'exp':('Integer(1)','Exp[0]'),
            'sqrt':('Pow(Integer(2),Rational(1,2))','Sqrt[2]'),
            're':('Integer(2)','Re[2+3*I]'), 'im':('Integer(3)','Im[2+3*I]'),
        }
        # Index fixtures by the independent native source head, never by the
        # table's mutable target. A wrong target must reach the comparison.
        case_keys=dict(zip(
            ('Plus','Times','Power','Equal','Unequal','Greater','Less','GreaterEqual','LessEqual',
             'List','Integrate','Inactive/Limit','Element','ForAll','Exists','And','Or','Not',
             'Piecewise','Rational','Min','Max','Abs','Sin','Cos','Tan','Sec','ArcCos','Log','Exp','Sqrt','Re','Im'),
            ('Add','Mul','Pow','Equality','Unequality','StrictGreaterThan','StrictLessThan','GreaterThan','LessThan',
             'Sequence','Integral','Limit','Contains','ForAll','Exists','And','Or','Not','Piecewise','Rational',
             'Min','Max','Abs','sin','cos','tan','sec','acos','log','exp','sqrt','re','im')))
        self.assertEqual({r.wl[0] for r in c.FUNCTION_TABLE},set(case_keys))
        scalar_identity=self.fixture('Integer(7)','7')[0][0]['zero_class']
        boolean_identity=self.fixture('true','True')[0][0]['samples'][0]['comparison']
        for row in c.FUNCTION_TABLE:
            case_key=case_keys[row.wl[0]]
            py,wl=cases[case_key]
            for alias in row.wl:
                original=row.wl[0]
                source=lambda name: ('Inactive['+name.split('/')[1]+']') if name.startswith('Inactive/') else name
                native=wl.replace(source(original),source(alias))
                with self.subTest(function=row.py,alias=alias):
                    before=self.fixture(py,native)[0][0]
                    if before['outcome']=='predicate':
                        self.assertTrue(all(r['comparison']==boolean_identity for r in before['samples']))
                    else:
                        self.assertEqual(before.get('zero_class'),scalar_identity)
                        if before['samples']: self.assertTrue(all(c.sample_formed(r['comparison']) for r in before['samples']))
                    # Repoint to a distinct head. Numeric constructors use an
                    # executable alternative; unsupported heads must change the
                    # computed coverage outcome, not merely the operand text.
                    target={'sqrt':'cbrt','Add':'Mul','Mul':'Add','sin':'cos','cos':'sin',
                            'tan':'sin','sec':'cos','acos':'sin','log':'exp','exp':'log',
                            'Abs':'re','Min':'Max','Max':'Min','re':'im','im':'re'}.get(case_key,'unbound_repoint_'+case_key)
                    altered=tuple(replace(r,py=target) if r==row else r for r in c.FUNCTION_TABLE)
                    # A swap keeps the table injective when the new target is
                    # another row's target, just as physical-name repoints do.
                    altered=tuple(replace(r,py=row.py) if r.py==target and r.wl!=row.wl else r for r in altered)
                    with patch.object(c,'FUNCTION_TABLE',altered): after=self.fixture(py,native)[0][0]
                    self.assertNotEqual(products(before),products(after))
                    with patch.object(c,'mapped_function',return_value=None):
                        disabled_before=self.fixture(py,native)[0]
                        with patch.object(c,'FUNCTION_TABLE',altered): disabled_after=self.fixture(py,native)[0]
                        with self.assertRaises(AssertionError): self.assertNotEqual(products(disabled_before),products(disabled_after))
        row=c.FUNCTION_TABLE[0]
        with self.assertRaises(c.InputError): c.function_bindings(c.FUNCTION_TABLE+(row,))

    def test_wrong_sqrt_target_fails_through_comparison(self):
        altered=tuple(replace(r,py='cbrt') if 'Sqrt' in r.wl else r for r in c.FUNCTION_TABLE)
        # The table is still injective and its source aliases are unchanged.
        # The native fixture must fail on the resulting mathematical residual.
        with patch.object(c,'FUNCTION_TABLE',altered), patch.object(self,'subTest',side_effect=lambda **kw:nullcontext()):
            c.function_bindings()
            with self.assertRaisesRegex(AssertionError,"nonzero.*zero"):
                self.test_every_function_binding_repoint()

    def test_ungrammatical_input(self):
        for py,wl in [('Integer(','7'),('Integer(7)','{7'),('evil.attribute()','7')]:
            with self.assertRaises(c.InputError): self.fixture(py,wl)

    def test_union_complement_coverage(self):
        result=self.fixture('Or(StrictGreaterThan(Symbol("GM"),Integer(0)),StrictLessThan(Symbol("GM"),Integer(0)))','GM>0 || GM<0')[0][0]
        self.assertTrue(any(s['piece']=='union complement' for s in result['strata']))
        self.assertTrue(any(s['coverage']=='sampled' for s in result['strata'] if s['piece']=='union complement'))

    def test_unresolved_strata_visible(self):
        result=self.fixture('Equality(Mystery(Symbol("GM")),Integer(0))','Mystery[GM]==0')[2]
        inventories=[r for r in result if r['kind']=='stratum_inventory']
        self.assertTrue(all(r['strata'] for r in inventories))


if __name__=='__main__':
    unittest.main(verbosity=2)
