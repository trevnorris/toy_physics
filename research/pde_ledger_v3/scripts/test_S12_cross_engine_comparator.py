#!/usr/bin/env python3
"""Synthetic native-stream controls. No measured stream is opened.

Each defect control changes a computed channel with the other raw operand
held fixed, then demonstrates that assertion ceases to distinguish the
mutation when its comparison channel is disabled. No measured residual target.
"""
import argparse
from dataclasses import replace
import io
import json
from pathlib import Path
import resource
import tempfile
import time
import unittest
from unittest.mock import patch
import S12_cross_engine_comparator as c


WORK=Path('/var/projects/toy_physics/_scratch/s12/cmp_work')
ARGS_P='t,x1,x2,x3,w'
ARGS_W='t,x1,x2,x3,w'
FACE_P='wminus(t,x1,x2,x3)'
FACE_W='wm[t,x1,x2,x3]'
OTHER_W='wp[t,x1,x2,x3]'
DOMAIN_TAG='fixture'
SYNTH=c.Correspondence(DOMAIN_TAG,'value','synthetic scalar','synthetic fixture','synthetic fixture')


def computed(obj,samples_only=False):
    """No operand echoes, operand digests, coordinates, names or sample values.

    Coverage reductions are tested using only computed evaluation residuals;
    symbolic subtraction cannot rescue reduced evaluation coverage.
    """
    if samples_only:
        residuals=[]
        def collect(x):
            if isinstance(x,dict):
                if 'samples' in x:
                    residuals.extend(sorted({json.dumps({k:v for k,v in row.items()
                                      if k in ('residual','unresolved_coverage')},sort_keys=True)
                                      for row in x['samples']}))
                else:
                    for v in x.values(): collect(v)
            elif isinstance(x,list):
                for v in x: collect(v)
        collect(obj)
        return sorted(set(residuals))
    if isinstance(obj,list):
        return [computed(x) for x in obj]
    if not isinstance(obj,dict):
        return obj
    keep={'as_computed_residual','on_supplied_domain_residual','text_difference','equal_text','first_different_character','not_compared','not_formed',
          'py_coverage','wl_coverage','unpaired_components','empty_container_comparison'}
    ans={k:computed(v) for k,v in obj.items() if k in keep}
    for k in ('components','distribution_parts'):
        if k in obj: ans[k]=computed(obj[k])
    return ans


def clear_evaluation_caches():
    c.test_expression.cache_clear()
    c.evaluate.cache_clear()


class Controls(unittest.TestCase):
    def fixture(self,py,wl,*,corr=SYNTH,names=c.NAME_TABLE,suffix_py='',suffix_wl='',sample_options=None):
        WORK.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='fixture_',dir=WORK) as folder:
            p,w=Path(folder)/'py.out',Path(folder)/'wl.out'
            p.write_text(f'{DOMAIN_TAG}: {py}\n'+suffix_py)
            w.write_text(f'{DOMAIN_TAG}: {wl}\n'+suffix_wl)
            out=io.StringIO()
            c.run(p,w,out,names=names,correspondences=(corr,),progress=False,
                  sample_options=sample_options,probe_metadata=False)
        records=[json.loads(line) for line in out.getvalue().splitlines()]
        comparisons=[r['comparison'] for r in records if r['kind']=='comparison']
        accounting=[r for r in records if r['kind']=='accounting']
        return comparisons,accounting,records

    def assert_moved(self,before,after,*,samples_only=False,accounting=False):
        old=[r for r in before[2] if r['kind']=='operands']
        new=[r for r in after[2] if r['kind']=='operands']
        self.assertEqual([r['py'] for r in old],[r['py'] for r in new])
        if accounting:
            def count_products(rows):
                return [{k:v for k,v in r.items() if k in ('accounting','parsed_leaves','actually_compared_leaves',
                                                          'uncompared_leaf_paths','reason')} for r in rows]
            self.assertNotEqual(count_products(before[1]),count_products(after[1]))
        else:
            self.assertNotEqual(computed(before[0],samples_only),computed(after[0],samples_only))

    def mutation(self,py,wl,changed,*,corr=SYNTH,channel='scalar',samples_only=False,
                 disabled=None,options=None):
        before=self.fixture(py,wl,corr=corr,sample_options=options)
        after=self.fixture(py,changed,corr=corr,sample_options=options)
        self.assert_moved(before,after,samples_only=samples_only)
        target=c.Comparator if channel=='scalar' else c
        if disabled is None: disabled={}
        with patch.object(target,channel,return_value=disabled):
            with self.assertRaises(AssertionError):
                self.assert_moved(self.fixture(py,wl,corr=corr,sample_options=options),
                                  self.fixture(py,changed,corr=corr,sample_options=options),samples_only=samples_only)
        return before,after

    def test_one_sided_operand_corruption(self):
        self.mutation('x1+7','x1+7','x1+11',channel='compare_values')

    def test_form_mutation(self):
        self.mutation('x1*x2','x1*x2','x1+x2^2',channel='compare_values')

    def test_coefficient_and_sign_mutation(self):
        self.mutation('3*x1','3*x1','-5*x1',channel='compare_values')

    def test_every_name_table_repoint(self):
        for row in c.NAME_TABLE:
            with self.subTest(name=row.py):
                names=c.NAME_TABLE
                other=next((r for r in names if r!=row and r.arity==row.arity),None)
                if other is None:
                    other=c.NameRow('synthetic_other_order','syntheticOtherOrder','different synthetic two-offset order',
                                    'synthetic fixture','synthetic fixture',row.arity)
                    names=names+(other,)
                altered=tuple(replace(r,wl=other.wl) if r==row else replace(r,wl=row.wl) if r==other else r
                              for r in names)
                args=','.join(str(i+1) for i in range(row.arity))
                py=f'{row.py}({args})' if row.arity else row.py
                wl=f'{row.wl}[{args}]' if row.arity else row.wl
                before=self.fixture(py,wl,names=names)
                after=self.fixture(py,wl,names=altered)
                self.assert_moved(before,after)
                def residuals(result):
                    return [r['as_computed_residual'] for r in result[2] if r['kind']=='residual_before_domain']
                self.assertNotEqual(residuals(before),residuals(after))
                with patch.object(c.Comparator,'scalar',return_value={}):
                    with self.assertRaises(AssertionError):
                        self.assert_moved(self.fixture(py,wl,names=names),self.fixture(py,wl,names=altered))

    def test_stripped_applied_function_argument(self):
        self.mutation(f'n({ARGS_P})',f'n[{ARGS_W}]','n')
        self.mutation(f'n({ARGS_P})',f'n[{ARGS_W}]','n[t,x1,x2,x3]')

    def test_live_function_replaced_by_constant(self):
        self.mutation(f'f({ARGS_P})',f'f[{ARGS_W}]','3',channel='compare_values')

    def test_delta_coefficient_moved_or_removed(self):
        py=f'(w+2)*DiracDelta(w-{FACE_P})'
        wl=f'(w+2)*DiracDelta[w-{FACE_W}]'
        disabled=({'regular/below':c.s.S.Zero,'regular/inside':c.s.S.Zero,'regular/above':c.s.S.Zero,
                   'delta/lower':c.s.S.Zero,'delta/upper':c.s.S.Zero},[])
        for changed in (f'(w+2)*DiracDelta[w-{OTHER_W}]','0'):
            with self.subTest(mutation=changed):
                self.mutation(py,wl,changed,channel='distribution_parts',disabled=disabled)

    def test_delta_normalization_in_bulk_coordinate(self):
        self.mutation(f'DiracDelta(2*(w-{FACE_P}))',f'DiracDelta[2*(w-{FACE_W})]',
                      f'DiracDelta[3*(w-{FACE_W})]',channel='compare_values')

    def test_only_one_region_regular_part(self):
        py=f'Heaviside(w-{FACE_P})*Heaviside(wplus(t,x1,x2,x3)-w)'
        wl=f'HeavisideTheta[w-{FACE_W}]-HeavisideTheta[w-{OTHER_W}]'
        self.mutation(py,wl,f'2*({wl})',channel='region_part',disabled=(c.s.S.Zero,[]))

    def test_difference_only_on_order_vanishing_stratum(self):
        py='0'
        wl='0'
        changed=f'Piecewise[{{{{t+2,f[t,x1,x2,x3,{FACE_W}]==0}}}},0]'
        self.mutation(py,wl,changed,samples_only=True,channel='sample_values',disabled=[])

    def test_difference_only_on_jump_vanishing_stratum(self):
        cond=f'nRegion1[t,x1,x2,x3,{FACE_W}]==nRegion2[t,x1,x2,x3,{FACE_W}]'
        changed=f'Piecewise[{{{{x1+2,{cond}}}}},0]'
        self.mutation('0','0',changed,samples_only=True,channel='sample_values',disabled=[])

    def test_changed_label(self):
        corr=replace(SYNTH,kind='text')
        self.mutation('"a supplied premise"','"a supplied premise"','"a different premise"',
                      corr=corr,channel='text_difference')

    def test_parsed_uncompared_leaf_count(self):
        corr=replace(SYNTH,kind='nested_ledger')
        before=self.fixture('[(x1,1)]','{{x1,1}}',corr=corr)
        after=self.fixture('[(x1,1)]','{{x1,1},{x2,2}}',corr=corr)
        self.assert_moved(before,after,accounting=True)
        with patch.object(c,'census',return_value=0),patch.object(c,'leaf_paths',return_value=set()):
            with self.assertRaises(AssertionError):
                self.assert_moved(self.fixture('[(x1,1)]','{{x1,1}}',corr=corr),
                                  self.fixture('[(x1,1)]','{{x1,1},{x2,2}}',corr=corr),accounting=True)

    def test_changed_derivative_structure(self):
        self.mutation(f'Derivative(n({ARGS_P}),t)',f'Derivative[1,0,0,0,0][n][{ARGS_W}]',
                      f'Derivative[0,1,0,0,0][n][{ARGS_W}]',channel='sample_values',disabled=[],samples_only=True)

    def test_derivative_at_composite_argument(self):
        py='Subs(Derivative(G(a,b),a),(a,b),(w-wminus(t,x1,x2,x3),wplus(t,x1,x2,x3)-w))'
        wl=f'Derivative[1,0][g][w-{FACE_W},{OTHER_W}-w]'
        changed=f'Derivative[0,1][g][w-{FACE_W},{OTHER_W}-w]'
        self.mutation(py,wl,changed,channel='sample_values',disabled=[],samples_only=True)

    def test_derivative_of_composite_argument_chain(self):
        py='Derivative(G(t**2,x1),t)'
        wl='2*t*Derivative[1,0][g][t^2,x1]'
        changed='Derivative[1,0][g][t^2,x1]'
        self.mutation(py,wl,changed,channel='sample_values',disabled=[],samples_only=True)

    def test_generic_profiles_separate_each_argument(self):
        for row in c.NAME_TABLE:
            if not row.arity: continue
            for form in c.FORMS:
                for index in range(3):
                    with self.subTest(name=row.py,form=form,point=index):
                        point=c.point_for(form,index,'inside')
                        actual=tuple(point[x] for x in c.COORDS[:row.arity])
                        if row.py=='G':
                            lo,hi=c.face_forms(form)
                            actual=(point[c.W]-lo.subs(point),hi.subs(point)-point[c.W])
                        env=dict(zip(c.FORMAL,actual))
                        formula=c.profile(row.py,row.arity,form,'generic')
                        derivatives=[c.s.diff(formula,x).subs(env) for x in c.FORMAL[:row.arity]]
                        for i,left in enumerate(derivatives):
                            for right in derivatives[i+1:]:
                                self.assertNotEqual(left,right)

    def test_vanishing_order_profiles_separate_arguments_at_face(self):
        for name,arity in (('f',5),('G',2)):
            for form in c.FORMS:
                for index in range(3):
                    for stratum in (x for x in c.STRATA if 'order_zero' in x):
                        loc='lower' if stratum.endswith('lower') else 'upper'
                        with self.subTest(name=name,form=form,stratum=stratum,point=index):
                            point=c.point_for(form,index,loc)
                            actual=tuple(point[x] for x in c.COORDS[:arity])
                            if name=='G':
                                lo,hi=c.face_forms(form)
                                actual=(point[c.W]-lo.subs(point),hi.subs(point)-point[c.W])
                            env=dict(zip(c.FORMAL,actual))
                            formula=c.profile(name,arity,form,stratum)
                            derivatives=[c.s.diff(formula,x).subs(env) for x in c.FORMAL[:arity]]
                            for i,left in enumerate(derivatives):
                                for right in derivatives[i+1:]:
                                    self.assertNotEqual(left,right)

    def test_nested_sibling_removed(self):
        schema=(('first',('a','b')),('second',('c','d')))
        self.mutation('((x1,x2),(x3,w))','{{x1,x2},{x3,w}}','{{x1,x2},{x3}}',
                      corr=replace(SYNTH,kind='tree',schema=schema),channel='shaped',disabled=False)

    def test_one_on_supplied_domain_factor_changed(self):
        py=f'x1*Heaviside(wplus(t,x1,x2,x3)-{FACE_P})'
        wl=f'x1*HeavisideTheta[{OTHER_W}-{FACE_W}]'
        before,after=self.mutation(py,wl,'x1',channel='compare_values')
        # This specifically observes the unconditioned subtraction channel.
        def raw(rows):
            return [r['as_computed_residual'] for r in rows[2] if r['kind']=='residual_before_domain']
        self.assertNotEqual(raw(before),raw(after))

    def test_component_correspondence_repoint(self):
        corr=replace(SYNTH,kind='tree',schema=('first coordinate','second coordinate'))
        swapped=replace(corr,paths=(('first coordinate',(0,),(1,)),('second coordinate',(1,),(0,))))
        before=self.fixture('(x1,x2)','{x1,x2}',corr=corr)
        after=self.fixture('(x1,x2)','{x1,x2}',corr=swapped)
        self.assert_moved(before,after)
        with patch.object(c.Comparator,'scalar',return_value={}):
            with self.assertRaises(AssertionError):
                self.assert_moved(self.fixture('(x1,x2)','{x1,x2}',corr=corr),
                                  self.fixture('(x1,x2)','{x1,x2}',corr=swapped))

    def test_joined_structure_mismatch_is_not_formed(self):
        corr=replace(SYNTH,kind='tree',schema=('first coordinate','second coordinate'))
        before,after=self.mutation('(x1,x2)','{x1,x2}','{{x1,x2}}',corr=corr,channel='shaped',disabled=False)
        self.assertTrue(c.has_key(after[0],'not_formed'))
        self.assertFalse(c.has_key(after[0],'as_computed_residual'))

    def test_single_point_coverage_reduction_misses_fixture(self):
        before=self.fixture('0','0')
        after=self.fixture('0','t-1/7')
        self.assert_moved(before,after,samples_only=True)
        options={'points':(0,)}
        with self.assertRaises(AssertionError):
            self.assert_moved(self.fixture('0','0',sample_options=options),
                              self.fixture('0','t-1/7',sample_options=options),samples_only=True)

    def test_single_form_coverage_reduction_misses_fixture(self):
        changed=f'Derivative[2,0,0,0,0][n][{ARGS_W}]'
        self.assert_moved(self.fixture('0','0'),self.fixture('0',changed),samples_only=True)
        options={'forms':('affine',)}
        with self.assertRaises(AssertionError):
            self.assert_moved(self.fixture('0','0',sample_options=options),
                              self.fixture('0',changed,sample_options=options),samples_only=True)

    def test_unseparated_arguments_coverage_reduction_misses_fixture(self):
        py='Subs(Derivative(G(a,b),a),(a,b),(1,2))'
        wl='Derivative[1,0][g][1,2]'
        changed='Derivative[0,1][g][1,2]'
        self.assert_moved(self.fixture(py,wl),self.fixture(py,changed),samples_only=True)
        original=c.profile
        def symmetric(name,arity,form,stratum):
            return sum(c.FORMAL[:arity]) if name=='G' else original(name,arity,form,stratum)
        clear_evaluation_caches()
        with patch.object(c,'profile',side_effect=symmetric):
            with self.assertRaises(AssertionError):
                self.assert_moved(self.fixture(py,wl),self.fixture(py,changed),samples_only=True)
        clear_evaluation_caches()

    def test_equivalent_delta_coefficient_away_from_support(self):
        py=f'(w+2)*DiracDelta(w-{FACE_P})'
        wl=f'(w+2)*DiracDelta[w-{FACE_W}]'
        changed=f'(w+2+(w-{FACE_W})^2*(1+x1))*DiracDelta[w-{FACE_W}]'
        before,after=self.fixture(py,wl),self.fixture(py,changed)
        def deltas(result):
            return [computed(part) for row in result[0] for part in row['distribution_parts']
                    if part['part'].startswith('delta/')]
        self.assertEqual(deltas(before),deltas(after))

    def test_duplicate_name_or_join_row_rejected(self):
        with self.assertRaises(c.InputError):
            self.fixture('x1','x1',names=c.NAME_TABLE+(c.NAME_TABLE[0],))
        with self.assertRaises(c.InputError):
            c.validate_tables(c.NAME_TABLE,(SYNTH,SYNTH))

    def test_repeated_tag_rejects_unique_join_and_is_accounted(self):
        before=self.fixture('x1','x1')
        after=self.fixture('x1','x1',suffix_wl='fixture: x2\n')
        self.assertNotEqual([r.get('accounting') for r in before[1]],
                            [r.get('accounting') for r in after[1]])
        self.assertTrue(all(r.get('reason') for r in after[1]))

    def test_unknown_singular_atom_is_coverage(self):
        self.mutation(f'DiracDelta(w-{FACE_P})',f'DiracDelta[w-{FACE_W}]',
                      'DiracDelta[w^2-1]',channel='scalar')

    def test_operational_grammar_error(self):
        with self.assertRaises(c.InputError):
            self.fixture('x1','x1 + [')


class ProgressResult(unittest.TextTestResult):
    def startTest(self,test):
        print('WORKING_ON: '+test.id(),flush=True)
        super().startTest(test)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    started=time.monotonic()
    with args.output.open('x') as out:
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(Controls)
        result=unittest.TextTestRunner(stream=out,verbosity=2,resultclass=ProgressResult).run(suite)
        out.write(json.dumps({'tests_run':result.testsRun,'errors':len(result.errors),'failures':len(result.failures),
                              'runtime_seconds':time.monotonic()-started,
                              'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})+'\n')
    return 0 if result.wasSuccessful() else 1


if __name__=='__main__':
    raise SystemExit(main())
