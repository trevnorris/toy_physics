#!/usr/bin/env python3
"""Synthetic controls only: never read either accepted transcript.

Mutation assertions use computed residuals, evaluations, coverage or census,
never operand echoes. Every mutation channel is also disabled to demonstrate
that the assertion loses its sensitivity. No assertion specifies a residual.
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
import S14a_1_cross_engine_comparator as c

WORK=Path('/var/projects/toy_physics/_scratch/s14a/cmp_work')
PARGS='x,y,z,w,t'
WARGS='x1,x2,x3,w,t'
PLO='General_face_lower(x,y,z,t)'
PHI='General_face_upper(x,y,z,t)'
WLO='wLower[x1,x2,x3,t]'
WHI='wUpper[x1,x2,x3,t]'
SYN=c.Binding('synthetic',c.End('fixture py'),c.End('fixture wl'),'synthetic object','synthetic','synthetic')


def products(value,samples=False):
    result=[]
    def walk(x):
        if isinstance(x,dict):
            if samples and 'samples' in x:
                for sample in x['samples']:
                    result.append({key:sample[key] for key in ('residual','unresolved_coverage','value_coverage') if key in sample})
            elif not samples:
                for key in ('as_computed_residual','on_supplied_domain_residual','text_difference','not_compared','not_formed',
                            'py_coverage','wl_coverage'):
                    if key in x:
                        result.append({key:x[key]})
            for key in ('comparison','components','distribution_parts','domain_context_diagnostic'):
                if key in x:
                    walk(x[key])
        elif isinstance(x,(tuple,list)):
            for y in x:
                walk(y)
    walk(value)
    return result


def clear_caches():
    c.profile.cache_clear()
    c.face_forms.cache_clear()
    c.point_for.cache_clear()
    c.test_expression.cache_clear()
    c.evaluate.cache_clear()


def serialize_shape(shape,engine,offset=0):
    coord='x' if engine=='py' else 'x1'
    if shape==c.V:
        return f'{coord}+{offset+2}'
    if shape=='equation':
        return f'Eq({coord}+{offset+2},t+{offset+13})' if engine=='py' else f'Inactive[Equal][{coord}+{offset+2},t+{offset+13}]'
    if shape=='conditioned':
        return (f'{{domain: Ne({coord}+{offset+2},0), object: {coord}+{offset+7}}}' if engine=='py' else
                f'ConditionalExpression[{coord}+{offset+7},Unequal[{coord}+{offset+2},0]]')
    values=[serialize_shape(sub,engine,offset+17*i) for i,sub in enumerate(shape)]
    return ('('+','.join(values)+(',' if len(values)==1 else '')+')') if engine=='py' else '{'+','.join(values)+'}'


def wrap_end(value,path,engine,root_shape=None):
    if not path:
        return value
    part=path[0]
    inner=wrap_end(value,path[1:],engine,root_shape[part] if isinstance(part,int) and root_shape is not None else None)
    if isinstance(part,int):
        siblings=[serialize_shape(sub,engine,111+i*5) for i,sub in enumerate(root_shape or (c.V,)*(part+1))]
        siblings[part]=inner
        return ('('+','.join(siblings)+(',' if len(siblings)==1 else '')+')') if engine=='py' else '{'+','.join(siblings)+'}'
    if part=='value':
        return f'{{domain: Ne(t,0), object: {inner}}}' if engine=='py' else f'ConditionalExpression[{inner},Unequal[t,0]]'
    if part=='condition':
        return f'{{domain: {inner}, object: x+7}}' if engine=='py' else f'ConditionalExpression[x1+7,{inner}]'
    raise AssertionError('synthetic wrapping missing selector '+str(part))


class Controls(unittest.TestCase):
    def fixture(self,py,wl,*,row=SYN,rows=None,names=c.NAME_TABLE,extra_py='',extra_wl='',options=None):
        WORK.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='synthetic_',dir=WORK) as directory:
            p,w=Path(directory)/'py.out',Path(directory)/'wl.out'
            p.write_text(row.py.tag+': '+py+'\n'+extra_py)
            w.write_text(row.wl.tag+': '+wl+'\n'+extra_wl)
            out=io.StringIO()
            c.run(p,w,out,names=names,bindings=(row,) if rows is None else rows,
                  progress=False,options=options,probe_metadata=False)
        records=[json.loads(line) for line in out.getvalue().splitlines()]
        return records

    def comparison(self,records):
        return [r for r in records if r['kind']=='comparison']

    def moved(self,before,after,*,samples=False,unchanged='py',accounting=False):
        def operands(records):
            return [r[unchanged] for r in records if r['kind']=='operands']
        self.assertEqual(operands(before),operands(after))
        if accounting:
            def counts(records):
                return [{key:r[key] for key in ('parsed_leaves','actually_compared_leaves') if key in r}
                        for r in records if r['kind'] in ('binding_accounting','tag_accounting')]
            self.assertNotEqual(counts(before),counts(after))
        else:
            self.assertNotEqual(products(self.comparison(before),samples),products(self.comparison(after),samples))

    def mutation(self,py,wl,changed,*,row=SYN,channel='compare_values',disabled=None,samples=False,options=None):
        before=self.fixture(py,wl,row=row,options=options)
        after=self.fixture(py,changed,row=row,options=options)
        self.moved(before,after,samples=samples)
        target=c.Comparator if channel in ('component','compare') else c
        with patch.object(target,channel,return_value={} if disabled is None else disabled):
            with self.assertRaises(AssertionError):
                self.moved(self.fixture(py,wl,row=row,options=options),self.fixture(py,changed,row=row,options=options),samples=samples)
        return before,after

    def test_operand_form_coefficient_sign(self):
        for py,wl,mutated in (('x+7','x1+7','x1+11'),('x*y','x1*x2','x1+x2^2'),('3*x','3*x1','-5*x1')):
            with self.subTest(mutated=mutated):
                self.mutation(py,wl,mutated)

    def test_repoint_every_binding(self):
        for row in c.BINDING_TABLE:
            with self.subTest(binding=row.key):
                print('WORKING_ON: binding '+row.key+' '+row.object,flush=True)
                text=all(part.kind=='text' for part in row.components)
                py='native scope' if text else serialize_shape(row.py_shape,'py')
                wl=json.dumps('native scope') if text else serialize_shape(row.wl_shape,'wl')
                other=json.dumps('changed scope') if text else serialize_shape(row.wl_shape,'wl',97)
                py=wrap_end(py,row.py.path,'py',c.ROOT_SCHEMAS.get(('py',row.py.tag)))
                wl=wrap_end(wl,row.wl.path,'wl',c.ROOT_SCHEMAS.get(('wl',row.wl.tag)))
                other=wrap_end(other,row.wl.path,'wl',c.ROOT_SCHEMAS.get(('wl',row.wl.tag)))
                alternate='alternate/comparison' if text else 'alternate'
                extra=alternate+': '+other+'\n'
                repointed=replace(row,wl=c.End(alternate,row.wl.path))
                def measure(r):
                    return self.fixture(py,wl,row=row,rows=(r,),extra_wl=extra)
                before,after=measure(row),measure(repointed)
                self.moved(before,after)
                if not text:
                    self.moved(before,after,samples=True)
                with patch.object(c.Comparator,'component',return_value={}):
                    with self.assertRaises(AssertionError):
                        self.moved(measure(row),measure(repointed))
                if not text:
                    with patch.object(c,'sample_values',return_value=[]):
                        with self.assertRaises(AssertionError):
                            self.moved(measure(row),measure(repointed),samples=True)

    def test_repoint_every_name_row(self):
        for row in c.NAME_TABLE:
            with self.subTest(name=row.py):
                print('WORKING_ON: name '+row.py,flush=True)
                other=next(r for r in c.NAME_TABLE if r!=row and r.arity==row.arity)
                altered=tuple(replace(r,wl=other.wl) if r==row else replace(r,wl=row.wl) if r==other else r for r in c.NAME_TABLE)
                args=','.join(str(i+1) for i in range(row.arity))
                py=row.py+'('+args+')' if row.arity else row.py
                wl=row.wl+'['+args+']' if row.arity else row.wl
                def measure(names):
                    return self.fixture(py,wl,names=names)
                before,after=measure(c.NAME_TABLE),measure(altered)
                self.moved(before,after)
                self.moved(before,after,samples=True)
                with patch.object(c,'compare_values',return_value={}):
                    with self.assertRaises(AssertionError):
                        self.moved(measure(c.NAME_TABLE),measure(altered))
                with patch.object(c,'sample_values',return_value=[]):
                    with self.assertRaises(AssertionError):
                        self.moved(measure(c.NAME_TABLE),measure(altered),samples=True)

    def test_each_production_component_correspondence(self):
        for row in c.BINDING_TABLE:
            if len(row.components)<2:
                continue
            for i,component in enumerate(row.components):
                with self.subTest(binding=row.key,component=component.object):
                    j=(i+1)%len(row.components)
                    swapped=tuple(replace(part,wl_path=row.components[j].wl_path) if k==i else
                                  replace(part,wl_path=component.wl_path) if k==j else part
                                  for k,part in enumerate(row.components))
                    changed=replace(row,components=swapped)
                    py=wrap_end(serialize_shape(row.py_shape,'py'),row.py.path,'py',c.ROOT_SCHEMAS.get(('py',row.py.tag)))
                    wl=wrap_end(serialize_shape(row.wl_shape,'wl'),row.wl.path,'wl',c.ROOT_SCHEMAS.get(('wl',row.wl.tag)))
                    def measure(binding):
                        return self.fixture(py,wl,row=row,rows=(binding,))
                    before,after=measure(row),measure(changed)
                    self.moved(before,after)
                    self.moved(before,after,samples=True)
                    with patch.object(c.Comparator,'component',return_value={}):
                        with self.assertRaises(AssertionError):
                            self.moved(measure(row),measure(changed))

    def test_full_source_container_schema_before_selection(self):
        row=next(r for r in c.BINDING_TABLE if r.wl.tag=='general projection/material/regional currents' and r.wl.path==(0,1))
        ps=c.ROOT_SCHEMAS[('py',row.py.tag)]
        ws=c.ROOT_SCHEMAS[('wl',row.wl.tag)]
        py=wrap_end('x+2',row.py.path,'py',ps)
        wl=wrap_end('x1+2',row.wl.path,'wl',ws)
        changed=wrap_end('x1+2',row.wl.path,'wl',ws[:-1])
        def check():
            before=self.fixture(py,wl,row=row)
            after=self.fixture(py,changed,row=row)
            self.moved(before,after)
            self.assertTrue(c.has_key(self.comparison(after),'not_formed'))
            self.assertFalse(c.has_key(self.comparison(after),'as_computed_residual'))
        check()
        original=c.shaped
        def no_root(ref,shape):
            return True if shape==ws else original(ref,shape)
        with patch.object(c,'shaped',side_effect=no_root):
            with self.assertRaises(AssertionError):
                check()

    def test_split_against_grouped_every_counterpart(self):
        # Same production get_ref/component path as currents, face lists and
        # dimension fluxes. Each grouped counterpart is independently changed.
        rows=tuple(replace(SYN,key='split'+str(i),py=c.End('part '+str(i)),wl=c.End('grouped',(i,))) for i in range(3))
        py='x+2'
        wl='{x1+2,x2+3,t+5}'
        extra='part 1: y+3\npart 2: t+5\n'
        def measure(group):
            return self.fixture(py,group,row=rows[0],rows=rows,extra_py=extra)
        for changed in ('{x1+7,x2+3,t+5}','{x1+2,x2+9,t+5}','{x1+2,x2+3,t+11}'):
            self.moved(measure(wl),measure(changed))
            self.moved(measure(wl),measure(changed),samples=True)
            with patch.object(c.Comparator,'component',return_value={}):
                with self.assertRaises(AssertionError):
                    self.moved(measure(wl),measure(changed))

    def test_distinct_functions_in_each_form_and_stratum(self):
        fields=[r for r in c.NAME_TABLE if r.arity]
        row=replace(SYN,py_shape=(c.V,)*len(fields),wl_shape=(c.V,)*len(fields),
                    components=c.vector_components(tuple(r.py for r in fields)))
        py='('+','.join('0' for r in fields)+')'
        wl='{'+','.join(r.wl+'['+','.join(WARGS.split(',')[:r.arity])+']' for r in fields)+'}'
        def measure():
            return self.fixture(py,wl,row=row)
        def distinguish(records):
            parts=self.comparison(records)[0]['comparison']['components']
            for form in c.FORMS:
                for stratum in c.STRATA:
                    signatures={}
                    for field,part in zip(fields,parts):
                        signature=tuple(s['residual'] for s in part['samples'] if s['form']==form and s['stratum']==stratum)
                        self.assertTrue(signature)
                        group=signatures.setdefault(field.arity,set())
                        self.assertNotIn(signature,group,(field.py,form,stratum))
                        group.add(signature)
        original=c.profile
        source=next(r for r in fields if r.arity==5)
        target=next(r for r in fields if r.arity==5 and r!=source)
        def collision(name,arity,form,stratum):
            return original(source.py if name==target.py else name,arity,form,stratum)
        before=measure()
        distinguish(before)
        clear_caches()
        try:
            with patch.object(c,'profile',side_effect=collision):
                after=measure()
                self.moved(before,after,samples=True)
                with self.assertRaises(AssertionError):
                    distinguish(after)
            with patch.object(c,'sample_values',return_value=[]):
                with self.assertRaises((AssertionError,KeyError)):
                    distinguish(measure())
        finally:
            clear_caches()

    def test_stripped_argument_and_constant_replacement(self):
        for mutated in ('n1','n1[x1,x2,x3,w]','3'):
            self.mutation('General_below_n('+PARGS+')','n1['+WARGS+']',mutated,channel='component')

    def test_delta_face_removal_normalization_and_regular_region(self):
        py=f'(w+2)*DiracDelta(w-{PLO})'
        wl=f'(w+2)*DiracDelta[w-{WLO}]'
        for mutated in (f'(w+2)*DiracDelta[w-{WHI}]','0',f'(w+2)*DiracDelta[2*(w-{WLO})]'):
            self.mutation(py,wl,mutated)
        self.mutation(f'Heaviside(w-{PLO})-Heaviside(w-{PHI})',
                      f'HeavisideTheta[w-{WLO}]-HeavisideTheta[w-{WHI}]',
                      f'2*(HeavisideTheta[w-{WLO}]-HeavisideTheta[w-{WHI}])')

    def test_stratum_only_differences(self):
        conditions=(f'chi[x1,x2,x3,{WLO},t]==0',
                    f'n1[x1,x2,x3,{WLO},t]==n2[x1,x2,x3,{WLO},t]')
        for condition in conditions:
            self.mutation('0','0',f'Piecewise[{{{{x1+2,{condition}}}}},0]',samples=True,channel='sample_values',disabled=[])

    def test_integral_range_variable_integrand(self):
        py='Integral(x+w,(w,0,1))'
        wl='Inactive[Integrate][x1+w,{w,0,1}]'
        for changed in ('Inactive[Integrate][x1+w,{w,0,2}]','Inactive[Integrate][x1+w,{q,0,1}]',
                        'Inactive[Integrate][x1+w^2,{w,0,1}]'):
            self.mutation(py,wl,changed,samples=True,channel='sample_values',disabled=[])
            self.mutation(py,wl,changed,channel='compare_values')

    def test_nested_integral_variable_identity(self):
        py='Integral(radius**2*Integral(radius*w,(w,0,1)),(radius,1,2))'
        wl='Inactive[Integrate][r^2*Inactive[Integrate][r*w,{w,0,1}],{r,1,2}]'
        changed='Inactive[Integrate][r^2*Inactive[Integrate][w^2,{w,0,1}],{r,1,2}]'
        before,after=self.mutation(py,wl,changed,samples=True,channel='compare_values')
        self.moved(before,after)
        # The other engine's inner-variable substitution must be visible too.
        py_changed='Integral(radius**2*Integral(w**2,(w,0,1)),(radius,1,2))'
        self.moved(before,self.fixture(py_changed,wl),samples=True,unchanged='wl')
        with patch.object(c,'compare_values',return_value={}):
            with self.assertRaises(AssertionError):
                self.moved(self.fixture(py,wl),self.fixture(py_changed,wl),samples=True,unchanged='wl')

        # Reintroduce the actual scope defect, not merely a disabled parser.
        def capture(value):
            if not value.args:
                return value
            rebuilt=value.func(*(capture(a) for a in value.args))
            if isinstance(rebuilt,c.s.Integral):
                return rebuilt.xreplace({v:c.s.Symbol('_int'+str(i),real=True)
                                         for i,v in enumerate(rebuilt.variables)})
            return rebuilt
        with patch.object(c,'alpha_integrals',side_effect=capture):
            with self.assertRaises(AssertionError):
                self.moved(self.fixture(py,wl),self.fixture(py,changed),samples=True)
        self.equivalent(py,wl,'Inactive[Integrate][q^2*Inactive[Integrate][q*p,{p,0,1}],{q,1,2}]')

    def test_every_live_first_partial_at_another_argument_point(self):
        # H2 requires one discriminating form/stratum. The quadratic probes
        # cover every derivative slot crossed with every evaluation argument;
        # production still samples all forms, strata and rational points.
        options={'forms':('quadratic',),'strata':('generic',)}
        for field in (r for r in c.NAME_TABLE if r.arity):
            with self.subTest(function=field.py):
                print('WORKING_ON: derivative points '+field.py,flush=True)
                spherical=field.py.startswith('SphericalSpecialization_')
                chart='SphericalSpecialization' if spherical else 'General'
                pa=['radius','w','t'] if spherical else ['x','y','z','w','t']
                wa=['r','w','t'] if spherical else ['x1','x2','x3','w','t']
                if '_face_' in field.py:
                    pa.remove('w')
                    wa.remove('w')
                py_parts=[]
                wl_parts=[]
                changed_parts=[]
                labels=[]
                for derivative_slot in range(field.arity):
                    orders=','.join(str(int(k==derivative_slot)) for k in range(field.arity))
                    for point_slot in range(field.arity):
                        py_parts.append('Derivative('+field.py+'('+','.join(pa)+'),'+pa[derivative_slot]+')')
                        prefix='Derivative['+orders+']['+field.wl+']'
                        wl_parts.append(prefix+'['+','.join(wa)+']')
                        shifted=list(wa)
                        shifted[point_slot]='('+shifted[point_slot]+'+1/29)'
                        changed_parts.append(prefix+'['+','.join(shifted)+']')
                        labels.append('slot '+str(derivative_slot)+' at argument '+str(point_slot))
                row=replace(SYN,py_shape=(c.V,)*len(labels),wl_shape=(c.V,)*len(labels),
                            components=c.vector_components(labels),chart=chart)
                py='('+','.join(py_parts)+')'
                wl='{'+','.join(wl_parts)+'}'
                changed='{'+','.join(changed_parts)+'}'
                def measure(value):
                    return self.fixture(py,value,row=row,options=options)
                def check_component(before,after,index):
                    self.assertEqual([r['py'] for r in before if r['kind']=='operands'],
                                     [r['py'] for r in after if r['kind']=='operands'])
                    a=self.comparison(before)[0]['comparison']['components'][index]
                    b=self.comparison(after)[0]['comparison']['components'][index]
                    self.assertTrue(a.get('samples'))
                    self.assertTrue(b.get('samples'))
                    self.assertTrue(all('residual' in sample for part in (a,b) for sample in part['samples']))
                    self.assertNotEqual(products(a,True),products(b,True))
                before,after=measure(wl),measure(changed)
                for i,label in enumerate(labels):
                    with self.subTest(pair=label):
                        check_component(before,after,i)
                with patch.object(c,'sample_values',return_value=[]):
                    disabled_before,disabled_after=measure(wl),measure(changed)
                    for i in range(len(labels)):
                        with self.assertRaises(AssertionError):
                            check_component(disabled_before,disabled_after,i)

    def test_face_derivative_point_inside_shell_integral(self):
        pface='SphericalSpecialization_face_lower(radius,t)'
        py=f'Integral(radius**2*SphericalSpecialization_inside_n(radius,{pface},t)*Derivative({pface},t),(radius,inner_radius,outer_radius))'
        wl='Inactive[Integrate][r^2*sn2[r,swLower[r,t],t]*Derivative[0,1][swLower][r,t],{r,ri,ro}]'
        changed='Inactive[Integrate][r^2*sn2[r,swLower[r,t],t]*Derivative[0,1][swLower][ri,t],{r,ri,ro}]'
        row=replace(SYN,chart='SphericalSpecialization')
        self.mutation(py,wl,changed,row=row,samples=True,channel='sample_values',disabled=[])
        clear_caches()
        try:
            with patch.object(c,'argument_coupling',return_value=c.s.S.Zero):
                with self.assertRaises(AssertionError):
                    self.moved(self.fixture(py,wl,row=row),self.fixture(py,changed,row=row),samples=True)
        finally:
            clear_caches()

    def test_unjoined_spaced_comparison_note(self):
        tag='transverse integration / comparison'
        note='native comparison note'
        row=replace(SYN,py=c.End('contract'),wl=c.End(tag),components=(c.Component('note',kind='text'),))
        with patch.object(c,'read_label',wraps=c.read_label) as reader:
            self.mutation(note,json.dumps(note),json.dumps(note+' revised'),row=row,channel='component')
            reader.assert_any_call(json.dumps(note),'wl',tag)
        def check():
            records=self.fixture('x','x1',extra_wl=tag+': '+json.dumps(note)+'\n')
            entry=next(r for r in records if r['kind']=='tag_accounting' and r['engine']=='wl' and r['tag']==tag)
            self.assertEqual(entry.get('classification'),'engine-local comparison note')
            self.assertTrue(entry.get('read_as_label'))
            self.assertIn('current-run accounting',entry.get('reason_basis',''))
        check()
        with patch.object(c,'unjoined_accounting',return_value={}):
            with self.assertRaises(AssertionError):
                check()

    def test_freed_component_bindings_have_current_accounting(self):
        # These fixtures exercise every component-ended production binding.
        # Samples are immaterial to this accounting control; value comparison
        # and its independent parse census still use the production path.
        for row in c.BINDING_TABLE:
            if not (row.py.path or row.wl.path):
                continue
            with self.subTest(binding=row.key):
                py=wrap_end(serialize_shape(row.py_shape,'py'),row.py.path,'py',c.ROOT_SCHEMAS.get(('py',row.py.tag)))
                wl=wrap_end(serialize_shape(row.wl_shape,'wl'),row.wl.path,'wl',c.ROOT_SCHEMAS.get(('wl',row.wl.tag)))
                def measure(rows):
                    return self.fixture(py,wl,row=row,rows=rows,options={'forms':()})
                before=measure((row,))
                def check(after):
                    tags={(r['engine'],r['tag']):r for r in after if r['kind']=='tag_accounting'}
                    for engine,other in (('py','wl'),('wl','py')):
                        entry=tags[engine,getattr(row,engine).tag]
                        self.assertEqual(entry['accounting'],'unjoined')
                        self.assertIn('curated_binding_candidates',entry)
                        candidate=next(x for x in entry['curated_binding_candidates'] if x['curated_binding']==row.key)
                        self.assertFalse(candidate['active_binding_rows'])
                        self.assertFalse(candidate['counterpart_tag_accounted_bindings'])
                        counterpart=tags[other,getattr(row,other).tag]
                        self.assertEqual(candidate['counterpart_occurrences'],
                                         sum(r['kind']=='tag_accounting' and r.get('engine')==other and r.get('tag')==getattr(row,other).tag for r in after))
                        self.assertEqual(candidate['counterpart_tag_compared_leaves'],counterpart['actually_compared_leaves'])
                        self.assertIn('computed',entry['reason_basis'])
                        prior=next(r for r in before if r['kind']=='tag_accounting' and r['engine']==engine and r['tag']==entry['tag'])
                        self.assertNotEqual(prior['accounting'],entry['accounting'])
                check(measure(()))
                with patch.object(c,'unjoined_accounting',return_value={'reason':'counterpart bound or absent'}):
                    with self.assertRaises(AssertionError):
                        check(measure(()))

    def test_condition_only_mutation_and_failure_stratum(self):
        row=next(r for r in c.BINDING_TABLE if r.py.tag=='SphericalSpecialization/outer projected velocity')
        py='{domain: Ne(Integral(SphericalSpecialization_inside_n(radius,w,t),(w,SphericalSpecialization_face_lower(radius,t),SphericalSpecialization_face_upper(radius,t))),0), object: radius+2}'
        wl='ConditionalExpression[r+2,Unequal[Inactive[Integrate][sn2[r,w,t],{w,swLower[r,t],swUpper[r,t]}],0]]'
        changed='ConditionalExpression[r+2,r>0]'
        before,after=self.mutation(py,wl,changed,row=row)
        self.moved(before,after,samples=True)
        def excluded(records):
            return [sample.get('conditions') for result in self.comparison(records) for part in result['comparison']['components']
                    for sample in part.get('samples',[]) if sample['stratum']=='density_zero']
        self.assertNotEqual(excluded(before),excluded(after))

    def test_dimension_value_domain_context(self):
        for row in c.BINDING_TABLE:
            if not row.condition_context:
                continue
            with self.subTest(binding=row.key):
                py=serialize_shape(row.py_shape,'py')
                vector=serialize_shape(row.wl_shape,'wl')
                wl='ConditionalExpression['+vector+',Unequal[r,0]]'
                changed='ConditionalExpression['+vector+',False]'
                self.mutation(py,wl,changed,row=row,samples=True,channel='compare_values')
                before=self.fixture(py,wl,row=row)
                after=self.fixture(py,changed,row=row)
                self.moved(before,after)

    def test_multiline_payload_changed_line(self):
        row=replace(SYN,py_shape=(c.V,c.V),wl_shape=(c.V,c.V),components=c.vector_components(('first','second')))
        before=self.fixture('(\n x,\n y\n)','{x1,x2}',row=row)
        after=self.fixture('(\n x,\n y+7\n)','{x1,x2}',row=row)
        self.moved(before,after,unchanged='wl')
        with patch.object(c.Comparator,'component',return_value={}):
            with self.assertRaises(AssertionError):
                self.moved(self.fixture('(\n x,\n y\n)','{x1,x2}',row=row),
                           self.fixture('(\n x,\n y+7\n)','{x1,x2}',row=row),unchanged='wl')

    def label_measure(self,py,wl,list_label=False):
        row=(replace(SYN,py=c.End('engine'),wl=c.End('engine'),py_shape=(c.V,c.V),wl_shape=(c.V,c.V),
                     components=c.vector_components(('runtime label','engine label'),'text')) if list_label else
             replace(SYN,py=c.End('contract'),wl=c.End('scope/comparison'),components=(c.Component('scope',kind='text'),)))
        p=('('+', '.join(py)+')') if list_label else py
        w=('{'+', '.join(json.dumps(x) for x in wl)+'}') if list_label else json.dumps(wl)
        return self.fixture(p,w,row=row)

    def label_control(self,list_label):
        value=('runtime (branch, date) [compiler]','library version') if list_label else 'supplied premise (live fields)'
        for side in ('py','wl'):
            for index in range(len(value) if list_label else 1):
                changed=tuple(x+' revised' if i==index else x for i,x in enumerate(value)) if list_label else value+' revised'
                before=self.label_measure(value,value,list_label)
                after=self.label_measure(changed,value,list_label) if side=='py' else self.label_measure(value,changed,list_label)
                self.moved(before,after,unchanged='wl' if side=='py' else 'py')

    def test_labels_in_both_native_serializations(self):
        for is_list in (False,True):
            self.label_control(is_list)
            for channel in ('text_difference','read_label'):
                original=c.read_label
                replacement=(lambda text,engine,tag:c.seq(c.atom('Text',''),c.atom('Text','')) if tag=='engine' else c.atom('Text',''))
                with patch.object(c,channel,**({'side_effect':replacement} if channel=='read_label' else {'return_value':{}})):
                    with self.assertRaises(AssertionError):
                        self.label_control(is_list)
            with patch.object(c.Comparator,'component',return_value={}):
                with self.assertRaises(AssertionError):
                    self.label_control(is_list)

    def test_parsed_uncompared_leaf_census(self):
        row=replace(SYN,py_shape=(c.V,c.V),wl_shape=(c.V,c.V),components=(c.Component('only first',(0,),(0,)),))
        before=self.fixture('(x,y)','{x1,x2}',row=row)
        after=self.fixture('(x,y)','{x1,x2+t}',row=row)
        self.moved(before,after,accounting=True)
        self.assertEqual(products(self.comparison(before)),products(self.comparison(after)))
        with patch.object(c,'census',return_value=0),patch.object(c,'leaf_paths',return_value=set()):
            with self.assertRaises(AssertionError):
                self.moved(self.fixture('(x,y)','{x1,x2}',row=row),self.fixture('(x,y)','{x1,x2+t}',row=row),accounting=True)

    def test_derivative_slot_order_point_and_composite(self):
        py=f'Derivative(General_below_n({PARGS}),x)'
        wl=f'Derivative[1,0,0,0,0][n1][{WARGS}]'
        for changed in (f'Derivative[0,1,0,0,0][n1][{WARGS}]',f'Derivative[2,0,0,0,0][n1][{WARGS}]',
                        'Derivative[1,0,0,0,0][n1][x1+1,x2,x3,w,t]'):
            self.mutation(py,wl,changed,samples=True,channel='sample_values',disabled=[])
        py='Subs(Derivative(General_below_n(a,b,c,d,e),a),(a,b,c,d,e),(x**2,y,z,w,t))'
        wl='Derivative[1,0,0,0,0][n1][x1^2,x2,x3,w,t]'
        self.mutation(py,wl,'Derivative[0,1,0,0,0][n1][x1^2,x2,x3,w,t]',samples=True,channel='sample_values',disabled=[])
        py=f'Derivative(General_below_n(x**2,y,z,w,t),x)'
        wl='2*x1*Derivative[1,0,0,0,0][n1][x1^2,x2,x3,w,t]'
        self.mutation(py,wl,'Derivative[1,0,0,0,0][n1][x1^2,x2,x3,w,t]',samples=True,channel='sample_values',disabled=[])

    def test_nested_sibling_and_structure_mismatch(self):
        row=replace(SYN,py_shape=((c.V,c.V),c.V),wl_shape=((c.V,c.V),c.V),
                    components=(c.Component('nested first',(0,0),(0,0)),c.Component('nested second',(0,1),(0,1)),c.Component('last',(1,),(1,))))
        before,after=self.mutation('((x,y),t)','{{x1,x2},t}','{{x1},t}',row=row,channel='shaped',disabled=False)
        self.assertTrue(c.has_key(self.comparison(after),'not_formed'))
        self.assertFalse(c.has_key(self.comparison(after),'as_computed_residual'))

    def test_domain_factor_and_supplied_equation_are_not_applied(self):
        self.mutation(f'Heaviside({PHI}-{PLO})*x',f'HeavisideTheta[{WHI}-{WLO}]*x1','x1')
        self.mutation('General_below_GammaB('+PARGS+')','gamma1['+WARGS+']',f'return1[{WARGS}]-drain1[{WARGS}]')
        # U1, committed balances, and the face factor are never input rewrite rules.
        self.mutation('General_inside_chi('+PARGS+')','chi['+WARGS+']','0')

    def test_coverage_one_point_only(self):
        changed='t-1/7'
        self.moved(self.fixture('0','0'),self.fixture('0',changed),samples=True)
        with self.assertRaises(AssertionError):
            self.moved(self.fixture('0','0',options={'points':(0,)}),self.fixture('0',changed,options={'points':(0,)}),samples=True)

    def test_coverage_one_form_only(self):
        changed=f'Derivative[2,0,0,0,0][n1][{WARGS}]'
        options={'strata':('generic',)}
        self.moved(self.fixture('0','0',options=options),self.fixture('0',changed,options=options),samples=True)
        options.update(forms=('affine',))
        with self.assertRaises(AssertionError):
            self.moved(self.fixture('0','0',options=options),self.fixture('0',changed,options=options),samples=True)

    def test_coverage_nonseparating_forms_only(self):
        py='Subs(Derivative(General_below_n(a,b,c,d,e),a),(a,b,c,d,e),(1,2,3,4,5))'
        wl='Derivative[1,0,0,0,0][n1][1,2,3,4,5]'
        changed='Derivative[0,1,0,0,0][n1][1,2,3,4,5]'
        options={'strata':('generic',)}
        self.moved(self.fixture(py,wl,options=options),self.fixture(py,changed,options=options),samples=True)
        def symmetric(form,args,code):
            total=sum(args)
            degree=c.FORMS.index(form)+1
            return (code+1)/7+(code+1)*total**degree/101
        clear_caches()
        try:
            with patch.object(c,'basic_form',side_effect=symmetric):
                with self.assertRaises(AssertionError):
                    self.moved(self.fixture(py,wl,options=options),self.fixture(py,changed,options=options),samples=True)
        finally:
            clear_caches()

    def equivalent(self,py,wl,changed,**kwargs):
        before=self.fixture(py,wl,**kwargs)
        after=self.fixture(py,changed,**kwargs)
        def check(left,right):
            exact=lambda r:[x for x in products(self.comparison(r)) if 'as_computed_residual' in x]
            self.assertTrue(exact(left))
            self.assertTrue(products(self.comparison(left),True))
            self.assertEqual(products(self.comparison(left)),products(self.comparison(right)))
            self.assertEqual(products(self.comparison(left),True),products(self.comparison(right),True))
        check(before,after)
        with patch.object(c,'compare_values',return_value={}):
            with self.assertRaises(AssertionError):
                check(self.fixture(py,wl,**kwargs),self.fixture(py,changed,**kwargs))
        return before,after

    def test_equivalence_delta_changed_away_from_face(self):
        py=f'(w+2)*DiracDelta(w-{PLO})'
        wl=f'(w+2)*DiracDelta[w-{WLO}]'
        self.equivalent(py,wl,f'(w+2+(w-{WLO})^2*(x1+1))*DiracDelta[w-{WLO}]')
        with patch.object(c,'distribution_parts',side_effect=lambda expr,chart:({'raw':expr},[])):
            with self.assertRaises(AssertionError):
                self.equivalent(py,wl,f'(w+2+(w-{WLO})^2*(x1+1))*DiracDelta[w-{WLO}]')

    def test_equivalent_native_derivative_integral_relation_condition_labels(self):
        # Compare equivalent notation against an independently serialized
        # reference pair, not against a hand-stated residual value.
        reference=self.fixture('x','x1')
        fixtures=(
            ('Subs(Derivative(General_below_n(a,b,c,d,e),a),(a,b,c,d,e),(x**2,y,z,w,t))',
             'Derivative[1,0,0,0,0][n1][x1^2,x2,x3,w,t]'),
            ('Integral(x+w,(w,0,1))','Inactive[Integrate][x1+w,{w,0,1}]'),
            ('x>y','x2<x1'),
            ('Eq(x,y)','Inactive[Equal][x2,x1]'),
        )
        def exact_products(records):
            return [r['as_computed_residual'] for r in records if r['kind']=='residual_before_domain']
        for py,wl in fixtures:
            with self.subTest(object=py):
                result=self.fixture(py,wl)
                ref=self.fixture('x<y','x1<x2') if py=='x>y' else reference
                self.assertTrue(exact_products(result))
                self.assertEqual(set(exact_products(ref)),set(exact_products(result)))
        self.label_control(False)
        self.label_control(True)
        list_reference=self.label_measure(('alpha','beta'),('alpha','beta'),True)
        self.assertTrue(all('text_difference' in item for item in self.comparison(list_reference)[0]['comparison']['components']))
        row=next(r for r in c.BINDING_TABLE if r.py.tag=='SphericalSpecialization/outer projected velocity')
        self.equivalent('{domain: Ne(radius,0), object: radius+2}',
                        'ConditionalExpression[r+2,r!=0]','ConditionalExpression[r+2,Unequal[r,0]]',row=row)

    def test_equivalence_integral_alpha_and_closed_endpoint(self):
        self.equivalent('Integral(x+w,(w,0,1))','Inactive[Integrate][x1+w,{w,0,1}]',
                        'Inactive[Integrate][x1+q,{q,0,1}]')
        self.equivalent('Integral(x+w,(w,0,radius),(radius,1,2))',
                        'Inactive[Integrate][x1+w,{w,0,r},{r,1,2}]',
                        'Inactive[Integrate][Inactive[Integrate][x1+w,{w,0,r}],{r,1,2}]')
        py=f'Integral((x+w)*DiracDelta(w-{PLO}),(w,{PLO},{PHI}))'
        wl=f'Inactive[Integrate][(x1+w)*DiracDelta[w-{WLO}],{{w,{WLO},{WHI}}}]'
        changed=f'x1+{WLO}'
        self.equivalent(py,wl,changed)
        original=c.delta_integral
        def half(integral,**kwargs):
            value,ledger=original(integral,**kwargs)
            return value/2,ledger
        with patch.object(c,'delta_integral',side_effect=half):
            with self.assertRaises(AssertionError):
                self.equivalent(py,wl,changed)

    def test_duplicate_end_name_repeat_tag_and_bad_grammar(self):
        with self.assertRaises(c.InputError):
            self.fixture('x','x1',names=c.NAME_TABLE+(c.NAME_TABLE[0],))
        with self.assertRaises(c.InputError):
            c.validate_tables(c.NAME_TABLE,(SYN,replace(SYN,key='other')))
        before=self.fixture('x','x1')
        after=self.fixture('x','x1',extra_wl='fixture wl: x2\n')
        self.moved(before,after)
        with self.assertRaises(c.InputError):
            self.fixture('x','x1 + [')
        with self.assertRaises(c.InputError):
            self.fixture('x','x1\nORPHAN')

    def test_unknown_singular_atom_remains_coverage(self):
        self.mutation(f'DiracDelta(w-{PLO})',f'DiracDelta[w-{WLO}]','DiracDelta[w^2-1]',channel='component')


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
        result=unittest.TextTestRunner(stream=out,verbosity=2,resultclass=ProgressResult,failfast=True).run(suite)
        out.write(json.dumps({'tests_run':result.testsRun,'errors':len(result.errors),'failures':len(result.failures),
                              'runtime_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})+'\n')
    return 0 if result.wasSuccessful() else 1


if __name__=='__main__':
    raise SystemExit(main())
