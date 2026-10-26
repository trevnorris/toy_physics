#!/usr/bin/env python3
"""Synthetic native-stream controls. No measured transcript is opened.

Every mutation asserts a change in computed residual, samples, coverage or
accounting. Targeted disabling of comparison must invalidate that assertion.
Coverage controls are additionally run with each named reduction on its own.
"""
import argparse
from dataclasses import replace
from contextlib import contextmanager
from functools import lru_cache
from io import StringIO
import inspect
import hashlib
import json
import multiprocessing
import os
import tempfile
import importlib.util
from pathlib import Path
import resource
import time
import traceback
import sympy as s
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational
import S15_cross_engine_comparator as c

ROOT=Path('/var/projects/toy_physics/_scratch/s15/cmp_work/fixtures')
ROOT.mkdir(parents=True,exist_ok=True)
FAST={'forms':c.FORMS[:1],'strata':('generic',),'points':(0,)}
COUNTS={}
ABLATIONS=[]
ACTIVE_CONTROL=''
x=s.symbols('x_1 x_2 x_3 t'); F=s.Function('W_1')(*x)
MAP={r.py:r.wl for r in c.NAME_TABLE}


def serial(v,engine):
    if isinstance(v,str):return v
    if isinstance(v,(int,float)):v=s.sympify(v)
    if isinstance(v,dict):
        if set(v)=={'lhs','rhs'}:
            a,b=serial(v['lhs'],engine),serial(v['rhs'],engine)
            return f'Eq({a}, {b})' if engine=='py' else f'{a} == {b}'
        return ('{' if engine=='py' else '<|')+', '.join((repr(k)+': '+serial(q,engine)) if engine=='py' else (json.dumps(k)+' -> '+serial(q,engine)) for k,q in v.items())+('}' if engine=='py' else '|>')
    if isinstance(v,(list,tuple)):
        return ('[' if engine=='py' else '{')+', '.join(serial(q,engine) for q in v)+(']' if engine=='py' else '}')
    if engine=='py':return s.sstr(v)
    if v==s.pi:return 'Pi'
    if v in (s.true,s.false):return 'True' if v==s.true else 'False'
    if isinstance(v,s.Symbol):return MAP.get(v.name,v.name)
    if isinstance(v,s.Number):return str(v)
    sig=c.partial_signature(v)
    if sig is not None:
        head,orders,points=sig
        return 'Derivative['+','.join(map(str,orders))+']['+MAP.get(head.__name__,head.__name__)+']['+','.join(serial(q,engine) for q in points)+']'
    if isinstance(v,AppliedUndef):return MAP.get(v.func.__name__,v.func.__name__)+'['+','.join(serial(q,engine) for q in v.args)+']'
    if isinstance(v,s.Integral):return 'Inactive[Integrate]['+serial(v.function,engine)+','+','.join(serial(tuple(lim),engine) for lim in v.limits)+']'
    if isinstance(v,s.Add):return '('+' + '.join(serial(q,engine) for q in v.args)+')'
    if isinstance(v,s.Mul):return '('+'*'.join(serial(q,engine) for q in v.args)+')'
    if isinstance(v,s.Pow):return '('+serial(v.base,engine)+')^('+serial(v.exp,engine)+')'
    if isinstance(v,Relational):return '('+serial(v.lhs,engine)+v.rel_op+serial(v.rhs,engine)+')'
    heads={s.sin:'Sin',s.cos:'Cos',s.exp:'Exp',s.Abs:'Abs',s.And:'And',s.Or:'Or'}
    if v.func in heads:return heads[v.func]+'['+','.join(serial(q,engine) for q in v.args)+']'
    raise ValueError('fixture serializer '+str(v))


def basic_row(kind='value'):
    return c.Binding('synthetic','synthetic',c.End('object'),c.End('object'),'synthetic object','fixture constructor','fixture constructor',kind)


def fixture_run(py,wl,rows=None,names=c.NAME_TABLE,options=None,workers=1,submission_order="forward"):
    rows=(basic_row(),) if rows is None else rows
    paths=[]
    for engine,objects in (('py',py),('wl',wl)):
        items=objects.items() if isinstance(objects,dict) else objects
        txt='\n'.join(k+': '+serial(v,engine) for k,v in items)+'\n'
        path=ROOT/(hashlib.sha256((engine+txt).encode()).hexdigest()+'.'+engine)
        if not path.exists():path.write_text(txt)
        paths.append(path)
    out=StringIO()
    with tempfile.TemporaryDirectory(prefix='schedule_',dir=ROOT) as directory:
        c.run(*paths,out,names=names,bindings=rows,options=FAST if options is None else options,progress=False,metadata=False,workers=workers,spool=Path(directory)/'bindings',submission_order=submission_order)
    result=Fixture(json.loads(line) for line in out.getvalue().splitlines())
    result.recipe=(tuple(paths),tuple(rows),tuple(names),dict(FAST if options is None else options),tuple(c.SLOT_TABLE),workers,submission_order)
    return result


class Fixture(list):
    """A native-stream result with the exact inputs needed for fault replay."""


# Faults are confined to one production operation and one binding. The source
# on disk is never modified. Each replacement must identify exactly one site.
# In particular, no fault uses disable_comparison or removes all comparisons.
FAULTS={
    'exact':('compare_value','residual=exact_difference(aa,bb)','residual=exact_difference(aa,aa)'),
    'sample_py':('compare_value','for v in (aa,bb)]','for v in (bb,bb)]'),
    'sample_wl':('compare_value','for v in (aa,bb)]','for v in (aa,aa)]'),
    'sample_offset':('compare_value',"rec.update(py=str(vals[0])", "vals[0]=s.Not(vals[0]) if vals[0] in (s.true,s.false) else vals[0]+1\n                rec.update(py=str(vals[0])"),
    'text':('compare_value',"if aa!=bb else []","if False else []"),
    'text_equality':('compare_value',"if aa!=bb else []","if True else []"),
    'dimension':('compare_value','ca,cb=Counter(aa),Counter(bb)','ca,cb=Counter(aa),Counter(aa)'),
    'dimension_sample':('compare_value',"'py':sorted(str(v.subs(env)) for v in aa)","'py':sorted(str(v.subs(env)) for v in bb)"),
    'dimension_offset':('compare_value',"'py':sorted(str(v.subs(env)) for v in aa)","'py':sorted(str(2*v.subs(env)) for v in aa)"),
    'formation':('compare_binding',"computed={'not_formed':str(exc)}","computed={}"),
    'domain':('compare_value','domains=[evaluate_cached(v,form,stratum,index,opts) for v in conditions]','domains=[s.true for v in conditions]'),
    'leaf_credit':('compare_binding',"if not options.get('omit_leaf_credit'):","if False:"),
}


@lru_cache(None)
def altered_function(fault):
    name,old,new=FAULTS[fault]
    source=inspect.getsource(getattr(c,name))
    assert source.count(old)==1, 'fault site absent or ambiguous: '+fault
    namespace={}
    exec(compile(source.replace(old,new),'<s15-test-fault:'+fault+'>','exec'),c.__dict__,namespace)
    return namespace[name]


@contextmanager
def operation_fault(fault,binding):
    name=FAULTS[fault][0];original=getattr(c,name);mutant=altered_function(fault)
    def selected(*args,**kwargs):
        row=args[2] if name=='compare_value' else args[0]
        return (mutant if row.key==binding else original)(*args,**kwargs)
    setattr(c,name,selected)
    try:yield
    finally:setattr(c,name,original)


def replay(records,options=None):
    paths,rows,names,opts,slots,workers,order=records.recipe
    out=StringIO();old=c.SLOT_TABLE
    try:
        c.SLOT_TABLE=slots
        with tempfile.TemporaryDirectory(prefix='ablation_',dir=ROOT) as directory:
            c.run(*paths,out,names=names,bindings=rows,options={**opts,**(options or {})},progress=False,metadata=False,workers=workers,spool=Path(directory)/'bindings',submission_order=order)
    finally:c.SLOT_TABLE=old
    return list(map(json.loads,out.getvalue().splitlines()))


def rejected_assertion(check,message):
    # Only the intended assertion is evidence. Grammar errors, missing keys,
    # an inactive baseline, or another assertion are test-harness bugs.
    try:check()
    except AssertionError as exc:
        assert str(exc)==message, 'fault rejected at the wrong assertion: '+str(exc)
    else:raise AssertionError('targeted control survived: '+message)


def require(value,message):
    assert value,message


def audit_fault(fault,binding,fixtures,check,message):
    check(*fixtures)
    with operation_fault(fault,binding):mutants=[replay(q) for q in fixtures]
    # Every non-target binding retains its full computation, not just operands.
    for ordinary,mutated in zip(fixtures,mutants):
        assert [q for q in comps(ordinary) if q['binding']!=binding]==[q for q in comps(mutated) if q['binding']!=binding], 'fault changed another binding'
        for a,b in zip(comps(ordinary),comps(mutated)):
            if a['binding']!=binding:continue
            if fault in ('exact','text','dimension'):
                assert a.get('samples')==b.get('samples') and a.get('dimension_evaluation')==b.get('dimension_evaluation'), 'residual fault disabled evaluation'
            if fault.startswith('sample') or fault in ('domain','dimension_sample','dimension_offset'):
                assert a.get('exact_residual')==b.get('exact_residual') and a.get('dimension_difference')==b.get('dimension_difference'), 'evaluation fault changed exact comparison'
            if fault in ('sample_py','sample_wl'):
                other='wl' if fault=='sample_py' else 'py'
                assert [v.get(other) for v in a['samples']]==[v.get(other) for v in b['samples']], 'sample fault changed the other operand'
    rejected_assertion(lambda:check(*mutants),message)
    ABLATIONS.append({'control':ACTIVE_CONTROL,'instance':len(ABLATIONS)+1,'binding':binding,'fault':fault,
        'site':FAULTS[fault][0],'removed':FAULTS[fault][1],'installed':FAULTS[fault][2],
        'assertion':message,'fixtures':[[str(p) for p in q.recipe[0]] for q in fixtures],
        'configurations':[{'options':q.recipe[3],
            'name_rows_changed':[c.asdict(r) for r in q.recipe[2] if r not in c.NAME_TABLE],
            'slot_rows':[c.asdict(r) for r in q.recipe[4]]} for q in fixtures],
        'bindings':[[c.asdict(r) for r in q.recipe[1] if r.key==binding] for q in fixtures]})


def comparison(records,binding):return next(q for q in comps(records) if q['binding']==binding)


def equivalent_samples(records,message='native-equivalent samples changed'):
    for row in comps(records):
        if not row.get('samples'):continue
        key=row['binding']
        def check(q):
            ss=comparison(q,key)['samples']
            require(bool(ss) and all(v.get('coverage')=='evaluated' and v['py']==v['wl'] for v in ss),message)
        audit_fault('sample_offset',key,(records,),check,message)


def custom_fault(name,description,fixtures,check,message,owner,attribute,mutant):
    check(*fixtures)
    original=getattr(owner,attribute)
    try:
        setattr(owner,attribute,mutant)
        altered=[replay(q) for q in fixtures]
    finally:setattr(owner,attribute,original)
    rejected_assertion(lambda:check(*altered),message)
    ABLATIONS.append({'control':ACTIVE_CONTROL,'instance':len(ABLATIONS)+1,'fault':name,
        'site':attribute,'installed':description,'assertion':message,
        'fixtures':[[str(p) for p in q.recipe[0]] for q in fixtures]})


def record_order_equivalent(first,reordered):
    message='record order changed comparison'
    def check(a,b):require(computed(a)==computed(b),message)
    original=c.Stream
    class OrdinalRecords(original):
        def __init__(self,*args,**kwargs):
            super().__init__(*args,**kwargs)
            if self.engine!='py':return
            for prefix in ('dimension_argument/name@','dimension_argument/operand_a@','dimension_argument/operand_b@','dimension_argument/residuals@','partition/unlabelled_coefficient@'):
                keys=[k for k in self.index if k.startswith(prefix)]
                values=[self.index[k] for k in keys]
                for key,value in zip(sorted(keys),values):self.index[key]=value
    custom_fault('record_ordinal','assign repeated records to sorted identities in stream order',
        (first,reordered),check,message,c,'Stream',OrdinalRecords)


def duplicate_control(invoke,raise_text):
    def check():
        try:invoke()
        except c.InputError:return
        raise AssertionError('duplicate accepted')
    check()
    original=c.validate_tables;source=inspect.getsource(original)
    assert source.count(raise_text)==1,'duplicate validator fault site absent'
    namespace={}
    exec(compile(source.replace(raise_text,'pass'),'<s15-duplicate-fault>','exec'),c.__dict__,namespace)
    try:
        c.validate_tables=namespace['validate_tables']
        rejected_assertion(check,'duplicate accepted')
    finally:c.validate_tables=original
    ABLATIONS.append({'control':ACTIVE_CONTROL,'instance':len(ABLATIONS)+1,'fault':'duplicate_validation',
        'site':'validate_tables','removed':raise_text,'installed':'pass','assertion':'duplicate accepted'})


def comps(records):return [q for q in records if q['kind']=='comparison']

def computed(records):
    # Original operands are excluded; only computation can detect a mutation.
    return [{k:v for k,v in q.items() if k not in ('py_value','wl_value','py_alias_resolved','wl_alias_resolved','kind')} for q in comps(records)]


def samples(records,engine='py'):
    return [q[engine] for r in comps(records) for q in r.get('samples',[]) if q.get('coverage')=='evaluated']


def exacts(records):return [q.get('exact_residual',q.get('dimension_difference',q.get('text_difference'))) for q in comps(records)]


def detects(base,changed,evaluation=False,engine='py',exact=False):
    targets=[(a,b) for a,b in zip(comps(base),comps(changed)) if computed([a])!=computed([b])]
    assert computed(base)!=computed(changed),'mutation invisible to comparison'
    assert all('not_compared' not in a and 'not_formed' not in a for a,b in targets), 'target baseline comparison inactive'
    other='wl' if engine=='py' else 'py'
    a=[q[other] for q in base if q['kind']=='operands'];b=[q[other] for q in changed if q['kind']=='operands']
    assert a==b,'unchanged operand changed'
    if samples(base,other) and samples(changed,other):assert samples(base,other)==samples(changed,other),'other operand evaluation changed'
    if evaluation:assert samples(base,engine) and samples(changed,engine) and samples(base,engine)!=samples(changed,engine),'mutation invisible to samples'
    if exact:assert exacts(base)!=exacts(changed),'mutation invisible to exact residual'
    exercised=0
    for a,b in targets:
        key=a['binding']
        if 'not_formed' in b:
            message='required formation diagnostic absent'
            def check(p,q):require(bool(comparison(q,key).get('not_formed')),message)
            audit_fault('formation',key,(base,changed),check,message);exercised+=1
            continue
        for field,fault in (('exact_residual','exact'),('text_difference','text'),('dimension_difference','dimension')):
            if a.get(field)==b.get(field):continue
            message='mutation invisible to '+field
            def check(p,q):require(comparison(p,key).get(field)!=comparison(q,key).get(field),message)
            audit_fault(fault,key,(base,changed),check,message);exercised+=1
        if a.get('dimension_evaluation')!=b.get('dimension_evaluation'):
            message='mutation invisible to dimension evaluation'
            def check(p,q):require(comparison(p,key).get('dimension_evaluation')!=comparison(q,key).get('dimension_evaluation'),message)
            audit_fault('dimension_sample',key,(base,changed),check,message);exercised+=1
        av=[q.get(engine) for q in a.get('samples',[]) if q.get('coverage')=='evaluated']
        bv=[q.get(engine) for q in b.get('samples',[]) if q.get('coverage')=='evaluated']
        if av and bv and av!=bv:
            message='mutation invisible to samples of '+engine
            def check(p,q):
                def values(records):return [v.get(engine) for v in comparison(records,key).get('samples',[]) if v.get('coverage')=='evaluated']
                require(values(p)!=values(q),message)
            audit_fault('sample_'+engine,key,(base,changed),check,message);exercised+=1
        if any(v.get('coverage')=='unresolved' for v in b.get('samples',[])):
            message='malformed domain was assumed satisfied'
            def check(p,q):require(any(v.get('coverage')=='unresolved' for v in comparison(q,key).get('samples',[])),message)
            audit_fault('domain',key,(base,changed),check,message);exercised+=1
    assert exercised,'mutation has no targeted comparison assertion'


def basic_mutations():
    base=fixture_run({'object':F},{'object':F})
    for value in (F+1,F+x[0]**3,-F,F.subs(x[0],x[0]+1),s.Integer(7)):
        detects(base,fixture_run({'object':value},{'object':F}),True,exact=True)
    detects(base,fixture_run({'object':'W_1(x_1,x_2,t)'},{'object':F}))


def labels_and_conditions():
    for tag in ('REGULARITY','pressure_expansion/order_and_condition','inner/pressure/component_1/first_order_domain','total_velocity/comparison','region/component_1/balance/reduction'):
        row=replace(basic_row('label'),py=c.End(tag),wl=c.End(tag))
        base=fixture_run({tag:'alpha condition'},{tag:'"alpha condition"'},(row,))
        detects(base,fixture_run({tag:'beta condition'},{tag:'"alpha condition"'},(row,)))
        detects(base,fixture_run({tag:'alpha condition'},{tag:'"beta condition"'},(row,)),engine='wl')
    label_row=replace(basic_row('label'),py=c.End('REGULARITY'),wl=c.End('REGULARITY'))
    literal=fixture_run({'REGULARITY':'7'},{'REGULARITY':'"7"'},(label_row,))
    detects(literal,fixture_run({'REGULARITY':'7'},{'REGULARITY':s.Integer(7)},(label_row,)),engine='wl')
    row=basic_row('condition')
    py='(r_in(t)>0) & (R(t)>r_in(t))';wl='rout[t]>rin[t]>0'
    base=fixture_run({'object':py},{'object':wl},(row,))
    detects(base,fixture_run({'object':'(r_in(t)>2) & (R(t)>r_in(t))'},{'object':wl},(row,)),True,exact=True)
    detects(base,fixture_run({'object':'Eq(r_in(t),0)'},{'object':wl},(row,)),True,exact=True)
    domain='DOMAIN_FOR_EVERY_OBJECT'
    valid_py={'object':F,domain:[s.Gt(x[0],0,evaluate=False),repr('domain label')]}
    valid_wl={'object':F,domain:[s.Gt(x[0],0,evaluate=False),'"domain label"']}
    valid=fixture_run(valid_py,valid_wl)
    malformed=fixture_run({**valid_py,domain:[repr('not a relation'),repr('domain label')]},valid_wl)
    detects(valid,malformed)
    assert any(q.get('coverage')=='unresolved' for r in comps(malformed) for q in r.get('samples',[])), 'malformed condition was assumed satisfied'


def equivalences():
    cases=[('Subs(Derivative(W_1(u,v,w,z),u),(u,v,w,z),(x_1+1,x_2,x_3,t))','Derivative[1,0,0,0][ww[1]][x+1,y,z,t]'),
           ('Derivative(p_br(x_1),x_1)',"pb'[x]"),
           ('(r_in(t)>0) & (R(t)>r_in(t)) & (X_1_1(t)>R(t))','xx[1,1][t]>rout[t]>rin[t]>0')]
    for py,wl in cases:
        row=basic_row('condition' if py.startswith('(r_') else 'value')
        records=fixture_run({'object':py},{'object':wl},(row,))
        equivalent_samples(records)
    # Same bindings with the alternate native representation must have the same
    # computed comparison as the independently serialized canonical expression.
    for payload in ([F,2*F],{'a':F,'b':2*F},{'lhs':F,'rhs':2*F},[[F,2*F],[3*F,4*F]]):
        if isinstance(payload,dict):paths=tuple((k,) for k in payload)
        elif isinstance(payload[0],list):paths=tuple((i,j) for i in range(2) for j in range(2))
        else:paths=((0,),(1,))
        rows=tuple(replace(basic_row(),key=str(path),py=c.End('object',path),wl=c.End('object',path)) for path in paths)
        records=fixture_run({'object':payload},{'object':payload},rows)
        equivalent_samples(records)
    row=replace(basic_row(),py=c.End('object',(1,)),wl=c.End('object',(1,)))
    records=fixture_run({'object':'Matrix([\n[W_1(x_1,x_2,x_3,t)],\n[2*W_1(x_1,x_2,x_3,t)]])'},{'object':[F,2*F]},(row,))
    equivalent_samples(records)
    row=basic_row('dimension');r=fixture_run({'object':[[1,2,-2]]},{'object':[[1,2,-2]]},(row,))
    message='dimension notation changed value'
    def check(q):require(all(v['dimension_evaluation']['py']==v['dimension_evaluation']['wl'] for v in comps(q)),message)
    audit_fault('dimension_offset','synthetic',(r,),check,message)
    row=basic_row('label');row=replace(row,py=c.End('REGULARITY'),wl=c.End('REGULARITY'))
    a=fixture_run({'REGULARITY':'same label'},{'REGULARITY':'"same label"'},(row,))
    b=fixture_run({'REGULARITY':'other label'},{'REGULARITY':'"other label"'},(row,))
    message='label notation changed comparison'
    def check(a,b):require([q['text_difference'] for q in comps(a)]==[q['text_difference'] for q in comps(b)],message)
    audit_fault('text_equality','synthetic',(a,b),check,message)
    th,ph=s.symbols('theta phi');r=s.Function('r_in')(x[3]);trig=s.sin(th)**2*s.cos(ph)**2+s.sin(th)**2*s.sin(ph)**2+s.cos(th)**2
    for integrand in (F, s.diff(F,x[1])):
        first=s.Integral(integrand*r**2*s.sin(th),(th,0,s.pi),(ph,0,2*s.pi))
        second=s.Integral(integrand*r**2*s.sin(th)*trig,(th,0,s.pi),(ph,0,2*s.pi))
        a=fixture_run({'object':first},{'object':first});b=fixture_run({'object':second},{'object':first})
        assert samples(a) and samples(a)==samples(b),'surface integral representation changed value'
        equivalent_samples(a,'surface integral representation changed value')
        equivalent_samples(b,'surface integral representation changed value')


def long_native_payload():
    py=' + '.join(f'{i}*W_1(x_1,x_2,x_3,t)' for i in range(1,1501))
    wl=' + '.join(f'{i}*ww[1][x,y,z,t]' for i in range(1,1501))
    base=fixture_run({'object':py},{'object':wl})
    equivalent_samples(base)
    detects(base,fixture_run({'object':py+' + x_1**2'},{'object':wl}),True,exact=True)


def derivative_points():
    # Each live function, every derivative slot and evaluation-point argument,
    # evaluated inside Subs, at a composite point, and within both integral types.
    th,ph,r=s.symbols('theta phi s')
    for name in c.NAME_TABLE:
        if not name.arity:continue
        n=name.arity;slots=c.FORMAL[:n]
        for derivative_slot in range(n):
            for point_slot in range(n):
                for context in ('substitution','composite','surface','volume'):
                    points=list(x if n==4 else (x[0],x[3]) if n==2 else (x[3],))
                    if context=='composite':points[point_slot]+=x[0]**2+s.Rational(2,9)
                    if context in ('surface','volume'):
                        radius=s.Function('r_in')(x[3]) if context=='surface' else r
                        if n>=2:points[0]=2+radius*s.sin(th)*s.cos(ph)
                    fun=s.Function(name.py)
                    d=s.Derivative(fun(*slots),slots[derivative_slot],evaluate=False)
                    original=s.Subs(d,slots,points)
                    shifted=points.copy();shifted[point_slot]+=s.Rational(1,7)
                    mutation=s.Subs(d,slots,shifted)
                    def wrap(v):
                        if context=='surface':return s.Integral(v*s.sin(th),(th,0,s.pi),(ph,0,2*s.pi))
                        if context=='volume':return s.Integral(v*r**2*s.sin(th),(r,1,2),(th,0,s.pi),(ph,0,2*s.pi))
                        return v
                    a,b=wrap(original),wrap(mutation)
                    base=fixture_run({'object':a},{'object':a})
                    detects(base,fixture_run({'object':b},{'object':a}),True,exact=True)
    baseline=s.diff(F,x[0]);base=fixture_run({'object':baseline},{'object':baseline})
    for mutation in (s.diff(F,x[1]),s.diff(F,x[0],2)):
        detects(base,fixture_run({'object':mutation},{'object':baseline}),True,exact=True)


def mixed_second_derivative_points():
    th,ph,r=s.symbols('theta phi s')
    for name,n in (('W_1',4),('rho_1',2),('u_1',2)):
        slots=c.FORMAL[:n];fun=s.Function(name)
        for first in range(n):
            for second in range(n):
                for changed_arg in range(n):
                    for context in ('substitution','composite','surface','volume'):
                        points=list(x if n==4 else (x[0],x[3]))
                        if context=='composite':points[changed_arg]+=x[0]**2+s.Rational(2,9)
                        if context in ('surface','volume'):points[0]=2+r*s.sin(th)*s.cos(ph) if context=='volume' else 2+s.sin(th)*s.cos(ph)
                        derivative=s.Derivative(fun(*slots),slots[first],slots[second],evaluate=False)
                        original=s.Subs(derivative,slots,points)
                        shifted=points.copy();shifted[changed_arg]+=s.Rational(1,7)
                        altered=s.Subs(derivative,slots,shifted)
                        def wrap(v):
                            if context=='surface':return s.Integral(v*s.sin(th),(th,0,s.pi),(ph,0,2*s.pi))
                            if context=='volume':return s.Integral(v*r**2*s.sin(th),(r,1,2),(th,0,s.pi),(ph,0,2*s.pi))
                            return v
                        a,b=wrap(original),wrap(altered)
                        options={**FAST,'forms':c.FORMS[2:]}
                        detects(fixture_run({'object':a},{'object':a},options=options),fixture_run({'object':b},{'object':a},options=options),True,exact=True)


def partial_time_semantics():
    # Independent numeric oracle: differentiate the printed test polynomial in
    # its formal time argument BEFORE inserting the moving chart. This operand
    # never passes through the live-function derivative translator being tested.
    th,ph,r=s.symbols('theta phi s');slots=c.FORMAL;tvar=x[3]
    polynomial=c.profile('W_1',4,c.FORMS[0],'generic',0)
    for context in ('composite','surface','volume'):
        radius=r if context=='volume' else s.Function('r_in')(tvar)
        moving=(s.Function('X_2_1')(tvar)+radius*s.sin(th)*s.cos(ph),
                s.Function('X_2_2')(tvar)+radius*s.sin(th)*s.sin(ph),
                s.Function('X_2_3')(tvar)+radius*s.cos(th),tvar)
        partial=s.Subs(s.diff(s.Function('W_1')(*slots),slots[3]),slots,moving)
        total=s.diff(s.Function('W_1')(*moving),tvar)
        expected=s.diff(polynomial,slots[3]).xreplace(dict(zip(slots,moving)))
        def wrap(v):
            if context=='surface':return s.Integral(v*radius**2*s.sin(th),(th,0,s.pi),(ph,0,2*s.pi))
            if context=='volume':return s.Integral(v*r**2*s.sin(th),(r,s.Function('r_in')(tvar),s.Function('R')(tvar)),(th,0,s.pi),(ph,0,2*s.pi))
            return v
        options={**FAST,'strata':c.STRATA,'points':tuple(range(len(c.POINTS)))}
        base=fixture_run({'object':wrap(partial)},{'object':wrap(expected)},options=options)
        equivalent_samples(base,'partial time became chart time')
        changed=fixture_run({'object':wrap(total)},{'object':wrap(expected)},options=options)
        detects(base,changed,True,exact=True)


def integral_changes():
    th,ph,r=s.symbols('theta phi s');t=x[3]
    for mode in ('surface','volume'):
        radius=s.Function('r_in')(t)
        limits=[(th,0,s.pi),(ph,0,2*s.pi)]
        if mode=='volume':limits.insert(0,(r,radius,s.Function('R')(t)))
        measure=(r**2 if mode=='volume' else radius**2)*s.sin(th)
        value=s.Integral(F*measure,*limits)
        base=fixture_run({'object':value},{'object':value})
        # Surface changed by area/radius, region changed by live upper bound;
        # integrand field is kept in the first mutation.
        moved=value.subs(s.Function('r_in')(t),s.Function('R')(t)) if mode=='surface' else s.Integral(F*measure,(r,radius,s.Function('R')(t)+1),*limits[1:])
        altered=s.Integral((F+x[0]**2)*measure,*limits)
        for mutation in (moved,altered):detects(base,fixture_run({'object':mutation},{'object':value}),True,exact=True)


def schema_payload(row,engine,value):
    end=getattr(row,engine);shape=c.ROOT_SHAPES.get((engine,end.tag))
    def build(sh):
        if sh=='predicate':return s.Gt(x[0],0,evaluate=False)
        if sh=='v':return F
        if isinstance(sh,dict):return {k:build(v) for k,v in sh.items()}
        return [build(v) for v in sh]
    if shape is not None:result=build(shape)
    else:
        result=value
        for key in reversed(end.path):result={key:result} if isinstance(key,str) else [F]*key+[result]
    if end.path:
        at=result
        for k in end.path[:-1]:at=at[k]
        at[end.path[-1]]=value
    elif shape is None or shape=='v':result=value
    return result


def binding_repoints():
    for rule,row in {r.rule:r for r in c.BINDING_TABLE}.items():
        if row.kind=='label':
            p,w,bad=(repr('alpha') if row.py.path else 'alpha'),'"alpha"',(repr('beta') if row.py.path else 'beta')
        elif row.kind=='basis':p,w,bad="('Density','Length','Time')",'{Dens,Length,Time}',"('Time','Length','Density')"
        elif row.kind=='condition':p,w,bad='x_1>0','x>0','x_1>3'
        elif row.kind=='dimension':p,w,bad=[[1,0,0]],[[1,0,0]],[[0,1,0]]
        elif row.kind=='dimension_field':p,w,bad=[1,0,0],[1,0,0],[0,1,0]
        else:p,w,bad=F,F,F+x[0]**2
        pp=schema_payload(row,'py',p);ww=schema_payload(row,'wl',w);alt=schema_payload(row,'py',bad)
        py={row.py.tag:pp,'alternate':alt};wl={row.wl.tag:ww}
        if row.kind=='label':
            # Native label path includes tag-dependent unquoted reading.
            alternate='alternate' if row.py.path else 'REGULARITY'
            if row.py.tag==alternate:alternate='CHOICES'
            py={row.py.tag:pp,alternate:alt}
        else:alternate='alternate'
        base=fixture_run(py,wl,(row,));changed=fixture_run(py,wl,(replace(row,py=c.End(alternate,row.py.path)),))
        detects(base,changed,evaluation=row.kind in ('value','condition','dimension_field'),exact=True)
        if row.kind=='dimension':
            assert [q.get('dimension_evaluation') for q in comps(base)]!=[q.get('dimension_evaluation') for q in comps(changed)]
    # Semantic mistakes in the enumerated surface, component and group axes.
    for source,target in [('inner/carried/component_1','outer/carried/component_1'),
                          ('inner/carried/component_1','inner/carried/component_2'),
                          ('inner/carried/component_1/labels_drain1','inner/carried/component_1/labels_drain2')]:
        row=next(r for r in c.BINDING_TABLE if r.py.tag==source)
        py={source:F,target:2*F};wl={row.wl.tag:F}
        detects(fixture_run(py,wl,(row,)),fixture_run(py,wl,(replace(row,py=c.End(target)),)),True,exact=True)


def name_repoints():
    for row in c.NAME_TABLE:
        candidates=[r for r in c.NAME_TABLE if r!=row and r.arity==row.arity]
        preferred={'Pi_o_1_2':'Pi_o_2_1','Pi_o_2_1':'Pi_o_1_2','u_1':'u_2','u_2':'u_1','rho_1':'rho_2','rho_2':'rho_1','r_in':'R','R':'r_in'}
        partner=next((r for r in candidates if r.py==preferred.get(row.py)),candidates[0])
        names=tuple(replace(r,wl=partner.wl) if r==row else replace(r,wl=row.wl) if r==partner else r for r in c.NAME_TABLE)
        args=x if row.arity==4 else (x[0],x[3]) if row.arity==2 else (x[3],)
        value=s.Function(row.py)(*args) if row.arity else s.Symbol(row.py)
        base=fixture_run({'object':value},{'object':value},options={**FAST,'points':tuple(range(len(c.POINTS)))})
        changed=fixture_run({'object':value},{'object':value},names=names,options={**FAST,'points':tuple(range(len(c.POINTS)))})
        detects(base,changed,True,engine='wl',exact=True)


def bound_name_repoints():
    original=c.SLOT_TABLE
    for row in original:
        partner=original[(row.py_index+1)%len(original)]
        value=s.diff(F,x[row.py_index])
        base=fixture_run({'object':value},{'object':value})
        try:
            c.SLOT_TABLE=tuple(replace(q,wl_index=partner.wl_index) if q==row else replace(q,wl_index=row.wl_index) if q==partner else q for q in original)
            changed=fixture_run({'object':value},{'object':value})
            detects(base,changed,True,engine='wl',exact=True)
        finally:c.SLOT_TABLE=original
    try:
        c.SLOT_TABLE=original+(original[0],)
        duplicate_control(lambda:c.validate_tables(c.NAME_TABLE,(basic_row(),)),
            "raise InputError('duplicate bound-name row '+field)")
    finally:c.SLOT_TABLE=original


def pieces_groups_and_components():
    for kind in ('piece','group'):
        rows=tuple(replace(basic_row(),key=kind+str(i),py=c.End(kind+str(i)),wl=c.End(kind+str(i))) for i in range(3))
        a={kind+'0':F,kind+'1':2*F,kind+'2':3*F}
        base=fixture_run(a,a,rows)
        moved={**a,kind+'0':F+1,kind+'1':2*F-1}
        changed=fixture_run(moved,a,rows);detects(base,changed,True)
        before,after=comps(base),comps(changed)
        assert before[0]['exact_residual']!=after[0]['exact_residual'] and before[1]['exact_residual']!=after[1]['exact_residual']
        assert before[2]['exact_residual']==after[2]['exact_residual'],'unchanged total moved'
        for i in (0,1):
            one={**a,kind+str(i):a[kind+str(i)]+x[0]**2}
            changed_one=fixture_run(one,a,rows)
            detects(base,changed_one,True,exact=True)
            assert comps(base)[i]['exact_residual']!=comps(changed_one)[i]['exact_residual'], 'targeted piece/group did not move'
            for j in range(3):
                if j!=i: assert comps(base)[j]['exact_residual']==comps(changed_one)[j]['exact_residual'], 'untouched piece/group moved'
    # Every source component-correspondence place (tag and semantic path) is
    # exercised, including matrices, equations, nested surface dictionaries and
    # field signature tuples split into two Wolfram associations.
    for row in c.BINDING_TABLE:
        if not row.py.path or row.kind not in ('value','dimension_field'):continue
        p=schema_payload(row,'py',F if row.kind=='value' else [1,0,0])
        w=schema_payload(row,'wl',F if row.kind=='value' else [1,0,0])
        q=schema_payload(row,'py',2*F if row.kind=='value' else [0,1,0])
        a=fixture_run({row.py.tag:p,'alternate':q},{row.wl.tag:w},(row,))
        b=fixture_run({row.py.tag:p,'alternate':q},{row.wl.tag:w},(replace(row,py=c.End('alternate',row.py.path)),))
        detects(a,b,True,exact=True)
    # Repoint selectors themselves in each extraction mechanism.
    fixtures=[([F,2*F],(0,),(1,)),({'a':F,'b':2*F},('a',),('b',)),
              ({'lhs':F,'rhs':2*F},('lhs',),('rhs',)),([[F],[2*F]],(0,0),(1,0)),
              ('Matrix([[W_1(x_1,x_2,x_3,t),2*W_1(x_1,x_2,x_3,t)],[3*W_1(x_1,x_2,x_3,t),4*W_1(x_1,x_2,x_3,t)]])',(0,0),(1,0))]
    for payload,path,alt in fixtures:
        wl=[[F,2*F],[3*F,4*F]] if isinstance(payload,str) else payload
        row=replace(basic_row(),py=c.End('object',path),wl=c.End('object',path))
        a=fixture_run({'object':payload},{'object':wl},(row,))
        b=fixture_run({'object':payload},{'object':wl},(replace(row,py=c.End('object',alt)),))
        detects(a,b,True,exact=True)


def packaging_and_accounting():
    rows=tuple(replace(basic_row(),key='split'+str(i),py=c.End('split'+str(i)),wl=c.End('group',(i,))) for i in range(2))
    py={'split0':F,'split1':2*F};wl={'group':[F,2*F]};base=fixture_run(py,wl,rows)
    for i in range(2):
        pp={**py,'split'+str(i):py['split'+str(i)]+1};ww={'group':[F,2*F]};ww['group'][i]+=1
        changed_py=fixture_run(pp,wl,rows);changed_wl=fixture_run(py,ww,rows)
        detects(base,changed_py,True,exact=True)
        detects(base,changed_wl,True,engine='wl',exact=True)
        for changed in (changed_py,changed_wl):
            assert comps(base)[i]['exact_residual']!=comps(changed)[i]['exact_residual'], 'targeted packaged component did not move'
            assert comps(base)[1-i]['exact_residual']==comps(changed)[1-i]['exact_residual'], 'untouched packaged component moved'
    row=next(r for r in c.BINDING_TABLE if r.py.tag=='velocity_profiles_and_remainder' and r.py.path==(0,0))
    intact=[[F,2*F,3*F],[4*F,5*F,6*F],[7*F,8*F,9*F]]
    removed=[[F,2*F,3*F],[4*F,5*F],[7*F,8*F,9*F]]
    a=fixture_run({row.py.tag:intact},{row.wl.tag:intact},(row,))
    b=fixture_run({row.py.tag:removed},{row.wl.tag:intact},(row,))
    detects(a,b);assert any('not_formed' in q for q in comps(b))
    a=fixture_run({'object':F},{'object':F});b=fixture_run({'object':[F]},{'object':F})
    detects(a,b);assert any('structural mismatch' in q.get('not_formed','') for q in comps(b))
    row=replace(basic_row(),py=c.End('object',(0,)),wl=c.End('object',(0,)))
    a=fixture_run({'object':[F]},{'object':[F]},(row,));b=fixture_run({'object':[F,2*F]},{'object':[F]},(row,))
    aa=next(q for q in a if q['kind']=='tag_accounting' and q['engine']=='py');bb=next(q for q in b if q['kind']=='tag_accounting' and q['engine']=='py')
    message='uncompared parsed leaf disappeared'
    def check(a,b):
        aa=next(q for q in a if q['kind']=='tag_accounting' and q['engine']=='py')
        bb=next(q for q in b if q['kind']=='tag_accounting' and q['engine']=='py')
        require(aa['parsed_leaves']!=bb['parsed_leaves'] and aa['actually_compared_leaves']==bb['actually_compared_leaves'] and bb['unpaired_components'],message)
    custom_fault('unpaired_leaf','omit unpaired-component enumeration',(a,b),check,message,c,'unpaired',lambda *args:[])
    custom_fault('parsed_census','replace parsed subtree leaf count with one',(a,b),check,message,c,'census',lambda node:1)
    bad=fixture_run({'object':[F]},{'object':[F]},(row,),options={**FAST,'omit_leaf_credit':True})
    message='omitted leaf credit invisible'
    def check(a,b):require([q for q in a if q['kind']=='binding_accounting']!=[q for q in b if q['kind']=='binding_accounting'],message)
    audit_fault('leaf_credit','synthetic',(a,bad),check,message)
    a=fixture_run({'object':'Matrix([\n[W_1(x_1,x_2,x_3,t)],\n[2*W_1(x_1,x_2,x_3,t)]])'},{'object':[F,2*F]},(row,))
    b=fixture_run({'object':'Matrix([\n[W_1(x_1,x_2,x_3,t)+1],\n[2*W_1(x_1,x_2,x_3,t)]])'},{'object':[F,2*F]},(row,))
    detects(a,b,True,exact=True)


def equation_relation_kind():
    row=next(r for r in c.BINDING_TABLE if r.py.tag=='region/component_1/balance_equation' and r.py.path==('lhs',))
    py={row.py.tag:{'lhs':F,'rhs':2*F}};wl={row.wl.tag:{'lhs':F,'rhs':2*F}}
    base=fixture_run(py,wl,(row,))
    changed=fixture_run(py,{row.wl.tag:'ww[1][x,y,z,t] -> 2*ww[1][x,y,z,t]'},(row,))
    detects(base,changed,engine='wl')
    assert any('equation relation is absent' in q.get('not_formed','') for q in comps(changed))


def repeated_records():
    def record(key,value):return [('dimension_argument/name',repr((key,0,'x_1'))),('dimension_argument/operand_a',[[1,0,0]]),('dimension_argument/operand_b',value),('dimension_argument/residuals',[[0,0,0]])]
    first=record('u_1',[[1,0,0]]);second=record('u_2',[[0,1,0]])
    key=json.dumps(c.data(c.py_parse(first[0][1])),separators=(',',':'))
    row=replace(basic_row('dimension'),py=c.End('dimension_argument/operand_b@'+key),wl=c.End('group',(0,)))
    wl={'group':[[[1,0,0]],[[0,1,0]]]}
    a=fixture_run(first+second,wl,(row,));b=fixture_run(record('u_1',[[0,1,0]])+second,wl,(row,))
    detects(a,b,exact=True)
    detects(a,fixture_run(second,wl,(row,)))
    reordered=fixture_run(second+first,wl,(row,))
    record_order_equivalent(a,reordered)
    def partition(tag,value):return [(tag+'/dimension_pair_residuals',[]),('partition/unlabelled_coefficient',value)]
    row=replace(basic_row(),py=c.End('partition/unlabelled_coefficient@inner/carried/component_1'),wl=c.End('object'))
    first=partition('inner/carried/component_1',F);second=partition('outer/carried/component_1',2*F)
    a=fixture_run(first+second,{'object':F},(row,));b=fixture_run(partition('inner/carried/component_1',F+1)+second,{'object':F},(row,))
    detects(a,b,True,exact=True);detects(a,fixture_run(second,{'object':F},(row,)))
    record_order_equivalent(a,fixture_run(second+first,{'object':F},(row,)))


def supplied_equation_and_alias():
    # Supplied G=rho*V applied only on one side changes the original operand.
    # Independent symbols here stand for the independently printed named terms.
    row=c.NameRow('alias_G','aliasG','synthetic printed alias of G','fixture','fixture',4,'alias_definition')
    names=c.NAME_TABLE+(row,)
    lhs='alias_G(x_1,x_2,x_3,t)';rhs='aliasG[x,y,z,t]'
    py={'object':lhs,'alias_definition':'Eq('+lhs+', W_1(x_1,x_2,x_3,t))'}
    wl={'object':rhs,'alias_definition':rhs+' == ww[1][x,y,z,t]'}
    a=fixture_run(py,wl,names=names)
    altered={**py,'alias_definition':'Eq('+lhs+', W_2(x_1,x_2,x_3,t))'}
    detects(a,fixture_run(altered,wl,names=names),True,exact=True)
    # S15 supplied equations, decision list 55-58: construct their live
    # Cartesian G and Pi on synthetic fields, then apply only the supplied
    # local balance to one operand. No equation is used by the comparator.
    xyz=s.Matrix(x[:3]);centres=[s.Matrix([s.Function(f'X_{k}_{j}')(x[3]) for j in (1,2,3)]) for k in (1,2)]
    distances=[s.sqrt((xyz-centre).dot(xyz-centre)) for centre in centres]
    density=sum(s.Function(f'rho_{k}') (distances[k-1],x[3]) for k in (1,2))+s.Function('varrho')(*x)
    velocity=sum((s.Function(f'u_{k}')(distances[k-1],x[3])*(xyz-centres[k-1])/distances[k-1] for k in (1,2)),s.zeros(3,1))+s.Matrix([s.Function(f'W_{j}')(*x) for j in (1,2,3)])
    gradient=velocity.jacobian(xyz)
    stress=-s.Function('p_br')(density)*s.eye(3)+s.Function('eta')(*x)*(gradient+gradient.T-s.Rational(2,3)*s.eye(3)*s.trace(gradient))+s.Function('zeta')(*x)*s.eye(3)*s.trace(gradient)
    current=density*velocity*velocity.T-stress+s.Matrix(3,3,lambda i,j:s.Function(f'Pi_o_{i+1}_{j+1}')(*x))
    original=s.diff(density*velocity[0],x[3])
    mutation=s.Function('f_1')(*x)-sum(s.diff(current[0,i],x[i]) for i in range(3))
    base=fixture_run({'object':original},{'object':original})
    detects(base,fixture_run({'object':mutation},{'object':original}),True,exact=True)



def coverage_controls():
    def witness(py,wl,reduction,full=None):
        options={**FAST,**(full or {})}
        a=fixture_run({'object':py},{'object':wl},options=options)
        b=fixture_run({'object':py},{'object':wl},options={**options,**reduction})
        message='full coverage missed witness'
        def check(q):require(samples(q) and samples(q)!=samples(q,'wl'),message)
        check(a)
        rejected_assertion(lambda:check(b),message)
        assert exacts(a)==exacts(b),'coverage reduction changed exact comparison'
        ABLATIONS.append({'control':ACTIVE_CONTROL,'instance':len(ABLATIONS)+1,'fault':'coverage_reduction',
            'site':'sample grid or profile options','installed':reduction,'assertion':message,
            'fixtures':[[str(p) for p in a.recipe[0]]]})
        assert samples(b) and samples(b)==samples(b,'wl'),'individual reduction did not miss its witness'
        return a,b
    witness((x[0]-c.POINTS[0][0])*F,s.S.Zero,{'points':(0,)},{'points':(0,1,2)})
    witness(s.diff(F,x[0],3),s.S.Zero,{'forms':c.FORMS[:1]},{'forms':c.FORMS})
    slots=c.FORMAL
    witness(s.Subs(s.diff(s.Function('W_1')(*slots),slots[0])-s.diff(s.Function('W_1')(*slots),slots[1]),slots,[s.Rational(1,3)]*4),s.S.Zero,{'symmetric_arguments':True})
    witness(s.diff(F,x[0],x[1]),s.S.Zero,{'uncoupled_arguments':True})
    for stratum,name in [('X2_velocity_zero','X_2_1'),('radii_velocity_zero','r_in'),('X2_and_radii_velocity_zero','X_2_1'),('X1_velocity_zero','X_1_1')]:
        speed=s.diff(c.profile(name,1,c.FORMS[0],'generic',0),c.FORMAL[0]).subs(c.FORMAL[0],c.POINTS[0][3])
        witness((s.diff(s.Function(name)(x[3]),x[3])-speed)*F,s.S.Zero,{'strata':('generic',)},{'strata':('generic',stratum)})
    th,ph,r=s.symbols('theta phi s')
    for integral in (s.Integral(F*s.sin(th),(th,0,s.pi),(ph,0,2*s.pi)),s.Integral(F*r**2*s.sin(th),(r,1,2),(th,0,s.pi),(ph,0,2*s.pi))):
        witness(integral,s.S.Zero,{'vanishing_form':True})
    witness(F-s.Function('W_2')(*x),s.S.Zero,{'duplicate_functions':True})


def test_functions():
    # Distinctness, dependence and domain are checked on every declared profile
    # in every form/stratum, not inferred from the profile generator's loops.
    for form in c.FORMS:
        for stratum in c.STRATA:
            for index in range(len(c.POINTS)):
                polys={r.py:c.profile(r.py,r.arity,form,stratum,index) for r in c.NAME_TABLE if r.arity}
                assert len(set(polys.values()))==len(polys),'live test functions collided'
                for row in c.NAME_TABLE:
                    if not row.arity:continue
                    p=polys[row.py]
                    for arg in c.FORMAL[:row.arity]:assert s.diff(p,arg)!=s.diff(s.S.One,arg),'constant test argument'
                    if form==c.FORMS[-1]:
                        for i in range(row.arity):
                            for j in range(row.arity):
                                for k in range(row.arity):
                                    assert s.diff(p,c.FORMAL[i],c.FORMAL[j],c.FORMAL[k])!=s.diff(s.S.One,c.FORMAL[i],c.FORMAL[j],c.FORMAL[k]), 'second derivative loses a point argument'
                    if row.arity>1:
                        values=[s.diff(p,arg).subs(dict(zip(c.FORMAL,c.POINTS[index]))) for arg in c.FORMAL[:row.arity]]
                        assert len(set(values))==len(values),'argument derivatives collided'
                        for i in range(row.arity):
                            for j in range(row.arity):
                                if i!=j:assert s.diff(p,c.FORMAL[i],c.FORMAL[j])!=s.diff(s.S.One,c.FORMAL[i],c.FORMAL[j]),'uncoupled derivative'
                def at(name):return polys[name].subs(c.FORMAL[0],c.POINTS[index][3])
                separation=sum((at(f'X_1_{i}')-at(f'X_2_{i}'))**2 for i in (1,2,3))
                assert separation>at('R')**2 and at('R')>at('r_in')>0,'geometry outside domain'
                point_radius_squared=sum((c.POINTS[index][j]-at(f'X_2_{j+1}'))**2 for j in range(3))
                assert at('r_in')**2<point_radius_squared<at('R')**2, 'Cartesian test point outside the shell'


def duplicates():
    row=basic_row()
    duplicate_control(lambda:c.validate_tables(c.NAME_TABLE,(row,replace(row,key='other'))),
        "raise InputError('duplicate binding end '+field+' '+str(end))")
    duplicate_control(lambda:c.validate_tables(c.NAME_TABLE+(c.NAME_TABLE[0],),(row,)),
        "raise InputError('duplicate name row '+field)")
    duplicate_control(lambda:c.validate_tables(c.NAME_TABLE,(row,),('same','same')),
        "raise InputError('duplicate binding rule')")
    original_rules=c.RULE_NAMES
    try:
        c.RULE_NAMES=original_rules+(original_rules[0],)
        duplicate_control(lambda:fixture_run({'object':F},{'object':F}),
            "raise InputError('duplicate binding rule')")
    finally:c.RULE_NAMES=original_rules
    a=fixture_run({'object':F},{'object':F})
    b=fixture_run([('object',F),('object',2*F)],{'object':F})
    detects(a,b)
    assert any('missing or repeated' in q.get('not_formed','') for q in comps(b))


def dag_sampling():
    # Compare against explicit polynomial substitution, through native streams.
    # This oracle does not use the point-first DAG evaluator to build its operand.
    th,ph,r=s.symbols('theta phi s'); slots=c.FORMAL
    field=s.Function('W_1')(*slots)
    moving=(s.Function('X_2_1')(x[3])+r*s.sin(th)*s.cos(ph),
            s.Function('X_2_2')(x[3])+r*s.sin(th)*s.sin(ph),
            s.Function('X_2_3')(x[3])+r*s.cos(th),x[3])
    d=s.Subs(s.diff(field,slots[0],slots[3]),slots,moving)
    objects=(s.diff(F,x[0],x[3])+s.Function('rho_1')(x[0]**2+1,x[3]),
             s.Integral(d*s.sin(th),(th,0,s.pi),(ph,0,2*s.pi)),
             s.Integral(d*r**2*s.sin(th),(r,s.Function('r_in')(x[3]),s.Function('R')(x[3])),(th,0,s.pi),(ph,0,2*s.pi)))
    for form in c.FORMS:
        for stratum in c.STRATA:
            for index in range(len(c.POINTS)):
                options={'forms':(form,),'strata':(stratum,),'points':(index,)}
                for obj in objects:
                    polynomial=c.substituted(obj,form,stratum,index)
                    base=fixture_run({'object':obj},{'object':polynomial},options=options)
                    equivalent_samples(base,'DAG sampling changed a substituted value')
                    detects(base,fixture_run({'object':obj+1},{'object':polynomial},options=options),True,exact=True)


def quadrature():
    r=s.Symbol('s');a=s.Integral(s.cos(r),(r,0,1));b=s.Integral(s.cos(r)+r**2,(r,0,1))
    base=fixture_run({'object':a},{'object':a})
    detects(base,fixture_run({'object':b},{'object':a}),True,exact=True)
    assert c.NUMERICAL_INTEGRALS and all(q.get('working_precisions') and 'error_estimate' in q for q in c.NUMERICAL_INTEGRALS.values())


def parallel_schedule():
    # Whole native streams, extraction, source-cited components, comparisons,
    # exact-first records and leaf accounting must survive worker repartition.
    th,ph,r=s.symbols('theta phi s')
    derivative=s.Subs(s.diff(F,x[0]),x[0],r*s.sin(th)*s.cos(ph))
    surface=s.Integral(derivative*s.sin(th),(th,0,s.pi),(ph,0,2*s.pi))
    volume=s.Integral(derivative*r**2*s.sin(th),(r,1,2),(th,0,s.pi),(ph,0,2*s.pi))
    objects={'scalar':F,'surface':surface,'volume':volume,'siblings':[F,2*F],
             'REGULARITY':'synthetic label','broken':[F,F]}
    wl={**objects,'REGULARITY':'"synthetic label"','broken':[F]}
    rows=tuple(replace(basic_row('label' if tag=='REGULARITY' else 'value'),key=tag,
                       py=c.End(tag),wl=c.End(tag)) for tag in ('scalar','surface','volume','REGULARITY','broken'))
    rows+=tuple(replace(basic_row(),key='sibling'+str(i),py=c.End('siblings',(i,)),wl=c.End('siblings',(i,))) for i in (0,1))
    options={'forms':c.FORMS,'strata':c.STRATA,'points':tuple(range(len(c.POINTS)))}
    def normalized(records):
        return [{k:v for k,v in q.items() if k!='elapsed_seconds'} for q in records if q['kind']!='resources']
    base=fixture_run(objects,wl,rows=rows,options=options)
    assert any(q.get('exact_residual') is not None for q in comps(base)), 'parallel comparison inactive'
    for workers,order in ((2,'forward'),(8,'forward'),(8,'reverse')):
        result=fixture_run(objects,wl,rows=rows,options=options,workers=workers,submission_order=order)
        message='scheduling changed mathematical output or accounting'
        def check(q):require(normalized(q)==normalized(base),message)
        audit_fault('sample_offset','scalar',(result,),check,message)
        changed=fixture_run({**objects,'scalar':F+1},wl,rows=rows,options=options,workers=workers,submission_order=order)
        detects(result,changed,True,exact=True)


def nested_radial_derivative():
    # Build the polynomial reference by cancelling the formal radial factor
    # before inserting a distance or differentiating. The production operand
    # differentiates the live radial function first, leaving roots inside sums.
    r,th,ph=s.symbols('s theta phi');time=x[3]
    dx=r*s.sin(th)*s.cos(ph)-time**2-1
    q=dx**2+(r*s.sin(th)*s.sin(ph)-2*time-3)**2+(r*s.cos(th)-time**3-5)**2
    velocity=dx*s.Function('u_1')(s.sqrt(q),time)/s.sqrt(q)
    derivative=s.diff(velocity,time)
    radial=(r,s.Function('r_in')(time),s.Function('R')(time))
    angles=((th,0,s.pi),(ph,0,2*s.pi));measure=r**2*s.sin(th)
    boundary=s.Function('W_1')(s.sin(th)*s.cos(ph),s.sin(th)*s.sin(ph),s.cos(th),time)*s.sin(th)
    nested=s.Integral(boundary+s.Integral(derivative*measure,radial),*angles)
    for form in c.FORMS:
        formal=c.profile('u_1',2,form,c.STRATA[0],0)
        radial_polynomial=s.Poly(s.cancel(formal/c.FORMAL[0]),c.FORMAL[0])
        assert all(powers[0]%2==0 for powers,_ in radial_polynomial.terms()),'reference requires the specified smooth radial test profile'
        polynomial=dx*s.Add(*(coefficient.subs(c.FORMAL[1],time)*q**(powers[0]//2) for powers,coefficient in radial_polynomial.terms()))
        reference=s.Integral(boundary,*angles)+s.Integral(s.diff(polynomial,time)*measure,radial,*angles)
        for stratum in c.STRATA:
            for index in range(len(c.POINTS)):
                options={'forms':(form,),'strata':(stratum,),'points':(index,)}
                base=fixture_run({'object':nested},{'object':reference},options=options)
                equivalent_samples(base,'nested radial derivative lost its exact polynomial value')
                mutated=s.Integral(boundary+s.Integral((derivative+1)*measure,radial),*angles)
                detects(base,fixture_run({'object':mutated},{'object':reference},options=options),True,exact=True)


def nested_integral_parameters():
    # The inner radial integral must retain the outer angular variables.
    # A separate native-stream operand supplies the flat iterated integral.
    r,th,ph=s.symbols('s theta phi');time=x[3]
    normal=(s.sin(th)*s.cos(ph),s.sin(th)*s.sin(ph),s.cos(th))
    points=tuple(s.Function(f'X_2_{i+1}')(time)+r*normal[i] for i in range(3))+(time,)
    formal=s.Function('W_1')(*c.FORMAL)
    derivative=s.Subs(s.diff(formal,c.FORMAL[0],c.FORMAL[3]),c.FORMAL,points)
    density=s.Function('rho_1')(r,time)
    integrand=(s.Function('W_1')(*points)+derivative+density)*r**2*s.sin(th)
    radial=(r,s.Function('r_in')(time),s.Function('R')(time))
    angles=((th,0,s.pi),(ph,0,2*s.pi))
    boundary=s.Function('W_2')(*(tuple(2*v for v in normal)+(time,)))*s.sin(th)
    nested=s.Integral(boundary+s.Integral(integrand,radial),*angles)
    flat=s.Integral(boundary,*angles)+s.Integral(integrand,radial,*angles)
    options={'forms':c.FORMS,'strata':c.STRATA,'points':tuple(range(len(c.POINTS)))}
    base=fixture_run({'object':nested},{'object':flat},options=options)
    assert samples(base) and len(samples(base))==len(c.FORMS)*len(c.STRATA)*len(c.POINTS)
    equivalent_samples(base,'nested integral lost outer parameters')
    for change in (integrand+s.sin(th),integrand.xreplace({derivative:s.Subs(s.diff(formal,c.FORMAL[0],c.FORMAL[3]),c.FORMAL,(points[0]+1,)+points[1:])})):
        mutated=s.Integral(boundary+s.Integral(change,radial),*angles)
        detects(base,fixture_run({'object':mutated},{'object':flat},options=options),True,exact=True)
    # Retain a free non-angular outer variable, and pi from an inner angular
    # moment, rather than prematurely demanding a closed numerical value.
    for nested,flat in (
        (s.Integral(r+s.Integral(s.Function('W_1')(r,s.cos(th),x[2],time)*s.sin(th),(th,0,s.pi)),(r,1,2)),
         s.Integral(r,(r,1,2))+s.Integral(s.Function('W_1')(r,s.cos(th),x[2],time)*s.sin(th),(th,0,s.pi),(r,1,2))),
        (s.Integral(s.Integral(s.Function('W_1')(r,s.cos(ph),x[2],time)*(1+s.cos(th)**2),(th,0,s.pi))+r,(ph,0,2*s.pi),(r,1,2)),
         s.Integral(s.Function('W_1')(r,s.cos(ph),x[2],time)*(1+s.cos(th)**2),(th,0,s.pi),(ph,0,2*s.pi),(r,1,2))+s.Integral(r,(ph,0,2*s.pi),(r,1,2)))):
        values=fixture_run({'object':nested},{'object':flat})
        equivalent_samples(values)
        detects(values,fixture_run({'object':nested+1},{'object':flat}),True,exact=True)


def split_integral_evaluation():
    # A rational split with polynomial sum. Test the same integral through
    # native serialization, extraction, joining and all evaluation strata.
    r=s.Symbol('s');field=s.Function('W_1')(r,x[1],x[2],x[3])
    first=s.Integral(field/(1+r**2),(r,0,1))
    second=s.Integral(r**2*field/(1+r**2),(r,0,1))
    whole=s.Integral(field,(r,0,1))
    options={'forms':c.FORMS,'strata':c.STRATA,'points':tuple(range(len(c.POINTS)))}
    base=fixture_run({'object':first+second},{'object':whole},options=options)
    equivalent_samples(base,'linear integral evaluation changed value')
    changed=fixture_run({'object':first+second+s.Integral(1,(r,0,1))},{'object':whole},options=options)
    detects(base,changed,True,exact=True)
    # Different bounds and a product stay distinct under the linear evaluator.
    # Their values are separately established by the native comparison path.
    for obj in (s.Integral(field,(r,0,2))+s.Integral(field,(r,0,1)),
                s.Integral(field,(r,0,2))*s.Integral(field,(r,0,1))):
        reference=c.substituted(obj,c.FORMS[0],c.STRATA[0],0)
        values=fixture_run({'object':obj},{'object':reference})
        equivalent_samples(values)
        detects(values,fixture_run({'object':obj+1},{'object':reference}),True,exact=True)


def sparse_integral_oracle():
    # Independent dense Fourier expansion on synthetic polynomials only.
    # The oracle integrates each exponential mode at its endpoints. It neither
    # calls angular_moment nor uses that function's beta/gamma moment formula.
    th,ph,r=s.symbols('theta phi s')
    sn,cs,sp,cp=s.sin(th),s.cos(th),s.sin(ph),s.cos(ph)
    @lru_cache(None)
    def fourier_moment(p,q,full):
        terms=[]
        for i in range(p+1):
            for j in range(q+1):
                frequency=p+q-2*(i+j)
                coefficient=s.binomial(p,i)*s.binomial(q,j)*(-1)**i/(2**(p+q)*s.I**p)
                endpoint=(2*s.pi if full else s.pi) if frequency==0 else (s.S.Zero if full else ((s.Integer(-1))**frequency-1)/(s.I*frequency))
                terms.append(coefficient*endpoint)
        return s.simplify(s.Add(*terms))
    def dense(integrand,radial=False):
        generators=(sn,cs,sp,cp,r) if radial else (sn,cs,sp,cp)
        polynomial=s.Poly(s.expand(integrand),*generators)
        terms=[]
        for powers,coefficient in polynomial.terms():
            value=coefficient*fourier_moment(*powers[:2],False)*fourier_moment(*powers[2:4],True)
            if radial:value*=(s.Rational(5,2)**(powers[4]+1)-s.Rational(2,3)**(powers[4]+1))/(powers[4]+1)
            terms.append(value)
        return s.Add(*terms)
    pressure=(1+(sn*cp+1)**2+(cs+2)**2)**4*sn
    volume=(1+(r*sn*cp+1)**2+(r*cs+2)**2)**4*sn*r**2
    denominator=1+sn+cs*sp
    numerator=s.expand(denominator*(2+sn+cp)**4)
    rational=s.Mul(numerator,s.Pow(denominator,-1,evaluate=False),evaluate=False)*sn
    for expression,polynomial,radial in ((pressure,pressure,False),(volume,volume,True),(rational,(2+sn+cp)**4*sn,False)):
        limits=((r,s.Rational(2,3),s.Rational(5,2)),(th,0,s.pi),(ph,0,2*s.pi)) if radial else ((th,0,s.pi),(ph,0,2*s.pi))
        integral=s.Integral(expression,*limits);oracle=dense(polynomial,radial)
        base=fixture_run({'object':integral},{'object':oracle})
        equivalent_samples(base,'sparse moments differ from independent Fourier expansion')
        # Alter the actual production moment, not the comparison's output.
        # Clear evaluator caches so a previously evaluated integral cannot hide
        # the defect; the already serialized Fourier oracle remains unchanged.
        original=c.angular_moment
        def clear():
            for fn in (c.evaluated_node,c.evaluate_cached,c.integrate_exact):fn.cache_clear()
        try:
            clear();c.angular_moment=lambda p,q,full:original(p,q,full)/(2 if full else 1)
            altered=replay(base)
        finally:c.angular_moment=original;clear()
        message='sparse moments differ from independent Fourier expansion'
        rejected_assertion(lambda:require(samples(altered) and samples(altered)==samples(altered,'wl'),message),message)
        assert samples(altered,'wl')==samples(base,'wl'),'moment defect changed the independent oracle'
        assert exacts(altered)==exacts(base),'moment defect changed symbolic subtraction'
        ABLATIONS.append({'control':ACTIVE_CONTROL,'instance':len(ABLATIONS)+1,'binding':'synthetic',
            'fault':'full_range_moment_halved','site':'angular_moment','removed':'full-range moment',
            'installed':'original moment / 2 for full=True only; evaluator caches cleared',
            'assertion':message,'fixtures':[[str(p) for p in base.recipe[0]]]})
        changed=fixture_run({'object':s.Integral(expression+sn,*limits)},{'object':oracle})
        detects(base,changed,True,exact=True)


def exclusive_guard_controls():
    # No systemd launch here: test the opt-in lock and limit predicates, then
    # the enclosing real guarded test run verifies their kernel/systemd values.
    path=Path('/var/projects/toy_physics/scripts/s11c_guarded_run.py')
    spec=importlib.util.spec_from_file_location('s15_guard_controls',path)
    guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
    cpus=sorted(os.sched_getaffinity(0))
    assert len(cpus)==8, 'eight-CPU guarded allocation absent'
    manifest={'pool':'s15-comparator','memoryMax':55*guard.GIB,'swapMax':0,'tasksMax':32,
              'poolMemoryMax':guard.POOL_MEMORY_MAX,'priorityPolicy':'desktop-managed',
              'affinityEnforcement':'systemd-and-inherited-process','cpu':cpus[-1],'cpus':cpus}
    actual={'nativeAddressSpaceBytes':[manifest['memoryMax']]*2,'systemd.CPUAffinity':' '.join(map(str,cpus)),
            'memory.max':str(manifest['memoryMax']),'memory.swap.max':'0','pids.max':'32',
            'affinity':cpus,'threads':{k:'1' for k in guard.THREADS}}
    assert guard.limits_match(actual,manifest)
    assert guard.limits_match({**actual,'systemd.CPUAffinity':'0-7'},manifest)
    for invalid in ('0-7 0','7-0','0-x','',None):
        assert guard.affinity_cpus(invalid) is None
    for key,value in (('affinity',cpus[:1]),('systemd.CPUAffinity',str(cpus[-1])),('memory.max',str(56*guard.GIB)),('memory.swap.max','1')):
        assert not guard.limits_match({**actual,key:value},manifest),'guard ignored allocation mutation'
    legacy={k:v for k,v in manifest.items() if k!='cpus'}
    assert guard.limits_match({**actual,'affinity':[legacy['cpu']],'systemd.CPUAffinity':str(legacy['cpu'])},legacy)
    with tempfile.TemporaryDirectory(prefix='guard_lock_',dir=ROOT) as directory:
        lockpath=Path(directory)/'admission.lock'
        held=guard.acquire_execution_lock(lockpath)
        try:
            for shared in (True,False):
                try:guard.acquire_execution_lock(lockpath,pooled=shared)
                except BlockingIOError:pass
                else:raise AssertionError('exclusive allocation admitted overlapping launcher')
        finally:held.close()
        guard.acquire_execution_lock(lockpath,pooled=True).close()


TESTS=(exclusive_guard_controls,nested_radial_derivative,nested_integral_parameters,split_integral_evaluation,sparse_integral_oracle,parallel_schedule,long_native_payload,basic_mutations,labels_and_conditions,equivalences,derivative_points,mixed_second_derivative_points,partial_time_semantics,integral_changes,binding_repoints,name_repoints,bound_name_repoints,pieces_groups_and_components,packaging_and_accounting,equation_relation_kind,repeated_records,supplied_equation_and_alias,coverage_controls,test_functions,duplicates,dag_sampling,quadrature)


def control_job(name):
    global ACTIVE_CONTROL,ABLATIONS,ROOT
    ACTIVE_CONTROL=name;ABLATIONS=[]
    ROOT=Path('/var/projects/toy_physics/_scratch/s15/cmp_work/fixtures')/name
    ROOT.mkdir(parents=True,exist_ok=True)
    print('WORKING_ON: control '+name,flush=True)
    try:
        globals()[name]()
        return name,ABLATIONS,None
    except Exception:return name,ABLATIONS,traceback.format_exc()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
    started=time.monotonic()
    global ACTIVE_CONTROL,ABLATIONS
    try:
        # Independent synthetic groups share no comparator state. The scheduling
        # control runs in the parent, separately, because it creates its own pool.
        jobs=[fn.__name__ for fn in TESTS if fn!=parallel_schedule]
        results={}
        with multiprocessing.get_context('fork').Pool(8) as pool:
            for name,ablations,error in pool.imap_unordered(control_job,jobs):
                results[name]=(ablations,error)
                print('WORKING_ON: completed '+name,flush=True)
        if not any(error for _,error in results.values()):
            name,ablations,error=control_job('parallel_schedule');results[name]=(ablations,error)
        ABLATIONS=[]
        errors=[]
        for fn in TESTS:
            if fn.__name__ not in results:continue
            ablations,error=results[fn.__name__]
            for entry in ablations:entry['instance']=len(ABLATIONS)+1;ABLATIONS.append(entry)
            if error:errors.append(fn.__name__+'\n'+error)
            else:COUNTS[fn.__name__]='completed'
        if errors:raise AssertionError('\n'.join(errors))
        for fn in TESTS:
            if fn in (exclusive_guard_controls,test_functions):continue
            assert any(q['control']==fn.__name__ for q in ABLATIONS),'control has no targeted ablation: '+fn.__name__
    except Exception:
        args.output.write_text(json.dumps({'completed':COUNTS,'targeted_ablations':ABLATIONS,'exception':traceback.format_exc(),'runtime_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2)+'\n')
        raise
    args.output.write_text(json.dumps({'controls':COUNTS,'targeted_ablations':ABLATIONS,'synthetic_only':True,'runtime_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2)+'\n')
    print(json.dumps({'control_groups':len(COUNTS),'targeted_ablations':len(ABLATIONS),'runtime_seconds':time.monotonic()-started}),flush=True)

if __name__=='__main__':main()
