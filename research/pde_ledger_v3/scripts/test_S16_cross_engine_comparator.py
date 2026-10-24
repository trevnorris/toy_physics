#!/usr/bin/env python3
"""Synthetic native fixtures only. Never opens a measured transcript."""
import argparse
from dataclasses import replace, asdict
import hashlib
from io import StringIO
import json
from pathlib import Path
import resource
import time
import traceback
import sympy as s
from sympy.core.function import AppliedUndef
import S16_cross_engine_comparator as c

ROOT=Path('/var/projects/toy_physics/_scratch/s16/cmp_work/fixtures')
ROOT.mkdir(parents=True,exist_ok=True)
FAST={'forms':c.FORMS[:1],'strata':('generic',),'points':(0,)}
COUNTS={}


def serial(v,e):
    if isinstance(v,str): return v
    if isinstance(v,dict):
        if set(v)=={'lhs','rhs'}:
            a,b=serial(v['lhs'],e),serial(v['rhs'],e)
            return f'Eq({a}, {b})' if e=='py' else f'{a} -> {b}'
        return ('{' if e=='py' else '<|')+', '.join((k+': '+serial(x,e)) if e=='py' else ('"'+k+'" -> '+serial(x,e)) for k,x in v.items())+('}' if e=='py' else '|>')
    if isinstance(v,(list,tuple)):
        return ('[' if e=='py' else '{')+', '.join(serial(x,e) for x in v)+(']' if e=='py' else '}')
    if e=='py': return s.sstr(v)
    if v==s.pi:return 'Pi'
    if v==s.true:return 'True'
    if v==s.false:return 'False'
    if isinstance(v,s.Symbol): return v.name
    if isinstance(v,s.Number): return str(v)
    if isinstance(v,AppliedUndef): return v.func.__name__+'['+', '.join(serial(x,e) for x in v.args)+']'
    if isinstance(v,s.Derivative) and isinstance(v.expr,AppliedUndef):
        f=v.expr; counts=dict(v.variable_count)
        return 'Derivative['+','.join(str(counts.get(x,0)) for x in f.args)+']['+f.func.__name__+']['+', '.join(serial(x,e) for x in f.args)+']'
    if isinstance(v,s.Add):return '('+' + '.join(serial(x,e) for x in v.args)+')'
    if isinstance(v,s.Mul):return '('+'*'.join(serial(x,e) for x in v.args)+')'
    if isinstance(v,s.Pow):return '('+serial(v.base,e)+')^('+serial(v.exp,e)+')'
    heads={s.sin:'Sin',s.cos:'Cos',s.exp:'Exp',s.sqrt:'Sqrt'}
    if v.func in heads:return heads[v.func]+'['+','.join(serial(x,e) for x in v.args)+']'
    raise ValueError('fixture serializer '+str(v))


def basic_row(kind='value',normalized=False):
    return c.Binding('synthetic','synthetic',c.End('object'),c.End('object'),'synthetic object','fixture','fixture',kind,normalized)

TEMPLATES={
'py':{'integral.Volume':'Integral(rho**2*q(Xc1(t)+rho*sin(theta)*cos(phi), Xc2(t)+rho*sin(theta)*sin(phi), Xc3(t)+rho*cos(theta), t)*sin(theta), (rho,r_in(t),R(t)), (theta,0,pi),(phi,0,2*pi))',
      'integral.Surface.0':'Integral(q(theta,phi,t)*r_in(t)**2*sin(theta),(theta,0,pi),(phi,0,2*pi))',
      'integral.Surface.1':'Integral(q(theta,phi,t)*R(t)**2*sin(theta),(theta,0,pi),(phi,0,2*pi))'},
'wl':{'integral.V':'Inactive[Integrate][rho^2*q[c1[t]+rho*Sin[theta]*Cos[phi],c2[t]+rho*Sin[theta]*Sin[phi],c3[t]+rho*Cos[theta],t]*Sin[theta],{rho,rin[t],rout[t]},{theta,0,Pi},{phi,0,2*Pi}]',
      'integral.B.1':'Inactive[Integrate][q[theta,phi,t]*rin[t]^2*Sin[theta],{theta,0,Pi},{phi,0,2*Pi}]',
      'integral.B.2':'Inactive[Integrate][q[theta,phi,t]*rout[t]^2*Sin[theta],{theta,0,Pi},{phi,0,2*Pi}]'}}
DOMAINS={'py':{'domain.unnormalized':'(r_in(t)>0) & (R(t)>r_in(t))','domain.normalized':'(r_in(t)>0) & (R(t)>r_in(t)) & Ne(N(t),0)'},
         'wl':{'domain.unnormalized':'0<rin[t]<rout[t] && 0<theta<Pi && Element[theta | phi | t, Reals]','domain.normalized':'0<rin[t]<rout[t] && 0<theta<Pi && Element[theta | phi | t, Reals] && count[t]!=0'}}


def fixture_run(py,wl,rows=None,names=c.NAME_TABLE,options=None,templates=True):
    rows=(basic_row(),) if rows is None else rows
    paths=[]
    for e,objects in (('py',py),('wl',wl)):
        data={**(TEMPLATES[e] if templates else {}),**DOMAINS[e],**objects}
        text='\n'.join(k+' = '+serial(v,e) for k,v in data.items())+'\n'
        path=ROOT/(hashlib.sha256((e+text).encode()).hexdigest()+'.'+e)
        if not path.exists():path.write_text(text)
        paths.append(path)
    output=StringIO()
    c.run(*paths,output,names=names,bindings=rows,options=FAST if options is None else options,progress=False,metadata=False)
    return [json.loads(line) for line in output.getvalue().splitlines()]


def computations(records):
    return [r for r in records if r['kind'] in ('comparison','binding_accounting','tag_accounting')]

def comparable(records):
    result=[]
    for r in records:
        if r['kind']=='comparison':
            result.append({k:v for k,v in r.items() if k not in ('py_value','wl_value','py_alias_resolved','wl_alias_resolved','kind')})
    return result


def samples(records,engine='py'):
    return [x.get(engine) for r in records if r['kind']=='comparison' for x in r.get('samples',[]) if x.get('coverage')=='evaluated']


def detects(base,mutated,evaluation=False):
    assert comparable(base)!=comparable(mutated), 'mutation invisible to computation'
    assert [r['wl'] for r in base if r['kind']=='operands']==[r['wl'] for r in mutated if r['kind']=='operands'], 'other operand changed'
    if samples(base,'wl') and samples(mutated,'wl'):
        assert samples(base,'wl')==samples(mutated,'wl'), 'other operand samples changed'
    if evaluation: assert samples(base)!=samples(mutated), 'mutation invisible to evaluation'


def test(name,fn):
    print('WORKING_ON = control '+name,flush=True)
    fn(); COUNTS[name]='completed'

x=s.symbols('x1 x2 x3 t'); F=s.Function('F1')(*x)


def mutation_control():
    base=fixture_run({'object':F},{'object':F})
    for v in (F+1,F+x[0]**3,-F,F.subs(x[0],x[0]+1),s.Integer(7),'F1(x1,x2,t)'):
        changed=fixture_run({'object':v},{'object':F})
        detects(base,changed,evaluation=not isinstance(v,str))
    disabled=fixture_run({'object':F},{'object':F},options={**FAST,'disable_comparison':True})
    disabled2=fixture_run({'object':F+1},{'object':F},options={**FAST,'disable_comparison':True})
    assert comparable(disabled)==comparable(disabled2)


def own_labels():
    row=basic_row('label')
    for py,wl in [('alpha','"alpha"'),('alpha','Missing["NoSecondRoute","alpha"]')]:
        # Native unquoted labels must be tagged as label payloads in sstr.
        r=replace(row,py=c.End('coverage.fluxDecomposition'),wl=c.End('coverage.fluxDecomposition'))
        base=fixture_run({r.py.tag:py},{r.wl.tag:wl},(r,))
        mut=fixture_run({r.py.tag:'beta'},{r.wl.tag:wl},(r,))
        detects(base,mut)
        other=fixture_run({r.py.tag:py},{r.wl.tag:wl.replace('alpha','beta')},(r,))
        assert comparable(base)!=comparable(other)


def equivalents():
    cases=[('Subs(Derivative(F1(u,v,w,z),u),(u,v,w,z),(x1+1,x2,x3,t))',
            'Derivative[1,0,0,0][F1][x1+1,x2,x3,t]'),
           ('Derivative(F1(x1,x2,x3,t),x2)','Derivative[0,1,0,0][F1][x1,x2,x3,t]'),
           ('Surface(0, sin(theta)**2*cos(phi)**2 + sin(theta)**2*sin(phi)**2 + cos(theta)**2)','B[1,1]'),
           ('Volume(F1(x1,x2,x3,t))','V[F1[x1,x2,x3,t]]')]
    for py,wl in cases:
        records=fixture_run({'object':py},{'object':wl})
        for r in records:
            if r['kind']=='comparison':
                for q in r.get('samples',[]):
                    assert q.get('coverage')=='evaluated',q
                    assert q['py']==q['wl'],'native notation samples differ'
    same=fixture_run({'object':'Volume(F1(x1,x2,x3,t))'},{'object':'V[F1[x1,x2,x3,t]]'})
    expanded=fixture_run({'object':'Volume(F1(x1,x2,x3,t)+x1-x1)'},{'object':'V[F1[x1,x2,x3,t]]'})
    assert comparable(same)==comparable(expanded)
    compositional=fixture_run(
        {'object':'Derivative(Xc1(t),t)*(Volume(F1(x1,x2,x3,t))+Surface(0,F2(x1,x2,x3,t)))**2'},
        {'object':'Derivative[1][c1][t]*(V[F1[x1,x2,x3,t]]^2+2*V[F1[x1,x2,x3,t]]*B[1,F2[x1,x2,x3,t]]+B[1,F2[x1,x2,x3,t]]^2)'},options={})
    assert samples(compositional)==samples(compositional,'wl')
    # Compare same objects under relations, domains, labels, containers, dimensions.
    r=replace(basic_row(),py=c.End('object',('rhs',)),wl=c.End('object',('rhs',)))
    records=fixture_run({'object':'Eq(N(t),F1(x1,x2,x3,t))'},{'object':'count[t] -> F1[x1,x2,x3,t]'},(r,))
    assert samples(records)==samples(records,'wl')
    for kind,py,wl in [('condition','(r_in(t)>0) & (R(t)>r_in(t))','0<rin[t]<rout[t]'),
                       ('dimension','[unitCount/unitTime]','{{0,-1,1,0}}')]:
        records=fixture_run({'object':py},{'object':wl},(basic_row(kind),))
        comparisons=[r for r in records if r['kind']=='comparison']
        assert comparisons and all('not_formed' not in r for r in comparisons)
        if kind=='condition':
            assert samples(records) and samples(records)==samples(records,'wl')
        else:
            for r in comparisons:
                assert r['dimension_evaluation']['py']==r['dimension_evaluation']['wl']


def derivative_points():
    for derivative_slot in range(4):
        var=('u','v','w','z')[derivative_slot]
        orders=[0]*4;orders[derivative_slot]=1
        for point_slot in range(4):
            for context in ('subs','composite','volume','surface'):
                pts=['x1','x2','x3','t'];wpts=pts.copy()
                if context=='composite': pts[0]='x1+2/9';wpts[0]='x1+2/9'
                if context=='surface':
                    pts=['Xc1(t)+r_in(t)*sin(theta)*cos(phi)','Xc2(t)+r_in(t)*sin(theta)*sin(phi)','Xc3(t)+r_in(t)*cos(theta)','t']
                    wpts=['c1[t]+rin[t]*Sin[theta]*Cos[phi]','c2[t]+rin[t]*Sin[theta]*Sin[phi]','c3[t]+rin[t]*Cos[theta]','t']
                shifted=pts.copy();shifted[point_slot]+='+1/7'
                py='Subs(Derivative(F1(u,v,w,z),'+var+'),(u,v,w,z),('+','.join(pts)+'))'
                changed='Subs(Derivative(F1(u,v,w,z),'+var+'),(u,v,w,z),('+','.join(shifted)+'))'
                wl='Derivative['+','.join(map(str,orders))+'][F1]['+','.join(wpts)+']'
                if context in ('volume','surface'):
                    ph='Volume(' if context=='volume' else 'Surface(0,'
                    wh='V[' if context=='volume' else 'B[1,'
                    py=ph+py+')';changed=ph+changed+')';wl=wh+wl+']'
                base=fixture_run({'object':py},{'object':wl})
                mut=fixture_run({'object':changed},{'object':wl})
                detects(base,mut,True)
    py='Subs(Derivative(Xc1(u),(u,2)),u,t)';wl='Derivative[2][c1][t]'
    detects(fixture_run({'object':py},{'object':wl}),fixture_run({'object':py.replace(',u,t)',',u,t+1/7)')},{'object':wl}),True)
    base=fixture_run({'object':'Derivative(F1(x1,x2,x3,t),x1)'},{'object':'Derivative[1,0,0,0][F1][x1,x2,x3,t]'})
    for bad in ('Derivative(F1(x1,x2,x3,t),x2)','Derivative(F1(x1,x2,x3,t),(x1,2))'):
        detects(base,fixture_run({'object':bad},{'object':'Derivative[1,0,0,0][F1][x1,x2,x3,t]'}),True)


def integral_and_domains():
    for p,w,bad in [('Surface(0,F1(x1,x2,x3,t))','B[1,F1[x1,x2,x3,t]]','Surface(1,F1(x1,x2,x3,t))'),
                    ('Volume(F1(x1,x2,x3,t))','V[F1[x1,x2,x3,t]]','Volume(F1(x1,x2,x3,t)+x1**2)'),
                    ('Surface(0,F1(x1,x2,x3,t))','B[1,F1[x1,x2,x3,t]]','Surface(0,F1(x1,x2,x3,t)+x1**2)')]:
        base=fixture_run({'object':p},{'object':w}); detects(base,fixture_run({'object':bad},{'object':w}),True)
    base=fixture_run({'object':'Volume(F1(x1,x2,x3,t))'},{'object':'V[F1[x1,x2,x3,t]]'})
    py={'object':'Volume(F1(x1,x2,x3,t))','integral.Volume':TEMPLATES['py']['integral.Volume'].replace('(rho,r_in(t),R(t))','(rho,r_in(t),R(t)+1)')}
    detects(base,fixture_run(py,{'object':'V[F1[x1,x2,x3,t]]'}),True)
    base=fixture_run({'object':F},{'object':F})
    mut=fixture_run({'object':F,'domain.unnormalized':'(r_in(t)>0) & (R(t)>r_in(t)) & Ne(N(t),0)'},{'object':F},options={**FAST,'strata':('generic','count_zero')})
    full=fixture_run({'object':F},{'object':F},options={**FAST,'strata':('generic','count_zero')})
    assert comparable(full)!=comparable(mut)
    real_domain=fixture_run({'object':'x1>0'},{'object':'x1>0 && Element[theta | phi | t,Reals]'},(basic_row('condition'),))
    assert samples(real_domain)
    a=fixture_run({'object':'x1>0'},{'object':'x1>0'},(basic_row('condition'),))
    b=fixture_run({'object':'x1>1'},{'object':'x1>0'},(basic_row('condition'),))
    detects(a,b,True)


def volume_partial_time():
    # The independent operand spells out evaluation AFTER differentiation,
    # bypassing the volume head whose Cartesian binding is under test.
    point='Xc1(t)+rho*sin(theta)*cos(phi),Xc2(t)+rho*sin(theta)*sin(phi),Xc3(t)+rho*cos(theta),t'
    explicit='Integral(rho**2*sin(theta)*Subs(Derivative(F1(u,v,w,z),z),(u,v,w,z),('+point+')),(rho,r_in(t),R(t)),(theta,0,pi),(phi,0,2*pi))'
    wexplicit='Inactive[Integrate][rho^2*Sin[theta]*Derivative[0,0,0,1][F1][c1[t]+rho*Sin[theta]*Cos[phi],c2[t]+rho*Sin[theta]*Sin[phi],c3[t]+rho*Cos[theta],t],{rho,rin[t],rout[t]},{theta,0,Pi},{phi,0,2*Pi}]'
    options={**FAST,'strata':c.STRATA}
    for py,wl in [('Volume(Derivative(F1(x1,x2,x3,t),t))',wexplicit),
                  (explicit,'V[Derivative[0,0,0,1][F1][x1,x2,x3,t]]')]:
        records=fixture_run({'object':py},{'object':wl},options=options)
        assert samples(records) and samples(records)==samples(records,'wl'), 'partial time became chart time'
    partial=fixture_run({'object':'Volume(Derivative(F1(x1,x2,x3,t),t))'},
                        {'object':'V[Derivative[0,0,0,1][F1][x1,x2,x3,t]]'},options=options)
    convective='+'.join('Derivative(Xc%d(t),t)*Derivative(F1(x1,x2,x3,t),x%d)'%(i,i) for i in (1,2,3))
    along=fixture_run({'object':'Volume(Derivative(F1(x1,x2,x3,t),t)+'+convective+')'},
                      {'object':'V[Derivative[0,0,0,1][F1][x1,x2,x3,t]]'},options=options)
    detects(partial,along,True)
    assert [r.get('resolved_exact_residual') for r in partial if r['kind']=='comparison'] != [r.get('resolved_exact_residual') for r in along if r['kind']=='comparison']


def count_zero_fields():
    options={'forms':c.FORMS,'strata':('count_zero',),'points':range(len(c.POINTS))}
    for py,wl in [('Volume(Sigma(x1,x2,x3,t))','0'),
                  ('0','V[Sigma[x1,x2,x3,t]]'),('N(t)','0')]:
        records=fixture_run({'object':py},{'object':wl},options=options)
        assert samples(records) and samples(records)==samples(records,'wl'), 'printed count does not vanish'
    live=fixture_run({'object':'Derivative(Sigma(x1,x2,x3,t),x1)+Derivative(Sigma(x1,x2,x3,t),t)'},
                     {'object':'0'},options=options)
    assert samples(live)!=samples(live,'wl'), 'density was frozen'
    # A witness in the printed count, rather than just the separate N symbol.
    generic=fixture_run({'object':'Volume(Sigma(x1,x2,x3,t))'},{'object':'0'})
    witness='(Volume(Sigma(x1,x2,x3,t))-('+samples(generic)[0]+'))*F1(x1,x2,x3,t)'
    full=fixture_run({'object':witness},{'object':'0'},options={**FAST,'strata':('generic','count_zero')})
    reduced=fixture_run({'object':witness},{'object':'0'},options={**FAST,'strata':('generic','count_zero'),'normalize_all':True})
    assert samples(full)!=samples(full,'wl') and samples(reduced)==samples(reduced,'wl')
    normalized=fixture_run({'object':'Volume(Sigma(x1,x2,x3,t))/N(t)'},{'object':'V[Sigma[x1,x2,x3,t]]/count[t]'},
                          (basic_row(normalized=True),),options=options)
    coverage=[v['coverage'] for r in normalized if r['kind']=='comparison' for v in r.get('samples',[])]
    assert coverage and set(coverage)=={'outside normalized domain'}


def wrapper_resolution():
    identity=fixture_run({'object':F},{'object':F})
    def resolved(records):
        return [r['resolved_exact_residual'] for r in records if r['kind']=='comparison']
    for context in ('plain','volume','surface'):
        py='Subs(u*Derivative(F1(u,v,w,z),z),(u,v,w,z),(x1,x2,x3,t))'
        wl='x1*Derivative[0,0,0,1][F1][x1,x2,x3,t]'
        wrap=lambda p,w:(p,w)
        if context=='volume': wrap=lambda p,w:('Volume('+p+')','V['+w+']')
        if context=='surface': wrap=lambda p,w:('Surface(0,'+p+')','B[1,'+w+']')
        a,b=wrap(py,wl)
        base=fixture_run({'object':a},{'object':b})
        assert resolved(base)==resolved(identity), 'feasible wrapper remains unresolved'
        assert all('exact_residual' in r for r in base if r['kind']=='comparison')
        for changed in (py.replace('),z)', '),u)'),py.replace('),z)', '),(z,2))'),
                        py.replace('(x1,x2,x3,t))','(x1,x2,x3,t+1/7))')):
            a,_=wrap(changed,wl)
            mutated=fixture_run({'object':a},{'object':b})
            detects(base,mutated,True)
            assert resolved(mutated)!=resolved(base), 'wrapper resolution lost derivative identity'
    # Composite points retain partial time, even when the point moves in time.
    base=fixture_run({'object':'Subs(u*Derivative(F1(u,v,w,z),z),(u,v,w,z),(Xc1(t),x2,x3,t))'},
                     {'object':'c1[t]*Derivative[0,0,0,1][F1][c1[t],x2,x3,t]'})
    assert resolved(base)==resolved(identity)


def repair_ablations():
    def clears():
        c.substituted.cache_clear(); c.evaluate_cached.cache_clear()
    def rejected(control):
        try: control()
        except AssertionError: return
        raise AssertionError('repair control survived its targeted ablation')
    original=c.Names.integral
    def old_volume(self,engine,head,args):
        if head not in ('Volume','V'): return original(self,engine,head,args)
        template=self.templates[engine,'volume']
        q=next(f for f in template.atoms(AppliedUndef) if f.func.__name__==self.name('q',engine))
        coords=tuple(s.Symbol(self.name(n,engine)) for n in ('x1','x2','x3','t'))
        return template.xreplace({q:args[-1].subs(dict(zip(coords,q.args)),simultaneous=True)})
    try:
        c.Names.integral=old_volume; clears(); rejected(volume_partial_time)
    finally: c.Names.integral=original; clears()
    original_profile=c.profile
    try:
        c.profile=lambda name,arity,form,stratum,point,options=None:original_profile(name,arity,form,'generic' if name=='Sigma' and stratum=='count_zero' else stratum,point,options)
        clears(); rejected(count_zero_fields)
    finally: c.profile=original_profile; clears()
    original_resolution=c.resolved_exact_difference
    try:
        c.resolved_exact_difference=lambda residual:residual
        rejected(wrapper_resolution)
    finally: c.resolved_exact_difference=original_resolution


def schema_fixture(row,e,value):
    end=getattr(row,e); shape=c.ROOT_SHAPES.get((e,end.tag))
    def build(shape):
        if shape=='v': return s.Function('F1')(*x)
        if isinstance(shape,dict):return {k:build(v) for k,v in shape.items()}
        return [build(v) for v in shape]
    if shape is not None: result=build(shape)
    else:
        result=value
        for k in reversed(end.path): result={k:result} if isinstance(k,str) else [s.Function('F1')(*x)]*k+[result]
    if end.path:
        obj=result
        for k in end.path[:-1]:obj=obj[k]
        obj[end.path[-1]]=value
    else:result=value if shape is None or shape=='v' else result
    return result


def binding_repoints():
    # Every explicit rule is exercised with native source-shaped endpoints.
    # Kinds that need special literals are tested below in their own paths.
    representatives={r.rule:r for r in c.BINDING_TABLE}
    for rule,row in representatives.items():
        if row.kind in ('dimension','dimension_field','basis','label','condition'):
            val=s.Integer(3);other=s.Integer(5)
            if row.kind=='dimension': py,wl='[unitLength]','{{1,0,0,0}}';bad='[unitTime]'
            elif row.kind=='dimension_field':py,wl='unitLength','{1,0,0,0}';bad='unitTime'
            elif row.kind=='basis':py,wl='[unitLength]','{"length"}';bad='[unitTime]'
            elif row.kind=='condition':py,wl='x1>0','x1>0';bad='x1>1'
            else:py,wl='"a"','"a"';bad='"b"'
            rr=replace(row,py=c.End('object'),wl=c.End('object'))
            a=fixture_run({'object':py,'alternate':bad},{'object':wl},(rr,))
            b=fixture_run({'object':py,'alternate':bad},{'object':wl},(replace(rr,py=c.End('alternate')),))
        else:
            py=schema_fixture(row,'py',F);wl=schema_fixture(row,'wl',F)
            alt=schema_fixture(row,'py',F+x[0]**2)
            rr=replace(row,py=c.End(row.py.tag,row.py.path))
            a=fixture_run({row.py.tag:py,'alternate':alt},{row.wl.tag:wl},(rr,))
            b=fixture_run({row.py.tag:py,'alternate':alt},{row.wl.tag:wl},(replace(rr,py=c.End('alternate',row.py.path)),))
        assert comparable(a)!=comparable(b),rule
        if row.kind=='value':assert samples(a)!=samples(b),rule


def name_repoints():
    for row in c.NAME_TABLE:
        candidates=[r for r in c.NAME_TABLE if r!=row and r.arity==row.arity and not r.alias]
        if not candidates:candidates=[r for r in c.NAME_TABLE if r.arity==4 and not r.alias]
        partner=next((r for r in candidates if r.py=={'Pi12':'Pi21','Pi21':'Pi12'}.get(row.py)),candidates[0])
        arity=4 if row.arity==-1 else row.arity
        args=x if arity==4 else (x[3],)
        py=s.Symbol(row.py) if arity==0 else s.Function(row.py)(*args)
        wl=row.wl if arity==0 else row.wl+'['+','.join(map(str,args))+']'
        if arity==0 and row.py not in ('x1','x2','x3','t','rho','theta','phi') and not row.py.startswith('unit'):continue
        # Swap canonical names on the Python side, preserving injectivity.
        names=tuple(replace(r,py=partner.py,arity=partner.arity) if r==row else replace(r,py=row.py,arity=row.arity) if r==partner else r for r in c.NAME_TABLE)
        pobjects={'object':py,'domain.unnormalized':'True'};wobjects={'object':wl,'domain.unnormalized':'True'}
        if row.alias:
            pf,pp=row.alias.split('.term.')
            wtag=dict(c.FAMILIES)[pf]+'.term.'+dict(c.PARTS)[pp]
            pobjects[row.alias]=f'Eq({row.py}(t),F1(x1,x2,x3,t))'
            wobjects[wtag]=row.wl+'[t] -> F1[x1,x2,x3,t]'
        a=fixture_run(pobjects,wobjects,templates=False)
        b=fixture_run(pobjects,wobjects,names=names,templates=False)
        assert comparable(a)!=comparable(b),row.py
        assert samples(a)!=samples(b) or samples(a,'wl')!=samples(b,'wl'),row.py


def piece_and_packaging():
    rows=tuple(replace(basic_row(),key='piece'+str(i),py=c.End('split'+str(i)),wl=c.End('group',(i,))) for i in range(2))
    base=fixture_run({'split0':F,'split1':2*F},{'group':[F,2*F]},rows)
    for i in range(2):
        mutated={'split0':F,'split1':2*F}; mutated['split'+str(i)]+=x[0]**2
        detects(base,fixture_run(mutated,{'group':[F,2*F]},rows),True)
    detects(base,fixture_run({'split0':F+1,'split1':2*F-1},{'group':[F,2*F]},rows),True)
    # Repoint each component-selector mechanism (sequence, dictionary, equation,
    # matrix, nested sequence) through the same extraction path.
    fixtures=[([F,2*F],(0,),(1,)),({'a':F,'b':2*F},('a',),('b',)),
              ({'lhs':F,'rhs':2*F},('lhs',),('rhs',)),([[F],[2*F]],(0,0),(1,0))]
    for payload,path,alt in fixtures:
        r=replace(basic_row(),py=c.End('object',path),wl=c.End('object',path))
        a=fixture_run({'object':payload},{'object':payload},(r,))
        b=fixture_run({'object':payload},{'object':payload},(replace(r,py=c.End('object',alt)),))
        detects(a,b,True)
    for i in range(2):
        payload=[F,2*F];payload[i]+=x[0]**2
        changed=fixture_run({'split0':F,'split1':2*F},{'group':payload},rows)
        assert comparable(base)!=comparable(changed)
        assert samples(base)==samples(changed)
        assert samples(base,'wl')!=samples(changed,'wl')
    for kind in ('dimension','dimension_field'):
        row=replace(basic_row(kind),py=c.End('object',(0,)),wl=c.End('object',(0,)))
        pp=['[unitLength]','[unitTime]'] if kind=='dimension' else ['unitLength','unitTime']
        ww=['{{1,0,0,0}}','{{0,1,0,0}}'] if kind=='dimension' else ['{1,0,0,0}','{0,1,0,0}']
        wtag='object'
        if kind=='dimension':
            # Use the source-defined three-component native dimension layout.
            wtag='supplied.momentum.dimensions'
            row=replace(row,wl=c.End(wtag,(0,)))
            pp.append(pp[0]);ww.append(ww[0])
        a=fixture_run({'object':pp},{wtag:ww},(row,))
        b=fixture_run({'object':pp},{wtag:ww},(replace(row,py=c.End('object',(1,))),))
        assert comparable(a)!=comparable(b)
    tag='supplied.momentum';row=next(r for r in c.BINDING_TABLE if r.py.tag==tag and r.py.path==(0,))
    a=fixture_run({tag:[F,2*F,3*F]},{tag:[F,2*F,3*F]},(row,))
    b=fixture_run({tag:[F,2*F]},{tag:[F,2*F,3*F]},(row,))
    assert comparable(a)!=comparable(b)
    assert any('not_formed' in q for q in b)
    nested='geometry.surfaces'
    nr=next(r for r in c.BINDING_TABLE if r.py.tag==nested and r.py.path==(0,0))
    intact=[[F,2*F,3*F],[4*F,5*F,6*F]]
    removed=[[F,2*F,3*F],[4*F,5*F]]
    a=fixture_run({nested:intact},{nr.wl.tag:intact},(nr,))
    b=fixture_run({nested:removed},{nr.wl.tag:intact},(nr,))
    assert comparable(a)!=comparable(b) and any('not_formed' in r for r in b)
    dr=next(r for r in c.BINDING_TABLE if r.kind=='dimension' and r.wl.tag=='surfaces.chart.dimensions' and r.wl.path==(0,0))
    dim='{{1,0,0,0}}'
    a=fixture_run({dr.py.tag:'[unitLength]'},{dr.wl.tag:[[dim,dim,dim],[dim,dim,dim]]},(dr,))
    b=fixture_run({dr.py.tag:'[unitLength]'},{dr.wl.tag:[[dim,dim,dim],[dim,dim]]},(dr,))
    assert comparable(a)!=comparable(b) and any('not_formed' in r for r in b)
    mismatch=fixture_run({'object':[F,2*F]},{'object':F})
    assert any('not_formed' in r for r in mismatch)
    lr=replace(basic_row(),py=c.End('object',(0,)),wl=c.End('object',(0,)))
    a=fixture_run({'object':[F]},{'object':[F]},(lr,))
    b=fixture_run({'object':[F,2*F]},{'object':[F]},(lr,))
    assert comparable(a)==comparable(b)
    la=next(r for r in a if r['kind']=='tag_accounting' and r['tag']=='object' and r['engine']=='py')
    lb=next(r for r in b if r['kind']=='tag_accounting' and r['tag']=='object' and r['engine']=='py')
    assert la['parsed_leaves']!=lb['parsed_leaves'] and la['actually_compared_leaves']==lb['actually_compared_leaves']
    assert lb['unpaired_components']
    # A multiline Matrix is not split into separate tags.
    r=replace(basic_row(),py=c.End('object',(0,)),wl=c.End('object',(0,)))
    a=fixture_run({'object':'Matrix([\n[F1(x1,x2,x3,t)],\n[2*F1(x1,x2,x3,t)]])'},{'object':[F,2*F]},(r,))
    b=fixture_run({'object':'Matrix([\n[F1(x1,x2,x3,t)+1],\n[2*F1(x1,x2,x3,t)]])'},{'object':[F,2*F]},(r,))
    detects(a,b,True)
    c0=fixture_run({'object':[F,2*F]},{'object':[F,2*F]},(r,),options={**FAST,'omit_leaf_credit':True})
    assert any(q['kind']=='binding_accounting' and q['parsed_leaves']!=q['actually_compared_leaves'] for q in c0)


def domain_and_balance_factor():
    for a,b in [('sqrt(r_in(t)**2)/r_in(t)*F1(x1,x2,x3,t)','F1[x1,x2,x3,t]'),
                ('Derivative(Sigma(x1,x2,x3,t),t)','Derivative[0,0,0,1][Sigma][x1,x2,x3,t]')]:
        baseline=fixture_run({'object':a},{'object':b})
        new='F1(x1,x2,x3,t)' if a.startswith('sqrt') else 'a_minus(x1,x2,x3,t)*n1(x1,x2,x3,t)-a_plus(x1,x2,x3,t)*n3(x1,x2,x3,t)-Derivative(F1(x1,x2,x3,t),x1)-Derivative(F2(x1,x2,x3,t),x2)-Derivative(F3(x1,x2,x3,t),x3)'
        altered=fixture_run({'object':new},{'object':b})
        assert [r.get('exact_residual') for r in baseline if r['kind']=='comparison'] != [r.get('exact_residual') for r in altered if r['kind']=='comparison']


def coverage_controls():
    # Each witness has no sample difference under the named reduction, and a
    # sample difference with the full coverage. No exact-residual-only assertion.
    def witness(py,wl,reduction,full=None):
        options={**FAST,**(full or {})}
        records=fixture_run({'object':py},{'object':wl},options=options)
        reduced=fixture_run({'object':py},{'object':wl},options={**options,**reduction})
        assert samples(records)!=samples(records,'wl'), 'full coverage missed witness'
        assert samples(reduced)==samples(reduced,'wl'), 'reduction did not remove witness'
    witness('(x1-1/5)*F1(x1,x2,x3,t)','0',{'points':(0,)},{'points':(0,1,2)})
    witness('Derivative(F1(x1,x2,x3,t),(x1,3))','0',{'forms':c.FORMS[:1]},{'forms':c.FORMS})
    # Under symmetric slots, first derivatives coincide at equal argument points.
    witness('Subs(Derivative(F1(u,v,w,z),u)-Derivative(F1(u,v,w,z),v),(u,v,w,z),(1/3,1/3,1/3,1/3))','0',{'symmetric_arguments':True})
    witness('Derivative(F1(x1,x2,x3,t),x1,x2)','0',{'uncoupled_arguments':True})
    speed=s.diff(c.profile('Xc1',1,c.FORMS[0],'generic',0),c.FORMAL[0]).subs(c.FORMAL[0],c.POINTS[0][3])
    witness(f'(Derivative(Xc1(t),t)-({speed}))*F1(x1,x2,x3,t)','0',{'strata':('generic',)},{'strata':('generic','center_velocity_zero')})
    witness(f'(N(t)-({s.Rational(c.CODES["N"],13)}))*F1(x1,x2,x3,t)','0',{'normalize_all':True},{'strata':('generic','count_zero')})
    # Also check that the zero-count coverage outcome itself changes.
    full=fixture_run({'object':F},{'object':F},options={**FAST,'strata':('count_zero',)})
    reduced=fixture_run({'object':F},{'object':F},options={**FAST,'strata':('count_zero',),'normalize_all':True})
    assert samples(full) and not samples(reduced)
    for expression in ('F1(x1,x2,x3,t)-2*F2(x1,x2,x3,t)+F3(x1,x2,x3,t)', 'Xc1(t)-2*Xc2(t)+Xc3(t)', 'Derivative(Xc1(t),(t,2))-2*Derivative(Xc2(t),(t,2))+Derivative(Xc3(t),(t,2))'):
        independent=fixture_run({'object':expression},{'object':'0'},options={**FAST,'forms':c.FORMS,'strata':('generic','all_velocities_zero')})
        assert samples(independent)!=samples(independent,'wl'), 'test profiles impose an affine component relation'
    duplicate=fixture_run({'object':'F1(x1,x2,x3,t)-F2(x1,x2,x3,t)'},{'object':'0'},options={**FAST,'duplicate_functions':True})
    distinct=fixture_run({'object':'F1(x1,x2,x3,t)-F2(x1,x2,x3,t)'},{'object':'0'})
    assert samples(duplicate)==samples(duplicate,'wl') and samples(distinct)!=samples(distinct,'wl')


def alias_mutation():
    row=c.NAME_TABLE[-1];pf,pp=row.alias.split('.term.')
    wl_tag=dict(c.FAMILIES)[pf]+'.term.'+dict(c.PARTS)[pp]
    py={'object':row.py+'(t)',row.alias:f'Eq({row.py}(t),F1(x1,x2,x3,t))'}
    wl={'object':row.wl+'[t]',wl_tag:row.wl+'[t] -> F1[x1,x2,x3,t]'}
    base=fixture_run(py,wl)
    changed={**py,row.alias:f'Eq({row.py}(t),F2(x1,x2,x3,t))'}
    detects(base,fixture_run(changed,wl),True)
    py[row.alias.replace('.term.','.termRate.')]=f'Eq(Derivative({row.py}(t),t),F2(x1,x2,x3,t))'
    wl[wl_tag.replace('.term.','.termRate.')]='Derivative[1]['+row.wl+'][t] -> F2[x1,x2,x3,t]'
    py['object']=f'Derivative({row.py}(t),t)';wl['object']='Derivative[1]['+row.wl+'][t]'
    a=fixture_run(py,wl)
    py['object']=f'Subs(Derivative({row.py}(u),u),u,t+1/7)'
    detects(a,fixture_run(py,wl),True)


def duplicates():
    def rejected(fn):
        try:fn()
        except c.InputError:return
        raise AssertionError('duplicate not rejected')
    r=basic_row()
    rejected(lambda:c.validate_tables(c.NAME_TABLE,(r,replace(r,key='other'))))
    rejected(lambda:c.validate_tables(c.NAME_TABLE+(c.NAME_TABLE[0],),(r,)))
    rejected(lambda:c.validate_tables(c.NAME_TABLE,(r,),('rule','rule')))
    dl=replace(r,kind='dimension',wl=c.End('dim',('lhs',)))
    di=replace(r,key='other',kind='dimension',py=c.End('other'),wl=c.End('dim',(0,)))
    rejected(lambda:c.validate_tables(c.NAME_TABLE,(dl,di)))
    a=fixture_run({'object':F},{'object':F})
    # Change the serialized fixture, not the already-parsed object.
    path=ROOT/'duplicate_tag.py';path.write_text('object = 1\nobject = 2\n')
    other=ROOT/'duplicate_tag.wl';other.write_text('object = 1\n')
    out=StringIO();c.run(path,other,out,bindings=(r,),options=FAST,progress=False,metadata=False)
    assert any('missing or repeated' in line for line in out.getvalue().splitlines())


def quadrature_path():
    original=c.s.integrate
    try:
        c.s.integrate=lambda expr,lim:s.Integral(expr,lim)
        records=fixture_run({'object':'Integral(cos(rho),(rho,0,1))'},
                            {'object':'Inactive[Integrate][Cos[rho],{rho,0,1}]'})
        changed=fixture_run({'object':'Integral(cos(rho)+rho**2,(rho,0,1))'},
                            {'object':'Inactive[Integrate][Cos[rho],{rho,0,1}]'})
        detects(records,changed,True)
    finally: c.s.integrate=original
    assert samples(records)==samples(records,'wl')
    assert c.NUMERICAL_INTEGRALS
    assert all('error_estimate' in r and 'working_precisions' in r for r in c.NUMERICAL_INTEGRALS.values())


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    start=time.monotonic()
    tests=[mutation_control,own_labels,equivalents,derivative_points,integral_and_domains,volume_partial_time,count_zero_fields,wrapper_resolution,repair_ablations,binding_repoints,name_repoints,piece_and_packaging,domain_and_balance_factor,coverage_controls,alias_mutation,duplicates,quadrature_path]
    try:
        for fn in tests:test(fn.__name__,fn)
    except Exception:
        args.output.write_text(json.dumps({'completed':COUNTS,'exception':traceback.format_exc()},indent=2)+'\n')
        raise
    args.output.write_text(json.dumps({'controls':COUNTS,'synthetic_only':True,'runtime_seconds':time.monotonic()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2)+'\n')
    print(json.dumps({'control_groups':len(COUNTS),'runtime_seconds':time.monotonic()-start}),flush=True)

if __name__=='__main__':main()
