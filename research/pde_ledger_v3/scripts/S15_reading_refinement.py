"""Target the initial comparison's precision flags and compound exponent branches."""
import argparse
import ast
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
from functools import lru_cache
import gzip
import itertools
import json
import multiprocessing
import os
import sys
from pathlib import Path
import mpmath as mp
import sympy as s

import S15_reading_math as m
import S15_reading_comparator as c


@lru_cache(maxsize=65536)
def jet(name,arity,orders,form,stratum,point):
    expr, origin = m.sample_jet(name, arity, orders, form, stratum, point)
    fn = s.lambdify(m.original.FORMAL[:arity], expr, modules='mpmath', cse=True, docstring_limit=0)
    if origin is None:
        return fn
    return lambda t: fn(t - rational(origin))


def rational(v):
    v=s.Rational(v)
    return mp.mpf(int(v.p))/int(v.q)


class PrecisionSample:
    """Scalar evaluation of closed forms on the actual native piecewise branches."""
    def __init__(self,form,stratum,point,exponents):
        self.form,self.stratum,self.point=form,stratum,point
        self.env=dict(zip(s.symbols('x_1 x_2 x_3 t'),map(rational,m.POINTS[point])))
        self.env.update(zip(m.EXP,map(rational,exponents)))
        t=self.env[s.Symbol('t')]
        vec=[jet(f'X_1_{j}',1,(0,),form,stratum,point)(t)-jet(f'X_2_{j}',1,(0,),form,stratum,point)(t) for j in (1,2,3)]
        d=mp.sqrt(mp.fsum(q*q for q in vec))
        self.env[s.Symbol('d')]=d
        self.env.update(zip(m.DIRECTION,(q/d for q in vec)))
        self.cache={}

    def value(self,v):
        if v not in self.cache:self.cache[v]=self.compute(v)
        return self.cache[v]

    def compute(self,v):
        ev=self.value
        if v is s.true or v is s.false:return v is s.true
        if v.is_Rational:return rational(v)
        if v==s.pi:return mp.pi
        if v.is_Symbol:return self.env[v]
        if v.is_Add:return mp.fsum(ev(q) for q in v.args)
        if v.is_Mul:return mp.fprod(ev(q) for q in v.args)
        if v.is_Pow:return mp.power(ev(v.base),ev(v.exp))
        if v.func in (s.sin,s.cos,s.log):return {s.sin:mp.sin,s.cos:mp.cos,s.log:mp.log}[v.func](ev(v.args[0]))
        if v.func==s.polylog:return mp.polylog(*(ev(q) for q in v.args))
        if v.func==s.appellf1:
            a,b,c,d,x,y=(ev(q) for q in v.args)
            if not d>a>0 or max(x,y)>=1:raise m.Unsupported('outside real Euler domain')
            return mp.quad(lambda u:u**(a-1)*(1-u)**(d-a-1)*(1-x*u)**(-b)*(1-y*u)**(-c),[0,1])/mp.beta(a,d-a)
        if isinstance(v,s.Piecewise):
            for val,cond in v.args:
                if ev(cond):return ev(val)
            raise m.Unsupported('missing branch')
        if isinstance(v,m.Relational):
            x,y=map(ev,v.args)
            return {'==':x==y,'!=':x!=y,'<':x<y,'<=':x<=y,'>':x>y,'>=':x>=y}[v.rel_op]
        if v.func in (s.And,s.Or):return (all if v.func==s.And else any)(ev(q) for q in v.args)
        if isinstance(v,m.AppliedUndef):sig=v.func,(0,)*len(v.args),v.args
        else:sig=m.original.partial_signature(v) if isinstance(v,(s.Subs,s.Derivative)) else None
        if sig:
            head,orders,args=sig
            return jet(head.__name__,len(args),orders,self.form,self.stratum,self.point)(*(ev(q) for q in args))
        raise m.Unsupported('no precision evaluator for '+str(v.func))


def refine_closed(a,b,sample):
    form,stratum,point,_=sample['case']
    values=[]
    for digits in (45,75):
        with mp.workdps(digits):
            ctx=PrecisionSample(form,stratum,point,sample['exponents'])
            x,y=ctx.value(a),ctx.value(b)
            if not all(mp.isfinite(v) and not mp.im(v) for v in (x,y)):raise m.Unsupported('non-real or nonfinite result')
            values.append([mp.nstr(q,digits) for q in (x,y,x-y)])
    with mp.workdps(80):
        x,y,diff=map(mp.mpf,values[-1])
        change=max(abs(mp.mpf(old)-mp.mpf(new)) for old,new in zip(*values))
        estimate=8*change+mp.mpf('1e-34')*max(abs(x),abs(y),1)
        status='within_precision_estimate' if abs(diff)<=estimate else 'persistent_difference'
    return dict(case=sample['case'],exponents=sample['exponents'],old_residual=sample['residual'],
                values_and_residual=values,estimate=mp.nstr(estimate,12),status=status)


def hyperplane(rel):
    poly=s.Poly(rel.lhs-rel.rhs,*m.EXP)
    if poly.total_degree()!=1:raise m.Unsupported('nonlinear boundary')
    return tuple(poly.coeff_monomial(q)/poly.LC() for q in (*m.EXP,1))


def branch_sets(raw,engine):
    terms=([m.parse(x,'py') for x in ast.literal_eval(raw)] if engine=='py' else list(m.parse(raw,'wl')))
    planes=set()
    joints=set()
    for q in terms:
        if q in (True,False):continue
        clauses=q.args if q.func==s.And else (q,)
        if not all(isinstance(rel,(s.Equality,s.Unequality)) for rel in clauses):raise m.Unsupported('unexpected Boolean exponent clause')
        planes.update(map(hyperplane,clauses))
        if len(clauses)>1:
            if any(rel.func!=s.Equality for rel in clauses):raise m.Unsupported('mixed compound branch')
            joints.add(tuple(sorted(set(map(hyperplane,clauses)))))
    return planes,joints


def points_for(planes,point):
    base=next(m.exponent_samples(point,()))[1]
    matrix=s.Matrix([list(q[:-1])+[-q[-1]] for q in planes])
    reduced,pivots=matrix.rref()
    if 4 in pivots:raise ValueError('inconsistent branch')
    vals=list(base)
    for r,j in enumerate(pivots):
        vals[j]=reduced[r,-1]-sum(reduced[r,k]*base[k] for k in range(4) if k not in pivots)
    return vals


def evaluate_extra(a,b,planes,label):
    out=[]
    for form,stratum,point in itertools.product(m.FORMS,m.STRATA,range(len(m.POINTS))):
        exp=points_for(planes,point)
        sample=dict(case=[form,stratum,point,label],exponents=list(map(str,exp)))
        try:
            values=[]
            for order in m.RULES:
                ctx=m.Sample(form,stratum,point,exp,order)
                values.append([float(ctx.value(v)) for v in (a,b)])
            x,y=values[-1]
            scale=max(1.,abs(x),abs(y))
            deltas=[max(abs(values[j][i]-values[j-1][i]) for i in (0,1)) for j in (1,2)]
            tol=8*max(deltas)+1.e-9*scale
            if not all(mp.isfinite(v) for pair in values for v in pair):raise ValueError('nonfinite')
            sample.update(values=values,residual=x-y,tolerance=tol,status='separated' if abs(x-y)>tol else
                          'poorly_resolved' if deltas[-1]>1.e-6*scale else 'within_estimate')
        except (ValueError,TypeError,ZeroDivisionError,m.Unsupported) as exc:sample.update(status='unresolved',reason=str(exc))
        out.append(sample)
    return out


def stationary_case(expr,case):
    _,stratum,_,_=case
    t=s.Symbol('t')
    heads=[]
    if stratum in ('X2_velocity_zero','X2_and_radii_velocity_zero'):heads += [f'X_2_{j}' for j in (1,2,3)]
    if stratum in ('radii_velocity_zero','X2_and_radii_velocity_zero'):heads += ['R','r_in']
    if stratum=='X1_velocity_zero':heads += [f'X_1_{j}' for j in (1,2,3)]
    # Exact substitutions of first derivatives only, after native derivatives have been bound.
    return expr.xreplace({s.Derivative(s.Function(name)(t),t):s.S.Zero for name in heads})


def work(idx):
    obj=json.load(gzip.open(base/'objects'/f'{idx:03d}.json.gz'))
    with ExitStack() as files:
        handles={k:files.enter_context(v.open('rb')) for k,v in c.PATHS.items()}
        cmp=c.ObjectComparison(handles)
        names={q['field']:q for q in obj['rows'] if q['role']=='cross'}
        row=names['item3/exponent_branch_inventory']
        extra=[]
        branch=None
        if row['status']=='unresolved':
            py,wl=[branch_sets(cmp.text(r),e) for r,e in zip((row['left'][0],row['right'][0]),('py','wl'))]
            loci=sorted(py[0]|wl[0])
            joints=sorted(py[1]|wl[1])
            branch=dict(same_elementary_boundaries=py[0]==wl[0],py_explicit_conjunctions=len(py[1]),
                        wl_explicit_conjunctions=len(wl[1]),single_planes=[list(map(str,q)) for q in loci],
                        joint_planes=[[list(map(str,p)) for p in group] for group in joints])
            extra=[('single_'+str(i),(p,)) for i,p in enumerate(loci)]+[('joint_'+str(i),group) for i,group in enumerate(joints)]
        refined=[]
        for f in ('item3/object','item3/Q_closed_operand','item3/family_integral','velocity_part_2/residual'):
            if f not in names:continue
            row=names[f]
            parent=obj['rows'][row['same_comparison_as']] if 'same_comparison_as' in row else row
            bad=[sample for sample in parent.get('samples',[]) if sample['status']=='separated']
            evaluate=extra if f.startswith('item3/') else []
            if not bad and not evaluate:continue
            a,b=map(cmp.operand,(row['left'][0],row['right'][0]))
            add=[x for label,planes in evaluate for x in evaluate_extra(a,b,planes,label)]
            bad+= [x for x in add if x['status']=='separated']
            checked=[]
            for sample in bad:
                if f=='velocity_part_2/residual':
                    xs=[stationary_case(v,sample['case']) for v in (a,b)]
                    residual=m.exact_residual(*xs)
                    checked.append(dict(case=sample['case'],old_residual=sample['residual'],
                                        status='exact_on_stratum' if residual==0 else 'unresolved',
                                        both_operands_zero=all(v==0 for v in xs),exact_residual='0' if residual==0 else 'unreduced'))
                elif not (a.has(s.Integral) or b.has(s.Integral)):
                    try:checked.append(refine_closed(a,b,sample))
                    except (m.Unsupported,ValueError,TypeError,ZeroDivisionError) as exc:
                        checked.append(dict(case=sample['case'],status='unresolved',reason=str(exc)))
                else:checked.append(dict(case=sample['case'],status='unresolved',reason='high precision integration unverified'))
            refined.append(dict(field=f,left=row['left'],right=row['right'],additional_samples=add,rechecks=checked))
    result=dict(object=obj['object'],build_id=obj['build_id'],branch=branch,fields=refined)
    c.save(output/f'{idx:03d}.json.gz',result)
    print('refined',idx,obj['object'],flush=True)


def main():
    global base,output
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('root',type=Path)
    opt=ap.parse_args()
    sys.setrecursionlimit(10000)
    base=opt.root
    output=base/'refinement'
    output.mkdir(exist_ok=False)
    config=json.loads((base/'configuration.json').read_text())
    c.PATHS={k:c.inputs.ROOT/v['path'] for k,v in config['inputs'].items()}
    c.LOCATOR={r['id']:(r['id'].split(':')[0],r['span'],r['payload'])
               for r in json.load(gzip.open(base/'coverage.json.gz')) if 'payload' in r}
    jobs=[]
    for p in sorted((base/'objects').glob('*.json.gz')):
        obj=json.load(gzip.open(p))
        if any(q['status'] in ('sample_difference','unresolved') for q in obj['rows']):jobs.append(obj['index'])
    with ProcessPoolExecutor(len(os.sched_getaffinity(0)),mp_context=multiprocessing.get_context('fork')) as pool:
        list(pool.map(work,jobs))
    tally=defaultdict(Counter)
    for p in output.glob('*.gz'):
        r=json.load(gzip.open(p))
        for field in r['fields']:
            for x in field['rechecks']:tally[field['field']][x['status']]+=1
            for x in field['additional_samples']:tally[field['field']+'/new'][x['status']]+=1
    c.save(output/'summary.json',dict(matched_objects=len(jobs),counts=tally))


if __name__=='__main__':
    main()
