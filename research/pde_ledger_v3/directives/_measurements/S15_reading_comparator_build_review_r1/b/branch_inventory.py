"""Audit native object/Q guard surfaces against actually exercised points.

For each nontrivial native guard extract its primitive affine equations,
retaining multi-equation branches, and test rational membership only on the
stored per-object sample points. Do not promote coverage to an identity.
"""
import ast,collections,sys
import sympy as s
from packet_io import *
import S15_reading_math as m
from S15_cross_engine_comparator import application,call
sys.setrecursionlimit(20000)
d=Records()
def guards(raw,engine):
    clauses=[]
    if engine=='py':
        tree=ast.parse(raw,mode='eval')
        for v in ast.walk(tree):
            if not isinstance(v,ast.Call) or getattr(v.func,'id',None)!='Piecewise':continue
            for branch in v.args:
                cond=branch.elts[1]
                if isinstance(cond,ast.Constant) and cond.value==True:continue
                clauses.append(ast.unparse(cond))
        return {m.parse(v,engine) for v in set(clauses)}
    seen=set()
    mapped=m.WL_SCALARS.sub(lambda mt:m.CHART[mt[0].replace(' ','')],raw)
    tree=m.original.wl_parse(mapped)
    def walk(n):
        if n in seen:return
        seen.add(n)
        if n.head=='Apply' and call(n,'Piecewise'):
            for entry in n.args[1].args:clauses.append(entry.args[1])
        for child in n.args:walk(child)
    walk(tree)
    return {m.translate(v,engine,()) for v in set(clauses)}
def locus(v):
    parts=v.args if v.func in (s.And,s.Or) else (v,)
    canon=[]
    for q in parts:
        if q in (True,False):continue
        if q.func not in (s.Eq,s.Ne):raise ValueError(q)
        eq=s.Poly(q.lhs-q.rhs,*m.EXP); assert eq.total_degree()==1
        canon.append(tuple(eq.coeff_monomial(u)/eq.LC() for u in (*m.EXP,1)))
    return tuple(sorted(canon))
def evaluate(case,plane):
    return sum(a*b for a,b in zip((*case,1),plane))==0
counts=collections.Counter();miss=[];printed_example=False
all_labels={'single':set(),'joint':set()}
for i,name in enumerate(d.cfg['required_objects']):
    obj=d.obj(i)
    if not any(r['field']=='item3/Q_closed_operand' for r in obj['rows']):continue
    ref_path=f'refinement/{i:03d}.json.gz';new=member(d.tar,ref_path) if ref_path in d.tar.getnames() else {'fields':[]}
    for field in ('item3/object','item3/Q_closed_operand'):
        row=d.row(name,field)
        orig=obj['rows'][row['same_comparison_as']] if 'same_comparison_as' in row else row
        if orig['status']=='exact_identity':continue
        raw=[d.raw(r) for r in d.ref(name,field)]
        gs=[guards(x,e) for x,e in zip(raw,('py','wl'))]
        exposed={locus(g) for s0 in gs for g in s0}; exposed.discard(())
        initial=orig['samples']
        extra=next((f['additional_samples'] for f in new['fields'] if f['field']==field),[])
        samples={tuple(map(s.Rational,v['exponents'])) for v in initial+extra}
        unhit=[l for l in exposed if not any(all(evaluate(x,p) for p in l) for x in samples)]
        if unhit:miss.append((name,field,unhit))
        counts[field]+=1; counts['explicit guard surfaces'+field]+=len(exposed)
        if exposed and len(exposed)>5 and not printed_example:
            print('guard example',name,field, 'count',len(exposed),'example',sorted(exposed)[:3], 'sample exponents',len(samples),flush=True);printed_example=True
        for l in exposed:all_labels['single' if len(l)==1 else 'joint'].add(l)
print('branches',dict(counts),'distinct',[(k,len(v)) for k,v in all_labels.items()],'unhit',len(miss))
for row in miss[:10]:print('missing',row)
