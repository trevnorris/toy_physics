"""Read the scoped comparator catalogue to compare endpoint conventions and argument checks.

Run through s11c_guarded_run.py after the object pass; never launches either engine.
"""
import argparse
import ast
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import sympy as s
import S15_reading_math as m
from S15_reading_comparator import save, InputError


class Catalogue:
    def __init__(self,root):
        with gzip.open(root/'coverage.json.gz','rt') as f:
            records=json.load(f)
        wanted={}
        self.tags=defaultdict(list)
        for record in records:
            stream=record['id'].split(':')[0]
            if not stream.endswith('_reading') or 'payload' not in record:
                continue
            tag=record['tag']
            needed=(tag.startswith('item3_kernel/') or 'exponent_divisor_loci' in tag
                    or tag in ('reading/dimension_argument_audit','field_argument_dimension_signatures',
                               'family/response_dimension_signatures','dimension_field_signatures'))
            if needed and not tag.endswith('/comparison'):
                self.tags[stream[:2],tag].append(record)
                wanted[record['payload']]=None
        with gzip.open(root/'operands.jsonl.gz','rt') as f:
            for line in f:
                row=json.loads(line)
                if row['sha256'] in wanted:
                    wanted[row['sha256']]=row['native_operand']
        self.blobs=wanted
        self.rows=[]
        self.covered=set()

    def get(self,engine,tag,index=0):
        r=self.tags[engine,tag][index]
        self.covered.add(r['id'])
        return r['id'],self.blobs[r['payload']]

    def kernel(self,tag,engine,index=0):
        rid,raw=self.get(engine,'item3_kernel/'+tag,index)
        parse=m.original.py_parse if engine=='py' else m.original.wl_parse
        if engine=='wl':
            raw=m.WL_SCALARS.sub(lambda match:m.CHART[match[0].replace(' ','')],raw)
        return rid,parse(raw)

    def note(self,field,py,wl,a,b,unit=False):
        raw=a-b
        reduced=s.factor(raw)
        row=dict(field=field,left=py,right=wl,translated_py=str(a),translated_wl=str(b),
                 as_computed=str(raw),factored=str(reduced))
        if unit and reduced!=0:
            # Eliminate one squared component only after printing the raw difference.
            num,den=s.fraction(s.cancel(reduced))
            direction=m.DIRECTION
            relation=sum(q*q for q in direction)-1
            reduced=s.cancel(s.rem(num,relation,direction[0])/den)
            row.update(unit_reduction=str(reduced),equation='sum N_i^2 = 1, profile_family/chart_condition')
        if field.startswith('log_primitive') and reduced.has(s.polylog):
            reduced=s.simplify(reduced.replace(lambda v:v.func==s.polylog and v.args[0]==1,
                                               lambda v:-s.log(1-v.args[1])))
            row.update(real_branch_reduction=str(reduced),equation='Li1(z) = -log(1-z) for the principal real branch at |z| < 1; native log variable is r/D.')
        row['status']='exact_identity' if reduced==0 else 'remaining_symbolic_difference'
        self.rows.append(row)


def kernel_checks(cat):
    py_names=('kernel_distance','integration_radius','squared_distance','mu','kernel_power','log_argument')
    wl_names=('distance','radial','squared','mu','power','logVariable')
    z=s.symbols('kernel_D kernel_r kernel_w kernel_mu kernel_p kernel_log',real=True)
    scopes=dict(py=tuple(zip(py_names,z)),wl=tuple(zip(wl_names,z)))
    def value(node,engine):
        return m.translate(node,engine,scopes[engine])
    def pair(a,b,pa=lambda v:v,pb=lambda v:v,index=0,label=None):
        i,x=cat.kernel(a,'py',index)
        j,y=cat.kernel(b,'wl',index)
        x,y=value(pa(x),'py'),value(pb(y),'wl')
        cat.note(label or a,i,j,x,y)
    tag='distance_substitution'
    pair('squared_distance',tag,pb=lambda n:n.args[0])
    pair('change_of_variable',tag,lambda n:n.args[0],lambda n:n.args[1].args[0].args[1],label='inverse angle')
    pair('change_of_variable',tag,lambda n:n.args[1],lambda n:n.args[2],label='positive distance Jacobian')
    for k in (0,1):
        pair('change_of_variable',tag,lambda n:n.args[2].args[k],lambda n:n.args[3].args[k],label=f'squared endpoint {k}')
        i,x=cat.kernel('distance_endpoints','py')
        j,y=cat.kernel(tag,'wl')
        cat.note(f'positive endpoint {k}, squared on D > r > 0',i,j,value(x.args[k],'py')**2,value(y.args[3].args[k],'wl'))
    pair('power_primitive_operand','power_primitive_operand')
    i,x=cat.kernel('power_primitive','py')
    j,y=cat.kernel('power_primitive_branches','wl')
    py=value(x,'py')
    # Explicitly evaluate guard selection. Pick the generic branch only after checking
    # the native opposing guards; otherwise a reversed inequality could pass unnoticed.
    if not (py.args[0].cond==m.original.native_relation('Unequality',z[4],-s.S.One)
            and value(y.args[0].args[1],'wl')==m.original.native_relation('Equality',z[4],-s.S.One)
            and py.args[1].cond==s.true and value(y.args[1].args[1],'wl')==True):
        raise InputError('unrecognized power primitive partition')
    cat.note('power primitive, exponent different from -1',i,j,py.args[0].expr,value(y.args[1].args[0],'wl'))
    cat.note('power primitive, exponent -1',i,j,py.args[1].expr,value(y.args[0].args[0],'wl'))
    for k in (0,1):
        for name in ('log_primitive_operand','log_primitive','log_primitive_derivative/operand_a',
                     'log_primitive_derivative/operand_b','log_primitive_derivative/residual'):
            pair(name,name,index=k,label=f'{name}, orientation {(-1,1)[k]}')
    for a,b in [('geometric_distance_derivative','geometric_distance_derivative'),
                ('separation_time_derivative','separation_rate')]:
        pair(a,b)
    for a,b in [('geometric_direction_derivative','geometric_direction_derivative'),
                ('direction_time_derivative','direction_rate')]:
        i,x=cat.kernel(a,'py')
        j,y=cat.kernel(b,'wl')
        xs=[v.args[0] for v in x.args[1].args] if m.original.call(x,'Matrix') else x.args
        for k,(u,v) in enumerate(zip(xs,y.args,strict=True)):
            cat.note(f'{a} component {k+1}',i,j,value(u,'py'),value(v,'wl'),unit=True)


def angular_moments(cat):
    """Recover the exponent tuple from the emitted SymPy operand before comparing its value."""
    axis=m.DIRECTION
    mu,az=s.symbols('kernel_mu kernel_az',real=True)
    py_scope=tuple(zip(('axis_1','axis_2','axis_3','mu','az'),(*axis,mu,az)))
    wl_scope=tuple(zip(('axisX','axisY','axisZ','mu'),(*axis,mu)))
    h=s.Matrix([0,0,1])-s.Matrix(axis)
    rotated=(s.eye(3)-h*h.T/(1-axis[2]))*s.Matrix(
        [s.sqrt(1-mu**2)*s.cos(az),s.sqrt(1-mu**2)*s.sin(az),mu])
    wtag='angular_moment'
    candidates={}
    for k in range(len(cat.tags['wl','item3_kernel/'+wtag])):
        rid,node=cat.kernel(wtag,'wl',k)
        powers=tuple(int(v.value) for v in node.args[0].args)
        if powers in candidates:raise InputError('duplicate Wolfram angular powers')
        candidates[powers]=(rid,m.translate(node.args[1],'wl',wl_scope))
    # Exact rational points identify a candidate efficiently. Confirm the full
    # native integrand symbolically before accepting its binding.
    numbered=[(tuple(s.Rational(n,3) for n in ax),s.Rational(3,5),s.Rational(sn,5),s.Rational(cs,5))
              for ax,sn,cs in [((2,1,2),3,4),((-2,2,1),-4,3)]]
    subs=[dict(zip((*axis,mu,s.sin(az),s.cos(az)),(*a,v,sn,cs))) for a,v,sn,cs in numbered]
    def key(v):return tuple(v.xreplace(row) for row in subs)
    options={}
    for p,(rid,value) in candidates.items():
        monomial=s.prod(arg**power for arg,power in zip(rotated,p))
        options.setdefault(key(monomial),[]).append((p,rid,value,monomial))
    assigned=set()
    for k in range(len(cat.tags['py','item3_kernel/azimuth_operand'])):
        input_id,integral=cat.kernel('azimuth_operand','py',k)
        body=integral.args[1]
        expr=m.translate(body,'py',py_scope)
        # Source emits operand and result as adjacent records with intervening /comparison labels.
        answer_id,ans=cat.kernel('azimuth_object','py',k)
        opts=options.get(key(expr),[])
        verified=[candidate for candidate in opts if s.cancel(expr-candidate[3])==0]
        if len(verified)!=1:
            raise InputError('angular binding failed '+input_id)
        powers,wl_id,wl_result,monomial=verified[0]
        if powers in assigned:
            raise InputError('unverified angular power tuple '+str(powers))
        assigned.add(powers)
        cat.note('angular moment '+str(powers),answer_id,wl_id,m.translate(ans,'py',py_scope),wl_result,unit=True)
        cat.rows[-1].update(binding_operand=input_id,integrand_minus_rotated_monomial='0',powers=powers)
    if assigned!=set(candidates):raise InputError('angular roster mismatch')


def audit_arguments(cat,root):
    with gzip.open(root/'context.json.gz','rt') as f: rows=json.load(f)
    tag='field_argument_dimension_signatures'
    rid,raw=cat.get('wl',tag)
    signatures={}
    for pair in m.original.wl_parse(raw).args:
        n,v=pair.args
        signatures[m.NAMES.name(m.original.key_text(n),'wl')]=[
            [int(m.translate(k,'wl',())) for k in vector.args] for vector in v.args]
    out=[]
    for r in rows:
        if r['field']!='reading/dimension_argument_audit':continue
        key=r['argument']
        a=key[0]
        row=dict(py_record=r['left'],argument=key,expected=r['expected'],observed=r['observed'])
        if len(key)==3 and isinstance(key[1],int):
            # The Wolfram sources register change[h] with argumentDims[h]. Record that
            # separate source-defined inheritance, not a second observed-argument audit.
            base=a.removesuffix('_change') if a not in signatures else a
            if base in signatures:
                expected=signatures[base][key[1]]
                row.update(wl_signature_record=rid,wl_head=base,wl_slot=key[1]+1,wl_expected=expected,
                           inherited_from_base=(a!=base),status='matches_declared_slot' if
                           [expected]==row['expected']==row['observed'] else 'slot_dimension_difference')
            else:row.update(status='no_wl_declared_slot')
        else:
            row.update(status='sympy_only_bound_or_trig_check')
        out.append(row)
    return out


def divisor_audit(cat):
    divisors=[]
    for (engine,tag),rows in list(cat.tags.items()):
        if engine!='wl' or not tag.endswith('/exponent_divisor_loci'):continue
        for idx,record in enumerate(rows):
            rid,raw=cat.get(engine,tag,idx)
            payload=m.original.wl_parse(raw)
            flags=[q.args[1] for q in payload.args]
            false=sum(q.head=='Boolean' and q.value=='False' for q in flags)
            divisors.append(dict(id=rid,tag=tag,branch_supports=len(flags),empty_returns=false,
                                 other_returns=len(flags)-false))
    return divisors


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('root',type=Path)
    ap.add_argument('--divisors-only',action='store_true',help='Recheck the native Boolean return tags without replaying completed comparisons.')
    args=ap.parse_args()
    cat=Catalogue(args.root)
    if args.divisors_only:
        path=args.root/'context_checks.json.gz'
        with gzip.open(path,'rt') as f: data=json.load(f)
        data['wolfram_divisor_loci']=divisor_audit(cat)
        save(path,data)
        ds=data['wolfram_divisor_loci']
        print('Wolfram divisor records',len(ds),'supports',sum(v['branch_supports'] for v in ds),
              'non-False returns',sum(v['other_returns'] for v in ds))
        return
    kernel_checks(cat)
    angular_moments(cat)
    arguments=audit_arguments(cat,args.root)
    divisors=divisor_audit(cat)
    remaining={}
    for (engine,tag),rows in cat.tags.items():
        if tag.startswith('item3_kernel/') and not tag.endswith('/comparison'):
            ids=[r['id'] for r in rows if r['id'] not in cat.covered]
            if ids:remaining[engine+' '+tag]=ids
    data=dict(kernels=cat.rows,argument_slots=arguments,wolfram_divisor_loci=divisors,
              unpaired_internal_kernel_occurrences=remaining,
              note='Diagnostic blocks have engine-specific definitions and multiplicities. Unpaired IDs refer to unchanged catalogue occurrences.')
    save(args.root/'context_checks.json.gz',data)
    print('typed kernel comparisons',dict(Counter(x['status'] for x in cat.rows)))
    print('argument slot comparisons',dict(Counter(x['status'] for x in arguments)))
    print('Wolfram divisor records',len(divisors),'supports',sum(x['branch_supports'] for x in divisors),
          'non-False returns',sum(x['other_returns'] for x in divisors))


if __name__=='__main__':
    main()
