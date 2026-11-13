from collections import Counter,defaultdict
import ast,math,mpmath as mp,numpy as np,re
import sympy as s
from packet_io import *
import S15_reading_comparator as c
import S15_reading_math as m
import S15_reading_refinement as r
from S15_reading_log_branches import LogBranchSample
import S15_reading_context_checks as ctx
d=Records()
# The packet copy lacks a .git root. Only the accepted git object locator needs
# the checkout: force exactly the git-show dependency expected by the candidate.
c.inputs.ROOT=R;c.inputs.LEDGER=R/'research/pde_ledger_v3'
look={e:defaultdict(list) for e in ('py','wl')}
for (stream,tag),ids in d.n.items():
    if stream.endswith('reading'): look[stream[:2]][tag]=ids
def context(e,tag,frm,to):
    key,=look[e][tag]
    def get(rid):
        v=d.raw(rid)
        if rid==key:
            assert frm in v,(tag,frm)
            return v.replace(frm,to,1)
        return v
    return list(c.family_context(look['py'],look['wl'],get))
for e,tag,frm,to in [
    ('py','family/centre_substitution','X_1_2(t):','X_2_2(t):'),
    ('py','family/centre_substitution',' + X_1_2(t),',' + X_2_2(t),'),
    ('wl','profile_family/separation_chart','xx[1, 2][t] ->','xx[2, 2][t] ->'),
    ('wl','family/separation_vector','xx[1, 2][t] - xx[2, 2][t]','xx[2, 2][t] - xx[1, 2][t]'),
    ('wl','profile_family/domain','Reals','Integers'),('wl','profile_family/domain',' | bb[1]','')]:
    rows=context(e,tag,frm,to);changes=[(v['field'],v['status']) for v in rows if v['status']!='same_values']
    print('native change',e,tag,frm,'=>',to,changes);assert changes
print('stored context argument check gaps',[(x['argument'],x['expected'],x['observed']) for x in member(d.tar,'context_checks.json.gz')['argument_slots'] if x['status']=='no_wl_declared_slot'])
class DirectKernels:
    def __init__(self,edit):
        self.tags=defaultdict(list)
        self.rows=[];self.edit=edit
        for (stream,tag),ids in d.n.items():
            if stream.endswith('_reading') and tag.startswith('item3_kernel/'):
                self.tags[stream[:2],tag]=ids
    def kernel(self,tag,engine,index=0):
        rid=self.tags[engine,'item3_kernel/'+tag][index];raw=d.raw(rid)
        if self.edit and self.edit[:3]==(engine,tag,index):
            old,new=self.edit[3:];assert old in raw;raw=raw.replace(old,new,1)
        if engine=='wl':
            raw=m.WL_SCALARS.sub(lambda mt:m.CHART[mt[0].replace(' ','')],raw)
        return rid,(m.original.py_parse if engine=='py' else m.original.wl_parse)(raw)
    note=ctx.Catalogue.note
for edit,func,field in [
    (('py','distance_endpoints',0, ', integration_radius + kernel_distance',', -integration_radius - kernel_distance'),ctx.kernel_checks,'positive endpoint 1'),
    (('wl','generator_operand',0,'2*Pi','Pi'),ctx.angular_moments,'angular moment'),
    (None,ctx.angular_moments,'angular moment')]:
    k=DirectKernels(edit);func(k);rows=[v for v in k.rows if v['field'].startswith(field)]
    print('kernel ablation',edit,'status',dict(Counter(v['status'] for v in rows)),'example',[(v['field'],v.get('domain_residual',v.get('wl_bounds'))) for v in rows[:1]])
# Native geometry must exercise only specified zeros, retain accelerations,
# and give distinct velocities to each of the two centres on other strata.
names=[f'X_{i}_{j}' for i in (1,2) for j in (1,2,3)]+['r_in','R']
for st in m.STRATA:
    zeros=[];n=0;amin=None
    for form in m.FORMS:
        for p in range(3):
            exp=next(m.exponent_samples(p,()))[1];t=m.POINTS[p][-1]
            z=[];first=[];second=[]
            for nm in names:
                for od,values in [(1,first),(2,second)]:
                    jet,origin=m.sample_jet(nm,1,(od,),form,st,p);v=jet.subs(m.original.FORMAL[0],t-origin)
                    v0=m.jet_function(nm,1,(od,),form,st,p)(float(t))
                    with mp.workdps(70):
                        v1=r.jet(nm,1,(od,),form,st,p)(r.rational(t))
                    assert abs(v0-float(v))<1e-13 and abs(float(v1)-float(v))<1e-13
                    if od==1 and v==0:z.append(nm);assert v0==v1==0;n+=1
                    values.append(str(v))
                    if od==2:
                        assert v;amin=min(float(abs(v)),amin) if amin is not None else float(abs(v))
            zeros.append(tuple(z))
    assert len(set(zeros))==1;print('stationary',st,'zero heads',zeros[0],'num exact zero checks in both backends',n,'minimum live acceleration',amin,'example first second',first,second)
# Native signed zero comparison and both selected centre projections.
for nm in ('region/component_1/momentum_rate','inner/viscous/component_1','outer/transport_centre/component_3'):
    for k in (1,2):
        print('velocity input',nm,'drain',k)
        pe, we=[],[]
        for fld in ('operand_a','operand_b','residual'):
            ids=d.ref(nm,f'velocity_part_{k}/'+fld)
            for rid,eng,store in zip(ids,('py','wl'),(pe,we)):
                raw=d.raw(rid);store.append(m.parse(raw,eng))
                if fld=='residual' and len(raw)<120: print('literal',rid,raw)
        tt=s.Symbol('t')
        drops={s.Derivative(s.Function(f'X_{k}_{j}')(tt),tt):s.S.Zero for j in range(1,4)}
        for eng,ops in zip(('py','wl'),(pe,we)):
            x,z,y=ops
            proofs=[m.exact_residual(x.xreplace(drops),z),m.exact_residual(x-z,y)]
            v=m.Sample(m.FORMS[0],m.STRATA[0],0,next(m.exponent_samples(0,()))[1],32)
            print(eng,[float(v.value(p)) for p in ops],'exact check', [p==0 for p in proofs]);assert all(p==0 for p in proofs)
# Check the one genuinely precision-flagged exceptional printed expression
# against its own family integral, with a direct 48-node spherical integral and
# a 85-digit evaluation of the opposite printed closed form.
name='inner/carried/component_2/labels_drain1+drain2'; o=member(d.tar,'refinement/019.json.gz')
z=next(f['rechecks'][0] for f in o['fields'] if f['rechecks']);form,st,p,_=z['case'];exps=z['exponents']
with mp.workdps(85):
    pc=r.PrecisionSample(form,st,p,exps)
    for fld in ('item3/object','item3/family_integral','item3/Q_closed_operand'):
        vals=[]
        for rid,eng in zip(d.ref(name,fld),('py','wl')):
            a=m.parse(d.raw(rid),eng)
            v=float(LogBranchSample(form,st,p,tuple(map(s.Rational,exps)),48).value(a)) if fld.endswith('family_integral') else pc.value(a)
            vals.append(v)
        print('exceptional grounded by family integral',fld,z['case'],exps,d.ref(name,fld),[mp.nstr(v,62) for v in vals],'difference',mp.nstr(vals[0]-vals[1],6))
# Difference between sample test fields at different inputs/derivative slots.
for form in m.FORMS:
    v=[]
    for name,od,x in [('Pi_o_2_1_change',(0,0,0,0),(1,3,2,1)),('Pi_o_1_2_change',(0,0,0,0),(1,3,2,1)),('W_2_change',(1,0,0,0),(1,2,3,1)),('W_2_change',(0,1,0,0),(1,2,3,1)),('W_2_change',(1,0,0,0),(3,2,1,1))]:
        jet,_=m.sample_jet(name,4,od,form,'generic',0);v.append(jet.subs(dict(zip(m.original.FORMAL,x))))
    print('slot/evaluation-point values',form,list(map(str,v)));assert len(set(v))==5
