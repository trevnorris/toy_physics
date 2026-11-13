"""Directly compare archived arguments, error estimates and within-engine family values."""
import tarfile,collections,itertools,sys
from streams import PACKET,detail
from native_controls import m
import sympy as s
def run():
    with tarfile.open(PACKET/'scripts/out/S15_reading_comparator_details.out') as tar:
        config=detail(tar,'configuration.json')
        raw_counts=collections.defaultdict(collections.Counter)
        own=collections.Counter()
        worst=(0,None)
        for i,name in enumerate(config['required_objects']):
            obj=detail(tar,f'objects/{i:03d}.json.gz')
            rows=obj['rows']
            fields={r['field']:rows[r.get('same_comparison_as',j)] for j,r in enumerate(rows) if r['role']=='cross'}
            for field in ('Q','E_geo','E_resp','separation_derivative','item3/family_integral','item3/object','velocity_part_1/residual','velocity_part_2/residual'):
                raw_counts[field][fields[field]['status']]+=1
            if 'item3/Q_closed_operand' not in fields:continue
            pair=[fields[f] for f in ['item3/family_integral','item3/object']]
            if not all('samples' in rr for rr in pair):continue
            objvals={tuple(v['case']):v for v in pair[1]['samples']}
            for src in pair[0]['samples']:
                c=tuple(src['case'])
                ret=objvals.get(c)
                if not (ret and 'values' in src and 'values' in ret):continue
                for j,e in enumerate(('py','wl')):
                    o=src['values'];v=ret['values'];x=o[-1][j];y=v[-1][j]
                    sd=abs(x-y);sc=max(1,abs(x),abs(y))
                    delta=max(abs(q[k][j]-q[k-1][j]) for q in (o,v) for k in (1,2))
                    ratio=sd/(8*delta+1.e-9*sc)
                    own[e,'archived comparisons']+=1
                    if ratio>1:own[e,'initial flags']+=1
                    if ratio>worst[0]:worst=(ratio,(name,e,c,x,y,sd,delta,src['status'],ret['status']))
        print('initial field counts reconstructed from complete separate object files')
        for k,v in raw_counts.items():print(k,dict(v))
        print('within-engine integral vs closed return where both sampled, prior to precision rechecks',dict(own))
        print('worst method-specific binary64 ratio (diagnostic; cancellation in closed form is possible)',worst)
    flags=collections.Counter();zeros=collections.Counter()
    t=s.Symbol('t')
    for form,stratum,p in itertools.product(m.FORMS,m.STRATA,range(len(m.POINTS))):
        h=[]
        if stratum in ('X2_velocity_zero','X2_and_radii_velocity_zero'):h+=['X_2_1','X_2_2','X_2_3']
        if stratum in ('radii_velocity_zero','X2_and_radii_velocity_zero'):h+=['r_in','R']
        if stratum=='X1_velocity_zero':h+=['X_1_1','X_1_2','X_1_3']
        if not h:continue
        ctx=m.Sample(form,stratum,p,(1,2,3,4),10)
        vals=[ctx.value(s.diff(s.Function(k)(t),t)) for k in h]
        if any(vals): flags[stratum]+=1
        else:zeros[stratum]+=1
        print('stationary offsets',form,stratum,p,vals)
    print('off/on intended stationary stratum (6 cases per stratum)',dict(flags),dict(zeros))
if __name__=='__main__':run()
