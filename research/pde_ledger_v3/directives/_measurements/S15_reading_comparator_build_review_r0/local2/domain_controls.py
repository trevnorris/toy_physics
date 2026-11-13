"""Ablate standalone context records while other operands stay byte-for-byte the same."""
import sys,tarfile,gzip,itertools,collections,json
from streams import PACKET, raw, scan, files, text_at, detail
from native_controls import m,c
from S15_reading_context_checks import Catalogue, kernel_checks, angular_moments
import sympy as s

def main():
    ix=scan()
    with tarfile.open(PACKET/'scripts/out/S15_reading_comparator_details.out') as tar:
        stored=detail(tar,'context_checks.json.gz')
        rows=stored['kernels']
        objs=detail(tar,'configuration.json')['required_objects']
        archive_result=next(q for q in rows if q['field'].startswith('positive endpoint 0'))
        print('accepted lower positive endpoint',archive_result)
        full=detail(tar,'context.json.gz')
        for v in full:
            if v['field'] in ('profile_family/chart_condition','family/centre_substitution'):print('accepted context',v)
        kidx=detail(tar,'coverage.json.gz')
    byid={r['id']:r for r in kidx}
    # Read all operands from the hash-verified native files, matched by id in the accepted run.
    class Cat(Catalogue):
        def __init__(self, changes):
            self.rows=[];self.covered=set();self.tags={};self.blobs={}
            for e in ('py','wl'):
                k=e+'_reading'
                for tag,rr in ix[k].items():
                    if not tag.startswith('item3_kernel/') or tag.endswith('/comparison'):continue
                    for i,ref in enumerate(rr):
                        rid=f'{k}:{ref[0]}'
                        rec=byid[rid];self.tags.setdefault((e,tag),[]).append(rec)
                        if (e,tag,i) in changes:
                            val=text_at(k,ref)
                            old,new=changes[e,tag,i]
                            assert old in val
                            self.blobs[rec['payload']]=val.replace(old,new,1)
                        elif rec['payload'] not in self.blobs:self.blobs[rec['payload']]=text_at(k,ref)
    cat=Cat({('py','item3_kernel/distance_endpoints',0):('-integration_radius + kernel_distance','integration_radius - kernel_distance')})
    kernel_checks(cat)
    print('mutated lower endpoint',cat.get('py','item3_kernel/distance_endpoints')[1],
          next(q for q in cat.rows if q['field'].startswith('positive endpoint 0')))
    cat=Cat({('py','item3_kernel/azimuth_operand',0):('(az, 0, 2*pi))','(az, 0, pi))')})
    angular_moments(cat)
    print('mutated first angular bound, original (az, 0, 2*pi) -> (az, 0, pi):',cat.rows[0]['field'],cat.rows[0]['status'], 'number of moments',len(cat.rows))
    py={t:[('py',t)] for t in ix['py_reading']}
    wl={t:[('wl',t)] for t in ix['wl_reading']}
    def family(changes):
        def get(rec):
            k,t=rec;rv=raw(ix,k+'_reading',t)
            if rec in changes:
                old,new=changes[rec];assert old in rv
                return rv.replace(old,new,1)
            return rv
        return list(c.family_context(py,wl,get))
    for tag,mutations in [
        ('profile_family/chart_condition',{('wl','profile_family/domain'):('Reals','Integers')}),
        ('family/centre_substitution',{('py','family/centre_substitution'):(' + X_1_1(t),',' + X_2_1(t),')})]:
        print('changed context',mutations,next(row for row in family(mutations) if row['field']==tag))
    flags=collections.Counter()
    maxima={}
    for form,stratum,point in itertools.product(m.FORMS,m.STRATA,range(len(m.POINTS))):
        targets=[]
        if stratum in ('X2_velocity_zero','X2_and_radii_velocity_zero'):targets += ['X_2_1','X_2_2','X_2_3']
        if stratum in ('radii_velocity_zero','X2_and_radii_velocity_zero'):targets+=['r_in','R']
        if stratum=='X1_velocity_zero': targets+=['X_1_1','X_1_2','X_1_3']
        for name in targets:
            val=m.Sample(form,stratum,point,(1,2,3,4),10).value(s.diff(s.Function(name)(s.Symbol('t')),s.Symbol('t')))
            flags['total']+=1
            if val!=0:
                flags['nonzero']+=1
                if stratum not in maxima or abs(val)>abs(maxima[stratum][-1]):maxima[stratum]=(name,form,point,val)
        if targets:flags['cases']+=1
    print('floating stationary substitutions',flags,maxima)
if __name__=='__main__':main()
