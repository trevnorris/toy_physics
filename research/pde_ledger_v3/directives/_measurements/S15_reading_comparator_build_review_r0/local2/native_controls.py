"""One-sided controls on native strings in private memory; no production replay."""
import sys, hashlib, io, itertools
from collections import Counter
from streams import PACKET,scan,raw,files
sys.path.insert(0,str(PACKET/'scripts'))
import S15_reading_comparator as c
import S15_reading_math as m
import sympy as s
def compare(raw1,raw2,field,item,loci=()):
    a,b=raw1.encode(),raw2.encode()
    c.LOCATOR={'left':('py_reading',{'offset':0,'length':len(a)+1},hashlib.sha256(a).hexdigest()),
               'right':('wl_reading',{'offset':0,'length':len(b)+1},hashlib.sha256(b).hexdigest())}
    ctx=c.ObjectComparison({'py_reading':io.BytesIO(a+b'\n'),'wl_reading':io.BytesIO(b+b'\n')})
    ctx.loci=loci
    ctx.compare(field,item,['left'],['right'],'cross')
    return ctx.rows[0]

def describe(label,row):
    print(label,row['status'],row.get('sample_counts'),flush=True)
    if 'samples' in row:
        for rec in [row['samples'][0], row['samples'][-1]]:
            print('case',rec['case'],'values',rec.get('values'),'residual',rec.get('residual'),'tolerance',rec.get('tolerance'),'status',rec['status'],flush=True)
    else:print('residual',row.get('residual'),flush=True)

def test():
    indices=scan()
    for q,field,old,new in [
        ('inner/other/component_1','Q','oo[1, 2]','oo[1, 1]'),
        ('outer/other/component_1','Q','oo[1, 2]','(-1)*oo[1, 2]'),
        ('inner/other/component_1','E_resp','change[oo[1, 2]]','change[oo[1, 1]]'),
        ('both_surfaces/sum_crossing/component_1','E_resp','change[oo[1, 2]]','change[oo[1, 1]]')]:
        tag='reading/'+q+'/'+field
        py,wl=(raw(indices,k+'_reading',tag) for k in ('py','wl'))
        print('mutation',q,field,old,'to',new,'matching original occurrences',wl.count(old),'replaced 1',flush=True)
        if q.startswith('both_surfaces'):
            # Repeat only the corrupted comparison. Baseline for this row was retrieved in streams.stdout.
            pass
        else:
            describe('baseline',compare(py,wl,field,1 if field=='Q' else 2))
        assert old in wl
        desc=compare(py,wl.replace(old,new,1),field,1 if field=='Q' else 2)
        describe('changed',desc)
    q='inner/other/component_1'
    for kind,mut in [('E_geo/dimension_reference_operand',lambda x:x.replace('-2','2')),
                     ('E_resp/dimension_pair_residuals',lambda x:'{{0,0,0}}')]:
        tag=f'reading/{q}/{kind}'
        a,b=(raw(indices,e+'_reading',tag) for e in ('py','wl'))
        print('dimensions',kind,a,b,mut(b),flush=True)
        describe('baseline',compare(a,b,kind,6))
        describe('changed',compare(a,mut(b),kind,6))
if __name__=='__main__':test()
