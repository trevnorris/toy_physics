"""Count joint exponent specializations of printed elementary additive loci, separately from explicit And syntax."""
from native_challenges import member,tf,obj,findrow,s
from collections import Counter
from functools import lru_cache
from itertools import combinations
from fractions import Fraction as F
@lru_cache(None)
def strata(planes):
    planes=tuple(tuple(map(s.Rational,n)) for n in planes)
    def canonical(equations):
        if not equations:return ()
        rr,piv=s.Matrix(equations).rref()
        return None if 4 in piv else tuple(tuple(rr[i,:]) for i in range(len(piv)))
    known={()}
    for q in planes:
        extended={canonical(h+(q,)) for h in known};known.update(extended-{None})
    result={}
    for h in known:
        sig=tuple(i for i,q in enumerate(planes) if canonical(h+(q,))==h)
        result[sig]=len(h)
    return result
def signature(planes,v):
    return tuple(i for i,p in enumerate(planes) if F(p[-1])+sum(F(a)*F(b) for a,b in zip(p,v))==0)
def main():
    hist=Counter();rank=Counter();tot=Counter()
    for i in range(336):
        o=obj(i); row=findrow(o,'item3/object')
        if row['status'] not in ('sample_agreement','sample_difference'): continue
        row=o['rows'][row['same_comparison_as']] if 'same_comparison_as' in row else row
        branch={}
        extra=[]
        key=f'refinement/{i:03}.json.gz'
        if key in tf.getnames():
            data=member(tf,key);branch=data['branch']
            ex=next((q for q in data['fields'] if q['field']=='item3/object'),{})
            extra=ex.get('additional_samples',[])
        p=tuple(tuple(k) for k in (branch['single_planes'] if branch else o['exponent_loci']))
        if not p:continue
        truth=strata(p)
        chosen={tuple(q['exponents']) for q in row['samples']+extra}
        hit={signature(p,v) for v in chosen}
        unseen=set(truth)-hit
        hist[len(p),len(truth),len(hit),len(unseen)]+=1
        for state in unseen: rank[truth[state]]+=1
        tot.update(objects=1,nonempty_gap=bool(unseen),missing=len(unseen),observed=len(hit))
    print('final item3/object sample cases with native hyperplanes; count is per object, syntactic joint branch and admissible intersection differ')
    print('planes feasible truth-vectors sampled absent -> number of objects',dict(sorted(hist.items())))
    print('missing by codimension',dict(rank),'totals',dict(tot))
if __name__=='__main__':main()
