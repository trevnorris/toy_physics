"""All fixed sorted modal actions on the reciprocal angular spring network.

Deleting springs to form a tree gives a gauge-independent lower bound on
every sorted squared frequency, because link phases on a tree can be
removed by site rephasings. Signed twists and site mixtures are unrestricted.
Eigenvalue comparisons use checked dyadic residuals and orthogonality.
"""
import argparse

import numpy as np
from mpmath import iv
from scipy import sparse
from scipy.sparse.csgraph import minimum_spanning_tree

from mouth_cavity_occupation import S
from mouth_collar_link_bounds import R, PI, rounded_bound
from mouth_collar_mobile_rim import prepare


SCALE=2**25


def dyad(x):
    k=round(float(x.mid)*SCALE)
    assert abs(x-R(k)/SCALE)<R(1)/SCALE
    return k


def stiffness(p,edges,tau,delta):
    n=len(p)
    h=np.zeros((2*n,2*n),dtype=np.int64)
    im=np.zeros_like(h)
    deg=np.zeros(n,dtype=int)
    for i,j in edges:
        deg[i]+=1
        deg[j]+=1
        h[i,j]=h[j,i]=h[n+i,n+j]=h[n+j,n+i]=dyad(-tau)
    for i,(x,y,z) in enumerate(p):
        h[i,i]=dyad(int(deg[i])*tau+delta*(1-z))
        h[i+n,i+n]=dyad(int(deg[i])*tau+delta*(1+z))
        h[i,n+i]=h[n+i,i]=dyad(-delta*x)
        im[i,n+i]=dyad(delta*y)
        im[n+i,i]=-im[i,n+i]
    # Row-sum bound on true H minus its dyadic representation.
    return h,im,R(int(max(deg))+3)/SCALE


def norm_majorant(re,im):
    absolute=np.abs(re)+np.abs(im)
    assert int(absolute.max())*max(absolute.shape)<2**62
    return max(int(absolute.sum(axis=0).max()),int(absolute.sum(axis=1).max()))


def sorted_enclosures(p,edges,tau,delta):
    h,b,material_error=stiffness(p,edges,tau,delta)
    H=(h+1j*b)/SCALE
    values,vec=np.linalg.eigh(H)
    D=np.rint(values*SCALE).astype(np.int64)
    X=np.rint(vec.real*SCALE).astype(np.int64)
    Y=np.rint(vec.imag*SCALE).astype(np.int64)
    # Intermediate products are exact signed 64-bit integers within these bounds.
    n=len(D)
    limit=max(int(np.abs(m).max()) for m in (h,b,D,X,Y))+1
    assert (2*n+2)*limit**2 < 2**62
    xt,yt=np.ascontiguousarray(X.T),np.ascontiguousarray(Y.T)
    gram_real=xt@X+yt@Y
    gram_real.flat[::n+1]-=SCALE**2
    gram_imag=xt@Y-yt@X
    g=R(norm_majorant(gram_real,gram_imag))/SCALE**2
    hr,hi=sparse.csr_matrix(h),sparse.csr_matrix(b)
    rr=hr@X-hi@Y-X*D
    ri=hi@X+hr@Y-Y*D
    residual=R(norm_majorant(rr,ri))/SCALE**2
    assert g<R(".01")
    # Let U be the unitary polar factor of (X+iY)/S. ||U-Q||<=g.
    # Weyl applied after comparison with U diag(D/S) U* gives all ranks.
    norms=R(norm_majorant(h,b)+int(np.abs(D).max()))/SCALE
    error=(residual+norms*g+material_error).b
    return [R(int(x))/SCALE for x in D],error


def trees(n,edges):
    e=np.array(edges)
    i,j=e.T
    for seed in (175,389):
        w=np.random.default_rng(seed).uniform(.5,1.5,len(i))
        graph=sparse.csr_matrix((w,(i,j)),shape=(n,n))
        tree=minimum_spanning_tree(graph+graph.T).tocoo()
        assert tree.nnz==n-1
        yield list(zip(tree.row,tree.col))


def run(level):
    cfg=S["parameters"]
    ref,L=[R(str(v)) for v in S["phase_A"]["q_star"]]
    c,mi,mo=[R(str(cfg[v])) for v in ("c_w","mu_in","mu_out")]
    axial=(c*PI/2)**2/L**3
    # Antipodal vertex permutation exchanges mouth labels and preserves the
    # full graph. Both mirrors therefore have identical sorted-sector bounds.
    p,edges,adj,nf=prepare(level,(".26","-.19",".17"),1)
    forests=list(trees(len(p),edges))
    for ratio in (".9","1","1.125","1.25","1.5","2"):
        scale=R(ratio)
        radius=scale*ref
        tau=R(".04")/scale**2
        delta=R(".16")
        base=mi**2+(c*PI/radius)**2+(c*PI/(2*L))**2
        full,err=sorted_enclosures(p,edges,tau,delta)
        lowers=[]
        for tree in forests:
            vals,tree_err=sorted_enclosures(p,tree,tau,delta)
            lowers.append([(v-tree_err).a for v in vals])
        # Independent triangle bound: 2L0 is the sum of all untwisted
        # triangular Laplacians, each bounded by 3 on its face. La>=0.
        universal=9*tau
        ratio_bound=iv.sqrt(base+universal)/(iv.sqrt(base+universal)+iv.sqrt(base))
        analytic=universal*ratio_bound/axial
        gains=[]
        for k,val in enumerate(full):
            if base+val-err >= mo**2:
                continue
            top=min((base+val+err).b,(mo**2).b)
            bot=(base+max(R(0).a,*(row[k] for row in lowers))).a
            wtop=iv.sqrt(top)
            wbot=iv.sqrt(bot)
            gains.append((wtop*(wtop-wbot)/axial).b)
        bound=min(analytic.b,max(gains))
        cost=PI/(32*scale)
        load=4*PI*R(str(cfg["sigma"]))*radius**2+4*PI*R(str(cfg["alpha_volume"]))*radius**3/3
        print(f"forest level={level} R/R0={ratio} candidate_ranks={len(gains)} "
              f"flat_spectrum_error_upper={rounded_bound(err,8,True):.8f} "
              f"all_occupations_gain_upper={rounded_bound(load*bound,8,True):.8f} "
              f"exterior_cost_lower={rounded_bound(cost,8):.8f} "
              f"analytic_only_gain_upper={rounded_bound(load*analytic,8,True):.8f} "
              f"any_occupation_exclusion_below_load={rounded_bound(cost/bound,6):.6f}",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--level",type=int,choices=[2,3],required=True)
    run(parser.parse_args().level)
