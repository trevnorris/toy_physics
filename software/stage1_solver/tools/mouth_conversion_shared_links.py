"""Regroup the conversion readouts on physical links of one closed cut.

A parcel visits exactly one triangle; each of its three sides is equally
likely to be the marked closing bond. A controlled edge changes the two
adjacent oriented loop angles oppositely. Energy units set hbar=gap=1.
"""
import numpy as np
from scipy.spatial.transform import Rotation
from numpy.polynomial.legendre import leggauss

from mouth_collar_order_preparation import predicted, calibration
from mouth_material_links import section, spherical_areas
from mouth_conversion_twist import read


RATE=.25
TAU=2*np.pi


def prepared(level,r):
    n,faces=section(level)
    n=Rotation.from_rotvec([.26,-.19,.17]).apply(n)
    weights=spherical_areas(n,faces)/(4*np.pi)
    left=faces.ravel()
    right=np.roll(faces,-1,axis=1).ravel()
    edges,inv=np.unique(np.sort(np.column_stack((left,right)),axis=1),
                        axis=0,return_inverse=True)
    rows=np.repeat(np.arange(len(faces)),3)
    sign=np.sign(right-left)
    # face and outward sign for each incidence, grouped by undirected edge
    order=inv.argsort()
    neighbor=rows[order].reshape(-1,2)
    orient=sign[order].reshape(-1,2)
    assert np.all(orient.sum(axis=1)==0) and np.isclose(weights.sum(),1)
    mid=n[edges].sum(axis=1)
    mid/=np.linalg.norm(mid,axis=1)[:,None]
    test_edge=np.argmax(mid@np.array([1.,2.,-1.]))
    mirrored=(predicted(-n,-1,r,6))[faces[:,::-1]]
    return [predicted(n,1,r,6)[faces],mirrored],weights,neighbor,orient,test_edge


def response(p,f):
    # Cyclic relabeling is gauge equivalent in the held tag readout. Its
    # derivative corresponds to moving a different marked bond. Average as
    # stipulated in the parcel law to avoid selecting a preferred corner.
    total_tag=0.
    total_imp=0.
    for k in range(3):
        tag,imp,entry,exit_work=read(np.roll(p,k,axis=1),f)
        assert np.max(abs(tag+entry+exit_work))<2e-12
        total_tag+=tag/3
        total_imp+=imp/3
    return total_tag,total_imp


def integrated(p,upper,n_nodes):
    root,w=leggauss(n_nodes)
    start=np.zeros_like(upper)
    for x,weight in zip(root,w):
        _,imp=response(p,.5*(1+x)*upper)
        start+=.5*weight*upper*imp
    return start


def shell(level,r):
    pairs,w,neighbors,orient,e=prepared(level,r)
    D=RATE*w
    tests=[]
    for eta,p in zip((1,-1),pairs):
        # Mirror rows are antipodal faces ordered by the reversed boundary.
        # Apply equal face twists to corresponding outward-oriented faces.
        tag0,imp0=response(p,0.)
        idx=neighbors[e]
        sign=orient[e]
        results=[]
        for f in (.08,.04,.02):
            positive,jp=response(p[idx],sign*f)
            negative,jm=response(p[idx],-sign*f)
            odd=.5*D[idx]@(positive-negative)
            mean_imp=.5*(D[idx]*sign)@(jp+jm)
            baseline=D@tag0
            plus=(D[idx]@(positive-tag0[idx]))
            minus=(D[idx]@(negative-tag0[idx]))
            print(f"edge level={level} eta={eta:+} edge={e} adjacent={idx.tolist()} F={f:g} "
                  f"ordinary_tag_rate={baseline:.12g} extra_plus={plus:.12g} extra_minus={minus:.12g} "
                  f"odd_tag_rate={odd:.12g} even_generalized_force={mean_imp:.12g}",flush=True)
            results.append([odd,mean_imp])
        # The sum of the unwrapped curls is zero. For a sector endpoint
        # with principal angles q*2pi*w, any one face may take the return turn.
        per32=integrated(p,np.full_like(w,TAU),32)
        period=integrated(p,np.full_like(w,TAU),48)
        period_error=max(abs(per32-period))
        assert period_error<1e-10
        mu=D*period
        cycle=np.sum(mu[neighbors]*orient,axis=1)
        print(f"period level={level} eta={eta:+} per_face_min={min(period):.12g} "
              f"per_face_max={max(period):.12g} max_abs_shared_cycle_work={max(abs(cycle)):.12g} "
              f"cycle_test_edge={cycle[e]:.12g} quadrature_error={period_error:.3g}",flush=True)
        paths=[]
        for q in (1,-1):
            smooth=D@integrated(p,q*TAU*w,8)
            work=smooth-q*mu
            # Independent direct integration of the long return on one face.
            f0=idx[0]
            direct=D[f0]*(integrated(p[[f0]],np.array([q*TAU*(w[f0]-1)]),32)[0]-
                          integrated(p[[f0]],np.array([q*TAU*w[f0]]),8)[0])+smooth
            error=abs(direct-work[f0])
            print(f"insert level={level} eta={eta:+} q={q:+} smooth_part={smooth:.12g} "
                  f"min_work={min(work):.12g} max_work={max(work):.12g} "
                  f"return_check_error={error:.3g}",flush=True)
            assert error<1e-12
            paths.append(work)
        tests.append((np.array(results),np.array(paths),cycle))
    edge_err=max(abs(tests[0][0]+tests[1][0]).ravel())
    path_err=np.max(abs(tests[0][1]-tests[1][1][::-1]))
    cyc_err=max(abs(tests[0][2]+tests[1][2]))
    print(f"mirror_check level={level} held_odd_error={edge_err:.3g} insertion_error={path_err:.3g} "
          f"cycle_error={cyc_err:.3g} max_solid_angle={TAU*2*w.max():.9g}")
    assert max(edge_err,path_err,cyc_err)<2e-12


if __name__=="__main__":
    retained=calibration()
    for mesh in (2,3):
        shell(mesh,retained)
