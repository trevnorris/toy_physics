"""Bound the energy a standing-wave preparation can release into free rim links.

Four real mass-normalized patch displacements are represented by a complex
doublet. At fixed radius the positive additional stiffness consists of
squared neighbor differences and local alignment springs. The axial force
in the unbiased preparation fixes sum(I_b/omega_b) if the axial dependence
is solely the shared D/N term. Actions are held during the flux comparison.
"""
import numpy as np
from mpmath import iv
from scipy import sparse
from scipy.sparse.linalg import eigsh, spsolve

from mouth_cavity_occupation import S
from mouth_collar_link_bounds import R, PI, dot, rounded_bound
from mouth_collar_mobile_rim import prepare


def as_interval(x):
    return R(str(x))


def spectrum_bounds(p,edges,adj,tau,delta):
    # Obtain trial vectors in double precision; only outward interval bounds
    # on their Rayleigh and positive-vector quotients enter the energy test.
    deg = np.array(adj.sum(axis=1)).ravel().astype(int)
    lap = sparse.diags(deg)-adj
    tau0,delta0 = float(tau.mid),float(delta.mid)
    xyz = np.array([[float(y.mid) for y in z] for z in p])
    x,y,z=xyz.T
    scalar = tau0*lap+sparse.diags(delta0*(1-np.linalg.norm(xyz,axis=1)))
    nz=delta0*(x-1j*y)
    spin = sparse.bmat([[tau0*lap+sparse.diags(delta0*(1-z)),sparse.diags(-nz)],
                        [sparse.diags(-nz.conj()),tau0*lap+sparse.diags(delta0*(1+z))]],format="csr")
    seed = np.linspace(.8,1.1,len(p))
    bottom=eigsh(scalar,k=1,which="SA",v0=seed,tol=1e-12,return_eigenvectors=False)[0]
    pos=spsolve((scalar-(bottom-1.e-8)*sparse.eye(len(p))).tocsc(),np.ones(len(p)))
    _,states=eigsh(spin,k=1,which="SA",v0=np.r_[seed,.3*seed],tol=1e-12)
    up,down=np.split(states[:,0],2)
    v=[R(float(a)) for a in pos]
    assert min(v)>0
    waves=[tuple(R(float(u)) for u in row) for row in zip(up.real,up.imag,down.real,down.imag)]
    neigh=[R(0) for _ in p]
    num,den=R(0),R(0)
    for i,j in edges:
        neigh[i]+=v[j]
        neigh[j]+=v[i]
        d=[a-b for a,b in zip(waves[i],waves[j])]
        num+=tau*dot(d,d)
    lb=[]
    for i,(a,b,c,d) in enumerate(waves):
        prob=dot(waves[i],waves[i])
        pol=(2*(a*c+b*d),2*(a*d-b*c),a*a+b*b-c*c-d*d)
        num+=delta*(prob-dot(pol,p[i]))
        den+=prob
        lb.append(tau*(int(deg[i])-neigh[i]/v[i])+delta*(1-iv.sqrt(dot(p[i],p[i]))))
    return max(R(0).a,min(s.a for s in lb)), (num/den).b, max(deg)


def run():
    cfg=S["parameters"]
    ref,L=map(as_interval,S["phase_A"]["q_star"])
    c,mi,mo=map(as_interval,[cfg["c_w"],cfg["mu_in"],cfg["mu_out"]])
    slope_axial=(c*PI/2)**2/L**3
    for level in (2,3):
        # Each mesh is a possible fixed set of angular patches. At the reference
        # radius this uniform per-bond squared-frequency increment is prescribed.
        for eta in (1,-1):
            p,edges,adj,_=prepare(level,(".26","-.19",".17"),eta)
            for multiple in (".9","1","1.125","1.25","1.5","2"):
                scale=R(multiple)
                radius=ref*scale
                base=mi**2+(c*PI/radius)**2+(c*PI/(2*L))**2
                tau=R(".04")/scale**2
                delta=R(".16")
                low,upper,max_deg=spectrum_bounds(p,edges,adj,tau,delta)
                w_low, w_high=iv.sqrt(base+low),iv.sqrt(base+upper)
                assert w_high<mo
                # F_L = (c*pi/2)^2/L^3 * sum(I_b/omega_b). For the
                # relaxed lowest mode this also gives I/(omega_0) directly.
                per_load=(upper-low)*w_high/((w_high+w_low)*slope_axial)
                # Ordered eigenvalue comparison for any fixed modal actions
                # in the initial bound spectrum. Unoccupied branches are free.
                d2=min((2*int(max_deg)*tau).b,(mo**2-base).b)
                multi_per_load=d2*mo/((mo+iv.sqrt(base))*slope_axial)
                cost=PI/(2*16*scale)
                load=4*PI*as_interval(cfg["sigma"])*radius**2+4*PI*as_interval(cfg["alpha_volume"])*radius**3/3
                # Illustrative free-end support, no curvature term. A measured
                # wave-only load could replace this comparison input directly.
                print(f"budget level={level} eta={eta:+} R/R0={multiple} "
                      f"tau={float(tau.mid):.8g} delta={float(delta.mid):.8g} "
                      f"lowest_frequency_at_zero_upper={rounded_bound(w_high,8,True):.8f} "
                      f"lowest_frequency_any_field_lower={rounded_bound(w_low,8):.8f} "
                      f"axial_support_load={float(load.mid):.9g} "
                      f"exterior_cost_lower={rounded_bound(cost,8):.8f} "
                      f"ground_saving_upper={rounded_bound(load*per_load,8,True):.8f} "
                      f"any_fixed_occupation_saving_upper={rounded_bound(load*multi_per_load,8,True):.8f} "
                      f"ground_exclusion_below_load={rounded_bound(cost/per_load,6):.6f} "
                      f"any_fixed_occupation_exclusion_below_load={rounded_bound(cost/multi_per_load,6):.6f}",
                      flush=True)
                assert per_load>0 and multi_per_load>=per_load and cost>0


if __name__=="__main__":
    run()
