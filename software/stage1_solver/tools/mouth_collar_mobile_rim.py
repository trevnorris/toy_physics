"""One mobile doublet shared by all rim sites coupled to the same free links.

Certified upper bound on the energy it can gain from ANY link field. The
six-event material orientations are held. A scalar hopping problem gives a
lower bound in arbitrary flux by the norm inequality on every bond. The
zero-field spinor trial and scalar positive-vector lower bound are evaluated
with interval arithmetic; eigensolvers only supply candidate test vectors.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, spsolve
from mpmath import iv

from mouth_collar_link_bounds import (
    R, PI, TAU, dot, mesh, rotated, rounded_bound, spring
)


def prepare(level, rotation, eta):
    vertices,faces = mesh(level)
    axes = rotated(vertices,rotation)
    t6 = iv.cos(iv.sqrt(R(2))*R(".33"))**12
    h = iv.sqrt(t6)
    p = []
    for n in axes:
        x,y,z = (eta*c for c in n)
        b = (1-h)*(1+h-h*z)
        p.append((b*x,b*y,h+b*z))
    edges = sorted({tuple(sorted((i,j))) for f in faces for i,j in zip(f,f[1:]+f[:1])})
    i,j = np.array(edges).T
    nv = len(vertices)
    half = sparse.csr_matrix((np.ones(len(i)),(i,j)),shape=(nv,nv))
    return p, edges, half+half.T, len(faces)


def witnesses(p,adj,j):
    x,y,z = np.array([[float(c.mid) for c in v] for v in p]).T
    density = np.sqrt(x*x+y*y+z*z)
    scalar = -adj-sparse.diags(j*density)
    diag = sparse.diags
    spin = sparse.bmat([[-adj-diag(j*z), diag(-j*(x-1j*y))],
                        [diag(-j*(x+1j*y)), -adj+diag(j*z)]], format="csr")
    n = len(x)
    guess = np.linspace(.9,1.1,n)
    eig = eigsh(scalar,k=1,which="SA",v0=guess,tol=1e-12,return_eigenvectors=False)[0]
    # A resolvent slightly below the scalar bottom supplies a strictly positive
    # comparison vector even where the true scalar eigenmode is extremely small.
    positive = spsolve((scalar-(eig-1e-7)*sparse.eye(n)).tocsc(),np.ones(n))
    assert min(positive)>0
    _,trial = eigsh(spin,k=1,which="SA",v0=np.r_[guess,.3*guess],tol=1e-12)
    upper,lower = np.split(trial[:,0],2)
    return positive, list(zip(upper.real,upper.imag,lower.real,lower.imag))


def certify(level,rotation,eta,j,geometry):
    p,edges,adj,nf = geometry
    positive, trial = witnesses(p,adj,j)
    v = [R(float(a)) for a in positive]
    neighbors = [R(0) for _ in p]
    spinors = [tuple(R(float(a)) for a in row) for row in trial]
    energy, normalization = R(0),R(0)
    for i,l in edges:
        neighbors[i] += v[l]
        neighbors[l] += v[i]
        energy -= 2*dot(spinors[i], spinors[l])
    floors = []
    for pol,amp,neighbor,weight in zip(p,spinors,neighbors,v):
        floors.append(-j*iv.sqrt(dot(pol,pol))-neighbor/weight)
        a,b,c,d = amp
        sigma = (2*(a*c+b*d),2*(a*d-b*c),a*a+b*b-c*c-d*d)
        energy -= j*dot(pol,sigma)
        normalization += dot(amp,amp)
    scalar_lower = min(s.a for s in floors)
    zero_upper = energy/normalization
    saving = (zero_upper-scalar_lower).b
    field_floor = min(spring(PI,4).a,(5*spring(TAU/5,4)).a,
                      (PI/32+7*nf*spring(TAU/nf,4)).a)
    print(f"rim level={level} rotation={','.join(rotation)} eta={eta:+} J_per_t={j} "
          f"absolute_energy_lower={float(scalar_lower):.10f} "
          f"zero_trial_energy={float(zero_upper.mid):.10f} "
          f"possible_saving_per_t_upper={rounded_bound(saving,8,True):.8f} "
          f"nonneutral_or_tie_cost_per_K_lower={rounded_bound(field_floor,8):.8f} "
          f"neutral_for_t_per_K_below={rounded_bound(field_floor/saving,6):.6f}",flush=True)
    assert saving > 0


if __name__ == "__main__":
    for level in (2,3):
        for rotation in ((".26","-.19",".17"),("-.13",".24",".31")):
            for eta in (1,-1):
                geometry = prepare(level,rotation,eta)
                for j in (1,4,16):
                    certify(level,rotation,eta,j,geometry)
