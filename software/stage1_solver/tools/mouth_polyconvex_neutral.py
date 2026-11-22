"""Homogeneous, wave-empty 4D crystal specimens at prescribed order and finite strain."""
import hashlib
import json
from itertools import permutations, product
from math import factorial
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


MU = 1.0
K = 1.0
LAMBDA = K - MU/2
ROOT = Path(__file__).resolve().parents[3]


def energy(F, chi):
    jac = np.linalg.det(F)
    if np.any(jac <= 0):
        raise ValueError("Material orientation changed")
    return chi*(MU/2*(np.sum(F*F, axis=(-2, -1)) - 4) - MU*np.log(jac)
                + LAMBDA/2*(jac-1)**2)


def stress(F, chi):
    """First Piola stress work-conjugate to F in the reference four-volume."""
    jac = np.linalg.det(F)
    beta = -MU + LAMBDA*jac*(jac-1)
    return chi*(MU*F + beta[..., None, None]*np.linalg.inv(F).swapaxes(-1,-2))


def emit(**data):
    print(json.dumps(data, allow_nan=False), flush=True)


def reference_mesh():
    # Subdivide 2^4 congruent hypercubes, each into 4! equal-volume simplices.
    # All vertices other than the center lie on a prescribed outer face.
    cell_nodes, gradients = [], []
    for corner in product(range(2), repeat=4):
        for perm in permutations(range(4)):
            vertices = [np.array(corner)]
            for axis in perm:
                vertices.append(vertices[-1] + np.eye(4,dtype=int)[axis])
            v = np.array(vertices)/2
            D = v[1:]-v[0]
            grads = np.linalg.inv(D).T
            grads = np.vstack((-grads.sum(axis=0), grads))
            mobile = np.flatnonzero(np.all(v==.5, axis=1))
            cell_nodes.append(v)
            gradients.append(grads[mobile[0]] if mobile.size else np.zeros(4))
    return np.array(cell_nodes), np.array(gradients)


def patch(F, vertex, cells):
    """Deformation gradients reconstructed from positions of the five simplex vertices."""
    coords = cells @ F.T + np.all(cells==.5,axis=-1)[:,:,None]*vertex
    return np.linalg.solve(cells[:,1:]-cells[:,:1], coords[:,1:]-coords[:,:1]).transpose(0,2,1)


def patch_hold(F, chi, cells, g):
    # Releasing an interior node supplies a genuine non-affine lattice degree of freedom.
    v = np.array([.002, -.001, .0015, .00125])*min(np.linalg.svd(F, compute_uv=False))
    energies = []
    for iteration in range(20):
        Fc = patch(F,v,cells)
        P = stress(Fc,chi)
        grad = np.mean((P @ g[:,:,None])[:,:,0],axis=0)
        val = float(energy(Fc,chi).mean())
        energies.append(val)
        if np.linalg.norm(grad/chi) < 2.e-10:
            break
        J = np.linalg.det(Fc)
        normal = np.linalg.solve(Fc.transpose(0,2,1),g[:,:,None])[:,:,0]
        h = (chi*MU*np.mean(np.sum(g*g,axis=1))*np.eye(4)
             + chi*np.einsum('s,si,sj->ij', MU+LAMBDA*J**2,normal,normal)/len(g))
        v -= np.linalg.solve(h,grad)
    else:
        raise ArithmeticError(f"Interior lattice vertex did not equilibrate: {grad}")
    if max(np.diff(energies),default=0) > 5.e-12*max(1., abs(val)):
        raise ArithmeticError("Relaxation increased the energy")
    return float(np.linalg.norm(v)), float(np.linalg.norm(grad)/chi), energies[0]-energies[-1]


def hydro(J, chi, cells):
    ident = np.eye(4)
    def read(V):
        scale = V**.25
        Fc = patch(scale*ident,np.zeros(4),cells)
        traction = stress(Fc,chi) @ Fc.transpose(0,2,1) / V
        pressure = -np.trace(traction, axis1=-2,axis2=-1).mean()/4
        return pressure, energy(Fc,chi).mean()
    p, _ = read(J)
    d = J*.0001
    pm, em = read(J-d)
    pp, ep = read(J+d)
    bulk = -J*(pp-pm)/(2*d)
    p_energy = -(ep-em)/(2*d)
    theory = chi*(MU/J-MU/np.sqrt(J)-LAMBDA*(J-1))
    norm = chi*(MU/J + MU/np.sqrt(J)+LAMBDA*(J+1))
    expect = chi*(MU/J-MU/(2*np.sqrt(J))+LAMBDA*J)
    return dict(volume_ratio=J, shear_order=chi, pressure_from_traction=p,
                pressure_from_work=p_energy, tangent_from_traction=bulk,
                predicted_tangent=expect,
                relative_pressure_error=max(abs(p-theory),abs(p-p_energy))/norm,
                relative_tangent_error=abs(bulk/expect-1), old_log_law_tangent=chi*K*(1-np.log(J))/J)


def acoustic(F, chi, cells):
    # Physical-plane diagonal and oblique pointwise probes: use centered traction
    # differences and an absolute stress derivative independent of the wave solver.
    n = np.array([1.,-2.,.5,.2])
    n /= np.linalg.norm(n)
    m = np.array([.5,1.,.6,-.3])
    m -= n.dot(m)*n
    m /= np.linalg.norm(m)
    mode = np.outer(m,n)@F
    J = np.linalg.det(F)
    transverse = chi*MU*np.linalg.norm(F.T@n)**2
    errors = []
    for polarization in (m,n,(m+n)/np.sqrt(2.)):
        H = np.outer(polarization,n)@F
        step = 1e-5
        pp = stress(patch(F+step*H,np.zeros(4),cells),chi).mean(axis=0)
        pm = stress(patch(F-step*H,np.zeros(4),cells),chi).mean(axis=0)
        measured = np.sum(H*(pp-pm))/(2*step)
        target = transverse + chi*(MU+LAMBDA*J**2)*(polarization@n)**2
        errors.append(abs(measured/target-1.))
    return dict(transverse_speed_squared=transverse/chi,
                max_rank_one_relative_error=max(errors))


def archive_envelope():
    folder = ROOT/"_scratch/s11c/finite_bore_graph01/fields"
    hashes = json.loads((ROOT/"research/pde_ledger_v3/_measurements/mouth_finite_bore_support.json"
                          ).read_bytes())["central_fields"]
    for radius, depth, side in product((1,2),("long","unit"),("ordered","receiver")):
        name = f"r{radius}_{depth}_g2_{side}.npz"
        blob = folder/name
        if hashlib.sha256(blob.read_bytes()).hexdigest()!=hashes[name]:
            raise ValueError("Archived deformation changed")
        with np.load(blob) as d:
            rs,rz,ws,wz,h = d["gradient"].T
            chi,J = d["chi"],d["J"]
            mer = np.stack((rs,rz,ws,wz),axis=-1).reshape(-1,2,2)
            slowest = np.minimum(np.linalg.svd(mer,compute_uv=False)[:,-1],h)
            eigen = chi*MU*slowest**2/J
            pure = chi*(MU/J-MU/(2*np.sqrt(J))+LAMBDA*J)
            if min(eigen)<=0 or min(pure)<=0:
                raise ArithmeticError("Vanishing stiffness in a finite positive-order specimen")
            emit(case="archived_strain_envelope",radius=radius,depth=depth,side=side,
                 points=len(J), minimum_current_acoustic_eigenvalue=min(eigen),
                 min_volume_ratio=min(J),max_volume_ratio=max(J),min_order=min(chi),max_order=max(chi))


def main():
    cells,g=reference_mesh()
    if cells.shape!=(2**4*factorial(4),5,4) or np.max(abs(g.mean(axis=0)))>1e-12:
        raise ArithmeticError("4D patch closure failed")
    for chi in (1.,.01,.0001):
        for J in (.0001,.01,.25,1.,4.,300.):
            row=hydro(J,chi,cells)
            emit(case="volumetric_hold", **row)
            if (row["relative_pressure_error"]>1e-7 or row["relative_tangent_error"]>1e-7
                    or row["tangent_from_traction"]<=0):
                raise ArithmeticError("Volume response mismatch")
        for name, F in (("oblique_compression",np.array([[.23,.31,0.,0.],[0.,1.2,.17,0.],
                                                          [0.,0.,.62,.09],[0.,0.,0.,.8]])),
                        ("mixed_extension",np.array([[2.7,1.,0.,.4],[0.,.3,.12,0.],
                                                      [0.,0.,1.7,.3],[0.,0.,0.,1.1]]))):
            v,grad,release=patch_hold(F,chi,cells,g)
            row=acoustic(F,chi,cells)
            emit(case="free_interior_node",deformation=name,chi=chi,volume_ratio=np.linalg.det(F),
                 vertex_residual=v, scaled_force_residual=grad, energy_released=release, **row)
            if v>1e-9 or grad>2e-10 or release<=0 or row["max_rank_one_relative_error"]>1e-7:
                raise ArithmeticError("Mixed finite-strain or acoustic control failed")
    # Axial engineering strain, one normal traction free, other two strains fixed.
    forces = []
    for e in (-.001,.001):
        def P(norm):
            return stress(np.diag([1.+e,1.,1.,norm]),1.)
        free = brentq(lambda w:P(w)[3,3],.8,1.2)
        forces.append(P(free)[0,0])
    axial=(forces[1]-forces[0])/.002
    emit(case="strip_and_bulk", strip_tangent=axial,zero_order_lattice_energy=
         float(energy(np.diag([1.1,.9,1.3,1.2]),0.)))
    if abs(axial-2.4)>1e-5:
        raise ArithmeticError("Neutral strip tangent changed")
    archive_envelope()


if __name__=="__main__":
    main()
