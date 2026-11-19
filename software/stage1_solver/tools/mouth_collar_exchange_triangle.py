"""A mobile internal doublet on one closed held triangle of prepared material.

Each loop is measured separately. No link angles or phases of incident paths
are registered. The net flux preferred by isolated triangles is a local
readout; neighboring triangles on a closed surface share their physical links.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.spatial.transform import Rotation

from mouth_collar_order_preparation import calibration, predicted
from mouth_collar_defect import states, face_angle
from mouth_material_links import section
from mouth_pair_constituent_response import PAULI


T = 1.
DIRECTIONS = np.arange(48)*2*np.pi/48
DELTA = 2*np.pi/len(DIRECTIONS)


def ring(p,j):
    H = np.zeros((6,6),complex)
    for i,v in enumerate(p):
        sl = slice(2*i,2*i+2)
        H[sl,sl] = -j*np.einsum("a,aij->ij",v,PAULI)
    for i in range(2):
        H[2*i:2*i+2,2*i+2:2*i+4] = -T*np.eye(2)
        H[2*i+2:2*i+4,2*i:2*i+2] = -T*np.eye(2)
    return H


def energy(H,f,occupancy):
    f = np.asarray(f)
    h = np.broadcast_to(H, f.shape+(6,6)).copy()
    hop = -T*np.exp(-1j*f)
    for s in (0,1):
        h[...,4+s,s] = hop
        h[...,s,4+s] = hop.conj()
    return np.linalg.eigvalsh(h)[...,:occupancy].sum(axis=-1)


def inspect(p,j,occupancy=1):
    h = ring(p,j)
    sampled = energy(h,DIRECTIONS,occupancy)
    width = np.ptp(sampled)
    if width < 1.e-11:
        return np.nan, width, 0.
    minima = np.flatnonzero((sampled<=np.roll(sampled,1)) & (sampled<=np.roll(sampled,-1)))
    minima_fit = []
    for idx in minima:
        fit = minimize_scalar(lambda u: energy(h,u,occupancy),
                              bounds=(DIRECTIONS[idx]-DELTA,DIRECTIONS[idx]+DELTA),
                              method="bounded",options={"xatol":1e-10})
        assert fit.success
        minima_fit.append((fit.fun,fit.x))
    value,phi = min(minima_fit)
    phase = np.angle(np.exp(1j*phi))
    eps = .001
    stiffness = (energy(h,phi+eps,occupancy)+energy(h,phi-eps,occupancy)-2*value)/eps**2
    return phase,width,stiffness


def held():
    for orientation in (1,-1):
        p = np.diag([1,orientation,1])
        target = face_angle(states(p)[None,:,:])[0]
        for j in (0.,1.,4.,16.,64.):
            for occupied in (1,2,6):
                f,width,stiffness = inspect(p,j,occupied)
                print(f"triangle orient={orientation:+} J={j:g} occupancy={occupied} "
                      f"preferred_F={f:.9g} projected_F={target:.9g} "
                      f"sampled_span={width:.9g} stiffness={stiffness:.9g}",flush=True)
                if occupied==6:
                    assert np.isnan(f)
    for scale in (1.,.5,0.):
        p = np.diag([1.,scale,1.])
        f,width,stiffness = inspect(p,4)
        print(f"middle_site J=4 polarization={scale:g} preferred_F={f:.9g} "
              f"sampled_span={width:.9g} stiffness={stiffness:.9g}",flush=True)


def collar(r):
    normals,faces = section(2)
    normals = Rotation.from_rotvec([.26,-.19,.17]).apply(normals)
    for eta in (1,-1):
        for k in (0,3,4,6):
            p = predicted(normals,eta,r,k)
            projected = sum(face_angle(states(p)[faces]))/(2*np.pi)
            for j in (1.,4.,16.):
                values = np.array([inspect(triangle,j) for triangle in p[faces]])
                print(f"prepared level=2 eta={eta:+} encounter={k} J={j:g} occupation=1 "
                      f"preferred_sum_turns={values[:,0].sum()/(2*np.pi):.9g} "
                      f"projected_ray_turns={projected:.9g} "
                      f"min_stiffness={min(values[:,2]):.9g} max_stiffness={max(values[:,2]):.9g}",
                      flush=True)
                assert np.isfinite(values).all() and min(values[:,2]) > 0


if __name__ == "__main__":
    held()
    collar(calibration())
