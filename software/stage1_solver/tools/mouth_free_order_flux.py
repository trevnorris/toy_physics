"""Independent flux integral across a saved diffuse half-order surface in the conduit."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_curved_phase import quadratic, tensor
from mouth_flared_receiver import EXIT, ROUNDING, FlaredDomain
from mouth_free_order_diffuse import emit


def evaluate(dom, field, eta, z):
    """Read the Q2 field and derivatives in (eta=r/R(z), z) at paired points."""
    pg, vg, dofs = dom.blocks[2]
    v = 1-eta
    i = np.clip(np.searchsorted(pg, z)-1, 0, len(pg)-2)
    j = np.clip(np.searchsorted(vg, v)-1, 0, len(vg)-2)
    hp, hv = np.diff(pg)[i]/2, np.diff(vg)[j]/2
    ap, dp = quadratic((z-pg[i])/hp-1)
    av, dv = quadratic((v-vg[j])/hv-1)
    values = field[dofs[i*(len(vg)-1)+j]]
    def product(basis):
        return np.einsum('ij,ij->i', basis, values)
    return product(tensor(ap,av)), -product(tensor(ap,dv))/hv, product(tensor(dp,av))/hp


def crossing(dom, chi, carrier, subdivisions):
    _, vg, _ = dom.blocks[2]
    edges = 1-vg[::-1]
    nodes, weights = leggauss(6)
    panel = (np.arange(subdivisions)[:,None]+(1+nodes)/2).ravel()/subdivisions
    eta = (edges[:-1,None] + np.diff(edges)[:,None]*panel).ravel()
    weights = (np.diff(edges)[:,None]*np.tile(weights,subdivisions)[None,:]/(2*subdivisions)).ravel()
    lo, hi = np.full_like(eta,ROUNDING), np.full_like(eta,EXIT)
    if np.any(evaluate(dom,chi,eta,lo)[0] <= .5) or np.any(evaluate(dom,chi,eta,hi)[0] >= .5):
        raise ValueError('midpoint surface leaves the conduit bracket')
    for _ in range(48):
        mid = .5*(lo+hi)
        above = evaluate(dom,chi,eta,mid)[0] > .5
        lo, hi = np.where(above,mid,lo), np.where(above,hi,mid)
    z = (lo+hi)/2
    level, c_eta, c_z = evaluate(dom,chi,eta,z)
    if np.any(c_z >= 0):
        raise ValueError('nonregular crossing')
    dz = -c_eta/c_z
    R, Rp = dom.radius(z)
    r, dr = eta*R, R+eta*Rp*dz
    f, f_eta, f_z = evaluate(dom,carrier,eta,z)
    jr = (f.conj()*f_eta/R).imag
    jz = (f.conj()*f_z).imag - eta*Rp*jr
    measure = jz*dr-jr*dz
    return dict(subdivisions=subdivisions, points=len(eta), integrated_j=weights @ (4*np.pi*r*r*measure),
                min_j_normal=np.min(measure/np.hypot(dr,dz)), minimum_radial_derivative=dr.min(),
                level_error=np.max(abs(level-.5)))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('orders', nargs='+', type=Path)
    args = parser.parse_args()
    import json
    state = json.loads((args.source/'stdout').read_text().splitlines()[0])
    dom = FlaredDomain(state['n'],state['profile'],state['roof'],state['feed_far'])
    with np.load(args.source/'field.npz') as saved:
        np.testing.assert_array_equal(dom.points,saved['points'])
        psi = saved['psi']
    for path in args.orders:
        with np.load(path) as saved:
            chi = saved['chi']
        emit(order=str(path), order_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
             source=str(args.source), target=state['bulk'],
             results=[crossing(dom,chi,psi,n) for n in (1,4,12)])
