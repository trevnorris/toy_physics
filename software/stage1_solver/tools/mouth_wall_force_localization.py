"""Localize the direct-pressure minus interior-momentum residual on saved fields.

Every row pairs the same smooth wall test displacement with its extension into
the fluid. The cutoff is one at the complete wall and zero at the reservoirs.
Radial distance to the exit edge is measured in the (r,w) meridian, and thus
the exit neighborhoods include the whole spherical junction. No field solve.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_inertia import WINDOWS, taper
from mouth_carrier_lateral_response import prepare, trace
from mouth_flared_receiver import EXIT, EXIT_RADIUS
from mouth_occupied_push_baths import load_bath


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def windows(r, w):
    a, b, c, d = WINDOWS[0]
    cr, dr = taper(r, a, b)
    cw, dw = taper(w, c, d)
    whole, wr, ww = cr*cw, dr*cw, cr*dw
    yield 'whole', whole, wr, ww
    s = np.hypot(r-EXIT_RADIUS, w-EXIT)
    for inner, outer in ((.08, .20), (.20, .40), (.4, .70), (.70, 1.)):
        g, dg = taper(s, inner, outer)
        u = dg / np.maximum(s, 1e-30)
        yield f'exit_{inner:g}_{outer:g}', g*whole, g*wr + whole*u*(r-EXIT_RADIUS), g*ww + whole*u*(w-EXIT)
    for lo, hi in ((.35,.65), (1.1,1.4), (2.3,2.6)):
        g, dg = taper(w, lo, hi)
        yield f'upstream_{lo:g}_{hi:g}', g*whole, g*wr, g*ww + whole*dg


def phasor(re, im):
    return re + 1j*im


def fields(r, f, fr, fw, omega):
    base, hc, hq = f.T
    br, rc, rq = fr.T
    bw, wc, wq = fw.T
    def coefficients(h, dh_r, dh_w):
        grad_x = dh_r + 2*h/r
        e = (br.conj()*dh_r + bw.conj()*dh_w + .5*(abs(base)**8-1)*base.conj()*h).real
        n_phase = (base.conj()*h).imag
        jx = (base.conj()*grad_x + h.conj()*br).imag
        pr = (br.conj()*dh_r + grad_x.conj()*br).real - e
        pw = (br.conj()*dh_w + grad_x.conj()*bw).real
        return np.array([e,n_phase,jx,pr,pw])
    energy, time, jx, pr, pw = phasor(coefficients(hc,rc,wc),coefficients(hq,rq,wq))
    # On the wall, relative carriage cancels the Neumann traction, leaving
    # pressure. In the volume tensor it cancels the translating-chart boost.
    actuator_pressure = energy + 1j*omega*time - (base.conj()*br).imag
    return actuator_pressure, 1j*omega*jx, -pr + 1j*omega*time, -pw


def measure(dom, nodes, omega, order):
    results = {}
    root, wt = leggauss(order)
    flat = lambda a: a.reshape((-1, a.shape[-1]))
    for block, (pg, vg, dofs) in dom.blocks.items():
        if block in (0, 2):
            half = np.diff(pg)[:,None]/2
            p = pg[:-1,None]+half*(1+root)
            r,w,_,_,wp,_,f,fr,fw = trace(dom,block,p,np.zeros_like(p),nodes)
            weight = (4*np.pi/3*r**2*wp*half*wt).ravel()
            pressure, *_ = fields(r.ravel(),flat(f),flat(fr),flat(fw),omega)
            for name, g, _, _ in windows(r.ravel(), w.ravel()):
                results.setdefault(name, np.zeros(2,dtype=complex))[0] += weight @ (g*pressure)
        u,v = np.meshgrid(root,root,indexing='ij')
        tensor_weights = np.outer(wt,wt).ravel()
        for start in range(0,len(dofs),128):
            ids = np.arange(start,min(start+128,len(dofs)))
            i,j = np.unravel_index(ids,(len(pg)-1,len(vg)-1))
            hp,hv = np.diff(pg)[i,None]/2,np.diff(vg)[j,None]/2
            p,q = pg[i,None]+hp*(1+u.ravel()),vg[j,None]+hv*(1+v.ravel())
            r,w,rp,rv,wp,wv,f,fr,fw = trace(dom,block,p,q,nodes)
            weight = (4*np.pi/3*r**2*abs(rp*wv-rv*wp)*hp*hv*tensor_weights).ravel()
            _, storage, pr, pw = fields(r.ravel(),flat(f),flat(fr),flat(fw),omega)
            for name, g, gr, gw in windows(r.ravel(),w.ravel()):
                results.setdefault(name, np.zeros(2,dtype=complex))[1] += weight @ (g*storage+gr*pr+gw*pw)
    return results


def normal_residuals(dom, nodes, order):
    """One-sided Neumann traces; omit the edge of zero face measure itself.

    At unit real translational speed the material wall obeys partial_n psi =
    i V n_x psi. The first angular harmonic includes the Y=x/r factor in
    both the prescribed derivative and the observed response.
    """
    root, weight = leggauss(order)
    for name,block,along,fixed,outward in (
            ('conduit',2,'p',0.,(1.,0.)), ('receiver_plane',4,'v',EXIT,(0.,-1.))):
        pg,vg,_ = dom.blocks[block]
        edges = pg if along=='p' else vg
        half = np.diff(edges)/2
        q = (edges[:-1,None]+half[:,None]*(1+root)).ravel()
        other = np.full_like(q,fixed)
        p,v = (q,other) if along=='p' else (other,q)
        r,w,rp,rv,wp,wv,f,fr,fw = trace(dom,block,p,v,nodes)
        tr,tw = (rp,wp) if along=='p' else (rv,wv)
        length = np.hypot(tr,tw)
        orient = np.sign(tw*outward[0]-tr*outward[1])
        nr,nw = orient*tw/length,-orient*tr/length
        normal = nr[:,None]*fr+nw[:,None]*fw
        normal[:,1] -= 1j*nr*f[:,0]
        area = 4*np.pi*r*r*length*(half[:,None]*weight).ravel()
        distance = np.hypot(r-EXIT_RADIUS,w-EXIT)
        for lo,hi in ((0.,.1),(.1,.25),(.25,.5),(.5,1.)):
            mask = (distance>=lo)&(distance<hi)
            if not np.any(mask):
                continue
            rms = np.sqrt(area[mask] @ (abs(normal[mask])**2) / area[mask].sum())
            print(json.dumps(dict(order=order,test='boundary_normal',face=name,
                                  distance=[lo,hi],quadrature_samples=int(mask.sum()),
                                  area=float(area[mask].sum()),
                                  base_rms=rms[0],cos_dipole_rms=rms[1],
                                  quadrature_dipole_rms=rms[2])),flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', required=True)
    ap.add_argument('--orders', type=int,nargs='+',default=(4,8))
    ap.add_argument('--normal-only',action='store_true')
    args=ap.parse_args()
    arc=Path('research/pde_ledger_v3/_measurements/mouth_rigid_drag_table.json')
    case=json.loads(arc.read_text())['runs'][args.run]
    cmd=case['command']
    field_path=Path(cmd[cmd.index('--save-response')+1])
    source_flag='--source' if '--source' in cmd else '--bath'
    folder=Path(cmd[cmd.index(source_flag)+1])
    if digest(field_path)!=case['response_field_sha256'] or digest(folder/'field.npz')!=case['source_field_sha256']:
        raise ValueError('Saved field differs from the original measurement')
    if source_flag=='--source':
        dom, base=prepare(argparse.Namespace(source=folder))
    else:
        dom, base=load_bath(folder)
    with np.load(field_path) as saved:
        np.testing.assert_array_equal(saved['points'],dom.points)
        omega=float(saved['omega'])
        hr,hi=saved['harmonic_real'],saved['harmonic_imag']
        nodes=np.column_stack([base, hr.real+1j*hi.real, hr.imag+1j*hi.imag])
        reference=complex(saved['impedance'])
    print(json.dumps(dict(run=args.run,field_sha=case['response_field_sha256'],omega=omega,
                          wall_legacy=[reference.real,reference.imag])),flush=True)
    for order in args.orders:
        if args.normal_only:
            normal_residuals(dom,nodes,order)
            continue
        results=measure(dom,nodes,omega,order)
        for name,(wall,volume) in results.items():
            print(json.dumps(dict(order=order,test=name,wall=[wall.real,wall.imag],
                                  momentum=[volume.real,volume.imag],
                                  residual=[(wall-volume).real,(wall-volume).imag])),flush=True)
        if order>=8:
            np.testing.assert_allclose(results['whole'][0],reference,atol=3e-10,rtol=3e-10)


if __name__=='__main__':
    main()
