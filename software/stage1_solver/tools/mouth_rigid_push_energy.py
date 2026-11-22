"""Cycle-mean excess Hamiltonian energy crossing the translating rigid carrier.

The nondimensional Hamiltonian is e=|grad psi|²/2+(n⁵-5n+4)/20. Its lab
current is -Re[(dt_lab psi)* grad psi]. Report powers per |V|² for unit real
center velocity and e^(+i omega t). The independent control volumes use no
wall-pressure trace. The driving force fixes the complex velocity afterwards.
"""
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_carrier_free_center import momentum_force_rows, motion_source, wall_row
from mouth_carrier_lateral_inertia import WINDOWS, taper
from mouth_carrier_lateral_response import Dipole, emit, parser, prepare, trace, weighted_gram
from mouth_flared_receiver import EXIT
from mouth_occupied_push_baths import digest, load_bath


def unit_velocity_response(dom, f, omega):
    """Linear drive in resting baths in the actual benchmark's translating box."""
    model = Dipole(dom, f)
    mass = weighted_gram(dom.N, dom.weight)
    skew = sparse.bmat([[None, -mass], [mass, None]], format='csr')
    dynamic = model.operator.astype(complex) - 1j * omega * skew
    source = motion_source(dom, f)
    free = model.free
    operator = dynamic[free][:, free].tocsc()
    factor = splu(operator)
    h = np.zeros(2 * len(f), complex)
    h[free] = factor.solve(source[free].astype(complex))
    residual = max(abs(operator @ h[free] - source[free]))
    static, temporal = momentum_force_rows(dom, f, source)
    impedance = (static + 1j * omega * temporal) @ h
    trace_static, trace_temporal, trace_boost = wall_row(dom, f)
    direct = -((trace_static + 1j * omega * trace_temporal) @ h + trace_boost)
    emit(test='unit_velocity', nodes=len(f), omega=omega, dynamic_residual=residual,
         actuator_impedance=[impedance.real, impedance.imag],
         force_method='volume_momentum_and_remote_stress',
         direct_trace_impedance=[direct.real, direct.imag],
         actuator_mean_power_per_speed_squared=.5 * impedance.real)
    if residual > 2e-8 or abs(impedance) == 0:
        raise ArithmeticError('Harmonic carrier solve or applied impedance unresolved')
    return h, impedance


def angular_current(r, f, fr, fw, omega):
    """Radial and w currents, scaled by 3 after the brane solid-angle average.

    Columns store stationary carrier, first-harmonic cosine, sine (each complex
    in carrier space). In the moving chart, dt_lab = dt_chart - V(t) d_x.
    At order V² the unknown second-order field enters J only via its chart time
    derivative times grad f; bounded periodic motion makes this term average
    to zero. The swept term -V e e_x uses the linear energy density.
    """
    base, hc, hs = np.moveaxis(f, -1, 0)
    base_r, hc_r, hs_r = np.moveaxis(fr, -1, 0)
    base_w, hc_w, hs_w = np.moveaxis(fw, -1, 0)
    de_c = (base_r.conj()*hc_r + base_w.conj()*hc_w +
            (.5 * (abs(base)**8-1)) * base.conj()*hc).real
    # Solid-angle averages of 3 d_x[h Y] and 3 Y² are (h_r+2h/r) and 1.
    dx_c = hc_r + 2*hc/r
    def direction(dc, ds, grad_base, sweep):
        time = -omega*(hs.conj()*dc - hc.conj()*ds).real
        convected_base = (base_r.conj()*dc).real
        convected_dipole = (dx_c.conj()*grad_base).real
        return .5*np.stack((time, convected_base, convected_dipole, sweep), axis=-1)
    return (direction(hc_r, hs_r, base_r, -de_c),
            direction(hc_w, hs_w, base_w, np.zeros_like(de_c)))


def section(dom, fields, omega, block, along, value, normal, clip=None):
    pg, vg, _ = dom.blocks[block]
    edges = pg if along == 'p' else vg
    if clip:
        lo, hi = clip
        edges = np.r_[lo, edges[(edges > lo) & (edges < hi)], hi]
    root, weights = leggauss(8)
    half = np.diff(edges) / 2
    q = (edges[:-1,None] + half[:,None]*(1+root)).ravel()
    dq = (half[:,None]*weights).ravel()
    fixed = np.full_like(q, value)
    p, v = (q, fixed) if along == 'p' else (fixed, q)
    r, _, rp, rv, wp, wv, f, fr, fw = trace(dom, block, p, v, fields)
    tr, tw = (rp, wp) if along == 'p' else (rv, wv)
    orient = tw*normal[0] - tr*normal[1]
    if np.any(abs(orient)<1e-12):
        raise ValueError('Ambiguous boundary normal')
    nr, nw = np.sign(orient)*tw, -np.sign(orient)*tr
    jr, jw = angular_current(r, f, fr, fw, omega)
    weight = 4*np.pi/3*r*r*dq
    return weight @ (jr*nr[:,None] + jw*nw[:,None])


def window_currents(dom, fields, omega):
    """Integrate the physical flux on a continuum of cuts by -grad(chi).Jrel."""
    root, weights = leggauss(6)
    u, v = [x.ravel() for x in np.meshgrid(root, root, indexing='ij')]
    weights = np.outer(weights, weights).ravel()
    results = np.zeros((len(WINDOWS), 2, 4))
    for block, (pg, vg, dofs) in dom.blocks.items():
        for start in range(0, len(dofs), 128):
            cells = np.arange(start, min(len(dofs), start+128))
            i, j = np.unravel_index(cells, (len(pg)-1, len(vg)-1))
            hp, hv = np.diff(pg)[i,None]/2, np.diff(vg)[j,None]/2
            p, s = pg[i,None]+hp*(1+u), vg[j,None]+hv*(1+v)
            r, w, rp, rv, wp, wv, f, fr, fw = trace(dom, block, p, s, fields)
            dv = 4*np.pi/3*r*r*hp*hv*abs(rp*wv-rv*wp)*weights
            jr, jw = angular_current(r, f, fr, fw, omega)
            for k, window in enumerate(WINDOWS):
                radial, dr = taper(r, *window[:2])
                axial, dw = taper(w, *window[2:])
                integrand = jr*(dr*axial)[...,None] + jw*(radial*dw)[...,None]
                # The two windows transition wholly outside the curved wall.
                results[k, int(block >= 3)] -= (dv[...,None]*integrand).sum(axis=(0,1))
    return results


def partition(rows):
    return {name:dict(terms=terms.tolist(), total=float(sum(terms)))
            for name,terms in rows.items()}


def run(args):
    if args.bath and args.source:
        raise ValueError('Choose a saved original or nested bath')
    if args.omega <= 0 or args.amplitude <= 0:
        raise ValueError('Expected positive frequency and reference push')
    if args.bath:
        dom, field = load_bath(args.bath)
        emit(test='nested_preparation', source=str(args.bath/'field.npz'), far=dom.far,
             roof=dom.roof, sha256=digest(args.bath/'field.npz'))
    else:
        dom, field = prepare(args)
    h, impedance = unit_velocity_response(dom, field, args.omega)
    hre, him = np.split(h, 2)
    fields = np.column_stack([field, hre.real+1j*him.real, -(hre.imag+1j*him.imag)])
    if args.save_response:
        source = args.source or args.bath
        with args.save_response.open('xb') as output:
            np.savez_compressed(output, points=dom.points, harmonic_real=hre,
                                harmonic_imag=him, omega=args.omega,
                                background_sha256=digest(source/'field.npz') if source else '',
                                velocity_per_force=1/impedance, impedance=impedance)
    def face(block, along, value, normal, clip=None):
        return section(dom, fields, args.omega, block, along, value, normal, clip)
    outer = dict(feed=face(1,'v',dom.far,(1,0)), bulk_side=face(4,'p',1.,(1,0)),
                 roof_inner=face(3,'v',dom.roof,(0,1)), roof_outer=face(4,'v',dom.roof,(0,1)))
    reflecting = dict(lip=face(0,'p',0.,(1,1)), neck_flare=face(2,'p',0.,(1,0)),
                      feed_plane=face(1,'p',0.,(0,1)), feed_floor=face(1,'p',1.,(0,-1)),
                      inner_floor=face(0,'p',1.,(0,-1),(-np.pi/2, dom.bottom_turn)),
                      bulk_plane=face(4,'v',EXIT,(0,-1)))
    windows = window_currents(dom, fields, args.omega)
    boundary = sum(v.sum() for v in outer.values())
    wall_trace = sum(v.sum() for v in reflecting.values())
    wall = .5*impedance.real
    emit(test='boundary_outward_power_per_speed_squared',
         terms=['field_time','advection_of_base','advection_of_dipole','moving_cut'],
         remote=partition(outer), reflecting_raw=partition(reflecting),
         remote_total=boundary, reflecting_raw_total=wall_trace,
         all_raw_boundary_residual=boundary+wall_trace, actuator_minus_remote=wall-boundary)
    for i, parts in enumerate(windows):
        emit(test='enclosing_outward_power_per_speed_squared', window=WINDOWS[i],
             parts=partition(dict(feed=parts[0], receiver=parts[1])),
             total=parts.sum(), actuator_minus_outward=wall-parts.sum(),
             outward_minus_remote=parts.sum()-boundary)
    for relative in [0.,1.,-1.,.5,-.5]:
        push=args.amplitude*relative
        velocity=push/impedance
        speed2=abs(velocity)**2
        emit(test='linear_response_quadratic_power', applied_force=push,
             velocity=[velocity.real,velocity.imag], actuator=.5*push*velocity.real,
             remote=speed2*boundary, windows=(speed2*windows.sum(axis=(1,2))).tolist())
    if not np.all(np.isfinite(windows)) or not np.isfinite(wall_trace+boundary):
        raise ArithmeticError('Nonfinite energy readout')


if __name__=='__main__':
    p=parser()
    p.add_argument('--bath',type=Path)
    p.add_argument('--omega',type=float,required=True)
    p.add_argument('--amplitude',type=float,default=1e-5)
    p.add_argument('--save-response',type=Path)
    run(p.parse_args())
