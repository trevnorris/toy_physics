"""Compare finite polar material energies to direct four-dimensional map composition."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_bonded_material_motion import Bonded
from mouth_internal_collar_check import read
from mouth_live_waist_order import Order
from mouth_polar_bonded import PolarBonded, StaticPolar, order_at, geometry, outer_crystal_source_delta
from mouth_sealed_shape_relax import DELTA
from mouth_supercore_front import check_guard


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def energy4(material, epsilon, quad):
    wt, j, K, hoop, d, s, v = material
    # Form full spatial gradients, then compute a volume ratio and the
    # squared Frobenius norm independently of the invariant expansion.
    F = np.zeros((len(j), 4, 4))
    F[:, :2, :2] = np.linalg.cholesky(K)
    F[:, 2, 2] = F[:, 3, 3] = hoop
    result = 0.
    for p, pa, w in zip(*quad):
        A = np.broadcast_to(np.eye(4), F.shape).copy()
        A[:, 0, 0] += epsilon*d*p
        A[:, 0, 1] += epsilon*v*p
        A[:, 0, 2] += epsilon*s*pa
        A[:, 2, 2] += epsilon*s*p
        A[:, 3, 3] += epsilon*s*p
        G = A @ F
        jac = np.linalg.det(G)
        W = (np.sum(G**2, axis=(1, 2))-4)/2-np.log(jac)+.25*(jac-1)**2
        result += w*(wt@W)
    return float(result)


def direct_order(data, epsilon, quad):
    wt, gr, gz, potential, d, s, v = data
    grad = np.column_stack((gr, gz, np.zeros_like(gr), np.zeros_like(gr)))
    res = 0.
    for p, pa, w in zip(*quad):
        A = np.broadcast_to(np.eye(4), (len(gr), 4, 4)).copy()
        A[:, 0, 0] += epsilon*d*p
        A[:, 0, 1] += epsilon*v*p
        A[:, 0, 2] += epsilon*s*pa
        A[:, 2, 2] += epsilon*s*p
        A[:, 3, 3] += epsilon*s*p
        pulled = np.linalg.solve(A.transpose(0, 2, 1), grad[..., None])[..., 0]
        res += w*(wt @ (np.linalg.det(A)*(DELTA**2*np.sum(pulled**2, axis=1)+potential)))
    return res


def run(args):
    check_guard(args.state)
    sh, _, inner, _ = read(args.state)
    b = PolarBonded(args.state, 5, 12)
    bare = Bonded(args.state, 5, 12, full_depth=True)
    stat = StaticPolar(args.state)
    rng = np.random.default_rng(58347)
    q = .002*rng.normal(size=b.mask.sum())
    direction = rng.normal(size=len(q))
    direction /= np.linalg.norm(direction)
    u = b.free.copy()
    u[b.mask] = q/b.scale
    nodes = u.reshape(-1, 2)[b.id]
    val = np.einsum('qi,qik->qk', b.basis, nodes)
    r, z = b.r+val[:, 0], b.w+val[:, 1]
    B = np.eye(2) + np.stack([np.einsum('qi,qik->qk', g, nodes) for g in (b.dr, b.dw)], axis=-1)
    chi = order_at(b, r, z)[0]
    body_sites = (b.measure*chi, b.j*np.linalg.det(B)*(r/b.r)**2,
                  B @ b.K @ B.transpose(0, 2, 1), b.h*r/b.r, *geometry(b.chart.shape, r, z, False))
    source_sha = hashlib.sha256(args.state.read_bytes()).hexdigest()
    val, grad = b.energy(q)
    en0, g0 = bare.energy(q)
    st0 = stat.terms(0)
    ine0 = direct_order(stat.inner, 0, b.quad)
    native = Order(sh.shape, len(inner[1])-1).energy_gradient(inner[-1])[0]
    emit(kind='baseline', source_sha256=source_sha, body=val, spherical=en0,
         gradient_error=np.max(abs(grad-g0)), static_delta=st0,
         direct_inner_order=ine0, native_inner_order=native)
    if val != en0 or np.any(grad != g0) or any(st0.values()) or abs(native-ine0)>1e-8:
        raise ArithmeticError('Spherical limit')
    refout = direct_order(stat.outer, 0, b.quad)
    ext0 = energy4(stat.elastic, 0, b.quad)
    for e in (-.025, .025):
        b.epsilon = e
        val, grad = b.energy(q)
        true_body = energy4(body_sites, e, b.quad)
        direct_ext = energy4(stat.elastic, e, b.quad)-ext0
        terms = stat.terms(e)
        shift, dw = stat.tail
        P, _, w = b.quad
        tail = sum(wi*(dw @ ((12.**3-(12.+e*shift*pi)**3)/3)) for pi, wi in zip(P, w))
        direct_in = direct_order(stat.inner, e, b.quad)-ine0
        direct_out = direct_order(stat.outer, e, b.quad)-refout+tail
        metric = b.energy(q, diagnose=True)
        emit(kind='composition', epsilon=e, body=val, direct_body=true_body, terms=terms,
             direct_outer_crystal=direct_ext, direct_inner_order=direct_in,
             direct_outer_order=direct_out, flat_continuation=tail,
             min_J=metric['min_J'], core_delta=metric['core_polar_change'],
             inner_collar_delta=metric['inner_collar_polar_change'])
        if max(abs(val-true_body), abs(direct_in-terms['inner_order']),
               abs(direct_out-terms['outer_order']), abs(direct_ext-terms['outer_crystal']))>1e-8:
            raise ArithmeticError('Full map and invariant contraction disagree')
        for h in (1e-5, 5e-6):
            fd = (b.energy(q+h*direction)[0]-b.energy(q-h*direction)[0])/(2*h)
            emit(kind='material_derivative', epsilon=e, step=h, finite_difference=fd,
                 gradient=float(grad@direction))
            if abs(fd-grad@direction)>5e-7:
                raise ArithmeticError('Material gradient after finite polar displacement')
    for h in (.002, .001):
        signs = []
        for sign in (-1, 1):
            b.epsilon = sign*h
            signs.append(dict(body=b.energy(q)[0]-en0, **stat.terms(sign*h)))
        slope = {k:(signs[1][k]-signs[0][k])/(2*h) for k in signs[0]}
        stiffness = {k:(signs[1][k]+signs[0][k])/h**2 for k in signs[0]}
        emit(kind='shape_stiffness', step=h, minus=signs[0], plus=signs[1],
             derivative=slope, curvature=stiffness, total=sum(stiffness.values()))


def check_receiver(args):
    """Envelope theorem on the constrained, re-solved exterior order field."""
    from mouth_receiver_order_equilibrium import receiver
    from mouth_poly_shallow_release import CUT

    check_guard(args.state)
    pair, labels, inner, outer = read(args.state)
    sh = pair.shape
    mesh = len(inner[1])-1
    base = pair.evaluate(labels, inner, outer)['outer_source']

    def solve(e):
        profile, info = receiver(sh, mesh, CUT, lambda wa, x: base,
                                 polar_shift=e, material_path=args.state)
        static = StaticPolar(args.state, outer_order=profile.side)
        value = info['energy']+info['static_work']+static.tail_energy(e)
        return value, static, info, profile

    h = 1.e-4
    for e in (-.025, .025):
        value, static, info, profile = solve(e)
        def held(shift):
            terms = static.terms(shift)
            return terms['outer_order']+terms['outer_crystal']
        wa, x, chi = profile.side
        source = outer_crystal_source_delta(pair, labels, wa, x, e)
        received = float(np.sum(source*chi))
        crystal = static.terms(e)['outer_crystal']
        slope = (solve(e+h)[0]-solve(e-h)[0])/(2*h)
        frozen = (held(e+h)-held(e-h))/(2*h)
        emit(kind='relaxed_receiver_envelope', epsilon=e, step=h, reduced_energy=value,
             residual=info['residual'], relaxed_slope=slope, frozen_slope=frozen,
             deposited_crystal_change=received, crystal_change=crystal,
             flat_continuation=static.tail_energy(e))
        if abs(received-crystal)>1e-9 or abs(slope-frozen)>2e-6:
            raise ArithmeticError('Native polar receiver stationarity or crystal adjoint fails')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state', type=Path)
    parser.add_argument('--receiver', action='store_true')
    args = parser.parse_args()
    (check_receiver if args.receiver else run)(args)
