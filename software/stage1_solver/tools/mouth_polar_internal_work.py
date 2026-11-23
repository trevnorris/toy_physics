"""An independent radial quadrupole of the crystal at fixed polar wall and Eulerian order.

A material site on the pre-polar meridional chart (r,z) goes to
R_e(r,z,theta)+s*x³(1−x)*h(z)*P2, x=r/b(z). Both materials share this map.
The wall, axis, axial cuts, front, carrier and receiver retain their positions.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_bonded_material_motion import Bonded
from mouth_joined_angular_work import bump
from mouth_polar_bonded import angles, geometry, order_at
from mouth_polar_joined import PolarCoupled
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_joined_front import file_sha, check_guard, emit


def accepted(state, point):
    # An immutable accepted snapshot may be examined during its enclosing
    # shape sweep. Check admission and its own completed, hash-identified row.
    context = point.parent
    limits = json.loads((context/'limit-validation.json').read_text())
    if not limits['verified'] or int(limits['actual']['memory.swap.max']):
        raise ValueError('Unverified shape run')
    rows = list(map(json.loads, (context/'stdout').read_text().splitlines()))
    prep = rows[0]
    sha = file_sha(point)
    hit = next(row for row in rows if row['kind']=='accepted_shape' and row['sha256']==sha)
    if prep['winding'] != 1 or prep['state'] != file_sha(state):
        raise ValueError('Snapshot and preparation differ')
    with np.load(point) as data:
        if int(data['iteration']) != hit['iteration'] or data['state_sha'].item() != prep['state']:
            raise ValueError('Snapshot content and accepted metadata differ')
        return (Meridian(4., data['shape_params'], extra=float(data['shape_extra'])), float(data['x'][5]),
                data['displacement'], data['field'], tuple(prep['field_grid']), sha)


def inspect(state, point):
    sh, e, disp, field, grid, sha = accepted(state, point)
    body = Bonded(state, (disp.shape[1]-1)//3, disp.shape[0]-1, shape=sh, full_depth=True)
    m = PolarCoupled(sh, body.grid, 1, grid, polar_shift=e, precondition=False)
    corr = field[2,m.ids[1,0]]
    depths = m.z[1]/sh.length
    support = (depths[:-1] < .62) & (depths[1:] > .06)
    front_overlap = max(abs(corr[:,:,:-1][:,:,support]).max(),abs(corr[:,:,1:][:,:,support]).max())
    if front_overlap != 0 or .62*sh.length >= sh.length-m.h:
        raise ArithmeticError('The material probe moves the reciprocal front load')
    nodes = disp.reshape(-1, 2)[body.id]
    val = np.einsum('qi,qik->qk', body.basis, nodes)
    B = np.eye(2) + np.stack([np.einsum('qi,qik->qk', g, nodes)
                              for g in (body.dr, body.dw)], axis=-1)
    r, z = body.r+val[:,0], body.w+val[:,1]
    mask = (z > .06*sh.length) & (z < .62*sh.length)
    j = (body.j*np.linalg.det(B)*(r/body.r)**2)[mask]
    hoop = (body.h*r/body.r)[mask]
    K = (B @ body.K @ B.transpose(0, 2, 1))[mask]
    rr, rz, zz = K[:,0,0], K[:,0,1], K[:,1,1]
    r, z, wt = r[mask], z[mask], body.measure[mask]
    b, bw = sh.radius(z)
    h, hw = bump(z, sh.length)
    x = r/b
    k = x**3*(1-x)
    kx = x**2*(3-4*x)
    phi = k*h
    derivative, angular, axial = kx*h/b, phi/r, k*hw-kx*x*bw/b*h
    chi, cr, _ = order_at(body, r, z)
    d, s, v = [e*a for a in geometry(sh, r, z, False)]
    quad = angles(10)
    steps = (-.0002, -.0001, .0001, .0002)
    delta = np.zeros(len(steps))
    strain = occupancy = 0.
    min_j = np.inf
    for p, pa, w in zip(*quad):
        a, g, c = 1+d*p, 1+s*p, v*p
        J = j*a*g*g
        trace = a*a*rr+2*a*c*rz+(1+c*c)*zz+hoop**2*((s*pa)**2+2*g*g)
        W = (trace-4)/2-np.log(J)+.25*(J-1)**2
        trhalf = a*derivative*p*rr+(derivative*p*c+a*axial*p)*rz+c*axial*p*zz
        trhalf += hoop**2*(s*angular*pa**2+2*g*angular*p)
        dlog = derivative*p/a+2*angular*p/g
        strain += w*(wt @ (chi*(trhalf+(-1+.5*J*(J-1))*dlog)))
        occupancy += w*(wt @ (W*cr*phi*p/a))
        for i, step in enumerate(steps):
            aa, gg, cc = a+step*derivative*p, g+step*angular*p, c+step*axial*p
            dest = r*gg
            # Invert the held Eulerian map at the moved material quadrature
            # point. Occupancy is not advected with this internal probe.
            pre = r.copy()
            for _ in range(5):
                pre -= (pre + e*(pre/b)**3*h*p-dest)/(1+3*e*pre**2*h*p/b**3)
            if np.max(abs(pre+e*(pre/b)**3*h*p-dest)) > 2e-14:
                raise ArithmeticError('Inverse polar chart residual')
            trial_chi, _, _ = order_at(body, pre, z)
            trial_J = j*aa*gg*gg
            min_j = min(min_j, trial_J.min())
            trial_trace = aa**2*rr+2*aa*cc*rz+(1+cc**2)*zz
            trial_trace += hoop**2*(((s+step*angular)*pa)**2+2*gg**2)
            trial_W = (trial_trace-4)/2-np.log(trial_J)+.25*(trial_J-1)**2
            delta[i] += w*(wt @ ((trial_chi-chi)*trial_W+chi*(trial_W-W)))
    emit(kind='fixed_wall_internal_P2', state=file_sha(state), accepted=sha, snapshot=str(point),
         wall_amplitude=e, field_grid=grid, material_cells=[disp.shape[1]-1,disp.shape[0]-1],
         front_overlap=float(front_overlap), minimum_trial_J=min_j, steps=steps, signed_changes=delta.tolist(),
         strain_derivative=float(strain), occupancy_derivative=float(occupancy),
         energy_derivative=float(strain+occupancy), generalized_force=float(-strain-occupancy),
         coarse_difference=(delta[3]-delta[0])/.0004, fine_difference=(delta[2]-delta[1])/.0002)
    if min_j<=0 or abs((delta[2]-delta[1])/.0002-strain-occupancy) > 1e-6:
        raise ArithmeticError('Independent angular load and signed energy disagree')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('state', type=Path)
    p.add_argument('snapshots', nargs='+', type=Path)
    a = p.parse_args()
    check_guard(a.state)
    for point in a.snapshots:
        inspect(a.state, point)
