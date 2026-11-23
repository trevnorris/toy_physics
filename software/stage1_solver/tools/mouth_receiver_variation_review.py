"""Signed full-energy work in a native receiver-order direction outside the mobile front."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_bonded_front_adjoint import Coupled
from mouth_internal_collar_check import read
from mouth_joined_angular_work import bump
from mouth_polar_bonded import order_sites
from mouth_poly_shallow_release import CUT
from mouth_receiver_order_equilibrium import receiver
from mouth_sealed_shape_relax import DELTA
from mouth_supercore_joined_front import emit, check_guard, file_sha


def run(state, joint):
    for path in (state, joint):
        check_guard(path)
    rows = list(map(json.loads, (joint.parent/'stdout').read_text().splitlines()))
    if (rows[-1]['sha256'] != file_sha(joint) or rows[-1]['kind'] != 'held_equilibrium'
        or rows[0]['source_sha256'] != file_sha(state) or rows[0]['rest'] or rows[0]['winding'] != 1):
        raise ValueError('Require the matched held driven R4 field')
    pair, labels, inner, outer = read(state)
    source = pair.evaluate(labels, inner, outer)['outer_source']
    sh = pair.shape
    profile, info = receiver(sh, len(inner[1])-1, CUT, lambda w, x: source)
    with np.load(joint) as data:
        u = data['field']
    m = Coupled(sh, inner[2], 1, tuple(rows[0]['grid']), terminal_profile=profile, precondition=False)
    m.material_load = np.zeros(m.size)  # The varied outer side has no inner material load.
    wa, x, c = profile.side
    h, _ = bump(wa, sh.length)
    # Vary only free cells of the receiver; vanishing at pegged nodes and where
    # chi=0 or 1 makes both signs feasible. Includes the carrier's coarse nodes.
    direction = (h[:, None]*c*(1-c))
    direction[:, [0, -1]] = 0
    direction[[0, -1], :] = 0
    direction /= np.max(abs(direction))
    b, _ = sh.radius(wa)
    interp = RegularGridInterpolator((wa, x), direction, bounds_error=False, fill_value=0.)
    xs, _, zs = np.meshgrid(m.rr[1], m.a, m.z[1], indexing='ij')
    b, _ = sh.radius(zs)
    nodes = np.zeros(m.size)
    nodes[m.ids[1, 1]] = interp((zs, xs*(28.-b)/(CUT-b)))
    ref0 = m.ref.copy()
    dref = np.column_stack([m.V@nodes, *[G@nodes for G in m.G]])
    if np.max(abs(dref[m.D[:, 8].astype(bool)])) > 1e-12:
        raise ValueError('Receiver trial overlaps mobile order')
    emit(kind='input', state=file_sha(state), joint=file_sha(joint), receiver_residual=info['residual'],
         native_peak=float(abs(direction).max()), carrier_grid_peak=float(abs(nodes).max()))
    for step in (2e-4, 1e-4):
        signs = []
        for sign in (-1., 1.):
            delta = sign*step
            cnew = c+delta*direction
            if cnew.min() < 0 or cnew.max() > 1:
                raise ValueError('Infeasible signed variation')
            (wt, gr, gz, well, *_), _ = order_sites(sh, (wa/sh.length, x, cnew), True)
            static = wt @ (DELTA**2*(gr**2+gz**2)+well)+np.sum(source*cnew)
            m.ref = ref0+delta*dref
            carrier = m.evaluate(u)[0]
            signs.append([static, carrier])
        grad = (np.array(signs[1])-signs[0])/(2*step)
        emit(kind='exterior_reference_variation', step=step, minus=signs[0], plus=signs[1],
             native_and_crystal=float(grad[0]), carrier_and_front=float(grad[1]), total=float(grad.sum()))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('state', 'joint'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    run(a.state, a.joint)
