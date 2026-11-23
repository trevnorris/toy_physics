"""Independent cap tractions and receiving-order traces at saved R4 endpoints."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_bonded_front_adjoint import angular_weights
from mouth_joined_front_audit import Front
from mouth_joined_front_trace import face
from mouth_poly_shallow_release import CUT
from mouth_receiver_order_control import side_source
from mouth_receiver_order_equilibrium import receiver
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_joined_front import check_guard, emit
from mouth_supercore_terminal_order import terminal


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def crystal_cap(base, m, mean, solid, quadrature):
    """One-sided anchored top edge of the displacement and inverse-label triangles."""
    sh = m.shape
    if np.any(solid[-1]) or np.any(base.labels[-1, :, 1] != 1.):
        raise ValueError('Direct recovery requires zero cap displacement and constant cap reference depth')
    nz, nr = np.subtract(solid.shape[:2], 1)
    n = nr//3
    p = np.exp(base.co.item())
    edge = np.exp(-p/2)
    material_x = np.r_[np.linspace(0., edge, 2*n+1)[:-1], np.linspace(edge, 1., n+1)]
    collar_e = np.linspace(0., .5, base.na+1)
    boundaries = np.unique(np.r_[m.rr[0], material_x, np.exp(-p*collar_e), .75, .95])
    boundaries = boundaries[np.r_[True, np.diff(boundaries) > 2e-14]]
    q, weights = leggauss(quadrature)
    x = (boundaries[:-1, None]+(1.+q)/2.*np.diff(boundaries)[:, None]).ravel()
    w = (weights/2.*np.diff(boundaries)[:, None]).ravel()
    cell = np.searchsorted(material_x, x)-1
    # Top-edge displacement vanishes identically. Its one-sided depth
    # derivative equals the top minus bottom difference on the triangle's east edge.
    stretch = 1.-solid[-2, cell+1, 1]*nz/sh.length
    outer = x >= edge
    kr = sh.length/base.ref.length
    qratio = sh.exit*edge/(base.ref.exit*np.exp(-.5*np.exp(base.co0.item())))
    zz = np.full(x.size, kr**2)
    jac = np.full(x.size, qratio**3*kr)
    # At the collar cap T_e = 0, F_wZ = L/T_t and det(Dlabels) = S_e*T_t.
    radius = sh.exit*x[outer]
    ecoord = -np.log(x[outer])/p
    col = np.searchsorted(collar_e, ecoord)-1+base.na
    Tt = (1.-base.labels[-2, col+1, 1])*base.nt
    Se = (base.labels[-1, col+1, 0]-base.labels[-1, col, 0])*(2*base.na)
    eta_nodes = np.linspace(-.5, .5, 2*base.na+1)
    S = sh.neck+np.interp(ecoord, eta_nodes, base.labels[-1, :, 0])
    jac[outer] = p*radius*sh.length*(radius/S)**2/(Se*Tt)
    zz[outer] = (sh.length/Tt)**2
    zz *= stretch**2
    jac *= stretch
    if min(jac) <= 0:
        raise ValueError('Inverted solid on the disk')
    stress = -np.interp(x, m.rr[0], mean)*(zz-1.+.5*jac*(jac-1.))/jac
    weight = w*4*np.pi*sh.exit**3*x**2
    ramp = np.clip((x-.75)/.2, 0., 1.)
    return dict(core=(weight*stress)[~outer].sum(), collar=(weight*stress)[outer].sum(),
                cap=weight@stress, compact=weight@(stress*(1.-3*ramp**2+2*ramp**3)))


def review(a):
    for p in (a.state, a.joint):
        check_guard(p)
    with np.load(a.joint) as archive:
        solid = archive['displacement']
    state_hash, joint_hash = sha(a.state), sha(a.joint)
    source_rows = [json.loads(line) for line in (a.joint.parent/'stdout').read_text().splitlines()]
    if (source_rows[0]['source_sha256'] != state_hash or source_rows[-1]['kind'] != 'held_equilibrium'
            or source_rows[-1]['sha256'] != joint_hash):
        raise ValueError('Joint source and shape are not the matched R4 endpoint')
    for fields, audit in a.endpoint:
        for p in (fields, audit):
            check_guard(p)
        rows = [json.loads(line) for line in audit.read_text().splitlines()]
        prep = rows[0]
        if (prep['source_sha256'] != state_hash or prep['joint_sha256'] != joint_hash
                or prep['fields_sha256'] != sha(fields)):
            raise ValueError('Check filenames against the physical endpoints used by the pillbox')
        grid, winding = tuple(prep['grid']), prep['winding']
        base = Front(a.state, grid, winding)
        model = base.joined(base.sh, False)
        with np.load(fields) as saved:
            pair = [saved[k] for k in ('fields', 'rest')]
        row = next(v for v in rows if v['kind'] == 'pillbox' and v['radial'] == 'full')
        for record, sign in zip(row['raw_and_weak'], (-1, 1)):
            refs = [face(model, base.c, u, sign, 12) for u in pair]
            stress = refs[0][0] @ (refs[0][-1]-refs[1][-1])
            difference = stress-sum(record['raw_cap'])
            emit(kind='independent_face', grid=grid, winding=winding, side=sign,
                 stress=stress, error=difference)
            if abs(difference) > 2e-8:
                raise ArithmeticError('Independent carrier and order cap stresses disagree')
        dc = (pair[0][2]-pair[1][2])[model.ids[1, 0][:, :, -1]]
        mean = dc@angular_weights(model.a)
        emit(kind='material_cap_trace', grid=grid, winding=winding, delta_chi_max=abs(dc).max(),
             delta_chi_rim_max=abs(dc[-1]).max(), angular_mean_peak=abs(mean).max())
        for quad in (6, 10):
            found = crystal_cap(base, model, mean, solid, quad)
            emit(kind='independent_crystal', winding=winding, grid=grid, quadrature=quad, **found)
            for radius, key in [('full', 'cap'), ('interior_to_.95', 'compact')]:
                ref = next(r for r in rows if r['kind'] == 'sharp_solid' and
                           r['quad'] == quad and r['radial'] == radius)
                if abs(ref['crystal_axial']-found[key]) > 2e-8:
                    raise ArithmeticError('Direct one-sided crystal stress differs across implementations')
    if a.receiver:
        base = Front(a.state, (13, 10, 14), 1)
        for shift in (0., .005):
            q = base.sh.params.copy()
            q[3] += shift
            sh = Meridian(4., q, extra=base.sh.extra)
            source = side_source(base, sh)
            def load(w, x):
                return source
            previous = terminal(sh, base.order_mesh, CUT, load)
            new, info = receiver(sh, base.order_mesh, CUT, load)
            t, x = np.meshgrid(np.linspace(0., 1., 91), np.linspace(0., 1., 91), indexing='ij')
            dome = abs(previous(x, t*.6)-new(x, t*.6)).max()
            diff = []
            for upper, z in [(False, t*sh.length), (True, sh.length+t*.6)]:
                b = sh.radius(z)[0]
                rr = b+x*(CUT-b)
                diff.append(abs(previous.outside(rr,z,upper)-new.outside(rr,z,upper)).max())
            emit(kind='native_receiver_traces', dq=shift, terminal=dome, side=diff[0],
                 receiving=diff[1], projected_residual=info['residual'])
            if max(dome, *diff) > 1e-5:
                raise ArithmeticError('New and native receiving-order traces differ')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state', type=Path)
    parser.add_argument('joint', type=Path)
    parser.add_argument('--endpoint', required=True, nargs=2, type=Path, action='append',
                        metavar=('FIELDS', 'AUDIT_STDOUT'))
    parser.add_argument('--receiver', action='store_true')
    review(parser.parse_args())
