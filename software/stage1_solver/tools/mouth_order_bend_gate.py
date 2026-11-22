"""Archive transfer test for a local curvature-dependent order mobility.

The force condition supplies a common affinity per constituent along each
regular order level. We change only its conductance on the measured geometry.
No trial current is fed back into the saved phase or geometry.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.polynomial import polyfit

from mouth_connected_order_current import depth
from mouth_flared_receiver import FlaredDomain
from mouth_free_order_diffuse import NU, SIGMA, read
from mouth_order_wave_pulse import ARCHIVE, ROOT, digest, emit


DELTA = .04
MOBILITY = NU/(6*SIGMA*DELTA)
LEVELS = (.25, .5, .75)
CUTOFFS = (.002, .004, .008, .016)
WIDTHS = (.04, .08, .12)
METHODS = ('level_height', 'order_normal')


def smooth_fraction(b, cutoff):
    t = np.clip(b/cutoff-1, 0., 1.)
    return t*t*(3-2*t)


def cap_sample(domain, order, eta, level):
    z = depth(domain, order, eta, level)
    a, ap = domain.radius(z)
    r = eta*a
    value, cr, cz = read(domain, order, r, z)
    g = np.hypot(cr, cz)
    along = cz + eta*ap*cr
    dz = -a*cr/along
    dr = a + eta*ap*dz
    if abs(value-level) > 2e-9 or along >= 0 or dr <= 0 or cz >= 0:
        raise ArithmeticError('Order level is not an upward oriented radial graph')
    return r, z, -cr/g, g, np.hypot(dr,dz)


def bend(reference, radii, width):
    r, z, nr = reference
    # At the axis height is even and radial normal is odd.
    r = np.r_[-r[:0:-1], r]
    z = np.r_[z[:0:-1], z]
    nr = np.r_[-nr[:0:-1], nr]
    out = []
    for center in radii:
        ij = slice(np.searchsorted(r,center-width),np.searchsorted(r,center+width,side='right'))
        coord = (r[ij]-center)/width
        # Independently reconstruct principal curvatures from the graph and
        # from the actual order normal. At the outer edge only inner data exist.
        height = polyfit(coord,z[ij]-z[ij][0],4)
        normal = polyfit(coord,nr[ij],3)
        slope = height[1]/width
        norm = np.hypot(1,slope)
        graph_mer = -2*height[2]/(width**2*norm**3)
        graph_ang = -slope/(center*norm)
        normal_mer = normal[1]/width
        normal_ang = normal[0]/center
        out.append([DELTA*np.sqrt(graph_mer**2+2*graph_ang**2),
                    DELTA*np.sqrt(normal_mer**2+2*normal_ang**2)])
    return np.asarray(out).T


def control():
    r = np.linspace(0,1,1001)
    points = np.linspace(.002, .999, 237)
    # A regular 3-surface of revolution in 4-space has two angular curvatures.
    for sphere in (False, True):
        radius = 2.4
        z = np.sqrt(radius**2-r*r) if sphere else np.zeros_like(r)
        nr = r/radius if sphere else np.zeros_like(r)
        expected = DELTA*np.sqrt(3)/radius if sphere else 0.
        for width in WIDTHS:
            estimate = bend((r,z,nr), points, width)
            residual = np.max(abs(estimate-expected),axis=1)
            emit(test='geometry_control',shape='sphere' if sphere else 'plane', width=width,
                 expected=expected, max_error_height=residual[0], max_error_normal=residual[1])
            if np.max(residual) > 8e-5:
                raise ArithmeticError('Geometric calibration failed')


def prepare(mesh):
    archive = json.loads(ARCHIVE.read_text())
    carrier_path = ROOT/f'mouth_flared_cubic_n{mesh}'/'field.npz'
    suffix = '_outer_residual' if mesh==48 else ''
    level_path = ROOT/f'mouth_free_order_cubic_n{mesh}_d04{suffix}'/'order_seed1.55.npz'
    chash, ohash = digest(carrier_path), digest(level_path)
    assert chash==archive['input_carriers_sha256'][str(carrier_path)]
    assert ohash==archive['diffuse_orders_sha256'][str(level_path)]
    config = json.loads((carrier_path.parent/'stdout').read_text().splitlines()[0])
    domain = FlaredDomain(mesh,config['profile'],config['roof'],config['feed_far'])
    with np.load(carrier_path) as data:
        np.testing.assert_array_equal(domain.points, data['points'])
        carrier = data['psi']
    with np.load(level_path) as data:
        order = data['chi']
    emit(test='preparation',mesh=mesh, order_sha256=ohash, carrier_sha256=chash,
         mobility=MOBILITY, delta=DELTA, cutoffs=CUTOFFS,widths=WIDTHS,methods=METHODS)
    return domain, order, carrier


def baseline(mesh, path):
    domain, order, carrier = prepare(mesh)
    control()
    etas = np.linspace(0.,1.,801)
    v = domain.blocks[2][1]
    edges = 1-v[::-1]
    roots, weights = leggauss(6)
    quad = (edges[:-1,None]+np.diff(edges)[:,None]*(roots+1)/2).ravel()
    weights = (np.diff(edges)[:,None]*weights/2).ravel()
    saved = {'mesh':mesh,'cutoffs':CUTOFFS,'levels':LEVELS,'widths':WIDTHS}
    for k, level in enumerate(LEVELS):
        ref = np.array([cap_sample(domain,order,e,level)[:3] for e in etas]).T
        r,z,nr,grad,ds = np.array([cap_sample(domain,order,e,level) for e in quad]).T
        density = np.array([abs(read(domain,carrier,x,y)[0])**2 for x,y in zip(r,z)])
        areas = weights*4*np.pi*r*r*ds
        conductance = areas*MOBILITY*density**2/grad
        b = np.array([bend(ref,r,width) for width in WIDTHS])
        emit(test='uniform',mesh=mesh,level=level,area=sum(areas),conductance=sum(conductance),
             flat_weight_per_three_area=MOBILITY*DELTA/(level*(1-level)),
             reference_radius=ref[0,-1],axis_depth=ref[1,0],rim_depth=ref[1,-1],
             bend_range_on_halfwidth_008=[b[1].min(axis=1).tolist(), b[1].max(axis=1).tolist()])
        saved.update({f'b{k}':b,f'area{k}':areas,f'weight{k}':conductance})
    np.savez_compressed(path, **saved)
    emit(test='derived_samples',path=str(path),sha256=digest(path))


def thresholds(path, guard_dir, cutoffs):
    guard = json.loads((guard_dir/'outcome.json').read_text())
    if not (guard['limitsVerified'] and guard['exitCode']==guard['childOutcome']['exitCode']==0):
        raise ValueError('Baseline guard was not resolved')
    rows = [json.loads(line) for line in (guard_dir/'stdout').read_text().splitlines()]
    last = rows[-1]
    assert last['test']=='derived_samples' and last['sha256']==digest(path)
    emit(test='using_geometry', source=str(path),sha256=last['sha256'],
         law='M/M0=t^2*(3-2t), t=clip(b/bc-1,0,1)')
    with np.load(path) as cache:
        for k,level in enumerate(cache['levels']):
            weights = cache[f'weight{k}']
            for i,width in enumerate(cache['widths']):
                for method, b in zip(METHODS,cache[f'b{k}'][i],strict=True):
                    fractions = [weights@smooth_fraction(b,c)/sum(weights) for c in cutoffs]
                    active = [sum(weights[b>c])/sum(weights) for c in cutoffs]
                    emit(test='conditional_cap',mesh=int(cache['mesh']),level=level,halfwidth=width,
                         method=method,cutoffs=cutoffs,conductance_fraction=fractions,
                         conductance_weight_above_onset=active)
    # Previously calibrated curved intact patches with tr K * delta = .03.
    for name, principal in [('flat',[0.,0.,0.]),('sphere',[.01,.01,.01]),
                            ('collar',[.015,.015,0.]),('saddle',[.015,-.015,0.])]:
        b = np.linalg.norm(principal)
        emit(test='unpunctured_control',shape=name,principal_times_delta=principal,b=b,
             relative_mobility=[smooth_fraction(b,c) for c in cutoffs])


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    base = sub.add_parser('baseline')
    base.add_argument('mesh',type=int,choices=(26,34,48))
    base.add_argument('path',type=Path)
    gate = sub.add_parser('gate')
    gate.add_argument('path',type=Path)
    gate.add_argument('guard_dir',type=Path)
    gate.add_argument('--cutoffs', type=float, nargs='+', default=CUTOFFS)
    args = parser.parse_args()
    if args.command=='baseline':
        baseline(args.mesh,args.path)
    else:
        thresholds(args.path,args.guard_dir,args.cutoffs)
