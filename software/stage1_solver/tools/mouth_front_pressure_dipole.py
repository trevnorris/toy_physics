"""Read companion-directed moments on the three-ball at a held flat mouth, R=1."""
import json
from pathlib import Path
import shlex
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator

from incident_flow_falloff import check_current, checked_grid
from lateral_source_rate import verified_packet


MANIFEST = Path('docs/_measurements/held_force_incident_flow.json')
VOLUME = 4*np.pi/3


def grid(extent, cells, depth):
    # Coordinates only. A new field solve or eigensystem is unnecessary.
    count = extent*cells
    dx = 1/cells
    r = (np.arange(count)+0.5)*dx
    return SimpleNamespace(x=(np.arange(2*count)+0.5-count)*dx, r=r,
                           dx=dx, dz=1/depth, dv=2*np.pi*r*dx**2)


def fields(path, pins, s, centers):
    with verified_packet(path, pins) as packet:
        checked_grid(s, packet)
        if 'centers' in packet:
            np.testing.assert_array_equal(packet['centers'], centers)
            np.testing.assert_array_equal(packet['rates'], np.full(len(centers), 0.25))
        elif int(packet['d']) != -centers[0]:
            raise ValueError('Wrong separation in archive')
        currents = packet['currents']
        if len(currents) != len(centers):
            raise ValueError('Wrong number of mouths')
        for i, center in enumerate(centers):
            check_current(s, center, currents[i:i+1])
        # Only the covering stencil about the target is needed after validation.
        ix = abs(s.x) < 1+3*s.dx
        ir = s.r < 1+3*s.dx
        layers = packet['layers' if 'layers' in packet else 'field'][-4:,ix,:][:,:,ir]
        current = currents[-1,ix,:][:,ir]
    x, r = s.x[ix], s.r[ir]
    face = layers[-1]
    n = abs(face)**2
    if not np.isfinite(layers).all() or np.min(abs(layers)) <= 0:
        raise ValueError('Invalid face data')
    weight = np.maximum(1-x[:,None]**2-r**2, 0)**3
    pressure = np.expm1(5*np.log(n))/5
    carriage = current**2/n
    # At interior lateral cells use the parent stationary equation and its
    # imposed normal trace, as in mouth_normal_load_return. The cropped outer
    # cell stencils are discarded; axis symmetry sets the inner radial flux.
    faces_r = np.zeros((len(x), len(r)+1), dtype=complex)
    faces_r[:,1:-1] = (r[:-1]+s.dx/2)*np.diff(face,axis=1)/s.dx
    lateral = np.diff(faces_r,axis=1)/(s.dx*r)
    lateral[1:-1] += np.diff(face,n=2,axis=0)/s.dx**2
    fww = np.expm1(4*np.log(n))*face/2 - lateral
    fw = (-weight*(1-1/np.sqrt(n))+1j*current/n)*face
    normal = pressure+(abs(fw)**2-(face.conj()*fww).real)/2
    # Every packet has four face-adjacent layers, allowing the same second-
    # order one-sided stencil for both the isolated and paired preparations.
    log_n = np.log(abs(layers[::-1])**2)
    normal_layers = pressure+carriage-n*np.einsum('i,ijk->jk',[2,-5,4,-1],log_n)/(4*s.dz**2)
    names = ('pressure','current','normal_speed','carriage','normal_pde','normal_layers')
    samples = np.stack((pressure,current,current/n,carriage,normal,normal_layers),axis=-1)
    return x[1:-1], r[:-1], names, samples[1:-1,:-1]


def ball_moments(x, r, values, order):
    """Integrate xi=-x times the bilinear face trace over the round ball."""
    radii = np.r_[0,r]
    interpolator = RegularGridInterpolator((x,radii),np.concatenate((values[:,:1],values),axis=1))
    # Split wherever either interpolation stencil or radial upper limit changes.
    turns = np.sqrt(1-r[r<1]**2)
    bounds = np.unique(np.r_[-1,x[abs(x)<1],-turns,turns,1])
    nodes, weights = leggauss(order)
    xx = (bounds[:-1,None]+bounds[1:,None]+np.diff(bounds)[:,None]*nodes)/2
    wx = np.diff(bounds)[:,None]*weights/2
    points, coeffs = [], []
    for axial, w in zip(xx.flat,wx.flat,strict=True):
        edge = np.sqrt(1-axial*axial)
        cuts = np.r_[0,r[r<edge],edge]
        mid, span = (cuts[1:]+cuts[:-1])/2, np.diff(cuts)/2
        # Two Gauss nodes exactly integrate s times the linear interpolant.
        rr = mid[:,None]+span[:,None]*np.array([-1,1])/np.sqrt(3)
        points.extend(np.column_stack((np.full(rr.size,axial),rr.ravel())))
        coeffs.extend((-axial*w*2*np.pi*rr*span[:,None]).ravel())
    points, coeffs = np.array(points), np.array(coeffs)
    np.testing.assert_allclose(coeffs @ (-points[:,0]), VOLUME/5, rtol=2e-14)
    return coeffs @ interpolator(points)


def measure(record):
    if record['guard_exit'] or record['worker_exit'] or not record['limits_verified']:
        raise ValueError('Input run failed verification')
    reference = json.loads(record['stdout'][-1])
    d = reference['d']
    s = grid(*reference['grid'])
    args = shlex.split(record['command'])
    pins = {Path(p):v for p,v in (record['input_sha256']|record['output_fields_sha256']).items()}
    paths = [Path(args[args.index(flag)+1]) for flag in ('--baseline','--pair' if '--pair' in args else '--output')]
    x,r,names,alone = fields(paths[0],pins,s,[0])
    xp,rp,npair,paired = fields(paths[1],pins,s,[-d,0])
    np.testing.assert_array_equal(x,xp)
    np.testing.assert_array_equal(r,rp)
    if names != npair or np.max(abs(alone-alone[::-1])) > 3e-10:
        raise ValueError('Isolated reflection or readout mismatch')
    difference = paired-alone
    moments = ball_moments(x,r,difference,12)
    coarse = ball_moments(x,r,difference,6)
    np.testing.assert_allclose(coarse,moments,rtol=1e-6,atol=1e-13)
    mask = x[:,None]**2+r**2 < 1
    midpoint = np.einsum('ij,ijk->k',(-x[:,None])*mask*(2*np.pi*r*s.dx**2),difference)
    vals = dict(zip(names,map(float,moments),strict=True))
    # The paired speed is an even target shape times a single normalization.
    # A finite reading here would challenge that boundary-law prediction.
    velocity = paired[:,:,names.index('normal_speed')]
    speed_asymmetry = float(np.max(abs(velocity-velocity[::-1])))
    if speed_asymmetry > 1e-12:
        raise ValueError('Unexpected sideways odd part of normal speed')
    print(json.dumps(dict(grid=reference['grid'],d=d,
        inputs={str(p):pins[p] for p in paths},moment=vals,
        midpoint=dict(zip(names,map(float,midpoint),strict=True)),
        quadrature_difference=float(np.max(abs(moments-coarse))),
        max_paired_speed_asymmetry=speed_asymmetry,
        equivalent_linear_pressure=5*vals['pressure']/VOLUME,
        flat_momentum=reference['force'],pressure_moment_over_flat=vals['pressure']/reference['force']),
        allow_nan=False),flush=True)


if __name__ == '__main__':
    rows=[r for r in json.loads(MANIFEST.read_text())['runs'] if r['name']!='compare']
    if {r['name'] for r in rows} != {f'L{L}C{C}_d{d}' for L,C in ((48,4),(96,4),(96,6)) for d in (8,16)} or len(rows)!=6:
        raise ValueError('Missing preparation')
    for row in rows:
        measure(row)
