"""Linear response of a curved stationary conversion layer to one localized site.

The monopole, dipole and quadrupole of a 4D Gaussian source decouple exactly
on the angularly symmetric empty mouth. Fields use ordinary Legendre coefficients,
with the polar axis pointing from the mouth center toward the site. All responses
are per unit h for interaction energy -h int f chi in physical volume.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import ive

from mouth_curved_phase import quadratic, tensor
from mouth_flared_receiver import FlaredDomain
from mouth_free_order_diffuse import NU, SIGMA, Order, positive_factor, read, emit


ROOT = Path('_scratch/s11c')
MEASURE = Path('research/pde_ledger_v3/_measurements')
WIDTH = .04


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(mesh):
    flow = ROOT/f'mouth_flared_cubic_n{mesh}'
    order = ROOT/f'mouth_free_order_cubic_n{mesh}_d04'
    if mesh == 48:
        order = order.with_name(order.name+'_outer_residual')
    for folder in (flow,order):
        out = json.loads((folder/'outcome.json').read_text())
        assert out['exitCode'] == out['childOutcome']['exitCode'] == 0 and out['limitsVerified']
    src, base = flow/'field.npz', order/'order_seed1.55.npz'
    setup = json.loads((flow/'stdout').read_text().splitlines()[0])
    meta = json.loads((order/'stdout').read_text().splitlines()[0])
    accepted = json.loads((MEASURE/'mouth_free_order_depth.json').read_text())
    assert sha(base) == accepted['diffuse_orders_sha256'][str(base)]
    assert sha(src) == meta['source_sha256'] and meta['delta'] == WIDTH and meta['sigma'] == SIGMA
    np.testing.assert_allclose(setup['bulk'], .25, rtol=0, atol=1e-8)
    dom = FlaredDomain(mesh,'cubic',setup['roof'],setup['feed_far'])
    with np.load(src) as a, np.load(base) as b:
        np.testing.assert_array_equal(dom.points,a['points'])
        psi, chi = a['psi'], b['chi']
    emit(test='preparation', mesh=mesh, field=str(src), field_sha256=sha(src),
         order=str(base), order_sha256=sha(base), intake=setup['bulk'], nu=NU,
         delta=WIDTH, mobility=NU/(6*SIGMA*WIDTH), dofs=len(chi))
    return dom, psi, chi


def quadrature_radius(dom):
    root = leggauss(4)[0]
    x,y = [q.ravel() for q in np.meshgrid(root,root,indexing='ij')]
    radii=[]
    for block,(pg,vg,_) in dom.blocks.items():
        i,j=[a.ravel() for a in np.indices((len(pg)-1,len(vg)-1))]
        p=pg[i,None]+np.diff(pg)[i,None]*(1+x)/2
        v=vg[j,None]+np.diff(vg)[j,None]*(1+y)/2
        radii.append(dom.chart(block,p,v)[0].ravel())
    return np.concatenate(radii)


def source(dom, psi, chi, depth, width, degree):
    """Source on the material side of the wall; resolve its size separately from the field."""
    radius = dom.radius(depth)[0]
    roots,weights = leggauss(degree)
    x,y = [q.ravel() for q in np.meshgrid(roots,roots,indexing='ij')]
    px,_ = quadratic(x)
    vy,dy = quadratic(y)
    basis = tensor(px,vy)
    pg,vg,cells = dom.blocks[2]
    # Only exponentially small axial Gaussian tails extend beyond 12 widths.
    assert depth-12*width>pg[0] and depth+12*width<pg[-1]
    i,j = [q.ravel() for q in np.indices((len(pg)-1,len(vg)-1))]
    keep = (pg[i]<depth+12*width) & (pg[i+1]>depth-12*width)
    i,j,cells = i[keep],j[keep],cells[keep]
    hp,hv = np.diff(pg)[i,None]/2,np.diff(vg)[j,None]/2
    p,v = pg[i,None]+hp*(1+x),vg[j,None]+hv*(1+y)
    a = dom.radius(p)[0]
    r = a*(1-v)
    measure = 4*np.pi*r*r*a*hp*hv*np.outer(weights,weights).ravel()
    alpha = radius*r/width**2
    gauss = np.exp(-((r-radius)**2+(p-depth)**2)/(2*width**2))/(2*np.pi*width**2)**2
    prefactor = gauss*np.sqrt(np.pi/(2*alpha))
    moments = [prefactor*((2*l+1)*ive(l+.5,alpha)) for l in range(3)]
    # Normalize in physical volume before multiplying by the integrating factor.
    norm = np.sum(measure*moments[0])
    phase = np.angle(psi[cells]@basis.T)
    physical = measure/norm
    weighted = physical*np.exp(-phase/NU)
    loads = []
    for f in moments:
        local = (weighted*f)@basis
        loads.append(np.bincount(cells.ravel(), weights=local.ravel(), minlength=len(psi)))
    cr = -(chi[cells]@tensor(px,dy).T)/(hv*a)
    source_grad = np.sum(physical*moments[1]*cr)/3
    mean_cos = np.sum(physical*moments[1])/3
    return np.array(loads), norm, source_grad, mean_cos


def front_read(dom, empty, fields, fraction):
    def location(z):
        a, ap = dom.radius(z)
        return fraction*a, fraction*ap
    def order(z):
        return read(dom,empty,location(z)[0],z)[0]
    z = brentq(lambda w: order(w)-.5, 1.0, 1.5, xtol=3e-14)
    r,rp = location(z)
    _, cr,cw = read(dom,empty,r,z)
    deriv = cw+rp*cr
    if deriv >= 0:
        raise ValueError('Wrong order crossing')
    shifts = np.array([-read(dom,val,r,z)[0]/deriv for val in fields])
    return z,r,shifts


def run(args):
    dom,psi,chi=prepare(args.mesh)
    model=Order(dom,psi,WIDTH)
    residual, physical=model.measure(chi)
    model_mobility = NU/(6*SIGMA*WIDTH)
    free = dom.free
    radius=quadrature_radius(dom)
    cuts={'entry':1.175,'half':front_read(dom,chi,[],1.)[0]}
    rhs={}
    saved=dict(points=dom.points,chi=chi)
    for label,z in cuts.items():
        a,n,g,cosine=source(dom,psi,chi,z,args.width,12)
        b,n2,g2,_=source(dom,psi,chi,z,args.width,9)
        err=np.max(np.sum(abs(a-b),axis=1)/np.sum(abs(a),axis=1))
        emit(test='site', label=label, depth=z, radius=float(dom.radius(z)[0]), width=args.width,
             physical_gaussian_fraction_inside=n, alternate_fraction=n2, quadrature_load_relative_l1=err,
             source_dipole_grad=g, alternate_source_dipole_grad=g2, source_mean_cosine=cosine)
        if err>5e-5 or abs(n-n2)>2e-5 or not .4<n<.6:
            raise ValueError('Unresolved Gaussian source')
        saved[f'source_dipole_grad_{label}']=g
        saved[f'depth_{label}']=z
        rhs[label]=model_mobility*a
    raw=abs(residual).max()
    assert raw<5e-9
    regular=dom.points[free,0]>0
    angular=model.mass_weight(NU*model.wt/radius**2)
    c_entry=read(dom,chi,dom.radius(cuts['entry'])[0],cuts['entry'])[0]
    curvature=1-6*c_entry*(1-c_entry)
    # Compare the same empty shape held by a background bias, replacing only
    # the tangent restoring curvature with its positive upstream value.
    control=model.K+model.mass_weight(model.scale*curvature*model.wt)
    if curvature<=0: raise ValueError('Expected a positive-curvature reference')
    for name, base in [('converting',physical),('uniform_restoring_curvature',control)]:
        fields={label:[] for label in cuts}
        for l in range(3):
            active=np.ones(regular.shape,dtype=bool) if l==0 else regular
            mat=(base+l*(l+1)*angular)[active][:,active]
            fac=positive_factor(mat)
            if fac is None: raise ValueError('Perturbation branch is unstable')
            loads=np.stack([v[l,free][active] for v in rhs.values()],axis=1)
            out=fac.solve(loads)
            err=np.linalg.norm(mat@out-loads,axis=0)/np.linalg.norm(loads,axis=0)
            emit(test='linear_solve', branch=name, l=l, relative_residual=err.tolist(),
                 empty_max_weak_residual=raw, control_scaled_curvature=curvature)
            if max(err)>2e-9: raise ArithmeticError('Linear residual failed')
            for col,label in enumerate(cuts):
                f=np.zeros_like(chi)
                f[np.flatnonzero(free)[active]]=out[:,col]
                fields[label].append(f)
                if name=='converting' and l<=1: saved[f'eta_{l}_{label}']=f
        for label, vals in fields.items():
            for frac in (0.,.25,.5,.75,.9,1.):
                z,r,shifts=front_read(dom,chi,vals,frac)
                emit(test='front_shift',branch=name,site=label,radial_fraction=frac,
                     empty_half_order_depth=z,radius=r,minimum_distance_to_site=
                     np.hypot(dom.radius(cuts[label])[0]-r,z-cuts[label]),
                     mean_depth_shift_per_h=shifts[0],p1_depth_shift_per_h=shifts[1],
                     p2_depth_shift_per_h=shifts[2],
                     angular_first_moment=shifts[1]/(3*shifts[0]),
                     angular_second_legendre_moment=shifts[2]/(5*shifts[0]))
    if args.fields:
        # Downstream stress workers consume these only after the guarded outcome exists.
        np.savez_compressed(args.fields,**saved)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',required=True,choices=(26,34,48),type=int)
    parser.add_argument('--width',type=float,choices=(.01,.02),default=.01)
    parser.add_argument('--fields',type=Path)
    run(parser.parse_args())
