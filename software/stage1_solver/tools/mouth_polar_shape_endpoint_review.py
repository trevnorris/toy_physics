"""Independently inspect the geometry and signed depth/polar work at one saved R4 shape."""
import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from mouth_bonded_front_adjoint import angular_weights
from mouth_disk_independent_review import crystal_cap
from mouth_disk_pillbox import Pillbox
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_polar_joined import PolarCoupled
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_joined_front import check_guard, emit, ORDER_INERTIA


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def wall(shape, eps, p, t):
    """Wall and d wall/dt in terms of the reference spline and polynomial coordinates."""
    q = shape.params
    t = np.asarray(t)
    u = t*(1.-t)
    base = shape.reference(shape.unit_length*t)
    neutral = shape.reference(shape.unit_length*t, 1)*shape.unit_length
    exit_change = 1.5*np.expm1(q[1])
    meridian = shape.neck*(base+exit_change*t+4*q[2]*u+4*q[3]*u*(2*t-1)+16*shape.extra*u*u)
    derivative = shape.neck*(neutral+exit_change+4*q[2]*(1-2*t)
                              +4*q[3]*(-6*t*t+6*t-1)+32*shape.extra*u*(1-2*t))
    v = np.clip((t-.06)/.56, 0., 1.)
    h = (4*v*(1-v))**2
    dh = 32*v*(1-v)*(1-2*v)/.56
    return meridian+eps*p*h, derivative+eps*p*dh


def extrema(sh, eps, p):
    grid = np.unique(np.r_[np.linspace(0.,1.,2001), .06, .62])
    deriv = wall(sh, eps, p, grid)[1]
    roots = [0., 1.]
    for i in np.flatnonzero(deriv[:-1]*deriv[1:]<0):
        roots.append(brentq(lambda t: wall(sh,eps,p,t)[1],grid[i],grid[i+1]))
    ts = np.array(roots)
    radius, _ = wall(sh,eps,p,ts)
    lo, hi = np.argmin(radius), np.argmax(radius)
    return dict(minimum=float(radius[lo]), maximum=float(radius[hi]), narrowest_depth=ts[lo]*sh.length,
                widest_depth=ts[hi]*sh.length)


def solid_and_currents(trial, d, u):
    """Independent singular-value and weak-flux reconstructions on saved coefficients."""
    from numpy.polynomial.legendre import leggauss
    from mouth_joined_angular_work import bump
    b, model = trial.body, trial.model
    nodal = d.reshape(-1,2)[b.id]
    value = np.einsum('qi,qik->qk', b.basis, nodal)
    r, z = b.r+value[:,0], b.w+value[:,1]
    B = np.stack([np.einsum('qi,qik->qk',g,nodal) for g in (b.dr,b.dw)],axis=-1)+np.eye(2)
    h, hz = bump(z,b.L)
    R, Rz = trial.shape.radius(z)
    q = trial.epsilon*(r/R)**3
    s = q*h/r
    v = q*(hz-3*h*Rz/R)
    # Cholesky gives a different reference representative with the same left
    # Cauchy–Green tensor. Compose all three maps, then take singular values.
    F = np.zeros((len(r),4,4))
    F[:,:2,:2] = B@np.linalg.cholesky(b.K)
    F[:,2,2] = F[:,3,3] = b.h*r/b.r
    mu = np.r_[0.,(leggauss(8)[0]+1.)/2.,1.]
    limits = [np.inf,-np.inf,np.inf,-np.inf]
    for costheta in mu:
        P = (3*costheta**2-1)/2
        Pa = -3*costheta*np.sqrt(1-costheta**2)
        M = np.broadcast_to(np.eye(4),F.shape).copy()
        M[:,0,0],M[:,0,1],M[:,0,2] = 1+3*s*P,v*P,s*Pa
        M[:,2,2] = M[:,3,3] = 1+s*P
        full = M@F
        singular = np.linalg.svd(full,compute_uv=False)
        J = np.linalg.det(full)
        limits = [min(limits[0],J.min()),max(limits[1],J.max()),
                  min(limits[2],singular[:,-1].min()),max(limits[3],singular[:,0].max())]
    n = (model.D[:,4]*(model.V@u[0]))**2
    p = np.stack([a@u[1] for a in model.G],axis=1)
    c = model.ref[:,1:]+np.stack([a@u[2] for a in model.G],axis=1)
    pc = (p*c).sum(axis=1)
    vel = p-c*(ORDER_INERTIA*pc/(n+ORDER_INERTIA*(c*c).sum(axis=1)))[:,None]
    current = n[:,None]*vel
    flux = []
    for ids, sign in ((model.feed,-1.),(model.bulk,1.)):
        test = np.zeros(model.size)
        test[ids]=1.
        flux.append(sign*sum((model.hq*(g@test))@current[:,i] for i,g in enumerate(model.G)))
    wt = model.D[:,9]
    cap = model.D[:,7] == -1
    mach = np.sqrt((vel*vel).sum(axis=1)+(model.winding/model.D[:,6])**2)/n**2
    dense = cap&(n>.5)&(mach<1.)
    jp = np.maximum(current[:,2],0.)
    return dict(min_J=limits[0],max_J=limits[1],min_stretch=limits[2],max_stretch=limits[3],
                feed_current=flux[0],bulk_current=flux[1],dense_total_subsonic_share=wt[dense]@jp[dense]/(wt[cap]@jp[cap]))


def review(args):
    state, joint, point = args.state, args.joint, args.checkpoint
    for path in (state, joint):
        check_guard(path)
    limits = json.loads((point.parent/'limit-validation.json').read_bytes())
    if not limits['verified'] or int(limits['actual']['memory.swap.max']):
        raise ValueError('Producer limits not verified')
    history = [json.loads(s) for s in (point.parent/'stdout').read_text().splitlines()
               if s.startswith('{') and s.endswith('}')]
    prep = next(r for r in history if r['kind']=='preparation')
    blob = point.read_bytes()
    fingerprint = hashlib.sha256(blob).hexdigest()
    found = next(r for r in history if r.get('sha256')==fingerprint and r['kind']=='accepted_shape')
    material = [json.loads(s) for s in (joint.parent/'stdout').read_text().splitlines()]
    if (material[-1]['kind'] != 'held_equilibrium' or material[-1]['sha256'] != prep['joint']
            or digest(joint) != prep['joint'] or digest(state) != prep['state']
            or material[0]['source_sha256'] != prep['state'] or material[0]['winding']!=1
            or material[0]['grid'] != prep['field_grid']):
        raise ValueError('Shapes, preparation and historical equilibrium must match')
    with np.load(io.BytesIO(blob)) as saved:
        x, d, u = saved['x'], saved['displacement'], saved['field']
        shape = Meridian(4.,saved['shape_params'],extra=saved['shape_extra'].item())
        energy = saved['energy'].item()
    if args.transfer_from:
        transfer = next(r for r in history if r['kind']=='mesh_transfer')
        with np.load(args.transfer_from) as old:
            old_grid = list(np.subtract(old['displacement'].shape[:2],1)[::-1])
            if (digest(args.transfer_from)!=transfer['state'] or not np.array_equal(old['x'],x)
                    or old_grid!=transfer['source_material_mesh'] or prep['material_cells']!=transfer['target_material_mesh']):
                raise ValueError('Endpoint does not match the coarse shape recorded in the transfer')
        emit(kind='transfer_source',coarse_sha256=transfer['state'],fine_sha256=fingerprint,
             source_material_mesh=transfer['source_material_mesh'],target_material_mesh=transfer['target_material_mesh'])
    base = Front(state,tuple(prep['field_grid']),1)
    eps = x[5]
    if not np.array_equal(shape.params,base.sh.params+x[:4]*[1.,.25,.25,.25]) or (
            shape.extra != base.sh.extra+.25*x[4] or not np.array_equal(x,found['coordinates'])):
        raise ValueError('Physical shape and normalized record disagree')
    for angle, p in [('pole',1.),('equator',-.5)]:
        probe = extrema(shape, eps, p)
        emit(kind='wall_geometry',checkpoint=fingerprint, angle=angle, depth=shape.length, exit=shape.exit,
             minimum_difference=probe['minimum']-found[angle+'_minimum'],
             widest_depth_difference=probe['widest_depth']-found['widest_depth_'+angle], **probe)
    mid = [float(wall(shape,eps,p,.34)[0]) for p in (1.,-.5)]
    emit(kind='maximum_polar_split',radii=mid,quadrupole=eps,depth=.34*shape.length,
         polar_minus_equatorial=mid[0]-mid[1])
    cells, nz = prep['material_cells']
    def trial(v):
        return Trial(base,state,v,cells//3,nz,False)
    center = trial(x)
    baseline = center.energy(d,u)
    emit(kind='held_energy',recorded=energy,reconstructed=baseline,difference=baseline-energy)
    if args.transfer_from:
        verified=solid_and_currents(center,d,u)
        reference={**found['crystal'],**found['flow']}
        discrepancy=max(abs(reference[k]-value) for k,value in verified.items())
        emit(kind='independent_strain_and_flow',maximum_difference=discrepancy,**verified)
        if discrepancy>5e-8:
            raise ArithmeticError('Independent solid and flow tensors disagree with the accepted endpoint')
    grad = -np.asarray(found['global_forces'])
    directions = [] if args.flow_only else [(1.,0.),(0.,1.),(.6,-.8)]
    for depth,polar in directions:
        direction = np.array([depth,0.,0.,0.,0.,polar])
        for h in (5.e-5,2.5e-5):
            samples = [trial(x+sgn*h*direction).energy(d,u) for sgn in (-1.,1.)]
            slope = (samples[1]-samples[0])/(2*h)
            emit(kind='signed_endpoint_work',direction=direction.tolist(),h=h,slope=slope,
                 recorded_projection=grad@direction,difference=slope-grad@direction)
    if args.controls:
        check_guard(args.controls)
        rows = [json.loads(s) for s in (args.controls.parent/'stdout').read_text().splitlines()]
        flow = next(r for r in rows if r['kind']=='matched_flow')
        if flow['sha256'] != digest(args.controls) or rows[0]['checkpoint_sha256'] != fingerprint:
            raise ValueError('Flow readout and accepted shape disagree')
        with np.load(args.controls) as archive:
            pair = [archive[k] for k in ('wound','rest')]
            unwound = archive['unwound']
        if not np.array_equal(pair[0],u):
            raise ValueError('Wrong wound checkpoint')
        m, b = center.model, center.body
        pill = Pillbox(m,base.c,b,d,pair)
        src = (u[2]-pair[1][2])[m.ids[1,0][:,:,-1]]@angular_weights(m.a)
        solid = crystal_cap(base,m,src,d,10)['cap']
        for side in (-1,1):
            pressure = pill.sides(1 if side<0 else 2,.12,False,8)[0].sum()
            total = pressure+(solid if side<0 else 0.)
            reported = next(r for r in rows if r['kind']=='disk' and r['side']==side)['total']
            emit(kind='released_cap_comparison',side=side,fluid_and_order=pressure,total=total,
                 difference=total-reported)
            if abs(total-reported)>2e-8:
                raise ArithmeticError('Released disk stress disagreement')
        un = PolarCoupled(m.shape,base.c,0,base.grid,polar_shift=eps,
                          terminal_profile=m.terminal_profile,precondition=False)
        un.material_load = b.get_material_load(d[b.mask]*b.scale,un)
        retained = m.measure(u)['feed_current']/un.measure(unwound)['feed_current']
        emit(kind='retention_comparison',ratio=retained,difference=retained-flow['retained_fraction'])


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('state','joint','checkpoint'):
        p.add_argument(key,type=Path)
    p.add_argument('--controls',type=Path)
    p.add_argument('--transfer-from',type=Path,help='Check the matched coarse source, full material stretches and weak currents')
    p.add_argument('--flow-only',action='store_true',help='Check controls for an endpoint already tested for signed work')
    options=p.parse_args()
    if options.flow_only and options.controls is None:
        p.error('--flow-only requires --controls')
    review(options)
