"""A conditional 0PN wall-only drive with two slipping, renewing material sheets.

On the round cubic wall the vector displacement of each sheet is
q*h*e_x + a*n_x*T + b*grad_S(n_x). T is the meridional tangent. Both sheets
follow the displaced surface normally. The four tangential fields are solved.
"""
import argparse
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_calibrated_flare_coating import read_inputs
from mouth_carrier_lateral_inertia import taper
from mouth_empty_bend_stiffness import INTERVALS
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_periodic_compression import Bilayer
from neutral_patch_sites import GAP, emit, morse


def registry(patch, scale):
    cell = Bilayer(1,1,True,patch['shell'])
    d = patch['pitch']*cell.vectors[cell.cross]
    r = np.sqrt(np.sum(d*d,axis=1)+GAP**2)
    depth, alpha, length = [v[cell.cross] for v in cell.params]
    _, first = morse(r,depth,alpha,length)
    y = np.exp(-alpha*(r-length))
    second = 2*depth*alpha**2*y*(2*y-1)
    return scale*np.sum(first/r+(second-first/r)*(d[:,0]/r)**2)/PITCH**2


def assemble(trace, patch, refine):
    scale = patch['shared_scale']*2*GAP*PITCH
    moduli = patch['homogeneous_modes']
    mu = scale*moduli['area_preserving_shear']/4
    c11 = scale*(moduli['area_preserving_shear']+moduli['dilation'])/4
    grip = registry(patch,scale)
    edges = np.unique(np.r_[trace.edges, *INTERVALS[1]])
    points = np.r_[np.concatenate([np.linspace(a,b,refine,endpoint=False)
                                   for a,b in zip(edges[:-1],edges[1:])]),edges[-1]]
    size = 4*len(points)+1
    # The last coordinate is the externally displaced flare. The remaining
    # coordinates at each station are the two slips for each sheet.
    membrane, between, lining, drag, drift = [np.zeros((size,size)) for _ in range(5)]
    shape_q = np.eye(9)[-1]
    roots, weights = leggauss(6)
    geo = trace.geometry
    for e,(lo,hi) in enumerate(zip(points[:-1],points[1:])):
        w = (lo+hi)/2+(hi-lo)*roots/2
        radius, slope = geo.radius(w)
        arc = np.hypot(1.,slope)
        cs, ns = slope/arc, 1/arc
        # Differentiate within each smooth chart; the collar and flare join at an edge.
        gap = min(1.e-6,.1*np.min(np.r_[w-lo,hi-w]))
        curvature = (geo.radius(w+gap)[1]-geo.radius(w-gap)[1])/(2*gap*arc**3)
        speed = trace.speed(w)*arc
        order = np.clip(trace.sample(trace.chi,w)[0],0.,1.)
        gate, derivative = taper(w,*INTERVALS[1])
        h, hs = 1-gate,-derivative/arc
        idx = np.r_[np.arange(4*e,4*e+8),size-1]
        ix = np.ix_(idx,idx)
        for j, root in enumerate(roots):
            weight = 4*np.pi*radius[j]**2/3*arc[j]*(hi-lo)*weights[j]/2
            val = np.column_stack((np.eye(4)*(1-root)/2,np.eye(4)*(1+root)/2,np.zeros(4)))
            der = np.column_stack((-np.eye(4),np.eye(4),np.zeros(4)))/((hi-lo)*arc[j])
            a0,b0,a1,b1 = val
            def squared(x): return np.outer(x,x)
            lining[ix] += weight*grip*(squared(a0)+2*squared(b0))
            between[ix] += weight*grip*order[j]*(squared(a0-a1)+2*squared(b0-b1))
            for i in (0,2):
                a, b = val[i:i+2]
                da, db = der[i:i+2]
                us = hs[j]*shape_q
                ss = cs[j]*us+da
                tt = (cs[j]*a-b)/radius[j]
                st = us+db+(a-cs[j]*b)/radius[j]
                # Each sheet carries half the calibrated affine membrane energy.
                membrane[ix] += .5*weight*(c11*(squared(ss)+squared(tt)) + mu*squared(st)
                                  +(c11-2*mu)*(np.outer(ss,tt)+np.outer(tt,ss)))
                # Angular vector has mean square 2/3, scalar n_x has 1/3.
                geom = h[j]*shape_q
                velocity = np.array([cs[j]*geom+a, ns[j]*geom, np.sqrt(2)*(geom+b)])
                spatial = np.array([cs[j]*us+da, ns[j]*us+curvature[j]*a, np.sqrt(2)*(us+db)])
                drag[ix] += weight*velocity.T@velocity
                drift[ix] += weight*speed[j]*velocity.T@spatial
    # Both incoming lateral fields are zero at the undriven upstream station.
    return points,grip,[sparse.csc_matrix(mat[4:,4:]) for mat in (membrane,between,lining,drag,drift)]


def trajectory(points, matrices, trace, grip, corrugated, tau, nsteps, mesh, refine):
    membrane,interlayer,wall,metric,transport = matrices
    transit = trace.transit(points[0],points[-1],8)
    zeta = grip*tau*transit
    stiffness = membrane + interlayer + corrugated*wall
    drag = zeta*metric
    flow = zeta*transport
    op = stiffness+flow
    def separate(m):
        return m[:-1,:-1].tocsc(),m[:-1,-1].toarray().ravel()
    d, dq = separate(drag)
    k, kq = separate(op)
    mass_solve = splu(d).solve
    amplitudes = np.array([1.e-5,-1.e-5,.5e-5,-.5e-5])
    equilibrium = splu(k).solve(-kq)
    steady = np.r_[equilibrium,1.]
    steady_force = (op@steady)[-1]
    steady_relative = np.max(abs((op@steady)[:-1]))
    if steady_relative>1.e-8: raise ArithmeticError('Unbalanced steady slip')
    # Smoothly impose the displacement and hold it while resident material leaves.
    dt = transit/nsteps
    ramp = transit/8
    def height(t):
        if t>=ramp: return 1.,0.
        return (1-np.cos(np.pi*t/ramp))/2,np.pi*np.sin(np.pi*t/ramp)/(2*ramp)
    left = splu((d/dt+k/2).tocsc()).solve
    right = d/dt-k/2
    state = np.zeros((len(kq),len(amplitudes)))
    radius,slope = trace.geometry.radius(points)
    cs = slope/np.hypot(1.,slope)
    h = 1-taper(points,*INTERVALS[1])[0]
    station = np.argmin(abs(points-1.6))
    elapsed = np.array([trace.transit(points[0],z,6) if z>points[0] else 0. for z in points])
    worst = 0.
    emit(test='fixture',mesh=mesh,refine=refine,steps_per_transit=nsteps,wall_registry=float(corrugated*grip),
         drag_per_sheet=zeta,registry_relaxation_over_transit=tau,
         depth_limits=[points[0],points[-1]],carrier_transit=transit,
         force_unit='material contribution to external flare force per displacement',
         joint_clamp_membrane_tangent=membrane[-1,-1],steady_coating_slope=steady_force,
         readout_depth=points[station],steady_free_equation_residual=steady_relative)
    for step in range(1,4*nsteps+1):
        t = step*dt
        q_old,_ = height(t-dt)
        q,rate = height(t)
        avg = (q+q_old)/2
        state = left(right@state - np.outer(kq*avg+dq*(q-q_old)/dt,amplitudes))
        if step not in (nsteps//8,nsteps//2,nsteps,2*nsteps,4*nsteps): continue
        y = np.vstack((state,q*amplitudes))
        vel = mass_solve(-k@state - np.outer(kq*q+dq*rate,amplitudes))
        dy = np.vstack((vel,rate*amplitudes))
        force = drag@dy+op@y
        worst = max(worst,float(np.max(abs(force[:-1]/amplitudes))))
        observed = force[-1]/amplitudes
        np.testing.assert_allclose(observed,observed[0],rtol=1e-9,atol=1e-11)
        slips = np.r_[np.zeros(4),state[:,0]/amplitudes[0]].reshape(-1,4)
        angular = slips[station,[1,3]] + q*h[station]
        meridional = slips[station,[0,2]] + cs[station]*q*h[station]
        emit(test='hold',wall='lattice' if corrugated else 'smooth',tau=tau,t_over_transit=t/transit,
             signed_displacements=amplitudes.tolist(),force_per_signed_displacement=observed.tolist(),
             elastic_force=(stiffness@y)[-1,0]/amplitudes[0],reservoir_drag_force=(drag@dy+flow@y)[-1,0]/amplitudes[0],
             wall_angular_dipole=q*h[station],sheet_angular_dipoles=angular.tolist(),sheet_meridional_dipoles=meridional.tolist(),
             fresh_front_depth=np.interp(min(t,transit),elapsed,points),all_original_passed=step>=nsteps,
             steady_difference=float(observed[0]-steady_force))
    if worst>2e-8: raise ArithmeticError(f'Unbalanced free sheet {worst}')
    emit(test='balance',max_interior_residual_per_step=worst)


def main(args):
    patch,baseline,_ = read_inputs(args.mesh)
    trace,intake = cubic_trace(args.mesh)
    points,grip,mats = assemble(trace,patch,args.refine)
    if not .24999999 < intake < .25000001: raise ValueError('Intake changed')
    diag = json.loads(open('research/pde_ledger_v3/_measurements/mouth_calibrated_flare_coating.json').read())
    lines = diag['runs'][f'_scratch/s11c/flare_coating_n{args.mesh}_01']['literal_stdout'].splitlines()
    clamp, = [row['material_slope'] for row in map(json.loads,lines) if row.get('test')=='clamp'
              and row['bend']==2 and row['engagement']=='persistent']
    np.testing.assert_allclose(mats[0][-1,-1],clamp,rtol=1.e-7)
    emit(test='reference',mesh=args.mesh,intake=intake,rigidly_attached_slope=clamp,
         static_carrier_and_wave_slope=baseline[2,2],local_registry=grip)
    for tau in args.tau:
        for lining in (False,True):
            trajectory(points,mats,trace,grip,lining,tau,args.steps,args.mesh,args.refine)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=(26,34),required=True)
    parser.add_argument('--refine',type=int,default=1)
    parser.add_argument('--steps',type=int,default=128)
    parser.add_argument('--tau',type=float,nargs='+',default=[1/32,1/8])
    main(parser.parse_args())
