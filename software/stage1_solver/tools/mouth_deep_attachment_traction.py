"""Independent radial-work test of converting central bonds as deep wall support."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import linprog

from mouth_carrier_lateral_response import prepare
from mouth_deep_support_calibration import band_work, quadrature
from mouth_flared_receiver import emit
from mouth_middle_tag_grip import cubic_trace
from neutral_patch_conversion import NEUTRAL, SEEDS, Sites
from neutral_patch_sites import BETWEEN, COUNT, GAP, STEP, TEMPERATURE, morse

LEDGER = Path('research/pde_ledger_v3/_measurements')
START, END = 1.175, 2.
BANDS = ((1.22,1.32), (1.34,1.42), (1.44,1.50), (1.52,1.59),
         (1.62,1.69), (1.71,1.78), (1.81,1.89), (1.91,1.99))
TRAIN = [0,3,7]


def patch_clock(trace, refinement):
    grid = trace.edges[(trace.edges > START) & (trace.edges < END)]
    edges = np.r_[START, grid, END]
    w = np.r_[np.concatenate([np.linspace(a,b,17,endpoint=False) for a,b in zip(edges[:-1],edges[1:])]), END]
    widths = np.diff(w)
    nodes, weights = np.polynomial.legendre.leggauss(6)
    samples = w[:-1,None]+widths[:,None]*(nodes+1)/2
    time = np.r_[0., np.cumsum(widths*(1/trace.speed(samples)@weights)/2)]
    np.testing.assert_allclose(time[-1], trace.transit(START,END,8), rtol=2e-8)
    t = np.linspace(0.,time[-1],refinement*int(np.ceil(time[-1]/STEP))+1)
    return t, PchipInterpolator(time,w)(t)


def pair_observables(delta):
    depth, alpha, rest = BETWEEN
    r = np.sqrt(GAP**2+np.sum(delta**2,axis=-1))
    _, tension = morse(r,*BETWEEN)
    normal = tension * GAP/r
    y = np.exp(-alpha*(r-rest))
    radial_k = 2*depth*alpha**2*y*(2*y-1)
    tangent_k = tension/r + (radial_k-tension/r)*(delta[...,0]/r)**2
    return normal.sum(axis=-1), tangent_k.sum(axis=-1)


def energy_check(delta):
    def energy(x, normal):
        shifted = delta.copy()
        shifted[...,0] += x
        r = np.sqrt((GAP+normal)**2+np.sum(shifted**2,axis=-1))
        return morse(r,*BETWEEN)[0].sum(axis=-1)
    normal, tangent = pair_observables(delta)
    for h in (1e-4,5e-5):
        np.testing.assert_allclose((energy(0,h)-energy(0,-h))/(2*h),normal,rtol=2e-6)
        np.testing.assert_allclose((energy(h,0)+energy(-h,0)-2*energy(0,0))/h**2,
                                   tangent,rtol=3e-5,atol=2e-5)


def passage(trace, refinement):
    t, w = patch_clock(trace,refinement)
    order = np.clip(trace.sample(trace.chi,w)[0],0.,1.)
    strength = order/order[0]
    a, _ = trace.geometry.radius(w)
    record = json.loads((LEDGER/'neutral_patch_conversion.json').read_text())
    ids = record['independent_id_extraction']['laws']['central']['strict_surviving_ids']
    frame = json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']['central']
    rest = np.array(frame['rested_xy'])
    rest -= rest.mean(axis=0)
    sites = Sites('central',ids)
    pos = np.broadcast_to(rest,(len(SEEDS),*rest.shape)).copy()
    patch = sites.patch
    i,j = patch.i[patch.cross],patch.j[patch.cross]
    assert len(i)==COUNT**2 and len(ids)==9
    dt = t[1]
    rng = [np.random.default_rng(seed) for seed in SEEDS]
    bridge = [np.random.default_rng(seed+10000) for seed in SEEDS]
    all_live = np.ones((len(SEEDS),len(ids)),bool)
    first_loss = None
    survivors = np.full(len(t),len(ids),dtype=int)
    # Store every eighth coarse step and the end. Geometry and clock use all steps.
    samples, normal, tangent = [], [], []
    for n in range(len(t)):
        grad, _, live, _ = sites.read(pos,strength[n])
        all_live &= live
        survivors[n] = all_live.sum(axis=1).min()
        if first_loss is None and not np.all(all_live):
            first_loss = w[n]
        if n%(8*refinement)==0 or n+1==len(t):
            delta = pos[:,j]-pos[:,i]
            values, stiffness = pair_observables(delta)
            if not samples or n+1==len(t):
                energy_check(delta)
            samples.append(n)
            normal.append(strength[n]*values)
            tangent.append(strength[n]*stiffness)
        if n+1==len(t):
            break
        if n%refinement==0:
            dw = np.sqrt(refinement*dt)*np.stack([g.normal(size=rest.shape) for g in rng])
            if refinement==2:
                other = np.sqrt(2*dt)*np.stack([g.normal(size=rest.shape) for g in bridge])
                remaining = (dw-other)/2
                dw = (dw+other)/2
        else:
            dw = remaining
        # The rigid wall provides the lower sites' unit tangential tethers.
        # Normal separation is constrained. All 144 cross bonds load this clamp.
        grad[:,:COUNT,0] += pos[:,:COUNT,0] - (a[n]/a[0])*rest[:COUNT,0]
        pos = pos*(a[n+1]/a[n])-dt*grad+np.sqrt(2*TEMPERATURE)*dw
    counts = [int(survivors[min(np.searchsorted(w,hi),len(w)-1)]) for _,hi in BANDS]
    return w[samples], np.array(normal), np.array(tangent), first_loss, counts


def best_possible(curvature, pressure, lower, upper):
    # Unknown pressure zero, area tension, nonnegative patch scale and Chebyshev
    # error. Each band can independently choose any traction inside its envelope.
    ones = np.ones_like(curvature)
    below = np.column_stack((ones,curvature,lower,-ones))
    above = np.column_stack((-ones,-curvature,-upper,-ones))
    matrix = np.vstack((below,above))
    rhs = np.r_[pressure,-pressure]
    cost = np.array([0.,0.,0.,1.])
    fit = linprog(cost,A_ub=matrix,b_ub=rhs,bounds=[(None,None),(None,None),(0,None),(0,None)],
                  method='highs')
    if not fit.success:
        raise ValueError(fit.message)
    # Primal and dual objectives certify the numerical minimax value.
    dual = fit.ineqlin.marginals
    gradient = cost-matrix.T@dual
    np.testing.assert_allclose(gradient[:2],0.,atol=2e-8)
    assert max(dual)<1e-10 and min(gradient[2:])>=-2e-8
    np.testing.assert_allclose(fit.fun,rhs@dual,atol=2e-9)
    assert np.max(matrix@fit.x-rhs)<2e-8
    return float(fit.fun)


def run(mesh, refinement):
    trace, intake = cubic_trace(mesh)
    assert abs(intake-.25)<1e-8
    w, normal, tangent, first_loss, counts = passage(trace,refinement)
    # The lower force bound needs registration only in the two upstream bands.
    # At >=60% well depth, r<=r_max; no lower bound is needed farther along.
    depth, alpha, rest = BETWEEN
    rmax = rest-np.log(1-np.sqrt(1-.6))/alpha
    lower_one = min(morse(x,*BETWEEN)[1] for x in (GAP,rmax))*GAP/rmax
    upper_one = depth*alpha/2
    emit(test='patch', mesh=mesh, refinement=refinement, intake=intake, start=START, end=END,
         first_sampled_loss=first_loss, surviving_through_band=counts, seeds=SEEDS, upper_per_bond=upper_one,
         registered_lower_per_bond=lower_one, fixed_normal_gap=GAP,
         clock='carrier-following patch', total_pairs=COUNT**2)
    carrier,field = prepare(SimpleNamespace(source=Path(f'_scratch/s11c/mouth_flared_cubic_n{mesh}')))
    np.testing.assert_array_equal(carrier.points,trace.geometry.points)
    entry = np.array([START])
    a0,s0 = trace.geometry.radius(entry)
    flux_factor = (a0**2*trace.speed(entry)*np.hypot(1.,s0)).item()
    entry_order = trace.sample(trace.chi,entry)[0].item()
    curves, actual, interior, means, grips, lower, upper = [],[],[],[],[],[],[]
    for band,(lo,hi) in enumerate(BANDS):
        cur, surf, bulk = band_work(carrier,field,lo,hi)
        z,dz = quadrature(trace.edges,lo,hi)
        radius,slope = trace.geometry.radius(z)
        x = (z-lo)/(hi-lo)
        weights = dz*radius**2*(x*(1-x))**2
        weights /= weights.sum()
        fixed = np.ones_like(z)
        passing = flux_factor/(radius**2*trace.speed(z)*np.hypot(1.,slope))
        scale = np.clip(trace.sample(trace.chi,z)[0],0.,1.)/entry_order
        means.append([[float(weights @ (area*np.interp(z,w,n))) for n in normal.T]
                      for area in (fixed,passing)])
        grips.append([[float(weights @ (area*np.interp(z,w,k))) for k in tangent.T]
                      for area in (fixed,passing)])
        lower.append(float(weights @ (np.minimum(fixed,passing)*scale))*
                     (counts[band]*lower_one if band<=1 else 0.))
        upper.append(float(weights @ (np.maximum(fixed,passing)*scale))*COUNT**2*upper_one)
        curves.append(cur); actual.append(surf); interior.append(bulk)
    cur, p, lo, hi = map(np.array,(curves,actual,lower,upper))
    means,grips = np.array(means),np.array(grips)
    assert np.all(means>=lo[:,None,None]-1e-9) and np.all(means<=hi[:,None,None]+1e-9)
    emit(test='band_measurements',bands=BANDS,curvature=curves,pressure=actual,volume_pressure=interior,
         contact_normal_per_reference_patch=means.tolist(),instant_tangent_per_reference_patch=grips.tolist(),
         density_models=['fixed per actual wall area','conserved patch flux per wall area'],
         normal_lower=lower,normal_upper=upper, sample_order=['band','density model','thermal history'])
    for coverage, data in zip(('fixed','passing'),means.mean(axis=2).T):
        design = np.column_stack((np.ones(len(cur)),cur,data))
        params = np.linalg.solve(design[TRAIN],p[TRAIN])
        residual = p-design@params
        emit(test='prediction',coverage=coverage,training=TRAIN,backing=params[0],tension=params[1],
             signed_pressure_per_patch_force=params[2],residual=residual.tolist())
    # Include the opposite normal orientation as a sensitivity test. Choosing
    # a pressure datum cannot change a residual because the intercept is free.
    for sign in (1,-1):
        L,U = (lo,hi) if sign==1 else (-hi,-lo)
        synthetic = -.001+.002*cur+.0001*(L+U)/2
        assert best_possible(cur,synthetic,L,U)<1e-9
        emit(test='any_admissible_profile',direction=sign,pressure_residual_floor=best_possible(cur,p,L,U),
             volume_residual_floor=best_possible(cur,np.array(interior),L,U))
    emit(test='complete',mesh=mesh,refinement=refinement)


if __name__=='__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--mesh',required=True,type=int,choices=(26,34))
    cli.add_argument('--refinement',required=True,type=int,choices=(1,2))
    args = cli.parse_args()
    run(args.mesh,args.refinement)
