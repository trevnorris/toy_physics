"""Read the same microscopic pair along the cubic carrier used for the middle push.

A weak displacement of the two tagged constituents is released immediately.
The stochastic material, support strip and eligibility test are inherited from
the neutral-patch experiment. Downstream motion of the patch is prescribed.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator

from mouth_bend_static_shear import archive_inputs
from mouth_contact_tracking import WallTrace
from mouth_order_wave_pulse import ROOT, digest, emit
from mouth_two_depth_load import DEPTHS
from neutral_patch_conversion import NEUTRAL, SEEDS, Sites, wall_line
from neutral_patch_sites import COUNT, TEMPERATURE, STEP


def cubic_trace(mesh):
    flow = ROOT / f'mouth_flared_cubic_n{mesh}/field.npz'
    order = ROOT / f'mouth_free_order_cubic_n{mesh}_d04/order_seed1.55.npz'
    geometry, chi, _ = archive_inputs(order, None)
    meta = json.loads((order.parent/'stdout').read_text().splitlines()[0])
    if digest(flow) != meta['source_sha256']:
        raise ValueError('Unmatched cubic carrier')
    trace = WallTrace.__new__(WallTrace)
    trace.geometry = geometry
    trace.edges, radial, cells = geometry.order_blocks[2]
    nodes = cells[::len(radial)-1, [0, 3, 6]]
    trace.chi = chi[nodes]
    with np.load(flow) as saved:
        np.testing.assert_array_equal(saved['points'], geometry.points)
        trace.psi = saved['psi'][nodes]
    return trace, meta['source_flux']


def passage(trace, refine):
    start = float(DEPTHS[1])
    levels = (.5, .25, .1, .01)
    gates = np.array([start, *(trace.crossing(c) for c in levels)])
    knots = np.unique(np.r_[gates, trace.edges[(trace.edges > start) & (trace.edges < gates[-1])]])
    w = np.r_[np.concatenate([np.linspace(a,b,17)[:-1] for a,b in zip(knots[:-1],knots[1:])]), knots[-1]]
    roots, weights = leggauss(8)
    width = np.diff(w)[:,None]
    tau = np.r_[0.,np.cumsum(np.sum(width*weights/(2*trace.speed(w[:-1,None]+width*(1+roots)/2)), axis=1))]
    times = tau[np.searchsorted(w,gates)]
    np.testing.assert_allclose(times, [trace.transit(start,stop,6) for stop in gates], atol=1e-9)
    counts = refine*np.ceil(np.diff(times)/STEP).astype(int)
    event = np.r_[0,np.cumsum(counts)]
    t = np.r_[np.concatenate([np.linspace(a,b,n+1)[:-1] for a,b,n in zip(times[:-1],times[1:],counts,strict=True)]), times[-1]]
    depth = PchipInterpolator(tau,w)(t)
    order = trace.sample(trace.chi,depth)[0]
    dilation = trace.geometry.radius(depth)[0]
    return t, depth, order, dilation/dilation[0], event


def run(mesh, refine):
    trace, intake = cubic_trace(mesh)
    t, depth, chi, scale, events = passage(trace, refine)
    k, mass, distal, left, mix, _ = wall_line(trace,80,depth,.0001)
    record = json.loads((Path('research/pde_ledger_v3/_measurements')/'neutral_patch_conversion.json').read_text())
    ids = record['independent_id_extraction']['laws']['central']['strict_surviving_ids']
    frame = json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']['central']
    tag = frame['sites']
    which = ids.index(tag)
    sites = Sites('central', ids)
    rest = np.array(frame['rested_xy'])
    rest -= rest.mean(axis=0)
    h = 1e-5
    jumps = h*np.array([0.,1.,-1.,.5,-.5])
    n, d = len(SEEDS),len(jumps)
    pos = np.broadcast_to(rest,(n*d,*rest.shape)).copy()
    # Equal/opposite shifts preserve the center and change relative slip by jumps.
    pos[:,tag[0],0] -= np.tile(jumps,n)/2
    pos[:,tag[1],0] += np.tile(jumps,n)/2
    u = np.zeros((n*d,len(mass)))
    v = np.zeros_like(u)
    live = np.ones((n,len(ids)),dtype=bool)
    onset = np.full(n,np.nan)
    below_since = np.full(n,np.nan)
    generators = [np.random.default_rng(s) for s in SEEDS]
    bridges = [np.random.default_rng(s+10000) for s in SEEDS]
    half_error = 0.
    # A first dip is sensitive to microscopic noise. Require 0.05 time
    # units continuously below half for the release-memory reading.
    emit(test='passage',intake=float(intake),mesh=mesh,refine=refine,law='central',tag=tag,
         pairs=ids,depth=depth[events].tolist(),order=chi[events].tolist(),
         elapsed=t[events].tolist(),waist_to_middle=float(trace.transit(.5,DEPTHS[1],8)),
         seeds=list(SEEDS),tap=float(h),sustained_half_window=.05)
    for j,time in enumerate(t):
        strength=chi[j]/chi[0]
        grad,pair,now,_=sites.read(pos,strength)
        live &= now[::d]
        load=pair[:,which].reshape(n,d)
        odd=(load[:,1]-load[:,2])/(2*h)
        half=(load[:,3]-load[:,4])/h
        if j==0:
            initial=odd.copy()
            if min(initial)<=0 or not live.all():
                raise ValueError('Invalid initial grip or labels')
        ratio=odd/initial
        half_error=max(half_error,float(np.max(np.abs((odd-half)/initial))))
        below_since=np.where(ratio>=.5,np.nan,np.where(np.isnan(below_since),time,below_since))
        reached=np.isnan(onset)&(time-below_since>=.05)
        onset[reached]=below_since[reached]
        if j in events:
            # Measure instantaneous tangent grip on new frozen copies along
            # the unperturbed trajectory, independent of the released tap.
            copies=pos[::d].copy()
            copies[:,tag[0],0]-=h/2
            copies[:,tag[1],0]+=h/2
            _,plus,_,_=sites.read(copies,strength)
            copies[:,tag[0],0]+=h
            copies[:,tag[1],0]-=h
            _,minus,_,_=sites.read(copies,strength)
            tangent=(plus[:,which]-minus[:,which])/(2*h)
            emit(test='reading',elapsed=float(time),depth=float(depth[j]),order=float(chi[j]),
                 survivor_counts=live.sum(axis=1).tolist(),original_tag_survives=live[:,which].tolist(),
                 tangent_over_entry=(tangent/initial).tolist(),released_over_entry=ratio.tolist())
        if j+1==len(t):
            break
        dt=t[j+1]-time
        if refine==1 or j%2==0:
            coarse=dt if refine==1 else t[j+2]-time
            wiener=np.sqrt(coarse)*np.stack([g.normal(size=rest.shape) for g in generators])
            if refine==2:
                bridge=np.sqrt(coarse)*np.stack([g.normal(size=rest.shape) for g in bridges])
                remaining=(wiener-bridge)/2
                wiener=(wiener+bridge)/2
        else:
            wiener=remaining
        a=left[j]
        b=mix[j]
        local=(1-b)*u[:,a]+b*u[:,a+1]
        grip=pos[:,:COUNT,0]-scale[j]*rest[:COUNT,0]-local[:,None]
        drift=-grad
        drift[:,:COUNT,0]-=grip
        pos=pos*(scale[j+1]/scale[j])+dt*drift+np.repeat(wiener,d,axis=0)*np.sqrt(2*TEMPERATURE)
        tension=np.diff(u,axis=1)*k
        wall=np.zeros_like(u)
        wall[:,1:]-=tension
        wall[:,:-1]+=tension
        wall[:,a]+=(1-b)*grip.sum(axis=1)
        wall[:,a+1]+=b*grip.sum(axis=1)
        wall[:,-1]-=distal*v[:,-1]
        v[:,1:]+=dt*wall[:,1:]/mass[1:]
        u+=dt*v
        if not np.isfinite(pos).all():
            raise ArithmeticError('Nonfinite patch motion')
    if half_error>1e-4:
        raise ArithmeticError('Weak probe changed with amplitude')
    emit(test='memory',half_onset=[None if np.isnan(x) else float(x) for x in onset],
         maximum_tap_halving_error_per_initial=half_error)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=(26,34),required=True)
    parser.add_argument('--refine',type=int,choices=(1,2),required=True)
    options=parser.parse_args()
    run(options.mesh,options.refine)
