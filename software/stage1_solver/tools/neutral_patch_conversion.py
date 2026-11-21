"""Follow previously registered sites along the native 0PN wall transit.

The conversion dependence, wall attachments and scalar transmission line are
trial material properties. The unperturbed carrier supplies the clock, local
angular dilation and spatial order. Only the weak lateral drive changes sign.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator

from mouth_contact_tracking import WallTrace
from mouth_flared_receiver import ROUNDING
from mouth_order_wave_pulse import emit, digest
from neutral_patch_sites import COUNT, GAP, PITCH_UNIT, Patch, TEMPERATURE, morse


LEVELS = (.95, .75, .5, .3, .1, .01)
SEEDS = (8401, 8402, 8403, 8404)
DRIVES = np.array([0., .02, -.02, .01, -.01])
NEUTRAL = Path("research/pde_ledger_v3/_measurements/neutral_patch_sites.json")


def clock(trace, refinement):
    gates = np.array([trace.crossing(c) for c in LEVELS])
    edges = np.unique(np.r_[gates, trace.edges[(trace.edges > gates[0]) & (trace.edges < gates[-1])]])
    z = np.r_[np.concatenate([np.linspace(a,b,17)[:-1] for a,b in zip(edges[:-1],edges[1:])]), edges[-1]]
    nodes, weights = leggauss(6)
    widths = np.diff(z)[:,None]
    samples = z[:-1,None]+widths*(1+nodes)/2
    elapsed = np.r_[0., np.cumsum(np.sum(widths*weights/(2*trace.speed(samples)),axis=1))]
    marks = elapsed[np.searchsorted(z,gates)]
    np.testing.assert_allclose(marks,[trace.transit(gates[0],w,8) for w in gates],atol=1e-10)
    counts = refinement*np.ceil(np.diff(marks)/.0005).astype(int)
    breaks = np.r_[0,np.cumsum(counts)]
    time = np.r_[np.concatenate([np.linspace(a,b,n+1)[:-1] for a,b,n in zip(marks[:-1],marks[1:],counts)]),marks[-1]]
    w = PchipInterpolator(elapsed,z)(time)
    radius = trace.geometry.radius(w)[0]
    return time, w, trace.sample(trace.chi,w)[0], radius/radius[0], breaks


def wall_line(trace, pieces, w, end_level):
    grid = np.linspace(ROUNDING,trace.crossing(end_level),pieces+1)
    x, weights = leggauss(6)
    mid = (grid[:-1]+grid[1:])/2
    half = np.diff(grid)/2
    points = mid[:,None]+half[:,None]*x
    ds = np.hypot(1,trace.geometry.radius(points)[1])
    length = half*(ds@weights)
    # Linear displacement basis in w, physical tension chi and mass per arc length 1.
    stiffness = ((trace.sample(trace.chi,points)[0]/ds)@weights)/(4*half)
    mass = np.zeros(pieces+1)
    mass[:-1] += half*(ds@((1-x)*weights/2))
    mass[1:] += half*(ds@((1+x)*weights/2))
    distal = np.sqrt(trace.sample(trace.chi,grid[-1])[0])
    left = np.searchsorted(grid,w,side='right')-1
    if left.min() < 1 or left.max() >= pieces or min(stiffness) <= 0:
        raise ValueError("Invalid wall or parcel overlap with a boundary")
    fraction = (w-grid[left])/(grid[left+1]-grid[left])
    return stiffness,mass,distal,left,fraction,float(length.sum())


class Sites:
    def __init__(self, law, identities):
        self.patch = p = Patch(law)
        self.incidence = np.zeros((len(p.i),2*COUNT))
        self.incidence[np.arange(len(p.i)),p.i] = -1
        self.incidence[np.arange(len(p.i)),p.j] = 1
        self.a, self.b = np.array(identities).T
        self.idx = np.array([np.flatnonzero((p.i==a)&(p.j==b))[0] for a,b in identities])

    def read(self, pos, scale):
        p = self.patch
        delta = pos[:,p.j]-pos[:,p.i]
        tang = np.einsum('bik,bik->bi',delta,delta)
        cross = tang[:,p.cross]
        radius2 = tang.copy()
        radius2[:,p.cross] = GAP**2+(cross if p.law=='central' else cross**2/PITCH_UNIT**2)
        radius = np.sqrt(radius2)
        energy, axial = morse(radius,p.depth,p.alpha,p.length)
        coeff = axial/radius
        coeff[:,p.cross] *= scale
        if p.law=='layered':
            coeff[:,p.cross] *= 2*cross/PITCH_UNIT**2
        bond_gradient = delta*coeff[:,:,None]
        nodal = (bond_gradient.transpose(0,2,1)@self.incidence).transpose(0,2,1)
        # For an interplane pair i<j, this is the +x physical force on lower i.
        tag_loads = bond_gradient[:,self.idx,0]
        d = cross.reshape(-1,COUNT,COUNT)
        mutual = ((d.argmin(axis=2)[:,self.a] == self.b-COUNT)
                  & (d.argmin(axis=1)[:,self.b-COUNT] == self.a))
        depth = -energy[:,self.idx]/p.depth[self.idx]
        return nodal,tag_loads,mutual & (depth > .6),mutual & (depth*scale > .6)


def experiment(args):
    prior = json.loads(NEUTRAL.read_text())
    archive = prior["saved_measurements"]["frames1321.json"]
    frame = archive["contents"][args.law]
    extraction = json.loads(args.cohorts.read_text())
    if extraction["archive_sha256"][str(NEUTRAL)] != digest(NEUTRAL):
        raise ValueError("Identity record refers to a different material preparation")
    cohort = extraction["laws"][args.law]
    identities = cohort["strict_surviving_ids"]
    tag = frame["sites"]
    which = identities.index(tag)
    if len({tuple(p) for p in identities}) != 9 or len(identities) != 9:
        raise ValueError("Wrong seed-1321 cohort")
    trace = WallTrace(args.mesh,.5)
    time,w,chi,dilation,gates = clock(trace,args.refine)
    k,mass,Z,station,mix,length = wall_line(trace,args.line,w,args.end_order)
    sites = Sites(args.law,identities)
    rest = np.array(frame['rested_xy'])
    if (tag != cohort['prior_tag'] or hashlib.sha256(rest.astype('<f8').tobytes()).hexdigest()
            != cohort['archive_rested_float64le_sha256']):
        raise ValueError("Neutral tag or rested coordinates changed")
    rest -= rest.mean(axis=0)
    seeds = SEEDS[args.first_history:args.first_history+args.histories]
    if len(seeds) != args.histories:
        raise ValueError("Requested histories outside the declared ensemble")
    S, D = len(seeds),len(DRIVES)
    B = S*D
    push = np.tile(DRIVES,S)
    pos = np.broadcast_to(rest,(B,)+rest.shape).copy()
    u = np.zeros((B,len(mass)))
    v = np.zeros_like(u)
    total_rim = np.zeros(B); total_port = np.zeros(B)
    total_end = np.zeros(B); total_bath = np.zeros(B)
    alive = np.ones((B,len(identities)),dtype=bool)
    absolute = alive.copy()
    first = np.full(alive.shape,np.nan)
    generators = [np.random.default_rng(s) for s in seeds]
    bridges = [np.random.default_rng(s+10000) for s in seeds]
    gate_lookup = {int(j):i for i,j in enumerate(gates)}
    tag_logs = []
    if args.history:
        store_xy = np.empty((len(time),S,2*COUNT,2))
        store_port = np.empty((len(time),D))
        store_rim = np.empty_like(store_port)

    emit(test="preparation",law=args.law,power=args.power,mesh=args.mesh,pieces=args.line,
         frame_sha256=archive["sha256"],cohort_sha256=digest(args.cohorts),pairs=identities,
         tag=tag,seeds=seeds,drive_order=DRIVES.tolist(),end_order=args.end_order,
         tether_k=args.tether,wall_arc_length=length,head=trace.head,intake=trace.intake,
         fields=trace.hashes,timing_gates=time[gates].tolist(),refine=args.refine)
    for j,t in enumerate(time):
        scale = float((chi[j]/chi[0])**args.power)
        grad, pair_force, now, now_abs = sites.read(pos,scale)
        newly_lost = alive & ~now
        first[newly_lost] = chi[j]
        alive &= now
        absolute &= now_abs
        a, b = station[j], mix[j]
        local = (1-b)*u[:,a]+b*u[:,a+1]
        spring = args.tether*(pos[:,:COUNT,0] - dilation[j]*rest[:COUNT,0] - local[:,None])
        port = spring.sum(axis=1)
        tension = np.diff(u,axis=1)*k
        rim = tension[:,0]
        terminal = Z*v[:,-1]
        if args.history:
            store_xy[j] = pos[::D]
            store_port[j] = port[:D]
            store_rim[j] = rim[:D]
        if j in gate_lookup:
            momentum = v@mass
            error = total_port-total_rim-total_end-momentum
            site_error = total_bath+total_port-push*t
            if max(np.max(abs(error)),np.max(abs(site_error))) > 3e-8:
                raise ArithmeticError("Wall or site impulse ledger failed")
            unforced = alive.reshape(S,D,-1)[:,0]
            continuous = alive.reshape(S,D,-1)
            event = dict(test="crossing",law=args.law,power=args.power,chi=float(chi[j]),w=float(w[j]),
                        t=float(t),dilation=float(dilation[j]),scale=scale,
                        strict_zero_ids=[np.flatnonzero(r).tolist() for r in unforced],
                        absolute_zero_count=absolute.reshape(S,D,-1)[:,0].sum(axis=1).tolist(),
                        loaded_tag_continuous=continuous[:,:,which].tolist(),
                        unforced_tag_current=now.reshape(S,D,-1)[:,0,which].tolist(),
                        max_wall_balance=float(max(abs(error))),max_site_balance=float(max(abs(site_error))))
            for name,arr in (("rim",rim),("rim_impulse",total_rim),("port",port),
                             ("port_impulse",total_port),("end_impulse",total_end),
                             ("line_momentum",momentum),("bath_impulse",total_bath),
                             ("tag_force",pair_force[:,which])):
                event[name] = arr.reshape(S,D).tolist()
            emit(**event)
        if j == len(time)-1:
            for r in range(S):
                tag_logs.append([None if np.isnan(x) else float(x) for x in first[r*D]])
            break

        dt = time[j+1]-t
        if args.refine==1 or j % 2 == 0:
            base_dt = dt if args.refine==1 else time[j+2]-t
            gaussian = np.stack([g.normal(size=rest.shape) for g in generators])
            wiener = np.sqrt(base_dt)*gaussian
            if args.refine==2:
                noise_bridge = np.sqrt(base_dt)*np.stack([g.normal(size=rest.shape) for g in bridges])
                leftover = (wiener-noise_bridge)/2
                wiener = (wiener+noise_bridge)/2
        else:
            wiener = leftover
        drift = -grad
        drift[:,:COUNT,0] -= spring
        drift[:,tag[1],0] += push
        random = np.repeat(wiener,D,axis=0)*np.sqrt(2*TEMPERATURE)
        affine = dilation[j+1]/dilation[j]-1
        change = drift*dt + random
        transported = (1+affine)*pos
        pos = transported + change
        # The net bath impulse is friction relative to the affine carrier,
        # less the random impulse supplied by that same bath (unit drag).
        total_bath += (pos[:,:,0]-transported[:,:,0]-random[:,:,0]).sum(axis=1)
        wall_force = np.zeros_like(u)
        wall_force[:,1:] -= tension
        wall_force[:,:-1] += tension
        wall_force[:,a] += (1-b)*port
        wall_force[:,a+1] += b*port
        wall_force[:,-1] -= terminal
        v[:,1:] += dt*wall_force[:,1:]/mass[1:]
        u += dt*v
        total_rim += dt*rim
        total_port += dt*port
        total_end += dt*terminal
    emit(test="first_loss_order",unforced_by_seed_and_original_pair=tag_logs,
         tag=tag, null_means_survived_to_end=True)
    if args.history:
        np.savez_compressed(args.history,t=time,w=w,chi=chi,xy_zero=store_xy,
                            station=station,mix=mix,k=k,mass=mass,Z=Z,port=store_port,rim=store_rim)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--law",choices=("central","layered"),required=True)
    p.add_argument("--power",type=int,choices=(0,1,2),required=True)
    p.add_argument("--cohorts",type=Path,required=True)
    p.add_argument("--mesh",type=int,choices=(34,48),default=48)
    p.add_argument("--line",type=int,default=80)
    p.add_argument("--refine",type=int,choices=(1,2),default=1)
    p.add_argument("--histories",type=int,choices=range(1,5),default=4)
    p.add_argument("--first-history",type=int,choices=range(4),default=0)
    p.add_argument("--end-order",type=float,default=.0001)
    p.add_argument("--tether",type=float,default=1.)
    p.add_argument("--history",type=Path)
    experiment(p.parse_args())
