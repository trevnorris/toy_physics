"""Open a converted patch to a finite material bath across its free normal face."""
import argparse
from copy import copy
from pathlib import Path

import numpy as np

from mouth_brace_grip import Material, clock, evolve_interval
from mouth_brace_bulk_shear import Strained, shear_work
from mouth_middle_tag_grip import cubic_trace
from neutral_patch_open_w import mix
from neutral_patch_sites import BETWEEN, TEMPERATURE, WITHIN, emit, morse
from neutral_patch_structural_endpoints import PITCH


def converted(mat, seed, times, order):
    pos = mat.flat.copy()
    perturb = np.random.default_rng(seed + 11).normal(size=pos.shape)
    pos += .001 * (perturb - perturb.mean(axis=0))
    rng = np.random.default_rng(seed)
    bridge = np.random.default_rng(seed + 7013)
    for start, end in zip(times[:-1], times[1:]):
        evolve_interval(mat, [pos], start, end, [order], rng, bridge, 1)
    return pos - pos.mean(axis=0)


def double_images(i, j, shifts, n, joined):
    ii, jj, offset = [i, i+n], [j, j+n], [shifts, shifts]
    if joined:
        a, b = np.indices((n, n)).reshape(2, -1)
        tiles = np.unique(shifts, axis=0)
        ii.append(np.tile(a, len(tiles)))
        jj.append(np.tile(b+n, len(tiles)))
        offset.append(np.repeat(tiles, len(a), axis=0))
    return (np.ascontiguousarray(np.concatenate(ii), dtype=np.int32),
            np.ascontiguousarray(np.concatenate(jj), dtype=np.int32), np.vstack(offset))


def material_pair(mat, joined):
    pair = copy(mat)
    pair.cell = copy(mat.cell)
    pair.n = pair.cell.n = 2*mat.n
    pair.load = pair.height = np.zeros(pair.n)  # fully converted, free external faces
    pair.cell.i, pair.cell.j, pair.cell.shifts = double_images(mat.cell.i, mat.cell.j, mat.cell.shifts, mat.n, joined)
    i,j,pair.shifts = double_images(*mat.indices, mat.shifts, mat.n, joined)
    pair.indices = [i,j]
    return pair


def observations(cell, pos, check=False):
    # Cut follows overall translation. Split each bond's tangential virial by
    # its straight segment's length fraction below the plane (equal for images).
    d = pos[cell.j] - pos[cell.i]
    d[:,:2] += cell.shifts - cell.box*np.rint(d[:,:2]/cell.box)
    r = np.linalg.norm(d, axis=1)
    keep = r < cell.cutoff
    i,j,d,r = cell.i[keep],cell.j[keep],d[keep],r[keep]
    f, df, _ = mix(d[:,2], cell.power)
    ea, ga = morse(r, *WITHIN)
    eb, gb = morse(r, *BETWEEN)
    radial = (ga+f*(gb-ga))/r
    virial = radial*(d[:,0]**2 - d[:,1]**2)/cell.area
    z = pos[:,2]-pos[:,2].mean()
    za,zb = z[i],z[j]
    crosses = (za<0) != (zb<0)
    low = (za<0).astype(float)
    low[crosses] = -np.minimum(za[crosses],zb[crosses])/abs(d[crosses,2])
    # Force on the lower material exerted by upper material, positive in tension.
    traction = np.sum((radial*d[:,2] + df*(eb-ea))[crosses]*np.sign(d[crosses,2]))/cell.area
    stress = virial.sum()
    if check:
        h = 2e-7
        numerical = (shear_work(cell,pos,h)[0]-shear_work(cell,pos,-h)[0])/(2*h*cell.area)
        np.testing.assert_allclose(numerical, stress, rtol=3e-6, atol=2e-5)
    return np.array([stress, np.dot(low,virial)]), traction


def pair_distances(pos, box):
    d = pos[:,None,:] - pos[None,:,:]
    d[:,:,:2] -= box*np.rint(d[:,:,:2]/box)
    return np.linalg.norm(d,axis=2)


class Tracking:
    """Persistent exchange through a band, and neighbor identities with two cutoffs."""
    def __init__(self, pos, box, joined):
        self.n = len(pos)//2
        self.possible = np.triu(np.ones((len(pos),len(pos)),dtype=bool),1)
        if not joined:
            self.possible[:self.n,self.n:] = False
        self.bonds = (pair_distances(pos,box)<1.15) & self.possible
        self.initial = self.bonds.copy()
        self.last = np.sign(pos[:,2]-pos[:,2].mean()).astype(int)
        self.origin = self.last.copy()
        self.ever = np.zeros(len(pos),dtype=bool)
        self.from_lower = np.arange(len(pos)) < self.n
        self.joined = joined
        self.peak_crossers = 0

    def update(self,pos,box):
        r = pair_distances(pos,box)
        self.bonds = ((self.bonds & (r<1.30)) | (r<1.15)) & self.possible
        z = pos[:,2]-pos[:,2].mean()
        side = np.where(z<-.20,-1,np.where(z>.20,1,self.last))
        if self.joined:
            self.ever |= (side != self.origin)
            self.peak_crossers = max(self.peak_crossers,int(np.sum(side!=self.origin)))
        self.last = side
        return dict(original_bonds=int(self.initial.sum()),lost_original=int((self.initial & ~self.bonds).sum()),
                    new_bonds=int((self.bonds & ~self.initial).sum()),cross_origin_bonds=int(self.bonds[:self.n,self.n:].sum()),
                    outside_starting_half_now=int((side!=self.origin).sum()) if self.joined else None,
                    distinct_reaching_other_half=int(self.ever.sum()) if self.joined else None,
                    peak_simultaneous=int(self.peak_crossers) if self.joined else None,
                    lower_origin_in_upper=int((self.from_lower & (z>.20)).sum()) if self.joined else None,
                    upper_origin_in_lower=int((~self.from_lower & (z<-.20)).sum()) if self.joined else None,
                    lower_population=int((z<0).sum()),height_span=float(np.ptp(z)))


def run(args):
    args.work.mkdir(parents=True,exist_ok=True)
    args.mesh,args.cutoff,args.side,args.brace,args.mobility,args.temperature = 34,4.5,4,'elastic',1.,TEMPERATURE
    depths,times,order,intake = clock(args.mesh)
    trace,_ = cubic_trace(args.mesh)
    transit = trace.transit(depths[-1],trace.edges[-1],8)
    mat = Material(args)
    bottom = converted(mat,args.seed,times,order)
    top = converted(mat,args.partner,times,order)
    # Start with a full pitch of clear normal space between the closest faces;
    # a free normal translation and an in-plane registration are preparations.
    top[:,2] += bottom[:,2].max() - top[:,2].min() + PITCH
    top[:,:2] += np.array([.5, np.sqrt(3)/6])*PITCH*args.registration
    initial = np.vstack([bottom,top]); initial[:,2] -= initial[:,2].mean()
    rng = np.random.default_rng(args.seed+23*args.partner)
    materials = [material_pair(mat,True),material_pair(mat,False)]
    zero = np.zeros((1,len(initial),3))
    for material in materials:
        trial = initial.copy()
        material.step(trial,1e-6,np.zeros(1),zero)
        force = -material.cell.evaluate(initial)[1]
        np.testing.assert_allclose((trial-initial)/1e-6,force,atol=1.e-8)
    e0 = sum(shear_work(mat.cell,p,0.)[0] for p in (bottom,top))
    s0 = sum(shear_work(mat.cell,p,0.)[1] for p in (bottom,top))
    np.testing.assert_allclose(shear_work(materials[1].cell,initial,0.),[e0,s0],atol=1e-10)
    emit(test='setup',seeds=[args.seed,args.partner],mix_power=args.power,registration=args.registration,
         samples=[mat.n,mat.n],area=mat.cell.area,temperature=mat.temperature,applied_outer_normal_traction=0.,
         applied_normal_brace=0.,intake=intake,transit=transit,prep_transits=1,hold_transits=1,
         committed_crossing_half_band=.2,projected_sites_per_patch=mat.n/mat.cell.area,dt_max=args.dt)
    phases = ['preparation','hold']
    positions = [initial.copy() for _ in materials]
    nread = 48
    for phase in phases:
        amplitudes = [0.] if phase=='preparation' else [0.,1e-4,-1e-4,5e-5,-5e-5]
        clones = [[Strained(m,a,p) for a in amplitudes] for m,p in zip(materials,positions)]
        trackers = [Tracking(p,m.cell.box,index==0) for index,(p,m) in enumerate(zip(positions,materials))]
        for frame in range(nread+1):
            if frame:
                steps = int(np.ceil(transit/nread/args.dt))
                dt = transit/nread/steps
                for start in range(0,steps,128):
                    count = min(steps-start,128)
                    noise = rng.normal(size=(count,len(initial),3))*np.sqrt(2*mat.temperature*dt)
                    orders = np.zeros(count)
                    for m,states,tracker in zip(materials,clones,trackers):
                        for state in states:
                            state.select()
                            m.step(state.pos,dt,orders,noise)
                        tracker.update(states[0].pos,states[0].box)
            for name,m,states,tracker in zip(['joined','disconnected'],materials,clones,trackers):
                measured = []
                for state in states:
                    state.select()
                    measured.append(observations(m.cell,state.pos,frame in (0,nread)))
                pos = states[0].pos
                samples = np.array([item[0] for item in measured])
                diag = tracker.update(pos,states[0].box)
                data = dict(test=phase,contact=name,fraction=frame/nread,bond_normal_traction=measured[0][1],
                            rest_stress=samples[0].tolist(),**diag)
                if phase=='hold':
                    data.update(full_derivative=((samples[1]-samples[2])/2e-4).tolist(),
                                half_derivative=((samples[3]-samples[4])/1e-4).tolist())
                emit(**data)
        positions = [states[0].pos for states in clones]
        for states in clones:
            states[0].select()


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--power',type=int,choices=(1,2),required=True)
    parser.add_argument('--seed',type=int,default=319)
    parser.add_argument('--partner',type=int,default=817)
    parser.add_argument('--registration',type=int,choices=(0,1),default=1)
    parser.add_argument('--dt',type=float,default=.001)
    run(parser.parse_args())
