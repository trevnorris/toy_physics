"""Hold one registered pair through the planar conversion band and read its zero load."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq
from scipy.special import expit, logit

from mouth_conversion_wave_patch import KAPPA, MOBILITY, settings
from neutral_patch_conversion import NEUTRAL, Sites
from neutral_patch_sites import COUNT, GAP, PITCH_UNIT, STEP, TEMPERATURE, Patch, emit, morse


def order_clock(refine):
    setup, _ = settings()
    hi, lo = setup['ordered_side'], setup['less_ordered_side']
    speed, width = setup['material_speed'], setup['width']
    length, head = setup['length'], 1.01*setup['head_reference']
    mesh = np.linspace(0., length, 65)
    c = expit(logit(hi)-mesh/width)
    def rhs(z, y):
        c, slope = y
        return [slope,(speed*slope/MOBILITY + c*(1-c)*(
            KAPPA/width**2*(1-2*c)+6*head))/KAPPA]
    front = solve_bvp(rhs,lambda a,b: [a[0]-hi,b[0]-lo],mesh,[c,-c*(1-c)/width],tol=1e-9)
    if not front.success:
        raise ArithmeticError(front.message)
    # Eight equal-age reading intervals. Each coarse step splits in two in the refinement.
    count = 8*refine*int(np.ceil(length/(8*speed*STEP)))
    t = np.linspace(0.,length/speed,count+1)
    chi = front.sol(speed*t)[0]
    assert np.all(np.diff(chi)<0.)
    return t, chi


def shifted(pos, tag, slip):
    result = pos.copy()
    i,j = tag
    delta = slip-(pos[j,0]-pos[i,0])
    result[i,0] -= delta/2
    result[j,0] += delta/2
    return result


def root_reading(patch, pos, tag, meter):
    i,j = tag
    held = pos[j,0]-pos[i,0]
    pair = (patch.i==i) & (patch.j==j)
    def reaction(s):
        points = shifted(pos,tag,s)
        if meter=='bond':
            _, grad, _ = patch.evaluate(points,removed=~pair)
        else:
            _, grad, _ = patch.evaluate(points)
        return (grad[j,0]-grad[i,0])/2
    # Bracket the incremental zero on both sides, including a flat cubic zero.
    # Report an unresolved local readout if a warm configuration loses the bracket.
    center = 0. if meter=='bond' else held
    for span in (.015,.06,.18,.35):
        grid = center+span*np.linspace(-1.,1.,10)
        force = np.array([reaction(s) for s in grid])
        crossings = np.flatnonzero((force[:-1]<=0)&(force[1:]>0))
        if len(crossings)==1:
            k=crossings[0]
            zero=brentq(reaction,grid[k],grid[k+1],xtol=1e-13)
            return float(zero),float(held-zero)
        if len(crossings)>1:
            break
    return None,None


def cross_energy(sites, positions):
    p = sites.patch
    diff = positions[:,p.j[p.cross]]-positions[:,p.i[p.cross]]
    lateral = np.sum(diff**2,axis=-1)
    effective = lateral if p.law=='central' else lateral**2/PITCH_UNIT**2
    raw, _ = morse(np.sqrt(GAP**2+effective),p.depth[p.cross],p.alpha[p.cross],p.length[p.cross])
    return raw.sum(axis=1)


def experiment(args, meter, times, orders, frame, tag):
    # Re-equilibrate the same archived constituents at the band's entrance.
    # The previous conversion law normalizes cross-layer well depth at chi=.95.
    probe = Patch(args.law)
    cross_depth = probe.depth[probe.cross].copy()
    scale = orders/.95
    probe.depth[probe.cross] = cross_depth*scale[0]
    rest = np.asarray(frame['rested_xy'])
    i,j = tag
    if meter=='bond':
        rest = rest.copy()
        rest[[i,j]] = rest[[i,j]].mean(axis=0)
        rest,_ = probe.relaxed(rest,fixed=tag)
    else:
        rest,_ = probe.relaxed(rest)
    origin,_ = root_reading(probe,rest,tag,meter)
    if origin is None:
        raise ArithmeticError('Entrance zero-load state unresolved')
    prestrains = np.array([0.,.002,-.002,.001,-.001])
    pos = np.array([shifted(rest,tag,origin+eps) for eps in prestrains])
    entry_energy = np.array([probe.evaluate(x)[0].sum() for x in pos])
    sites = Sites(args.law,[tag])
    rng = np.random.default_rng(args.seed)
    split = np.random.default_rng(args.seed+10000)
    dt = times[1]
    temperature = TEMPERATURE if args.warm else 0.
    changing = args.path=='conversion'
    if not changing:
        scale = np.full_like(scale,scale[0])
    positive_work = np.zeros(len(pos))
    dissipation = np.zeros(len(pos))
    alive = np.ones(len(pos),dtype=bool)
    entry_clamp = pos[:,[i,j]].copy()
    before = cross_energy(sites,pos)
    for k, s in enumerate(scale):
        grad, _, live, _ = sites.read(pos,s)
        alive &= live[:,0]
        grad[:,[i,j]] = 0.
        if k%(int((len(times)-1)//8))==0:
            probe.depth[probe.cross] = cross_depth*s
            readings = [root_reading(probe,x,tag,meter) for x in pos]
            zeros,extensions = zip(*readings)
            def gain(a,b):
                if not all(alive[[a,b]]) or extensions[a] is None or extensions[b] is None:
                    return None
                return (extensions[a]-extensions[b])/(prestrains[a]-prestrains[b])
            end_energy = np.array([probe.evaluate(x)[0].sum() for x in pos])
            energy_change = end_energy-entry_energy
            emit(test='tag', law=args.law, meter=meter, path=args.path, tag=tag,
                 temperature=temperature, seed=args.seed, refinement=args.refine,
                 time=times[k], chi=.95*s, held_slips=(origin+prestrains).tolist(),
                 input_prestrains=prestrains.tolist(), zero_force_positions=zeros,
                 elastic_extensions=extensions, registered_at_all_steps=alive.tolist(),
                 full_odd_ratio=gain(1,2),half_odd_ratio=gain(3,4),
                 order_work_into_patch_per_site=(positive_work/(2*COUNT)).tolist(),
                 delta_material_energy_per_site=(energy_change/(2*COUNT)).tolist(),
                 heat_to_bath_per_site=((positive_work-energy_change)/(2*COUNT)).tolist(),
                 cold_drag_integral_per_site=None if temperature else (dissipation/(2*COUNT)).tolist())
            np.testing.assert_array_equal(pos[:,[i,j]],entry_clamp)
            if not temperature:
                np.testing.assert_allclose(positive_work-energy_change,dissipation,atol=1e-8,rtol=.001)
            if k==len(times)-1:
                break
        if temperature:
            if k%args.refine==0:
                coarse = np.sqrt(2*temperature*dt*args.refine)*rng.normal(size=(2*COUNT,2))
                if args.refine==2:
                    bridge = np.sqrt(4*temperature*dt)*split.normal(size=coarse.shape)
                    noise,remaining = (coarse+bridge)/2,(coarse-bridge)/2
                else:
                    noise = coarse
            else:
                noise = remaining
            noise[[i,j]] = 0.
        else:
            noise = 0.
        trial = pos+noise-dt*grad
        other, _, _, _ = sites.read(trial,scale[k+1])
        other[:,[i,j]] = 0.
        new = pos+noise-.5*dt*(grad+other)
        after = cross_energy(sites,new)
        positive_work += (scale[k+1]-s)*(before+after)/2
        if not temperature:
            dissipation += .5*dt*(np.sum(grad**2,axis=(1,2))+np.sum(other**2,axis=(1,2)))
        pos,before = new,after


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--law',choices=('central','layered'),required=True)
    parser.add_argument('--path',choices=('conversion','fixed'),required=True)
    parser.add_argument('--refine',choices=(1,2),type=int,required=True)
    parser.add_argument('--warm',action='store_true')
    parser.add_argument('--seed',type=int,default=8401)
    args = parser.parse_args()
    frames = json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']
    cohort = json.loads(Path('research/pde_ledger_v3/_measurements/neutral_patch_conversion.json').read_text())
    frame = frames[args.law]
    tag = frame['sites']
    assert tag in cohort['independent_id_extraction']['laws'][args.law]['strict_surviving_ids']
    times,chi = order_clock(args.refine)
    for meter in ('bond','clamp'):
        experiment(args,meter,times,chi,frame,tag)
