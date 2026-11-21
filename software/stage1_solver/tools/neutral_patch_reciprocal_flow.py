"""Two-way local transverse momentum exchange during a matched 0PN transit.

Sites are mobile in both tangent coordinates. The solvent resolves the transverse
coordinate y and the layer normal z and averages along the drive direction x.
The same unit drag acts at every site. Only the remote upper edge sees the inlet
perturbation. Net site reaction is returned at its interpolation stencil.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_contact_tracking import WallTrace
from mouth_order_wave_pulse import digest, emit
from neutral_patch_companion_rim import CASES, driving, preparation
from neutral_patch_conversion import LEVELS, NEUTRAL, SEEDS, clock, wall_line
from neutral_patch_shear_channel import ShearChannel
from neutral_patch_sites import COUNT, TEMPERATURE

FIELD = Path('_scratch/s11c/companion_patch_fields_v1.npz')
PRIOR = Path('research/pde_ledger_v3/_measurements/neutral_patch_companion_rim.json')
OBSERVABLES = ('top_speed','fluid_lower','fluid_upper','material_lower','material_upper',
               'tag_fluid','tag_material','tag_drag','lower_drag_cv','lower_drag_cv2',
               'lower_surface_shear','lower_fluid_storage','lower_drag_mechanical','upper_drag_mechanical',
               'lower_drag_slip','grip','rim','strip_storage','distal')


def measurements(total, duration, seeds, cases=CASES):
    by_case = (total/duration).reshape(len(seeds),len(cases),len(OBSERVABLES))
    differential = by_case-by_case[:,:1]
    mean = differential.mean(axis=0)
    return {name:dict(zip(OBSERVABLES,values.tolist(),strict=True)) for name,values in zip(cases,mean,strict=True)}, differential


def noise_steps(times, refine, shape, seeds, case_count=len(CASES)):
    base = [np.random.default_rng(s) for s in seeds]
    fine = [np.random.default_rng(s+10000) for s in seeds]
    for j,dt in enumerate(np.diff(times)):
        if refine==1 or j%2==0:
            width = dt if refine==1 else times[j+2]-times[j]
            step = np.sqrt(width)*np.stack([g.normal(size=shape) for g in base])
            if refine==2:
                bridge = np.sqrt(width)*np.stack([g.normal(size=shape) for g in fine])
                other = (step-bridge)/2
                step = (step+bridge)/2
        else:
            step = other
        yield np.repeat(step,case_count,axis=0)*np.sqrt(2*TEMPERATURE)


def experiment(args, probe=None):
    tr = WallTrace(48,.5)
    time,w,chi,dilation,gates = clock(tr,args.refine)
    k,mass,Z,station,mix,_ = wall_line(tr,args.pieces,w,.0001)
    if probe is None:
        pinned = json.loads(PRIOR.read_text())['output_profiles_sha256'][str(FIELD)]
        if digest(FIELD)!=pinned:
            raise ValueError('Wrong companion input')
        speed = driving(tr,FIELD,w,.01)
        initial_speed = speed[0]
        cases = CASES
    else:
        pinned = None
        cases = probe.cases
        speed,initial_speed = probe.prepare(time)
    rest,sites,_ = preparation(args.law)
    tag = json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents'][args.law]['sites'][0]
    if not 0<=tag<COUNT:
        raise ValueError('Expected the archived wall-facing tag')
    seeds = SEEDS[:args.histories]
    copies = len(seeds)*len(cases)
    # All copies enter unstressed at the grips, at the same thermal site frame.
    pos = np.broadcast_to(rest,(copies,*rest.shape)).copy()
    strip = np.zeros((copies,len(mass)))
    strip_speed = np.zeros_like(strip)
    fluid = ShearChannel(args.grid,args.back_gap,args.backing,args.nu,1.,copies,
                         np.tile(initial_speed,len(seeds)),not args.cold,back_impedance=args.back_impedance)
    summed = np.zeros((copies,len(OBSERVABLES)))
    window = np.zeros_like(summed)
    gate_number = {int(g):i for i,g in enumerate(gates)}
    max_bulk_balance = 0.
    emit(test='setup',law=args.law,backing=args.backing,grid=args.grid,refine=args.refine,
         strip_pieces=args.pieces,back_gap=args.back_gap,back_impedance=args.back_impedance,
         nu=args.nu,rho=1.,gamma=1.,grip_k=1.,
         initialization=('specified_per_probe' if probe is not None else ('cold' if args.cold else 'steady_empty_channel')),seeds=seeds,tag=int(tag),
         length_unit_R=.01,field_sha256=pinned,fields=tr.hashes,intake=tr.intake,cases=cases,
         labels=OBSERVABLES,transit=time[-1],layer_order='lower wall-facing first, upper brane-facing second')
    for j,thermal in enumerate(noise_steps(time,args.refine,rest.shape,seeds,len(cases))):
        if probe is not None:
            probe.sample(j,fluid)
        duration = time[j+1]-time[j]
        top = np.tile(speed[j],len(seeds))
        local,stencil = fluid.sample(pos[:,:,1])
        mechanical = -sites.read(pos,chi[j]/chi[0])[0]
        a = station[j]
        b = mix[j]
        along = (1-b)*strip[:,a]+b*strip[:,a+1]
        grip = pos[:,:COUNT,0]-dilation[j]*rest[:COUNT,0]-along[:,None]
        mechanical[:,:COUNT,0] -= grip
        port = grip.sum(axis=1)
        # Peculiar material motion in the dilating patch. Thermal noise is
        # shared across paired and isolated copies and cancels in differences.
        drift = mechanical.copy()
        drift[:,:,0] += local
        transported = pos*(dilation[j+1]/dilation[j])
        pos = transported+(duration*drift+thermal)
        dx = (pos[:,:,0]-transported[:,:,0])
        lower = slice(None,COUNT)
        upper = slice(COUNT,None)
        # At unit site friction the deterministic internal/grip force is also
        # the net reaction on the carrier: friction and the Brownian reaction
        # sum to it for a massless site. Fluid momentum is evolved explicitly.
        cv1,cv2,bulk_error = fluid.step(duration,top,stencil,mechanical[:,:,0])
        max_bulk_balance = max(max_bulk_balance,float(np.max(abs(bulk_error))))
        tension = np.diff(strip,axis=1)*k
        rim = tension[:,0]
        end = Z*strip_speed[:,-1]
        load = np.zeros_like(strip)
        load[:,1:] -= tension
        load[:,:-1] += tension
        load[:,a] += (1-b)*port
        load[:,a+1] += b*port
        load[:,-1] -= end
        change_v = duration*load[:,1:]/mass[1:]
        strip_speed[:,1:] += change_v
        strip += duration*strip_speed
        imp = np.stack((duration*top,duration*local[:,lower].mean(axis=1),duration*local[:,upper].mean(axis=1),
                 dx[:,lower].mean(axis=1),dx[:,upper].mean(axis=1),duration*local[:,tag],dx[:,tag],
                 -duration*mechanical[:,tag,0],cv1,cv2,fluid.lower_shear_impulse,fluid.lower_storage_change,
                 -duration*mechanical[:,lower,0].sum(axis=1),-duration*mechanical[:,upper,0].sum(axis=1),
                 (duration*local[:,lower]-dx[:,lower]+thermal[:,lower,0]).sum(axis=1),duration*port,
                 duration*rim,change_v@mass[1:],duration*end),axis=1)
        window += imp
        summed += imp
        if j+1 in gate_number:
            i = gate_number[j+1]
            elapsed = time[j+1]-time[gates[i-1]]
            rows,_ = measurements(window,elapsed,seeds,cases)
            emit(test='interval',order=[LEVELS[i-1],LEVELS[i]],duration=elapsed,**{cases[1]:rows[cases[1]]})
            window.fill(0.)
    rows,diff = measurements(summed,time[-1]-time[0],seeds,cases)
    absolute = {key:arr for key,arr in zip(OBSERVABLES,summed.T,strict=True)}
    target = absolute['lower_drag_mechanical']
    checks = {key:float(np.max(abs(absolute[key]-target))) for key in ('lower_drag_cv','lower_drag_cv2','lower_drag_slip')}
    checks['strip'] = float(max(abs(absolute['grip']-absolute['rim']-absolute['strip_storage']-absolute['distal'])))
    checks['raw_step_fluid_momentum'] = max_bulk_balance
    checks['site_pair_cancellation'] = float(max(abs(target+absolute['upper_drag_mechanical']-absolute['grip'])))
    suffix = '8' if probe is None else 'pulse'
    emit(test='result',cases=rows,ledger_errors=checks,
         **{'seed_lower_force_'+suffix:diff[:,1,OBSERVABLES.index('lower_drag_cv')].tolist(),
            'seed_rim_force_'+suffix:diff[:,1,OBSERVABLES.index('rim')].tolist()})
    if max(checks.values())>2e-9 or not np.isfinite(summed).all():
        raise ArithmeticError('Fluid, site or strip momentum balance failed')
    if probe is not None:
        probe.sample(len(time)-1,fluid)
        probe.finish(time,diff)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--law',choices=('central','layered'),required=True)
    p.add_argument('--backing',choices=('rest','free','robin'),required=True)
    p.add_argument('--back-impedance',type=float)
    p.add_argument('--grid',type=int,choices=(4,8),default=8)
    p.add_argument('--refine',type=int,choices=(1,2),default=2)
    p.add_argument('--pieces',type=int,default=160)
    p.add_argument('--back-gap',type=float,default=1.5)
    p.add_argument('--nu',type=float,default=1.)
    p.add_argument('--histories',type=int,choices=(1,2,3,4),default=4)
    p.add_argument('--cold',action='store_true')
    experiment(p.parse_args())
