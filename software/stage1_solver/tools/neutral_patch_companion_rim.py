"""Companion and isolated shear-bath experiments on the same converting patch.

The six copies share the microscopic initial state and every thermal impulse.
The lateral carrier trace changes the bath velocity in the overdamped material
equation. The held waist reads the elastic tension at the end of a moving-load
line. Distances, time, and force are in the units specified in the preparation.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_contact_tracking import WallTrace
from mouth_order_wave_pulse import digest, emit
from neutral_patch_conversion import NEUTRAL, SEEDS, Sites, clock, wall_line
from neutral_patch_sites import COUNT, TEMPERATURE

PRIOR = Path('research/pde_ledger_v3/_measurements/neutral_patch_conversion.json')
CASES = ('alone', 'toward8', 'mirror8', 'toward16', 'mirror16', 'half8')


def preparation(law):
    archive = json.loads(PRIOR.read_text())
    neutral = json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents'][law]
    controls = {}
    for suffix in ('fine_01', 'fine_rest_01'):
        run = archive['verified_runs'][f'converting_{law}_{suffix}']
        if (run['guard_exit'] or run['child_exit'] or not run['limits_verified']):
            raise ValueError('Unverified isolated control')
        lines = list(map(json.loads, run['literal_stdout']))
        header = lines[0]
        for i, seed in enumerate(header['seeds']):
            controls[seed] = [event['rim_impulse'][i][0] for event in lines if event['test']=='crossing']
    if header['frame_sha256'] != json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['sha256']:
        raise ValueError('Neutral frame archive mismatch')
    rest = np.array(neutral['rested_xy'])
    return rest-rest.mean(axis=0), Sites(law,header['pairs']), controls


def driving(trace, path, w, length_ratio):
    with np.load(path,allow_pickle=False) as source:
        np.testing.assert_array_equal(source['wall_edges'], trace.edges)
        expected_w = trace.geometry.points[trace.nodes,1]
        np.testing.assert_allclose(source['w_nodes'],expected_w,rtol=0,atol=1e-10)
        speeds = [trace.sample(source[name],w)[0] for name in ('velocity8','velocity16')]
    if length_ratio <= 0 or not np.all(np.isfinite(speeds)):
        raise ValueError('Invalid length scale or velocity')
    u8,u16 = speeds
    return np.column_stack([np.zeros_like(w),u8,-u8,u16,-u16,u8/2])/length_ratio


def run(args, observe_step=None):
    trace = WallTrace(48,.5)
    t,w,order,dilation,gates = clock(trace,args.refine)
    stiffness,mass,termination,station,weight,_ = wall_line(trace,args.pieces,w,.0001)
    speed = driving(trace,args.profile,w,args.length_ratio)
    rest,sites,controls = preparation(args.law)
    seeds = SEEDS[:args.histories]
    cases = len(CASES)
    copies = len(seeds)*cases
    pos = np.broadcast_to(rest,(copies,*rest.shape)).copy()
    displ = np.zeros((copies,len(mass)))
    wall_speed = np.zeros_like(displ)
    rim_impulse,port_impulse,end_impulse,bath_impulse = np.zeros((4,copies))
    generators = [np.random.default_rng(s) for s in seeds]
    bridges = [np.random.default_rng(s+10000) for s in seeds]
    exposed = slice(COUNT,None) if args.exposure=='upper' else slice(None)
    exposed_count = COUNT if args.exposure=='upper' else 2*COUNT
    gate_number = {int(j):i for i,j in enumerate(gates)}
    emit(test='patch_setup',law=args.law,seeds=seeds,cases=CASES,intake=trace.intake,
         fields=trace.hashes,flow_profile_sha256=digest(args.profile),order_head=trace.head,
         histories_source_sha256=digest(PRIOR),length_per_patch_unit_R=args.length_ratio,
         patch_force_unit='within-plane well energy / preferred same-plane pair distance',
         exposure=args.exposure,grip=args.grip,pieces=args.pieces,refine=args.refine,
         entry_velocity_patch_units=speed[0].tolist(),local_span_in_R=(np.ptp(rest,axis=0)*args.length_ratio).tolist())
    for j,time in enumerate(t):
        gradient = sites.read(pos,order[j]/order[0])[0]
        a,b = station[j],weight[j]
        wall_local = (1-b)*displ[:,a]+b*displ[:,a+1]
        grip = args.grip*(pos[:,:COUNT,0]-dilation[j]*rest[None,:COUNT,0]-wall_local[:,None])
        port = grip.sum(axis=1)
        tension = np.diff(displ,axis=1)*stiffness
        rim = tension[:,0]
        terminal = termination*wall_speed[:,-1]
        if j in gate_number:
            momentum = wall_speed@mass
            wall_error = np.max(abs(port_impulse-rim_impulse-end_impulse-momentum))
            bath_error = np.max(abs(bath_impulse-port_impulse))
            signals = {}
            for name,values in [('rim_force',rim),('rim_impulse',rim_impulse),
                                ('grip_impulse',port_impulse),('strip_momentum',momentum),
                                ('distal_impulse',end_impulse),('bath_on_sites',bath_impulse)]:
                signals[name] = values.reshape(-1,cases).tolist()
            index = gate_number[j]
            reference_error = None
            if (args.refine,args.pieces,args.grip)==(2,160,1.):
                reference = [controls[s][index] for s in seeds]
                reference_error = float(np.max(abs(np.array(signals['rim_impulse'])[:,0]-reference)))
            emit(test='patch_gate',index=index,t=float(time),w=float(w[j]),order=float(order[j]),
                 carrier_velocity_R_per_time=(speed[j]*args.length_ratio).tolist(),
                 wall_balance_max=float(wall_error),bath_balance_max=float(bath_error),
                 prior_isolated_impulse_max_error=reference_error,**signals)
            if max(wall_error,bath_error,reference_error or 0.)>3e-9:
                raise ArithmeticError('Isolated replay or impulse balance failed')
            if j==len(t)-1:
                break
        dt = t[j+1]-time
        if args.refine==1 or j%2==0:
            noise_dt = dt if args.refine==1 else t[j+2]-time
            noise = np.sqrt(noise_dt)*np.stack([g.normal(size=rest.shape) for g in generators])
            if args.refine==2:
                bridge = np.sqrt(noise_dt)*np.stack([g.normal(size=rest.shape) for g in bridges])
                second_half = (noise-bridge)/2
                noise = (noise+bridge)/2
        else:
            noise = second_half
        thermal = np.repeat(noise,cases,axis=0)*np.sqrt(2*TEMPERATURE)
        local_bath = np.tile(speed[j],len(seeds))
        drift = -gradient
        drift[:,:COUNT,0] -= grip
        drift[:,exposed,0] += local_bath[:,None]
        transported = (1+(dilation[j+1]/dilation[j]-1))*pos
        pos = transported+(dt*drift+thermal)
        if observe_step is not None:
            observe_step(dt, gate_number.get(j+1), local_bath, pos[:,:,0]-transported[:,:,0],
                         gradient[:,:,0], port)
        # Actual bath reaction after advection and thermal impulses (unit drag).
        bath_impulse += (dt*exposed_count*local_bath
                         -(pos[:,:,0]-transported[:,:,0]-thermal[:,:,0]).sum(axis=1))
        net = np.zeros_like(displ)
        net[:,1:] -= tension
        net[:,:-1] += tension
        net[:,a] += (1-b)*port
        net[:,a+1] += b*port
        net[:,-1] -= terminal
        wall_speed[:,1:] += dt*net[:,1:]/mass[1:]
        displ += dt*wall_speed
        rim_impulse += dt*rim
        port_impulse += dt*port
        end_impulse += dt*terminal

    elapsed = t[-1]-t[0]
    mean_force = np.array(signals['rim_impulse'])/elapsed
    differential = mean_force-mean_force[:,:1]
    delta = differential.mean(axis=0)
    record = dict(test='patch_result',cases=CASES,delta_force_by_seed=differential.tolist(),
                  absolute_mean_force=mean_force.mean(axis=0).tolist(),delta_mean_force=delta.tolist(),
                  min_seed_delta=differential.min(axis=0).tolist(),max_seed_delta=differential.max(axis=0).tolist(),
                  elapsed=float(elapsed),average='time integral over one conversion transit, equal seed weights')
    if delta[1]!=0.:
        record.update(ratio16_to_8=delta[3]/delta[1],double_half8_to_full=2*delta[5]/delta[1],
                      mirror_sum_over_full=(delta[2]+delta[1])/delta[1])
    emit(**record)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile',required=True,type=Path)
    parser.add_argument('--law',required=True,choices=('central','layered'))
    parser.add_argument('--exposure',choices=('upper','both'),default='upper')
    parser.add_argument('--length-ratio',type=float,default=.01)
    parser.add_argument('--grip',type=float,default=1.)
    parser.add_argument('--pieces',type=int,default=160)
    parser.add_argument('--refine',type=int,choices=(1,2),default=2)
    parser.add_argument('--histories',type=int,choices=(1,2,3,4),default=4)
    run(parser.parse_args())
