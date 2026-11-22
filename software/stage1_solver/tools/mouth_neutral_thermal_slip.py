"""Tangential impulse of a dilute neutral Maxwell reservoir on a held periodic sheet.

Native units have m=k_B=1. The Cartesian engine traces 2 tangents and 1
normal. A complete incident crossing of z=5.05 is the denominator. The
zero-slip trace set also integrates the displaced Maxwell densities exactly
in velocity space by importance weighting (up to sampling quadrature).
"""
import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc, t


HEIGHT = 5.05
TEMPERATURE = 1.024
DENSITY = .001
MACH_DRAWS = (.125, -.125, .0625, -.0625)  # slip / one-component thermal deviation
SEEDS = tuple(range(61391, 61399))
CROSS_MORSE = (.6, 10., .94)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def emit(label, **data):
    print(json.dumps(dict(test=label, **data), allow_nan=False), flush=True)


def estimate(values):
    values = np.asarray(values, dtype=float)
    average = values.mean(axis=0)
    half = t.ppf(.975, len(values)-1)*values.std(ddof=1, axis=0)/np.sqrt(len(values))
    return dict(mean=average.tolist(), ci95=np.stack((average-half, average+half), axis=-1).tolist())


def prepare(seeds, temperature, drives, power):
    """Uniform primitive registry, Gaussian in-plane velocities, flux-Rayleigh normal speed."""
    rows, index = [], []
    for seed in seeds:
        points = qmc.Sobol(5, scramble=True, seed=seed).random_base2(power)
        for drift in drives:
            for k, (a,b,qx,qy,qn) in enumerate(points):
                s = math.sqrt(temperature)
                vx = s*(ndtri(qx)+drift)
                vy = s*ndtri(qy)
                vz = -s*math.sqrt(-2*math.log1p(-qn))
                index.append((seed, drift, k))
                rows.append((a,b,vx,vy,vz))
    return np.array(rows), np.array(index)


def simulate(args, label, chi, temperature, time, dt, power, drives, pitch, seeds=SEEDS):
    samples, index = prepare(seeds, temperature, drives, power)
    input_path = args.out/(label+'.input')
    raw = args.out/(label+'.jsonl')
    with input_path.open('x') as f:
        for i, (a,b,vx,vy,vz) in enumerate(samples):
            print(i, a,b,HEIGHT,vx,vy,vz,chi, file=f)
    steps = int(round(time/dt))
    command = [str(args.engine.resolve()), '--host', 'rigid', '--side', '1', '--pitch',
               str(pitch), '--mixed', '1', '--dt', str(dt), '--steps', str(steps),
               '--input', str(input_path.resolve())]
    with raw.open('x') as fd:
        subprocess.run(command, stdout=fd, check=True)
    rows = [json.loads(line) for line in raw.read_text().splitlines()]
    if len(rows) != len(samples) or any(row['id'] != k for k,row in enumerate(rows)):
        raise ArithmeticError('Incomplete Cartesian trajectories')
    def col(key):
        return np.array([row[key] for row in rows])
    x, y, z = col('vx_end'), col('vy_end'), col('vz_end')
    vx, vy, vz = samples[:,2:].T
    ix, iy = col('support_impulse_x'), col('support_impulse_y')
    refl = (col('z_end') >= HEIGHT) & (z > 0)
    tran = (col('z_end') <= -HEIGHT) & (z < 0)
    pending = ~(refl|tran)
    no_visit = pending & (col('entry20') < 0)
    crossed = col('plane_crossing_time') >= 0
    eerr = col('max_energy_error')
    emit('traces', case=label, temperature=temperature, chi=chi, native_duration=time,
         dt=dt, seeds=list(seeds), power=power, signed_slip_in_thermal_units=list(drives),
         traces=len(rows), command=command, input_sha256=sha(input_path), output_sha256=sha(raw),
         max_abs_energy_error=float(eerr.max()), min_incoming_energy=float(col('energy_tag_initial').min()),
         max_tangent_vector_momentum_residual=float(col('max_momentum_error').max()),
         pending=int(pending.sum()), reflected=int(refl.sum()), transmitted=int(tran.sum()),
         plane_crossed=int(crossed.sum()), pending_never_within_two_lengths=int(no_visit.sum()))
    if not np.all(np.isfinite(np.c_[x,y,z,ix,iy,eerr])) or np.any(col('energy_tag_initial') <= 0):
        raise ArithmeticError('Invalid trajectory energy or final state')
    # Each independent scramble yields a separate integral. The symmetry image
    # (-x,-y,vz) at inverted registry has (-Ix,-Iy,Iz). Its even observables
    # have the same value. Average that pair analytically only in the
    # reweighted estimator. The explicit moved launches use their own measured impulse.
    summary = {}
    for drift in drives:
        zero = drift==0
        cloud = []
        for seed in seeds:
            select = (index[:,0]==seed)&(index[:,1]==drift)
            if zero:
                per_speed = []
                for scaled in MACH_DRAWS:
                    u = scaled*np.sqrt(temperature)
                    odd_weight = np.exp(-scaled**2/2)*np.sinh(u*vx[select]/temperature)
                    even_weight = np.exp(-scaled**2/2)*np.cosh(u*vx[select]/temperature)
                    # If a pending trajectory later exits at the far plane, its
                    # asymptotic tangential speed is bounded by total incoming
                    # speed. The smooth potential at the read plane is negligible.
                    vin = np.linalg.norm(samples[select,2:], axis=1)
                    late_bound = np.mean(np.abs(odd_weight)*pending[select]*
                                         (abs(x[select])+vin))/abs(u)
                    per_speed.append([np.mean(odd_weight*ix[select])/u,
                                      np.mean(odd_weight*iy[select])/u,
                                      np.mean(even_weight*refl[select]), np.mean(even_weight*tran[select]),
                                      np.mean(even_weight*pending[select]),np.mean(even_weight), late_bound])
                score = np.mean(vx[select]*ix[select])/temperature
                reciprocity = np.mean(ix[select]**2)/(2*temperature)
                cloud.append((per_speed,score,reciprocity))
            else:
                cloud.append([np.mean(ix[select]),np.mean(iy[select])])
        key = str(drift)
        if zero:
            weighted, score, reciprocity = zip(*cloud)
            summary['weighted'] = np.asarray(weighted)
            emit('weighted', case=label, signed_slip_in_thermal_units=MACH_DRAWS,
                 components=['Ix/u','Iy/u','return_flux','transmitted_flux','pending_flux',
                             'integrated_probability','eventual_exit_impulse_slope_correction_bound'],
                 by_speed=estimate(weighted), zero_slip_slope=estimate(score),
                 squared_impulse_over_2T=estimate(reciprocity),
                 zero_signed_impulse_under_inversion=[0.,0.])
        else:
            summary[key] = np.asarray(cloud)
        parts = []
        for seed in seeds:
            sel=(index[:,0]==seed)&(index[:,1]==drift)
            # Measured incoming momentum partitions into returned and passing
            # outward momentum, still-present tracer momentum, and sheet impulse.
            parts.append([vx[sel].mean(),(x*refl)[sel].mean(),(x*tran)[sel].mean(),
                          (x*pending)[sel].mean(),ix[sel].mean(),iy[sel].mean(),
                          (x+ix-vx)[sel].mean(),(y+iy-vy)[sel].mean(),
                          (vz-z)[sel].mean(),refl[sel].mean(),tran[sel].mean(),pending[sel].mean()])
        emit('direct', case=label, signed_slip_in_thermal_units=drift,
             components=['incoming_x','reflected_x','transmitted_x','pending_x','support_x',
                         'support_y','ledger_residual_x','ledger_residual_y','tag_in_minus_out_normal',
                         'reflected_fraction','transmitted_fraction','pending_fraction'],
             reading=estimate(parts))
    for d in (.125,.0625):
        if d in drives and -d in drives:
            slope=(summary[str(d)]-summary[str(-d)])/(2*d*np.sqrt(temperature))
            emit('shifted_launch_slope',case=label, relative_slip=d, xy=estimate(slope))
    return summary['weighted']


def run(args):
    record=json.loads(args.pitch_record.read_text())
    source=record['guarded_runs']['primary']
    if not source['limits_verified'] or source['guard_exit'] or source['child_exit']:
        raise ValueError('Pitch source did not complete')
    sample, = [s for s in map(json.loads,source['literal_stdout'].splitlines()) if s['test']=='sample']
    pitch=float(sample['pitch'])
    # The discrete engine tapers pair energies between r=4.5 and 5. The exact
    # lateral mean of a radially tapered sum has Fz = 2 pi rho z V_taper(z).
    # It equals the engine's untapered mean force for z<=4.5. Bound the
    # difference in mean potential at any height by the missing far tail.
    depth, alpha, rest=CROSS_MORSE
    area=math.sqrt(3)*pitch*pitch/2
    tail=4*math.pi*depth/area*math.exp(-alpha*(4.5-rest))*(4.5/alpha+1/alpha**2)
    force_tail=4*math.pi*depth/area*4.5*math.exp(-alpha*(4.5-rest))
    emit('preparation', engine_sha256=sha(args.engine), pitch=pitch, pitch_record_sha256=sha(args.pitch_record),
         lateral_primitive_area=area, height=HEIGHT, mass=1., gas_number_density=DENSITY,
         half_space='incident Maxwell reservoir, positive z, slip along primitive vector 1',
         surface='periodic rigid triangular sheet, chi=1, or uniform lateral mean, chi=0',
         pair_morse=CROSS_MORSE, mean_potential_error_upper_bound=tail, mean_force_error_upper_bound=force_tail,
         temperatures=[dict(T=T, kinetic_pressure=DENSITY*T, incident_number_flux=DENSITY*np.sqrt(T/(2*np.pi)))
                       for T in (TEMPERATURE, TEMPERATURE/2)],
         unit_conversion='All forces and times are local microscopic native units. No mouth-derived coupling.')
    if args.independent_sampling:
        # Fresh scrambles and a larger phase-space sample check noisy differences
        # of directly shifted launches at both speeds without reusing trajectories.
        simulate(args,'ordered_repeat',1.,TEMPERATURE,48.,.004,8,(0.,*MACH_DRAWS),pitch,
                 seeds=range(72107,72123))
        emit('complete', cases=['ordered_repeat'])
        return
    results={}
    for name, chi in (('ordered',1.),('smooth',0.)):
        results[name]=simulate(args,name,chi,TEMPERATURE,48.,.004,6,(0.,*MACH_DRAWS),pitch)
    for label,time,dt,temp in [('long',96.,.004,TEMPERATURE),('fine',48.,.002,TEMPERATURE),
                               ('cool',48.,.004,TEMPERATURE/2)]:
        results[label]=simulate(args,label,1.,temp,time,dt,6,(0.,),pitch)
    for name in ('smooth','long','fine','cool'):
        delta=results[name][:,0,0]-results['ordered'][:,0,0]
        emit('paired_change', comparison=name+' minus ordered', Ix_over_u=estimate(delta))
    emit('complete', cases=list(results))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--engine',type=Path,required=True)
    p.add_argument('--pitch-record',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--independent-sampling',action='store_true')
    run(p.parse_args())
