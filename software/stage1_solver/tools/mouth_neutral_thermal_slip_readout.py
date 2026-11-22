"""Reduce the saved Cartesian neutral-interface branch ledgers, in incident-gas units."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import t


def estimate(x):
    x=np.asarray(x,dtype=float)
    m=x.mean(axis=0)
    h=t.ppf(.975,len(x)-1)*x.std(axis=0,ddof=1)/math.sqrt(len(x))
    return {'mean':m.tolist(),'ci95':np.stack((m-h,m+h),axis=-1).tolist()}


def main(path):
    rows=list(map(json.loads,(path/'guard/stdout').read_text().splitlines()))
    setup,=[v for v in rows if v['test']=='preparation']
    finished,=[v for v in rows if v['test']=='complete']
    for name in finished['cases']:
        case,=[v for v in rows if v['test']=='traces' and v['case']==name]
        temperature=case['temperature']
        speed=math.sqrt(temperature)
        flux=setup['gas_number_density']*speed/math.sqrt(2*math.pi)
        pressure=setup['gas_number_density']*temperature
        data=list(map(json.loads,(path/(name+'.jsonl')).read_text().splitlines()))
        n=2**case['power']
        slips=case['signed_slip_in_thermal_units']
        nseed=len(case['seeds'])
        def col(key):
            return np.array([d[key] for d in data]).reshape(nseed,len(slips),n)
        vx, final_vx, vz, final_vz=map(col,['vx0','vx_end','vz0','vz_end'])
        impulse=col('support_impulse_x')
        equilibrium=slips.index(0.)
        x0=vx[:,equilibrium,:]
        i0=impulse[:,equilibrium,:]
        # From the actual sampling weights of the two Maxwell populations.
        # Inverting velocities and registry pairs removes the even impulse.
        u=.125*speed
        w=np.sinh(u*x0/temperature)*math.exp(-u*u/(2*temperature))
        full=(w*i0).mean(axis=1)/u
        half=(w[:,:n//2]*i0[:,:n//2]).mean(axis=1)/u
        normal=(final_vz[:,equilibrium,:]-vz[:,equilibrium,:]).mean(axis=1)*flux/pressure
        results=dict(test='traction',case=name,temperature=temperature,kinetic_pressure=pressure,
                     incident_number_flux=flux,weak_slip=u,zero_slip_normal_load_over_pressure=estimate(normal),
                     support_impulse_per_incident_tag=estimate(u*full),shear_traction=estimate(u*flux*full),
                     shear_traction_over_kinetic_pressure=estimate(u*flux*full/pressure),
                     drag_coefficient_per_area=estimate(flux*full),slope_half_sample=estimate(half),
                     paired_slope_change_full_minus_half_sample=estimate(full-half))
        if .125 in slips:
            for mag in (.125,.0625):
                d=mag*speed
                wp=math.exp(-mag*mag/2)*np.sinh(d*x0/temperature)
                weighted=(wp*i0).mean(axis=1)/d
                shifted=(impulse[:,slips.index(mag),:]-impulse[:,slips.index(-mag),:]).mean(axis=1)/(2*d)
                results[f'shifted_minus_weighted_slope_{mag}']=estimate(shifted-weighted)
        if not np.allclose(vx-final_vx,impulse,rtol=0.,atol=1.e-10):
            raise ArithmeticError('Tangential momentum reduction disagrees with support meter')
        print(json.dumps(results,allow_nan=False),flush=True)


def mean_force_audit(path):
    """Bound the engine's additional lateral-mean force error from tabulating W(r²).

    Cubic Hermite interpolation error is bounded by dq⁴ max|d⁴W/dq⁴|/384.
    Integrating the radial force over a plane gives 2 pi z W(z²) / cell area,
    so the energy interpolation error enters the normal force. Each
    exponential term in the fourth-derivative bound decreases with r.
    """
    rows=list(map(json.loads,(path/'guard/stdout').read_text().splitlines()))
    setup,=[v for v in rows if v['test']=='preparation']
    depth,alpha,rest=setup['pair_morse']
    dq=1./2048
    lower=.25  # actual engine table in squared distance
    radii=np.sqrt(lower+np.arange(int((32.-lower)/dq))*dq)
    def fourth(b):
        r=radii
        return np.exp(-b*(r-rest))*(b**4/r**4+6*b**3/r**5+15*b*b/r**6+15*b/r**7)/16
    energy_bound=depth*dq**4/384*(fourth(2*alpha)+2*fourth(alpha))
    zmax=np.sqrt(radii**2+dq)
    force_bound=2*np.pi*zmax*energy_bound/setup['lateral_primitive_area']
    # Below the first tabulated q and above the final one the engine uses
    # the original Morse expression. The independent bound on taper remains.
    taper=setup['mean_force_error_upper_bound']
    minima=[]
    for row in rows:
        if row['test']=='traces' and row['chi']==1.:
            with (path/(row['case']+'.jsonl')).open() as f:
                minima.append(min(json.loads(line)['min_z'] for line in f))
    measured_min=min(minima)
    sampled=force_bound[zmax>=measured_min]
    print(json.dumps(dict(test='area_mean_normal_force_bound', smallest_sampled_height=measured_min,
                          maximum_potential_interpolation_error=float(energy_bound.max()),
                          tail_taper_contribution=taper,
                          maximum_force_interpolation_plus_taper=float(force_bound.max()+taper),
                          sampled_heights_force_interpolation_plus_taper=float(sampled.max()+taper),
                          read_height_1p1_bound=float(force_bound[zmax>=1.1].max()+taper))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output_directory',type=Path)
    p.add_argument('--mean-force-audit',action='store_true')
    args=p.parse_args()
    (mean_force_audit if args.mean_force_audit else main)(args.output_directory)
