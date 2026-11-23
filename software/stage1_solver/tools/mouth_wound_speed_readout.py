"""Reduce two-sided stationary circulation samples, preserving mass and angular momentum flux."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_poly_shallow_release import ROOT, emit


def audited(path):
    outcome=json.loads((path/'outcome.json').read_text())
    limits=json.loads((path/'limit-validation.json').read_text())
    if outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not outcome['limitsVerified'] or not limits['verified'] or int(limits['actual']['memory.swap.max']):
        raise ValueError(f'Unverified run {path}')
    raw=(path/'stdout').read_bytes()
    return list(map(json.loads,raw.splitlines())),hashlib.sha256(raw).hexdigest()


def quantiles(values,weights,fractions):
    sort=values.argsort()
    cum=np.cumsum(weights[sort])/weights.sum()
    return np.interp(fractions,cum,values[sort]).tolist()


def run(folder):
    cases={}
    for guard in sorted(folder.glob('a*/outcome.json')):
        path=guard.parent
        rows,digest=audited(path)
        prep=rows[0]
        radius,winding,n=prep['radius'],prep['winding'],prep['mesh']
        samples=[r for r in rows if r['kind']=='geometry']
        if prep['mode']=='signed':
            for load in [r for r in rows if r['kind']=='signed_work']:
                cases[radius,n,winding,load['half_offset']]=load['force_on_exit']
            continue
        field=path/'fields.npz'
        if not field.exists(): continue
        with np.load(field,allow_pickle=False) as save:
            weights,_,_,s,j,n0,ur,up,swirl=save['disk'].T
            if prep['material_sha256']!=save['material_sha256'].item():
                raise ValueError('Stale disk field')
        normal=(j/n0)/n0**2
        swirl_mach=abs(swirl)/n0**2
        total=np.sqrt(normal**2+(ur/n0**2)**2+(up/n0**2)**2+swirl_mach**2)
        flow=weights@j
        cap_normal=np.array([.25,.5,.75,.95])
        angular=weights@(j*s*swirl)
        emit(kind='disk',radius=radius,mesh=n,winding=winding,source_sha256=digest,
             field_sha256=hashlib.sha256(field.read_bytes()).hexdigest(),throughput=flow,
             angular_momentum_output=angular,angular_per_mass=angular/flow,
             normal_mach_above_one_area=weights[abs(normal)>1].sum()/weights.sum(),
             normal_mach_above_one_current=(weights*j)[abs(normal)>1].sum()/flow,
             azimuthal_mach_above_one_area=weights[swirl_mach>1].sum()/weights.sum(),
             total_mach_above_one_area=weights[total>1].sum()/weights.sum(),
             total_mach_above_one_current=(weights*j)[total>1].sum()/flow,
             min_total_mach=total.min(),minimum_abs_normal_mach=abs(normal).min(),
             rim_equator=dict(s=s.max(),density=n0[s.argmax()],normal_mach=normal[s.argmax()],
                              meridional_speed=np.sqrt((j/n0)**2+ur**2+up**2)[s.argmax()],
                              azimuthal_speed=swirl[s.argmax()],total_speed=(total*n0**2)[s.argmax()],
                              sound_speed=n0[s.argmax()]**2,total_mach=total[s.argmax()]),
             area_quantile_fractions=cap_normal.tolist(),normal_mach_quartiles=quantiles(normal,weights,cap_normal),
             total_mach_quartiles=quantiles(total,weights,cap_normal),
             angular_flux_error=abs(angular-winding*flow))
        if abs(flow-samples[0]['total_current'])>2e-7 or abs(angular-winding*flow)>2e-10:
            raise ArithmeticError('Conserved cap flux failed')
    data=json.loads((ROOT/'research/pde_ledger_v3/_measurements/mouth_breathing_exit_support.json').read_bytes())
    for (radius,n,m,h),force in sorted(cases.items()):
        if not m:continue
        label=f'r{radius}n72'+('_direction' if radius==1 else '')
        source=data['runs'][f'mouth_breathing_matched01/{label}']
        archive=list(map(json.loads,source['literal_stdout'].splitlines()))
        load,=[r['current_closing_load'] for r in archive if r['kind']=='signed' and r['index']==0 and r['exit_step']==h]
        zero=cases[radius,n,0,h]
        emit(kind='matched_work',radius=radius,mesh=n,winding=m,half_offset=h,
             wound_carrier_force=force,unwound_carrier_force=zero,extra_force=force-zero,
             zero_wind_closing_load=load,net_outward=force-zero-load,
             closing_source_sha256=source['stdout_sha256'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder',type=Path)
    run(p.parse_args().folder)
