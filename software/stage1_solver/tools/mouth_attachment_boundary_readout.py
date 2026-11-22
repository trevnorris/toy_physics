"""Reduce signed pole early-capture measures after continuous clearance integration."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, trapezoid

from mouth_attachment_incidence import HalfOrderBand
from mouth_attachment_capture_boundary import flux_integral
from mouth_attachment_impact_trial import sha
from neutral_patch_sites import emit

def read(path):
    outcome=json.loads((path/'outcome.json').read_text())
    assert outcome['limitsVerified'] and outcome['exitCode']==outcome['childOutcome']['exitCode']==0
    lines=[json.loads(r) for r in (path/'stdout').read_text().splitlines()]
    setup,=[r for r in lines if r['test']=='boundary_setup']
    audit,=[r for r in lines if r['test']=='boundary_audit']
    assert sha(Path(audit['path']))==audit['sha256']
    return setup,[r for r in lines if r['test']=='clearance_integral'],audit

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('blocks',nargs='+',type=Path)
parser.add_argument('--gated-time',action='store_true',
                    help='Integrate truncated gated joint time on the paired clearance meshes')
parser.add_argument('--controls',action='store_true',
                    help='Compare the opposite poles across a single matched set of launch rays')
args=parser.parse_args()
blocks=[read(p) for p in args.blocks]
setups,audits=[b[0] for b in blocks],[b[2] for b in blocks]
if args.controls:
    common=min(s['fine_cells'] for s in setups)
    seed=setups[0]['seed']
    band=HalfOrderBand.load(Path(setups[0]['section']))
    u=.002
    for path,(setup,rows,audit) in zip(args.blocks,blocks,strict=True):
        if (setup['seed']!=seed or setup['probe_count']!=setups[0]['probe_count']
                or setup['section_sha256']!=setups[0]['section_sha256']
                or setup['engine_sha256']!=setups[0]['engine_sha256']):
            raise ValueError('Controls must share launch rays, section and engine')
        r={(v['probe'],v['drive']):v for v in rows if v['cells']==common}
        idx=range(setup['probe_count'])
        z,plus,minus=(np.array([r[k,w]['mass_over_undriven_incident_rate'] for k in idx])
                      for w in (0.,u,-u))
        c0,cp,cm=z.mean(),plus.mean(),minus.mean()
        signed=(cp-cm)/(2*u*c0)
        fixed=np.mean([flux_integral(band,u,r[k,0.]['intervals'])-
                       flux_integral(band,-u,r[k,0.]['intervals']) for k in idx])/(2*u*c0)
        end_error=np.mean([r[k,u]['localized_boundary_mass_uncertainty']+
                           r[k,-u]['localized_boundary_mass_uncertainty'] for k in idx])/(2*u*c0)
        emit(test='paired_control',path=str(path),source_sha256=sha(path/'stdout'),seed=seed,cells=common,
             dt=setup['dt'],command=setup['command'],min_initial_energy=audit['min_incoming_tag_energy'],
             max_total_energy_error=audit['max_total_energy_error'],zero_capture=c0,
             near_pole_capture_supply=cp,far_pole_capture_supply=cm,odd_capture_supply=signed,
             frozen_interval_slope=fixed,change_from_interval_motion=signed-fixed,
             bracket_mass_error_over_2U_c0=end_error,
             probe_zero_captures=z.tolist(),probe_signed_raw_derivatives=((plus-minus)/(2*u)).tolist(),
             probe_interval_counts=[[len(r[k,w]['intervals']) for w in (0.,u,-u)] for k in idx])
    raise SystemExit(0)
cells=setups[0]['fine_cells']
drives=setups[0]['drives']
for setup in setups:
    for key in ('offset','depth_gate','native_time','dt','cosine','fine_cells','probe_count','engine_sha256','section_sha256','drives'):
        assert setup[key]==setups[0][key]
band=HalfOrderBand.load(Path(setups[0]['section']))
incoming={u:flux_integral(band,u,[(0.,1.)]) for u in drives}
records={}
for setup,rr,_ in blocks:
    for r in rr:
        records[setup['seed'],r['probe'],r['drive'],r['cells']]=r
probes=[(s['seed'],k) for s in setups for k in range(s['probe_count'])]
def vector(key,u,level,field):
    return np.array([records[*p,u,level][field] for p in key])
def report(key,level,scope):
    def avg(u):
        return np.mean(vector(key,u,level,'mass_over_undriven_incident_rate'))
    c0=avg(0.)
    for u in (.002,.001):
        cp,cm=avg(u),avg(-u)
        b=(cp-cm)/(2*u*c0) if c0>0 else None
        if c0>0:
            fixed=[records[*p,0.,level]['intervals'] for p in key]
            supply=np.mean([flux_integral(band,u,segments)-flux_integral(band,-u,segments)
                            for segments in fixed])/(2*u*c0)
            bp=np.mean(vector(key,u,level,'localized_boundary_mass_uncertainty'))
            bm=np.mean(vector(key,-u,level,'localized_boundary_mass_uncertainty'))
            b0=np.mean(vector(key,0.,level,'localized_boundary_mass_uncertainty'))
            # Enclose either sign using the measured root-bracket widths.
            # Entirely unseen intervals between identically labelled endpoints
            # are outside this numerical bound.
            bounds=[(cp-cm+sgn*(bp+bm))/(2*u*(c0+sgn0*b0))
                    for sgn in (-1,1) for sgn0 in (-1,1)]
            bracket=[min(bounds),max(bounds)]
        else:
            bracket=None
            supply=None
        emit(test='capture_response',scope=scope,cells=level,signed_magnitude=u,
             fixed_rays=len(key),zero_capture=c0,near_pole_capture_supply=cp,far_pole_capture_supply=cm,
             near_prob_per_incident=cp/incoming[u],far_prob_per_incident=cm/incoming[-u],
             raw_odd_capture_derivative=(cp-cm)/(2*u),fractional_odd_supply=b,
             fractional_slope_frozen_capture_intervals=supply,
             additional_slope_from_changed_intervals=None if b is None else b-supply,
             localized_bracket_only_beta_range=bracket)

emit(test='inputs',blocks=[s['seed'] for s in setups],audits=audits,
     inputs={str(p):sha(p/'stdout') for p in args.blocks})
if args.gated_time:
    times={}
    for setup,_,audit in blocks:
        with Path(audit['path']).open() as stream:
            for line in stream:
                r=json.loads(line)
                times[setup['seed'],r['probe'],r['drive'],r['clearance']]=r['gated_joint_residence']
    for n in (cells//4,cells//2,cells):
        s=np.linspace(0.,1.,n+1)**2
        for method,integral in (('trapezoid',trapezoid),('simpson',simpson)):
            means={}
            for seed,k in probes:
                for u in drives:
                    tau=np.array([times[seed,k,u,v] for v in s])
                    means[seed,k,u]=integral(tau*band.launch(u,1.,s)['weight'],x=s)
            groups=[(str(a['seed']),[p for p in probes if p[0]==a['seed']]) for a in setups]
            groups.append(('whole',probes))
            for scope,key in groups:
                c0=np.mean([means[*p,0.] for p in key])
                for u in (.002,.001):
                    cp,cm=[np.mean([means[*p,v] for p in key]) for v in (u,-u)]
                    emit(test='gated_time_response',scope=scope,rule=method,cells=n,
                         signed_magnitude=u,probes=len(key),undriven_mean_native_time=c0,
                         near_pole_time_supply=cp,far_pole_time_supply=cm,
                         raw_odd_time_derivative=(cp-cm)/(2*u),
                         fractional_odd_time_supply=None if c0<=0 else (cp-cm)/(2*u*c0))
    raise SystemExit(0)
for n in (cells//2,cells):
    for setup in setups:
        report([p for p in probes if p[0]==setup['seed']],n,str(setup['seed']))
    report(probes,n,'whole')
for lo,hi in ((0.,.25),(.25,.5),(.5,.75),(.75,1.)):
    subset=[p for p in probes if lo<=records[*p,0.,cells]['draw'][3]<hi]
    if subset: report(subset,cells,f'thermal_quantile_{lo}_{hi}')
active=[p for p in probes if records[*p,0.,cells]['intervals']]
zero=vector(probes,0.,cells,'mass_over_undriven_incident_rate')
unequal=[]
multis=0
for p in probes:
    for u in drives:
        lo=records[*p,u,cells//2]
        hi=records[*p,u,cells]
        if len(lo['intervals'])!=len(hi['intervals']):
            unequal.append([*p,u,len(lo['intervals']),len(hi['intervals'])])
        multis+=int(len(hi['intervals'])>1)
emit(test='coverage',fixed_rays=len(probes),rays_with_zero_drive_interval=len(active),
     zero_mass_probe_range=[min(zero),max(zero)],resolved_disconnected_cases=multis,
     cases_with_changed_interval_count=unequal,
     max_absolute_change_between_nested_flux_estimates=max(
         abs(records[*p,u,cells]['mass_over_undriven_incident_rate']-
             records[*p,u,cells//2]['mass_over_undriven_incident_rate']) for p in probes for u in drives))
