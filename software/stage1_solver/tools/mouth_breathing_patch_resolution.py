"""Validate the two-sided planar response with finer resolution; adopt existing zero-current readings."""
import argparse
import hashlib
import json

from mouth_breathing_mobile_patch import ROOT, preparation, response, emit


def run(radius, matched):
    case, occupied_frequency, source=preparation(radius)
    j=case['current']/case['cap_three_area']
    if matched:
        path=ROOT/f'_scratch/s11c/mouth_breathing_matched01/r{radius}n72'
        if radius==1:
            path=path.with_name(path.name+'_direction')
        meta=json.loads((path/'outcome.json').read_bytes())
        if meta['childOutcome']['exitCode'] or meta['exitCode'] or not meta['limitsVerified']:
            raise ValueError('Unverified empty-background frequency')
        raw=(path/'stdout').read_bytes()
        rows=list(map(json.loads,raw.splitlines()))
        chosen,=[r for r in rows if r.get('index')==0 and r.get('exit_step')==.006 and 'omega' in r]
        omega=chosen['omega']
        source=hashlib.sha256(raw).hexdigest()
    else:
        omega=occupied_frequency
        path=ROOT/f'_scratch/s11c/mouth_breathing_mobile_patch01/r{radius}'
        meta=json.loads((path/'outcome.json').read_bytes())
        raw=(path/'stdout').read_bytes()
        previous=list(map(json.loads,raw.splitlines()))
        if (meta['exitCode']!=1 or not meta['limitsVerified'] or len(previous)!=3
                or any(r['normal_current']!=0 or r['frequency']!=omega for r in previous[1:])
                or 'Window or spatial resolution affected response' not in (path/'stderr').read_text()):
            raise ValueError('Invalid previously completed zero-current resolutions')
        emit(kind='adopt_zero',source_sha=hashlib.sha256(raw).hexdigest(),
             prior_refinements=[r['mesh_refinement'] for r in previous[1:]])
    emit(kind='patch_input',radius=radius, matched_empty=matched, frequency=omega, source_sha=source,
         shape_hash=case['endpoint_sha256'], measured_mean_flux=j)
    readings={}
    for current in (0.,j,-j):
        for end,n in [(16.,4),(21.,6)]:
            row=response(current,omega,end,n)
            emit(kind='two_sided_patch',**row)
            if min(row['toward_feed'], row['toward_receiver'])<=0 or row['acoustic_work_relative_error']>2e-5:
                raise ArithmeticError('Unresolved actuator and two-reservoir power')
            readings[current,n]=row
        a,b=readings[current,4],readings[current,6]
        keys=('normal_impedance_re','normal_impedance_im','toward_receiver','toward_feed')
        errors={k:abs(a[k]/b[k]-1) for k in keys}
        emit(kind='resolution',normal_current=current,relative_errors=errors)
        if max(errors.values())>3e-4:
            raise ArithmeticError('Normal response is still resolution sensitive')
    pos,neg=readings[j,6],readings[-j,6]
    errors=[abs(pos[a]/neg[b]-1) for a,b in [('toward_feed','toward_receiver'),('toward_receiver','toward_feed')]]
    emit(kind='reversal',signed_flux=j,port_errors=errors)
    if max(errors)>2e-5:
        raise ArithmeticError('Doppler port reversal failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',choices=(1,2),type=int,required=True)
    p.add_argument('--matched-empty',action='store_true')
    args=p.parse_args()
    run(args.radius,args.matched_empty)
