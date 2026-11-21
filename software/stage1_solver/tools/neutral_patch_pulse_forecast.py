"""Use single-mouth measurements to forecast previously unseen two-drain loads.

The four inputs span cubic polynomials of absolute age during one specified
conversion transit. Their measured responses can be added to first order even
if material properties evolve during conversion. A separate release measures
the sensitivity to the specified, initially established empty-channel profile.
The forecast reads no backing labels or companion load measurements.
"""
import argparse
from math import comb
import json
from pathlib import Path

import numpy as np

from mouth_contact_tracking import WallTrace
from mouth_order_wave_pulse import digest, emit
from neutral_patch_companion_rim import driving
from neutral_patch_conversion import clock

FIELD = Path('_scratch/s11c/companion_patch_fields_v1.npz')
FIELD_HASH = '859fc720afeff889e864a8873bc514110ff998555ff631045a2e93dfb9cbe192'
FOLDER = Path('_scratch/s11c/pulse_probe01')
TRUTH = Path('research/pde_ledger_v3/_measurements/neutral_patch_reciprocal_flow.json')
PRIMARY = ('lower_drag_cv','rim')
READINGS = (*PRIMARY,'fluid_lower','material_lower','lower_surface_shear','lower_fluid_storage')


def bernstein(p,degree):
    return np.stack([comb(degree,k)*p**k*(1-p)**(degree-k) for k in range(degree+1)],axis=1)


def waveform_fit(t,velocity,degree):
    # Use the same left endpoints the numerical experiment advances.
    durations = np.diff(t)
    p = (t[:-1]-t[0])/(t[-1]-t[0])
    basis = bernstein(p,degree)
    root = np.sqrt(durations/durations.sum())
    c,_,rank,_ = np.linalg.lstsq(root[:,None]*basis,root*velocity[:-1],rcond=None)
    if rank != degree+1:
        raise ArithmeticError('Dependent input pulses')
    error = basis@c-velocity[:-1]
    # Degree elevation represents every lower-order fit in the measured basis.
    for d in range(degree+1,4):
        fractions = np.arange(1,d)/d
        c = np.r_[c[0],fractions*c[:-1]+(1-fractions)*c[1:],c[-1]]
    np.testing.assert_allclose(bernstein(p,3)@c,velocity[:-1]+error,rtol=3e-14,atol=1e-17)
    return dict(coefficients=c.tolist(),degree=degree,initial=float(velocity[0]),
                relative_rms=float(np.linalg.norm(root*error)/np.linalg.norm(root*velocity[:-1])),
                peak_normalized_max=float(max(abs(error))/max(abs(velocity))))


def lookups(packet,keys,field):
    raw = [str(s) for s in packet[field].tolist()]
    if len(raw)!=len(set(raw)):
        raise ValueError(f'Duplicate {field} label')
    return [raw.index(key) for key in keys]


def run_forecast(args):
    if digest(FIELD)!=FIELD_HASH:
        raise ValueError('External velocity history differs from the pinned file')
    tr = WallTrace(48,.5)
    t,w,_,_,_ = clock(tr,2)
    speed = driving(tr,FIELD,w,.01)[:,[1,3]]
    fits = {}
    for k,r in enumerate((8,16)):
        fits[r] = [waveform_fit(t,speed[:,k],degree) for degree in (1,2,3)]
        emit(test='drive_approximation',separation=r,models=fits[r])

    rows = []
    inputs = {str(FIELD):FIELD_HASH}
    for code in 'ABCD':
        for law in ('central','layered'):
            filename = args.directory/f'{code}_{law}_fine.npz'
            inputs[str(filename)] = digest(filename)
            # Only guard status and the final artifact hash enter here. Opening
            # the acquisition command or setup would reveal the blind condition.
            guard = args.directory/f'run_{code}_{law}_fine'
            outcome = json.loads((guard/'outcome.json').read_text())
            limits = json.loads((guard/'limit-validation.json').read_text())
            sealed = json.loads((guard/'stdout').read_text().splitlines()[-1])
            if (outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not limits['verified']
                    or sealed.get('test')!='probe_saved' or sealed['sha256']!=inputs[str(filename)]):
                raise ValueError('Pulse measurement failed provenance checks')
            with np.load(filename,allow_pickle=False) as packet:
                amplitude = float(packet['amplitude'])
                if amplitude!=.001:
                    raise ValueError('Unexpected calibration amplitude')
                np.testing.assert_allclose([packet['t'][0],packet['t'][-1]],t[[0,-1]],rtol=1e-12)
                col = lookups(packet,READINGS,'obs')
                coef = lookups(packet,('b0','b1','b2','b3','release'),'case_names')
                cases = packet['response']
                zero, = lookups(packet,('zero',),'case_names')
                np.testing.assert_allclose(cases[zero,col],0.,atol=1e-14,rtol=0.)
                gains = cases[np.ix_(coef,col)]/amplitude
                if not np.isfinite(gains).all():
                    raise ValueError('Nonfinite single-mouth response')
            for r,models in fits.items():
                for model in models:
                    profile_result = np.array(model['coefficients'])@gains[:4]
                    for state in ('cold','warm'):
                        load = profile_result.copy()
                        if state=='warm':
                            load += model['initial']*gains[4]
                        row = dict(code=code,law=law,separation=r,preparation=state,degree=model['degree'],
                                   predicted=dict(zip(READINGS,load.tolist(),strict=True)))
                        rows.append(row)
                        if model['degree']==3:
                            emit(test='blind_load',**row)
    payload = dict(method='Least squares on incident velocity only, then single-mouth response superposition',
                   inputs_sha256=inputs,approximation=fits,forecast=rows)
    with args.output.open('x') as stream:
        json.dump(payload,stream,allow_nan=False,indent=2)
        stream.write('\n')
    emit(test='forecast_freeze',path=str(args.output),sha256=digest(args.output),rows=len(rows))


def validate(args):
    if digest(args.forecast)!=args.sha256:
        raise ValueError('Forecast differs from the frozen digest')
    data = json.loads(args.forecast.read_text())
    mapping = json.loads(args.mapping.read_text())
    if set(mapping)!=set('ABCD') or set(mapping.values())!={'free','rest','z05','z4'}:
        raise ValueError('Expected blind codes for the four tested backings')
    truth = {}
    for run in json.loads(TRUTH.read_text())['runs'].values():
        if run['guard_exit'] or run['child_exit'] or not run['limits_verified']:
            raise ValueError('Unverified target archive')
        rows = list(map(json.loads,run['literal_stdout']))
        h = rows[0]
        if h.get('test')!='setup' or h['grid']!=8 or h['refine']!=2:
            continue
        label = h['backing']
        if label=='robin':
            label = {.5:'z05',4.:'z4'}[h['back_impedance']]
        state = {'cold':'cold','steady_empty_channel':'warm'}[h['initialization']]
        if (h['rho'],h['nu'],h['gamma'],h['grip_k'],h['back_gap'])!=(1,1,1,1,1.5):
            raise ValueError('Different material or grip')
        for r in (8,16):
            key=(label,h['law'],state,r)
            if key in truth:
                raise ValueError('Duplicate comparison')
            truth[key] = rows[-1]['cases'][f'toward{r}']
    if len(truth)!=24:
        raise ValueError('Incomplete comparison set')
    groups={}
    for row in data['forecast']:
        state,law,r,degree = (row[k] for k in ('preparation','law','separation','degree'))
        key=(mapping[row['code']],law,state,r)
        if key not in truth:
            if state!='cold' or key[0] not in ('z05','z4'):
                raise ValueError('Unexpected missing comparison')
            continue
        target=np.array([truth[key][k] for k in PRIMARY])
        pred=np.array([row['predicted'][k] for k in PRIMARY])
        if np.any(target==0):
            raise ValueError('Relative error undefined for zero target')
        errors=abs(pred-target)
        fractions=errors/abs(target)
        signs=(np.sign(target)==np.sign(pred)).astype(int)
        groups.setdefault(degree,[]).append((errors,fractions,signs))
        if degree==3:
            emit(test='holdout',code=row['code'],backing=key[0],law=law,preparation=state,separation=r,
                 ordering=PRIMARY,target=target.tolist(),predicted=pred.tolist(),relative_error=fractions.tolist())
    for degree,group in sorted(groups.items()):
        errors,fractions,signs = map(np.stack,zip(*group,strict=True))
        if len(group)!=len(truth):
            raise ArithmeticError('Incomplete predictions for degree')
        emit(test='holdout_summary',degree=degree,count=len(group),ordering=PRIMARY,
             max_absolute_error=errors.max(axis=0).tolist(),max_relative_error=fractions.max(axis=0).tolist(),
             median_relative_error=np.median(fractions,axis=0).tolist(),matched_signs=signs.sum(axis=0).tolist())
    emit(test='compared_files',forecast_sha256=args.sha256,target_sha256=digest(TRUTH),
         blind_key_sha256=digest(args.mapping),missing_comparators='cold finite-resistance cases')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_subparsers(dest='mode',required=True)
    predict=modes.add_parser('forecast')
    predict.add_argument('--directory',type=Path,default=FOLDER)
    predict.add_argument('--output',type=Path,required=True)
    compare=modes.add_parser('compare')
    compare.add_argument('--forecast',type=Path,required=True)
    compare.add_argument('--sha256',required=True)
    compare.add_argument('--mapping',type=Path,required=True)
    args=parser.parse_args()
    if args.mode=='forecast':
        run_forecast(args)
    else:
        validate(args)
