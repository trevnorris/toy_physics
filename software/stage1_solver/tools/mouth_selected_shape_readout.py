"""Reduce signed controlled carrier fields and compare their inner and outer windows."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def checked_run(root, radius, mode, step):
    name = f'r{radius}_{mode}_{round(step*1000):+d}'
    if mode == 'depth' and step == -.05:
        name += '_corner'
    folder = root/name
    outcome = read(folder/'outcome.json')
    child = read(folder/'child-outcome.json')
    limits = read(folder/'limit-validation.json')
    duration = read(folder/'duration-validation.json')
    invocation = read(folder/'invocation.json')
    if (outcome['exitCode'] or child['exitCode'] or not outcome['limitsVerified']
            or not limits['verified'] or int(limits['actual']['memory.swap.max']) != 0
            or not duration['verified'] or duration['wallDeadlineSeconds'] is not None):
        raise ValueError(f'Guard acceptance failed at {folder}')
    lines = (folder/'stdout').read_text().splitlines()
    records = [json.loads(line) for line in lines]
    body = [r for r in records if r['test']=='shape']
    windows = [r for r in records if r['test']=='shape_response']
    if len(body)!=1 or len(windows)!=2 or windows[-1]['window'][1]<=windows[0]['window'][1]:
        raise ValueError(f'Missing shape or response window at {folder}')
    shape = body[0]
    if (shape['radius'],shape['mode'],shape['step']) != (radius,mode,step):
        raise ValueError(f'Wrong signed preparation at {folder}')
    return dict(name=name, worker=invocation['command'], guard_exit=outcome['exitCode'],
                child_exit=child['exitCode'], limits_verified=limits['verified'],
                swap_bytes=int(limits['actual']['memory.swap.max']),
                stdout_sha256=digest(folder/'stdout'), literal_stdout=lines,
                shape=shape, windows=windows)


def reduce(root):
    out = {'runs':{}, 'signed_pairs':[]}
    for radius in (1,2):
        for mode in ('depth','bulge'):
            for h in (.05,.025):
                minus,plus = [checked_run(root,radius,mode,s) for s in (-h,h)]
                for case in (minus,plus):
                    out['runs'][case['name']] = {k:v for k,v in case.items()
                                                  if k not in ('shape','windows')}
                logscale = math.log((1+h)/(1-h))
                window_slopes = []
                for wm,wp in zip(minus['windows'],plus['windows'],strict=True):
                    slopes = {key:math.log(wp[key]/wm[key])/logscale
                              for key in ('D_over_R3','M_over_R3','D_over_M','C_over_M')}
                    if abs(slopes['D_over_M']-(slopes['D_over_R3']-slopes['M_over_R3']))>1e-11:
                        raise ArithmeticError('Intake and mass slopes do not give D/M slope')
                    window_slopes.append(slopes)
                derivatives = {key:(plus['shape'][key]-minus['shape'][key])/(2*h)
                               for key in ('finite_phase_energy','material_wall_area')}
                out['signed_pairs'].append(dict(radius=radius,mode=mode,amplitude=h,
                                                log_slopes=window_slopes,derivatives=derivatives))
    # Quantities used for cross-checking radius-one and radius-two reference
    # fields come from the earlier guarded wide-window solutions.
    archive = Path('research/pde_ledger_v3/_measurements/mouth_radius_pivot.json')
    saved = [r for r in read(archive)['cases'] if r['name'] in ('r1_n30_box28','r2_n30_box28')]
    if len(saved)!=2:
        raise ValueError('Two archived baselines are required')
    out['reference'] = {r['name']:{k:r[k] for k in ('D','M','D_over_M','C_over_M')} for r in saved}
    names = ('mouth_selected_shape_sensitivity','mouth_selected_shape_readout','mouth_radius_pivot',
             'mouth_carrier_weak_radius','mouth_carrier_radius_family','mouth_flared_receiver',
             'mouth_curved_phase','mouth_carrier_matched_bias','mouth_carrier_lateral_inertia',
             'mouth_carrier_lateral_response')
    tools = Path(__file__).parent
    sources = [tools/(name+'.py') for name in names] + [archive,Path('scripts/s11c_guarded_run.py')]
    original = root/'first_solver.py'
    for name, record in out['runs'].items():
        solver = tools/'mouth_selected_shape_sensitivity.py' if name.endswith('_corner') else original
        record['as_run_main_sha256'] = digest(solver)
    out['excluded'] = {}
    for radius in (1,2):
        case = root/f'r{radius}_depth_-50'
        fail = read(case/'outcome.json')
        if fail['exitCode'] != 1 or fail['childOutcome']['exitCode'] != 1:
            raise ValueError('Expected failed seam classification on initial negative-depth case')
        err = (case/'stderr').read_text()
        if 'assert_array_equal(self.bulk' not in err or (case/'stdout').read_text().strip():
            raise ValueError('Initial negative-depth case exited outside initial geometry classification')
        out['excluded'][case.name] = dict(reason='Rounded seam failed bath-membership assertion before a field solve',
                                           guard_exit=fail['exitCode'],traceback_sha256=digest(case/'stderr'),
                                           worker=read(case/'invocation.json')['command'],
                                           as_run_main_sha256=digest(original))
    out['sha256'] = {str(p):digest(p) for p in sources}
    out['design'] = read(root/'preparation.json')
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input_root',type=Path)
    print(json.dumps(reduce(p.parse_args().input_root),indent=2,allow_nan=False))
