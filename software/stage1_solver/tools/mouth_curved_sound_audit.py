"""Independently reduce the native radius-two resonance rows without solving them again."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def pairs(x):
    a = np.asarray(x, dtype=float)
    return a[..., 0] + 1j*a[..., 1]


def read(case):
    root = case.parent
    raw = (root/'stdout').read_bytes()
    guard = json.loads(case.read_bytes())
    child = json.loads((root/'child-outcome.json').read_bytes())
    if guard['exitCode'] or child['exitCode'] or not guard['limitsVerified']:
        raise ArithmeticError(f'Unvalidated input {root}')
    data = [json.loads(s) for s in raw.splitlines()]
    if [r['kind'] for r in data] != ['input', 'mode', 'traces']:
        raise ArithmeticError(f'Incomplete input {root}')
    pre, mode, traces = data
    if pre['material_sha256'] != '954a5d6d2f4a679e4bb84a437c671458ec9a79c5abf7290c1015d49380eb86be':
        raise ArithmeticError('Different radius-two material')
    args = pre['args']
    setting = tuple(args[k] for k in ('mesh', 'box', 'pml', 'stretch', 'front'))
    omega, decay = mode['pole'][0], -mode['pole'][1]
    Q = omega/(2*decay)
    fraction = -np.expm1(-2*np.pi/Q)
    power = mode['real_frequency_power']
    feed, bulk = power['feed_cut'], power['bulk_cut']
    work = mode['real_frequency_external_work']
    integrated = np.array(traces['weights'])
    area = 4*np.pi*pre['physical_exit']**3/3
    vector = pairs(mode['mass_unit_vector'])
    mass_norm = (vector.conj()@np.array(pre['mass'])@vector).real
    cut_total = feed + bulk
    E = mode['peak_unit_stored_energy']
    weight_errors = {}
    for key, target in [('relative_current','disk_relative_current'), ('absolute_current','disk_absolute_current'), ('front_force','disk_front_reaction')]:
        midpoint = integrated@pairs(traces[key])
        stored = pairs(mode.get('midpoint_'+target, mode[target]))
        weight_errors[key] = abs(midpoint - stored)
    speed = pairs(mode['disk_sweep'])/area
    density = (pairs(traces['absolute_current'])-pairs(traces['relative_current']))/speed
    # Cap velocity is uniform for the saved three shape coordinates. The two
    # frames must differ at every radial sample by precisely the local rho*V.
    errs = (abs(Q-mode['outgoing_Q'])/Q,abs(fraction-mode['pole_cycle_loss']),
            abs(integrated.sum()/area-1),abs(mass_norm-1),np.max(abs(density.imag)),
            max(weight_errors.values()))
    if max(errs)>1.e-8 or density.real.min()<=0 or not 0<decay<omega:
        raise ArithmeticError(f'Bad normalization, reference frame or quadrature in {root}: {errs}')
    if not np.isclose(cut_total,work,rtol=1.e-6) or min(feed,bulk)<=0:
        raise ArithmeticError(f'Cut power imbalance in {root}')
    per_cycle = cut_total*2*np.pi/(omega*E)
    return (setting,args['offset']), dict(
        kind='center_or_offset', source=str(root), stdout_sha256=hashlib.sha256(raw).hexdigest(),
        setting=setting, offset=args['offset'], omega=omega, decay=decay, Q=Q,
        radius=pre['physical_exit'], pole_loss=fraction, energy_efold_periods=Q/(2*np.pi),
        reservoir_cut_energy_per_fixed_amplitude_period=per_cycle,
        feed_fraction_of_cut_power=feed/cut_total,
        bulk_fraction_of_cut_power=bulk/cut_total,
        loss_from_real_frequency_energy_estimate=-np.expm1(-per_cycle),
        work_minus_cut_fraction=(work-cut_total)/work,
        effective_density_range_at_cap=[density.real.min(),density.real.max()],
        min_mass_overlap=min(s[3] for s in mode['tracking']),
        unit_effective_energy=E, coordinate=mode['mass_unit_vector'])


def write_stroke_spec(path, measured, setting):
    from scipy.optimize import minimize_scalar
    from mouth_curved_sound_stroke import finite_shapes
    from mouth_curved_steady import reference
    from mouth_poly_shallow_release import ROOT
    h = .006
    center = measured[setting, 0.]
    p, m = measured[setting, h], measured[setting, -h]
    slope = (m['omega']-p['omega'])/(p['radius']-m['radius'])
    shape, _, _, scale, outer, rows, digest = reference()
    signed, = [r for r in rows if r['kind']=='signed' and r['exit_step']==h and r['index']==0]
    if slope<=0:
        raise ValueError('There is no positive conditional balancing action')
    J = signed['current_closing_load']/slope
    v = pairs(center['coordinate'])*np.sqrt(J*center['omega']/center['unit_effective_energy'])
    def outer_strain(angle):
        step = (v*np.exp(1j*angle)).real
        return finite_shapes(shape,scale,outer,step)['outermost_registered']['peak_absolute_infinitesimal_strain']
    angles = np.arange(32)*np.pi/32
    peak = angles[np.argmax([outer_strain(x) for x in angles])]
    fit = minimize_scalar(lambda x:-outer_strain(x),bounds=(peak-np.pi/32,peak+np.pi/32),
                          method='bounded',options=dict(xatol=1.e-8))
    if not fit.success:
        raise ArithmeticError('Phase search failed')
    phase = fit.x%np.pi
    phase -= np.pi if phase>np.pi/2 else 0.
    spec = dict(source=str((Path(center['source'])/'stdout').resolve().relative_to(ROOT)),
                source_sha256=center['stdout_sha256'],material_sha256=digest,
                label='conditional mean-front signed pole-frequency slope',relative_exit_step=h,
                frequency_slope_force=slope, formal_action=J,
                signed_readout_sha256=[m['stdout_sha256'],p['stdout_sha256']],phase_at_stroke=phase,
                peak_outer_infinitesimal_strain=outer_strain(phase),
                coordinate_peak=(v*np.exp(1j*phase)).real.tolist())
    path.write_text(json.dumps(spec,indent=2)+'\n')


def run(cases, stroke_spec):
    measured = {}
    def emit(**data):
        print(json.dumps(data, allow_nan=False), flush=True)
    for case in cases:
        key,row = read(case)
        if key in measured:
            raise ValueError(f'Duplicate numerical settings {key}')
        measured[key] = row
        emit(**row)
    for setting,offset in sorted(measured):
        if offset!=0.:
            continue
        center = measured[setting,offset]
        for h in (.012,.006):
            if (setting,h) in measured and (setting,-h) in measured:
                positive = measured[setting,h]; negative = measured[setting,-h]
                force = (negative['omega']-positive['omega'])/(positive['radius']-negative['radius'])
                def dynamic_stiffness(row, name):
                    mode = json.loads((Path(row['source'])/'stdout').read_text().splitlines()[1])
                    A = pairs(mode[name])
                    return (A + A.conj().T)/2
                difference = dynamic_stiffness(positive,'anchor_stiffness')-dynamic_stiffness(negative,'anchor_stiffness')
                v = pairs(center['coordinate'])
                restoring = (v.conj()@difference@v).real/(positive['radius']-negative['radius'])
                normalization = (v.conj()@dynamic_stiffness(center,'anchor_slope')@v).real
                emit(kind='signed_frequency_derivative',setting=setting,relative_step=h,
                     effective_outward_force_per_action=force,ratio_to_old_local_force=force/1.07730,
                     reactive_virtual_work_at_real_anchor=restoring/normalization,
                     Q_required_for_one_percent_loss=2*np.pi/-np.log1p(-.01),
                     measured_Q=center['Q'],may_proceed_on_loss_alone=bool(center['pole_loss']<.01))
    if stroke_spec is not None:
        settings = [s for s,offset in measured if offset==0 and s[-1]=='mean' and
                    all((s,h) in measured for h in (-.006,.006))]
        write_stroke_spec(stroke_spec,measured,max(settings))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('outcomes',nargs='+',type=Path)
    parser.add_argument('--stroke-spec',type=Path)
    args=parser.parse_args()
    run(args.outcomes,args.stroke_spec)
