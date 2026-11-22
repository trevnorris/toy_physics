"""Validate the signed open resonances used in the fixed-action stiffness test."""
import argparse
import hashlib
import json
from pathlib import Path
import math

TOOLS=Path(__file__).parent
ROOT=Path('_scratch/s11c/exchange_wave01')
GEOMETRY_RECORD=Path('research/pde_ledger_v3/_measurements/mouth_conversion_shape.json')


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_case(name):
    d=ROOT/name
    def read(f): return json.loads((d/f).read_text())
    g,child,inv,limits,duration=(read(p) for p in ('outcome.json','child-outcome.json','invocation.json',
                                                  'limit-validation.json','duration-validation.json'))
    if (g['exitCode'] or child['exitCode'] or not g['limitsVerified'] or not limits['verified']
            or limits['actual']['memory.swap.max']!='0' or duration['wallDeadlineSeconds'] is not None):
        raise ValueError(f'Unverified calculation {name}')
    output=(d/'stdout').read_text()
    return (dict(command=['python3','scripts/s11c_guarded_run.py','--pool',inv['pool'],'--memory-gib',
                          str(inv['memoryMax']//(1024**3)),'--log-directory',str(d),'--',*inv['command']],
                 guard_exit=0,worker_exit=0,limits_verified=True,swap_bytes=0,
                 wall_deadline=None,stdout_sha256=digest(d/'stdout'),literal_stdout=output),
            list(map(json.loads,output.splitlines())))


def run():
    cases=('n16_p3_f0','n24_p3_f0','n24_p4_f0','n24_p3_f0.1','interior_n16','interior_n24')
    records={};summary={}
    shape_record=json.loads(GEOMETRY_RECORD.read_text())['primary']['exchange']
    geometry=json.loads(shape_record['literal_stdout'])
    baseline,=(row for row in geometry if row['waist']==1. and row['exit_to_waist']==1.5 and row['cosines']==10)
    for name in cases:
        records[name],rows=read_case(name)
        if name.startswith('interior'):
            record,=(r for r in rows if r['test']=='second_direction')
            assert record['frequency_first']>0 and len(rows)==4
            summary[name]=record
            continue
        signed=[r for r in rows if r['test']=='signed']
        modes=[r for r in rows if r['test']=='mode']
        if len(signed)!=2 or {r['step'] for r in modes}!=set([-.01,-.005,0.,.005,.01]):
            raise ValueError('Missing signed pair or occupied reference')
        for r in modes:
            if (r['pole_residual']>3e-10 or r['overlap_with_zero_on_material_nodes']<.999
                    or abs(r['outward_power_over_interior_energy']/(2*r['amplitude_loss'])-1)>.0001):
                raise ArithmeticError('Open pole, motion or flux validation failed')
        if any(abs(s['frequency_second']-s['energy_work_curvature'])>.00015 for s in signed):
            raise ArithmeticError('Energy work does not match open pole curvature')
        K=signed[0]['frequency_second']
        base,wide=next(m for m in modes if m['step']==0),next(m for m in modes if m['step']==.01)
        ws=[]
        for r in (base,wide):
            rate=2*r['amplitude_loss']
            ws.append(dict(step=r['step'],energy_decay_rate=rate,energy_efold=1/rate,
                          loss_fraction_at_200=-math.expm1(-200*rate),
                          loss_fraction_at_400pi=-math.expm1(-400*math.pi*rate)))
        summary[name]=dict(wave_only_required_action_over_tension=-baseline['area_stiffness']/K,
                           first_per_action=signed[0]['frequency_first'],second_per_action=K,
                           widest_over_reference_loss_rate=wide['amplitude_loss']/base['amplitude_loss'],
                           retention=ws)
    src={str(p):digest(p) for p in [GEOMETRY_RECORD,*(TOOLS/n for n in (
        'mouth_exchange_geometry.py','mouth_exchange_interior.py','mouth_conversion_shape_mode.py',
        'mouth_exchange_shear.py','mouth_exchange_shear_interior.py','mouth_exchange_shear_summary.py'))]}
    prior=ROOT/'m14p3'
    previous=json.loads((prior/'outcome.json').read_text())
    if previous['exitCode']!=1 or 'Mode correspondence failed' not in (prior/'stderr').read_text():
        raise ValueError('Initial rejected-result metadata changed')
    return dict(base='017b95e2',scope='Unit radius and stored-action secular work of the degree-four material-rest-frame torsion. No physical action assigned; no radius-3 run.',
                source_sha256=src,runs=records,summary=summary,
                exclusion=dict(source_sha256=digest(ROOT/'initial_physical_sort.py'),
                               stdout_sha256=digest(prior/'stdout'),stderr_sha256=digest(prior/'stderr'),
                               stage='Displaced eigenvector overlap used a physical-coordinate ordering which changes under deformation; excluded before readout.'))


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
