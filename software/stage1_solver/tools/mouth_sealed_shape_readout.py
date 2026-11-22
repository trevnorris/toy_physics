"""Reduce existing guarded duct relaxations without replaying any fields."""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / '_scratch/s11c/route_a_duct_relax01'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metadata(name, expected=0):
    path = DATA/name
    invocation = json.loads((path/'invocation.json').read_text())
    outcome = json.loads((path/'outcome.json').read_text())
    limit = json.loads((path/'limit-validation.json').read_text())
    if (outcome['exitCode']!=expected or outcome['childOutcome']['exitCode']!=expected
            or not outcome['limitsVerified'] or not limit['verified']
            or limit['actual']['memory.swap.max']!='0' or invocation['wallDeadlineSeconds'] is not None):
        raise ValueError(f'Guard outcome differs at {name}')
    samples = [json.loads(s) for s in (path/'resource-samples.jsonl').read_text().splitlines() if s]
    literal = (path/'stdout').read_text()
    return dict(command=invocation['command'],cwd=invocation['cwd'],pool=invocation['pool'],
                started_utc=invocation['startedUtc'],wall_seconds=outcome['wallSeconds'],
                last_resource_sample_utc=samples[-1]['utc'],guard_exit=outcome['exitCode'],
                worker_exit=outcome['childOutcome']['exitCode'],memory_bytes=invocation['memoryMax'],
                peak_sample_bytes=max(map(int,(s.get('memory.peak',0) for s in samples))),
                swap_bytes=0, wall_deadline=None,stdout_path=str((path/'stdout').relative_to(ROOT)),
                stdout_sha256=digest(path/'stdout'),stdout_size=len(literal.encode()),literal_stdout=literal,
                stderr_sha256=digest(path/'stderr'),diagnostic=(path/'stderr').read_text().splitlines()[-1:])


def rows(run):
    return [json.loads(line) for line in run['literal_stdout'].splitlines()]


def inputs():
    files = ['software/stage1_solver/tools/'+name+'.py' for name in (
        'mouth_sealed_shape_relax','mouth_exchange_shear','mouth_exchange_carrier',
        'mouth_exchange_carrier_density','mouth_exchange_geometry','mouth_exchange_interior',
        'mouth_conversion_shape_mode','mouth_sealed_shape_readout','mouth_sealed_shape_transfer')]
    files.append('research/pde_ledger_v3/_measurements/mouth_conversion_shape.json')
    return {f:digest(ROOT/f) for f in files}


def main():
    evidence, training = {}, {}
    for r in (1,2):
        prefix=f'r{r}'
        entries={'relax':prefix+'m10s','coarse_extension':prefix+('extension' if r==1 else 'continuation'),
                 'starting_control':prefix+'m16start','last_control':prefix+'m16last',
                 'fine_extension':prefix+'m16tracked'}
        evidence[prefix+'m10']=metadata(prefix+'m10',expected=1)
        for name in entries.values():
            evidence[name]=metadata(name,expected=1 if name=='r2m16tracked' else 0)
        records={k:rows(evidence[v]) for k,v in entries.items()}
        original=records['relax']; extension=records['coarse_extension']
        assert original[-1]['reason']=='thin_front'
        actual=[v for v in original if v['test']=='state']
        extra=[v for v in extension if v['test']=='state']
        fine=[v for v in records['fine_extension'] if v['test']=='state']
        for path in (original,extension,records['fine_extension']):
            states=[v for v in path if v['test']=='state']
            assert all(b['total']<a['total'] for a,b in zip(states,states[1:]))
            assert all(t['overlap']>.9 for t in path if t['test']=='accepted')
        assert extension[-1]['reason']=='waist comparable to planar half-to-90-percent width'
        if r==1:
            assert records['fine_extension'][-1]['reason']=='waist comparable to planar half-to-90-percent width'
        else:
            assert 'current fails' in evidence['r2m16tracked']['diagnostic'][0]
            evidence['r2m16transfer']=metadata('r2m16transfer')
            fixed_path=[v for v in rows(evidence['r2m16transfer']) if v['test']=='independent_mesh']
            assert fixed_path[-1]['coordinates']==extra[-1]['coordinates']
            assert all(b['total']<a['total'] for a,b in zip(fixed_path,fixed_path[1:]))
        start,last=[records[k][-1] for k in ('starting_control','last_control')]
        assert abs(start['coordinates'][0])<1e-12 and last['coordinates']==actual[-2]['coordinates']
        assert start['half_signed'][0]['total']>0 and last['half_signed'][0]['total']>0
        max_rate=max(v['feed_energy_flux_over_stored'] for v in actual)
        max_frac = max(1-v['wave_pocket_fraction'] for v in actual)
        final_coarse=extra[-1]; final_fine=fine[-1] if r==1 else fixed_path[-1]
        training[r]=dict(states=actual,extension_states=extra,diagnostic_fine_descent_states=fine,
                         fine_endpoint_method='descent' if r==1 else 'transfer along accepted coarse shapes',
                         previous_guarded_stage=actual[-1]['coordinates'],
                         initial_and_last_fine_forces=[[t['total'] for t in k['half_signed']] for k in (start,last)],
                         maximum_loss_percent_one_slow_cycle=-100*math.expm1(-max_rate*400*math.pi),
                         maximum_loss_percent_over_200=-100*math.expm1(-max_rate*200),
                         maximum_exterior_kinetic_fraction=max_frac,
                         last_fine_log_exit_slope_without_cap=last['half_signed'][1]['total']-last['half_signed'][1]['cap'],
                         last_fine_unused_quartic_slope=last['half_signed'][4]['total'],
                         energy_at_order_scale=[final_coarse['total'],final_fine['total']],
                         radius_over_interface_width_at_order_scale=[v['minimum_radius']/.04 for v in (final_coarse,final_fine)])
    evidence['r2extension']=metadata('r2extension',expected=1)
    versions={}
    for v in ('field_snapshot','extension_snapshot'):
        versions[v]={p.name:digest(p) for p in (DATA/v).glob('*.py')}
    print(json.dumps(dict(base='223de658', source_sha256=inputs(), archived_run_source_versions=versions,
                          protocol='Fixed R=1 or 2, J=1, sigma=1, cT=1, delta=.04; sharp seams carry no independent energy; final disk surface energy included. Last coordinates within 0.1 use a conventional weak-curvature reporting window. The extension is diagnostic.',
                          training=training, raw_runs=evidence),indent=2,allow_nan=False))


if __name__=='__main__':
    main()
