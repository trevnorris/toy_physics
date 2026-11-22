"""Reduce the two partially shear-free radius-one force and reordering checks."""
from datetime import datetime, timezone
import json
from pathlib import Path

from compare import digest, read_guarded, signed_slope

BASE = Path('_scratch/s11c')
REF = BASE / 'waist_coupled_ref01'
MAIN = BASE / 'live_waist_dynamics01'
OUTER = BASE / 'waist_outer01'
SCRIPTS = Path('software/stage1_solver/tools')


def archive(folder):
    _, fingerprint = read_guarded(folder)
    capture = folder / 'stdout'
    invocation = json.loads((folder / 'invocation.json').read_text())
    peak = max(int(json.loads(line)['memory.peak'])
               for line in (folder / 'resource-samples.jsonl').read_text().splitlines())
    return dict(command=invocation['command'], start=invocation['startedUtc'], worker_exit=0,
                finish_record_utc=datetime.fromtimestamp(
                    (folder/'outcome.json').stat().st_mtime,timezone.utc).isoformat(),
                peak_memory_bytes=peak, output=str(capture), output_bytes=capture.stat().st_size,
                sha256=fingerprint, literal_results=capture.read_text().splitlines())


def compare(name):
    fe = REF / f'p1{name}72'
    halves = REF / f'halves72{name}'
    main = MAIN / f'r1_{name}_72_support'
    outer = OUTER / {'mid': 'first_r1mid70', 'tight': 'r1tight70'}[name]
    pe, _ = read_guarded(fe)
    probes, _ = read_guarded(halves)
    reference, refhash = read_guarded(main)
    exterior, exthash = read_guarded(outer)
    assert pe[0]['sha256'] == probes[0]['sha256'] == digest(pe[0]['source'])
    assert max(pe[0]['overlap']) > .99
    directions = {}
    for work in pe[1:]:
        dim = {'radial': 'waist', 'axial': 'depth'}[work['coord']]
        m1, m2 = [signed_slope(reference, dim, h, 'total_inside_and_drives') for h in (.001, .002)]
        outer_force = next(row['slope'] for row in exterior if row['axis'] == work['coord'] and row.get('h') == .001)
        assert min(work['total']+outer_force, m1+outer_force, m2+outer_force) > 0
        directions[dim] = dict(primary_interior=[m1,m2], independent_interior=work['total'],
                               exterior=outer_force, primary_total=m1+outer_force,
                               independent_plus_primary_exterior=work['total']+outer_force)
    assert all(row['energy_slope'] > 0 for test in probes[1:] for row in test['slopes'])
    return dict(signed=directions, primary_sha256=refhash, exterior_sha256=exthash,
                independent_closed_force=archive(fe), half_patch_reopening=archive(halves))


if __name__ == '__main__':
    base = Path('research/pde_ledger_v3/_measurements/mouth_live_waist_independent.json')
    sources = {str(p): digest(p) for p in [
        base, *(SCRIPTS / s for s in [
            'mouth_live_waist_checks/compare.py', 'mouth_live_waist_checks/compare_bounded.py',
            'mouth_live_waist_checks/reorder.py', 'mouth_live_waist_checks/bounded_shape.py',
            'mouth_live_waist_fields.py', 'mouth_live_waist_order.py', 'mouth_live_waist_active_wave.py',
            'mouth_live_waist_wave.py', 'mouth_exchange_carrier.py', 'mouth_exchange_carrier_density.py',
            'mouth_sealed_shape_relax.py', 'mouth_neck_order_check.py'])]}
    print(json.dumps(dict(
        base='bca2e63c', inputs=sources,
        original_readout=str(base),
        method='P1 eigenvectors on independent straight triangles, exact-null mass removal and overlap with the measured Q1 mode. Fixed transported-order quadrature gives the virtual energy derivative, with a separate weak-carrier check. Finite local order reintroduction uses the original outgoing operator in an independent driver, resolving the fields after separate depth and radial releases. Applied steps were 0.00005 and 0.0001.',
        scope='Signed work on the supplied partially empty fields at R1. An upward cost on halved patches and a projected bounded residual are directional tests, insufficient by themselves to certify local material minima.',
        cases={name: compare(name) for name in ['mid','tight']},
    ), indent=2))
