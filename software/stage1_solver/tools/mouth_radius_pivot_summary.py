"""Reduce guarded fixed-depth radius tests and compare to authenticated archived fields."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path('research/pde_ledger_v3/_measurements/mouth_carrier_0pn.json')
DEPTH = 3.05
CASES = [
    (radius, 30, 28) for radius in (1, 2, 4, 6, 8)
] + [(6, 40, 28), (8, 40, 28)] + [
    (radius, mesh, 42) for mesh in (30, 34) for radius in (4, 6, 8)
]


def read_run(folder, expected_hash=None):
    outcome = json.loads((folder / 'outcome.json').read_text())
    stdout = (folder / 'stdout').read_bytes()
    digest = hashlib.sha256(stdout).hexdigest()
    if expected_hash and expected_hash != digest:
        raise ValueError(f'{folder} differs from archived manifest')
    if not outcome['limitsVerified']:
        raise ValueError(f'{folder} lacks guarded limits')
    exit_codes = [outcome['exitCode'], outcome['childOutcome']['exitCode']]
    invocation = json.loads((folder / 'invocation.json').read_text())
    limits = json.loads((folder / 'limit-validation.json').read_text())
    if not limits['verified'] or limits['actual']['memory.swap.max'] != '0':
        raise ValueError(f'{folder} has absent bounds or permits swap')
    return dict(exit_codes=exit_codes, stdout_sha256=digest,
                stdout_path=str(folder / 'stdout'), worker_command=invocation['command'],
                verified_limits=True, memory_bytes=invocation['memoryMax'],
                swap_bytes=invocation['swapMax'], wall_deadline=invocation['wallDeadlineSeconds']), [
        json.loads(line) for line in stdout.splitlines() if line.startswith(b'{')
    ]


def archive():
    manifest = json.loads(BASE.read_text())
    rows = []
    for label in ('one', 'two', 'four'):
        record = manifest['runs'][f'carrier_action_w{label}_01']
        folder = Path(record['stdout_path']).parent
        run, data = read_run(folder, record['stdout_sha256'])
        if run['exit_codes'] != [0, 0]:
            raise ValueError(f'Archived {folder} has a failed worker')
        wall = [row for row in data if row['test'] == 'matched_radius']
        outer = wall[-1]
        geom = data[0]
        D, C, M = outer['rate'], outer['C'], outer['M']
        radius = outer['radius']
        assert len(wall) == 2 and (geom['outer'], geom['roof'], geom['depth']) == (28., 22., 1.)
        rows.append(dict(radius=radius, mesh=geom['mesh'], radius_over_depth=radius/DEPTH,
                         D=D, C=C, M=M, D_over_R3=D/radius**3, M_over_R3=M/radius**3,
                         C_over_D=C/D, C_over_M=C/M, D_over_M=D/M, **run))
    factors = {key: rows[-1][key]/rows[0][key] for key in
               ('D_over_R3', 'M_over_R3', 'C_over_D', 'C_over_M', 'D_over_M')}
    product = factors['D_over_R3']*factors['C_over_D']/factors['M_over_R3']
    if abs(product-factors['C_over_M']) > 1.e-13:
        raise ArithmeticError('Radius-4 factor attribution failed')
    return rows, factors


def compare(root):
    rows, failures = [], []
    for radius, mesh, box in CASES:
        name = f'r{radius}_n{mesh}_box{box}'
        folder = root / name
        run, data = read_run(folder)
        if run['exit_codes'] != [0, 0]:
            partial = [row for row in data if row['test'] in ('pivot_preparation', 'regional_integrals')]
            failures.append(dict(name=name, error=(folder/'stderr').read_text().splitlines()[-1:],
                                 partial_diagnostics_only=partial, **run))
            continue
        prep = next(row for row in data if row['test'] == 'pivot_preparation')
        reg = next(row for row in data if row['test'] == 'regional_integrals')
        windows = [row for row in data if row['test'] == 'finite_radius']
        assert len(windows) == 2
        a, b = windows
        assert radius == b['radius'] and mesh == prep['mesh'] and box == prep['outer']
        rcost = reg['finite_flow_resistance']
        mzero = reg['zero_intake_added_mass']
        R3 = radius**3
        rows.append(dict(name=name, radius=radius, radius_over_depth=radius/DEPTH, mesh=mesh,
                         outer=box, roof=prep['roof'], D=b['intake'], C=b['C'], M=b['M'],
                         C_over_D=b['C_over_D'], C_over_M=b['C_over_M'], D_over_M=b['D_over_M'],
                         D_over_R3=b['normalized_rate'], M_over_R3=b['normalized_mass'],
                         alternate_window_C_over_M=a['C_over_M'], windows=windows,
                         residual_checks=[r for r in data if r['test'] == 'refined_dipole_residual'],
                         region_names=reg['region_names'], scaled_regional_phase_resistances=[v*R3 for v in rcost],
                         regional_phase_resistance_share=[v/sum(rcost) for v in rcost],
                         regional_zero_intake_mass_over_R3=[v/R3 for v in mzero],
                         zero_mass_over_finite=b['weak_M_over_finite_M'],
                         regional_phase_balance=sum(rcost)/reg['fixed_phase_drop_over_current']-1,
                         regional_mass_balance=sum(mzero)/reg['zero_intake_wall_work']-1,
                         **run))
    return rows, failures


def plot(rows, saved, destination):
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.75), sharex=True)
    styles = [(saved, 'Archived, 28 × 22', 'k', 'x', 'none'),
              ([r for r in rows if r['outer'] == 28 and r['mesh'] == 30],
               '28 × 22, mesh 30', '#084c72', 'o', '-'),
              ([r for r in rows if r['outer'] == 42], '42 × 33, mesh 34', '#b44c24', 's', '--'),
              ([r for r in rows if r['mesh'] == 40], '28 × 22, mesh 40', '#428678', '+', 'none')]
    for ax, column, title in zip(axes, ('D_over_R3', 'M_over_R3'),
                                 ('Intake divided by radius cubed', 'Inertia divided by radius cubed')):
        for data, label, color, marker, line in styles:
            if data:
                ax.plot([r['radius_over_depth'] for r in data], [r[column] for r in data],
                        label=label, color=color, marker=marker, linestyle=line,
                        markerfacecolor='none', markersize=6, linewidth=1.3)
        ax.set(xlabel='Mouth radius / axial depth', title=title)
        ax.set_ylim(bottom=0)
        ax.grid(axis='y', alpha=.24)
    axes[0].set_ylabel(r'$D/R^3$  (solver units)')
    axes[1].set_ylabel(r'$M/R^3$  (solver units)')
    axes[1].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path)
    p.add_argument('--plot', type=Path)
    args = p.parse_args()
    saved, factors = archive()
    result = dict(archive=saved, factors_R4_over_R1=factors)
    if args.root:
        rows, failures = compare(args.root)
        result.update(cases=rows, failed_cases=failures,
                      predeclaration=json.loads((args.root/'prediction.json').read_text()),
                      predeclaration_sha256=hashlib.sha256((args.root/'prediction.json').read_bytes()).hexdigest())
        source_names = ('mouth_radius_pivot', 'mouth_radius_pivot_summary', 'mouth_carrier_weak_radius',
                        'mouth_carrier_matched_bias', 'mouth_carrier_lateral_response',
                        'mouth_carrier_lateral_inertia', 'mouth_carrier_radius_family', 'mouth_curved_phase',
                        'mouth_curved_phase_stability', 'mouth_flared_receiver', 'mouth_free_contact')
        sources = [Path('software/stage1_solver/tools')/(s+'.py') for s in source_names]
        sources += [args.root/'source_initial.py', BASE, Path('scripts/s11c_guarded_run.py')]
        result['source_sha256'] = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
        if args.plot:
            plot(rows, saved, args.plot)
            result['figure'] = dict(path=str(args.plot), sha256=hashlib.sha256(args.plot.read_bytes()).hexdigest())
    print(json.dumps(result, indent=2, allow_nan=False))
