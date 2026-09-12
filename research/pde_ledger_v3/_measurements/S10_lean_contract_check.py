#!/usr/bin/env python3
"""Reproduce the bounded S10 action-fidelity controls; no CAS translation.

Canonical sources are never edited. Each source mutation keeps the original
proofs and must fail within an identified mathematical declaration. The original
module must first compile. Existing spectrum/stratum and Q6/Q7 controls are
referenced by the coverage contract rather than duplicated here.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import time

BASE = Path(__file__).resolve().parents[1]
LEAN = BASE / 'lean'
SCRATCH = LEAN / 's10/_scratch/contract'
CASES = (
    ('curl_normalization', 'S10Pilot/Action.lean',
     '(1 / 2 : ℝ) * ∑ i, ∑ j, antisym J i j ^ 2',
     '(1 : ℝ) * ∑ i, ∑ j, antisym J i j ^ 2',
     {'lagrangian_increment', 'modalAction_eq', 'antisym_mode_pair'}),
    ('fullgrad_replaced_by_divergence', 'S10Controls/Action.lean',
     'def fullGradientStiffness (J : Jet D) : ℝ := ∑ i : Fin D, ∑ j : Fin D, J i.succ j ^ 2',
     'def fullGradientStiffness (J : Jet D) : ℝ := (∑ i : Fin D, J i.succ i) ^ 2',
     {'lagrangian_increment', 'fullGradient_mode_pair', 'stiffnessPair_self'}),
    ('divergence_sum_of_squares', 'S10Controls/Action.lean',
     'def divergenceOnlyStiffness (J : Jet D) : ℝ := divergence J ^ 2',
     'def divergenceOnlyStiffness (J : Jet D) : ℝ := ∑ i : Fin D, J i.succ i ^ 2',
     {'lagrangian_increment', 'stiffnessPair_self'}),
    ('anisotropy_all_inertias', 'S10Anisotropic/Action.lean',
     'def kinetic (e : Fin D) (sigma : ℝ) (v : Vec D) : ℝ :=\n  ∑ i, (if i = e then sigma else 1) * v i ^ 2',
     'def kinetic (e : Fin D) (sigma : ℝ) (v : Vec D) : ℝ :=\n  ∑ i, sigma * v i ^ 2',
     {'kinetic_eq'}),
    ('ignore_scalar_control', 'S10Controls/Scalar.lean',
     'rho / 2 * normSq (J 0) - c * mu / 2 * stiffness J',
     'rho / 2 * normSq (J 0) - mu / 2 * stiffness J',
     {'lagrangian_eq', 'lagrangian_variation', 'modalAction_eq'}),
    ('phase_average_missing_half', 'S10Pilot/PhaseAverage.lean',
     'phaseAverage rho mu omega k a = (1 / 2 : ℝ) * modalAction rho mu omega k a := by',
     'phaseAverage rho mu omega k a = modalAction rho mu omega k a := by',
     {'phaseAverage_eq'}),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(path, outcome, required=()):
    source = path.read_text()
    rel = str(path.relative_to(LEAN))
    command = ['lake', 'env', 'lean', '-DwarningAsError=true', rel]
    started = time.monotonic()
    result = subprocess.run(command, cwd=LEAN, text=True, capture_output=True,
                            env={**os.environ, 'LAKE_CACHE_DIR': '.lake/cache'}, timeout=180)
    output = result.stdout + result.stderr
    declarations = [(source[:m.start()].count('\n') + 1, m[1]) for m in
                    re.finditer(r'^\s*(?:@\[[^\n]*\]\s*)?theorem\s+(\w+)', source, re.M)]
    errors = [int(n) for n in re.findall(re.escape(rel) + r':(\d+):\d+: error:', output)]
    failed = sorted({next((name for line, name in reversed(declarations) if line <= n), '')
                     for n in errors} - {''})
    name = path.stem if '_scratch' in path.parts else 'control_' + path.parent.name + '_' + path.stem
    record = {'name': name, 'expected_outcome': outcome, 'command': command,
              'source_sha256': digest(path), 'exit_status': result.returncode,
              'wall_seconds': round(time.monotonic() - started, 3),
              'failed_declarations': failed, 'required_failure_in': sorted(required),
              'output': output, 'output_sha256': hashlib.sha256(output.encode()).hexdigest()}
    if outcome == 'PASS':
        valid = result.returncode == 0 and not re.search(r'\b(?:error|warning):', output)
    else:
        valid = (result.returncode == 1 and bool(set(failed) & set(required)) and
                 bool(re.search(r'unsolved goals|Type mismatch|type mismatch|failed to close|proved that the proposition', output)) and
                 not re.search(r'unknown (?:module|namespace|identifier)|unexpected token|maximum (?:recursion|heartbeats)', output))
    record['outcome'] = outcome if valid else 'UNEXPECTED'
    print(record['name'] + ': ' + record['outcome'], flush=True)
    assert valid, record
    return record


def main():
    SCRATCH.mkdir(parents=True, exist_ok=True)
    sources = sorted(p for step in ('s9', 's10') for p in (LEAN / step).rglob('*.lean')
                     if '_scratch' not in p.parts)
    before = {str(p.relative_to(BASE)): digest(p) for p in sources}
    records, controls = [], set()
    for name, rel, old, new, required in CASES:
        canonical = LEAN / 's10' / rel
        source = canonical.read_text()
        if rel not in controls:
            records.append(run(canonical, 'PASS'))
            controls.add(rel)
        assert source.count(old) == 1, (name, source.count(old))
        path = SCRATCH / (name + '.lean')
        path.write_text(source.replace(old, new))
        record = run(path, 'REJECTED', required)
        record.update(canonical_source=str(canonical.relative_to(BASE)),
                      canonical_sha256=digest(canonical), replacement={'old': old, 'new': new})
        records.append(record)

    scalar = '''import S10Controls.Scalar
open S10Pilot S10ScalarControls
noncomputable section
theorem scalar_sign_and_scale :
    S10ScalarControls.coneValue (-1) 1 1 (![1, 0, 0] : Vec 3) SIGN 0 ∧
    S10ScalarControls.coneValue 2 1 1 (![1, 0, 0] : Vec 3) = 2 := by
  norm_num [S10ScalarControls.coneValue, normSq, dot, Fin.sum_univ_succ]
'''
    for name, sign, outcome in (('scalar_admissible_control', '<', 'PASS'),
                                ('negative_branch_reported_positive', '>', 'REJECTED')):
        path = SCRATCH / (name + '.lean')
        path.write_text(scalar.replace('SIGN', sign))
        record = run(path, outcome, {'scalar_sign_and_scale'} if outcome == 'REJECTED' else ())
        record['source'] = path.read_text()
        records.append(record)
    after = {str(p.relative_to(BASE)): digest(p) for p in sources}
    assert before == after
    result = {'status': 'PASS', 'scope': 'S10 action fidelity and normalization controls; no CAS bridge expansion',
              'instrument_sha256': digest(Path(__file__)), 'checks': records,
              'canonical_sha256_before': before, 'canonical_sha256_after': after}
    (BASE / '_measurements/S10_lean_contract_checks.json').write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    print('PASS: action/normalization controls; canonical sources unchanged.', flush=True)


if __name__ == '__main__':
    main()
