#!/usr/bin/env python3
"""Targeted source mutations; canonical Lean modules are never edited."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
LEAN = ROOT / 'lean'
SCRATCH = LEAN / 's10' / '_scratch' / 'q6_q7'
SCRATCH.mkdir(parents=True, exist_ok=True)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_mutation(name, module, original, replacement, required_theorems):
    source = LEAN / 's10' / 'S10Audit' / module
    before = digest(source)
    body = source.read_text()
    assert body.count(original) == 1, name
    mutated = body.replace(original, replacement, 1)
    target = SCRATCH / module
    target.write_text(mutated)
    command = ['lake', 'env', 'lean', '-DwarningAsError=true', str(target.relative_to(LEAN))]
    started = time.monotonic()
    run = subprocess.run(command, cwd=LEAN, text=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, timeout=180,
                         env={**os.environ, 'LAKE_CACHE_DIR': '.lake/cache'})
    elapsed = time.monotonic() - started
    log = SCRATCH / (name + '.log')
    log.write_text(run.stdout)
    assert run.returncode == 1, (name, run.returncode, run.stdout)
    assert digest(source) == before, name
    error_lines = [int(n) for n in re.findall(r'\.lean:(\d+):\d+: error:', run.stdout)]
    declarations = [(m.group(1), mutated[:m.start()].count('\n') + 1)
                    for m in re.finditer(r'^theorem (\w+)', mutated, re.M)]
    failed = []
    for index, (theorem, line) in enumerate(declarations):
        end = declarations[index + 1][1] if index + 1 < len(declarations) else len(mutated.splitlines()) + 1
        if any(line <= error < end for error in error_lines):
            failed.append(theorem)
    assert set(required_theorems) <= set(failed), (name, failed, run.stdout)
    return {'name': name, 'outcome': 'REJECTED', 'exit_status': run.returncode,
            'wall_seconds': round(elapsed, 3), 'source': str(source.relative_to(ROOT)),
            'canonical_sha256_before': before, 'canonical_sha256_after': digest(source),
            'mutation_sha256': digest(target), 'failed_theorems': failed,
            'diagnostic_error_lines': error_lines, 'command': command,
            'log_sha256': digest(log)}


def main():
    checks = [run_mutation(
        'double_epsilon_contraction', 'Curl.lean',
        '(leviCivitaSymbol ![i, j, k] : ℝ) * J j.succ k',
        '2 * (leviCivitaSymbol ![i, j, k] : ℝ) * J j.succ k',
        ['epsilonCurl_eq_jetCurl']),
        run_mutation(
            'drop_root_wavevector', 'RootDimensions.lean',
            '.mul (.div (coefficientTree p) (.atom (.coefficient .rho))) normKTree',
            '.mul (.div (coefficientTree p) (.atom (.coefficient .rho))) (.scalar 1)',
            ['ordinaryRootTree_eval', 'ordinaryRootTree_hasDim'])]
    output = ROOT / '_measurements' / 'S10_lean_q6_q7_checks.json'
    output.write_text(json.dumps({'status': 'PASS', 'scope': 'two isolated source mutations',
                                 'checks': checks}, indent=2) + '\n')
    print('PASS: both mutations rejected by the intended theorems; canonical hashes unchanged.')


if __name__ == '__main__':
    main()
