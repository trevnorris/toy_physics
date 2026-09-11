#!/usr/bin/env python3
"""Read stored branch payloads and source passages; no spectral computation."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import sympy as sp
import ledger_fold


def pin(path):
    content = (ROOT / path).read_bytes()
    return {'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()}


def main():
    paths = [
        'directives/S11b_SHARED_PHYSICS.md',
        *('directives/S11c_' + stage + '_SHARED_PHYSICS.md'
          for stage in ('a', 'b', 'c1', 'c2', 'd')),
        'directives/S11c_d_sympy_build_PROGRAM_BRIEF.md',
        'directives/S11c_d_sympy_build_directive.md',
        'scripts/S11b_interface_coupling_law_sympy_audit.py',
        'scripts/S11c_a_interface_geometry_sympy_audit.py',
        'scripts/S11c_b_brane_operator_sympy_audit.py',
        'scripts/S11c_c1_bulk_closure_sympy_audit.py',
        'scripts/S11c_c2_selfenergy_fold_sympy_audit.py',
        'scripts/S11c_d_mixing_scattering_sympy_audit.py',
        'scripts/S11c_d_sheet_continuation_probe.py', 'scripts/ledger_fold.py',
        'scripts/S11b_exports.py',
        *('scripts/S11c_' + stage + '_exports.py' for stage in ('a', 'b', 'c1', 'c2')),
        '_measurements/S11c_d_channel_preflight_input.json',
    ]
    record = {
        'command': 'python3 _measurements/S11c_d_sheet_scope_lookup.py',
        'cwd': str(ROOT),
        'checkpoint': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'files': {path: pin(path) for path in paths},
    }
    # Exact source retrieval. These excerpts do not decide mathematical validity.
    excerpts = {
        'directives/S11b_SHARED_PHYSICS.md': [(84, 94), (113, 158)],
        'directives/S11c_c1_SHARED_PHYSICS.md': [(88, 122)],
        'directives/S11c_d_SHARED_PHYSICS.md': [(354, 380), (401, 414), (495, 524)],
        'scripts/S11c_d_mixing_scattering_sympy_audit.py':
            [(715, 732), (1626, 1630), (1778, 1796), (1948, 1956)],
        'scripts/S11c_c2_selfenergy_fold_sympy_audit.py': [(437, 450), (751, 758)],
        'scripts/S11b_interface_coupling_law_sympy_audit.py':
            [(1867, 1888), (1953, 1964)],
    }
    record['source_excerpts'] = []
    for path, ranges in excerpts.items():
        lines = (ROOT / path).read_text().splitlines()
        for first, last in ranges:
            record['source_excerpts'].append({
                'path': path, 'first_line': first, 'last_line': last,
                'lines': lines[first-1:last],
            })
    fold, audit = ledger_fold.load_model(*(
        ROOT / ('scripts/S11c_' + stage + '_exports.py')
        for stage in ('b', 'c1', 'c2')))
    record['fold_audit'] = audit
    record['closed_rows'] = []
    for key in ('s11cc2ClosedSlabOperator', 's11cc2ClosedCouplingKernel'):
        row = fold[key]
        cases = []
        for case, payload in row['value']:
            slots = {str(k): v for k, v in payload}
            bindings = slots['COMPUTED_BRANCH_BINDINGS']
            cases.append({
                'case': list(map(str, case)), 'slot_names': list(slots),
                'branch_bindings': [sp.srepr(v) for v in bindings],
                'branch_symbol_assumptions': {
                    s.name: s.assumptions0 for s in sorted(
                        bindings.atoms(sp.Symbol), key=str)},
            })
        record['closed_rows'].append({
            'key': key, 'step': row['step'], 'class': row['class'], 'cases': cases})
    # The roots are stored objects, retrieved independently of their producer code.
    inherited = ledger_fold._ledger_from_path(ROOT / 'scripts/S11b_exports.py')
    record['s11b_stored_root_rows'] = {
        k: {'step': inherited[k]['step'], 'value_srepr': sp.srepr(inherited[k]['value'])}
        for k in ('roots', 'sheet_of_each_root')
    }
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
