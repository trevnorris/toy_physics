"""Mechanical artifact and serialized-row census for the inertia rebuild."""
from pathlib import Path
import argparse
import ast
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def export_data(path):
    assignments = {}
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            assignments[node.targets[0].id] = node.value
    ledger = assignments['_LEDGER']
    values = {}
    for key, row in zip(ledger.keys, ledger.values):
        fields = {ast.literal_eval(k): v for k, v in zip(row.keys, row.values)}
        value = fields['value']
        values[ast.literal_eval(key)] = ast.literal_eval(value.args[0])
    digests = ast.literal_eval(assignments['BUILD_INPUT_DIGESTS'].args[0])
    imports = ast.literal_eval(assignments['IMPORT_KEYS']) if 'IMPORT_KEYS' in assignments else ()
    return values, digests, imports


def locate(name):
    options = (ROOT / name, ROOT / 'scripts' / name, ROOT / 'directives' / name)
    return next(p for p in options if p.is_file())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('b', 'c1', 'c2'))
    parser.add_argument('--baseline', type=Path, default=Path('/tmp/s11c-inertia-repair-baseline'))
    args = parser.parse_args()
    relative = Path('scripts') / ('S11c_' + args.stage + '_exports.py')
    new, pins, imports = export_data(ROOT / relative)
    old, old_pins, _ = export_data(args.baseline / relative)
    changed = sorted(key for key in new.keys() & old.keys() if new[key] != old[key])
    report = {
        'stage': args.stage,
        'export_sha256': hashlib.sha256((ROOT / relative).read_bytes()).hexdigest(),
        'export_bytes': (ROOT / relative).stat().st_size,
        'rows': len(new), 'added_keys': sorted(new.keys() - old.keys()),
        'removed_keys': sorted(old.keys() - new.keys()), 'changed_value_serializations': changed,
        'pin_mismatches': {name: {'pinned': digest, 'actual': hashlib.sha256(locate(name).read_bytes()).hexdigest()}
                           for name, digest in pins.items()
                           if hashlib.sha256(locate(name).read_bytes()).hexdigest() != digest},
        'pins': pins,
    }
    if args.stage == 'b':
        c1_old, _, c1_imports = export_data(args.baseline / 'scripts/S11c_c1_exports.py')
        report['changed_c1_direct_input_serializations'] = sorted(set(c1_imports) & set(changed))
    print(json.dumps(report, indent=2), flush=True)
