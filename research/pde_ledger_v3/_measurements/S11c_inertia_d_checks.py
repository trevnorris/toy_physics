"""Inventory the completed prototype transcript and its recorded input bindings."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import sys
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from ledger_fold import _restore

path = ROOT / 'scripts/out/S11c_d_mixing_scattering_sympy_audit.out'
selected = {
    'BUILD_INPUT_DIGESTS', 'REDUCED_DIMENSION_CONSTRAINT_RESIDUALS',
    'REDUCED_DIMENSION_UNRESOLVED', 'PENCIL_DIMENSION_CONSTRAINT_RESIDUALS',
    'IMPLEMENTATION_CHECKPOINT', 'OUTSTANDING_CONSTRUCTIONS',
    'RESOURCE_MEASUREMENTS', 'PROCESS_COMPLETION',
}
tags, records = [], {}
for line in path.open():
    if not line.startswith('PY_S11CD_'):
        continue
    tag, payload = line.split(': ', 1)
    tags.append(tag)
    if tag.removeprefix('PY_S11CD_') in selected:
        records[tag.removeprefix('PY_S11CD_')] = _restore(payload)
pins = records['BUILD_INPUT_DIGESTS']
result = {
    'tag_count': len(tags), 'duplicate_tags': [k for k,v in Counter(tags).items() if v > 1],
    'missing_final_records': sorted(selected - records.keys()),
    'reduced_case_records': [tag.removeprefix('PY_S11CD_REDUCED_ACTION_CENSUS_') for tag in tags
                             if tag.startswith('PY_S11CD_REDUCED_ACTION_CENSUS_')],
    'spectral_pit_sample_count': sum(tag.startswith('PY_S11CD_SPECTRAL_PIT_CARRIER_VALUES_') for tag in tags),
    'dimension_records': {key: sp.sstr(records[key]) for key in selected if 'DIMENSION' in key},
    'outstanding_constructions': list(map(str, records['OUTSTANDING_CONSTRUCTIONS'])),
    'resource_measurements': sp.sstr(records['RESOURCE_MEASUREMENTS']),
    'input_pin_mismatches': [str(p) for p,digest in pins
                           if hashlib.sha256((ROOT / str(p)).read_bytes()).hexdigest() != str(digest)],
    'output_bytes': path.stat().st_size,
    'output_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
}
print(json.dumps(result, indent=2))
