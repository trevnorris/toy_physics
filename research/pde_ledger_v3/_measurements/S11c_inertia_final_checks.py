"""Summarize actual emitted diagnostic values; no physical acceptance target."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import sys
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from ledger_fold import _restore

path = ROOT / 'scripts/out/S11c_d_transverse_sign_diagnostic.out'
records = {}
tags = []
for line in path.read_text().splitlines():
    if not line.startswith('PY_S11CD_SIGN_'):
        continue
    tag, payload = line.split(': ', 1)
    tags.append(tag)
    records[tag] = _restore(payload)


def named(obj, key):
    return next(v for k, v in obj if str(k) == key)


groups = {}
for prefix in ('ROW_MINUS_ACTION', 'CLOSED_MINUS_OPEN', 'CLOSED_SCALAR_ROWS',
               'THICKNESS_ROW_MINUS_ACTION', 'STANDING_WAVE_ENERGY_RATES',
               'S11B_ACTION_PURE_CURL_RESIDUAL', 'S11B_STORED_PURE_CURL_RESIDUAL',
               'S11CB_MINUS_S11B_STIFFNESS'):
    values = [named(obj, 'VALUE') for tag, obj in records.items()
              if tag.startswith('PY_S11CD_SIGN_' + prefix + '_')]
    groups[prefix] = {'count': len(values), 'values': sorted({sp.sstr(v) for v in values})}
frequency_pairs = [named(obj, 'VALUE') for tag, obj in records.items()
                   if tag.startswith('PY_S11CD_SIGN_FREQUENCY_SQUARED_')]
groups['FREQUENCY_SQUARED_PAIR_DIFFERENCES'] = {
    'count': len(frequency_pairs),
    'values': sorted({sp.sstr(sp.expand(v[0] - v[1])) for v in frequency_pairs})}
metadata = [named(obj, 'METADATA') for obj in records.values()
            if isinstance(obj, sp.Tuple) and len(obj) == 2
            and isinstance(obj[0], sp.Tuple) and str(obj[0][0]) == 'VALUE']
pins = records['PY_S11CD_SIGN_INPUT_DIGESTS']
result = {
    'tag_count': len(tags), 'duplicate_tags': [k for k, v in Counter(tags).items() if v > 1],
    'physical_record_count': len(metadata), 'groups': groups,
    'missing_leaf_units': sum(named(data, 'DIMENSION_L_T_M') is None for obj in metadata for _, data in obj),
    'dimension_constraints': sp.sstr(records['PY_S11CD_SIGN_DIMENSION_CONSTRAINT_RESIDUALS']),
    'dimension_unresolved': sp.sstr(records['PY_S11CD_SIGN_DIMENSION_UNRESOLVED']),
    'cache_flag': sp.sstr(records['PY_S11CD_SIGN_CACHED_DEVELOPMENT']),
    'input_pin_mismatches': [str(p) for p, digest in pins
                           if hashlib.sha256((ROOT / str(p)).read_bytes()).hexdigest() != str(digest)],
    'resource_measurements': sp.sstr(records['PY_S11CD_SIGN_RESOURCE_MEASUREMENTS']),
    'output_bytes': path.stat().st_size, 'output_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
}
print(json.dumps(result, indent=2))
