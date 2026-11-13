# S15 record decision list measurements

Reviewed baseline. The [setup lookups](S15_record_setup_lookups.md) give the input fingerprints, native expressions and measured comparison statuses. The fingerprint command below applies to the files as preserved at `5387b080`. Run from the repository root.

```sh
python3 - <<'PY'
from hashlib import sha256
from pathlib import Path
from subprocess import check_output
p = Path('research/pde_ledger_v3/directives')
for name in ('S15_record_decision_list.md', '_measurements/S15_record_setup_lookups.md'):
    print(name, sha256(check_output(['git', 'show', f'5387b080:{p / name}'])).hexdigest())
PY
```

```text
S15_record_decision_list.md 4d16e5346522b84c13187282c983d9d920601d046ded862a6a36f2e3c554c944
_measurements/S15_record_setup_lookups.md 3c39ab92b79d1c3160edae7f98c58003798ba704654dbb55e02ca7133f09dbdf
```

## Net shell quantities and final dispositions

The signed totals, flagged rechecks and elementary exponent boundaries can be retrieved directly. See the [original archive readout and the fold](S15_record_list_review_p1/README.md).

```sh
python3 - <<'PY'
import gzip, json, tarfile
from pathlib import Path
base = Path('research/pde_ledger_v3')
statuses = json.loads((base/'directives/_measurements/S15_separation_departure_wolfram_run_manifest.json').read_text())['item4_statuses']
with tarfile.open(base/'scripts/out/S15_reading_comparator_details.out') as archive:
    for i in (0, 1, 2, 19, 109, 288):
        def get(folder):
            return json.load(gzip.open(archive.extractfile(f'{folder}/{i:03d}.json.gz')))
        obj = get('objects')
        rows = {r['field']: r for r in obj['rows'] if r['role'] == 'cross'}
        if i < 3:
            print(i, rows['item3/object']['status'], len(rows['E_resp']['samples']), statuses[obj['object']])
        elif i == 288:
            print(i, rows['item3/exponent_branch_inventory']['status'], get('refinement')['branch']['same_elementary_boundaries'])
        else:
            field = next(f for f in get('refinement')['fields'] if f['field'] == 'item3/object')
            print(i, rows['item3/object']['status'], len(field['rechecks']), sorted({q['status'] for q in field['rechecks']}))
PY
```

```text
0 sample_agreement 30 "Not applicable to the general-field family integral."
1 sample_agreement 30 "Not applicable to the general-field family integral."
2 sample_agreement 30 "Not applicable to the general-field family integral."
19 sample_difference 40 ['within_precision_estimate']
109 sample_difference 12 ['within_precision_estimate']
288 unresolved True
```
