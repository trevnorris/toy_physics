# S15 record decision list, independent local reader A

**Verdict. Fold two clarifications before handing the S15 record decision list to the author for items 1–3, 5 and 6.** This verdict covers the list and companion lookups on `41adb784`. It gives no verdict on the future record or register pass. Candidate SHA-256 prefixes are `4d16e5346522b84c` (list) and `3c39ab92b79d1c31` (lookups). I used local source retrieval only.

1. **Identify the actual coverage of the proposed S21 observable.** The list at lines 11–17 correctly defines the family test and requests both a closed native term and the signed total. The example at lookup lines 91–97 is one inner, drain-2 centre-transport term. All three **total** crossing components instead retain a general-field integral on the power-law family, have only sampled cross-engine agreement for that item-3 result, and fall outside Wolfram’s item-4 classification domain. Completing the deferred classification of eligible closed pieces alone would still provide no classification of the signed total. Ask the author to state this explicitly in the record and S21 entry. Conditions on a closed piece and cancellation among pieces can be recorded at their measured scope. The test on the total allows both geometric and field-response contributions. Sources: `directives/S15_reading_decision_list.md:65–86,103–113`; committed Wolfram manifest, `item4_statuses`; paired native and archive results below.
2. **Make the 336-object index recover the final dispositions as well as the initial ones.** List line 15 and lookup lines 65–75 could induce the author to inventory the `status` in each `objects/NNN` row. That field is the initial status. For objects 019 and 109 the field reads `sample_difference`; all 40 and 12 flagged points, respectively, were within the later precision estimates. Object 288 is one of six whose elementary exponent boundaries were matched in a later pass. The committed exceptional-selection inventory covers 48 sampled emitted family results, not all 336 objects and fields. Require the index to associate each affected `(object, field)` with its refinement and log-branch evidence where applicable, preserving the original observations and the final qualified outcome. Sources: `scripts/S15_reading_report.py:187–245`; `scripts/out/S15_reading_comparator.out:103–120`; `directives/_measurements/S15_reading_comparator_coverage.json:4–10`. This is retrieval from the accepted calculation.

**Other checks.** Both lookup commands reproduced their complete displayed blocks (five and thirteen output lines). The shell signs, live fields, pressure domains, velocity-difference definition, regularity qualifications, dimensional-evidence distinctions and S16 face/core/inertia handoff are supported by the original decisions, accepted reading list and comparator dispositions. For later comparison with a companion-induced free acceleration, the register also distinguishes the companion’s field and the isolated baseline (`SUBSTRATE_REQUIREMENTS.md:1011–1026,2074–2082`). The finished record may cite the stopped Singular trial from its current November 12 summary (`directives/_measurements/S15_separation_departure_singular_run_record.md:3–15`, 13 of 14 trials returned conditions; survey never started). That source’s later November 10 details describe an earlier stop.

## Retrieval

Run from the repository root. Archive indices 000–002 are `both_surfaces/sum_crossing/component_1` through `_3`; 019, 109 and 288 are identified by the corresponding source rows above.

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

```sh
python3 - <<'PY'
import re, subprocess
from pathlib import Path
s = Path('research/pde_ledger_v3/directives/_measurements/S15_record_setup_lookups.md').read_text()
for cmd, expect in zip(re.findall(r'```sh\n(.*?)\n```', s, re.S), re.findall(r'```text\n(.*?)\n```', s, re.S), strict=True):
    out = subprocess.check_output(['bash', '-c', cmd], text=True)
    print(len(out.splitlines()), out == expect + '\n')
PY
```

```text
5 True
13 True
```
