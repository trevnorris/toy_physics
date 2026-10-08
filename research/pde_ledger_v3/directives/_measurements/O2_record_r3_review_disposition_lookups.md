# Measurements — O2 record r3 review dispositions (generated 2026-10-08 09:19)

Generator: `_scratch/s9b_build/gen/o2_record_r3_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). Leg evidence is quoted verbatim from the Claude leg's filed stdout, copied from `/tmp/o2_r3_review_claude/` to `_scratch/s9b_build/o2_record_review_r3_claude_evidence/`.

```
$ sha256sum -c _scratch/s9b_build/o2_record_review_baseline_r3.sha256
research/pde_ledger_v3/steps/O2_steady_brane_balance.md: OK
research/pde_ledger_v3/steps/_measurements/O2_record_measurements.md: OK
research/pde_ledger_v3/scripts/O2_record_measurements.py: OK
research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md: OK
```

## R4-1 — hold_normal's repeated OPEN_free entries: what the comparator compares, and what the record says
The comparator unions object fields across same-role entries:
```
$ sed -n '999p;1087,1093p' research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
OBJECT_FIELDS = frozenset(('head','role','orientation','named_OPEN_operands','live_arguments'))
def merge_fields(current,fields):
    for field,values in fields.items():
        inventory=current.setdefault(field,Counter())
        if field in OBJECT_FIELDS:
            inventory.update({key:1 for key in values if key not in inventory})
        else:
            inventory.update(values)
```

Stored entries (balance_entries, line 653) versus grouped differences (balance_comparison, line 654):
```
$ printf 'line 653 "role":"OPEN_free" count: '; sed -n 653p research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o '"role":"OPEN_free"' | wc -l; printf 'line 653 "open_free":true count: '; sed -n 653p research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o '"open_free":true' | wc -l; printf 'line 654 "OPEN_free" difference keys: '; sed -n 654p research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o '"OPEN_free":{' | wc -l
line 653 "role":"OPEN_free" count: 6
line 653 "open_free":true count: 6
line 654 "OPEN_free" difference keys: 1
```

```
$ sed -n 654p research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o '"OPEN_free":{[^}]*}'
"OPEN_free":{"head":[],"live_arguments":[],"named_OPEN_operands":[],"orientation":[],"role":[]}
```

The printed scope and the acceptance's N4 condition:
```
$ grep -n -o 'how many times an object occurs' research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out
1:how many times an object occurs
```

```
$ sed -n '69p;88p' research/pde_ledger_v3/directives/_measurements/O2_comparator_build_r5_review_disposition.md
| N4 | Inventories are unioned across same-role occurrences within a row. Component-indexed roles and per-component balances limit the effect. (Claude N4) | **Note; carried to the record.** This is a form of the stated multiplicity limit. The record does not read an empty difference on a row with repeated same-role occurrences a
claimed from its output, provided the record carries the limits above (N1, N4, the generalized-rates scope) together
```

The record (r3):
```
$ sed -n '201p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
| `hold_normal` / 654 | 32 / 32 | `[]` at balance-entry level; three `OPEN_free` inventory matches. | Named differences on all 26 OPEN role groups; live differences on all twelve momentum-flux roles. |
```

```
$ grep -n -F 'printed role/orientation/OPEN-free multiset match' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
455:| `ℬ_hold^live` in-plane/bulk/graph-normal components and `ℬ_E^steady` | Conditional named accounting. The printed role/orientation/OPEN-free multiset match holds **at the balance-entry classification level only**, alongside residual `0` at the six balance closed parts. Every inventory difference persists: §4 includes t
```

```
$ grep -n -F 'Inventories union repeated same-role' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
347:on the measured sign paths; no broader validation follows. Inventories union repeated same-role
```

History: the per-entry wording since r1, the 'printed' attribution since r2; the r0 disposition's own wording:
```
$ for c in 69f0ed3d 9eff4ae5; do printf '%s ' $c; git show $c:research/pde_ledger_v3/steps/O2_steady_brane_balance.md | grep -n 'hold_normal` / 654' | cut -c1-140; done
69f0ed3d 194:| `hold_normal` / 654 | 32 / 32 | All empty; matching three `OPEN_free` entries. | Named differences on all 26 OPEN role groups; live di
9eff4ae5 194:| `hold_normal` / 654 | 32 / 32 | All empty; matching three `OPEN_free` entries. | Named differences on all 26 OPEN role groups; live di
```

```
$ for c in 69f0ed3d 9eff4ae5; do printf '%s ' $c; git show $c:research/pde_ledger_v3/steps/O2_steady_brane_balance.md | grep -c -F 'printed role/orientation/OPEN-free multiset'; done
69f0ed3d 0
9eff4ae5 1
```

```
$ grep -n -o 'These are printed comparator fields[^.]*' research/pde_ledger_v3/directives/_measurements/O2_record_r0_review_disposition.md
41:These are printed comparator fields, which the record may report by retrieval
```

Claude leg stdout (verbatim):
```
$ cat _scratch/s9b_build/o2_record_review_r3_claude_evidence/hold_normal_entries.stdout
line 643 hold_inplane[0]
  py: entries=8 distinct=8 repeated={}
  wl: entries=8 distinct=8 repeated={}
  grouped entry_differences keys: 8 OPEN_free key present: True
line 645 hold_inplane[1]
  py: entries=8 distinct=8 repeated={}
  wl: entries=8 distinct=8 repeated={}
  grouped entry_differences keys: 8 OPEN_free key present: True
line 647 hold_inplane[2]
  py: entries=8 distinct=8 repeated={}
  wl: entries=8 distinct=8 repeated={}
  grouped entry_differences keys: 8 OPEN_free key present: True
line 650 hold_bulk[]
  py: entries=8 distinct=8 repeated={}
  wl: entries=8 distinct=8 repeated={}
  grouped entry_differences keys: 8 OPEN_free key present: False
line 653 hold_normal[]
  py: entries=32 distinct=27 repeated={('"OPEN_free"', True): 3, ('"[\\"OPEN_SumOverAllNativeFaces\\"]"', False): 4}
  wl: entries=32 distinct=27 repeated={('"[\\"OPEN_SumOverAllNativeFaces\\"]"', False): 4, ('"OPEN_free"', True): 3}
  grouped entry_differences keys: 27 OPEN_free key present: True
line 656 energy_balance[]
  py: entries=5 distinct=5 repeated={}
  wl: entries=5 distinct=5 repeated={}
  grouped entry_differences keys: 5 OPEN_free key present: False
```

```
$ grep -n 'OPEN_free' _scratch/s9b_build/o2_record_review_r3_claude_evidence/open_free_entries.stdout | head -4
1:653 py entry 0 role OPEN_free orientation -1 inventories {"named_OPEN_operands": [], "live_arguments": ["[\"LiveProfile\",\"V_r\",[],[[\"Pow\",\"\",[],[[\"Add\",\"\",[],[[\"Pow\",\"\",[],[[\"Name\",\"x1\",[],[]],[\"Number\",\"2\",[],[]]]],[\"Pow\",\"\",[],[[\"Name\",\"x2\",[],[]],[\"Number\",\"2\",[],[]]]],[\"Pow\",\"\",[],[[\
2:653 py entry 8 role OPEN_free orientation -1 inventories {"named_OPEN_operands": [], "live_arguments": ["[\"LiveProfile\",\"V_r\",[],[[\"Pow\",\"\",[],[[\"Add\",\"\",[],[[\"Pow\",\"\",[],[[\"Name\",\"x1\",[],[]],[\"Number\",\"2\",[],[]]]],[\"Pow\",\"\",[],[[\"Name\",\"x2\",[],[]],[\"Number\",\"2\",[],[]]]],[\"Pow\",\"\",[],[[\
3:653 py entry 16 role OPEN_free orientation -1 inventories {"named_OPEN_operands": [], "live_arguments": ["[\"LiveProfile\",\"V_r\",[],[[\"Pow\",\"\",[],[[\"Add\",\"\",[],[[\"Pow\",\"\",[],[[\"Name\",\"x1\",[],[]],[\"Number\",\"2\",[],[]]]],[\"Pow\",\"\",[],[[\"Name\",\"x2\",[],[]],[\"Number\",\"2\",[],[]]]],[\"Pow\",\"\",[],[[
4:653 wl entry 14 role OPEN_free orientation -1 inventories {"named_OPEN_operands": [], "live_arguments": ["[\"LiveProfile\",\"V_r\",[],[[\"Pow\",\"\",[],[[\"Add\",\"\",[],[[\"Pow\",\"\",[],[[\"Name\",\"x1\",[],[]],[\"Number\",\"2\",[],[]]]],[\"Pow\",\"\",[],[[\"Name\",\"x2\",[],[]],[\"Number\",\"2\",[],[]]]],[\"Pow\",\"\",[],[[
```

## Note (below the filter) — an attribution to the acceptance whose text the measurements file does not carry
```
$ sed -n '353,354p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
for `xi_w''` carries §4's undischarged/returned status. `wl::D`/`wl::Map` raw-census spellings
remain diagnostics, as the acceptance states. No trace string, equal bare name, or unformable
```

```
$ grep -n -F 'SOURCE directives/_measurements/O2_comparator_build_r5_review_disposition.md' research/pde_ledger_v3/steps/_measurements/O2_record_measurements.md
1680:SOURCE directives/_measurements/O2_comparator_build_r5_review_disposition.md
```

```
$ sed -n '1680,1800p' research/pde_ledger_v3/steps/_measurements/O2_record_measurements.md | grep -o '^[0-9]*:' | tr '\n' ' '
1: 2: 3: 4: 5: 6: 7: 8: 9: 10: 11: 12: 62: 63: 64: 65: 66: 67: 68: 69: 70: 71: 72: 73: 74: 75: 76: 77: 78: 79: 80: 81: 82: 83: 84: 85: 86: 87: 88: 89: 90: 91: 92: 1: 2: 3: 4: 5: 6: 7: 8: 9: 10: 11: 12: 13: 14: 15: 16: 17: 18: 19: 20: 21: 22: 23: 24: 25: 26: 50: 51: 52: 53: 54: 55: 56: 57: 58: 59: 60: 1100: 1101: 1102: 1103: 1104
```

```
$ grep -n -F 'W1 resolved' research/pde_ledger_v3/directives/_measurements/O2_comparator_build_r5_review_disposition.md
44:**W1 resolved.** In the production output, `"wl::6"` and `"wl::3,8"` now occur on 0 lines (r4: 6 each).
```

