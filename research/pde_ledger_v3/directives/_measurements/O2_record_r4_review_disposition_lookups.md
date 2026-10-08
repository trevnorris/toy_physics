# Measurements — O2 record r4 review dispositions (generated 2026-10-08 10:06)

Generator: `_scratch/s9b_build/gen/o2_record_r4_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters).

```
$ sha256sum -c _scratch/s9b_build/o2_record_review_baseline_r4.sha256
research/pde_ledger_v3/steps/O2_steady_brane_balance.md: OK
research/pde_ledger_v3/steps/_measurements/O2_record_measurements.md: OK
research/pde_ledger_v3/scripts/O2_record_measurements.py: OK
research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md: OK
```

## Verdicts
```
$ grep -n -o '\*\*Verdict: [A-Z]*' _scratch/s9b_build/o2_record_review_r4_grok.txt
3:**Verdict: CLEAR
```

```
$ grep -n -m1 'Verdict' _scratch/s9b_build/o2_record_review_r4_claude.md
1:**Verdict: CLEAR**
```

## R4-1 (r3) — the unioned OPEN_free group and the multiset provenance in r4
```
$ grep -n 'hold_normal` / 654' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
208:| `hold_normal` / 654 | 32 / 32 | `[]` at this unioned level; **one** all-empty `OPEN_free` group comparison. M7 separately counts **3** stored `OPEN_free` tuples per engine. | Named differences on all 26 OPEN role groups; live differences on all twelve momentum-flux roles. |
```

```
$ grep -c -F 'printed role/orientation/OPEN-free multiset' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
0
```

```
$ sed -n '54p;192,193p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
(M6), literal counts of joint stored role/orientation/OPEN-free classification tuples, separately
For **each** balance component, the record's literal counts of the joint stored
`(role, orientation, open_free)` tuples match between PY and WL **at the level of those stored
```

## Note (r3) — the acceptance's W1 text now carried in the measurements file
```
$ grep -n -F 'W1 resolved' research/pde_ledger_v3/steps/_measurements/O2_record_measurements.md
1693:44: **W1 resolved.** In the production output, `"wl::6"` and `"wl::3,8"` now occur on 0 lines (r4: 6 each).
```

## Claude leg notes below the filter — what the record claims about them
```
$ sed -n '178,182p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
The complete printed not-formed reason list, with literal counts in the joined residual rows (M4/M8),
is **“held, OPEN or unsupported applied head” (92)**, **“named operand has no supplied scalar value”
(34)**, **“container structure differs” (31)**, **“nested sibling absent” (27)**, **“text, native
boolean or name is not a subtractable operand” (21)**, and **“bindings alone are not value evidence”
(7)**. Unformed coordinate/basis leaves and all relation/container results remain in M8. The
```

```
$ grep -n -F 'every named operand and complete live-object dependence printed by either engine' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
469:Part D must retain **every named operand and complete live-object dependence printed by either engine**,
```

```
$ grep -n -F 'undischarged' research/pde_ledger_v3/steps/O2_steady_brane_balance.md | head -4
316:**The interpretation remains an undischarged sub-step-7 obligation, returned to Claude (orchestrator)
364:for `xi_w''` carries §4's undischarged/returned status. `wl::D`/`wl::Map` raw-census spellings
376:| Transport calculus | `momentum_storage`, `momentum_transport` and the three hold representations retain OPEN content and unformed full residuals. Explicit profile/velocity gradients print exact scalar-leaf residuals `0`; the six balance closed-part residuals are `0` at that level only (§3; M8). The WL-only derivative cont
466:| `ℬ_hold^live` in-plane/bulk/graph-normal components and `ℬ_E^steady` | Conditional named accounting. The record's literal counts match **at the level of joint stored `(role, orientation, open_free)` classification tuple counts only** (M7). Separately, the comparator prints empty head/role/orientation differences **at t
```

