# Measurements — O2 record r2 review dispositions (generated 2026-10-08 00:01)

Generator: `_scratch/s9b_build/gen/o2_record_r2_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). Leg evidence is quoted verbatim from the Claude leg's filed stdout, copied from `/tmp/o2_record_review_r2/` to `_scratch/s9b_build/o2_record_review_r2_claude_evidence/`. In a structure line each role entry prints its role, then its unpaired reason, then its difference fields, in that order (layout shown under F1).

```
$ sha256sum -c _scratch/s9b_build/o2_record_review_baseline_r2.sha256
research/pde_ledger_v3/steps/O2_steady_brane_balance.md: OK
research/pde_ledger_v3/steps/_measurements/O2_record_measurements.md: OK
research/pde_ledger_v3/scripts/O2_record_measurements.py: OK
research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md: OK
```

## R3-1 — paired-role orientation differences beneath the energy balance
Token layout of a role entry (first entries of stream line 207):
```
$ sed -n 207p research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o '{"role":"[^"]*"\|"orientation":\[\[[^]]*\]\]\|"orientation":\[\]\|"unpaired_reason":null\|"unpaired_reason":"[^"]*"' | head -6
{"role":"OPEN_ApplyNativeFaceReduction"
"unpaired_reason":null
"orientation":[]
{"role":"OPEN_FaceApplicationVelocity_0"
"unpaired_reason":null
"orientation":[]
```

Every paired (reason null) role entry whose orientation difference is nonempty, over the whole stream, counted by stream line:
```
$ grep -n -o '{"role":"[^"]*"\|"orientation":\[\[[^]]*\]\]\|"orientation":\[\]\|"unpaired_reason":null\|"unpaired_reason":"[^"]*"' research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -A1 ':"unpaired_reason":null' | grep 'orientation":\[\[' | cut -d: -f1 | sort | uniq -c
      4 207
```

```
$ sed -n 207p research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o '{"role":"[^"]*"\|"orientation":\[\[[^]]*\]\]\|"orientation":\[\]\|"unpaired_reason":null\|"unpaired_reason":"[^"]*"' | grep -B2 'orientation":\[\[' | grep -A2 '{"role"' | grep -A1 -B1 'reason":null'
{"role":"OPEN_MaterialEnergyDensity"
"unpaired_reason":null
"orientation":[["argument",1]]
{"role":"OPEN_MaterialEnergyFlux_0"
"unpaired_reason":null
"orientation":[["argument",1]]
{"role":"OPEN_MaterialEnergyFlux_1"
"unpaired_reason":null
"orientation":[["argument",1]]
{"role":"OPEN_MaterialEnergyFlux_2"
"unpaired_reason":null
"orientation":[["argument",1]]
```

The record (r2):
```
$ sed -n '188,189p;195p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
| Component / comparison stream line | PY / WL entries | Role and orientation differences | Remaining inventory differences |
|---|---:|---|---|
| `energy_balance` / 657 | 5 / 5 | All empty; no `OPEN_free` entry. | Named differences on all five role groups; live differences on all three energy-flux roles. |
```

```
$ sed -n '235,238p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
These are differences in the action inventories beneath the energy balance entries. They remain
OPEN with no owner named; the balance-entry role/orientation/OPEN-free multiset match in §3 does
not establish a match of the actions beneath those entries. No occurrence-level pairing, complete
work equality or duplicate-power conclusion is formed from either engine's one-sided content.
```

```
$ sed -n '432,434p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
In particular, the duplicate-power condition inherits the one-sided `OPEN_MaterialCompatibility`
actions in all four energy rows and the one-sided material-energy density/flux actions in
`energy_power` (§4; M6). A match of balance-entry classifications does not settle whether the complete
```

```
$ for t in 'JointPowerAccounting' '"argument"' 'argument",1' 'occurs twice'; do printf '%s: ' "$t"; grep -c -F "$t" research/pde_ledger_v3/steps/O2_steady_brane_balance.md; done
JointPowerAccounting: 0
"argument": 0
argument",1: 0
occurs twice: 0
```

The inheritance sentence is new in r2; repair brief 2 named the one-sided class's roles:
```
$ git show 69f0ed3d:research/pde_ledger_v3/steps/O2_steady_brane_balance.md | grep -c -F 'duplicate-power condition inherits'
0
```

```
$ grep -n -F 'OPEN_MaterialEnergyDensity' _scratch/s9b_build/o2_record_repair2_prompt.md
19:   in `energy_power`, `OPEN_MaterialEnergyDensity` and `OPEN_MaterialEnergyFlux_0/1/2`. The record has no match for
```

Claude leg stdout (verbatim):
```
$ cat _scratch/s9b_build/o2_record_review_r2_claude_evidence/orientation_diffs.stdout | grep 'field=orientation'
structure line=207 row=energy_balance role=OPEN_MaterialEnergyDensity field=orientation delta(py-wl)=[['argument', 1]] occ py=[{'argument': 1}, {'1': 1}] wl=[{'1': 1}]
structure line=207 row=energy_balance role=OPEN_MaterialEnergyFlux_0 field=orientation delta(py-wl)=[['argument', 1]] occ py=[{'argument': 1}, {'1': 1}] wl=[{'1': 1}]
structure line=207 row=energy_balance role=OPEN_MaterialEnergyFlux_1 field=orientation delta(py-wl)=[['argument', 1]] occ py=[{'argument': 1}, {'1': 1}] wl=[{'1': 1}]
structure line=207 row=energy_balance role=OPEN_MaterialEnergyFlux_2 field=orientation delta(py-wl)=[['argument', 1]] occ py=[{'argument': 1}, {'1': 1}] wl=[{'1': 1}]
```

```
$ grep -A12 'line=206 row=energy_balance engine=py' _scratch/s9b_build/o2_record_review_r2_claude_evidence/energy_nesting.stdout
line=206 row=energy_balance engine=py
   1x OPEN_MaterialCompatibility  enclosed-by: OPEN_JointPowerAccounting
   1x OPEN_MaterialCompatibility  enclosed-by: OPEN_JointPowerAccounting > OPEN_MaterialEnergyDensity
   1x OPEN_MaterialCompatibility  enclosed-by: OPEN_JointPowerAccounting > OPEN_MaterialEnergyFlux_0
   1x OPEN_MaterialCompatibility  enclosed-by: OPEN_JointPowerAccounting > OPEN_MaterialEnergyFlux_1
   1x OPEN_MaterialCompatibility  enclosed-by: OPEN_JointPowerAccounting > OPEN_MaterialEnergyFlux_2
   1x OPEN_MaterialCompatibility  enclosed-by: OPEN_MaterialEnergyDensity
   1x OPEN_MaterialCompatibility  enclosed-by: OPEN_MaterialEnergyFlux_0
   1x OPEN_MaterialCompatibility  enclosed-by: OPEN_MaterialEnergyFlux_1
   1x OPEN_MaterialCompatibility  enclosed-by: OPEN_MaterialEnergyFlux_2
   1x OPEN_MaterialEnergyDensity  enclosed-by: OPEN_JointPowerAccounting
   1x OPEN_MaterialEnergyDensity  enclosed-by: TOP
   1x OPEN_MaterialEnergyFlux_0  enclosed-by: OPEN_JointPowerAccounting
```

## R3-2 — the WL-only xi_w'' difference: the source's assignment and the record's owner statement
```
$ sed -n '81,83p' research/pde_ledger_v3/directives/_measurements/O2_comparator_build_r5_review_disposition.md
Its FORM ablation that forgets derivative order makes this difference disappear (135 leaf paths) and fails two
tests. This is a cross-engine difference in OPEN content. Interpreting it belongs to the O2 record (sub-step 7),
under M1: it is preserved, never designed away.
```

```
$ grep -n -F 'Interpreting it belongs to the O2 record' research/pde_ledger_v3/steps/_measurements/O2_record_measurements.md
1713:82: tests. This is a cross-engine difference in OPEN content. Interpreting it belongs to the O2 record (sub-step 7),
```

```
$ sed -n '256,260p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
These are OPEN cross-engine content differences, with no owner named for adjudication. The comparator
acceptance specifically calls the WL-only `xi_w''(r)` a **“cross-engine difference in OPEN content”**,
with filed reviewer measurements `{{1, 43}, {2, 12}}` for WL XiW derivatives versus `{1: 105}` for
PY `OPEN_MomentumFlux_*` occurrences (M1, acceptance's measurement section). These are quotations
of that evidence, not a new count or computation. **OPEN; no owner named for reconciliation.**
```

```
$ sed -n '450,454p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
The engine acceptance describes **“Claude's seven named representational differences”** and routes
them to the comparator. The comparator acceptance describes the second-derivative result as
**“a cross-engine difference in OPEN content.”** Both statements and scopes are retained (M1);
the former cannot settle all later printed inventory differences as merely representational.
Their complete physical reconciliation remains OPEN under M1, with no adjudication owner named.
```

```
$ grep -c -F 'Interpreting it belongs' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
0
```

Present since r0 (243a530f) and in r1 (69f0ed3d):
```
$ for c in 243a530f 69f0ed3d; do printf '%s ' $c; git show $c:research/pde_ledger_v3/steps/O2_steady_brane_balance.md | grep -n -F 'no owner named for reconciliation'; done
243a530f 185:of that evidence, not a new count or computation. **OPEN; no owner named for reconciliation.**
69f0ed3d 243:of that evidence, not a new count or computation. **OPEN; no owner named for reconciliation.**
```

The record directive's retrieval bound:
```
$ sed -n '54,56p' research/pde_ledger_v3/directives/O2_record_directive.md
   - Only retrieval is allowed: existence, verbatim retrieval, literal-match counts, and the shape of a named
     stored object.
   - A question that needs computation beyond retrieval is listed as open, not computed.
```

Claude leg stdout (verbatim):
```
$ cat _scratch/s9b_build/o2_record_review_r2_claude_evidence/xiw_sanity.stdout
('balance_comparison', 'order2') 27
('balance_entries', 'order1') 264
('balance_entries', 'order2') 81
('structure', 'order1') 2028
```

