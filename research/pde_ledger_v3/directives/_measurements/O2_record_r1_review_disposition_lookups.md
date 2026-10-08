# Measurements — O2 record r1 review dispositions (generated 2026-10-07 23:28)

Generator: `_scratch/s9b_build/gen/o2_record_r1_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). Leg evidence is quoted verbatim from the Claude leg's filed stdout, copied from `/tmp/o2_record_review_r1/` to `_scratch/s9b_build/o2_record_review_r1_claude_evidence/`.

```
$ sha256sum -c _scratch/s9b_build/o2_record_review_baseline_r1.sha256
research/pde_ledger_v3/steps/O2_steady_brane_balance.md: OK
research/pde_ledger_v3/steps/_measurements/O2_record_measurements.md: OK
research/pde_ledger_v3/scripts/O2_record_measurements.py: OK
research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md: OK
```

The Claude leg did not read the baseline file (its wrapper restricted _scratch reads); the shas it measured itself, against the baseline:
```
$ grep -A6 '| Artifact | sha256 |' _scratch/s9b_build/o2_record_review_r1_claude.md
| Artifact | sha256 |
|---|---|
| record `steps/O2_steady_brane_balance.md` | `8d64a476…` |
| `steps/_measurements/O2_record_measurements.md` | `37293a05…` |
| generator `scripts/O2_record_measurements.py` | `7be41b28…` |
| `SUBSTRATE_REQUIREMENTS.md` | `b0b3a690…` |

```

```
$ cut -c1-8,65- _scratch/s9b_build/o2_record_review_baseline_r1.sha256
8d64a476  research/pde_ledger_v3/steps/O2_steady_brane_balance.md
37293a05  research/pde_ledger_v3/steps/_measurements/O2_record_measurements.md
7be41b28  research/pde_ledger_v3/scripts/O2_record_measurements.py
b0b3a690  research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
```

Grok's report on the same eight entries:
```
$ grep -n -o 'OPEN action records 430.*\|.\{0,40\}8 one-sided declared roles.\{0,200\}' _scratch/s9b_build/o2_record_review_r1_grok.txt
26:OPEN action records 430 = 234 both-sides content deltas
51:_policy`. The 188 unbound roles and the 8 one-sided declared roles (`energy_storage` / `energy_transport` / `energy_power` / `energy_balance`, missing `OPEN_MaterialCompatibility` or energy-flux roles) stay OPEN in §4, with no owner assigned. Nothing in that set is 
```

```
$ grep -n 'declared-role missing one side' _scratch/s9b_build/o2_record_review_r1_grok.txt
27:  + 188 representation-specific unpaired + 8 declared-role missing one side
```

## Finding 1 — the comparator's one-sided declared-role entries; what the record says about them
Where the reason is assigned in the comparator:
```
$ sed -n '1124,1126p' research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
                       'unpaired_reason':None if left and right else
                           ('no cross-engine role binding; representation-specific action' if '::' in role else
                            'no occurrence of this declared role in the other operand'),
```

Rows (stream lines) that carry the literal, and the count per line:
```
$ grep -n -F 'no occurrence of this declared role in the other operand' research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | cut -d: -f1
195
199
203
207
```

```
$ for n in 195 199 203 207; do printf 'line %s row %s literal count %s; empty-WL-side-then-literal count %s\n' $n "$(sed -n "${n}p" research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o '"row":"[^"]*"' | head -1)" "$(sed -n "${n}p" research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o -F 'no occurrence of this declared role in the other operand' | wc -l)" "$(sed -n "${n}p" research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o -F '"wl":[]},"unpaired_reason":"no occurrence of this declared role in the other operand' | wc -l)"; done
line 195 row "row":"energy_storage" literal count 1; empty-WL-side-then-literal count 1
line 199 row "row":"energy_transport" literal count 1; empty-WL-side-then-literal count 1
line 203 row "row":"energy_power" literal count 5; empty-WL-side-then-literal count 5
line 207 row "row":"energy_balance" literal count 1; empty-WL-side-then-literal count 1
```

Role key preceding each literal, in stream order:
```
$ for n in 195 199 203 207; do echo "line $n:"; sed -n "${n}p" research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o '"role":"[^"]*"\|"unpaired_reason":"no occurrence of this declared role in the other operand"' | grep -B1 'no occurrence' | grep role; done
line 195:
"role":"OPEN_MaterialCompatibility"
"unpaired_reason":"no occurrence of this declared role in the other operand"
line 199:
"role":"OPEN_MaterialCompatibility"
"unpaired_reason":"no occurrence of this declared role in the other operand"
line 203:
"role":"OPEN_MaterialCompatibility"
"unpaired_reason":"no occurrence of this declared role in the other operand"
"role":"OPEN_MaterialEnergyDensity"
"unpaired_reason":"no occurrence of this declared role in the other operand"
"role":"OPEN_MaterialEnergyFlux_0"
"unpaired_reason":"no occurrence of this declared role in the other operand"
"role":"OPEN_MaterialEnergyFlux_1"
"unpaired_reason":"no occurrence of this declared role in the other operand"
"role":"OPEN_MaterialEnergyFlux_2"
"unpaired_reason":"no occurrence of this declared role in the other operand"
line 207:
"role":"OPEN_MaterialCompatibility"
"unpaired_reason":"no occurrence of this declared role in the other operand"
```

```
$ for n in 203 207; do printf 'line %s OPEN_JointPowerAccounting occurrences: ' $n; sed -n "${n}p" research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out | grep -o 'OPEN_JointPowerAccounting' | wc -l; done
line 203 OPEN_JointPowerAccounting occurrences: 12
line 207 OPEN_JointPowerAccounting occurrences: 12
```

The record (r1) on this class and these roles:
```
$ for t in 'no occurrence of this declared role in the other operand' 'OPEN_MaterialCompatibility' 'MaterialEnergyDensity' 'MaterialEnergyFlux' 'JointPowerAccounting' 'missing-role'; do printf '%s: ' "$t"; grep -c -F "$t" research/pde_ledger_v3/steps/O2_steady_brane_balance.md; done
no occurrence of this declared role in the other operand: 0
OPEN_MaterialCompatibility: 0
MaterialEnergyDensity: 0
MaterialEnergyFlux: 1
JointPowerAccounting: 0
missing-role: 1
```

```
$ sed -n '219p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
| Every nonempty matched-role head/role/orientation or balance-entry inventory difference | **OPEN cross-engine difference; no owner named for adjudication.** M6 preserves every key and sign, including role-only/missing-role occurrences and their reasons. Repeated same-role unions are not per-occurrence comparisons. |
```

```
$ grep -n -F 'matching printed role/orientation/OPEN-free entries' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
377:| `ℬ_hold^live` in-plane/bulk/graph-normal components and `ℬ_E^steady` | Conditional named accounting with matching printed role/orientation/OPEN-free entries and six compared closed-part residuals. All OPEN inventory differences persist; no complete momentum, normal, work or energy balance is established equal or solved
```

```
$ grep -n -F 'avoid duplicate stress/exchange/' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
398:with compatible native maps, retain untruncated unknown terms/grades, avoid duplicate stress/exchange/
```

```
$ sed -n '397,399p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
Part D must keep coordinate versus graph-normal content distinct, use the same material and measure
with compatible native maps, retain untruncated unknown terms/grades, avoid duplicate stress/exchange/
power content, and name the supplier/budget if a closure requires net power. Closed-part zeros and
```

The r0 baseline (243a530f): the missing-role row existed; the entry-match handoff wording did not:
```
$ git show 243a530f:research/pde_ledger_v3/steps/O2_steady_brane_balance.md | grep -n -F 'missing-role'
176:| Every nonempty matched-role head/role/orientation or balance-entry inventory difference not settled as representation-specific syntax above | **OPEN cross-engine difference; no owner named for adjudication.** M6 preserves every key and sign, including role-only/missing-role occurrences and their reasons. Repeated same-role
```

```
$ git show 243a530f:research/pde_ledger_v3/steps/O2_steady_brane_balance.md | grep -c -F 'matching printed role'
0
```

Claude leg stdout (verbatim):
```
$ cat _scratch/s9b_build/o2_record_review_r1_claude_evidence/action_entry_census.stdout
234 ('PAIRED', 'nonempty')
188 ('no cross-engine role binding; representation-specific action', 'nonempty')
8 ('no occurrence of this declared role in the other operand', 'nonempty')
total entries 430
declared-role-absent entries (row, role, py_occ, wl_occ):
   ('energy_balance', 'OPEN_MaterialCompatibility', 9, 0)
   ('energy_power', 'OPEN_MaterialCompatibility', 5, 0)
   ('energy_power', 'OPEN_MaterialEnergyDensity', 1, 0)
   ('energy_power', 'OPEN_MaterialEnergyFlux_0', 1, 0)
   ('energy_power', 'OPEN_MaterialEnergyFlux_1', 1, 0)
   ('energy_power', 'OPEN_MaterialEnergyFlux_2', 1, 0)
   ('energy_storage', 'OPEN_MaterialCompatibility', 1, 0)
   ('energy_transport', 'OPEN_MaterialCompatibility', 3, 0)
```

```
$ cat _scratch/s9b_build/o2_record_review_r1_claude_evidence/py_energy_power_shape.stdout
== PY ENERGY_POWER top head: OPEN_JointPowerAccounting
    OPEN_JointPowerAccounting > OPEN_MaterialEnergyDensity
    OPEN_JointPowerAccounting > OPEN_MaterialEnergyFlux_0
    OPEN_JointPowerAccounting > OPEN_MaterialEnergyFlux_1
    OPEN_JointPowerAccounting > OPEN_MaterialEnergyFlux_2
    OPEN_JointPowerAccounting > OPEN_MaterialCompatibility
== WL B_E_STEADY/Entries/2/Object top head: OpenAction
```

```
$ grep -n 'energy_balance\|MaterialEnergyDensity\|MaterialEnergyFlux_0' _scratch/s9b_build/o2_record_review_r1_claude_evidence/energy_role_occurrences.stdout | head -12
17:   role=OPEN_MaterialEnergyDensity                    py_occ=1 wl_occ=1 unpaired=None
22:   role=OPEN_MaterialEnergyFlux_0                     py_occ=1 wl_occ=1 unpaired=None
43:   role=OPEN_MaterialEnergyDensity                    py_occ=1 wl_occ=0 unpaired=no occurrence of this declared role in the other operand
44:   role=OPEN_MaterialEnergyFlux_0                     py_occ=1 wl_occ=0 unpaired=no occurrence of this declared role in the other operand
64:STREAM_LINE 207 row=energy_balance
80:   role=OPEN_MaterialEnergyDensity                    py_occ=2 wl_occ=1 unpaired=None
81:   role=OPEN_MaterialEnergyFlux_0                     py_occ=2 wl_occ=1 unpaired=None
```

## Finding 2 — the register's rest-on test against its R-S8-06 widening; the record's parallel S8 deferrals
```
$ sed -n '414,416p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
- **requirement** — `u` as the material displacement of the stuff whose density is `ρ_br`, and the
  quadratic kinetic form of that material and the thickness degree of freedom (B's `μ_W`) on the
  original homogeneous/slab domain; for O2, the same material's OPEN flowing/embedded inertial response.
```

```
$ sed -n '425,430p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
- **O2 scope** — The same material displaced by `u` has live velocity `V`; O2's flowing/embedded
  inertial response `ℐ_br^live` and momentum map remain OPEN. The homogeneous quadratic anchor above
  is not a live kinetic law and supplies no identification of momentum density with `ρ_br V`.
  The material identity and live inertial response are merged here by object, with the original slab
  and thickness domain retained; no quadratic live form is required. Without them the formal O2
  storage/transport actions cannot become a physically supplied flowing response
```

```
$ sed -n '675,679p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md

O1/O3–O7, A13, live stress/energy/reference forms, missing grades/scales, core selection, source/response
`GM` matching and S21's integration/sort remain the record's OPEN handoff. They are not entered merely
because a later closed profile calculation would use them. The banked claim is not a solution that
already rests on those deliveries. In particular, the emitted mass residual is unsolved, no induced-
```

```
$ sed -n '695,697p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
this register says “a **postulate** is a value or form assumed *here* and not derived” and “A postulate
with a named retirement condition generates a requirement.” O2's premise provenance is retained as
conditional adopted input; open closure needs alone source no new entry under the rest-on test.
```

```
$ sed -n '656,658p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
interface reservoir obligation, with no live owner named. The distinct O2 net-supplier/budget
requirement is `R-O2-01`, with no owner named for its delivery; spec §6 expressly requires that object
to accompany any relation requiring net supply. It is entered on that rest-on condition, not simply
```

```
$ grep -n '^| O1 \|^| O7 \|live inertia$' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
320:| O1 `ℳ_⊥` | Constitutive input to closure; optical speed/inertia ratio supplies neither stiffness response nor full stress. | S8, substrate reduction, S22 completion. |
325:| O7 `ℛ_br` and missing grades/scales | Density/ordering input to closure and truncation. Live `ρ_br` is retained without eliminating it using the bulk EOS or assigning grades that remove terms. | S8, substrate reduction, S22 completion for density; missing order/derivative contracts stay unassigned where sources name non
328:energy-reference/improvement convention including `C_ref` (S1.5/S8 as applicable); live inertia
```

```
$ grep -n -F '(S8 with S1.5 antecedents)' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
329:(S8 with S1.5 antecedents), full stress and normal material response (S5–S8 ingredients, S8/S22
```

The widening was present in the r0 baseline and absent before the O2 pass:
```
$ git show 243a530f:research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | grep -n -F 'for O2, the same material'
414:  original homogeneous/slab domain; for O2, the same material's OPEN flowing/embedded inertial response.
```

```
$ git show c9db665d:research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | grep -c -F 'for O2, the same material'
0
```

