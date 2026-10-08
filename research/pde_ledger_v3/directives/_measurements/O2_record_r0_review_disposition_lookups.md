# Measurements — O2 record r0 review dispositions (generated 2026-10-07 22:40)

Generator: `_scratch/s9b_build/gen/o2_record_r0_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). Leg evidence is quoted verbatim from the Claude leg's filed stdout, copied from `/tmp/o2_record_review/` to `_scratch/s9b_build/o2_record_review_r0_claude_evidence/`.

```
$ sha256sum -c _scratch/s9b_build/o2_record_review_baseline_r0.sha256
research/pde_ledger_v3/steps/O2_steady_brane_balance.md: OK
research/pde_ledger_v3/steps/_measurements/O2_record_measurements.md: OK
research/pde_ledger_v3/scripts/O2_record_measurements.py: OK
research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md: OK
```

## F1 — R-S1-03's target qualifier (committed c9db665d vs working tree); the record on time-reversibility
```
$ git show c9db665d:research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | grep -n 'register inference; no owner named by the records'
373:  (**register inference; no owner named by the records**) · **status** OPEN
390:  (**register inference; no owner named by the records**) · **status** OPEN
```

```
$ grep -n 'O2 explicitly names S8\*\*' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
394:  (**register inference for the original S11/S11b records; O2 explicitly names S8**) · **status** OPEN
411:  (**register inference for the original S11/S11b records; O2 explicitly names S8**) · **status** OPEN
```

```
$ grep -n 'R-S1-03' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | head -4
69:| **S1** *(register inference)* | `R-S1-03` | whether the substructure's microdynamics is time-reversible | S11b-B / unified S11b |
390:### R-S1-03 — the substructure's microscopic time-reversibility
520:S1 `R-S1-03`; S8 `R-S8-06`; S12 `R-S12-01`/`02`.
573:sources none here. Onsager–Casimir is a candidate only; `R-S1-03` instead rests on B's
```

```
$ grep -c -i 'revers\|onsager' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
0
```

## F2 — the record's classification of unbound roles; where the comparator's label comes from
```
$ sed -n '173p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
| Unpaired role/head syntax explicitly labelled “no cross-engine role binding; representation-specific action” | **Representational role-binding difference only to that stated extent** (M6). Includes PY domain/history/chart/reduction/rotational-work actions and WL native measure/tangent/unit-normal/hold actions. The accompan
```

```
$ sed -n '1108p;1125,1126p' research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
            role = maps[engine].get(key,engine+'::'+key)
                           ('no cross-engine role binding; representation-specific action' if '::' in role else
                            'no occurrence of this declared role in the other operand'),
```

## F3 — the comparator's balance structure, OPEN difference counts and residual totals (Claude leg stdout)
```
$ grep -h 'identical_role_orientation_multisets\|nonempty_fields_by_role_count' _scratch/s9b_build/o2_record_review_r0_claude_evidence/balance_role_orientation.stdout | head -10
line 643 balance_entries row=hold_inplane component=[0] n_py=8 n_wl=8 identical_role_orientation_multisets=True
line 644 balance_comparison row=hold_inplane component=[0] roles=8 nonempty_fields_by_role_count={'named_OPEN_operands': 7, 'live_arguments': 3}
line 645 balance_entries row=hold_inplane component=[1] n_py=8 n_wl=8 identical_role_orientation_multisets=True
line 646 balance_comparison row=hold_inplane component=[1] roles=8 nonempty_fields_by_role_count={'named_OPEN_operands': 7, 'live_arguments': 3}
line 647 balance_entries row=hold_inplane component=[2] n_py=8 n_wl=8 identical_role_orientation_multisets=True
line 648 balance_comparison row=hold_inplane component=[2] roles=8 nonempty_fields_by_role_count={'named_OPEN_operands': 7, 'live_arguments': 3}
line 650 balance_entries row=hold_bulk component=[] n_py=8 n_wl=8 identical_role_orientation_multisets=True
line 651 balance_comparison row=hold_bulk component=[] roles=8 nonempty_fields_by_role_count={'named_OPEN_operands': 8, 'live_arguments': 3}
line 653 balance_entries row=hold_normal component=[] n_py=32 n_wl=32 identical_role_orientation_multisets=True
line 654 balance_comparison row=hold_normal component=[] roles=27 nonempty_fields_by_role_count={'named_OPEN_operands': 26, 'live_arguments': 12}
```

```
$ grep -h 'paired role entries\|balance role groups' _scratch/s9b_build/o2_record_review_r0_claude_evidence/open_difference_counts.stdout
paired role entries with nonempty differences: 234  with all-empty differences: 0
unpaired role entries: 196
total balance role groups with differences: 60
```

```
$ grep -h 'leaf outcome totals\|NONZERO' _scratch/s9b_build/o2_record_review_r0_claude_evidence/enumerate_comparator.stdout
leaf outcome totals: {'not_formed': 212, 'exact_zero': 136}
NONZERO exact residuals: 0
```

```
$ sed -n '11,12p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
preserving differences in OPEN content. An empty OPEN difference means only the printed inventories
match. Sources: `O2_SHARED_PHYSICS.md`, §§1–10; engine acceptance `O2_build_r3_review_disposition.md`;
```

## F4 — the record's WL-only derivative sentence; the balance components the keys appear in (Claude leg stdout)
```
$ sed -n '179,181p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
In particular, M6's in-plane `MomentumCurrent`/`OPEN_MomentumFlux_i_j` balance entries print
WL-only first derivatives of `V_r`, `δ`, `f`, `h`, `j_n`, `μ_⊥`, `ρ_br`, and `xi_w''(r)`;
the `MaterialEnergyCurrent` entries likewise print differing live-derivative content. The comparator
```

```
$ sed -n '156,173p' _scratch/s9b_build/o2_record_review_r0_claude_evidence/enumerate_distinct_differences.stdout
=== 3. BALANCE entry differences, distinct (field, key, py_minus_wl) ===
distinct: 45
 live_arguments -1 ProfileDerivative<V_r>(Sequence(1),Sequence(Pow(Add(Pow(x1,2),Pow(x2,2),Pow(x3,2)),1/2)))
      balances=[('energy_balance', ()), ('hold_bulk', ()), ('hold_inplane', (0,)), ('hold_inplane', (1,)), ('hold_inplane', (2,)), ('hold_normal', ())]; roles=15
 live_arguments -1 ProfileDerivative<delta>(Sequence(1),Sequence(Pow(Add(Pow(x1,2),Pow(x2,2),Pow(x3,2)),1/2)))
      balances=[('energy_balance', ()), ('hold_bulk', ()), ('hold_inplane', (0,)), ('hold_inplane', (1,)), ('hold_inplane', (2,)), ('hold_normal', ())]; roles=15
 live_arguments -1 ProfileDerivative<f>(Sequence(1),Sequence(Pow(Add(Pow(x1,2),Pow(x2,2),Pow(x3,2)),1/2)))
      balances=[('energy_balance', ()), ('hold_bulk', ()), ('hold_inplane', (0,)), ('hold_inplane', (1,)), ('hold_inplane', (2,)), ('hold_normal', ())]; roles=15
 live_arguments -1 ProfileDerivative<h>(Sequence(1),Sequence(Pow(Add(Pow(x1,2),Pow(x2,2),Pow(x3,2)),1/2)))
      balances=[('energy_balance', ()), ('hold_bulk', ()), ('hold_inplane', (0,)), ('hold_inplane', (1,)), ('hold_inplane', (2,)), ('hold_normal', ())]; roles=15
 live_arguments -1 ProfileDerivative<j_n>(Sequence(1),Sequence(Pow(Add(Pow(x1,2),Pow(x2,2),Pow(x3,2)),1/2)))
      balances=[('energy_balance', ()), ('hold_bulk', ()), ('hold_inplane', (0,)), ('hold_inplane', (1,)), ('hold_inplane', (2,)), ('hold_normal', ())]; roles=15
 live_arguments -1 ProfileDerivative<mu_perp>(Sequence(1),Sequence(Pow(Add(Pow(x1,2),Pow(x2,2),Pow(x3,2)),1/2)))
      balances=[('energy_balance', ()), ('hold_bulk', ()), ('hold_inplane', (0,)), ('hold_inplane', (1,)), ('hold_inplane', (2,)), ('hold_normal', ())]; roles=15
 live_arguments -1 ProfileDerivative<o2_rho_br_live>(Sequence(1),Sequence(Pow(Add(Pow(x1,2),Pow(x2,2),Pow(x3,2)),1/2)))
      balances=[('energy_balance', ()), ('hold_bulk', ()), ('hold_inplane', (0,)), ('hold_inplane', (1,)), ('hold_inplane', (2,)), ('hold_normal', ())]; roles=15
 live_arguments -1 ProfileDerivative<xi_w>(Sequence(2),Sequence(Pow(Add(Pow(x1,2),Pow(x2,2),Pow(x3,2)),1/2)))
      balances=[('energy_balance', ()), ('hold_bulk', ()), ('hold_inplane', (0,)), ('hold_inplane', (1,)), ('hold_inplane', (2,)), ('hold_normal', ())]; roles=15
```

## F5 — the Part D handoff; spec admissible dependences and the mass-law qualification
```
$ sed -n '302,309p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md

A later S9b Part D calculation may use the supplied live radial setting, coordinate mass law,
graph/field identity, optical ratio/anchoring/counting and user-adopted premises 1–4, and may refer
to the accepted emitted geometry, material velocity, in-plane carried momentum, named `ℬ_hold^live`
component accounting and `ℬ_E^steady` with their full OPEN inputs and differences. It must keep
coordinate versus graph-normal content distinct, use the same material/measure and compatible native
maps, retain untruncated unknown terms/grades, avoid duplicate stress/exchange/power content, and
name the supplier/budget if a closure requires net power. It may not take closed-part comparator
```

```
$ sed -n '124,125p' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md
instantaneous response, tensor realization, constitutive family or derivative cutoff. Live fields,
gradients and material history remain admissible dependences (C §1).
```

```
$ sed -n '330,336p' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md
`GM` is the independent slow-test-matter orbital parameter, not a source, profile, mouth or drain
amplitude. Only `μ_⊥/ρ_br` inherits the speed-change grade. Keep the full density factor and
derivatives in the mass input. The recorded induced-metric qualification is a relative `O(ε)`
correction to `j_n`, not a live O6 law. O2 uses the law on the declared coordinate measure (§1) and
supplies no induced-metric mass balance. The qualification is carried as a limit on claims, not as a
term in the object: a claim that reads this `j_n` or `ρ_br` as a density per induced measure, or
compares it with an induced-metric mass balance, carries it. Fixed `ℓ` transfers the live `h` slope
```

## F6 — R-S12-01 in the working tree; the sources' owners for the O2 net supplier/budget
```
$ grep -n 'for O2, a named physical supplier\|merged here rather than creating a second' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
438:  coupling; for O2, a named physical supplier and its budget if a closure requires net supply.
```

```
$ sed -n '282p' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md
| `𝒮_E,net`, `𝒫_E,supply` | Identity of the physical supplier of any net power and its stated power budget. Both stay general OPEN inputs. Naming the drain does not specify available energy or a budget. | O2 requirement; source/holder/supplier forms retain their owners. |
```

```
$ sed -n '311,314p' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md

No non-passive interface response is selected. If a later closure adopts one, it must name its
reservoir and state its power budget. This conditional obligation is separate from S12's
non-variational source partners and from premise 1's relaxation power. The historical S11b-C → S11c
```

```
$ sed -n '457,461p' research/pde_ledger_v3/directives/O2_input_contract.md
This contract chooses neither `C_ref=0` nor an extension of this `n=5` expression to v9's symbolic `n`.

The existing non-passive-interface condition remains separate: if a later model adopts such a response,
it must name a reservoir and state a power budget (`steps/S11b_interface_coupling_law.md:57–63`;
`steps/S11bB_interface_assembly.md:195–197`). This does not select that response or settle relaxation
```

