# Measurements — S9b repair decision-list review dispositions (generated 2026-10-08 10:23)

Generator: `_scratch/s9b_build/gen/s9b_repair_dl_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). The reviewed version is `_scratch/s9b_build/S9b_repair_decision_list_reviewed_v0.md`.

```
$ cut -c1-64 _scratch/s9b_build/s9b_repair_dl_review_baseline.sha256; sha256sum _scratch/s9b_build/S9b_repair_decision_list_reviewed_v0.md | cut -c1-64
f4aaea5f0e4cb7548a15740ea3a20a3af750d60b5bd5ec330aadfc58afef64dc
f4aaea5f0e4cb7548a15740ea3a20a3af750d60b5bd5ec330aadfc58afef64dc
```

```
$ grep -n -o 'Verdict:[* ]*[A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_dl_review_codex_final.txt _scratch/s9b_build/s9b_repair_dl_review_grok.txt
_scratch/s9b_build/s9b_repair_dl_review_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_repair_dl_review_grok.txt:1:Verdict:** NEEDS REVISION
```

## C1 / G2 — the returned obligation, the stress operand and O1
```
$ grep -n 'Owner.\*\* The interpretation' -A1 _scratch/s9b_build/S9b_repair_decision_list_reviewed_v0.md
71:- **Owner.** The interpretation stays OPEN with the stress-law owner (S8; O1). The S9b record enters it in the
72-  register.
```

```
$ grep -n 'Part D supplies the in-plane momentum flux' _scratch/s9b_build/S9b_repair_decision_list_reviewed_v0.md
67:- **In-plane momentum flux.** Part D supplies the in-plane momentum flux from P5 and P6. It therefore takes no
```

```
$ sed -n '527,529p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
(M1). This is an undischarged sub-step-7 interpretation obligation returned to the orchestrator
under M1, who must route the work needed beyond retrieval; it is not an unowned difference.
No specific adjudication owner is named for the remaining cross-engine differences.
```

```
$ sed -n '440,444p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
The first route merges shared objects rather than turning each OPEN closure input into a new
requirement. Its one rest-on criterion enters only sourced conditions on which the reported
conditional object rests: its stated material/phase identifications and explicitly carried
accounting/admissibility obligations. Unselected response forms needed only to close a later state
are OPEN handoffs, even when S8 is their named future owner. The register pass applies this criterion
```

```
$ grep -n 'all eight WL-only derivative keys at every location' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
470:including all eight WL-only derivative keys at every location in §4 and all PY-only native/chart,
```

```
$ grep -n '^| `𝒯_br^live`' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | cut -c1-200
130:| `𝒯_br^live`; **OPEN form**, with character fixed by **adopted premise 1**, C §3 | Full live in-plane and normal stress, including any relaxation content. It enters the internal material-forc
```

```
$ grep -n '^| `ℳ_⊥` (O1)' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | cut -c1-260
137:| `ℳ_⊥` (O1); **OPEN**, C §6 | General stiffness response, including frequency regime and loading/material history, bulk/flow/embedding/thickness/projection dependence and gradients. It is a constitutive input to the linked optical/material descriptio
```

## G1 — the selected objects for P6 and the density link
```
$ grep -n 'asymptotic state' _scratch/s9b_build/S9b_repair_decision_list_reviewed_v0.md
27:| **P6** (2026-10-08) | In steady flow the brane's in-plane stress is an isotropic pressure `p_br(ρ_br)`, because shear relaxes under steady load (P1). Its compressional speed `c_comp`, with `c_comp² = dp_br/dρ_br` at the asymptotic state, is a live symbol. In the optical regime, light still sees the elastic transverse sti
28:| **Density link** (2026-10-06) | v8's `c_γ² = μ_⊥/ρ_br` is kept, with `μ_⊥` a function of `ρ_br`. Its logarithmic slope `α ≡ d ln μ_⊥ / d ln ρ_br`, at the asymptotic state, is a live symbol. | — |
```

```
$ grep -n 'μ_⊥(x)' <(git show c2f1cf2b:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md)
27:  c_γ(x)² ≡ μ_⊥(x)/ρ_br(x) .
31:  S9's basis writes it as `μ_R` (`steps/S9_light_requires_shear.md:78–79`). Both `μ_⊥(x)` and `ρ_br(x)`
```

## G3 — two order policies in D4
```
$ grep -n 'relative `O(ε)`' _scratch/s9b_build/S9b_repair_decision_list_reviewed_v0.md
53:- **v8's sentence.** v8 says the induced-metric mass balance differs from the flat form by "a relative `O(ε)`
58:  relative `O(ε)` for induced-measure readings (O2 record §8).
```

## C2 / G4 — the build commit's findings against D6
```
$ git show -s --format=%B bb94b885 | sed -n '5,16p'
- S9b_exports.py binds ansatz symbols w, q, L onto unrelated upstream rows
  (bulk normal coordinate, wave-norm coordinate, half-interval size) via a
  bare-symbol F9B_EQUAL; rows retagged as corroborated S9b KNOBs (both legs).
- Every-b requirement printed as an unevaluated ConditionSet/ForAll; the
  conditions and implied j_n are not reduced (Claude leg).
- SymPy Part C rows conjoin domain predicates built from the pre-response
  amplitude d (Grok leg; d appears 7x in the delta=0 row).
- Radar ln(1/b^2) coefficient from a hand rule, not from TIME_RADAR; path
  independence computed on a placeholder; typed demo refusal (Claude leg).
- Tag sets are not parallel across engines (directive item 4).
Part A agrees with independent derivations and exact rays in both engines.
Do not add S9b_exports.py to any fold until a repaired delta is accepted.
```

```
$ grep -n 'F9B_EQUAL' research/pde_ledger_v3/scripts/S9b_exports.py | cut -c1-200
23:'L': {'value': _restore("Symbol('L', positive=True)"), 'display': 'L', 'value_kind': 'COMPUTED_OBJECT', 'class': 'KNOB', 'step': 'S9b', 'f9_operands': _restore("Tuple(Symbol('L', positive=True), Sy
48:'q': {'value': _restore("Symbol('q', positive=True)"), 'display': 'q', 'value_kind': 'COMPUTED_OBJECT', 'class': 'KNOB', 'step': 'S9b', 'f9_operands': _restore("Tuple(Symbol('q', positive=True), Sy
54:'w': {'value': _restore("Symbol('w', real=True)"), 'display': 'w', 'value_kind': 'COMPUTED_OBJECT', 'class': 'KNOB', 'step': 'S9b', 'f9_operands': _restore("Tuple(Symbol('w', real=True), Symbol('w'
```

```
$ grep -n '^- \*\*B[0-9]' _scratch/s9b_build/S9b_repair_decision_list_reviewed_v0.md
77:- **B1.** The exports bind `w`, `q` and `L` onto unrelated upstream rows by bare-symbol equality.
78:- **B2.** The every-`b` requirement is restated, not reduced. No condition is solved relative to `GM` per stratum,
80:- **B3.** The radar `ln(1/b²)` coefficient comes from a hand rule, not from Part A's round trip.
81:- **B4.** Path independence is computed on a placeholder.
82:- **B5.** The demo refusal is a typed literal.
83:- **B6.** The Part C domain predicates use the pre-response amplitude `d`.
84:- **B7.** The tag sets are not parallel across engines.
```

## Process — the Grok leg read files outside the repository
```
$ grep -n -o '[^ ]*project_[a-z0-9_]*\.md' _scratch/s9b_build/s9b_repair_dl_review_grok.txt
29:.../project_premise6_relaxed_shear_pressure.md
32:.../project_brane_density_links_light_speed_and_flow.md
```

