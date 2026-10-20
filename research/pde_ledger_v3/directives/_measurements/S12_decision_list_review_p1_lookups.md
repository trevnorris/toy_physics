# Measurements: S12 decision list, review pass 1 (generated 2026-10-19 22:29)

Generator: `_scratch/s12/gen/dl_p1_lookups.sh`. Literal retrieval only. The list is read from the working tree;
H1 shows its hash equals the reviewed baseline. The reviewed list is preserved in the commit that adds this file.

## H1. The pass: prompt, baseline, exits, verdict lines
```
$ sha256sum /var/projects/toy_physics/_scratch/s12/dl_review_prompt_p1.md | cut -c1-64
509c33d8e57b498a84a098a25f5c55c40e02ce193ff454febe421f28ff9be1a4
```

```
$ sha256sum -c /var/projects/toy_physics/_scratch/s12/dl_baseline_p1.sha256 | sed 's#/var/projects/toy_physics/##'
research/pde_ledger_v3/directives/S12_decision_list.md: OK
research/pde_ledger_v3/directives/_measurements/S12_setup_lookups.md: OK
```

```
$ cat /var/projects/_review_scratch/s12-dl-p1-codex/DONE /var/projects/_review_scratch/s12-dl-p1-grok/DONE
EXIT=0
EXIT=0
```

```
$ grep -n '^FINDINGS\|^\*\*FINDINGS\*\*\|^CLEAR\|^[0-9]\. [A-Z]*\b —\|^[0-9]\. \*\*' /var/projects/_review_scratch/s12-dl-p1-codex/REVIEW.txt /var/projects/_review_scratch/s12-dl-p1-grok/grok_review.txt | cut -c1-160
/var/projects/_review_scratch/s12-dl-p1-codex/REVIEW.txt:61:FINDINGS
/var/projects/_review_scratch/s12-dl-p1-codex/REVIEW.txt:63:1. HIGH — The claimed U1 encoding has no stated admissibility relation; its
/var/projects/_review_scratch/s12-dl-p1-codex/REVIEW.txt:123:2. MEDIUM — Vanishing deviation does not by itself give the asserted source limit.
/var/projects/_review_scratch/s12-dl-p1-codex/REVIEW.txt:172:3. MEDIUM — The conversion-partner inventory drops the OPEN reaction and energy
/var/projects/_review_scratch/s12-dl-p1-codex/REVIEW.txt:231:4. MEDIUM — The OPEN list collapses the native O6/S12 handoff into the S14a bridge.
/var/projects/_review_scratch/s12-dl-p1-grok/grok_review.txt:9:**FINDINGS**
/var/projects/_review_scratch/s12-dl-p1-grok/grok_review.txt:11:1. **Medium.** The conversion-partner inventory treats O2 premise 3 as the momentum statement an
/var/projects/_review_scratch/s12-dl-p1-grok/grok_review.txt:56:2. **Medium.** The boundary inventory replaces the throat’s IR/return conditions with one oute
/var/projects/_review_scratch/s12-dl-p1-grok/grok_review.txt:83:3. **Low.** The same inventory calls the S9b/O2 domain a resting brane. In those steps the mass 
```

## H2. The U1 encoding (Codex 1)
```
$ grep -n 'would\* have is set\|plus a live deviation\|This encoding is wrong if' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
49:**Identification, flagged.** U1 is encoded as follows: the order a point *would* have is set by where it sits
51:slab's thickness `W(x,t)`. The actual `χ_B` is `F` plus a live deviation from it. How fast `χ_B` approaches `F` is
53:general case. It is not posed. This encoding is wrong if order depends on something other than position relative to
```

```
$ grep -n 'Material holds shear only\|de-orders it\. Bulk' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
38:  be one key way of explaining why it loses the shear as it falls into the throat." Material holds shear only
40:  de-orders it. Bulk material that settles into the slab re-orders. The user confirmed this one rule.
```

## H3. The limit claim (Codex 2)
```
$ grep -n 'reached as a limit\|follows as the limit case\|at any timescale' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
52:**not** chosen. The case with no deviation, where order follows the slab instantly, is reached as a limit of the
61:   not enter. The same object with the deviation set to zero follows as the limit case.
78:- **Regime:** none chosen. The objects are local, with every field live, at any timescale.
```

## H4. Conversion partners and power obligations (Codex 3, Grok 1)
```
$ grep -n 'Converted material carries the local brane velocity' -A2 /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
74:- **Momentum and energy partners of conversion.** Converted material carries the local brane velocity `V` (O2
75-  premise 3, `directives/O2_premise_decision_list.md`). The energy of ordering, and who supplies it, stay OPEN (`R-O2-01`, `R-S12-01`).
76-
```

```
$ sed -n '32,35p' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_premise_decision_list.md
3. **Exchange momentum: the brane's local velocity.** Converted material carries the brane material's local
   velocity `V`. This closes the transported momentum only.
   - Any additional non-variational momentum partners, and the system that carries their reaction, stay OPEN
     with S12.
```

```
$ grep -n 'Premise 3 fixes transported momentum, not the carried total energy' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md
280:| `𝒫_convert/exchange^live` | Order-conversion work, energy carried with exchanged material, and additional non-variational energy partners, with their reaction/supply systems. Premise 3 fixes transported momentum, not the carried total energy. | S12. |
```

```
$ grep -n 'It does not' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | grep -n 'supply the reduced bulk-direction'
1:243:It does not supply the reduced bulk-direction component `(Π_n^carry)^w` as `j_n` times one
```

```
$ grep -n 'none is assigned\|from S12.s additional non-variational source partners' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
521:  from S12's additional non-variational source partners, premise 1's relaxation power, and the
568:  object. S12 source forms and Q2/S22 holder forms keep their separate assignments; none is assigned
```

## H5. The native O6 handoff versus the S14a bridge (Codex 4)
```
$ grep -n 'which tie .Γ_B.' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
87:- the projection onto the slab and the far-field flux bridge, which tie `Γ_B` to the brane's exchange `j_n`
```

```
$ grep -n '^| O6' /var/projects/toy_physics/research/pde_ledger_v3/steps/O2_steady_brane_balance.md
399:| O6 `𝒥_map` | Exchange/projection and bulk-boundary input; common exchanged material and compatible mass/momentum/traction/work geometry and measures remain required but uncomputed. | S12 dynamical conversion and distinct native source/boundary data; S14a's projected order-loss/far-field bridge remains distinct. |
```

```
$ grep -n 'O6 identifications, and additional non-variational' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md
605:| **S12** | Conversion/source/controller inventory, distinct boundary/domain inventory, O6 identifications, and additional non-variational momentum/energy partners with reaction/supply systems. O2's local-velocity transported momentum does not determine those partners, a return law, branch choice or energy supply. |
```

```
$ grep -n 'the projected order-loss source · the definition' -A1 /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
636:⇒ **Derive, from the dynamical system:** the projected order-loss source · the definition of `J` · the
637-profile-dependent `J → Q` map · the controlled `ω → 0` return law.
```

## H6. The boundary inventory (Grok 2)
```
$ grep -n 'Galaxy-scale behaviour\|The outer edge stays open\|the outer edge and' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
16:follow in their own steps. S21 compares them, and a mismatch reopens S9b (user, 2026-10-19). Galaxy-scale behaviour
71:- **Boundary and domain data.** The mouth, the collar, the outer edge and the `𝔅` branch selection. The outer
```

```
$ sed -n '612p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
- **boundary / domain data** — the mouth, collar and IR conditions, and the `𝔅` branch selection.
```

```
$ sed -n '48p' /var/projects/toy_physics/research/pde_ledger_v2/notes/stage045_nonvariational_block_prep.md
- **The `𝔅` sleeve-contact + mouth-Robin/Neumann + collar-Kirchhoff + two-sided-IR-return BCs** (the boundary taxonomy —
```

```
$ grep -n '^| IR / return scheme' /var/projects/toy_physics/docs/two_throat_simulation_handoff_spec.md
340:| IR / return scheme | One- vs two-sided domain, return/absorber law, far-field truncation/extrapolation, treatment of zero modes/moduli, and (for moving/radiative rungs) the outgoing/absorbing prescription. |
```

## H7. The S9b/O2 domain (Grok 3)
```
$ grep -n 'resting' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
72:  edge stays open: a body sits inside a galaxy's flow inside an expanding brane. "One body alone in a resting
```

```
$ grep -n 'mass at rest with its drain flowing\|brane.s far-field rest frame' /var/projects/toy_physics/research/pde_ledger_v3/steps/O2_steady_brane_balance.md
74:The supplied setting is one isolated, spherically symmetric mass at rest with its drain flowing:
75:far field, Eulerian steady material profiles, lab time in the brane's far-field rest frame, linear
```

```
$ grep -n 'one isolated spherical mass at rest' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
84:The domain is one isolated spherical mass at rest, Eulerian steady draining profiles, far field,
```

