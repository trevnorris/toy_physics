# D4 quadratic invariant completeness contract

Authorized 2026-09-16 after the completed D3 bulk checkpoint `43794555`.
Governed by [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).
**Complete within D4.1–D4.4.** Local verification passed and both independent
fidelity reviews returned CLEAR. See [D4_FIDELITY_REVIEW.md](D4_FIDELITY_REVIEW.md)
for the fixed reviewed revision, reviewer limits and finding dispositions.

## Claim and domain

Classify all real quadratic forms on all real 4×4 matrices under the full
conjugation action `G ↦ R G Rᵀ`. As in native Q9, `G_ij = ∂_i u_j`, with
row-major entries. This is a pointwise density classification before any
Euler–Lagrange or total-divergence quotient. Dimension four is supplied.
Coefficients are arbitrary real constants; G need not be symmetric. There are
no positivity, field-equation, nonzero-coefficient or polarization assumptions.

Use the three trace forms `(tr G)²`, `tr(G²)`, `tr(G Gᵀ)` and the normalized
orientation form

```
P(G) = (G01-G10)(G23-G32) - (G02-G20)(G13-G31)
       + (G03-G30)(G12-G21).
```

Here indices run from 0 to 3 and the orientation is `(0,1,2,3)`. The compact
native check establishes `P_D = P` and the fully summed
`epsilon_ijkl Gij Gkl = 2P`, with `epsilon_0123 = +1`. All four basis forms
are homogeneous of degree two in gradient entries; units are inherited from
the supplied gradient and coefficients, not a separate dimensional-analysis
formalization.

| Item | Required result |
|---|---|
| D4.1 — object and completeness | Every full SO(4)-invariant quadratic form has a unique expression `a (tr G)² + b tr(G²) + c tr(G Gᵀ) + d P(G)`. Prove sufficiency for every proper rotation, including the determinant transformation of P; finite tests alone are not sufficiency. |
| D4.2 — reflection and census | The full O(4) invariant submodule is exactly `d=0`, the reflection-odd submodule inside SO(4) is exactly `a=b=c=0`, and the actual submodule dimensions are 4/3/1. Include zero forms and all coefficient choices. |
| D4.3 — compact native identification | Compare the complete native Q9 V1/V2/V6 spans with these forms, including the exact P_D normalization and reflection action. Dimensions alone do not identify the object. A bounded helper call and source inspection suffice; no production run or export replacement. |
| D4.4 — sensitivity and fidelity | Standard-axiom audit, mathematical mutations of completeness, reflection, counts and normalization, admissible positive controls, and two independent non-author reviews of a fixed packet. |

## Method and evidence boundary

Reuse the established finite-coordinate method for necessity: justify the
complete 136-coefficient representation of quadratic forms on 16 entries;
select rational constraints from a small explicit collection of proper
rotations; have Lean check each constraint and its coefficient reconstruction.
The generator chooses proof certificates and is not a trusted theorem oracle.
Trace identities and a determinant identity prove full-group sufficiency.
The necessary equations may be split into sequential compilation blocks;
splitting does not remove any equation or change the reconstruction.
Preserve all completed H/I/E/J/K sources and historical evidence.

The compact native check identifies entire polynomial spans and conventions,
with a same-count/wrong-span negative control. It must distinguish actual
SymPy helper execution from Wolfram source inspection. No correspondence of
unexamined outputs is claimed. Instrument errors and timeouts cannot count as
rejected mathematical mutations.

## Evidence map

All theorem names below belong to `S11D4Invariants`. The thirteen modules and
audit root passed verification run4; run5 reused those builds after checking
their full local source closures, pins, generator and input/output object
hashes. All fourteen mathematical mutations and sixteen positives ran afresh
and passed their intended checks. The 49-declaration axiom audit uses only
standard axioms. See [D4_VERIFICATION.txt](D4_VERIFICATION.txt) and
`_measurements/S11_lean_d4_validation.json` for validated provenance.

| Obligation | Formal or native evidence | Sensitivity controls |
|---|---|---|
| D4.1: all quadratic forms, necessity and sufficiency | `quadratic_representation`; `ConstraintBlock0–3` and `invariant_polynomial`; `orientation_conjugate`, `invariantForm_SO`, `SO_classification`, `SO_unique` | Omission of each of the four forms; a single-entry square falsely claimed invariant; trace cross-coefficient and sign mutations. |
| D4.2: complete reflection split and dimensions | `O_classification`, `odd_classification`, `so_dimension`, `o_dimension`, `odd_dimension`, `even_odd_disjoint`, `even_odd_span`, `census` | Wrong dimensions 3/4/0; a nonzero odd form falsely claimed O-invariant; an even contribution falsely claimed reflection-odd. |
| D4.3: actual native object | `S11_lean_d4_source_checks.json`: complete V1/V2/V6 row spaces, reflection operator, `P_D=P`, epsilon contraction `2P` | Same-count/wrong-span and historical generator-orientation controls, evaluated only in memory. |
| D4.4: normalization, nonvacuity and review | `orientationForm_apply`, `orientation_conjugate`, audit root and `S11_lean_d4_contract_checks.json`; fixed packet manifest `S11_lean_d4_review_packet.json` | Wrong P value and reflection sign; paired true statements; nonzero odd, zero, negative-coefficient and unique-coefficient positives. Claude and Grok independently returned CLEAR. |

The reflection-even and reflection-odd spaces overlap at zero; the proved
disjointness means their submodule intersection is exactly zero. Their sum
equals the whole SO-invariant space, so the decomposition is unique.

Claude and Grok reviewed the fixed 34-file packet after the user's explicit
“Proceed” authorization. Both returned substantive terminal CLEAR reports
with no blockers. Their optional observations are dispositioned in
[D4_FIDELITY_REVIEW.md](D4_FIDELITY_REVIEW.md); no proof changes were required.

## Completion and exclusions

D4.1–D4.4 are proved, the compact native identity is established, controls and
source/object/pin checks pass, and both independent fidelity reviews clear
with findings resolved. The finite completion criteria are satisfied; stop
adding coverage to this classification. The reviewed snapshot and historical
local evidence retain their original bytes and pending-review labels.

The D4 orientation term's divergence/zero bulk variation is the next separate
increment, not part of this density classification. This contract excludes
D5, variable coefficients, general null-Lagrangian classification, new
spectrum/stability/interface work, all S11c files/calculations and pinned
exports, production reruns and systematic CAS bridging. Use sequential,
memory-conscious jobs and silent local completion/error hooks.
