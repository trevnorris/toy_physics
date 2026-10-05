# D3 quadratic invariant completeness contract

Authorized 2026-09-16 after E1–E4 checkpoint `9e66a534`. Governed by
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). **J1–J4 is complete:**
the theorem, coverage, native identification and controls passed verification,
and both independent fidelity reviews returned CLEAR. See the
[review and closure record](D3_FIDELITY_REVIEW.md).

## Claim and domain

For all real quadratic forms on **all** real 3×3 matrices, with the full group
action `G ↦ R G Rᵀ`, prove the following finite contract. As in native Q9,
`G_ij = ∂_i u_j`, with row-major entries. This is a classification of pointwise
quadratic densities before any integration-by-parts or Euler–Lagrange quotient.
No positivity, field equations, symmetry of G, or transverse restriction is
assumed. Dimension three is supplied, not derived.

| Item | Required result |
|---|---|
| J1 — object and completeness | Every SO(3)-invariant quadratic form has a unique expression `a (tr G)² + b tr(G²) + c tr(G Gᵀ)`, and every such expression is invariant under all O(3). |
| J2 — coverage and counts | The actual SO(3) and O(3) invariant submodules coincide and have dimension 3. The reflection-minus submodule inside the SO(3) invariants is zero and has dimension 0. This includes the zero form and all real coefficient choices. |
| J3 — compact native identification | Compare the entire spans of actual native Q9 V1, V2 and V6 against the proved forms, including the empty odd basis and zero P_D; document coordinates and normalization. Matching dimensions alone is insufficient. |
| J4 — sensitivity and fidelity | Standard-axiom audit; mathematical mutations of completeness/count/oddness and form normalization with positive controls; two independent non-author reviews of a fixed packet and disposition of findings. |

The basis normalization is fixed by the displayed formula. In particular,
`tr(G²)` and `(tr G)²` are distinct forms; no equivalence or normalization is
inferred from sector names.

## Proof method and fidelity boundary

Use the associated bilinear form to justify an exhaustive finite coordinate
presentation of arbitrary quadratic forms. A small explicit collection of
proper rotations may impose necessary coefficient constraints, but sufficiency
must be proved for the **full** group, using trace identities. Coordinate
certificates are mathematical proof machinery, not a translation of CAS audit
outputs. Any generated proof is checked by Lean; its generator is not trusted
as a theorem oracle.

Reuse existing census patterns and pinned mathlib. Preserve completed D2 and
homogeneous proofs and their historical evidence. A compact native check may
call the existing Q9 helper for D3 without invoking production drivers or
regenerating exports. State precisely whether Wolfram evidence is source
inspection or execution. Instrument errors and timeouts are not mutations.

## Verified evidence and completed review

| Obligation | Evidence |
|---|---|
| J1 | `quadratic_representation`, `invariant_polynomial`, `SO_classification`, `SO_unique`, `invariantForm_O`. The finite certificate proves necessary constraints; trace identities prove full-group sufficiency. |
| J2 | `soEquiv`, `so_eq_o`, `odd_classification`, `odd_eq_bot`, `so_dimension`, `o_dimension`, `odd_dimension`, `census`. Counts are dimensions of the actual submodules. |
| J3 | `_measurements/S11_lean_d3_source_checks.json`: exact native V1/V2/V6 spans, zero odd P_D, two deliberate wrong-span controls. Wolfram source inspection only. |
| J4 local | Six modules and audit root, 40 standard-axiom audits, ten mathematical rejections, twelve positives; current source/object and package-pin validation. See [D3_VERIFICATION.txt](D3_VERIFICATION.txt). |
| J4 review | Claude and Grok independently returned CLEAR on the approved fixed 27-file packet. No blocking findings; optional points and review limits are recorded in [D3_FIDELITY_REVIEW.md](D3_FIDELITY_REVIEW.md). |

## Completion and exclusions

Stop after J1–J4: the full-group theorem, exact census, compact native span
connection, meaningful controls, source/object provenance and two cleared
fidelity reviews. These deliverables are complete. The user explicitly approved
this fixed packet with “yes go ahead”; both reviews are terminal and substantive.
The reviewed snapshot and archive remain unchanged. Only closure documentation
changed afterward; no proof, native source or instrument repair was required.

No D4/D5 census, general EL/divergence classification, dynamics or complete
XFORM_EXTRA spectrum, interface physics, S11c files/calculations, pinned export
changes, production reruns or systematic CAS bridge. No new commit requested.
Keep one Lean worker and sequential, memory-conscious jobs with durable logs
and silent completion/error hooks.
