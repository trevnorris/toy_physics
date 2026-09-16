# D3 bulk variation and equivalence contract

Authorized after J1–J4 checkpoint `0ecf30a8`: the user accepted the proposed
D3 bulk-variation increment and asked to continue. Governed by
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). Status: **complete for K1–K4**; both independent reviews CLEAR.
See [D3_BULK_FIDELITY_REVIEW.md](D3_BULK_FIDELITY_REVIEW.md).

The gap is S11 Move 2 / Q9 V5: the three classified invariant densities need
not give three independent bulk equations. Reuse the completed D3 classification
and S10 smooth-field, compact-test and integration-by-parts machinery.

| Item | Finite deliverable |
|---|---|
| K1 — density and variation | For arbitrary constant real coefficients `(a,b,c)`, identify the supplied stiffness density `Q(G)=a(tr G)^2+b tr(G^2)+c tr(G G^T)` from J1–J4 and use `L=-Q/2`, with `G_ij=partial_i u_j`. Derive derivative-defined momenta, the finite relative-action first variation and its local Euler–Lagrange expression on smooth fields on R^(3+1), tested against smooth compact variations. |
| K2 — exhaustive bulk equivalence | Prove that this three-parameter family has zero first variation on every smooth background exactly when `c=0` and `a+b=0`. Exhibit a divergence current for `(tr G)^2-tr(G^2)`. Classify equality of bulk operators by equality of `(c,a+b)`; the null family is exactly the one-dimensional span of `(1,-1,0)` and the independent bulk response has two parameters. Include zero, negative coefficients and admissible non-null examples. |
| K3 — physical operator map | Identify the actual density-derived plane-wave operator for every real k, including k=0, and its longitudinal and transverse action. Connect it to the existing homogeneous operator with `mu=c`, `B=a+b+c` and zero inertia. This is an operator identity, not a new spectrum calculation. |
| K4 — compact fidelity and closure | Compare native D3 Q9 V5 for the actual computed V1 basis with this general operator, retaining native basis changes and the overall variational sign. Keep sign/factor, false-null, false-equivalence and wrong-stiffness-map mutations with positives. Verify builds and standard axioms with one worker, freeze a compact packet, obtain two independent non-author fidelity reviews and resolve findings. |

Conventions: coordinates `(t,x1,x2,x3)`; `J j i=partial_j u_i`; spatial row i
is `J i.succ`. Phase `k·x-omega*t`. Lean variational derivative is
`-sum_j partial_j(dL/dJ_ji)`. The native helper uses the opposite overall sign;
record that explicitly. Coefficients are constant real stiffness coefficients,
with no positivity or nonzero assumptions. Dimension three is supplied.
Smooth backgrounds need not have finite total action: only the density change
under a compact variation is integrated. The divergence identity does not
discard boundary physics when variations reach a physical boundary.

Stop after K1–K4. This classifies nullness only within the already proved D3
invariant family, not general null Lagrangians. No D4/D5, variable coefficients,
new roots or stability classification, interfaces, S11c calculations/files,
pinned exports, production reruns or systematic CAS bridge. Preserve completed
H/I/E/J sources and historical records and unrelated remote-compute notes.
Use sequential jobs, one Lean worker, and silent local completion/error hooks.

No commit requested for this increment. The user explicitly approved the fixed
41-file K1–K4 packet; Claude and Grok independently reviewed it and returned
CLEAR with no blockers. The approved packet remains frozen.

## Verified evidence and closure

All five modules and the audit root pass. The current record validates all 49
new theorem roots with standard axioms, twelve intended mathematical rejections
and eleven positive controls. The native D3 V5 connection also passes, including
basis changes and sign/factor sensitivity. See
[D3_BULK_VERIFICATION.txt](D3_BULK_VERIFICATION.txt) and the source/compiled-object
hashes in `../../_measurements/S11_lean_d3_bulk_validation.json`.

| Obligation | Principal checked evidence |
|---|---|
| K1 | `density_identity`, `all_invariant_densities`, `momentum_eq`, `relativeAction_hasDerivAt`, `relativeAction_deriv_eq_eulerLagrange`, `actionStationary_iff_eulerLagrange`, `eulerLagrange_eq` |
| K2 | `boundary_identity`, `null_density_is_divergence`, `bulkEquivalent_iff`, `variationallyNull_iff`, `null_parameterization`, `nullSpace_identification`, `bulk_response_dimension`, `null_dimension`, `exists_nonzero_firstVariation` |
| K3 | `momentum_contraction`, `eulerLagrange_planeWave`, `homogeneous_operator`, `modal_longitudinal`, `modal_transverse`, `modal_zero_wavevector` |
| K4 | `S11_lean_d3_bulk_contract_checks.json`, `S11_lean_d3_bulk_source_checks.json`, `S11_lean_d3_bulk_validation.json`, and both CLEAR reports with dispositions in `D3_BULK_FIDELITY_REVIEW.md` |

The equivalence classes are indexed exhaustively by all `(c,a+b) in R^2`;
`response_surjective` proves every class parameter is attainable. The zero
response class is precisely the null family; every other class has a smooth
background and compact variation with nonzero first variation. This is equality
of operators on all smooth fields, not a generic-wavevector sampling claim.

Both substantive independent reviews are complete with no unresolved blockers.
Optional observations do not require additional formalization. Current source,
dependency, instrument and object correspondence is validated in
`../../_measurements/S11_lean_d3_bulk_closure_validation.json`. All K1–K4
deliverables are complete; stop here unless a new bounded scope is agreed.
