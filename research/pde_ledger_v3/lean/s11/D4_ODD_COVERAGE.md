# D4 odd-density divergence and bulk variation contract

Authorized after the D4 classification checkpoint `d6119f75`: the user asked
to commit the completed classification and proceed to its stated next increment.
Governed by [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).
**Complete for D4B.1–D4B.4.** Local verification and current source/object
correspondence pass; Claude and Grok independently returned CLEAR with no
required corrections. See [D4_ODD_FIDELITY_REVIEW.md](D4_ODD_FIDELITY_REVIEW.md).

The gap is the bulk meaning of the unique D4 reflection-odd density already
classified in D4.1–D4.4. Reuse that classification and S10's smooth-field,
compact-variation and integration-by-parts machinery. This increment concerns
only the odd family, not the full four-parameter D4 bulk classification.

| Item | Finite deliverable |
|---|---|
| D4B.1 — supplied density and actual variation | Use the reviewed P with `P_D=P` and `L_beta(J)=-(beta/2) P(G)`, `G_ij=J_(i+1),j`. Identify all SO-invariant reflection-odd densities through the existing classification. Derive momenta from actual derivatives of L and the finite relative-action first variation. |
| D4B.2 — explicit divergence | Give a current K with `sum_i partial_i K_i=P(grad u)` on every smooth field. Retain its exact factor/sign and show P can be nonzero; a divergence is not pointwise zero. |
| D4B.3 — exhaustive odd-family bulk conclusion | For every real constant beta, every smooth background and every smooth compact variation, prove the local Euler–Lagrange expression and the actual first variation vanish. Include beta=0, negative beta and arbitrary time dependence. No generic-wave or polarization restriction is used. |
| D4B.4 — compact fidelity and closure | Reuse the reviewed D4 Q9 identification and compare the actual native odd density and its derivative-defined momentum/local operator. Check native V5 on that odd combination only, with sign/factor controls on nonzero momenta and a wrong non-null density control. Verify standard axioms, meaningful mathematical mutations and positives, freeze a compact packet, and obtain two independent non-author fidelity reviews. |

## Domain and conventions

Coordinates are `(t,x1,x2,x3,x4)`, with `G_ij=partial_i u_j` and spatial
indices 0–3 corresponding to the four spatial coordinates. Backgrounds are
smooth real fields on R^(4+1). Variations are smooth with compact support in
spacetime. Only the change in density is integrated, so backgrounds need not
have finite total action. Beta is an arbitrary constant real coefficient with
the units needed by the supplied action; physical units are not separately
formalized. The inertia and other even densities are not added here.

Write `F_ij=G_ij-G_ji` and
`P=F01 F23-F02 F13+F03 F12`. The complementary antisymmetric matrix is

```
M = [[ 0,    F23, -F13,  F12],
     [-F23,  0,    F03, -F02],
     [ F13, -F03,  0,    F01],
     [-F12,  F02, -F01,  0  ]].
```

The proved current is `K_i=(1/2) sum_j u_j M_ij`. The contraction
`sum_ij G_ij M_ij=2P` and the mixed-derivative cancellation
`sum_i partial_i M_ij=0` establish the current identity. The Lagrangian momenta
are `-(beta/2) M_ij`, with zero time row. Lean's variational expression is
`-sum_j partial_j(dL/dJ_ji)`; the native helper uses the opposite overall sign.
Since the final bulk expression is zero, sign/factor sensitivity must also be
checked on the nonzero momenta, not inferred from zero-equals-zero.

The local divergence identity retains a possible boundary contribution.
Compact support justifies the bulk first-variation conclusion; this contract
does not claim absence of boundary effects or state interface conditions.

## Evidence and stopping rule

Recorded run2 passed four new modules and their audit root, with fourteen
unchanged imports reused under full source/pin/generator/input/output object
guards. All 41 selected declarations use only standard axioms; twelve
mathematical mutations were rejected and twelve positives passed. The compact
native check and current source/object correspondence also pass. See
[D4_ODD_VERIFICATION.txt](D4_ODD_VERIFICATION.txt) and
[D4_ODD_FIDELITY.md](D4_ODD_FIDELITY.md). D4B.1–D4B.3 are proved and the
independent reviews in D4B.4 are complete with no unresolved blockers.

Existing evidence: D4.1–D4.4 at `d6119f75`, including `orientationForm_apply`,
`odd_classification`, exact native `P_D=P`, and epsilon contraction `2P`.
S10 provides actual coordinate derivatives and compact integration by parts.
The controls detect a wrong density/current factor, wrong momentum sign,
a failed mixed-derivative cancellation, nonzero bulk claimed for the canonical
odd term, and the confusion between zero bulk variation and zero density.
Paired true statements and admissible nonzero, zero and negative coefficient
examples pass. Instrument errors and timeouts are not mathematical rejections.

The user approved the fixed 38-file D4B packet with “Yes you can”. Both
independent reviews cleared that same packet. The compact native connection,
controls and current source/object/pin correspondence pass. Optional reviewer
suggestions do not close missing contract items and are dispositioned in the
review record. D4B.1–D4B.4 are complete; the stopping rule now applies.

| Closed obligation | Principal evidence |
|---|---|
| D4B.1 | `all_odd_densities`, `density_identity`, derivative-defined `momentum_eq`, `lagrangian_variation`, `relative_density_integrable`, `relativeAction_hasDerivAt` |
| D4B.2 | `dualCurl_contraction`, `dualCurl_divergence_zero`, `boundary_identity`, `lagrangian_is_divergence`, `nonzero_density`, `current_normalization` |
| D4B.3 | `relativeAction_deriv_eq_eulerLagrange`, `eulerLagrange_zero`, `firstVariation_zero`, `every_background_stationary` |
| D4B.4 | Compact native and formal check reports, twelve mathematical rejections, twelve positives, both CLEAR fidelity reports and `S11_lean_d4_odd_closure_validation.json` |

Preserve all completed H/I/E/J/K and D4 classification sources and historical
records. No D5, full D4 even-family bulk census, general null-Lagrangian
classification, variable beta, spectra/stability/interfaces, S11c calculations
or pinned exports, production reruns or systematic CAS bridge. Use one Lean
worker, memory-conscious sequential jobs and silent completion/error hooks.
