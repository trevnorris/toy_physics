# D5 constant-coefficient bulk contract (D5B.1–D5B.4)

Status: **bounded D5B.1–D5B.4 complete**. Local verification passed and both
independent fidelity reviews returned CLEAR with no required correction. Governed by `../FORMALIZATION_POLICY.md`. This increment
reuses the completed D5 density classification at `ec27ecf0` and the existing
S10 variation machinery, following the completed D3/D4 bulk contracts.

## Claim and object

The entire classified D5 family is
`Q=a(tr G)^2+b tr(G^2)+c tr(G G^T)`, `L=-Q/2`, with arbitrary constant real
`a,b,c`, including zero and negative values. `G_ij=partial_i u_j` has all 25
independent real entries; fields have five components and five spatial
coordinates, plus time. Lean's row zero is time and row `i.succ` is spatial
coordinate `x_(i+1)`. There is no kinetic term or time momentum. This is a
conditional theorem about the supplied action, not a derivation of dimension
or physical coefficients.

| Item | Finite completion obligation |
|---|---|
| D5B.1 | Identify `L` with the reviewed D5 invariant form; use actual jet derivatives for momentum and prove the derivative of the finite relative action on smooth backgrounds with smooth compact variations. Stationarity is equivalent to the pointwise local EL equation. |
| D5B.2 | Prove `EL=(a+b)grad(div u)+c Delta u`; universally equal bulk responses iff `c` and `a+b` agree. Prove the null family is precisely `(t,-t,0)`, the response dimension is two and kernel dimension one. |
| D5B.3 | Prove the explicit current `J_i=sum_j(u_i partial_j u_j-u_j partial_j u_i)` satisfies `div J=(div u)^2-tr(G^2)`; hence the null Lagrangian is `-t div J/2`. Establish the compact plane-wave/homogeneous map `rho=0,mu=c,B=a+b+c`, including zero wavevector. |
| D5B.4 | Identify all three native Q9/V5 basis responses, actual action momentum and sign/factor conventions using selected existing helpers; retain nonzero density/momentum/current witnesses, independent response witnesses and paired mathematical mutation controls with positives. |

Bulk equivalence quantifies over every smooth background and point; nullness
quantifies over every smooth background and smooth compact variation. The
full coefficient space is retained. Null/non-null cases are complementary;
no generic witness substitutes for the universal equivalence theorem. The
native helper uses `+div(momentum)`, opposite the Lean EL convention for
`L=-Q/2`; its V5 of `Q` is twice that Lean EL. The compact source check must
verify those conventions on the actual three-dimensional native span.

## Evidence and review

Verified modules: Action, Variation, Bulk, Census and Controls, plus
an audit root. Reuse S10 analytic lemmas and the completed D5 classification;
do not change their reviewed source. The already dimension-general coordinate identities from `S11D4Odd.Calculus`
are reused directly. No new invariant
certificate is generated. The unchanged D5 generator is checked with `--check`.

The source comparison is a tested translation outside Lean's kernel. No
production driver or S11c source is executed. Controls must reject for an
explicit false mathematical statement, never a timeout, warning, import or
syntax error. Actual nonzero current/momentum witnesses protect normalization
which a zero null-family bulk expression alone cannot detect.

The user authorized transfer of this D5B packet to Claude and Grok on
2026-09-18, and a commit after both approve. Those conditions are satisfied.
The fixed packet and both substantive CLEAR reports are documented in
[D5_BULK_FIDELITY_REVIEW.md](D5_BULK_FIDELITY_REVIEW.md), including the
optional-note dispositions and evidence limits. Strict builds, the standard-axiom
audit, positive/negative controls and source/object correspondence all passed.

## Exclusions and stopping rule

No variable coefficients, general-dimensional theorem, general null-Lagrangian
classification, interface theory, new spectrum/stability/scattering/pole work,
S11c calculations/exports, production reruns or systematic CAS bridge. A null
bulk response does not imply pointwise zero density or absence of boundary
effects. Existing historical proofs/reports and protected objects remain
historical. Stop when D5B.1–D5B.4 and their two reviews are complete.
