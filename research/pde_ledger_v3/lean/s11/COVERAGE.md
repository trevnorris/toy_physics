# S11 homogeneous compression: bounded Lean contract

Started 2026-09-16 UTC, following the user's approval of the first S11 increment.
Status: **H1–H4 complete**, 2026-09-16 UTC. Both independent fidelity reviews
returned CLEAR; see [FIDELITY_REVIEW.md](FIDELITY_REVIEW.md).
Governed by [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).

This contract concerns the selected homogeneous curl-plus-compression action.
It does not encompass the whole S11 series. Reuse the S9/S10 density, calculus
and subspace results; add only the combination and its classification.

| Item | Finite deliverable |
|---|---|
| H1 — Action and identification | For supplied constant coefficients, identify the actual density `rho/2 |u_t|² - mu/2 S_curl(du) - B/2 (div u)²`, with the same antisymmetric double-sum normalization as S10. Reuse the curl and divergence-only variational proofs to derive the combined modal operator and connect it to integrated stationarity. Check the compact native SymPy/Wolfram MAIN constructor and normalization correspondence. |
| H2 — Exhaustive modes | On the physical D=3 domain `rho,mu,B > 0`, `k ≠ 0`, prove the entire amplitude kernel: the transverse plane at squared frequency `mu |k|²/rho` and longitudinal line at `B |k|²/rho` when `B ≠ mu`; the full three-space at the common root when `B = mu`; zero at every other real frequency. Count full subspaces, not witness vectors. Retain both frequency signs and positive-frequency existence. Explicitly prove the `B=0` recovery of S10 and treat `k=0` separately from the nonzero-wavevector census. General-D lemmas may be reused internally; D=3 is the completion claim. |
| H3 — Kinematic bulk threshold | With supplied `c_s > 0`, classify the sign of `k_w² = |k|² (B/(rho c_s²) - 1)` as negative, zero or positive according as `B` is below, equal to or above `rho c_s²`. Identify this with phase matching on the longitudinal branch. This is not a bound-state, coupling, radiation or leakage theorem. |
| H4 — Controls and fidelity | Compile sequentially, audit the load-bearing theorems without admissions/custom physics axioms, and reject isolated mathematical mutations with passing controls. Cover the compression sign/normalization, transverse independence, longitudinal coefficient, coincidence/full-kernel count, off-root exclusion, zero-compression/zero-wavevector limits and threshold sign/domain. Obtain two independent non-author reviews of a fixed packet and resolve findings. |

The analytic setting remains smooth real vector fields on flat whole spacetime,
with smooth compactly supported variations; finite total background action is
not required. Coordinates are `(t,x1,...,xD)`, `J j i = partial_j u_i`, and the
phase is `k·x - omega t`. The physical action and the material meaning of `u`
remain supplied premises. `rho` maps to `rho_br`/`rhoBr`, `mu` to `mu_R`/`muR`,
`B` to `B_comp`/`bComp`, and `c_s` to `c_s0`/`cs0`. Source spellings and modal
route scalar factors have been checked against the actual constructors at the
level recorded in [FIDELITY.md](FIDELITY.md).

The physical mode contract assumes positive compression. The `B=0` theorem is
a separately identified boundary limit, not an admissible positive-compression
sample. The whole homogeneous classification includes the coincidence locus;
no generic chart may remove it. The geometrical longitudinal/transverse spaces
remain meaningful at coincidence even though a basis of the full eigenspace
may mix them. Units are the inherited S10 inertia/stiffness units, with `B` in
the same units as `mu`; `B/rho` and `c_s²` have squared-speed units.

The native source connection is an inspected/tested translation boundary, not
a Lean-certified interpreter or whole-export certificate. Existing S11 CAS
stratum/comparator debts remain separate. No systematic emission bridge is
planned. The incorrect historical phrase “a longitudinal wave is pure trace”
has been corrected in the step record: its gradient is symmetric but can contain
both trace and symmetric-traceless parts. Curl-freeness removes the shear
contribution; compression acts through the divergence.

[VERIFICATION.txt](VERIFICATION.txt) records the passing builds, axiom audit,
mathematical mutations and positive controls, with unchanged verified sources.
Both independent non-author reviews of the fixed packet are complete, with no
mathematical blockers or required control gaps. Their nonblocking suggestions
and dispositions are recorded in [FIDELITY_REVIEW.md](FIDELITY_REVIEW.md).

Excluded: SO(D)/O(D) invariant completeness and exceptional dimensions; energy
bases modulo total divergences; other S11 packages; S11b interface/passivity;
S11c operators, calculations or files; variable coefficients, finite boundaries,
nonlinear physics, confinement, observability and general PDE/Fourier completeness.
The invariant-count repairs and later interface/scattering results need their
own separately agreed contracts. Stop when H1–H4 close.
