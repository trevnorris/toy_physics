# Review: Stage 002 — Breathing reduction

**Batch:** 1 — Geometry Lift & Coupling
**Status:** Hardened and verified (2× PASS, 2026-04-21)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage002_breathing_reduction.md`
- **SymPy:** `scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`

## Review Checklist

- [x] Equation-level correctness (signs, factors, indices, limits)
- [x] Logical flow from prior stage(s)
- [x] Assumptions stated and justified
- [x] Notation consistent with prior stages
- [x] Physical interpretation sensible
- [x] SymPy script faithfully implements notes
- [x] Mathematica script faithfully implements notes
- [x] Script runs without error
- [x] Script output matches notes claims
- [x] No missing edge cases or branches

## Hardening Summary

The earlier Stage `002` audits were too generic. They checked:

1. the `Y_00` normalization bridge,
2. a generic quadratic Euler-Lagrange matrix identity,
3. and the arithmetic specialization `ell(ell+1)|_{ell=2} = 6`.

Those checks were not enough to support the actual stage claim, which is that
the Stage `001` wall action reduces to the boxed overlap integrals for the
`(a,L)` truncation and that the grouped real `P_2` sector is genuinely
degenerate on the isotropic reference branch.

That gap is now closed.

## What the hardened scripts now prove

### 1. The `(a,L)` overlap matrices are derived from the Stage `001` action

Both CAS layers now start from the Stage `001` densitized quadratic wall action
on the axisymmetric branch and insert the Stage `002` ansatz

`eta_0 = 2 sqrt(pi) [alpha_a(w) delta a + alpha_L(w) delta L] Y_00`.

They integrate over the sphere and recover the exact reduced Lagrangian density
with the boxed overlap-integral matrices

- `M_AB = 4 pi ∫ dw mu_eta alpha_A alpha_B`,
- `K_AB = 4 pi ∫ dw [T_w alpha_A' alpha_B' + K_0 alpha_A alpha_B]`.

This is now a real action reduction, not a generic quadratic-form placeholder.

### 2. The conservative matrix system is checked with the actual overlap coefficients

After the overlap matrices are derived, both CAS layers build the reduced
time-domain Lagrangian using those same formal integral coefficients and verify
that the Euler-Lagrange equations are exactly

`M_AB Q¨^B + K_AB Q^B = 0`.

So the matrix equation is now tied to the derived overlaps rather than to
abstract `M_aa`, `K_aL`, etc.

### 3. The grouped real `P_2` degeneracy is supported by harmonic structure

The SymPy audit proves the full real-`P_2` Gram and angular-stiffness matrices:

- norm matrix = identity,
- angular matrix = `6 I`.

The Mathematica mirror uses a lighter but still substantive route:

- exact phase-shift identities link the sine harmonics to the cosine
  representatives,
- representative real `P_2` modes `Y20`, `Y21c`, and `Y22c` each have norm `1`,
  angular energy `6`, and satisfy `-Delta_{S^2} Y = 6 Y`,
- and each representative gives the same one-mode reduced density
  with stiffness shift `K_eta + 6 T_Omega`.

Together, that closes the isotropic grouped-real `P_2` degeneracy claim in both
CAS layers without importing any hidden constants.

## Verification Outcome

Both hardened scripts passed on `2026-04-21`:

- `python3 research/pde_ledger/scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py`
- `math -script research/pde_ledger/mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`

The Mathematica `P_2` section needed a lighter proof path than the full five-by-five
packet matrix, but the final mirror still proves the intended degeneracy claim.

## Verdict

**PASS**

Stage `002` is no longer relying on a generic quadratic-Lagrangian surrogate.
The dual-CAS audits now derive the overlap-integral closure from the Stage `001`
wall action itself and explicitly support the grouped-real `P_2` degeneracy in
the isotropic uncoupled limit.
