# Review: Stage 039 — Transport source asymmetry

**Batch:** 6 — Support & Threshold Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage039_transport_source_asymmetry.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage039_transport_source_asymmetry_sympy_audit.py`

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Script runs without error
- [ ] Script output matches notes claims
- [ ] No missing edge cases or branches

## Agent Reviews

<!-- Agents: append your review below this line using the template:

### Agent: [Model Name] — [Date]
**Verdict:** [PASS | MINOR | ISSUE | BLOCK]

**Notes Derivation Review:**

**Script Review:**

**Issues Found:**

**Questions:**

### Agent: Claude Opus 4.6 — 2026-04-02
**Verdict:** PASS

**Notes Derivation Review:**

1. **Equation-level correctness.** Drift-diffusion J = -D_sigma sigma' + v_sigma sigma with J=0 gives sigma = C exp(Pe s/L), Pe = v_sigma L/D_sigma. Normalization C = Pe/(exp(Pe)-1) verified. Identification alpha = Pe exact. Covariance identity dOmega_Pe/dPe = Cov_Pe(chi_0,s)/I_W verified (script independently computes derivative and covariance from three integrals). Monotonicity dOmega/dPe > 0 from FKG (covariance of co-monotone functions). Endpoints: Omega(0)=1, Omega(∞)=pi/2. Small-Pe coefficient (4-pi)/(2pi). Large-Pe: pi/2 - pi^3/(8Pe^2).

2. **Logical flow.** Clean: transport law → stationary branch → exponential profile → identification with Stage-36 family → covariance monotonicity → asymptotics.

3. **Physical interpretation.** Sound: converts abstract deformation parameter alpha into physical Peclet number. Pe > 0 = constructive (drift toward antinode).

**Script Review:** Zero-flux residual from ODE solution. Normalization by integration. Omega from actual sp.integrate vs closed form. Covariance identity: sp.diff(Omega,Pe) vs independently computed Cov from three separate integrals — strongest non-tautological check in Batch 6. Both asymptotics by sp.series. All pass (exit code 0). Exemplary verification depth.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

- The transport-origin derivation is correct. The stationary zero-flux solution of
  `partial_t sigma + partial_s J = 0`, `J = -D_sigma sigma' + v_sigma sigma`
  gives the normalized exponential branch exactly, and the identification `Pe = v_sigma L / D_sigma` is clean.
- I independently rechecked the overlap formula against the D/N lowest mode. The closed form for `Omega_Pe` matches the integral result exactly, and the endpoint limits `Omega_Pe(0)=1` and `lim_(Pe->oo) Omega_Pe = pi/2` are correct.
- The covariance identity is the key nontrivial claim here, and it holds: differentiating the exact overlap boost gives the covariance over `I_W` as stated. That makes the monotonicity argument on the constructive branch credible rather than rhetorical.
- Physically, the stage does what it should: it converts the earlier abstract source-shape asymmetry into a concrete axial transport Peclet number, so Stage 36’s exponential family is no longer a free knob.

**Script Review:**

- The audit script is faithful to the note and the cached output matches it line for line on the checks that matter: zero-flux residual, normalization, exact overlap, covariance identity, and asymptotics.
- I did not find a coding issue or a tautological proof. The use of separate integrals for the covariance check is the right level of independence.

**Issues Found:**

None.

**Questions:**

None.

---
