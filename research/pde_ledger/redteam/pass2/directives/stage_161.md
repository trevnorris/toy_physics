---
unit_id: 161
batch: IV.6
created_at: 2026-06-08T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-08T17:50:06Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 161

Apply the finding below. After applying, append an `## Applied: F1` block under it with:
`files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a
question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Do NOT change the
set of asserted identities or their expected-zero targets — only the Mathematica derivation *route* so
the second engine is genuinely independent of the SymPy one.

After editing, RUN the script (`math -script /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl`)
and iterate until it exits 0 with every `expectZero` printing `PASS:` and the final
`Stage 161 Mathematica audit passed.` line. The printed numeric `(1+r_F1^2)/9 = 0.46236233468786880…`
and the `Delta_Q in (dThat,dS)` coefficients (`5.352238871696225`, `10.70447774339245`) must be
unchanged from the committed transcript so engine agreement with SymPy is preserved.

Do NOT touch paper.tex, notes/, or the SymPy script.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl:26-135`

**Issue:** The `.wl` is a line-by-line port of the `.py`: identical symbol definitions
(wl:31-32 ↔ py:40-41), the same `eps -> 0` derivative perturbation trick (wl:40-41 ↔ py:49-50), the
same closed forms `epsKExact`/`epsGExact` and the same `12 lW^2 -> Pi^2 a^2 (1+rc)` branch replacement
(wl:50,54,67,70), the same `upsilonPi`/`deltaQ` construction so the "collapsed" assertion cancels by
identical construction (wl:96-110 ↔ py:100-110), the same re-typed Stage-159 transport literals
(wl:112,115 ↔ py:117,122), and even a verbatim-ported inline comment (wl:72 "On the branch d gamma_0
= (1+r_c)/9 * d ln gamma_0." ↔ py:81). Only the `PolynomialRemainder` branch reduction (wl:57-61,77-87)
genuinely differs from SymPy's `.subs(poly,0)`. The second engine echoes the first engine's algebra
rather than re-deriving the result, which defeats the dual-engine independence guarantee.

**Required change (route only — keep every asserted identity and its expected zero identical):**

1. **Linearizations (`dBW`, and reuse downstream):** replace the `eps -> 0` derivative trick at
   wl:40-41 with a `Series`/`Normal` first-order expansion. E.g. introduce the perturbation
   `epsKappa -> eps*depsK`, `epsGamma -> eps*depsG`, `rc -> rcStar + eps*drc` and take
   `Normal[Series[bW, {eps, 0, 1}]]`, then read the coefficient of `eps` as `dBW`. Assert the same
   `expectZero["linearized slippage law", dBW - (1 + rcStar)*(depsG - depsK)/9]`.

2. **D/N-tube even/odd defects (wl:50-76):** reach the branch-reduced differentials via explicit
   logarithmic derivatives rather than re-typing `epsKExact`/`epsGExact` and the polynomial
   substitution. Derive `d eps_kappa` from `D[Log[12*lW^2/(Pi^2*a^2*(1+rc))], {…}]`-style log
   differentials (so the `2 dLW/lW - 2 da/a - drc/(1+rc)` target emerges from `dln`), and `d eps_gamma`
   from the logarithmic variation of `gamma0Sym/(1+rc)`. Keep the SAME asserted targets:
   `expectZero["d eps_kappa identity", …]` reducing to the target `2*dLW/lW - 2*da/a - drc/(1+rc)`, and
   `expectZero["d eps_gamma = d ln gamma0 - d ln(1+r_c)", …]`. You may keep `PolynomialRemainder` or
   use `Eliminate`/`Reduce` against `12 lW^2 == Pi^2 a^2 (1+rc)` — just do not re-type the SymPy
   `epsk_exact`/`epsg_exact` derivation choreography step-for-step.

3. **Collapse (wl:96-110):** build `upsilonPi = (1 + rcStar)*(xiGamma - 2*xiL)/9` and the Stage-160
   prefactor `9*sigmaStar/((1-sigmaStar)*(1+rcStar))` as separate named quantities, form
   `deltaQ = -prefactor * upsilonPi * dPiTan`, and assert the SAME
   `expectZero["collapsed Delta_Q law", deltaQ + sigmaStar*(xiGamma - 2*xiL)*dPiTan/(1 - sigmaStar)]`
   and the `N_Q-1` mirror. (This is the same identity; it just makes the prefactor an explicit factor
   rather than an inline literal, slightly hardening the cancellation too.)

4. **Remove the ported comment** at wl:72 (or rewrite it in the engine's own words); do not echo the
   SymPy comment text.

5. Leave the carry-forward `Print` block (wl:126-132), the `D/N SIMILARITY PRESERVATION` substitution
   checks (wl:118-120), and the printed `(1+r_F1^2)/9` / `Delta_Q in (dThat,dS)` results unchanged in
   VALUE.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 161` and confirm the script exits
0 with all `expectZero` printing `PASS:`, the printed prefactor/coefficients unchanged (engine
agreement with SymPy preserved), and the `.wl` no longer mirroring the SymPy line structure (no
`eps -> 0` derivative trick, no verbatim "On the branch d gamma_0 …" comment).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl`
- summary: Reworked the Mathematica audit to use series linearization, logarithmic defect variations, and an explicit Stage-160 prefactor while preserving the asserted identities and printed values.
- deviation: none
</content>
