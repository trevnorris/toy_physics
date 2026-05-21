---
unit_id: 013
batch: I.2
created_at: 2026-05-21T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-21T12:49:15-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 013

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl` (create new file)

**Issue:** Unit 013 has no Mathematica counterpart. Stages 010, 011, 012 each carry one; the policy requires both engines for non-status-only, non-checkpoint units. The new `.wl` must independently verify the same six physical claims the SymPy script verifies, without transliterating the SymPy code. In particular, the chain-rule expansions `z0, z2, z4, n0, n2, n4` must be derived from a master primitive (`Q/(Delta - S2*ell^2 - Hport*ell^4)` for the z-tower; `P/(Delta - S2*ell^2 - Gw*ell^4)` for the n-tower) using Mathematica's `Series` rather than written down as literal polynomials.

**Required change:**

Create the file at the Target path with a script that performs the following claim manifest. Use Mathematica idioms (`Series`, `Coefficient`, `FullSimplify`, `Assuming`); do not import or echo the literal `z0/z2/z4/n0/n2/n4` polynomials from the `.py` file. Each claim should produce a `Print` of `OK <label>` on success and call `Exit[1]` on failure (use a helper like `assertZero[label_, expr_] := If[FullSimplify[expr] =!= 0, Print["FAIL ", label, ": ", expr]; Exit[1], Print["OK ", label]]` and a matching `assertNonzero`).

**Claim manifest:**

M1. **One-sided Taylor first moment (W = exp(-u)):** With `W = Exp[-u]` and `X = X0 + ell*u*X1 + ell^2*u^2*X2/2`, compute `Xproj = Integrate[W*X, {u, 0, Infinity}]` under `Assumptions -> {Re[u] > 0}` or use `Limit`/`Series` as needed. Assert that `Normal[Series[Xproj, {ell, 0, 1}]] == X0 + ell*X1`. (Equivalently, `Coefficient[Xproj, ell, 0] - X0 == 0` and `Coefficient[Xproj, ell, 1] - X1 == 0`.)

M2. **Second-moment closed form for W2 = u*exp(-u):** Assert `Integrate[u^2*Exp[-u], {u, 0, Infinity}] - 2 == 0` and the corresponding first-moment recovery: with `Xproj2 = Integrate[u*Exp[-u]*X, {u, 0, Infinity}]`, assert `Normal[Series[Xproj2, {ell, 0, 1}]] - (X0 + 2*ell*X1) == 0`.

M3. **Chain-rule derivation of the z-tower (M3 is the substantive replacement for hardcoded z0, z2, z4):** Define a master `Zsource[ell_] := Q[t] / (Delta[t] - S2[t]*ell^2 - Hport[t]*ell^4)` (or whatever local form matches the master notes; for this audit, use `Q[t]/(Delta[t] - S2[t]*ell^2 - Hport[t]*ell^4)` since that reproduces the polynomials in the SymPy script). Compute `Zexpansion = Series[Zsource[ell], {ell, 0, 4}]`. Take time-derivatives at `t -> 0` of the order-`ell^0`, `ell^2`, `ell^4` coefficients, then specialize the leading values to `Q[0] -> Q, Delta[0] -> Delta, S2[0] -> S2, Hport[0] -> Hport` and the t-derivatives to `Q'[0] -> q1, Delta'[0] -> d1, S2'[0] -> s1, Hport'[0] -> h1`. Assert that the resulting expressions equal the literal SymPy forms

  z0 expected: `(Delta*q1 - Q*d1)/Delta^2`
  z2 expected: `(-Delta^2*h1 + Delta*(Hport*d1 + Q*s1 + S2*q1) - 2*Q*S2*d1)/Delta^3`
  z4 expected: `(-Delta^2*Hport*s1 - Delta^2*S2*h1 - Delta^2*q1 + 2*Delta*Hport*S2*d1 + 2*Delta*Q*S2*s1 + 2*Delta*Q*d1 + Delta*S2^2*q1 - 3*Q*S2^2*d1)/Delta^4`

(Use `FullSimplify[derived - expected] == 0`.) If the master form does not actually reproduce these polynomials, do NOT silently change the expected values; instead `Exit[1]` and the verifier will treat it as a finding back to the auditor.

M4. **Chain-rule derivation of the n-tower:** Analogous to M3 but for `Nsource[ell_] := 2*P[t]^2 / (Delta[t] - S2[t]*ell^2 - Gw[t]*ell^4) - 2*P[t]^2/Delta[t]` (or the appropriate `P`-shift primitive — derive by inspection from the structure of n0, n2, n4 in the SymPy script and assert against the SymPy polynomials):

  n0 expected: `2*P*(Delta*p1 - P*d1)/Delta^3`
  n2 expected: `-(2*Delta^2*(Gw*p1 + P*g1) - 2*Delta*P*(2*Gw*d1 + P*s1 + 2*S2*p1) + 6*P^2*S2*d1)/Delta^4`
  n4 expected: `2*(Delta^3*Gw*g1 - Delta^2*Gw^2*d1 - 2*Delta^2*Gw*P*s1 - 2*Delta^2*Gw*S2*p1 - 2*Delta^2*P*S2*g1 - 2*Delta^2*P*p1 + 6*Delta*Gw*P*S2*d1 + 3*Delta*P^2*S2*s1 + 3*Delta*P^2*d1 + 3*Delta*P*S2^2*p1 - 6*P^2*S2^2*d1)/Delta^5`

If you cannot find a master primitive whose `Series` reproduces these, fall back to asserting the polynomial identities directly and append `## Blocked: F1.M4` to this directive explaining what primitive you tried.

M5. **Substantive coefficient checks (no transliteration):** Define `Xi, K1, H_even, deltaP2_der, deltaP4_der` from the derived z/n towers and the same combinations the SymPy script uses:

  `Xi = (2*p1/P - 2*d1/Delta + q1/(D0*Delta) - Q*d1/(D0*Delta^2)) /. {q1 -> mu1*Qx, ..., d1 -> mu1*Dx, ...} / mu1`
  `K1 = -(z2 + z0/9) /. (same subs) / mu1`
  `Heven = (-z4 + (2/3)*z2 - z0/27) /. (same subs) / mu1`

Assert:
  - `D[Xi, Px] - 2/P == 0` (matches SymPy line 91)
  - `D[deltaP2_der, Gx] + 2*P/(D0*Delta^2) == 0` (matches SymPy line 109)
  - `FullSimplify[D[deltaP4_der, Gx]] =!= 0` (matches SymPy lines 110-111)

M6. **Mechanism sieve:** With the same K1, H_even from M5, define the 2x2 source/denominator submatrix at `Sx -> 0, Hx -> 0` (derivatives w.r.t. `{Qx, Dx}`) and the 2x2 spectral submatrix at `Qx -> 0, Dx -> 0` (derivatives w.r.t. `{Sx, Hx}`). Assert both determinants are `FullSimplify`-nonzero, and that `Solve[{K1 /. {Sx -> 0, Hx -> 0} == 0, Heven /. {Sx -> 0, Hx -> 0} == 0}, {Qx, Dx}]` returns only the trivial solution `{Qx -> 0, Dx -> 0}` (and the analogous spectral solve returns only `{Sx -> 0, Hx -> 0}`).

End the script with a final `Print["STAGE 013 MATHEMATICA AUDIT: PASS"]` once all assertions pass.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 013` and confirm (a) the file exists at the Target path, (b) the script exits 0, (c) the transcript contains `OK` lines for M1-M6 and the final PASS print.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl`
- summary: Added the missing Mathematica audit witness for the Taylor moments, chain-rule towers, coefficient gates, and mechanism sieve.
- deviation: Used the stage-012 master primitives that reproduce the cited polynomials, `(Q-Hport*ell^2)/(Delta-S2*ell^2+ell^4)` and `(P-Gw*ell^2)^2/(Delta-S2*ell^2+ell^4)^2`, because the directive's denominator-only primitive does not reproduce the expected towers.

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:86-90`

**Issue:** The five assertions

```
for sym in (Sx, Hx, Gx):
    assert_zero(f"Xi independence from {sym}", sp.diff(Xi, sym))
for sym in (Px, Gx):
    assert_zero(f"even gates independence from {sym}", sp.diff(K1, sym))
    assert_zero(f"even gates independence from {sym}", sp.diff(H_even, sym))
```

each compute `sp.diff(EXPR, VAR)` where `VAR` is not a free symbol of `EXPR`. SymPy returns 0 unconditionally; `assert_zero` cannot fail. These five lines purchase "PASS" cheaply.

**Required change:**

Delete lines 86 through 90 inclusive. The two `for sym in ...:` loops and their bodies — total of five logical assertions across the loop iterations — must be removed.

Before:

```python
    for sym in (Sx, Hx, Gx):
        assert_zero(f"Xi independence from {sym}", sp.diff(Xi, sym))
    for sym in (Px, Gx):
        assert_zero(f"even gates independence from {sym}", sp.diff(K1, sym))
        assert_zero(f"even gates independence from {sym}", sp.diff(H_even, sym))
    assert_zero("dXi/dPprime", sp.diff(Xi, Px) - 2 / P)
```

After:

```python
    assert_zero("dXi/dPprime", sp.diff(Xi, Px) - 2 / P)
```

(The substantive `dXi/dPprime` line at line 91 is preserved verbatim.) Do not add replacement assertions; that is scope expansion. Do not modify any other line.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 013` and confirm (a) lines 86-90 (the two independence loops) are absent, (b) the `dXi/dPprime` assertion remains, (c) the script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py`
- summary: Removed the tautological independence assertion loops while preserving the substantive `dXi/dPprime` check.
- deviation: none

## F3 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:49-80`

**Issue:** The chain-rule expansion expressions `z0, z2, z4, n0, n2, n4` are written as literal symbolic polynomials with no in-script derivation and no precise upstream anchor. The docstring at line 2 cites a notes-file name but does not name an upstream script + assertion label that *verifies* those polynomials. The downstream assertions (sieve, determinant, coefficient checks) all sit on top of these hardcoded forms and cannot detect an error inside z0..n4.

**Required change:**

Insert a comment block immediately above line 49 (the `z0 = ...` line). The block must name (a) the master-primitive form that generates the z- and n-towers and (b) either an upstream stage script + line/label that already verifies those polynomials, OR a forward reference to the M3/M4 Mathematica derivation introduced by F1.

Concretely, insert the following block at line 49 (pushing the existing `z0 = ...` down):

```python
    # z0, z2, z4 are the order-ell^0, ell^2, ell^4 coefficients of the
    # bottleneck chain-rule expansion Q(t)/(Delta(t) - S2(t)*ell^2 -
    # Hport(t)*ell^4), evaluated at t=0 with leading values
    # {Q, Delta, S2, Hport} and t-derivatives {q1, d1, s1, h1}. The
    # polynomial forms below are independently derived and asserted
    # against the master primitive in
    # mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl
    # (claims M3 and M4); a discrepancy there is the canonical place to
    # discover an error in these literals.
```

Insert a parallel comment block immediately above the `n0 = ...` line (currently line 62) noting that n0, n2, n4 are the analogous coefficients of `2*P(t)^2/(Delta(t) - S2(t)*ell^2 - Gw(t)*ell^4) - 2*P(t)^2/Delta(t)` (or the matching primitive used in F1.M4), again citing the same `.wl` file.

If F1 is `## Blocked` (no Mathematica file lands), this anchor is invalid and Codex must instead substitute a citation to whatever upstream `.py` stage in `scripts/` already verifies these polynomials. If no such upstream verification exists, append `## Blocked: F3` to this directive with the question "no upstream verification found for z0/z2/z4/n0/n2/n4 — please specify the master notes derivation reference for the comment block."

**Verification command:**
After Codex applies, the verifier will read the script and confirm the two comment blocks are present immediately above the z0 and n0 definitions and cite the `.wl` file by name. The script's executable behavior is unchanged; `redteam exec-sympy 013` must still exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py`
- summary: Added upstream Mathematica derivation anchors immediately above the hardcoded z- and n-tower literals.
- deviation: The z-tower comment names the primitive used by the new Mathematica derivation rather than the directive's denominator-only form, which does not reproduce the listed z-polynomials.
