---
unit_id: 009
batch: I.1
created_at: 2026-05-21T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-21T11:12:11-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 009

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script

**Target:** create `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl`

**Issue:** The manifest declares unit 009 non-status-only and non-checkpoint, so both SymPy and Mathematica verification scripts are required. Only the SymPy script exists. A second engine must independently re-derive the projection identities (integration by parts on the half-line mouth kernel, even-kernel Taylor moments, effective-parameter ratios, Gaussian-localizer asymptotic) so that no claim depends solely on SymPy's `integrate`/`series` behavior.

**Required change:**
Create the file with this structure. Use Mathematica idioms throughout; do NOT line-by-line transliterate the SymPy code. Use `Integrate`, `Series`, `Normal`, `FullSimplify`, `Assuming`, and explicit `If[FullSimplify[lhs - rhs] =!= 0, Print["FAIL: …"]; Exit[1]]` guards. The script must terminate with `Print["STATUS: PASS"]; Exit[0]`.

Recommended skeleton (the implementer should choose Mathematica-idiomatic forms rather than copying SymPy expressions):

```mathematica
(* moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl *)
(* Independent Mathematica verification of unit 009. *)

assertZero[label_String, expr_] := If[
  FullSimplify[expr] =!= 0,
  Print["FAIL: ", label, " :: ", FullSimplify[expr]]; Exit[1],
  Print["OK: ", label]
];

(* --- M1: half-line integration-by-parts recombination --- *)
ClearAll[w, ell, q0, q1, q2, q3, q4];
Qw  = q0 + q1 w + q2 w^2/2 + q3 w^3/6 + q4 w^4/24;
Wel = Exp[-w/ell]/ell;
$Assumptions = ell > 0;
avgQ   = Integrate[Wel*Qw, {w, 0, Infinity}];
avgDQ  = Integrate[Wel*D[Qw, w], {w, 0, Infinity}];
bdry   = (Limit[Wel*Qw, w -> Infinity]) - (Wel*Qw /. w -> 0);
minusWp = -Integrate[D[Wel, w]*Qw, {w, 0, Infinity}];
assertZero["M1a half-line IBP recombination", bdry + minusWp - avgDQ];
assertZero["M1b half-line dQ expansion",
  avgDQ - (q1 + ell*q2 + ell^2*q3 + ell^3*q4)];
assertZero["M1c half-line Q expansion",
  avgQ - (q0 + ell*q1 + ell^2*q2 + ell^3*q3 + ell^4*q4)];

(* --- M2: even-kernel Taylor moments (unit Gaussian) --- *)
ClearAll[u, sigma, q0, q1, q2, q3, q4, q5];
Qser = q0 + q1*(sigma u) + q2*(sigma u)^2/2 + q3*(sigma u)^3/6
          + q4*(sigma u)^4/24 + q5*(sigma u)^5/120;
Wu = Exp[-u^2/2]/Sqrt[2 Pi];
avgQgauss  = Integrate[Wu * (Qser /. sigma -> 1), {u, -Infinity, Infinity}];
avgDQgauss = Integrate[Wu * D[(Qser /. sigma -> 1), u], {u, -Infinity, Infinity}];
assertZero["M2a Gaussian even-kernel Q",  avgQgauss  - (q0 + q2/2 + q4/8)];
assertZero["M2b Gaussian even-kernel dQ", avgDQgauss - (q1 + q3/2 + q5/8)];

(* --- M3: zero-mode effective-parameter series, mouth case --- *)
ClearAll[ell, mu0, xi, s0, s1, s2, z0, z1, z2, h0, h1, h2];
muHalf = mu0*(s0 + ell*s1 + ell^2*s2)/(z0 + ell*z1 + ell^2*z2);
xiHalf = xi *(z0 + ell*z1 + ell^2*z2)/(h0 + ell*h1 + ell^2*h2);
muHalfSer = Normal[Series[muHalf, {ell, 0, 1}]];
xiHalfSer = Normal[Series[xiHalf, {ell, 0, 1}]];
assertZero["M3a mouth mu_eff to O(ell)",
  muHalfSer - (mu0*s0/z0 + ell*(mu0*s1/z0 - mu0*s0*z1/z0^2))];
assertZero["M3b mouth xi_eff to O(ell)",
  xiHalfSer - (xi*z0/h0 + ell*(xi*z1/h0 - xi*z0*h1/h0^2))];

(* --- M4: zero-mode effective-parameter series, symmetric case --- *)
ClearAll[sigma, m2, mu0, xi, s0, s2, z0, z2, h0, h2];
muSym = mu0*(s0 + m2*sigma^2*s2/2)/(z0 + m2*sigma^2*z2/2);
xiSym = xi *(z0 + m2*sigma^2*z2/2)/(h0 + m2*sigma^2*h2/2);
muSymSer = Normal[Series[muSym, {sigma, 0, 2}]];
xiSymSer = Normal[Series[xiSym, {sigma, 0, 2}]];
assertZero["M4a symmetric mu_eff to O(sigma^2)",
  muSymSer - (mu0*s0/z0 + (m2*sigma^2/2)*(mu0*s2/z0 - mu0*s0*z2/z0^2))];
assertZero["M4b symmetric xi_eff to O(sigma^2)",
  xiSymSer - (xi*z0/h0 + (m2*sigma^2/2)*(xi*z2/h0 - xi*z0*h2/h0^2))];

(* --- M5: Gaussian-localizer asymptotic fingerprints --- *)
ClearAll[w, sigma, ell, lam];
$Assumptions = sigma > 0 && ell > 0 && lam > 0;
Ws = Exp[-w^2/(2 sigma^2)]/(Sqrt[2 Pi] sigma);
Wm = Exp[-w/ell]/ell;
Zg = Exp[-w^2/lam^2];
ISym   = FullSimplify[Integrate[Ws*Zg, {w, -Infinity, Infinity}]];
ISymSer = Normal[Series[ISym, {sigma, 0, 4}]];
assertZero["M5a symmetric Gaussian asymptotic",
  ISymSer - (1 - sigma^2/lam^2 + 3*sigma^4/(2 lam^4))];
IMou  = FullSimplify[Integrate[Wm*Zg, {w, 0, Infinity}]];
(* derive series in ell from the evaluated integral via r = 1/ell at r → ∞ *)
ClearAll[r];
IMouR = FullSimplify[IMou /. ell -> 1/r];
IMouSer = FullSimplify[(Normal[Series[IMouR, {r, Infinity, 7}]]) /. r -> 1/ell];
assertZero["M5b mouth Gaussian asymptotic",
  IMouSer - (1 - 2*ell^2/lam^2 + 12*ell^4/lam^4 - 120*ell^6/lam^6)];

Print["STATUS: PASS"];
Exit[0];
```

The script must NOT structurally mirror the SymPy script's variable choreography. Use Mathematica's `Assuming`, `Series[…, {x, Infinity, n}]`, `Normal`, and `FullSimplify` rather than re-naming `sp.series(...).removeO()` patterns one-for-one. Re-derive the half-line integral and the mouth-Gaussian asymptotic with Mathematica's `Integrate` directly (do not retype the erfc closed form).

If Mathematica's `Integrate` yields an `Erfc` form for the mouth Gaussian, that is fine — what matters is that the form comes from `Integrate`, not from a literal.

**Claim manifest** (the new script must independently verify):
- **M1.** On the half-line `w ∈ [0, ∞)` with `W_ell(w) = exp(−w/ell)/ell`, the integration-by-parts identity holds:
  `[W Q]_0^∞ − ∫₀^∞ W'(w) Q(w) dw = ∫₀^∞ W(w) Q'(w) dw`,
  and the resulting expansions are `<Q>_ell = q₀ + ell·q₁ + ell²·q₂ + ell³·q₃ + ell⁴·q₄`, `<∂_w Q>_ell = q₁ + ell·q₂ + ell²·q₃ + ell³·q₄`.
- **M2.** For the unit-variance Gaussian kernel on `ℝ`, `<Q_series>_W = q₀ + q₂/2 + q₄/8` and `<∂_u Q_series>_W = q₁ + q₃/2 + q₅/8`.
- **M3.** The mouth-kernel zero-mode effective parameters expand as
  `μ_eff = μ₀ s₀/z₀ + ell·μ₀·(s₁/z₀ − s₀ z₁/z₀²) + O(ell²)`,
  `ξ_eff = ξ z₀/h₀ + ell·ξ·(z₁/h₀ − z₀ h₁/h₀²) + O(ell²)`.
- **M4.** The symmetric-kernel zero-mode effective parameters expand as
  `μ_eff = μ₀ s₀/z₀ + (m₂ σ²/2)·μ₀·(s₂/z₀ − s₀ z₂/z₀²) + O(σ⁴)`,
  `ξ_eff = ξ z₀/h₀ + (m₂ σ²/2)·ξ·(z₂/h₀ − z₀ h₂/h₀²) + O(σ⁴)`.
- **M5.** With `Z(w) = exp(−w²/λ²)`,
  - symmetric: `<Z>_σ = 1 − σ²/λ² + 3σ⁴/(2λ⁴) + O(σ⁶)`,
  - mouth: `<Z>_ell = 1 − 2 ell²/λ² + 12 ell⁴/λ⁴ − 120 ell⁶/λ⁶ + O(ell⁸)`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 009` and confirm the new check appears AND the script exits 0. The output must include `OK: ` lines for labels M1a–M1c, M2a–M2b, M3a–M3b, M4a–M4b, M5a–M5b, and a final `STATUS: PASS`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl`
- summary: Added an independent Mathematica audit script covering half-line IBP, Gaussian moments, effective-parameter series, and Gaussian-localizer asymptotics.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py:181–199`

**Issue:** Lines 188–195 declare unevaluated symbolic Integrals `IH_conc = sp.Integral(W_conc*H_sym, ...)` and `IS_conc = sp.Integral(W_conc*S_sym, ...)`, then substitute `H_sym → Z_conc` (resp. `S_sym → C*Z_conc`) and assert the result equals the evaluated `IZ_conc` (resp. `C*IZ_conc`). After the substitution, both sides are the same integrand integrated over the same interval, so the assertion is guaranteed by symbol substitution and integral linearity. Likewise the ratio assertions `xi*IZ_conc/IH_conc.subs(H_sym, Z_conc) = xi` and `mu0*IS_conc.subs(S_sym, C*Z_conc)/IZ_conc = C*mu0` reduce to `xi*IZ_conc/IZ_conc = xi` and `mu0*C*IZ_conc/IZ_conc = C*mu0`. None of the four checks can fail for any choice of `Z_conc`. They give a false signal of verification.

**Required change:**
Replace lines 181–199 (the block beginning `# exact cancellations for H=Z and S=C Z`, up to and including the closing `print("  if H = Z globally, then ξ_eff =", ...)`) with a perturbative check that has discriminating power.

Before (current lines 181–199):
```python
# exact cancellations for H=Z and S=C Z
C = sp.symbols("C")
mu_eff_proportional = sp.simplify(mu0 * (C*(z0 + ell*z1 + ell**2*z2)) / (z0 + ell*z1 + ell**2*z2))
xi_eff_exact_HZ = sp.simplify(xi * (z0 + ell*z1 + ell**2*z2) / (z0 + ell*z1 + ell**2*z2))
W_conc = sp.exp(-w / ell) / ell
Z_conc = z0 + z1*w + z2*w**2 / 2
IZ_conc = sp.integrate(W_conc * Z_conc, (w, 0, sp.oo))
H_sym = sp.Function("H")(w)
S_sym = sp.Function("S")(w)
IH_conc = sp.Integral(W_conc * H_sym, (w, 0, sp.oo))
IS_conc = sp.Integral(W_conc * S_sym, (w, 0, sp.oo))
assert_zero("H=Z before effective-gauge cancellation", IH_conc.subs(H_sym, Z_conc) - IZ_conc)
assert_zero("S=CZ before effective-coupling cancellation", IS_conc.subs(S_sym, C * Z_conc) - C * IZ_conc)
assert_zero("H=Z effective gauge", sp.simplify((xi * IZ_conc / IH_conc).subs(H_sym, Z_conc)) - xi)
assert_zero("S=CZ effective coupling", sp.simplify((mu0 * IS_conc / IZ_conc).subs(S_sym, C * Z_conc)) - C * mu0)
print()
print("Exact profile-alignment checks:")
print("  if S = C Z globally, then μ_eff =", sp.sstr(mu_eff_proportional))
print("  if H = Z globally, then ξ_eff =", sp.sstr(xi_eff_exact_HZ))
```

After:
```python
# perturbative profile-alignment: H = Z + eps*Δh, S = C*Z + eps*Δs.
C, eps = sp.symbols("C epsilon", real=True)
W_conc = sp.exp(-w / ell) / ell
Z_conc = z0 + z1*w + z2*w**2 / 2
H_pert = Z_conc + eps * (h0 + h1*w + h2*w**2/2)
S_pert = C * Z_conc + eps * (s0 + s1*w + s2*w**2/2)
IZ_conc = sp.integrate(W_conc * Z_conc, (w, 0, sp.oo))
IH_pert = sp.integrate(W_conc * H_pert, (w, 0, sp.oo))
IS_pert = sp.integrate(W_conc * S_pert, (w, 0, sp.oo))

# the cancellation must be exact at eps = 0 (genuine cancellation, not symbol substitution)
xi_eff_HZ_zero  = sp.simplify((xi  * IZ_conc / IH_pert).subs(eps, 0))
mu_eff_SCZ_zero = sp.simplify((mu0 * IS_pert / IZ_conc).subs(eps, 0))
assert_zero("H=Z effective gauge (eps=0 cancellation)",  xi_eff_HZ_zero - xi)
assert_zero("S=CZ effective coupling (eps=0 cancellation)", mu_eff_SCZ_zero - C * mu0)

# leading correction must equal -xi * <W Δh> / <W Z> (gauge) and mu0 * <W Δs> / <W Z> (coupling)
xi_eff_pert  = sp.series(sp.simplify(xi  * IZ_conc / IH_pert), eps, 0, 2).removeO()
mu_eff_pert  = sp.series(sp.simplify(mu0 * IS_pert / IZ_conc), eps, 0, 2).removeO()
IDh = sp.integrate(W_conc * (h0 + h1*w + h2*w**2/2), (w, 0, sp.oo))
IDs = sp.integrate(W_conc * (s0 + s1*w + s2*w**2/2), (w, 0, sp.oo))
assert_zero("xi_eff first-order correction in eps",
            xi_eff_pert - (xi - eps * xi * IDh / IZ_conc))
assert_zero("mu_eff first-order correction in eps",
            mu_eff_pert - (C * mu0 + eps * mu0 * IDs / IZ_conc))

print()
print("Perturbative profile-alignment checks (H = Z + eps Δh, S = C Z + eps Δs):")
print("  ξ_eff |_{eps=0} =", sp.sstr(xi_eff_HZ_zero))
print("  μ_eff |_{eps=0} =", sp.sstr(mu_eff_SCZ_zero))
print("  ξ_eff to O(eps) =", sp.sstr(xi_eff_pert))
print("  μ_eff to O(eps) =", sp.sstr(mu_eff_pert))
```

Notes for the implementer:
- Keep the `C = sp.symbols("C")` declaration but add `eps` to it as shown.
- Do NOT remove `IZ_conc` — it is used both above and (after this directive's F3 / F4 patches) below.
- The new `assert_zero` calls use the existing `assert_zero` helper unchanged.
- If `sp.series` returns the leading-order in a form that `assert_zero` cannot reduce to 0, replace `xi_eff_pert - (xi - eps * xi * IDh / IZ_conc)` with `sp.expand(xi_eff_pert - (xi - eps * xi * IDh / IZ_conc))` before `assert_zero`.

**Verification:**
Verifier reruns SymPy. The output must contain the four labels "H=Z effective gauge (eps=0 cancellation)", "S=CZ effective coupling (eps=0 cancellation)", "xi_eff first-order correction in eps", "mu_eff first-order correction in eps". The labels "H=Z before effective-gauge cancellation", "S=CZ before effective-coupling cancellation", "H=Z effective gauge", "S=CZ effective coupling" must be gone. Exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py`
- summary: Replaced tautological profile-alignment checks with perturbative epsilon checks for exact cancellation and first-order corrections.
- deviation: none

## F3 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py:110, 149–156, 207–208`

**Issue:** Six quantitative results are printed but never asserted:
- `avg_Q_half` (line 110) — the `<Q>_ell` expansion.
- `mu_eff_sym_series` (line 149) — symmetric μ_eff to O(σ²).
- `xi_eff_sym_series` (line 150) — symmetric ξ_eff to O(σ²).
- `mu_eff_half_series` (line 155) — mouth μ_eff to O(ell).
- `xi_eff_half_series` (line 156) — mouth ξ_eff to O(ell).
- `IWZ_sym_gauss_series` (line 208) — symmetric Gaussian asymptotic.

These series carry the section's headline claims (output text lines 46, 87–92, 116–117). Each must be anchored to a literal target expansion.

**Required change:**
Add six new `assert_zero(...)` calls. Each must appear immediately after the assignment that produces the series (so the verifier can confirm the new lines are present).

1. After line 110 (i.e. as a new line between the existing `avg_Q_half = …` and `avg_dQ_half = …`), insert:
```python
assert_zero("half-line Q expansion", avg_Q_half - (q0 + ell*q1 + ell**2*q2 + ell**3*q3 + ell**4*q4))
```

2. After line 149 (`mu_eff_sym_series = sp.series(...)`), insert:
```python
assert_zero("symmetric mu_eff series",
            mu_eff_sym_series - (mu0*s0/z0 + (m2*sigma**2/2)*(mu0*s2/z0 - mu0*s0*z2/z0**2)))
```

3. After line 150 (`xi_eff_sym_series = sp.series(...)`), insert:
```python
assert_zero("symmetric xi_eff series",
            xi_eff_sym_series - (xi*z0/h0 + (m2*sigma**2/2)*(xi*z2/h0 - xi*z0*h2/h0**2)))
```

4. After line 155 (`mu_eff_half_series = sp.series(...)`), insert:
```python
assert_zero("mouth mu_eff series",
            mu_eff_half_series - (mu0*s0/z0 + ell*(mu0*s1/z0 - mu0*s0*z1/z0**2)))
```

5. After line 156 (`xi_eff_half_series = sp.series(...)`), insert:
```python
assert_zero("mouth xi_eff series",
            xi_eff_half_series - (xi*z0/h0 + ell*(xi*z1/h0 - xi*z0*h1/h0**2)))
```

6. After line 208 (`IWZ_sym_gauss_series = sp.series(...)`), insert:
```python
assert_zero("symmetric Gaussian asymptotic literal",
            IWZ_sym_gauss_series - (1 - sigma**2/lam**2 + 3*sigma**4/(2*lam**4)))
```

If `assert_zero` cannot collapse a residual to 0 because of `sp.series` ordering, wrap the difference in `sp.expand(...)` before passing it to `assert_zero`.

**Verification:**
Verifier reruns SymPy. The output must contain six new label lines: "half-line Q expansion", "symmetric mu_eff series", "symmetric xi_eff series", "mouth mu_eff series", "mouth xi_eff series", "symmetric Gaussian asymptotic literal". Exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py`
- summary: Added literal residual assertions for the half-line Q expansion, effective-parameter series, and symmetric Gaussian asymptotic.
- deviation: none

## F4 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py:213–217`

**Issue:** The mouth-Gaussian asymptotic series is derived from a hand-typed `erfc` closed form `IWZ_mouth_gauss_r = sqrt(pi)*lam*r*erfc(lam*r/2)*exp(lam²r²/4)/2` (line 214) rather than from the SymPy-evaluated integral `IWZ_mouth_gauss` (line 212). The SymPy integral is computed and printed but never used in any assertion. A typo in line 214 would propagate through `IWZ_mouth_gauss_series` undetected by A9 (which compares two literals) and would only be caught by A10 (Taylor-then-integrate) if the typo also broke that comparison. The hand-typed form is effectively a pre-baked answer; the script's structure should derive it.

**Required change:**
Replace lines 213–217 to derive the asymptotic series from `IWZ_mouth_gauss` itself, and add an explicit check that `IWZ_mouth_gauss` agrees with the erfc closed form.

Before (lines 213–217):
```python
r = sp.symbols("r", positive=True)
IWZ_mouth_gauss_r = sp.sqrt(sp.pi) * lam * r * sp.erfc(lam * r / 2) * sp.exp(lam**2 * r**2 / 4) / 2
IWZ_mouth_gauss_series = sp.simplify(
    sp.series(IWZ_mouth_gauss_r, r, sp.oo, 8).removeO().subs(r, 1 / ell)
)
```

After:
```python
r = sp.symbols("r", positive=True)
# derive the asymptotic series directly from the SymPy-evaluated integral via ell = 1/r at r → ∞
IWZ_mouth_gauss_r = sp.simplify(IWZ_mouth_gauss.subs(ell, 1/r))
IWZ_mouth_gauss_series = sp.simplify(
    sp.series(IWZ_mouth_gauss_r, r, sp.oo, 8).removeO().subs(r, 1 / ell)
)
# guard: the SymPy integral really equals the erfc closed form
IWZ_mouth_gauss_erfc = sp.sqrt(sp.pi) * lam * sp.erfc(lam / (2*ell)) * sp.exp(lam**2 / (4*ell**2)) / (2*ell)
assert_zero("mouth Gaussian integral equals erfc closed form",
            IWZ_mouth_gauss - IWZ_mouth_gauss_erfc)
```

Notes for the implementer:
- Keep lines 218–223 (the Taylor-integration cross-check and the existing A9, A10) unchanged. They still apply.
- If `sp.simplify(IWZ_mouth_gauss - IWZ_mouth_gauss_erfc)` does not reduce to 0 because SymPy returns the integral in a different equivalent form (e.g., with `erf` instead of `erfc`), replace the guard's body with `sp.simplify(sp.expand(IWZ_mouth_gauss - IWZ_mouth_gauss_erfc))`. If that still fails, use `sp.simplify((IWZ_mouth_gauss - IWZ_mouth_gauss_erfc).rewrite(sp.erfc))`. Do not skip the guard.

**Verification:**
Verifier reruns SymPy. The output must contain the new label "mouth Gaussian integral equals erfc closed form" before the existing two mouth-Gaussian asymptotic checks. The two existing checks on lines 222–223 ("mouth Gaussian asymptotic from erfc closed form" and "mouth Gaussian asymptotic from Taylor integration") must still pass. Exit 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py`
- summary: Derived the mouth-Gaussian asymptotic from the SymPy-evaluated integral and added an erfc-form equality guard.
- deviation: Rewrote the evaluated integral to `erfc` before the infinity-series expansion because SymPy cannot expand its raw `erf` form at infinity.
