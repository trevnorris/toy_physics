---
unit_id: 148
batch: IV.5
created_at: 2026-05-27T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 148

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.
After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — insufficient_verification (load-bearing D4 check is tautological / un-asserted)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py:87-89`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl:83-85`

**Issue:**
The Stage 228 broadening fraction `ξ_* = (-37√3 - 5π² + 2√(4107 - 168π²))/(5(8 - π²))` quoted in the notes (notes line 132) is the load-bearing comparison target for `1 - λ_{Π,0}`. The Mathematica script computes the right-hand side from the same `gMinus` already used to construct `λ_{Π,0}`, making the `expectZero` check tautological (the two sides are algebraically identical by construction — see report self-test). The SymPy script uses a hardcoded 15-digit float and only `print`s the residual without `assert`.

Replace both `ξ_*` constructions with the closed-form symbolic expression from the notes so the consistency check exercises a genuine algebraic identity between (a) the linear-response neutral point and (b) the independent Stage 228 closed form.

**Required change:**

SymPy — replace lines 86-89 (current text):
```python
# Stage 127 consistency value xi_* (exact positive-family compensation broadening)
xi_star = sp.N(sp.Float("0.183918405511538"), 30)
print("xi_* (Stage 127) =", xi_star)
print("(1-lambda_(Pi,0)) - xi_* =", sp.N((1-lam_Pi_zero) - xi_star, 20))
```
with:
```python
# Stage 228 positive-family compensation closed form (see notes section 3)
xi_star_closed = (-37*sp.sqrt(3) - 5*sp.pi**2 + 2*sp.sqrt(4107 - 168*sp.pi**2)) / (5*(8 - sp.pi**2))
xi_star = sp.N(xi_star_closed, 30)
print("xi_* (Stage 228 closed form) =", xi_star)
residual = sp.N((1 - lam_Pi_zero) - xi_star, 30)
print("(1-lambda_(Pi,0)) - xi_* =", residual)
assert abs(residual) < sp.Float("1e-25"), f"Stage 148 D4 consistency failed: residual = {residual}"
```

Mathematica — replace lines 83-85 (current text):
```mathematica
xiStar = FullSimplify[((Pi/4) - gMinus)/((Pi/4) - 2/Pi), Assumptions -> True];
Print["xi_* (Stage 126) = ", fmt[N[xiStar, 30]]];
expectZero["(1-lambda_(Pi,0)) - xi_*", (1 - lamPiZero) - xiStar];
```
with:
```mathematica
xiStarClosed = (-37*Sqrt[3] - 5*Pi^2 + 2*Sqrt[4107 - 168*Pi^2]) / (5*(8 - Pi^2));
Print["xi_* (Stage 228 closed form) = ", fmt[N[xiStarClosed, 30]]];
expectZero["(1-lambda_(Pi,0)) - xi_*", FullSimplify[(1 - lamPiZero) - xiStarClosed]];
```

The `FullSimplify` in the new Mathematica `expectZero` must reduce the difference of two algebraically distinct closed forms to exact zero. If it does not, the assertion will fail and reveal that the notes' closed form does not actually match the linear-response neutral point — which is precisely the bridge the paper claims.

**Claim manifest:**

M1: `1 - λ_{Π,0} = ξ_*` where `λ_{Π,0}` is the convex-interpolation bias-neutral point computed from the canonical lower-compensated branch, and `ξ_* = (-37·√3 - 5·π² + 2·√(4107 - 168·π²))/(5·(8 - π²))` is the Stage 228 closed-form positive-family compensation broadening fraction.

**Verification command:**
`redteam exec-sympy 148` and `redteam exec-mathematica 148` — both must exit 0 with the new assertion present in the output transcripts.

## F2 — hardcoded_result (sympy `xi_star`)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py:87`

**Issue:**
The literal `sp.Float("0.183918405511538")` is a hardcoded 15-digit numeric constant standing in for the Stage 228 closed form. The notes give the closed form explicitly.

**Required change:**
Subsumed by F1. After applying F1, F2 is resolved automatically (the literal is replaced by the symbolic closed form `xi_star_closed`). No additional edits required; record `## Applied: F2` as `files_changed: (subsumed by F1)`, `summary: subsumed by F1 fix`, `deviation: none`.

**Verification command:** Same as F1.

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl:26, 43-72`

**Issue:**
The Mathematica script duplicates the SymPy script's algebraic choreography line-for-line (`aT`, `bT` coefficients are bit-for-bit ports of the sympy `AT`, `BT` definitions; the convex-interpolation block reuses the same `(gLam - gStar)/gPrimeStar` form). The banner string on line 26 still reads `STAGE 131`, indicating the file was templated from another stage and the algebra was ported wholesale rather than re-derived.

**Required change:**

(a) Banner fix — line 26, replace
```mathematica
banner["STAGE 131 — REPRESENTATIVE NON-EXPONENTIAL POSITIVE FAMILIES"];
```
with
```mathematica
banner["STAGE 148 — REPRESENTATIVE NON-EXPONENTIAL POSITIVE FAMILIES"];
```

(b) Restructure the `dT` linear-response derivation in Mathematica to use a symbolically distinct path. Instead of pre-collecting `aT`, `bT` and then writing `dTU = aT*(gU - gStar) + bT*(sU - sStar)` etc. (the sympy form), express `dT` via the intermediate `dSigma`:

Replace lines 43-47:
```mathematica
aT = N[
  -(9/(40*tStar))*(1/(gPrimeStar*(1 - sStar/4)) + pStar*sPrimeStar/(4*gPrimeStar*(1 - sStar/4)^2)),
  30
];
bT = N[(9/(40*tStar))*pStar/(4*(1 - sStar/4)^2), 30];
```
with the following alternative formulation (derives `dT` via `dSigma`, matching the symbolic structure `Σ = Π/(1 - S/4)` so `dΣ = dΠ/(1 - S/4) + Π·dS/(4·(1 - S/4)²)`, and `T = √(9·Σ/20)` so `dT = (9/(40·T))·dΣ`):
```mathematica
(* dSigma in terms of dPi (= -(g_target - gStar)/gPrimeStar) and dS = (s_target - sStar) *)
dSigmaOfDeltas[dPi_, dS_] := dPi/(1 - sStar/4) + pStar*dS/(4*(1 - sStar/4)^2);
dTOfDeltas[dG_, dS_] := Module[{dPi},
  dPi = -dG/gPrimeStar;
  N[(9/(40*tStar))*dSigmaOfDeltas[dPi, dS], 30]
];
```

Then replace each downstream use accordingly:
- Line 52: `dTU = aT*(gU - gStar) + bT*(sU - sStar);` → `dTU = dTOfDeltas[gU - gStar, sU - sStar];`
- Line 62: `dTD = aT*(gD - gStar) + bT*(sD - sStar);` → `dTD = dTOfDeltas[gD - gStar, sD - sStar];`
- Line 72: `dTLam = FullSimplify[aT*(gLam - gStar) + bT*(sLam - sStar)];` → `dTLam = FullSimplify[dTOfDeltas[gLam - gStar, sLam - sStar]];`

This routes the Mathematica algebra through the intermediate `dSigma` step (the algebraically natural derivation) while sympy continues to use the pre-collected `AT`, `BT` form. The two paths converge symbolically only via the algebraic identity
`aT*dG + bT*dS == (9/(40*tStar))*(-(dG)/(gPrimeStar*(1 - sStar/4)) + pStar*dS/(4*(1 - sStar/4)^2))`
which the engines now establish independently. Remove (or comment out) the now-unused `aT`/`bT` definitions, or keep them as a `Print`-only sanity check.

**Verification command:**
`redteam exec-mathematica 148` — output must show `STAGE 148` in the banner, and the numerical values of `dTU`, `dTD`, `dTLam` must agree with the SymPy output to the printed precision (28-30 digits). If `dTLam` no longer factors to the same `0.5087... - 0.6257...·lam` literal form, the new derivation has diverged; investigate before declaring done.

## Self-test notes (auditor)

- For F1's new assertion, the necessary identity is
  `((π/4) - gMinus)/((π/4) - 2/π) == (-37√3 - 5π² + 2√(4107 - 168π²))/(5(8 - π²))`
  with `gMinus = rF1 - √(1+rF1²)/2`, `rF1 = √(12·(37/20)²/π² - 1)`.
  Both engines must reduce the difference to exact 0 under `FullSimplify` / `simplify`; the notes assert this agreement to 15 digits numerically (line 132).
- For F3, after the restructure, `dTLam` should still expand to `0.508756302... - 0.625700104...·lam` — i.e., the *numerics* are preserved while the *symbolic path* differs.
- Codex must not edit the notes or paper card; the Stage 250 / Stage 127 / Stage 126 / STAGE 131 label inconsistencies in those documents are out-of-scope.
