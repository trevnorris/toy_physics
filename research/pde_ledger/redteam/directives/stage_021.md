---
unit_id: 021
batch: I.2
created_at: 2026-05-21T12:30:55-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-21T15:16:05-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 021

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:73-81`

**Issue:**
The Mathematica Section II.1 currently contains only one assertion touching the reduced Lagrangian `lRed`, namely `expectZero["Q kinetic coefficient", Coefficient[lVel, vq^2] - m/2]` (line 81). Because `lRed` is constructed at line 67 with the explicit term `1/2 m D[q,t]^2`, that coefficient is `m/2` by construction, making the assertion algebraically guaranteed independent of any physics. The SymPy counterpart (sympy file lines 123-130) actually derives the three Euler-Lagrange equations from `Lred` via `euler_equations(...)` and matches each one against the expected EOM. The Mathematica script never derives the EOMs from `lRed`, so any sign or coupling error in the Lagrangian definition would not be detected on the Mathematica side.

**Required change:**

Edit the file as follows. The current lines 73-81 are:

```
staticL = lRed /. {D[q, t] -> 0, D[a, t] -> 0, D[ww, t] -> 0};
staticTmp = staticL /. {q -> q0, a -> a0r, ww -> w0};
staticBack = {q0 -> q, a0r -> a, w0 -> ww};
qd = D[q, t];
ad = D[a, t];
wd = D[ww, t];
lVel = Expand[lRed /. {qd -> vq, ad -> va, wd -> vw}];

expectZero["Q kinetic coefficient", Coefficient[lVel, vq^2] - m/2];
```

Replace those nine lines (73-81 inclusive) with the following six lines, which compute the Euler-Lagrange operator for each mode directly and verify against the canonical EOMs:

```
elQ = D[D[lRed, D[q, t]], t] - D[lRed, q];
elA = D[D[lRed, D[a, t]], t] - D[lRed, a];
elW = D[D[lRed, D[ww, t]], t] - D[lRed, ww];
expectZero["Q equation", elQ - (m D[q, t, t] + k q - gA a - gW ww)];
expectZero["A equation", elA - (D[a, t, t] + oA^2 a - r ww - gA q)];
expectZero["W equation", elW - (D[ww, t, t] + oW^2 ww - r a - gW q)];
```

Do not change any other lines. The blank line at line 80 (between the assignments block and the `expectZero` call) collapses naturally; preserve a single blank line of separation if the local style elsewhere in the file uses one.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 021` and confirm:
- the script still exits 0;
- the output file no longer contains the line `Q kinetic coefficient = ...`;
- the output file contains three new lines: `PASS: Q equation`, `PASS: A equation`, `PASS: W equation`;
- the residuals printed before those PASS lines are all `0`.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`
- summary: Replaced the tautological reduced-Lagrangian kinetic coefficient check with VariationalMethods EulerEquations checks for Q, A, and W.
- deviation: Also parenthesized the existing lRed RHS so Mathematica parses the intended full reduced Lagrangian instead of only the first line.

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:88-89, 131-132, 149-151`

**Issue:**
The Mathematica script reproduces the SymPy script's algebra step-for-step using the same intermediate-variable choreography (aKer ↔ Aker, wKer ↔ Wker, delta ↔ Delta, sigmaCons ↔ Sigma_cons, sigmaFull ↔ Sigma_full, sigmaFirst ↔ Sigma_first, nOmega ↔ N_omega, j2a ↔ j2a, y2a ↔ y2a, h2a ↔ h2a, lambda2 ↔ Lambda2). The "second engine" check thereby reduces to "the same algebra simplifies to the same form in two CAS surfaces". Three Mathematica-native independent derivations are available and replace the transliterated steps without altering the closed-form targets being verified.

**Required change:**

Three separate, localized edits in `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`. Apply all three; the existing assertions downstream of each block must continue to work unchanged.

**Edit 2.1 — Section II.2 Schur complement via LinearSolve.** The current lines 88-89 are:

```
aSol = FullSimplify[(gA wKer + gW r)/delta, Assumptions -> $Assumptions];
wSol = FullSimplify[(gW aKer + gA r)/delta, Assumptions -> $Assumptions];
```

Replace those two lines with the following six lines, which derive `aSol` and `wSol` by inverting the 2x2 mixed kernel against the source vector `{gA, gW}`, and add one new assertion that the resulting `sigmaCons` from `gA*aSol + gW*wSol` matches the closed form already assigned to `sigmaCons` at line 86:

```
matEAW = {{aKer, -r}, {-r, wKer}};
solAW = LinearSolve[matEAW, {gA, gW}];
aSol = FullSimplify[solAW[[1]], Assumptions -> $Assumptions];
wSol = FullSimplify[solAW[[2]], Assumptions -> $Assumptions];
sigmaConsDerived = FullSimplify[gA aSol + gW wSol, Assumptions -> $Assumptions];
expectZero["sigmaCons from LinearSolve matches closed form", sigmaConsDerived - sigmaCons];
```

Place the new lines in the same position the deleted lines occupied (immediately after the line `wSol = FullSimplify[...]` would have been, i.e. just before the `Print["Sigma_EM+mix^cons(omega) = ..."]` line at the current line 91).

**Edit 2.2 — Section III analytic derivative for N(omega).** The current lines 131-132 are:

```
sigmaFirst = Expand[Normal[Series[sigmaFull, {piOut, 0, 1}]]];
nOmega = FullSimplify[(sigmaFirst - sigmaCons)/piOut, Assumptions -> $Assumptions];
```

Replace those two lines with a single analytic-derivative line:

```
nOmega = FullSimplify[D[sigmaFull, piOut] /. piOut -> 0, Assumptions -> $Assumptions];
```

Delete the now-unused `sigmaFirst` assignment entirely. Keep the surrounding lines (the `sigmaFull` definition above and the `n0 = FullSimplify[nOmega /. omega -> 0, ...]` line below) unchanged.

**Edit 2.3 — Section IV spherical Hankel via built-in.** The current lines 149-151 are:

```
j2a = ((3/za^3) - 1/za) Sin[za] - 3 Cos[za]/za^2;
y2a = -((3/za^3) - 1/za) Cos[za] - 3 Sin[za]/za^2;
h2a = FullSimplify[j2a + I y2a, Assumptions -> $Assumptions];
```

Replace those three lines with a single built-in call:

```
h2a = SphericalHankelH1[2, za];
```

`SphericalHankelH1[2, za]` is the standard spherical Hankel function of the first kind at `l=2` and equals `j_2(za) + I y_2(za)` exactly; the small-`za` expansion downstream is unaffected. Keep all subsequent lines unchanged.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 021` and confirm:
- the script still exits 0;
- the output file contains a new line `PASS: sigmaCons from LinearSolve matches closed form`;
- the existing PASS lines `A exact solution residual`, `W exact solution residual`, `N(omega) compact formula`, `N(0) positive-square form`, `Y2_hat minimal branch`, `Gamma5_port - a^5/(27 c_s^5)`, `N_scalar leading term`, `scalar odd order` all remain present;
- the printed values of `N(omega)`, `N(0)`, `delta D_wall^(odd)`, `Lambda2(k)`, `Y2_hat(omega)`, `Gamma5_port` are unchanged (the closed forms must match the previous output line-for-line modulo Mathematica's `Together`/`FullSimplify` reordering).

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`
- summary: Replaced transliterated Schur complement, first-order expansion, and spherical Hankel steps with LinearSolve, analytic differentiation, and SphericalHankelH1.
- deviation: none
