---
unit_id: 021
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T15:25:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 021 (iteration 2)

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed (iter 2):**
Codex applied the delta directive's preferred option (a) — built-in
`EulerEquations` from the `VariationalMethods`` package — and additionally
parenthesized the `lRed` multi-line assignment (a second defect Codex caught
on its own and recorded under `## Applied: F1` as a deviation).

Concrete edits in
`mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`:

- Line 33 (new): `Needs["VariationalMethods``"];` inserted in the global
  preamble before the first `banner[...]` call.
- Lines 69-74: `lRed` RHS is now wrapped in parentheses:
  ```
  lRed = (
    1/2 m D[q, t]^2 - 1/2 k q^2
    + 1/2 D[a, t]^2 - 1/2 oA^2 a^2
    + 1/2 D[ww, t]^2 - 1/2 oW^2 ww^2
    + r a ww + gA q a + gW q ww
  );
  ```
  The pre-iter-2 form omitted the parens and ended with `+ gW q ww;` on its
  own line; Mathematica's `;` line continuation in that layout would only
  bind the first line of the RHS to `lRed` and silently discard the rest
  (the same defect we saw in batch I.1 stage 003). The parentheses force
  full-RHS binding.
- Lines 76-79: the manual EL operator from iter 1 is replaced by
  ```
  elList = EulerEquations[lRed, {qFun[t], aFun[t], wFun[t]}, t];
  expectZero["Q equation", (elList[[1]] /. Equal[lhs_, rhs_] :> lhs - rhs) + (m qFun''[t] + k qFun[t] - gA aFun[t] - gW wFun[t])];
  expectZero["A equation", (elList[[2]] /. Equal[lhs_, rhs_] :> lhs - rhs) + (aFun''[t] + oA^2 aFun[t] - r wFun[t] - gA qFun[t])];
  expectZero["W equation", (elList[[3]] /. Equal[lhs_, rhs_] :> lhs - rhs) + (wFun''[t] + oW^2 wFun[t] - r aFun[t] - gW qFun[t])];
  ```
  Note: Codex used `+` (not `-`) before the canonical-RHS expression. This
  is correct because `EulerEquations` returns the d'Alembert form
  `D[L,q] - D[D[L,q'],t] == 0`, which for our Lagrangian evaluates to
  `-(m q'' + k q - gA a - gW w)`; adding `+(m q'' + k q - gA a - gW w)`
  cancels to 0.

**Assessment:**
The fresh output
`mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt`
(mtime 2026-05-21 15:16:05, newer than .wl mtime 15:14:57) reports:

```
Q equation = 0
PASS: Q equation
A equation = 0
PASS: A equation
W equation = 0
PASS: W equation
```

These are genuine zeros printed by `expectZero` (which `FullSimplify`s the
expression and prints `fmt[res]` before declaring PASS), not just bare PASS
labels. The residuals are non-tautologically zero: each subtracts the
EulerEquations-derived EOM from the independently written canonical EOM,
and the simplification collapses to 0. The exit code is 0.

The parenthesization fix is independently confirmed by these PASS results:
if `lRed` had been bound to only its first line (just the q kinetic and
mass terms), `EulerEquations` would produce `m q'' + k q` for Q (no
couplings) and trivial equations for A and W, and the residuals against the
full canonical EOMs (which include `-gA a - gW ww` couplings) would
absolutely not be zero. The fact that all three residuals are exactly 0
proves lRed now binds the entire RHS as intended.

The Codex deviation (parenthesizing lRed) was unsolicited but correct and
necessary; without it, iter 2 would have failed at the Q assertion with a
different residual than iter 1. Recording it under `## Applied: F1` with
`deviation: <text>` is the right disclosure pattern.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
All three F2 edits from iter 1 remain in place and now actually execute
(iter 1's F1 failure had blocked them from running):

- Section II.2 (lines 86-91): the 2x2 Schur complement via `LinearSolve`.
  Produces `aSol = solAW[[1]]`, `wSol = solAW[[2]]`, then verifies
  `sigmaConsDerived = gA aSol + gW wSol` matches the closed-form
  `sigmaCons` via `expectZero["sigmaCons from LinearSolve matches closed form", sigmaConsDerived - sigmaCons]`.
- Section III (line 133): `nOmega = FullSimplify[D[sigmaFull, piOut] /. piOut -> 0, ...]`
  — direct analytic differentiation, no Series.
- Section IV (line 150): `h2a = SphericalHankelH1[2, za]` — built-in
  spherical Hankel function of the first kind at l=2.

**Assessment:**
Output confirms all three F2 edits execute and pass:

- Line 31-32 of output: `sigmaCons from LinearSolve matches closed form = 0; PASS: sigmaCons from LinearSolve matches closed form`. The
  printed residual is 0, meaning the matrix-inverse route algebraically
  matches the hand-written closed form `sigmaCons = (gA^2 wKer + 2 gA gW r + gW^2 aKer)/delta`.
- Lines 62-63: `N(omega) compact formula = 0; PASS: N(omega) compact formula`.
  This means the analytic-derivative-based `nOmega` collapses to
  `(aKer gW + r gA)^2/delta^2`, the closed form previously obtained from
  the first-order series coefficient. Independent route, same answer.
- Lines 73-74: `Y2_hat minimal branch = 0; PASS: Y2_hat minimal branch`,
  and lines 75-76: `Gamma5_port - a^5/(27 c_s^5) = 0; PASS`. The
  `SphericalHankelH1[2, za]` substitution produces the same `Lambda2(k)`
  expansion, `Y2_hat(omega)` rational expansion, and `Gamma5_port = a^5/(27 c_s^5)`.

Comparison to the iter-1 verifier's quoted "previously published"
Section II.2/III/IV/V values:

- `Sigma_EM+mix^cons(omega) = (gW^2*(-oA^2 + omega^2) + gA^2*(omega - oW)*(omega + oW) - 2*gA*gW*r)/((oA - omega)*(oA + omega)*(omega - oW)*(omega + oW) + r^2)` — matches (line 33 of output).
- `D_cons(omega) = k - m*omega^2 + (gA^2*(-omega^2 + oW^2) + gW*(gW*(oA - omega)*(oA + omega) + 2*gA*r))/((oA - omega)*(oA + omega)*(omega - oW)*(omega + oW) + r^2)` — matches (line 34).
- z0/z2/z4 EM+mix expressions — match (lines 45-47).
- `Sigma_full(omega) = (gW^2*(-oA^2 + omega^2) + gA^2*(omega^2 - oW^2 + piOut) - 2*gA*gW*r)/((oA - omega)*(oA + omega)*(omega^2 - oW^2 + piOut) + r^2)` — matches (line 58).
- `N(omega) = (gW*(oA - omega)*(oA + omega) + gA*r)^2/((oA - omega)*(oA + omega)*(omega - oW)*(omega + oW) + r^2)^2` — matches (line 59).
- `N(0) = (gW*oA^2 + gA*r)^2/(-(oA^2*oW^2) + r^2)^2` — matches (line 60).
- `delta D_wall^(odd) = ((-I)*gammaPort*omega^5*(gW*oA^2 + gA*r)^2)/(-(oA^2*oW^2) + r^2)^2` — matches (line 61).
- `Lambda2(k) = -3/radius + (kWave^2*radius)/3 + (kWave^4*radius^3)/9 + (I/9)*kWave^5*radius^4 - (2*kWave^6*radius^5)/27` — matches (line 70).
- `Y2_hat(omega) = 1 + (omega^2*radius^2)/(9*cS^2) + (4*omega^4*radius^4)/(81*cS^4) + ((I/27)*omega^5*radius^5)/cS^5` — matches (line 71).
- `Gamma5_port = radius^5/(27*cS^5)` — matches (line 72).
- `N_scalar(omega) = (eta^2*oA^4*omega^2)/(-(oA^2*oW^2) + r^2)^2` — matches (line 81).
- `Pi0_out * N_scalar = (I*eta^2*gamma1*oA^4*omega^3)/(-(oA^2*oW^2) + r^2)^2` — matches (line 82).

All Section II.2/III/IV/V values are byte-identical to the pre-iter-1
output the verifier noted as the published baseline. Downstream numbers
are unchanged.

## Exec log assessment

**SymPy:** exit=0. No SymPy edits were applied this iteration. The log
shows the full FINAL STAGE-4 LEDGER closing with the expected verified
quantities and `# exit_code: 0`. Notable lines:

```
scalar odd order = 0
  • the compact outgoing l=2 branch has the normalized fingerprint
    Y2_hat = 1 + a^2 omega^2/(9 c_s^2) + 4 a^4 omega^4/(81 c_s^4) + i a^5 omega^5/(27 c_s^5) + ...;
# exit_code: 0
```

**Mathematica:** exit=0. No dedicated `stage_021_mathematica.log` was
captured in `redteam/exec_logs/` (only `stage_021_diff.patch` and
`stage_021_sympy.log` present). The verifier read the saved
`mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt`
as instructed (mtime 2026-05-21 15:16:05, fresh — newer than .wl mtime
15:14:57). Header confirms `# Exit code: 0`, `# Status: PASS`, and the
trailing `EXIT_CODE: 0`. Notable PASS lines (all 19 of them):

```
PASS: E_w gauge variation
PASS: C_a gauge variation
PASS: Q equation                                  (new — iter 2 F1)
PASS: A equation                                  (new — iter 2 F1)
PASS: W equation                                  (new — iter 2 F1)
PASS: sigmaCons from LinearSolve matches closed form   (new — F2.1)
PASS: A exact solution residual
PASS: W exact solution residual
PASS: z0 formula
PASS: z2 formula
PASS: z4 formula
PASS: Sigma z0
PASS: Sigma z2
PASS: Sigma z4
PASS: N(omega) compact formula                    (now via D[]/piOut — F2.2)
PASS: N(0) positive-square form
PASS: Y2_hat minimal branch                       (now via SphericalHankelH1 — F2.3)
PASS: Gamma5_port - a^5/(27 c_s^5)
PASS: N_scalar leading term
PASS: scalar odd order
```

The "Q kinetic coefficient" line is correctly absent (F1 deletion target).

**Output freshness:** confirmed. .wl mtime 2026-05-21 15:14:57; saved .txt
mtime 2026-05-21 15:16:05; iter-1 .txt would have been written 15:02:41 so
this is genuinely the iter-2 re-run.

## Material-change assessment

`material_change`: false.

The published closed forms for `Sigma_EM+mix^cons`, `D_cons`, the z0/z2/z4
low-frequency coefficients, `Sigma_full`, `N(omega)`, `N(0)`, `delta
D_wall^(odd)`, `Lambda2(k)`, `Y2_hat(omega)`, `Gamma5_port`,
`N_scalar(omega)`, and `Pi0_out * N_scalar` are byte-identical to the
pre-edit output (verified line by line above). All downstream units (022+)
that consume these quantities are unaffected. The iter-2 changes are
purely verification-method changes (EulerEquations instead of a
tautological kinetic-coefficient check; LinearSolve instead of a
hand-written matrix inverse; analytic derivative instead of Series
truncation; built-in SphericalHankelH1 instead of a hand-built combination
of `j_2` and `y_2`), with the same algebraic targets.

## Side observations (non-blocking)

- No dedicated `redteam/exec_logs/stage_021_mathematica.log` is present
  this iteration either; only the .txt in `mathematica/output/`. Same
  observation as iter 1 — the verifier scope permits reading the .txt and
  the freshness check passes, so this does not block verification, but the
  orchestrator may want to know.
- The Codex `## Applied: F1` block records the parenthesization fix as a
  `deviation`. This is the correct disclosure pattern: the directive did
  not request it, but it was necessary for the requested F1 fix to actually
  exercise the full Lagrangian. Without it, EulerEquations would have been
  computed against a truncated lRed and produced different residuals.
- The shadowing of `gA`, `gW` in Section V (line 174-175 setting `gA = 0`
  and `gW = eta omega`) is intentional and matches the SymPy file, as noted
  in iter 1. Not new.
- The sign convention in F1's EL residual check (using `+ (m qFun''[t] + ...)`
  rather than `- (m qFun''[t] + ...)`) is the correct adaptation for
  `EulerEquations` (d'Alembert form). The fact that all three Q/A/W
  residuals collapse to exactly 0 confirms this is the right sign; if it
  were wrong by a factor of -1, the residuals would be 2x the EOM rather
  than 0.

## Verdict justification

Both F1 and F2 are now fully resolved. F1's `EulerEquations` adoption plus
the (unsolicited but correct) `lRed` parenthesization fix means the three
EL checks compare the EulerEquations-derived EOMs against independently
written canonical EOMs and obtain residuals that `FullSimplify` reduces to
0. F2's LinearSolve / analytic-derivative / SphericalHankelH1 substitutions
all execute and pass, with their downstream `Sigma_EM+mix^cons`, `N(omega)`,
`N(0)`, `Lambda2(k)`, `Y2_hat(omega)`, `Gamma5_port`, `N_scalar` values
byte-identical to the pre-edit baseline. Mathematica exits 0, SymPy exits
0, no regressions in the diff, all 19 PASS lines present in the fresh
output. Verified.
