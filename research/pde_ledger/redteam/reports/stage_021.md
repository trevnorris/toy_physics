---
unit_id: 021
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T12:30:55-06:00
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 021 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt`

## What the script claims to verify

The two scripts purport to verify a five-part chain for the moving-throat reduced one-port normal form (Stage 4 / Maxwell + mixed-sector reduction): (I) the mixed 4+1 Maxwell fields `E_w = F_w0` and `C_a = F_aw` are gauge invariant; (II) a reduced wall (Q) + brane-like Maxwell (A) + mixed mode (W) Lagrangian gives the Schur-complement self-energy `Sigma_cons = (g_A^2 W_ker + 2 g_A g_W R + g_W^2 A_ker)/Delta` with `Delta = A_ker W_ker - R^2`, and its low-frequency `omega^0/omega^2/omega^4` coefficients match a generic toy rational; (III) dressing the mixed block by a retarded port `Pi_out` transfers the odd part to the wall with coefficient `N(0) = (Omega_A^2 g_W + R g_A)^2 / (Omega_A^2 Omega_W^2 - R^2)^2 >= 0`; (IV) the compact outgoing l=2 normalized response equals `1 + a^2 omega^2/(9 c_s^2) + 4 a^4 omega^4/(81 c_s^4) + i a^5 omega^5/(27 c_s^5)`, so the leading outgoing imaginary coefficient is `Gamma_5_port = a^5/(27 c_s^5)`; (V) a derivative-coupled scalar mixed outlet with `g_A = 0`, `g_W = eta omega` and `Pi_0 = i gamma_1 omega` converts the naive `i omega` port law into an `i omega^3` wall correction with coefficient `gamma_1 eta^2 Omega_A^4 / (Omega_A^2 Omega_W^2 - R^2)^2`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 88 | `simplify(Ewp - Ew) == 0` | yes |
| A2 | sympy | 89 | `simplify(Cap - Ca) == 0` | yes |
| A3 | sympy | 128 | `EQ_Q.lhs + M Q'' + K Q - gA A - gW W == 0` (EL for Q) | yes |
| A4 | sympy | 129 | `EQ_A.lhs + A'' + Omega_A^2 A - R W - gA Q == 0` (EL for A) | yes |
| A5 | sympy | 130 | `EQ_W.lhs + W'' + Omega_W^2 W - R A - gW Q == 0` (EL for W) | yes |
| A6 | sympy | 148 | `Aker A_sol - R W_sol - gA == 0` (matrix inverse residual) | yes |
| A7 | sympy | 149 | `Wker W_sol - R A_sol - gW == 0` (matrix inverse residual) | yes |
| A8-10 | sympy | 161-163 | `z0/z2/z4` of generic toy rational match closed form | yes |
| A11-13 | sympy | 175-177 | `Sigma z0/z2/z4` under subs match generic formulas | yes |
| A14 | sympy | 223-225 | `N_omega - (Aker gW + R gA)^2 / Delta^2 == 0` | yes |
| A15 | sympy | 228-231 | `N(0) - (OA^2 gW + R gA)^2 / (OA^2 OW^2 - R^2)^2 == 0` | yes |
| A16 | sympy | 278-287 | `Y2_hat == 1 + a^2 omega^2/(9 c_s^2) + ... + i a^5 omega^5/(27 c_s^5)` | yes |
| A17 | sympy | 290 | `Gamma5_port - a^5/(27 c_s^5) == 0` | yes |
| A18 | sympy | 323-326 | `N_scalar - eta^2 OA^4 omega^2 / Delta_0^2 == 0` | yes |
| A19 | sympy | 332-335 | `Pi0 N_scalar - i gamma1 eta^2 OA^4 omega^3 / Delta_0^2 == 0` | yes |
| B1 | mathematica | 53 | `eWp - eW == 0` | yes |
| B2 | mathematica | 54 | `cAp - cA == 0` | yes |
| B3 | mathematica | 81 | `Coefficient[lVel, vq^2] - m/2 == 0` | **no — tautological** |
| B4 | mathematica | 93 | `aKer aSol - r wSol - gA == 0` | yes |
| B5 | mathematica | 94 | `wKer wSol - r aSol - gW == 0` | yes |
| B6-8 | mathematica | 102-104 | toy `z0/z2/z4` formulas match closed form | yes |
| B9-11 | mathematica | 115-117 | `Sigma z0/z2/z4` under subs match generic | yes |
| B12 | mathematica | 140 | `nOmega - (aKer gW + r gA)^2 / delta^2 == 0` | yes |
| B13 | mathematica | 141 | `n0 - (oA^2 gW + r gA)^2/(oA^2 oW^2 - r^2)^2 == 0` | yes |
| B14 | mathematica | 163-166 | `Y2_hat - (...closed form...) == 0` | yes |
| B15 | mathematica | 167 | `gamma5Port - radius^5/(27 cS^5) == 0` | yes |
| B16 | mathematica | 187 | `nSeries - eta^2 oA^4 omega^2/(oA^2 oW^2 - r^2)^2 == 0` | yes |
| B17 | mathematica | 188 | `deltaD0 - I gamma1 eta^2 oA^4 omega^3/(oA^2 oW^2 - r^2)^2 == 0` | yes |

Notable mismatch: the SymPy script verifies all three Euler-Lagrange equations (A3-A5) of the reduced Lagrangian via `euler_equations(Lred, ...)`. The Mathematica script has no corresponding check — only the tautological B3 — so it never derives the EOMs from `lRed` at all.

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:67-81`

**What's wrong:**

In Section II of the Mathematica script the reduced Lagrangian `lRed = 1/2 m D[q,t]^2 - 1/2 k q^2 + 1/2 D[a,t]^2 - 1/2 oA^2 a^2 + 1/2 D[ww,t]^2 - 1/2 oW^2 ww^2 + r a ww + gA q a + gW q ww` is defined at lines 67-72. The only assertion that touches this Lagrangian is:

```
expectZero["Q kinetic coefficient", Coefficient[lVel, vq^2] - m/2];   (* line 81 *)
```

where `lVel` is `lRed` with the time derivatives replaced by symbolic placeholders `vq, va, vw`. Because the Lagrangian was *written* with the literal coefficient `1/2 m D[q,t]^2`, `Coefficient[lVel, vq^2]` is `m/2` by construction; the residual is identically zero independent of any physics. This is a `tautological_check` posing as the Lagrangian-sector test.

By contrast the SymPy script (lines 123-130) actually computes the three Euler-Lagrange equations via `euler_equations(Lred, Q, [t])` etc. and verifies each one against the expected EOM:

```
EQ_Q = euler_equations(Lred, Q, [t])[0]
expect_zero("Q equation", EQ_Q.lhs + M*sp.diff(Q,t,2) + K*Q - gA*A - gW*W)
expect_zero("A equation", EQ_A.lhs + sp.diff(A,t,2) + OA**2*A - R*W - gA*Q)
expect_zero("W equation", EQ_W.lhs + sp.diff(W,t,2) + OW**2*W - R*A - gW*Q)
```

The Mathematica script omits the EL derivation entirely. Consequently the only piece of the Lagrangian → EOM map exercised on the Mathematica side is the kinetic-term coefficient, which is trivially built into the source. If a sign error or missing cross-coupling were introduced into `lRed`, the SymPy `Q/A/W equation` checks would catch it; the Mathematica side would silently pass.

**Why this matters:**

Section II is the algebraic heart of unit 021 — every downstream object (`sigmaCons`, `sigmaFull`, `nOmega`, `n0`) is supposed to follow from the Lagrangian. If Mathematica never verifies that the Lagrangian gives the claimed EOMs, the "second-engine" status of Section II is reduced to a check of the Schur complement of an *asserted* 2x2 matrix, with no audit of whether that matrix came from the stated Lagrangian. A typo in `lRed` (e.g. `gA q a` written as `gA q^2 a` or a missing factor of 1/2) would not be detected by Mathematica.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`, replace the tautological assertion at line 81 with three real Euler-Lagrange checks. The EL operator for a Lagrangian `L[q[t], q'[t], a[t], a'[t], ww[t], ww'[t]]` and variable `q` is `D[D[L, D[q,t]], t] - D[L, q]`. Add (after the `lRed` definition, before the Section II.2 `aKer`/`wKer` block):

```
elQ = D[D[lRed, D[q,t]], t] - D[lRed, q];
elA = D[D[lRed, D[a,t]], t] - D[lRed, a];
elW = D[D[lRed, D[ww,t]], t] - D[lRed, ww];
expectZero["Q equation", elQ - (m D[q,t,t] + k q - gA a - gW ww)];
expectZero["A equation", elA - (D[a,t,t] + oA^2 a - r ww - gA q)];
expectZero["W equation", elW - (D[ww,t,t] + oW^2 ww - r a - gW q)];
```

Delete the existing line 81 (`expectZero["Q kinetic coefficient", Coefficient[lVel, vq^2] - m/2]`) and the no-longer-needed `staticL`, `staticTmp`, `staticBack`, `qd`, `ad`, `wd`, `lVel` assignments at lines 73-79 (they are dead once the EL checks are added).

**Verification:**

After the fix the new output should contain three `PASS` lines: `PASS: Q equation`, `PASS: A equation`, `PASS: W equation`, and should no longer contain `Q kinetic coefficient`. If `lRed` is correct each residual will simplify to 0. Trivial-case sanity: at `q=a=ww=0` every term in both `elX` and the canonical EOM is 0, so the residual is 0; with nonzero coordinates the residual is also 0 only when both sides agree term-by-term, which exercises the Lagrangian.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:83-141`

**What's wrong:**

The `.wl` script does not derive the Stage-4 algebra independently; it transliterates the `.py` step-by-step. The intermediate-variable choreography is essentially identical between the two engines. Three corresponding excerpts (SymPy → Mathematica):

Section II.2, self-energy:
```
# sympy lines 132-136
Aker = OA**2 - omega**2
Wker = OW**2 - omega**2
Delta = sp.simplify(Aker * Wker - R**2)
Sigma_cons = sp.simplify((gA**2 * Wker + 2 * gA * gW * R + gW**2 * Aker) / Delta)
A_sol = sp.simplify((gA * Wker + gW * R) / Delta)
W_sol = sp.simplify((gW * Aker + gA * R) / Delta)
```
```
(* mathematica lines 83-89 *)
aKer = oA^2 - omega^2;
wKer = oW^2 - omega^2;
delta = FullSimplify[aKer wKer - r^2, ...];
sigmaCons = FullSimplify[(gA^2 wKer + 2 gA gW r + gW^2 aKer)/delta, ...];
aSol = FullSimplify[(gA wKer + gW r)/delta, ...];
wSol = FullSimplify[(gW aKer + gA r)/delta, ...];
```

Section III, port dressing:
```
# sympy lines 211-214
Sigma_full = sp.simplify((gA**2 * (Wker - Pi) + 2 * gA * gW * R + gW**2 * Aker) / (Aker * (Wker - Pi) - R**2))
Sigma_first = sp.expand(sp.series(Sigma_full, Pi, 0, 2).removeO())
N_omega = sp.simplify((Sigma_first - Sigma_cons) / Pi)
```
```
(* mathematica lines 130-132 *)
sigmaFull = FullSimplify[(gA^2 (wKer - piOut) + 2 gA gW r + gW^2 aKer)/(aKer (wKer - piOut) - r^2), ...];
sigmaFirst = Expand[Normal[Series[sigmaFull, {piOut, 0, 1}]]];
nOmega = FullSimplify[(sigmaFirst - sigmaCons)/piOut, ...];
```

Section IV, l=2 fingerprint:
```
# sympy lines 261-269
j2a = (3/za**3 - 1/za)*sin(za) - 3*cos(za)/za**2
y2a = -(3/za**3 - 1/za)*cos(za) - 3*sin(za)/za**2
h2a = sp.simplify(j2a + I*y2a)
Lambda2 = sp.simplify((k * sp.diff(h2a, za) / h2a).subs(za, k * a))
Lambda2_series = sp.series(Lambda2, k, 0, 7).removeO()
Y2 = sp.simplify(sp.series(1/Lambda2_series, k, 0, 6).removeO())
Y2_static = sp.simplify(Y2.subs(k, 0))
Y2_hat = sp.simplify(sp.expand(Y2/Y2_static))
```
```
(* mathematica lines 149-156 *)
j2a = ((3/za^3) - 1/za) Sin[za] - 3 Cos[za]/za^2;
y2a = -((3/za^3) - 1/za) Cos[za] - 3 Sin[za]/za^2;
h2a = FullSimplify[j2a + I y2a, ...];
lambda2 = FullSimplify[(kWave D[h2a, za]/h2a) /. za -> kWave radius, ...];
lambda2Series = Normal[Series[lambda2, {kWave, 0, 6}]];
y2 = Normal[Series[1/lambda2Series, {kWave, 0, 5}]] // FullSimplify;
y2Static = FullSimplify[y2 /. kWave -> 0, ...];
y2Hat = Expand[y2/y2Static];
```

Same variable names, same order of operations, same `Series → take coefficient → divide by lower-order term` pattern. The Mathematica script could (and should) at minimum derive `aSol/wSol` via `LinearSolve` of the 2x2 mixed system, get `nOmega` via `D[sigmaFull, piOut] /. piOut -> 0` (analytic derivative, not Series), and use the built-in `SphericalHankelH1[2, z]` instead of constructing `j2 + I y2` by hand. Each replacement uses different Mathematica machinery than the SymPy path and so constitutes a genuine independent re-derivation.

**Why this matters:**

The "second engine" policy is supposed to guarantee that if both scripts agree, the result is not an artifact of one CAS's symbolic-simplification quirks. If the second script just re-walks the first's algebra in a different surface syntax, that guarantee degrades to "the same algebra simplifies to the same form" — a much weaker statement. In particular, an error in the *derivation path* shared by both engines (e.g., wrong Schur complement, wrong order of the Pi expansion, wrong spherical-Hankel convention) would propagate to both and be missed by their agreement.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`, replace three derivation steps with their Mathematica-native independent computations. Do **not** delete the existing assertions — the closed-form targets remain the same; only the route to them changes.

1. Section II.2 (lines 88-89): replace the hand-written
   ```
   aSol = FullSimplify[(gA wKer + gW r)/delta, ...];
   wSol = FullSimplify[(gW aKer + gA r)/delta, ...];
   ```
   with a matrix-inversion derivation:
   ```
   matEAW = {{aKer, -r}, {-r, wKer}};
   solAW = LinearSolve[matEAW, {gA, gW}];
   aSol = FullSimplify[solAW[[1]], Assumptions -> $Assumptions];
   wSol = FullSimplify[solAW[[2]], Assumptions -> $Assumptions];
   sigmaConsDerived = FullSimplify[gA aSol + gW wSol, Assumptions -> $Assumptions];
   expectZero["sigmaCons from LinearSolve matches closed form", sigmaConsDerived - sigmaCons];
   ```
   (`matEAW` is the conservative kernel matrix for the (A, W) subsystem at fixed omega; solving against the source vector `{gA, gW}` gives `aSol, wSol` directly without writing the Schur-complement inverse by hand.)

2. Section III (line 132): replace
   ```
   sigmaFirst = Expand[Normal[Series[sigmaFull, {piOut, 0, 1}]]];
   nOmega = FullSimplify[(sigmaFirst - sigmaCons)/piOut, ...];
   ```
   with an analytic-derivative computation:
   ```
   nOmega = FullSimplify[D[sigmaFull, piOut] /. piOut -> 0, Assumptions -> $Assumptions];
   ```
   (i.e. take the first-order Pi coefficient via differentiation rather than via series expansion).

3. Section IV (lines 149-151): replace the hand-built spherical Bessel functions
   ```
   j2a = ((3/za^3) - 1/za) Sin[za] - 3 Cos[za]/za^2;
   y2a = -((3/za^3) - 1/za) Cos[za] - 3 Sin[za]/za^2;
   h2a = FullSimplify[j2a + I y2a, Assumptions -> $Assumptions];
   ```
   with the built-in spherical Hankel function of the first kind:
   ```
   h2a = SphericalHankelH1[2, za];
   ```
   Leave the subsequent `lambda2 = (kWave D[h2a, za]/h2a) /. za -> kWave radius` line and the rest of the section unchanged; `SphericalHankelH1[2, za]` is mathematically identical to `j2 + I y2` at l=2 but is derived through Mathematica's native special-function machinery, not by hand-coded power-series formulas that mirror the `.py`.

**Verification:**

After the fix the script must still exit 0, the existing assertions (`A exact solution residual`, `W exact solution residual`, `N(omega) compact formula`, `N(0) positive-square form`, `Y2_hat minimal branch`, `Gamma5_port - a^5/(27 c_s^5)`) must continue to pass, and there must be one *new* PASS line `PASS: sigmaCons from LinearSolve matches closed form` confirming the matrix-inversion route reproduces the closed-form `sigmaCons`. The Mathematica output's `Lambda2(k)` line should still report `-3/radius + (kWave^2 radius)/3 + (kWave^4 radius^3)/9 + (I/9) kWave^5 radius^4 - (2 kWave^6 radius^5)/27` (since `SphericalHankelH1[2, .]` and `j2 + I y2` agree exactly). The `Gamma5_port` line should still print `radius^5/(27 cS^5)`.

## Independent-derivation check (Mathematica)

A `.wl` file is present, but as documented in F2 above it is a step-by-step transliteration of the `.py` rather than an independent re-derivation. The variable names (aKer ↔ Aker, wKer ↔ Wker, delta ↔ Delta, sigmaCons ↔ Sigma_cons, sigmaFull ↔ Sigma_full, sigmaFirst ↔ Sigma_first, nOmega ↔ N_omega, n0 ↔ N0, j2a ↔ j2a, y2a ↔ y2a, h2a ↔ h2a, lambda2 ↔ Lambda2, y2 ↔ Y2, y2Hat ↔ Y2_hat) and the order of operations match line-for-line. F2 records this and gives the minimal set of Mathematica-native replacements (LinearSolve, analytic derivative, SphericalHankelH1) that would make the second engine genuinely independent.

## Engine cross-check

Where both engines do verify the same closed-form targets, they agree:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `Sigma_cons numerator` | `gA^2(OW^2-ω^2) + 2 gA gW R + gW^2(OA^2-ω^2)` | `gW^2(-oA^2+ω^2) + gA^2(ω-oW)(ω+oW) - 2 gA gW r` (equivalent up to sign factoring; `Sigma_cons` denominator also differs only by overall sign) |
| `z0^(EM+mix)` | `(OA^2 gW^2 + OW^2 gA^2 + 2 R gA gW)/(OA^2 OW^2 - R^2)` | `(gW^2 oA^2 + gA^2 oW^2 + 2 gA gW r)/(oA^2 oW^2 - r^2)` |
| `N(0)` | `(OA^2 gW^2 + 2 OA^2 R gA gW + R^2 gA^2)/(OA^2 OW^2 - R^2)^2` (expanded square) | `(gW oA^2 + gA r)^2/(-(oA^2 oW^2) + r^2)^2` (squared form; identical after expansion since the denominator squared loses the sign) |
| `Y2_hat` | `1 + a^2 ω^2/(9 c_s^2) + 4 a^4 ω^4/(81 c_s^4) + i a^5 ω^5/(27 c_s^5)` | `1 + ω^2 radius^2/(9 cS^2) + 4 ω^4 radius^4/(81 cS^4) + (I/27) ω^5 radius^5/cS^5` |
| `Gamma5_port` | `a^5/(27 c_s^5)` | `radius^5/(27 cS^5)` |
| `N_scalar leading` | `OA^4 eta^2 ω^2 / (OA^4 OW^4 - 2 OA^2 OW^2 R^2 + R^4)` | `eta^2 oA^4 ω^2/(-(oA^2 oW^2) + r^2)^2` |
| `Pi_0 N_scalar` | `i OA^4 eta^2 gamma_1 ω^3 / Delta_0^2` | `(I eta^2 gamma_1 oA^4 ω^3)/(-(oA^2 oW^2) + r^2)^2` |

All paired entries are algebraically identical (after recognizing `(−x)^2 = x^2` for the squared denominators and the symbol renaming `OA ↔ oA`, `OW ↔ oW`, `R ↔ r`, `a ↔ radius`, `c_s ↔ cS`). No engine_disagreement finding.

## Verdict justification

The SymPy script is substantive across all five sections: gauge invariance is checked symbolically, Euler-Lagrange equations are derived from the Lagrangian and checked against the expected EOMs, the Schur-complement self-energy is verified by an independent matrix-inverse residual check, the low-frequency `z0/z2/z4` coefficients are matched against a generic toy rational *and* against the EM/mixed substitution dictionary, the port dressing's first-order coefficient `N(omega)` is reduced to a manifestly nonnegative squared form, the l=2 outgoing fingerprint is built from the standard spherical Hankel function and the leading imaginary coefficient is identified, and the derivative-coupled scalar variant correctly demotes the odd port law from `i omega` to `i omega^3`. I attacked: a possible sign error in `Sigma_cons = -toy` (no — `Delta = D0 − S2 ω^2 + ω^4` carries a `+1` ω^4 coefficient so Sigma_cons = toy, no sign flip); a possible `Pi` series-vs-derivative discrepancy in Section III (no — the linear term in Series equals the derivative); a possible spherical-Hankel sign-convention mismatch (no — `h^{(1)} = j + i y` is consistent with the `+i omega^5` outgoing convention asserted); a possible parity error in the `gA=0, gW=eta omega` substitution of Section V (no — `(Aker gW)^2 ∝ ω^2`, the resulting wall correction goes as `ω·ω^2 = ω^3`). The SymPy side holds up. The Mathematica side reproduces the same closed forms but (a) replaces the entire Euler-Lagrange derivation with a tautological kinetic-coefficient check (F1) and (b) transliterates the rest of the algebra rather than re-deriving it via Mathematica-native machinery (F2). Verdict: `findings`, not stop-cold — the math is correct, but the second engine needs to be made genuinely second.

## Self-test notes

I checked variable-independence for the proposed F1 EL operators: `lRed` depends on `q[t], a[t], ww[t]` and their time derivatives, so `D[lRed, q]`, `D[lRed, D[q,t]]`, etc. are non-degenerate and the EL combinations match the canonical EOMs term-by-term, with the trivial-case substitution `q=a=ww=0` giving `0 - 0 = 0` (consistent). I checked parity in Section V — the scalar `(Aker gW)^2 = (Aker eta omega)^2` is even in omega, so the leading `omega^2` term is real and the subsequent `Pi_0 = i gamma_1 omega` multiplication gives `i omega^3`, consistent with the claim. I confirmed that `SphericalHankelH1[2, z]` and `j2(z) + i y2(z)` agree (both engines' small-z expansion `Lambda2 = -3/a + (k^2 a)/3 + (k^4 a^3)/9 + ...` matches exactly), so the F2 substitution preserves all downstream checks. I confirmed that `D[sigmaFull, piOut] /. piOut -> 0` equals the first-order Series coefficient algebraically, so the F2 derivative substitution preserves the existing `N(omega) compact formula` and `N(0) positive-square form` assertions.
