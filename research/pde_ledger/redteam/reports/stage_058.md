---
unit_id: 058
batch: III.2
auditor_model: claude-opus-4-7
audit_date: 2026-05-26
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage058_coupled_support_source_operator.md
  paper_appendix: present
---

# Audit unit 058 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_058.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage058_coupled_support_source_operator.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 94)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage058_coupled_support_source_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` reads: "Fixed-point equation \eqref{eq:app-stage058-Pe-fixed} and bracket \eqref{eq:app-stage058-bracket}." The bottom-line claims are (i) the self-consistent closure `Pe_* = Xi Delta(Pe_*; kappa, eta)` and (ii) the exact endpoint bracket `Xi Delta_0(kappa,eta) <= Pe_* <= Xi Delta_inf(kappa,eta)`, with the closed-form endpoint functions

- `Delta_0 = eta (cosh(alpha) - 1) / (alpha^2 (alpha sinh(alpha) + eta cosh(alpha)))`,
- `Delta_inf = (cosh(alpha) + (eta/alpha) sinh(alpha) - 1) / (alpha sinh(alpha) + eta cosh(alpha))`,

where `alpha := sqrt(kappa)`. The notes (sections 3-6) add several intermediate deliverables: the exact Green-kernel difference `K_(kappa,eta)(x) = G(1,x) - G(0,x)` that arises from the BVP `-Phi_xx + kappa Phi = Sigma, Phi_x(0) = eta Phi(0), Phi_x(1) = 0`; the strict kernel monotonicity `dK/dx > 0`; the covariance-representation monotonicity `dDelta/dPe = Cov_Pe(K, x) > 0` on the constructive branch; the weak-coupling branch law `Pe_* = Xi Delta_0 + O(Xi^2)`; and the IVT-based proof of the bracket via `F(Xi Delta_0) <= 0` and `F(Xi Delta_inf) >= 0` for `F(Pe) := Pe - Xi Delta(Pe; kappa, eta)`. The part-03 appendix row (line 94) restates the same fixed-point equation and bracket as the unit's status.

## What the script claims to verify

Both scripts verify a parallel chain of identities: (1) the kernel `K_(alpha,eta)(x)` is written in closed form and its derivative `dK/dx` matches the stated identity; (2) the Stage-39 family `Sigma_Pe` normalizes to one on `[0,1]`; (3) the antiderivatives `Fc, Fs` reproduce the integrands `exp(Pe x) cosh(alpha x)` and `exp(Pe x) sinh(alpha x)`; (4) the `Pe = 0` limit of `Delta` matches the closed-form `Delta_0`, and that closed form equals `int_0^1 K dx`; (5) `Delta_inf := K(1)` matches the stated closed form, and equals `lim_{Pe -> oo} Delta`; (6) the bracket gap `Delta_inf - Delta_0` has a stated closed-form, and is positive at nine numeric sample points `(alpha, eta) in {1/10, 1, 3} x {1/10, 1, 10}`; (7) the small-`Pe` series has constant term `Delta_0` and a nonvanishing first-order coefficient at `alpha = eta = 1`. The Mathematica script adds one independent check that `Integrate[kernel * sigmaPe, {x,0,1}]` reduces to the combination form built from `Ic, Is`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Fixed-point equation `Pe_* = Xi Delta(Pe_*; kappa, eta)` (eq:app-stage058-Pe-fixed) | Script defines `Pe_lo, Pe_hi` from `Xi*Delta_0`, `Xi*Delta_inf` and prints them; no fixed-point existence check, no `F(Pe_lo) <= 0`, no `F(Pe_hi) >= 0` | partial |
| Bracket `Xi Delta_0 <= Pe_* <= Xi Delta_inf` (eq:app-stage058-bracket) | Bracket *gap* `Delta_inf - Delta_0` is checked at 9 numeric sample points only; the bracket *theorem* (IVT/F-sign argument) is never exercised | partial |
| `Delta_0` closed form (eq:app-stage058-Delta0) | `expect_zero("Delta0 formula", Delta0 - Delta0_expected)` and integral identity check | match |
| `Delta_inf` closed form (eq:app-stage058-Deltainfty) | `expect_zero("Delta_inf direct substitution (sanity)", ...)` and `Delta_inf as Pe -> oo limit` | match |
| Green-kernel construction `K = G(1,x) - G(0,x)` satisfying the BVP with BCs (notes §3) | Kernel asserted as a closed-form ansatz; BVP solution not verified | missing |
| Kernel monotonicity `dK/dx > 0` (notes §4) | Closed-form identity check on `dK/dx`; positivity of the numerator (`alpha sinh(alpha x) + eta cosh(alpha x) + alpha sinh(alpha(1-x))`) under `alpha, eta > 0, x in [0,1]` not asserted | partial |
| Covariance monotonicity `dDelta/dPe = Cov_Pe(K, x) > 0` (notes §4) | Not checked | missing |
| Weak-coupling branch law `Pe_* = Xi Delta_0 + O(Xi^2)` (notes §6) | `Delta(Pe)` constant term equals `Delta_0`; first-order coefficient is shown nonzero at one point. Translation to `Pe_*` via the fixed-point equation not exercised. | partial |

`paper_alignment: partial` — the closed-form expressions for `Delta_0` and `Delta_inf` match the paper exactly, but the central paper deliverables (the fixed-point equation and the IVT-based bracket existence) are not actually verified by either engine.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 45-49 | `expect_zero` on `Kprime - (numerator)/W` | Notes §4 monotonicity (algebraic side only; sign not checked) | partial |
| A2 | sympy | 53 | `expect_zero` on `int_0^1 Sigma dx - 1` | Notes §3 `Sigma_Pe` definition | yes |
| A3 | sympy | 58-59 | antiderivative regression on `Fc, Fs` | Auxiliary `Ic, Is` (notes §3) | yes |
| A4 | sympy | 74 | `Delta0 - Delta0_expected` | Paper eq:app-stage058-Delta0 | yes |
| A5 | sympy | 75 | `Delta0 - int_0^1 K dx` | Notes §4 uniform-source identity | yes |
| A6 | sympy | 80 | `Delta_inf - Delta_inf_expected` | Paper eq:app-stage058-Deltainfty | yes |
| A7 | sympy | 94 | `bracket_gap - bracket_gap_expected` | Algebraic rewrite of `Delta_inf - Delta_0`; not a sign proof | partial |
| A8 | sympy | 95-102 | numeric loop `val > 0` at 9 points | Bracket non-emptiness, but only sampled | partial |
| A9 | sympy | 107 | `Delta_inf_limit - Delta_inf_expected` (as `Pe -> oo`) | Notes §4 sharp-bottom-source endpoint | yes |
| A10 | sympy | 112 | constant term of small-`Pe` series equals `Delta_0` | Notes §6 weak-coupling at `Pe=0` | partial |
| A11 | sympy | 115-119 | first-order coefficient nonzero at `alpha=eta=1` | Notes §6 (existence of `O(Xi^1)` correction, sampled at one point) | partial |
| M1-M11 | mathematica | parallel | Same chain as A1-A11, plus `delta independent integral matches combination form` (line 73) | Independent integration check is the one genuinely cross-engine derivation | partial |

There is no assertion that exercises the fixed-point equation `Pe_* = Xi Delta(Pe_*; kappa, eta)` itself, nor the IVT bracket existence proof, nor the BVP that defines the Green kernel, nor the covariance-form monotonicity.

## Findings

### F1 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:82-102`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:96-114`

**What's wrong:**
The paper's `\stagefield{Output}` is two things: (i) the fixed-point equation `Pe_* = Xi Delta(Pe_*; kappa, eta)` and (ii) the bracket `Xi Delta_0 <= Pe_* <= Xi Delta_inf`. The notes (§5) prove the bracket via continuity and the sign of `F(Pe) := Pe - Xi Delta(Pe;kappa,eta)` at the two endpoints (`F(Xi Delta_0) <= 0` and `F(Xi Delta_inf) >= 0`). Neither script verifies any of this. Both scripts only:

- compute `Pe_lo = Xi*Delta_0` and `Pe_hi = Xi*Delta_inf` and print them (sympy line 83-86; mathematica line 96-99),
- assert the algebraic rewrite `bracket_gap = Delta_inf - Delta_0` equals a hand-written closed form (sympy line 89-94; mathematica line 101-106) — this is a definitional rewrite, not a sign statement,
- sweep `bracket_gap > 0` at 9 numeric points `(alpha, eta) in {1/10, 1, 3} x {1/10, 1, 10}` (sympy line 95-102; mathematica line 107-114).

A 9-point numeric sweep is not a symbolic positivity proof; missed regimes (e.g., alpha very large with small eta, or eta dominating kappa) are not exercised. More importantly, the bracket *theorem* requires `Delta_0 <= Delta(Pe;kappa,eta) <= Delta_inf` for all `Pe >= 0`, plus continuity of `F` — neither inequality nor continuity is asserted. The fixed-point existence argument (root in the interval by IVT on `F`) is never executed.

**Why this matters:**
The paper's headline result — that the branch point `Pe_*` is trapped in `[Xi Delta_0, Xi Delta_inf]` — has no script-side verification. The scripts verify only the *formulas* `Delta_0`, `Delta_inf` and that `Delta_inf - Delta_0 > 0` at nine sample points. A regression that perturbed `Delta`'s definition would not necessarily be caught.

**Required change:**
Add symbolic / numeric verification of the bracket-existence chain:
1. For each `(alpha, eta) in {1/10, 1, 3} x {1/10, 1, 10}` AND each `Pe in {1/2, 1, 3, 10}`, evaluate `Delta(Pe; alpha, eta) - Delta_0(alpha, eta)` and `Delta_inf(alpha, eta) - Delta(Pe; alpha, eta)` numerically and `assert val >= 0` (within tolerance).
2. For each `(alpha, eta, Xi)` over the same `(alpha, eta)` grid with `Xi in {1/2, 1, 2}`, define `F(Pe) := Pe - Xi * Delta(Pe; alpha, eta)` and assert `F(Xi * Delta_0(alpha, eta)) <= 0` and `F(Xi * Delta_inf(alpha, eta)) >= 0` numerically.
3. Expand the `(alpha, eta)` grid to include one large-alpha point (e.g., `alpha = 10`) and one small-alpha point (e.g., `alpha = 1/100`) so the positivity sweep covers the asymptotic regimes.

**Verification:**
New `expect_close` / `If[val < 0, Exit[1]]` blocks appear after the existing bracket-gap sweep; the transcript shows new `PASS` lines for monotonicity-of-Delta and F-sign IVT checks.

### F2 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:38-44`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:37-44`

**What's wrong:**
The notes §3 derive the support-drop kernel `K_(kappa,eta)(x) = G(1,x) - G(0,x)` as the difference of Green's-function values for the BVP

```
-Phi_xx + kappa Phi = Sigma,  Phi_x(0) = eta Phi(0),  Phi_x(1) = 0.
```

Both scripts simply *write down* the closed-form

```
K = (cosh(alpha x) + (eta/alpha) sinh(alpha x) - cosh(alpha(1-x))) / (alpha sinh(alpha) + eta cosh(alpha))
```

(sympy line 39-41; mathematica line 38-41) and then verify only its derivative identity. Neither engine verifies that the Phi-field built from this kernel and the source `Sigma_Pe` actually satisfies the BVP. The entire Green-kernel construction — the load-bearing intermediate result the notes use to derive `Delta` — is asserted, not derived.

**Why this matters:**
If the kernel were transcribed with a sign error, a wrong factor of `eta/alpha`, or a Robin/Neumann BC swap, `Delta_0` and `Delta_inf` would still satisfy their algebraic identities (because the same wrong `K` defines both), but the resulting fixed-point equation would not correspond to the operator the paper claims to solve. The script cannot catch a kernel-side error.

**Required change:**
Add a BVP-verification block in both scripts. Solve the BVP `-Phi_xx + alpha^2 Phi = Sigma_Pe(x), Phi_x(0) = eta Phi(0), Phi_x(1) = 0` symbolically using `sp.dsolve` (SymPy) or `DSolve[]` (Mathematica) with `Sigma_Pe(x) = Pe exp(Pe x)/(exp(Pe)-1)` as the inhomogeneity. Call the resulting solution `Phi_solved(x)`. Then assert:
1. `expect_zero("BVP end-to-end drop matches Delta", Phi_solved.subs(x, 1) - Phi_solved.subs(x, 0) - Delta)` (where `Delta` is the existing closed-form `Delta(Pe;alpha,eta)` already computed in the script).

This single assertion is sufficient: it confirms that the kernel ansatz produces the same end-to-end drop as the actual BVP solution, without requiring a full symbolic Green's-function expansion. The Robin BC and Neumann BC are encoded in `dsolve`/`DSolve`, so any sign error in the kernel's BC handling would surface as a mismatch.

If `dsolve` / `DSolve` cannot produce a closed form with the symbolic `Sigma_Pe`, fall back to a numeric BVP solve at the same `(alpha, eta, Pe)` sample grid used in F1.

**Verification:**
New assertion `expect_zero("Phi(1) - Phi(0) from BVP matches Delta", ...)` (or Mathematica equivalent) appears at a fresh line in each script; transcript shows new `PASS`.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:42-49`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:42-49`

**What's wrong:**
The notes §4 state `dK/dx > 0` for `alpha > 0, eta > 0, x in [0,1]`. Both scripts assert the algebraic identity `dK/dx == (alpha sinh(alpha x) + eta cosh(alpha x) + alpha sinh(alpha(1-x))) / W` — an algebraic rewrite, not a positivity proof. The positivity of the numerator (which under the script's `x in [0,1]` and `alpha, eta > 0` assumptions follows from `sinh(alpha x) >= 0`, `cosh(alpha x) > 0`, `sinh(alpha(1-x)) >= 0`) is never asserted. A sign flip in any one of the three numerator terms would still pass the algebraic identity but break monotonicity.

**Why this matters:**
Kernel monotonicity is a load-bearing premise for the bracket argument (it implies, together with the covariance representation, `Delta_0 <= Delta(Pe) <= Delta_inf` on the constructive branch). Without a positivity check, the monotonicity is a stated claim, not a verified one.

**Required change:**
After the `Kprime identity` check, add numeric positivity sweep of the kernel numerator over the same `(alpha, eta)` grid used in F1 and `x in {0, 1/4, 1/2, 3/4, 1}`. Assert each value is `>= 0`. Same in Mathematica with `Table[...] // Flatten` followed by `AnyTrue[#, # < 0 &]`.

**Verification:**
A new `kernel numerator positivity sweep = PASS` line appears in the transcript.

### F4 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:109-120`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:123-132`

**What's wrong:**
Notes §6 states the weak-coupling branch law as a statement about the fixed-point root `Pe_*(Xi)`: `Pe_*(Xi) = Xi Delta_0(kappa, eta) + O(Xi^2)`. The scripts verify (i) the constant term of the small-`Pe` series of `Delta(Pe)` (not `Pe_*(Xi)`) equals `Delta_0`, and (ii) the first-order Pe-coefficient is nonzero at a single point `alpha = eta = 1`. Neither check is the weak-coupling branch law. The first is the consistency `Delta(0) = Delta_0` (already verified independently at sympy line 74); the second is a different statement entirely (sensitivity of `Delta` to `Pe`, not sensitivity of `Pe_*` to `Xi`).

**Why this matters:**
Conflating the small-`Pe` expansion of `Delta` with the small-`Xi` expansion of the fixed-point solution `Pe_*(Xi)` is a category error. The two are related (via the implicit function theorem at `Xi = 0`, `Pe_* = 0`), but the paper's claim is about `Pe_*(Xi)`, and that claim has no current verification.

**Required change:**
Add a symbolic implicit-function-theorem check at `Xi = 0`:
1. Define `F(Pe, Xi; alpha, eta) := Pe - Xi * Delta(Pe; alpha, eta)`. Verify `F(0, 0; alpha, eta) == 0` (trivial; `Delta(0) = Delta_0` is finite, so `F(0, 0) = 0`).
2. Compute `dPe_*/dXi|_{Xi=0}` as `-dF/dXi / dF/dPe`, evaluated at `(Pe, Xi) = (0, 0)`. This gives `dPe_*/dXi|_{Xi=0} = Delta(0; alpha, eta) / 1 = Delta_0(alpha, eta)`.
3. `expect_zero` on `(dPe_*/dXi)|_{Xi=0} - Delta_0(alpha, eta)`.

**Verification:**
A new `expect_zero("weak-coupling branch slope = Delta_0", ...)` line appears in the transcript with `= 0`.

## Independent-derivation check (Mathematica)

The `.wl` file is structurally parallel to the `.py` file:

- both define `W = alpha*Sinh[alpha] + eta*Cosh[alpha]` (sympy line 38; mathematica line 37);
- both define `K`/`kernel` via the same closed-form ratio (sympy line 39-41; mathematica line 38-41);
- both check `kernelPrime - (numerator)/W` (sympy line 45-49; mathematica line 46-49);
- both define antiderivatives `Fc, Fs`/`fc, fs` with the *same* structural form, though sympy writes the antiderivative explicitly at line 56-57 while mathematica uses `Integrate[Exp[Pe*x]*Cosh[alpha*x], x]` at line 55-56. That `Integrate` call is one genuine independent derivation step (the engine chooses the antiderivative rather than echoing a hand-written one).
- both numerical sweeps use exactly the same grid `(alpha, eta) in {1/10, 1, 3} x {1/10, 1, 10}` (sympy line 95-101; mathematica line 107-114) — identical sample points;
- both use the same series order `Series[Delta, {Pe, 0, 1}]` and the same `alpha = eta = 1` check point.

The Mathematica script does add one genuinely independent derivation step at line 65-73: `delta = Integrate[kernel * sigmaPe, {x, 0, 1}]` is computed directly and then verified against `deltaCombination` (the Ic/Is form). Given that independence point, plus the engine-chosen `Integrate` for `Fc/Fs`, I would not file a `mathematica_transliteration` finding — the script does exercise alternate paths — but the parallelism is otherwise quite strong and the numerical sweep grids being identical is worth noting for the future.

## Engine cross-check

Both scripts pass all stated assertions in their saved transcripts. The closed forms reported for `Delta_0`, `Delta_inf`, `bracket_gap`, and `Delta(Pe -> oo)` are identical between the two engines (up to trivial sign normalization, e.g., SymPy's denominator `(Pe^2 - alpha^2)` vs Mathematica's `(alpha^2 - Pe^2)` in `Ic`, which is compensated by a sign in the numerator). No engine disagreement.

## Verdict justification

`paper_alignment: partial` and `verdict: findings`. The script does verify the *closed forms* of `Delta_0` and `Delta_inf` against the paper exactly, and the engines agree. What it does not verify is the load-bearing paper deliverables that those formulas plug into: the BVP that the Green kernel solves (F2), the IVT-based bracket-existence argument (F1), the monotonicity inequalities `Delta_0 <= Delta(Pe) <= Delta_inf` (F1, F3), and the implicit-function-theorem weak-coupling branch law (F4). These are four `insufficient_verification` findings, each at the level of "the assertion exercises an algebraic side-identity instead of the physical claim." No `paper_misalignment` was found: the paper card and the script's targets are aligned; the script just under-verifies them. No `stop_cold`: the framework can proceed via Codex-applied script additions.

Attacks tried and failed: I attempted to spot a sign error in `Ic, Is` (the antiderivative regression catches it at sympy line 58-59 / mathematica line 57-58); I attempted a `kappa = alpha^2` substitution mismatch (the scripts and notes both use `alpha` consistently as `sqrt(kappa)`, and the paper card uses the same convention at line 25); I attempted a Robin-BC sign-error spot-check (the kernel's `eta/alpha` factor lines up with `Phi_x(0) = eta Phi(0)` — the linear-BC sign convention is consistent across paper, notes, and scripts); I attempted to find a value mismatch in the closed forms between paper, notes, and script (none — `Delta_0` and `Delta_inf` match exactly). The audit issues are about depth of verification, not target correctness.

## Self-test notes

For F1, I confirmed the proposed `F(Pe) := Pe - Xi*Delta(Pe; alpha, eta)` endpoint-sign check is non-trivial under the existing assumptions: the residual at `Pe = Xi*Delta_0` is nonzero in general (`F = Xi*Delta_0 - Xi*Delta(Xi*Delta_0)`), and only vanishes at the fixed point itself, so the inequality `F(Xi*Delta_0) <= 0` is exercising the monotonicity-in-Pe of `Delta`. For F2, the BVP residual `-Phi_xx + alpha^2 Phi - Sigma` involves `x` non-trivially (because `Sigma` and `Phi` both depend on `x`), so the second derivative is not identically zero — this is not a `sp.diff(EXPR, VAR)` where `VAR` is missing from `EXPR`. For F4, the implicit function theorem evaluation at `(Pe, Xi) = (0, 0)` requires `dF/dPe|_{(0,0)} = 1 != 0`, which is satisfied (independent of `Delta`). No symmetric-integral parity traps: the integration domain `[0,1]` is asymmetric, so no odd/even cancellation risk.
