---
unit_id: 083
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 083 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.txt`

## What the script claims to verify

The scripts purport to verify three things for the "Family-1 direct operator window" stage:
(1) the implicit-differentiation formula `dPe_*/dXi = Delta / (1 - Xi*dDelta/dPe)` for the fixed-point branch `Pe = Xi*Delta(Pe)`;
(2) the Family-1 wall constants `Delta_0`, `Delta_inf`, and the auxiliary constant `A_F1 = (kappa + pi^2/4)/(kappa + y^2)` where `y*tan(y) = eta`;
(3) the operator-selected `Pe_+-^(chi)`, `Pe_+-^(J)`, the transport window values `zeta` evaluated at these four Pe points (plus `zeta_max = A_F1*pi^2/4`), and the resulting `Pi/C_mix = 1 + zeta` window values.
The ledger-comment at the end further claims monotonicity of `zeta_phys(F1)` and `Pi/C_mix` in Xi on the stable branch, and that the values reproduce Stage-61/63/64 windows. The actual assertions in the scripts are: one symbolic `expect_zero` (sympy) / `expectZero` (mathematica) on the fixed-point derivative, plus eleven `expectApprox` numerical equality checks (mathematica only).

## Assertion inventory

| #  | Script       | Line       | Form                                                                                                | Anchored to claim? |
|----|--------------|------------|-----------------------------------------------------------------------------------------------------|--------------------|
| A1 | sympy        | 47         | `expect_zero("implicit fixed-point derivative", dPe_formula - dPe_expected)`                       | partial (tautological algebra) |
| A2 | mathematica  | 49         | `expectZero["implicit fixed-point derivative", dpeFormula - dpeExpected]`                          | partial (tautological algebra) |
| A3 | mathematica  | 103        | `expectApprox["Delta_0(F1) numeric check", delta0F1, 1.73302079021525*10^-4, 10^-16]`              | no (self-check)    |
| A4 | mathematica  | 104        | `expectApprox["Delta_inf(F1) numeric check", deltaInfF1, 2.01447565540522*10^-2, 10^-15]`          | no (self-check)    |
| A5 | mathematica  | 105        | `expectApprox["Pe_-^(chi) numeric check", peMinusChi, 96.5285247264385, 10^-10]`                   | no (self-check)    |
| A6 | mathematica  | 106        | `expectApprox["Pe_+^(chi) numeric check", pePlusChi, 11220.5441626259, 10^-7]`                     | no (self-check)    |
| A7 | mathematica  | 107        | `expectApprox["Pe_-^(J) numeric check", peMinusJ, 22.0062226330754, 10^-10]`                       | no (self-check)    |
| A8 | mathematica  | 108        | `expectApprox["Pe_+^(J) numeric check", pePlusJ, 2558.01892349205, 10^-8]`                         | no (self-check)    |
| A9 | mathematica  | 109-113    | five `expectApprox[...zeta_*...]` numeric checks                                                   | no (self-check)    |
| (sympy has no `assert` on items 2/3 — only `print`) | | | | |

A1/A2: confirm SymPy's `solve` (or Mathematica's `Solve`) correctly inverts a linear equation in `Pe'(Xi)` — algebraic identity, not a physics check.
A3..A9: compare a freshly-computed value against a stored decimal that came from the same formula in a previous run — they cannot fail unless the code is edited.
Sympy script (script asserts only A1): the entire claim-set (2)+(3) is reported via `print` with no executable check.

## Findings

### F1 — mathematica_transliteration

**Severity:** high
**Files:**
- `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:51-95`
- `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py:53-104`

**What's wrong:**
The `.wl` script is a near-line-by-line transliteration of the `.py` script. Side-by-side correspondence:

- Constants: sympy lines 53-55 `kappa_F1 = sp.Rational(12321, 5); eta_F1 = sp.Integer(37); alpha_F1 = sp.sqrt(kappa_F1)` vs mathematica lines 51-53 `kappaF1 = 12321/5; etaF1 = 37; alphaF1 = Sqrt[kappaF1]`.
- `Delta_0` closed form (sympy 57-60) vs (mathematica 55-57): both literally write `eta*(cosh(alpha)-1)/(alpha^2*(alpha*sinh(alpha)+eta*cosh(alpha)))`.
- `Delta_inf` closed form (sympy 61-64) vs (mathematica 58-60): both literally write `(cosh(alpha)+(eta/alpha)*sinh(alpha)-1)/(alpha*sinh(alpha)+eta*cosh(alpha))`.
- `y` from `y*tan(y)=eta`, identical initial guess `1.53` (sympy 67, mathematica 62).
- `A_F1` formula `(kappa + pi^2/4)/(kappa + y^2)` (sympy 68, mathematica 63).
- Hardcoded selectors `Theta_chi_coeff = 4.06863235008162`, `Theta_J_coeff = 0.927552032539308`, and the literal `136900` factor (sympy 75-78, mathematica 70-73).
- `Omega` formula `pi*Pe*(2*Pe*exp(Pe)+pi)/((4*Pe^2+pi^2)*(exp(Pe)-1))` (sympy 97, mathematica 88).
- `zeta = A_F1 * Omega^2` (sympy 98, mathematica 89).

Both engines start from identical pre-baked closed forms and only differ in numerical rendering. There is no independent derivation: the Mathematica script does not, e.g., solve the BVP that defines `Delta_0` / `Delta_inf`, does not derive `Omega` from a transport integral, and does not solve any cross-eigenvalue problem to obtain `Theta_chi` / `Theta_J`. As a result, an algebraic error in the SymPy script would be replicated, not caught, by the Mathematica engine.

**Why this matters:**
The whole point of a second engine is to expose copy-time mistakes in the closed forms. A line-for-line port provides no independent attestation; the "engine agreement" reported in the output is uninformative.

**Required change:**
In the Mathematica script, replace at least one of the closed-form constants with an independent derivation from the underlying definition. Concretely, derive `Delta_0` and `Delta_inf` by solving the linear two-point BVP that they come from (see Required change in directive F1), and re-derive `Omega` symbolically from its defining integral (or its defining ratio) rather than by typing in the same closed form. Then compare these independently-derived expressions to the pre-baked closed forms via `expectZero`.

**Verification:**
Verifier sees new `expectZero` lines in the `.wl` titled "Delta_0(F1) independent BVP derivation = 0", "Delta_inf(F1) independent BVP derivation = 0", and "Omega(Pe) independent derivation = 0", each printing residual `= 0` and `PASS:`. The pre-existing `expectApprox` numeric self-checks may stay but are now redundant.

### F2 — hardcoded_result

**Severity:** medium
**Files:**
- `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py:75-78`
- `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:70-73`

**What's wrong:**
The "operator selectors" `Theta_chi_coeff = 4.06863235008162` and `Theta_J_coeff = 0.927552032539308`, plus the constant `136900` used as `Xi_chi = 136900 * Theta_chi_coeff`, are pure literal decimals with no in-script derivation, no in-script citation, and no in-script algebraic constraint they must satisfy. They are simply asserted by typing them in. The same is true on both engines.

The downstream values `Pe_{+/-}^{chi}`, `Pe_{+/-}^{J}` are then computed as `Xi_* * Delta_*` and `Pi/C_mix` from those Pe values — so the entire downstream window is anchored on two unverified decimals. The `expectApprox` checks on lines 105-108 compare the script's own product `136900*Theta_chi_coeff*Delta_0` against a pre-stored decimal (e.g. `96.5285247264385`) that is exactly that same product, so they confirm nothing about the selectors.

**Why this matters:**
Any typo in `Theta_chi_coeff` or `Theta_J_coeff` (e.g. a transposed digit) would silently propagate and the script would still emit `PASS:` on every numeric line, because each check is a self-check. The window the stage claims to "reproduce directly" rests entirely on values the script does not anchor.

**Required change:**
Add a comment block above the four lines (sympy 75-78, mathematica 70-73) naming the upstream source for the three constants: which earlier stage / verified script computes `Theta_chi`, `Theta_J`, and where the `136900 = ?` factor comes from. Then, if those upstream scripts emit a checkable form for the selectors (e.g. a root of a known equation `f(Theta) = 0`), add an `expect_zero("Theta_chi selector residual", f(Theta_chi_coeff))` / `expectZero[...]` in each script. If the only upstream "source" is a paper/notes prose value, leave a `# SOURCE: stage_NN, paper section X.Y` comment so a future auditor can trace it. Do not invent equations; if no upstream script anchors them, state that explicitly: `# WARNING: no upstream-script anchor for Theta_chi`.

**Verification:**
Verifier sees either (a) a new `expect_zero` / `expectZero` block whose `residual = 0` and `PASS:` appears in the saved output, or (b) explicit `# SOURCE:` / `(* SOURCE: *)` comments naming the upstream verifying script for each constant, or (c) an explicit `# WARNING: no upstream-script anchor` if the audit cannot find one.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py:48-122` (entire body after item 1)
- `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:54-119`

**What's wrong:**
The SymPy script asserts exactly one thing (the algebraic rearrangement at line 47). All of items 2 and 3 — the wall constants, the `Pe` values, the `zeta` window, the `Pi/C_mix` window — are reported via `print` with no `assert`. The Mathematica script asserts the same item plus eleven `expectApprox` numeric self-checks, each of the form "value the script just computed agrees with a pre-stored decimal that equals the value the script just computed" (see F4 below for the tautology issue with these specific checks).

Additionally, the ledger comment

  > `zeta_phys(F1) and Pi/C_mix windows are monotone in Xi on the stable branch`

at sympy line 120 / mathematica line 124 is asserted in prose but not exercised: the scripts evaluate `zeta` at four discrete Pe values and at `Pe -> infinity` (`zeta_max`) — they never check the sign of `d zeta / d Pe` between `Pe_-^(chi)` and `Pe_+^(chi)` (or between the J-branch endpoints), so a non-monotone `Omega` would not be caught.

The ledger claim "reproduces the Stage-61/63/64 windows directly" is similarly only checked by self-comparison.

**Why this matters:**
The current SymPy script gives the verifier no executable check on the most physically interesting outputs of the stage. The monotonicity claim being a ledger comment with no verification means an algebra slip in `Omega` could yield a non-monotone window and still PASS.

**Required change:**
In `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py`:

1. After line 73 add an explicit residual check on `Delta_0_F1` and `Delta_inf_F1` against a numeric value that is **independently** recomputed (e.g. by mpmath at 50-digit precision from the closed form — equivalent precision but a different code path) so that a transcription typo in the closed form would surface:
   - `expect_zero("Delta_0_F1 vs mpmath", sp.N(Delta0_F1, 40) - mpmath_value)` with `mpmath_value` obtained by importing mpmath and re-evaluating the same formula via `mpmath.cosh` etc. at 50-digit precision (different evaluator, same formula — catches sympy-vs-mpmath bugs, not formula-vs-truth, but at least gives one executable check).
2. After line 110, add an explicit monotonicity check between `Pe_minus_chi` and `Pe_plus_chi`:
   - `dzeta_dPe = sp.diff(zeta_F1, Pe_sym)`
   - For a small sample of intermediate Pe values, assert `sp.N(dzeta_dPe.subs(Pe_sym, Pe_test), 30) > 0`. (zeta is monotone increasing in Pe over the relevant range; if the comment is "monotone in Xi on the stable branch" and Pe = Xi*Delta is monotone in Xi, this implies monotone in Pe.)
3. In the Mathematica script after line 113, add the same sign check via `Sign[D[zetaF1[pp], pp] /. pp -> peTest]` for the same sample points; assert `Sign == 1`.

If the user prefers not to introduce mpmath, change item 1 to an `expect_zero` on the identity `eta_F1*(cosh(alpha_F1)-1) - alpha_F1**2 * (alpha_F1*sinh(alpha_F1) + eta_F1*cosh(alpha_F1)) * Delta0_F1 == 0` — i.e., verify `Delta_0` is the root of the defining linear equation it claims to solve, not just an algebraic shorthand the script typed in.

**Verification:**
Verifier sees new `expect_zero` / `expectZero` lines in the SymPy output for `Delta_0_F1 defining equation = 0`, `Delta_inf_F1 defining equation = 0`, and `zeta_F1 monotone at Pe=...` (one line per test point), each printing `= 0` (or `Sign == 1`) and `PASS:`. Saved output grows correspondingly.

### F4 — tautological_check

**Severity:** medium
**Files:**
- `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:103-113`

**What's wrong:**
The eleven `expectApprox` checks compare a value the script just computed (`delta0F1`, `deltaInfF1`, `peMinusChi`, etc.) against a literal decimal target that is exactly what that computation produces. For example, on line 103:

  `expectApprox["Delta_0(F1) numeric check", delta0F1, 1.73302079021525*10^-4, 10^-16];`

The target `1.73302079021525*10^-4` is `N[etaF1*(Cosh[alphaF1]-1)/(alphaF1^2*(alphaF1*Sinh[alphaF1]+etaF1*Cosh[alphaF1])), 15]` — the script computes the symbolic form, then numerically compares it against its own numerical value to ~15 digits. The check is mathematically incapable of failing unless someone edits the formula and forgets to re-edit the target decimal. It is not testing the physics; it is testing that the file hasn't been edited in a particular pattern.

The same applies to lines 104-113.

The `expectZero` on line 49 is also borderline: it confirms `Solve` correctly inverts the linear equation `dpeFormula - dpeExpected == 0`, but both forms are algebraic rearrangements of the same `D[fixedPoint, xi] == 0` identity. It checks that Mathematica's `Solve` works as documented; it does not check a non-trivial physical claim. Left alone (low severity) since "algebra cross-check" is a useful sanity test even if not strictly a physics check.

**Why this matters:**
A reader of the output sees eleven `PASS:` lines and infers that eleven independent numeric facts were verified. None of them were. The signal-to-noise ratio of the output is poor and false confidence accumulates.

**Required change:**
Either (a) delete the eleven self-check `expectApprox` lines 103-113, since they don't add information after F1's independent-derivation checks land, or (b) replace each numeric target with an **independently** evaluated value — for example, compute `Delta_0(F1)` symbolically here (current line 55) AND additionally via the BVP-solution route added under F1, then compare the two values to `10^-30`. The latter is preferred because it forces the comparison to be against a different derivation, not the same formula.

If Codex chooses (a), it must verify F1's directive added the BVP/independent-derivation `expectZero` lines first; otherwise the unit loses all numeric coverage. Therefore (b) is the safer mechanical fix.

**Verification:**
Verifier sees either no more `expectApprox` self-checks (if (a)), or `expectApprox` lines whose targets are computed expressions referencing an independent variable created earlier in the script (if (b)). Saved output reflects the change.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script, as detailed in F1. Three corresponding-section quotes:

- Delta_0 closed form
  - sympy 57-60: `Delta0_F1 = sp.simplify(eta_F1 * (sp.cosh(alpha_F1) - 1) / (alpha_F1**2 * (alpha_F1 * sp.sinh(alpha_F1) + eta_F1 * sp.cosh(alpha_F1))))`
  - mathematica 55-57: `delta0F1 = FullSimplify[etaF1*(Cosh[alphaF1] - 1)/(alphaF1^2*(alphaF1*Sinh[alphaF1] + etaF1*Cosh[alphaF1]))]`

- Omega formula
  - sympy 97: `Omega = sp.simplify(pi * Pe_sym * (2 * Pe_sym * sp.exp(Pe_sym) + pi) / ((4 * Pe_sym**2 + pi**2) * (sp.exp(Pe_sym) - 1)))`
  - mathematica 88: `omega[pp_] := Pi*pp*(2*pp*Exp[pp] + Pi)/((4*pp^2 + Pi^2)*(Exp[pp] - 1));`

- A_F1
  - sympy 68: `A_F1 = sp.simplify((kappa_F1 + pi**2 / 4) / (kappa_F1 + y_F1**2))`
  - mathematica 63: `aF1 = N[(kappaF1 + Pi^2/4)/(kappaF1 + yRoot^2), 40];`

Verdict: the `.wl` is not an independent derivation. Same algebraic choreography, different syntax. This is the `mathematica_transliteration` finding F1.

## Engine cross-check

Both engines emit the same numerical values to displayed precision:

| Quantity | SymPy | Mathematica |
|----------|-------|-------------|
| Delta_0(F1)   | 0.000173302079021525149057156196550 | 0.0001733020790215251490571561965499... |
| Delta_inf(F1) | 0.0201447565540521594271032956099 | 0.0201447565540521594271032956099177... |
| y_F1          | 1.52948248371469963657887092268 | 1.5294824837146996499271076224024... |
| A_F1          | 1.00005192880219520364707725466 | 1.0000519288021953286593340837139... |
| Pi_max^(F1)/C_mix | 3.46752922945601192487420445402 | 3.467529229456012233329584501570... |

Numerical agreement to ~14-15 digits (mostly limited by SymPy's `Float`/`nsolve` precision). This agreement is **not informative** because both engines compute identical formulas; an error in either closed form would not be caught.

## Verdict justification

`findings`. The sole substantive symbolic assertion (item 1, fixed-point derivative) is an algebraic identity that confirms `Solve`/`solve` correctly inverts a linear equation in `Pe'(Xi)`; the eleven Mathematica `expectApprox` checks are pure self-checks against pre-baked decimals; the SymPy script asserts nothing about items 2 and 3 at all; both engines run identical algebra so "agreement" is uninformative; the monotonicity claim in the ledger is not exercised. None of the findings rise to `UNFIXABLE` or `CRITICAL_DOWNSTREAM` — the constants `Theta_chi`, `Theta_J`, the factor `136900`, the closed forms for `Delta_0`/`Delta_inf`/`Omega`, and the wall data may be correct (and probably are; the agreement with later stages is suggestive), but the scripts as written do not verify them.

Attacks tried that did *not* surface additional findings:
- Symbol-domain check: SymPy's `Xi`, `Pe_sym` declared `positive=True, real=True`; Mathematica's `xi > 0` and `kappa, eta, alpha` real/positive. After the symbols are reassigned to concrete numbers (`kappa=12321/5` etc.), the assumption space is irrelevant; no `simplify` is invoked under a contradictory assumption.
- Branch-cut hunt: `sp.cosh`, `sp.sinh` are entire; `sqrt(kappa)` with `kappa = 12321/5 > 0` is unambiguous; `sp.exp(Pe)` is entire. No hidden branch issues.
- Stale-output: sympy script Apr 1 12:39 < output May 11 12:44; mathematica script Apr 21 17:04 < output May 11 13:00. Both outputs are fresh — no `stale_output` finding.
- Engine disagreement: numerical outputs match to the precision SymPy carries; no disagreement.

## Self-test notes

Checked: variable-independence (the symbolic check at line 47/49 is over `Xi` with `Pe = pe(Xi)`; the derivative is not identically zero, the equation is well-posed). Parity: no integrals over unbounded symmetric domains in this script. Trivial-case: substituted `eta -> 0` mentally into `Delta_0` formula — gives `0`, consistent (Delta_0 vanishes when the no-flux source `eta` is off). Path specification: F1/F2/F3/F4 all edit existing files; no `missing_verification_script` finding, so no new-path issue. The monotonicity-check additions in F3 are designed so the derivative `dzeta/dPe` is over `Pe_sym` (which Omega genuinely depends on, line 97) and the assertion is checked at concrete numeric Pe-values in `(0, infinity)` where the analytic form is positive — non-tautological and non-trivial.
