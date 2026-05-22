---
unit_id: 036
batch: II.1
created_at: 2026-05-21T00:00:00Z
findings_count: 6
stop_cold: null
applied: true
applied_at: 2026-05-21T17:38:15-06:00
findings_applied: 6
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 036

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:72`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:59`

**Issue:**
SymPy line 57 defines `R_target = sp.simplify(NQ * A / (beta0 * (sp.Rational(8) / sp.pi**2)))`, which simplifies to `pi**2 * A * NQ / (8 * beta0)`. Line 72 then asserts that the difference between `R_target` and `sp.pi**2 * A * NQ / (8 * beta0)` is zero — this is `X - X = 0`. The Mathematica script repeats the same construction at lines 50, 59.

**Required change:**

SymPy script, delete line 72 entirely:

Before:
```python
expect_zero("R_target - pi^2 A NQ/(8 beta0)", R_target - sp.pi**2 * A * NQ / (8 * beta0))
```

After:
```
(line deleted; no replacement)
```

Mathematica script, delete line 59 entirely:

Before:
```wolfram
expectZero["R_target - Pi^2 A NQ/(8 beta0)", rTarget - Pi^2*A*NQ/(8*beta0)];
```

After:
```
(line deleted; no replacement)
```

Do not remove the definitions at SymPy line 57 or Mathematica line 50; `R_target` / `rTarget` are still printed for context (SymPy line 67, Mathematica line 56).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 036` and `redteam exec-mathematica 036` and confirm both scripts exit 0 AND that the strings `R_target - pi^2 A NQ/(8 beta0) = 0` (SymPy) and `R_target - Pi^2 A NQ/(8 beta0) = 0` / `PASS: R_target - Pi^2 A NQ/(8 beta0)` (Mathematica) no longer appear in the saved outputs.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- summary: Removed the tautological R_target closed-form assertions while leaving the contextual R_target prints intact.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:68-71,88-91`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:57,60-63,90-93`

**Issue:**
In the SymPy script, both `a_req = 9*pi**2*A*xi*(xi+delta)/(8*(9*delta+11*xi))` (line 60) and `G = 9*xi*(xi+delta)/(9*delta+11*xi)` (line 54) are written down as closed forms. Computing `(pi**2*A/8)*G` mechanically yields `a_req` because the two were typed in that way. Likewise `(pi**2*A/8)*Mmix_expr = alpha_mix` is mechanical. Therefore the assertion at lines 68-71 of the form
`gBreq_sq_over_varpi2 - (pi**2*A/8)*(G - Mmix_expr) == 0`
is `(a_req - alpha_mix) - (a_req - alpha_mix) = 0` by construction — it cannot fail no matter what is assigned to `a_req` and `G`, so long as their relationship is preserved at typing time.

Lines 88-91 substitute `xi -> xi_req` and re-check the same identity; same problem.

The Mathematica script mirrors this at lines 41-43, 49, 51, 60-63, 90-93.

The fix is to add a symbolic kappa-based derivation that arrives at one of these quantities through a different algebraic route, so that the lines 68-71 / 60-63 check stops being a self-restatement and becomes a confirmation that the kappa-derived expression factors as `(pi**2 A/8)*(G - M_mix)`.

**Required change:**

(Part A — SymPy.) Insert the following block immediately after SymPy line 120 (after the existing `R_target_host` check, before the `expect_true("admissible sample: M_mix <= G(xi_req,delta)", ...)` call at line 121):

```python
# Symbolic kappa-based cross-check of the support-feasibility content.
# Builds R_target_sym from the Stage-18 microscopic kappa expansion
# symbolically in (xi, delta, A_sym, beta0_sym) and confirms it equals F.
A_sym, beta0_sym = sp.symbols("A_sym beta0_sym", positive=True, real=True)
x_sym = A_sym * xi
deltaK_sym = A_sym * delta
N_sym = (
    beta0_sym
    * (KAPPA0_SQ * (x_sym + deltaK_sym) + KAPPA1_SQ * x_sym) ** 4
    / (
        KAPPA0_SQ
        * (A_sym - x_sym)
        * (KAPPA0_SQ * (x_sym + deltaK_sym) ** 2 + KAPPA1_SQ * x_sym ** 2) ** 2
    )
)
R_target_sym = sp.simplify(N_sym * A_sym / (beta0_sym * KAPPA0_SQ))
expect_zero(
    "symbolic kappa derivation: F(xi,delta) - R_target_sym",
    sp.simplify(F - R_target_sym),
)
```

Do NOT remove the existing assertions at lines 68-71 or 88-91; once the kappa derivation anchors F (and through it, the closed form of G via the same algebra), those assertions become valid bookkeeping confirmations of the factorization rather than self-restatements.

(Part B — Mathematica.) Insert the following block immediately after Mathematica line 114 (after the existing `expectZero["admissible sample: F(xi_req,delta) - R_target(host)", ...]`), before line 115's `expectTrue["admissible sample: M_mix <= G(xi_req,delta)", ...]`:

```wolfram
(* Symbolic kappa-based cross-check: derive R_target_sym symbolically and confirm = F. *)
Clear[Asym, beta0Sym];
$Assumptions =
  Element[{A, delta, xiReq, Chi, OmegaU, Delta0, beta0, NQ, Asym, beta0Sym}, Reals] &&
  A > 0 && delta > 0 && 0 <= xiReq < 1 && OmegaU > 0 && Delta0 > 0 && beta0 > 0 && NQ > 0 &&
  Asym > 0 && beta0Sym > 0;
xSym = Asym * xi;
deltaKSym = Asym * delta;
nSym = beta0Sym*(kappa0Sq*(xSym + deltaKSym) + kappa1Sq*xSym)^4 /
  (kappa0Sq*(Asym - xSym)*(kappa0Sq*(xSym + deltaKSym)^2 + kappa1Sq*xSym^2)^2);
rTargetSym = FullSimplify[nSym*Asym/(beta0Sym*kappa0Sq), Assumptions -> $Assumptions];
expectZero[
  "symbolic kappa derivation: F(xi,delta) - R_target_sym",
  FullSimplify[Together[Expand[fTarget - rTargetSym]], Assumptions -> $Assumptions]
];
```

Do NOT remove the existing assertions at Mathematica lines 57, 60-63, 90-93.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 036` and `redteam exec-mathematica 036`, confirm exit 0, and confirm both saved outputs contain a new PASS line of the form `symbolic kappa derivation: F(xi,delta) - R_target_sym = 0` (SymPy) and `PASS: symbolic kappa derivation: F(xi,delta) - R_target_sym` (Mathematica).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- summary: Added symbolic kappa-derived R_target cross-checks in both audit scripts.
- deviation: none

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:79,83`

**Issue:**
Line 79 defines `alphaCrit = FullSimplify[9*Pi^2*A*(1 + delta)/(8*(11 + 9*delta)), ...]` and line 83 asserts `(Pi^2*A/8)*gMaxTarget - alphaCrit == 0`. Since `gMaxTarget` (line 78) is `9*(1+delta)/(9*delta+11)`, `(Pi^2*A/8)*gMaxTarget` mechanically equals `alphaCrit`. The assertion is `X - X = 0` by construction.

**Required change:**

Delete Mathematica lines 79 and 83.

Before (lines 78-83):
```wolfram
gMaxTarget = FullSimplify[9*(1 + delta)/(9*delta + 11), Assumptions -> delta > 0];
alphaCrit = FullSimplify[9*Pi^2*A*(1 + delta)/(8*(11 + 9*delta)), Assumptions -> A > 0 && delta > 0];

Print["G_max(delta) = ", fmt[gMax]];
expectZero["G_max - closed form", gMax - gMaxTarget];
expectZero["(Pi^2 A / 8) G_max - alpha_crit", (Pi^2*A/8)*gMaxTarget - alphaCrit];
```

After:
```wolfram
gMaxTarget = FullSimplify[9*(1 + delta)/(9*delta + 11), Assumptions -> delta > 0];

Print["G_max(delta) = ", fmt[gMax]];
expectZero["G_max - closed form", gMax - gMaxTarget];
```

(I.e., remove the `alphaCrit = ...` line and the trailing `expectZero["(Pi^2 A / 8) G_max - alpha_crit", ...]` line.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 036`, confirm exit 0, and confirm the saved output no longer contains the line `(Pi^2 A / 8) G_max - alpha_crit = 0` or `PASS: (Pi^2 A / 8) G_max - alpha_crit`.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- summary: Removed the alphaCrit definition and the tautological alpha_crit assertion from the Mathematica audit.
- deviation: none

## F4 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:126`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:120`

**Issue:**
SymPy line 126's assertion is `bool(sp.Rational(9, 10) < 1)`, which is `True` regardless of any physics. No `R_target` is constructed from any parameter setting. The label "inadmissible sample: R_target < 1 blocks the branch" is misleading because no sample exists — only the bare arithmetic `9/10 < 1` is being asserted. Same issue at Mathematica line 120 (`9/10 < 1`).

**Required change:**

SymPy script, delete line 126:

Before:
```python
expect_true("inadmissible sample: R_target < 1 blocks the branch", bool(sp.Rational(9, 10) < 1), "R_target=9/10")
```

After:
```
(line deleted; no replacement)
```

Mathematica script, delete line 120:

Before:
```wolfram
expectTrue["inadmissible sample: R_target < 1 blocks the branch", 9/10 < 1, "R_target=9/10"];
```

After:
```
(line deleted; no replacement)
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 036` and `redteam exec-mathematica 036`, confirm exit 0, and confirm the lines `inadmissible sample: R_target < 1 blocks the branch = R_target=9/10` and the corresponding `PASS:` line are absent from both saved outputs.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- summary: Removed the bare arithmetic inadmissible R_target branch assertions from both audit scripts.
- deviation: none

## F5 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:43,66-69,78,128-129`

**Issue:**
The Mathematica script declares the same closed-form targets that SymPy declares, rather than deriving them inside Mathematica. Specifically:

- Line 43: `gTarget = FullSimplify[9*xi*(xi + delta)/(9*delta + 11*xi), ...]` — written by hand, identical to SymPy line 54.
- Lines 66-69: `dGTarget = FullSimplify[9*(9*delta^2 + 18*delta*xi + 11*xi^2)/(9*delta + 11*xi)^2, ...]` — written by hand, identical to SymPy lines 75-76.
- Line 78: `gMaxTarget = FullSimplify[9*(1 + delta)/(9*delta + 11), ...]` — written by hand, identical to SymPy line 81.
- Line 129: `gSeriesTarget = FullSimplify[xi - 2*xi^2/(9*delta), ...]` — written by hand, identical to SymPy line 135.

Because the targets are typed in, both engines would write down the same wrong target if the SymPy author got the closed form wrong. The fix is to derive each Mathematica target through a Mathematica-native algebraic step so that, if the SymPy target were wrong, the Mathematica derivation would produce a different expression and the cross-check would fail.

**Required change:**

(a) Replace lines 66-71 of the Mathematica script:

Before:
```wolfram
dG = FullSimplify[D[gTarget, xi], Assumptions -> $Assumptions];
dGTarget = FullSimplify[
  9*(9*delta^2 + 18*delta*xi + 11*xi^2)/(9*delta + 11*xi)^2,
  Assumptions -> $Assumptions
];
Print["dG/dxi = ", fmt[dG]];
expectZero["dG/dxi - manifestly positive form", dG - dGTarget];
```

After:
```wolfram
(* Derive dG/dxi from gTarget; multiply through by the squared denominator
   so the result is a polynomial in xi and delta. Mathematica must produce
   11*xi^2 + 18*delta*xi + 9*delta^2 on its own; the closed-form coefficients
   are not declared up front. *)
dG = FullSimplify[D[gTarget, xi], Assumptions -> $Assumptions];
dGPolynomial = FullSimplify[Expand[dG*(9*delta + 11*xi)^2/9], Assumptions -> $Assumptions];
Print["dG/dxi = ", fmt[dG]];
Print["9 * dG/dxi * (9 delta + 11 xi)^2 / 81 (polynomial) = ", fmt[dGPolynomial]];
expectZero[
  "dG/dxi positivity polynomial: 9 dG/dxi (9d+11xi)^2 / 9 == 11 xi^2 + 18 delta xi + 9 delta^2",
  dGPolynomial - (11*xi^2 + 18*delta*xi + 9*delta^2)
];
(* Also confirm the manifest non-negativity of the numerator polynomial
   for delta, xi >= 0: discriminant <= 0 in xi. *)
disc = Discriminant[11*xi^2 + 18*delta*xi + 9*delta^2, xi];
discSimplified = FullSimplify[disc, Assumptions -> delta > 0];
Print["discriminant (in xi) = ", fmt[discSimplified]];
expectZero["dG/dxi numerator discriminant should be 0 (perfect-square-like)", discSimplified - 0];
```

Note: hand-check — discriminant of `11 xi^2 + 18 delta xi + 9 delta^2` in `xi` is `(18 delta)^2 - 4*11*9*delta^2 = 324 delta^2 - 396 delta^2 = -72 delta^2`. So `disc = -72 delta^2`, NOT zero. **Reject the discriminant assertion above.** The correct assertion for nonzero discriminant is `disc + 72*delta^2 == 0`. Use this corrected form:

```wolfram
expectZero["dG/dxi numerator discriminant equals -72 delta^2", discSimplified + 72*delta^2];
```

(b) Replace line 78 of the Mathematica script:

Before:
```wolfram
gMaxTarget = FullSimplify[9*(1 + delta)/(9*delta + 11), Assumptions -> delta > 0];
```

After:
```wolfram
(* Derive gMaxTarget from gMax via substitution rather than declaring it. *)
gMaxTarget = FullSimplify[gMax, Assumptions -> delta > 0];
```

(After F3's removal of line 83, this leaves the `expectZero["G_max - closed form", gMax - gMaxTarget]` at the renumbered next line. Since `gMaxTarget = gMax` after substitution, this assertion becomes trivially `0 = 0`. To restore content, **replace** the existing `expectZero["G_max - closed form", gMax - gMaxTarget]` with the form:

```wolfram
expectZero["G_max - 9(1+delta)/(9delta+11)", gMax - 9*(1 + delta)/(9*delta + 11)];
```

This puts the closed-form `9*(1+delta)/(9*delta+11)` on the test side, where Mathematica's `gMax` (computed via Limit) is being compared against the literal closed form. Now the closed form is the target being verified, not a variable that was declared equal to it.

(c) Replace lines 128-131 of the Mathematica script:

Before:
```wolfram
gSeries = FullSimplify[Normal[Series[gTarget, {xi, 0, 2}]], Assumptions -> delta > 0];
gSeriesTarget = FullSimplify[xi - 2*xi^2/(9*delta), Assumptions -> delta > 0];
Print["G(xi,delta) near xi=0 = ", fmt[gSeries]];
expectZero["G near-onset series through O(xi^2)", gSeries - gSeriesTarget];
```

After:
```wolfram
(* Read the series coefficients directly out of gSeries; do not declare a target. *)
gSeries = FullSimplify[Normal[Series[gTarget, {xi, 0, 2}]], Assumptions -> delta > 0];
c0 = FullSimplify[Coefficient[gSeries, xi, 0], Assumptions -> delta > 0];
c1 = FullSimplify[Coefficient[gSeries, xi, 1], Assumptions -> delta > 0];
c2 = FullSimplify[Coefficient[gSeries, xi, 2], Assumptions -> delta > 0];
Print["G(xi,delta) near xi=0 = ", fmt[gSeries]];
Print["near-onset coefficients: c0=", fmt[c0], ", c1=", fmt[c1], ", c2=", fmt[c2]];
expectZero["near-onset c0 = 0", c0];
expectZero["near-onset c1 = 1", c1 - 1];
expectZero["near-onset c2 = -2/(9 delta)", c2 + 2/(9*delta)];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 036`, confirm exit 0, and confirm the saved output contains the new assertion labels `dG/dxi positivity polynomial: ...`, `dG/dxi numerator discriminant equals -72 delta^2`, `G_max - 9(1+delta)/(9delta+11)`, `near-onset c0 = 0`, `near-onset c1 = 1`, and `near-onset c2 = -2/(9 delta)`, all of which must PASS.

## Applied: F5

- files_changed:
  - `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- summary: Replaced the Mathematica closed-form derivative, endpoint, and near-onset target checks with derived polynomial, limit, and coefficient checks.
- deviation: none

## F6 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:120` (insertion point)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:114` (insertion point)

**Issue:**
The kappa-based identity `F(xi,delta) = R_target_sym(xi,delta,A,beta0)` is only checked at the single numeric sample `(A=3, beta0=5, xi=1/2, delta=1)`. The docstring claims "Exact dimensionless support-feasibility function"; "exact" should mean symbolic across the parameter range, not just at one point.

**Required change:**
This change is fulfilled by the same edits required under F2 (Parts A and B). No additional change required — the F2 insertion adds a symbolic-in-`(xi, delta, A_sym, beta0_sym)` version of the kappa derivation that asserts `F - R_target_sym == 0` for general symbols, thereby promoting the single-sample numeric witness to a symbolic identity.

If Codex has already applied F2, mark F6 as `Applied: F6` with `summary: covered by F2`. If F2 was blocked, F6 is also blocked.

**Verification command:**
After Codex applies, the verifier will confirm the new PASS line `symbolic kappa derivation: F(xi,delta) - R_target_sym = 0` appears in both saved outputs (same check as F2).

## Applied: F6

- files_changed:
  - `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- summary: covered by F2
- deviation: none
