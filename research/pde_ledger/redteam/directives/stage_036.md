---
unit_id: 036
batch: II.1
created_at: 2026-05-25T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-25T23:52:03-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 036 (v2)

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

This is a v2 directive that supersedes the previously-applied v1 directive (whose Applied blocks are above-the-fold in the prior file state). The v1 fixes are already in the current scripts; the two findings here are additional defects caught in the paper-grounded re-audit.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:96-97,140-149`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:108-109,139-148`

**Issue:**
The "admissible / inadmissible M_mix witness" checks for the paper-card feasibility leg `M_mix <= G(xi_req,delta)` are constructed so that they cannot fail regardless of what `G_sample` is. The sympy script (lines 96-97, 140-149) defines `Mmix_good = G_sample - 1/10` and `Mmix_bad = G_sample + 1/10` and then asserts `Mmix_good <= G_sample` and `Mmix_bad > G_sample`. These are arithmetic identities on the rationals: `x - 1/10 <= x` and `x + 1/10 > x` are unconditionally true. They never touch `M_mix_expr` (the actual paper definition `8 Chi^2/(pi^2 A Omega_U^2 Delta_0)`), and they would still pass if the closed forms of `G` or `M_mix` were corrupted. The Mathematica script mirrors the same construction at lines 108-109 and 139-148 (`mMixGood = gSample - 1/10`, `mMixBad = gSample + 1/10`).

The fix is to construct the sample `M_mix` from independent microscopic parameters `(Chi, Omega_U, Delta_0, A)` — i.e. by substituting into the existing `Mmix_expr` / `mMix` — and then numerically compare to `G_sample`. At `(delta=1, xi=1/2)`, `G_sample = 27/58 ≈ 0.4655`. Picking `(Chi, Omega_U, Delta_0, A) = (1, 1, 1, 29)` gives `Mmix = 8/(29 pi^2) ≈ 0.0280`, an admissible witness. Picking `(1, 1, 1, 1)` gives `Mmix = 8/pi^2 ≈ 0.8106`, an inadmissible witness.

**Required change:**

(Part A — SymPy.) Replace lines 96-97 and 140-149 of `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`.

Before (line 96-97):
```python
Mmix_good = sp.simplify(G_sample - sp.Rational(1, 10))
Mmix_bad = sp.simplify(G_sample + sp.Rational(1, 10))
```

After (line 96-97 replacement — define independent-parameter witnesses derived from Mmix_expr, not from G_sample):
```python
Mmix_admissible = sp.N(Mmix_expr.subs({Chi: 1, OmegaU: 1, Delta0: 1, A: 29}))
Mmix_inadmissible = sp.N(Mmix_expr.subs({Chi: 1, OmegaU: 1, Delta0: 1, A: 1}))
G_sample_n = sp.N(G_sample)
```

Then, before (lines 140-149):
```python
expect_true(
    "admissible sample: M_mix <= G(xi_req,delta)",
    bool(Mmix_good <= G_sample),
    f"M_mix={Mmix_good}, G={G_sample}",
)
expect_true(
    "inadmissible sample: support deficit blocks the branch",
    bool(Mmix_bad > G_sample),
    f"M_mix={Mmix_bad}, G={G_sample}",
)
```

After (lines 140-149 replacement — compare numerically against the independent-parameter Mmix evaluations):
```python
expect_true(
    "admissible sample: M_mix < G(xi_req,delta)",
    bool(Mmix_admissible < G_sample_n),
    f"M_mix={Mmix_admissible}, G={G_sample_n}",
)
expect_true(
    "inadmissible sample: support deficit blocks the branch",
    bool(Mmix_inadmissible > G_sample_n),
    f"M_mix={Mmix_inadmissible}, G={G_sample_n}",
)
```

Note: the two `expect_true` calls in the source are separated by other code (line 116-139 are the F-related kappa checks). Edit only the two `expect_true` blocks at lines 140-144 and 145-149, and the two `Mmix_*` definitions at lines 96-97. Do not edit the kappa block in between.

(Part B — Mathematica.) Replace lines 108-109 and 139-148 of `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`.

Before (line 108-109):
```wolfram
mMixGood = FullSimplify[gSample - 1/10];
mMixBad = FullSimplify[gSample + 1/10];
```

After (line 108-109 replacement):
```wolfram
mMixAdmissible = N[mMix /. {Chi -> 1, OmegaU -> 1, Delta0 -> 1, A -> 29}];
mMixInadmissible = N[mMix /. {Chi -> 1, OmegaU -> 1, Delta0 -> 1, A -> 1}];
gSampleN = N[gSample];
```

Then, before (lines 139-148):
```wolfram
expectTrue[
  "admissible sample: M_mix <= G(xi_req,delta)",
  mMixGood <= gSample,
  "M_mix=" <> fmt[mMixGood] <> ", G=" <> fmt[gSample]
];
expectTrue[
  "inadmissible sample: support deficit blocks the branch",
  mMixBad > gSample,
  "M_mix=" <> fmt[mMixBad] <> ", G=" <> fmt[gSample]
];
```

After (lines 139-148 replacement):
```wolfram
expectTrue[
  "admissible sample: M_mix < G(xi_req,delta)",
  mMixAdmissible < gSampleN,
  "M_mix=" <> fmt[mMixAdmissible] <> ", G=" <> fmt[gSampleN]
];
expectTrue[
  "inadmissible sample: support deficit blocks the branch",
  mMixInadmissible > gSampleN,
  "M_mix=" <> fmt[mMixInadmissible] <> ", G=" <> fmt[gSampleN]
];
```

Self-test (already performed by auditor; reproduced here so Codex can spot-check):
- `Mmix_expr = 8 Chi^2 / (pi^2 A Omega_U^2 Delta_0)`; at `(Chi=1, OmegaU=1, Delta0=1, A=29)`, `Mmix = 8/(29 pi^2) ≈ 0.02797`. `G_sample = 27/58 ≈ 0.46552`. So `0.02797 < 0.46552`, admissible witness passes.
- At `(Chi=1, OmegaU=1, Delta0=1, A=1)`, `Mmix = 8/pi^2 ≈ 0.81057`. `0.81057 > 0.46552`, inadmissible witness passes.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 036` and `redteam exec-mathematica 036`, confirm both scripts exit 0, and confirm both saved outputs contain admissible / inadmissible lines whose `M_mix=` value is approximately `0.02797` / `0.81057` (not the old `53/145` / `82/145`).

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- summary: Replaced the G-derived M_mix witness values with independent M_mix evaluations from microscopic parameters.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:67,86` (insertion points for comments)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:58,98` (insertion points for comments)

**Issue:**
The assertion at sympy lines 68-71 (`gBreq_sq_over_varpi2 - (pi^2 A/8)(G - Mmix_expr) == 0`) and its `xi -> xi_req` echo at lines 87-90 are identically zero by construction: both `a_req` and `G` are hardcoded closed forms (lines 54, 60) typed in such a way that `(pi^2 A/8) G = a_req` and `(pi^2 A/8) Mmix_expr = alpha_mix`, making the residual `(a_req - alpha_mix) - (a_req - alpha_mix) = 0`. The Mathematica script has the same pattern at lines 59-62 and 99-102.

The closed form of `G` is genuinely anchored elsewhere — by the symbolic kappa derivation at sympy:123-139 / mathematica:124-138, which proves `F = R_target_sym` and thereby fixes the algebraic structure of the F-G branch. So the local tautology is not load-bearing, but it should be labelled as such so future readers (and the next audit pass) do not double-count it as an independent verification.

**Required change:**

(Part A — SymPy.) Insert a comment block above the assertion at line 68 and above the assertion at line 87 of `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`.

Before (line 67-71):
```python
print("R_target =", R_target)
expect_zero(
    "g_B,req^2/varpi^2 - (pi^2 A / 8) (G - M_mix)",
    gBreq_sq_over_varpi2 - (sp.pi**2 * A / 8) * (G - Mmix_expr),
)
```

After (line 67-71 replacement; insert a 4-line comment between `print("R_target =", R_target)` and `expect_zero(`):
```python
print("R_target =", R_target)
# Definitional self-consistency: a_req and G are both hardcoded closed forms,
# so this just confirms a_req = (pi^2 A / 8) G after Mmix_expr cancels.
# The genuine anchor for F (and hence the closed form of G via the same algebra)
# is the symbolic kappa derivation below ("symbolic kappa derivation: F(xi,delta) - R_target_sym").
expect_zero(
    "g_B,req^2/varpi^2 - (pi^2 A / 8) (G - M_mix)",
    gBreq_sq_over_varpi2 - (sp.pi**2 * A / 8) * (G - Mmix_expr),
)
```

Before (line 86-90):
```python
print("Parametric frontier: xi -> (F(xi,delta), G(xi,delta))")
expect_zero(
    "final-test support inequality <-> nonnegative required support loading",
    (sp.pi**2 * A / 8) * (G.subs(xi, xi_req) - Mmix_expr) - gBreq_sq_over_varpi2.subs(xi, xi_req),
)
```

After (insert a 1-line comment):
```python
print("Parametric frontier: xi -> (F(xi,delta), G(xi,delta))")
# Same definitional identity as above, with xi -> xi_req. Definitional, not load-bearing.
expect_zero(
    "final-test support inequality <-> nonnegative required support loading",
    (sp.pi**2 * A / 8) * (G.subs(xi, xi_req) - Mmix_expr) - gBreq_sq_over_varpi2.subs(xi, xi_req),
)
```

(Part B — Mathematica.) Insert analogous comments above the assertions at line 59 and line 99 of `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`.

Before (lines 58-62):
```wolfram
expectZero["G - closed form", g - gTarget];
expectZero[
  "g_B,req^2/varpi^2 - (Pi^2 A / 8) (G - M_mix)",
  gBReqSqOverVarpi2 - (Pi^2*A/8)*(gTarget - mMix)
];
```

After (lines 58-62; insert a 4-line comment between the two expectZero calls):
```wolfram
expectZero["G - closed form", g - gTarget];
(* Definitional self-consistency: alphaReq and gTarget are both hardcoded
   closed forms, so this just confirms alphaReq = (Pi^2 A / 8) gTarget after mMix cancels.
   The genuine anchor for F (and hence the closed form of G via the same algebra)
   is the symbolic kappa derivation below (rTargetSym). *)
expectZero[
  "g_B,req^2/varpi^2 - (Pi^2 A / 8) (G - M_mix)",
  gBReqSqOverVarpi2 - (Pi^2*A/8)*(gTarget - mMix)
];
```

Before (lines 98-102):
```wolfram
  A > 0 && delta > 0 && 0 <= xiReq < 1 && OmegaU > 0 && Delta0 > 0 && beta0 > 0 && NQ > 0;
expectZero[
  "final-test support inequality <-> nonnegative required support loading",
  (Pi^2*A/8)*((gTarget /. xi -> xiReq) - mMix) - (gBReqSqOverVarpi2 /. xi -> xiReq)
];
```

After (insert a 1-line comment):
```wolfram
  A > 0 && delta > 0 && 0 <= xiReq < 1 && OmegaU > 0 && Delta0 > 0 && beta0 > 0 && NQ > 0;
(* Same definitional identity as above, with xi -> xiReq. Definitional, not load-bearing. *)
expectZero[
  "final-test support inequality <-> nonnegative required support loading",
  (Pi^2*A/8)*((gTarget /. xi -> xiReq) - mMix) - (gBReqSqOverVarpi2 /. xi -> xiReq)
];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 036` and `redteam exec-mathematica 036`, confirm exit 0, and confirm both saved outputs continue to print residuals of 0 for the four assertions. The comments themselves don't appear in the output transcripts; the verifier instead spot-reads the source files at the named line ranges to confirm the comments landed verbatim.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- summary: Added source comments labelling the two support-loading residual checks as definitional consistency checks.
- deviation: none
