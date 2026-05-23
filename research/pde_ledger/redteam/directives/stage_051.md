---
unit_id: 051
batch: III.2
created_at: 2026-05-22T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-22T17:11:48-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 051

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py:127-131`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:75-80`

**Issue:**
The `M_mix(Z_W^(twin,req)) - G_tr/2 == 0` assertion is tautological by construction. The script defines `Mmix_from_ZW` as the algebraic inverse of the linear coefficient that produced `ZW_twin_req` from `G_tr/2`, so `Mmix_from_ZW.subs(ZW, ZW_twin_req)` reduces to `f^{-1}(f(G_tr/2)) = G_tr/2` independent of whether the coefficient `pi^2 (1-eps_eta)(1-eps) / [8 (1+chi0)^2]` is the correct one. The docstring (sympy L20-23) claims the assertion checks the Stage 047/030 `M_mix <-> Z_W` bridge explicitly, but the current check does not.

**Required change:**

In the SymPy script, at `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py:127-131`, replace the block:

```python
# Verify Z_W threshold against the Stage-30 coherent map.
Mmix_from_ZW = sp.simplify(8 * ZW * (1 + chi0) ** 2 / (pi**2 * (1 - eps_eta) * (1 - eps)))
expect_zero(
    "M_mix(Z_W^(twin,req)) - G_tr/2",
    Mmix_from_ZW.subs(ZW, ZW_twin_req) - Gtr / 2,
)
```

with:

```python
# Stage 047/030 coherent forward map: Z_W = pi^2 (1-eps_eta)(1-eps) M_mix / [8 (1+chi0)^2].
ZW_from_Mmix = sp.simplify(pi**2 * (1 - eps_eta) * (1 - eps) * Mmix / (8 * (1 + chi0) ** 2))
# Apply the forward map to the M_mix threshold M_mix = G_tr/2 and compare to ZW_twin_req.
ZW_threshold_via_map = sp.simplify(ZW_from_Mmix.subs(Mmix, Gtr / 2))
expect_zero(
    "Z_W^(twin,req) - forward-map(M_mix=G_tr/2)",
    ZW_twin_req - ZW_threshold_via_map,
)
```

The symbol `Mmix` is already declared at sympy L51 — reuse it directly. `Gtr` is at L55, `ZW_twin_req` at L117. Do not change any other line.

In the Mathematica script, at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:75-80`, replace the block:

```mathematica
mMixFromZW = FullSimplify[8 zW (1 + chi0)^2/(Pi^2 (1 - epsEta) (1 - eps)), Assumptions -> $Assumptions];

Print["Lambda_twin,req = ", fmt[lambdaTwinReq]];
Print["M_mix^(twin,req) = ", fmt[mMixTwinReq]];
Print["Z_W^(twin,req) = ", fmt[zWTwinReq]];
expectZero["M_mix(Z_W^(twin,req)) - G_tr/2", (mMixFromZW /. zW -> zWTwinReq) - gTr/2];
```

with:

```mathematica
(* Stage 047/030 coherent forward map: Z_W = pi^2 (1-eps_eta)(1-eps) M_mix / [8 (1+chi0)^2]. *)
zWFromMmix = FullSimplify[Pi^2 (1 - epsEta) (1 - eps) mMix/(8 (1 + chi0)^2), Assumptions -> $Assumptions];
zWThresholdViaMap = FullSimplify[zWFromMmix /. mMix -> gTr/2, Assumptions -> $Assumptions];

Print["Lambda_twin,req = ", fmt[lambdaTwinReq]];
Print["M_mix^(twin,req) = ", fmt[mMixTwinReq]];
Print["Z_W^(twin,req) = ", fmt[zWTwinReq]];
expectZero["Z_W^(twin,req) - forward-map(M_mix=G_tr/2)", zWTwinReq - zWThresholdViaMap];
```

The symbol `mMix` is already in the `Clear[...]` and `$Assumptions` list (wl L28-32). Do not change any other line.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 051` and `redteam exec-mathematica 051` and confirm the new check line `Z_W^(twin,req) - forward-map(M_mix=G_tr/2) = 0` appears in both output files AND both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl`
- summary: Replaced the inverse-map threshold check with the explicit Stage 047/030 forward-map comparison in both audit scripts.
- deviation: none

## F2 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:40-50`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:82-87`

**Issue:**
The Mathematica script supplies the closed-form `piExpected` (wl L41-45) and the closed-form root `xi2x` (wl L82-85) verbatim from the SymPy script, then checks `piTr - piExpected == 0` and `gTr(xi2x) - 2 mMix == 0` using the same `FullSimplify` engine. This violates the checkpoint-stage second-engine policy: a typo or sign error replicated in both closed forms would not be detected by either. The Mathematica script must take an independent algebraic route for at least these two checks: derive `xi_(2x)` via `Solve` and verify `Pi_tr` cancellation via `Factor`/`Together` rather than `FullSimplify` of the difference.

**Required change:**

In the Mathematica script, at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:40-50`, replace the block:

```mathematica
piTr = FullSimplify[fTr gTr, Assumptions -> $Assumptions];
piExpected = FullSimplify[
  xi (xi + delta) (9 delta + (9 + 2 r) xi)^2 (9 delta + (9 + 2 r^2) xi)/
    (9 (1 - xi) (9 delta^2 + 18 delta xi + (9 + 2 r^2) xi^2)^2),
  Assumptions -> $Assumptions
];

Print["G_tr = ", fmt[gTr]];
Print["F_tr = ", fmt[fTr]];
Print["Pi_tr = ", fmt[piTr]];
expectZero["Pi_tr - expected closed form", piTr - piExpected];
```

with:

```mathematica
piTr = FullSimplify[fTr gTr, Assumptions -> $Assumptions];
(* Independent canonicalization via Factor/Together (no FullSimplify wrapper around the diff). *)
piTrFactored = Factor[Together[fTr gTr]];
piExpected = xi (xi + delta) (9 delta + (9 + 2 r) xi)^2 (9 delta + (9 + 2 r^2) xi)/
    (9 (1 - xi) (9 delta^2 + 18 delta xi + (9 + 2 r^2) xi^2)^2);
piExpectedFactored = Factor[Together[piExpected]];

Print["G_tr = ", fmt[gTr]];
Print["F_tr = ", fmt[fTr]];
Print["Pi_tr = ", fmt[piTr]];
Print["Pi_tr (Factor/Together) = ", fmt[piTrFactored]];
Print["Pi_tr (claim, Factor/Together) = ", fmt[piExpectedFactored]];
expectZero["Pi_tr - expected closed form", piTrFactored - piExpectedFactored];
```

Then, at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:82-87`, replace the block:

```mathematica
xi2x = FullSimplify[
  (2 mMix (9 + 2 r^2) - 9 delta + Sqrt[(2 mMix (9 + 2 r^2) - 9 delta)^2 + 648 mMix delta])/18,
  Assumptions -> $Assumptions
];
Print["xi_(2x) = ", fmt[xi2x]];
expectZero["G_tr(xi_(2x)) - 2 M_mix", (gTr /. xi -> xi2x) - 2 mMix];
```

with:

```mathematica
(* Independently derive the positive root of gTr == 2 mMix via Solve. *)
xi2xRoots = xi /. Solve[gTr == 2 mMix, xi];
xi2xDerived = FullSimplify[
  SelectFirst[xi2xRoots, TrueQ[Simplify[# > 0, Assumptions -> $Assumptions]] &],
  Assumptions -> $Assumptions
];
(* Closed-form claim (the docstring's answer). *)
xi2xClaim = FullSimplify[
  (2 mMix (9 + 2 r^2) - 9 delta + Sqrt[(2 mMix (9 + 2 r^2) - 9 delta)^2 + 648 mMix delta])/18,
  Assumptions -> $Assumptions
];
Print["xi_(2x) (Solve) = ", fmt[xi2xDerived]];
Print["xi_(2x) (claim) = ", fmt[xi2xClaim]];
expectZero["xi_(2x): Solve vs claim", xi2xDerived - xi2xClaim];
expectZero["G_tr(xi_(2x)) - 2 M_mix", (gTr /. xi -> xi2xDerived) - 2 mMix];
```

Do not modify any line outside the ranges named above. Do not change the SymPy script for this finding.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 051` and confirm the new output lines appear in the `.txt` file:
- `Pi_tr (Factor/Together) = ...`
- `Pi_tr (claim, Factor/Together) = ...`
- `xi_(2x) (Solve) = ...`
- `xi_(2x) (claim) = ...`
- `xi_(2x): Solve vs claim = 0` (PASS)

All `expectZero` assertions must PASS and the script must exit 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl`
- summary: Switched the Mathematica product check to Factor/Together canonicalization and derived xi_(2x) through Solve before comparing with the closed-form claim.
- deviation: none
