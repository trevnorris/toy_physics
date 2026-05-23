---
unit_id: 059
batch: III.2
created_at: 2026-05-22T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-22T17:51:13-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 059

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:63`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl:62`

**Issue:**
The check `expect_zero("residual bracket center identity", R_hi - R_lo - (zeta_hi - zeta_lo))` (sympy:63) and its Mathematica twin (wl:62) are algebraic identities given that the *previous* two lines define `R_lo = zeta_req - zeta_hi` and `R_hi = zeta_req - zeta_lo` (sympy:59-60, wl:51-52). The assertion cannot fail no matter what Omega, A_K, or any physical quantity equals.

**Required change:**
Delete the offending lines.

SymPy — remove line 63 in
`scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`:
```
expect_zero("residual bracket center identity", R_hi - R_lo - (zeta_hi - zeta_lo))
```

Mathematica — remove line 62 in
`mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`:
```
expectZero["residual bracket center identity", rHi - rLo - (zetaHi - zetaLo)];
```

Leave the `print`/`Print` of `R_-`, `R_+`, `zeta_-`, `zeta_+` in place (they are informational, not assertions). Do not introduce a replacement check in this directive — the finding's primary requirement is removal of the tautology.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 059` and `redteam exec-mathematica 059`. The lines `residual bracket center identity = 0` and `PASS: residual bracket center identity` must no longer appear in the saved outputs.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`
- summary: Removed the tautological residual bracket center identity assertions while preserving the informational bracket prints.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:74-79`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl:68-70`

**Issue:**
The three `expect_positive` checks on the ordered branch (sympy 74-79; wl 68-70) reduce to ratios of declared-positive symbols (`Pe_req*delta_gap/(Delta0*(Delta0+delta_gap))`, `Pe_req*delta_gap/(Delta0+delta_gap)`, `Pe_req*delta_gap/Delta0`). Sign positivity follows directly from `positive=True` declarations in SymPy (sympy 47-49, 69) and the `$Assumptions` block in Mathematica (wl 41-44, 64). No threshold physics is tested.

**Required change:**
Delete the offending assertion lines but keep the upstream definitions (`Xi_fail_ordered`, `Xi_suff_ordered`, `delta_gap`, `DeltaInf_ordered`, `zeta_req_branch`) intact, since they are used in the surrounding informational structure.

SymPy — remove lines 74-79 (inclusive) in
`scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`:
```
expect_positive("Xi_suff - Xi_fail on the ordered branch", Xi_suff_ordered - Xi_fail_ordered)
expect_positive("Pe_req - Xi_fail Delta_0", sp.simplify(Pe_req - Xi_fail_ordered * Delta0))
expect_positive(
    "Xi_suff Delta_inf - Pe_req",
    sp.simplify(Xi_suff_ordered * DeltaInf_ordered - Pe_req),
)
```
Keep lines 69-73 (the `delta_gap` symbol, `DeltaInf_ordered`, `Xi_fail_ordered`, `Xi_suff_ordered`, `zeta_req_branch` definitions) as-is.

Mathematica — remove lines 68-70 (inclusive) in
`mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`:
```
expectPositive["Xi_suff - Xi_fail on the ordered branch", xiSuffOrdered - xiFailOrdered];
expectPositive["Pe_req - Xi_fail Delta_0", FullSimplify[peReq - xiFailOrdered*delta0, Assumptions -> $Assumptions]];
expectPositive["Xi_suff Delta_inf - Pe_req", FullSimplify[xiSuffOrdered*deltaInfOrdered - peReq, Assumptions -> $Assumptions]];
```
Keep lines 63-67 (the `deltaGap` assumption update, `deltaInfOrdered`, `xiFailOrdered`, `xiSuffOrdered`) as-is.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 059` and `redteam exec-mathematica 059`. The output lines `Xi_suff - Xi_fail on the ordered branch = ...`, `Pe_req - Xi_fail Delta_0 = ...`, `Xi_suff Delta_inf - Pe_req = ...` must no longer appear (and the corresponding `PASS:` lines in the Mathematica output must also disappear).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`
- summary: Removed the ordered-branch positivity assertions while leaving the threshold definitions intact.
- deviation: none

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:94-118`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl:77-101`

**Issue:**
The numerical "saturation" checks define `zeta_req_probe := A_K_probe*Omega_Pe(Pe_req_probe)^2` (sympy:102, wl:85), then solve `A_K_probe*Omega_Pe(Xi*Delta)^2 - zeta_req_probe = 0` for `Xi` seeded at `Pe_req_probe/Delta` (sympy 105-116, wl 88-97), and assert the root equals `Pe_req_probe/Delta` (sympy 117-118, wl 100-101). By construction `Xi*Delta = Pe_req_probe` solves the equation exactly, and the seed *is* that answer; the check is tautological and would pass for any choice of `Pe_req_probe`/`Delta` regardless of whether the threshold formulas `Xi_fail = Pe_req/DeltaInf` and `Xi_suff = Pe_req/Delta0` have correct physical meaning. The two additional Mathematica-only `expectPositive` lines on `xiFailSolved` and `xiSuffSolved` (wl 98-99) only check that a positive seed remains positive.

**Required change:**
Break the circularity by giving `zeta_req_probe` a literal numerical value independent of `Omega_Pe`, then solve for `Pe_star` independently, then assert `Xi_fail*DeltaInf == Pe_star` and `Xi_suff*Delta0 == Pe_star`.

SymPy — replace lines 94-118 in
`scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`
with:
```
Xi_num = sp.symbols("Xi_num", real=True)
Pe_num = sp.symbols("Pe_num", real=True)
kappa_probe = sp.Integer(2)
y_probe = sp.Integer(1)
Delta0_probe = sp.Rational(3, 5)
delta_gap_probe = sp.Rational(2, 5)
DeltaInf_probe = Delta0_probe + delta_gap_probe
A_K_probe = sp.N(A_K.subs({kappa: kappa_probe, y: y_probe}), 80)
zeta_req_probe = sp.Rational(2, 5)  # independent target, NOT constructed from Omega_Pe
Pe_star = sp.nsolve(
    sp.N(A_K_probe * Omega_Pe.subs(Pe, Pe_num) ** 2 - zeta_req_probe, 80),
    Pe_num,
    sp.Rational(1, 2),
    prec=70,
)
Pe_req_probe = Pe_star  # operator-branch threshold scale derived from the independent zeta_req
Xi_fail_expected = sp.N(Pe_req_probe / DeltaInf_probe, 80)
Xi_suff_expected = sp.N(Pe_req_probe / Delta0_probe, 80)
Xi_fail_solved = sp.nsolve(
    sp.N(A_K_probe * Omega_Pe.subs(Pe, Xi_num * DeltaInf_probe) ** 2 - zeta_req_probe, 80),
    Xi_num,
    Xi_fail_expected,
    prec=70,
)
Xi_suff_solved = sp.nsolve(
    sp.N(A_K_probe * Omega_Pe.subs(Pe, Xi_num * Delta0_probe) ** 2 - zeta_req_probe, 80),
    Xi_num,
    Xi_suff_expected,
    prec=70,
)
expect_close("Xi_fail*DeltaInf saturates at Pe_star", Xi_fail_solved * DeltaInf_probe, Pe_star)
expect_close("Xi_suff*Delta0 saturates at Pe_star", Xi_suff_solved * Delta0_probe, Pe_star)
```

Note: this still uses `Pe_star` as the nsolve seed, but `Pe_star` is now an *independently-computed* numerical inversion of the operator-branch equation against a freely-chosen `zeta_req_probe = 2/5`; the assertions then test that `Xi_fail` and `Xi_suff`, multiplied by their respective Deltas, recover that same independently-derived `Pe_star`. The threshold *relations* `Xi_fail = Pe_req/DeltaInf` and `Xi_suff = Pe_req/Delta0` are what is now under test.

Mathematica — replace lines 77-101 in
`mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`
with:
```
Clear[xiNum, peNum];
kappaProbe = 2;
yProbe = 1;
delta0Probe = 3/5;
deltaGapProbe = 2/5;
deltaInfProbe = delta0Probe + deltaGapProbe;
aKProbe = N[(kappaProbe + Pi^2/4)/(kappaProbe + yProbe^2), 80];
zetaReqProbe = 2/5;  (* independent target, NOT constructed from omegaPe *)
peStar = peNum /. FindRoot[
  N[aKProbe*(omegaPe /. pe -> peNum)^2 - zetaReqProbe, 80],
  {peNum, 1/2},
  WorkingPrecision -> 70
];
peReqProbe = peStar;
xiFailExpected = N[peReqProbe/deltaInfProbe, 80];
xiSuffExpected = N[peReqProbe/delta0Probe, 80];
xiFailSolved = xiNum /. FindRoot[
  N[aKProbe*(omegaPe /. pe -> xiNum*deltaInfProbe)^2 - zetaReqProbe, 80],
  {xiNum, xiFailExpected},
  WorkingPrecision -> 70
];
xiSuffSolved = xiNum /. FindRoot[
  N[aKProbe*(omegaPe /. pe -> xiNum*delta0Probe)^2 - zetaReqProbe, 80],
  {xiNum, xiSuffExpected},
  WorkingPrecision -> 70
];
expectApprox["Xi_fail*DeltaInf saturates at Pe_star", xiFailSolved*deltaInfProbe, peStar, 10^-40];
expectApprox["Xi_suff*Delta0 saturates at Pe_star", xiSuffSolved*delta0Probe, peStar, 10^-40];
```

The two no-physics `expectPositive` calls at wl:98-99 are subsumed into the new `expectApprox` checks and should not be re-introduced.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 059` and `redteam exec-mathematica 059`. The output must now contain `Xi_fail*DeltaInf saturates at Pe_star diff = <something < 1e-40>` and `Xi_suff*Delta0 saturates at Pe_star diff = <something < 1e-40>`, and must no longer contain `nsolve Xi_fail saturation diff` / `nsolve Xi_suff saturation diff` / `FindRoot Xi_fail saturation solver` / `FindRoot Xi_suff saturation solver` / `nsolve-style Xi_fail root stayed positive` / `nsolve-style Xi_suff root stayed positive`. Both scripts must still exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`
- summary: Replaced the circular saturation probes with an independent numerical zeta target, a solved Pe_star, and checks that Xi thresholds recover Pe_star.
- deviation: none

## F4 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl:76,106`

**Issue:**
The Mathematica script's small-Pe series check at wl:76 (`Normal[Series[omegaPe^2, {pe, 0, 1}]]`) followed by `Coefficient[..., pe, 1]` at wl:106 mirrors the SymPy `sp.series(Omega_Pe**2, Pe, 0, 2).removeO()` + `.coeff(Pe, 1)` choreography (sympy:87,91). Combined with the identical closed-form `omegaPe` (wl:72-75 vs sympy:83-86), identical probe numerics (wl:78-83 vs sympy:95-100), and identical solver seeding (wl:88-97 vs sympy:105-116), the `.wl` is a line-by-line port. The two-engine policy requires the second engine to follow a structurally different path for at least one substantive claim.

**Required change:**
Break independence on the substantive (4-Pi)/Pi linear-coefficient claim by computing the linear coefficient via differentiation/limit rather than `Series`/`Coefficient`.

Mathematica — modify lines 76 and 106 in
`mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`.

Before (wl:76):
```
omegaSqSeries = FullSimplify[Normal[Series[omegaPe^2, {pe, 0, 1}]], Assumptions -> pe > 0];
```
After (wl:76):
```
omegaSqLinearCoeff = FullSimplify[Limit[D[omegaPe^2, pe], pe -> 0], Assumptions -> pe > 0];
```

Before (wl:106):
```
expectZero["Omega^2 linear coefficient", Coefficient[omegaSqSeries, pe, 1] - (4 - Pi)/Pi];
```
After (wl:106):
```
expectZero["Omega^2 linear coefficient", omegaSqLinearCoeff - (4 - Pi)/Pi];
```

Also update the corresponding `Print` at wl:104:

Before (wl:104):
```
Print["Omega_Pe^2 small-Pe series = ", fmt[omegaSqSeries]];
```
After (wl:104):
```
Print["Omega_Pe^2 linear coefficient = ", fmt[omegaSqLinearCoeff]];
```

This swaps the Mathematica path from a Taylor-series expansion to a derivative-at-zero limit, which exercises a different internal Mathematica code path (`D[]` plus `Limit[]`) than SymPy's `series().coeff()`. The asserted value `(4 - Pi)/Pi` is unchanged; only the route to it changes.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 059`. The output must show `Omega_Pe^2 linear coefficient = 4/Pi - 1` (or equivalent simplified form of `(4 - Pi)/Pi`) and `PASS: Omega^2 linear coefficient`. The Mathematica script must still exit 0 and agree with SymPy's coefficient value.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`
- summary: Switched the Mathematica Omega_Pe squared linear coefficient check from series extraction to a derivative-at-zero limit.
- deviation: none
