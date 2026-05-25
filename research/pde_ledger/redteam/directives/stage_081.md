---
unit_id: 081
batch: III.4
created_at: 2026-05-22T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-23T05:39:00Z
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 081

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — hardcoded_result

**Target:** `mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:41`

**Issue:** The Mathematica script derives `piOfZeta` via `Solve` on the premise, but then defines `qq` as a hand-written closed form rather than from `piOfZeta/cMix`. All downstream assertions exercise the hand-written form against itself, so the Solve-based inversion is never assertively verified.

**Required change:**
At line 41, change the bare hardcoded definition of `qq` into a derivation from `piOfZeta`, then add an immediate `expectZero` that ties the derived `qq` back to the same closed form the script previously asserted.

Replace:
```
qq = FullSimplify[(1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta), Assumptions -> $Assumptions];
```

with:
```
qq = FullSimplify[piOfZeta / cMix, Assumptions -> $Assumptions];
expectZero["Q matches closed form",
  qq - (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)];
```

Do not remove the `piOfZeta` definition on line 40; it is now consumed by the new `qq` definition.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 081` and confirm the new `PASS: Q matches closed form` line appears AND the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl`
- summary: Derived `qq` from `piOfZeta / cMix` and added the closed-form residual check.
- deviation: none

## F2 — tautological_check

**Target:** `mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:42,48`

**Issue:** The `dqq` definition (line 42) and its `expectZero` assertion (line 48) are pure calculus of the hardcoded `qq` and verify no physical content. With F1 applied, `qq` is now derived from `piOfZeta`, so the `dQ/dzeta` calculus identity adds no further coverage.

**Required change:**
Delete line 42 entirely:
```
dqq = FullSimplify[D[qq, zeta], Assumptions -> $Assumptions];
```

Delete line 48 entirely:
```
expectZero["dQ/dzeta exact formula", dqq - (1 - epsBlk)/(1 - epsBlk*zeta)^2];
```

If the printer that references `dqq` (none in current script — line 42 was only consumed by line 48), no further edits are needed.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 081` and confirm the output transcript no longer contains `dQ/dzeta exact formula` AND the script exits 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl`
- summary: Removed the `dqq` derivative helper and its tautological `dQ/dzeta exact formula` assertion.
- deviation: none

## F3 — tautological_check

**Target:** `mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:68-72`

**Issue:** Each `expectApprox` compares `qq /. zeta -> zeta_* /. epsBlk -> 0` (which is `1 + zeta_*` by the closed form) against a literal target that *is* `1 + zeta_*` to typed precision. The check is tautological. After F1, the closed form is verified, so the right assertion is the symbolic identity Q(zeta; 0) = 1 + zeta evaluated at each threshold, not a numeric match against a self-derived constant.

**Required change:**
Replace lines 68-72:
```
expectApprox["Pi_suff^(chi)/C_mix at eps=0", N[piSuffChiOverC /. epsBlk -> 0, 40], ToExpression["3.4662229134784601214`25"], 10^-14];
expectApprox["Pi_fail^(chi)/C_mix at eps=0", N[piFailChiOverC /. epsBlk -> 0, 40], ToExpression["3.4675291327386998930`25"], 10^-14];
expectApprox["Pi_suff^(J)/C_mix at eps=0", N[piSuffJOverC /. epsBlk -> 0, 40], ToExpression["3.4425757147717899187`25"], 10^-14];
expectApprox["Pi_fail^(J)/C_mix at eps=0", N[piFailJOverC /. epsBlk -> 0, 40], ToExpression["3.4675273685505798582`25"], 10^-14];
expectApprox["Pi_max^(F1)/C_mix at eps=0", N[piMaxOverC /. epsBlk -> 0, 40], ToExpression["3.4675292294560100537`25"], 10^-14];
```

with:
```
expectApprox["Pi_suff^(chi)/C_mix at eps=0 matches 1+zeta", N[(piSuffChiOverC - (1 + zetaSuffChi1)) /. epsBlk -> 0, 40], 0, 10^-14];
expectApprox["Pi_fail^(chi)/C_mix at eps=0 matches 1+zeta", N[(piFailChiOverC - (1 + zetaFailChi1)) /. epsBlk -> 0, 40], 0, 10^-14];
expectApprox["Pi_suff^(J)/C_mix at eps=0 matches 1+zeta", N[(piSuffJOverC - (1 + zetaSuffJ1)) /. epsBlk -> 0, 40], 0, 10^-14];
expectApprox["Pi_fail^(J)/C_mix at eps=0 matches 1+zeta", N[(piFailJOverC - (1 + zetaFailJ1)) /. epsBlk -> 0, 40], 0, 10^-14];
expectApprox["Pi_max^(F1)/C_mix at eps=0 matches 1+zeta", N[(piMaxOverC - (1 + zetaMaxF1)) /. epsBlk -> 0, 40], 0, 10^-14];
```

This converts each comparison from "value vs hardcoded value" into "value minus its expected functional form vs 0", exercising `qq`'s structure rather than recomputing the same constant on both sides.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 081` and confirm the five `... matches 1+zeta` PASS lines appear AND the script exits 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl`
- summary: Replaced the five epsilon-zero numeric literal checks with residual checks against the `1 + zeta` functional form.
- deviation: none

## F4 — tautological_check

**Target:** `mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:74-76`

**Issue:** The numeric blocking-ceiling check compares `epsCeiling = N[1/zetaMaxF1, 40]` against its own decimal expansion `0.40526368971137149977`. The assertion `1/zetaMaxF1 ≈ digits(1/zetaMaxF1)` cannot fail in a meaningful way; it should instead verify the reciprocal relationship symbolically.

**Required change:**
Replace lines 74-76:
```
epsCeiling = N[1/zetaMaxF1, 40];
Print["Blocking ceiling eps_blk < ", fmt[epsCeiling]];
expectApprox["blocking ceiling numeric check", epsCeiling, ToExpression["0.40526368971137149977`25"], 10^-14];
```

with:
```
epsCeiling = N[1/zetaMaxF1, 40];
Print["Blocking ceiling eps_blk < ", fmt[epsCeiling]];
expectApprox["blocking ceiling reciprocal", N[epsCeiling*zetaMaxF1 - 1, 40], 0, 10^-14];
```

This checks `epsCeiling * zetaMaxF1 == 1` (the defining identity of the reciprocal) instead of comparing decimal expansions.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 081` and confirm `PASS: blocking ceiling reciprocal` appears AND the script exits 0.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl`
- summary: Replaced the blocking ceiling decimal comparison with a reciprocal identity residual check.
- deviation: none
