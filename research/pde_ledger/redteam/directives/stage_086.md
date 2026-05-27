---
unit_id: 086
batch: III.5
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
verification_status: scripts_pass_pending_verifier
needs_user_resolution: false
---

# Codex directive — unit 086

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage086_family1_loading_ratio_window_sympy_audit.py:46-49`

**Issue:**
Lines 46–49 assert `rho_X - (1 + zeta_X) == 0` for the four carry-forward zetas, but `rho_X` is defined at lines 36–39 as `Q.subs({zeta: zeta_X, eps: 0})` and the script already verified at line 28 that `Q(zeta;eps=0) - (1+zeta) == 0` symbolically. The four numeric checks are therefore tautologically zero and do not anchor the script's `zeta` literals (lines 31–34) to the paper-stated 14-digit values. We replace each check with a direct numeric anchor of `rho_X` against the paper-stated rho-window constant.

**Required change:**

Step 1. Add a numeric helper. Immediately after the existing `expect_zero` helper (after line 17), append:

```python
def expect_close(name: str, value: sp.Expr, target: sp.Expr, tol: sp.Expr) -> None:
    diff = sp.Abs(sp.N(value - target, 40))
    print(f"{name} diff = {diff}")
    if diff > tol:
        raise AssertionError(f"{name} exceeds tolerance {tol}: diff={diff}")
```

Step 2. Replace lines 46–49 in their entirety with the following four lines, anchoring each computed rho to the paper-stated 14-digit constant from `paper/stages/stage_086.tex` lines 17–28 and notes section 2:

```python
expect_close("rho_suff^(chi) vs paper", rho_suff_chi, sp.Float("3.46622291347846", 30), sp.Float("1e-13", 30))
expect_close("rho_fail^(chi) vs paper", rho_fail_chi, sp.Float("3.46752913273870", 30), sp.Float("1e-13", 30))
expect_close("rho_suff^(J)   vs paper", rho_suff_J,   sp.Float("3.44257571477179", 30), sp.Float("1e-13", 30))
expect_close("rho_max        vs paper", rho_max,      sp.Float("3.46752922945601", 30), sp.Float("1e-13", 30))
```

Do not touch lines 31–44 or lines 51–63. Do not change `rho_X` definitions.

**Verification:**
The verifier will run `redteam exec-sympy 086` and confirm (a) exit code 0, (b) the transcript contains four new `rho_… vs paper diff = …` lines, (c) each printed diff is below `1e-13`. A negative control: perturbing any one of the four `zeta_*_chi` / `zeta_max` literals at lines 31–34 by `+1e-12` must cause exactly one of the new checks to fail.

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage086_family1_loading_ratio_window_mathematica_audit.wl:48-51`

**Issue:**
The four `expectApprox` checks at lines 63–66 compare each computed `rho` against an extended-precision target that is just `1 + zeta_literal` extended in precision. Because `qMap /. epsBlk -> 0 = 1 + zeta` is verified symbolically at line 46, the comparison reduces to `Abs[(1 + zetaSuffChi) - (1 + 2.46622291347846…)]`, which equals `Abs[zetaSuffChi - 2.46622291347846…]`. This is acceptable as an indirect anchor of `zetaSuffChi` but only at the same precision as the stored literal; in particular, mistyping `zetaSuffChi` would shift both sides and the check might still pass within tolerance. We add four explicit anchors comparing the four `zeta` literals directly to the paper-stated 14-digit values.

**Required change:**

Insert the following block immediately after line 51 (i.e., after the `zetaMaxNum = ToExpression["2.46752922945601`20"];` line, before the `rhoSuffChi = ...` line on line 53). Maintain a blank line before the inserted block for readability:

```
expectApprox["zeta_suff^(chi) vs paper", zetaSuffChi, 2.46622291347846, 10^-14];
expectApprox["zeta_fail^(chi) vs paper", zetaFailChi, 2.46752913273870, 10^-14];
expectApprox["zeta_suff^(J) vs paper",   zetaSuffJ,   2.44257571477179, 10^-14];
expectApprox["zeta_max^(F1) vs paper",   zetaMaxNum,  2.46752922945601, 10^-14];
```

Do not modify lines 48–51 themselves; do not modify the existing `expectApprox` rho checks at lines 63–66.

**Verification:**
The verifier will run `redteam exec-mathematica 086` and confirm (a) exit code 0, (b) the transcript contains four new `zeta_… vs paper diff = …` PASS lines preceding the `rho_…` numeric checks, (c) the existing rho-numeric checks still pass. Negative control: perturbing any one of `zetaSuffChi`, `zetaFailChi`, `zetaSuffJ`, `zetaMaxNum` by `+10^-13` must cause exactly one of the new `zeta_… vs paper` checks to fail with exit code 1.
