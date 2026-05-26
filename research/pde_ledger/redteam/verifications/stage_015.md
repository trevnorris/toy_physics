---
unit_id: 015
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: true
---

# Verification — unit 015

## Per-finding outcomes

### F1 — paper_misalignment

**Classification:** resolved

**What changed:**
Per orchestrator note, the user selected Q3 direction (b) (TRIM). Codex removed the wall-only / Y20 / grouped trace blocks from both engines:

- `scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`: removed the `gaunt` import, the `real_y20_square_ratio` helper, the `grouped_trace_anomaly` helper, and the wall-only / Y20 / grouped trace block (former lines ~103-208). Updated the final `print` summary to drop the wall-only/grouped wording (lines 110-113). Script is now 117 lines.
- `mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl`: removed M4-M9 (wall-only K1/H_even specializations, Gaussian overlap closed forms, perturbed-solve, Y20 ratios, grouped trace). Wall-only block at former lines 104-196 fully excised.

**Assessment:**
The trim correctly retains the K_eta exact-quadratic-recovery and IBP/boundary checks (the only items the stage card actually exports) and removes the orphaned `step_13_*_notes.md`-derived material. Stage 017 covers the trimmed content per orchestrator context. The sympy docstring still references `step_13_parent_throat_action_master_notes.md` at line 2 — minor cosmetic carry-over but the script body no longer pretends to verify those blocks. Not a defect for verification purposes since the orchestrator's resolution explicitly authorized direction (b) with the script trim.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
Both engines gained an asymmetric concrete IBP probe immediately after the original Gaussian baseline:

- sympy lines 56-79: `A_concrete_asym = w_ibp * exp(-w_ibp**2)` (odd), `eta_concrete_asym = exp(-w_ibp**2 / 2)` (even); asserts cross and bulk are nonzero individually, boundary discharge is zero, and the IBP identity `cross - (boundary + bulk) = 0` holds.
- mathematica lines 83-100: parallel `aConcreteAsym = w*Exp[-w^2]`, `etaConcreteAsym = Exp[-w^2/2]`, same four `expectNonzero` / `expectZero` checks.

**Assessment:**
Mathematica output confirms the residuals are non-tautological: `M2 asymmetric IBP cross nontrivial residual = Sqrt[Pi/2]/4` and `M2 asymmetric IBP bulk nontrivial residual = Sqrt[Pi/2]/4` — both individually `sqrt(pi/2)/4 ≠ 0` (matching the directive's predicted value), with their difference vanishing (true IBP cancellation, not 0=0). The original Gaussian baseline was preserved as a parity check, per directive. Matches the directive's required-change block character-for-character.

### F3 — tautological_check

**Classification:** blocked_legitimate (resolved by F1 trim)

**What changed:**
Codex correctly skipped F3 in iter1 per directive. The trim from F1 (direction b) physically removed the wall-only K1/H_even specialization asserts (former sympy 126-127 and mathematica 112-113) along with the rest of the wall-only block.

**Assessment:**
The tautological assertions no longer exist. The blocked classification is legitimate — F3 was contingent on F1's resolution direction; F1's TRIM outcome dissolves F3 mechanically (nothing left to assert against). This is `blocked_legitimate` in the prompt's taxonomy, but the underlying defect is removed, so for rollup purposes it counts as resolved.

### F4 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The K_eta portion of F4 was applied via iter2 (iter1 introduced a `Dt[..., Constants -> ...]` bug producing nonzero residual; iter2 reworked to use ordinary `D` with an explicit temporary `twR[w]` profile carrying the w-dependence of `Tw_R0`). The new M3 block (mathematica lines 102-145) computes K_eta by:

1. Defining `LDensity[R, Rt, Rw, gO, w]` with R/Rt/Rw/gO as independent slot variables.
2. Applying the Euler-Lagrange operator symbolically: `dLdRSlot - D[dLdRtSlot, t] - D[dLdRwSlot, w]` with the on-shell field substitution `rSlot -> R0[w] + eps*eta`, etc., taken AFTER the slot derivatives.
3. Reading off the O(eps) coefficient, then the eta-coefficient, then identifying it with `-K_eta`.
4. Collapsing `R0p * twR'[w] -> dTwRR0p - TwR0*R0pp` (the IBP product-rule rewrite, an identification of bookkeeping symbols, not a derivation shortcut).

The wall-only portion of F4 correctly skipped because the wall-only block was trimmed by F1.

**Assessment — is the iter2 EL derivation actually independent of SymPy's Series approach, or just retransliterated?**

It IS substantively independent. The SymPy script (lines 93-104) builds the quadratic Lagrangian `L2_raw = diff(L, eps, 2)|_{eps=0}/2`, peels the cross term `-TwR0*R0p*eta*eta_w`, IBPs it manually to `dTwRR0p*eta^2/2`, then compares the resulting `L2_after_ibp_derived` against a hand-written `canonical_L2 = mu0*eta_t^2/2 - Tw0*eta_w^2/2 - TO0*grad2/2 - K_eta*eta^2/2`. The Mathematica iter2 path does NOT touch a quadratic Lagrangian, does NOT do a manual IBP peel, does NOT compare against a hand-written canonical L2. Instead it: applies the EL operator to the FULL nonlinear `LDensity`, linearizes the resulting equation in eps, and extracts the mass coefficient from the field equation. The two paths arrive at K_eta from genuinely different intermediate representations (Lagrangian quadratic form vs. linearized field equation). The collapse rule `R0p * twR'[w] -> dTwRR0p - TwR0*R0pp` is just symbol identification (`dTwRR0p ≡ d/dw(TwR0 R0')` by definition), not an algebraic shortcut. A search-and-replace `s/sp.diff/D/`, `s/sp.expand/Expand/` transliteration would not produce this structure. The mutation guard returns `-2*dTwRR0p`, the correct expected diagnostic. The output line `M3 K_eta via EL linearization matches IBP form residual = 0` confirms agreement.

## Exec log assessment

**SymPy:** exit_code per saved output = PASS (canonical `scripts/output/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.txt` shows `STATUS: PASS`). The orchestrator-captured `redteam/exec_logs/stage_015_sympy.log` is stale (dated 2026-05-21, references removed wall-only print lines) — using canonical output per orchestrator instruction. The saved txt (mtime 2026-05-25 22:02) post-dates the script (mtime 2026-05-25 22:00), confirming freshness. Notable lines:
```
STEP 13 PARENT THROAT ACTION MASTER AUDIT
Checked promoted-action quadratic limit, concrete boundary discharge, and K_eta formula.
Boundary operator nonzero sanity check = PASS
STATUS: PASS
```

**Mathematica:** exit_code per saved output = PASS (canonical `mathematica/output/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.txt` shows `STATUS: PASS` on the final line with all per-check PASS lines and residuals). No stage_015_mathematica.log file exists in `exec_logs/`; using canonical output. Saved txt mtime 2026-05-25 22:06 > script mtime 22:05. Notable lines:
```
M2 asymmetric IBP cross nontrivial residual = Sqrt[Pi/2]/4
M2 asymmetric IBP bulk nontrivial residual = Sqrt[Pi/2]/4
M2 asymmetric IBP cross equals bulk residual = 0
M3 K_eta via EL linearization matches IBP form residual = 0
M3 K_eta via EL dTwRR0p sign mutation residual = -2*dTwRR0p
STATUS: PASS
```

**Output freshness:** confirmed (both sympy and mathematica output `.txt` files have mtimes post-dating their respective scripts, by ~2 and ~1 minutes respectively, indicating fresh regeneration after the iter2 fix).

## Material-change assessment

`material_change`: true.

Stage 015 lost ~half of its prior assertion count (wall-only K1/H_even gates, Jacobian determinant 1/27, Gaussian overlap closed forms, perturbed-solve diagnostics, real-Y20 overlap ratios, grouped trace identities). Per orchestrator context, that content was relocated to / now lives in stage 017. Downstream units that previously cited stage 015 assertions for wall-only K1/H_even, Y20 ratios, or grouped trace identities should be re-validated against stage 017 instead. Specific concern: any downstream stage that imports the `wall_only_specialization` block, `K1_wall`/`H_even_wall` symbols, the `delta_TO` 6→5 mutation guard, the `2*eps/3` determinant perturbation, or the `xbar = x0` / `bx = 3*ax` grouped identities now needs to look at stage 017 (which the orchestrator confirms already verifies the trimmed content). The orchestrator will mark all units > 015 as `upstream_stale: true`; recommend narrow re-audit of stages that referenced wall-only / Y20 / grouped material rather than a broad sweep.

## Side observations (non-blocking)

- The sympy docstring at line 2 still reads `"""Master-note audit for step_13_parent_throat_action_master_notes.md."""` — directive Q3=b mentions updating the docstring "to a stage 015 description that does not reference `step_13_*_notes.md`". Not a math defect; cosmetic and does not affect any assertion. Flagging only for cleanup.
- Mathematica `M1 mutated IBP boundary sign residual` prints a long expression instead of simplifying to a recognizable scalar — this is fine (the `expectNonzero` correctly returns PASS), just less polished display than M2's `Sqrt[Pi/2]/4`.

## Verdict justification

All four findings are addressed. F1 is resolved by the user-authorized trim to direction (b); the wall-only / Y20 / grouped blocks are gone from both engines and (per orchestrator) covered downstream in stage 017. F2 is resolved by the asymmetric IBP probe in both engines, with the predicted `Sqrt[Pi/2]/4` residual confirming the new check is non-trivial. F3 is mechanically dissolved by the F1 trim — the tautological asserts no longer exist to be tautological. F4's K_eta portion is resolved by an Euler-Lagrange linearization path that is genuinely distinct from the SymPy series-coefficient + IBP-peel choreography (different intermediate representations, no hand-written `canonical_L2` to match against); the wall-only portion of F4 is dissolved by the F1 trim. Both engines exit 0, outputs are fresh, and the iter2 fix correctly replaced the iter1 `Dt[..., Constants -> ...]` bug with explicit ordinary `D` plus a temporary `twR[w]` profile. Verdict: `verified`. `material_change: true` because the trim removes assertions that downstream units may have implicitly depended on (now covered by stage 017).
