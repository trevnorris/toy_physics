---
unit_id: 063
batch: III.3
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T20:15:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 063

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:102-122` — added three new `expect_zero` calls. A local solver symbol `gphi_sq = sp.symbols("gphi_sq_solve", ...)` is introduced; `G_micro` is rewritten in terms of `gphi_sq` and `sp.solve(G_micro_gphi - G_fail, gphi_sq)` (and the analogous suff call) is used to *derive* the threshold, then compared against the hand-rearranged `g_fail_sq` / `g_suff_sq`. A `Cauchy saturation` check substitutes `Osp**2 -> Nss*Npp` into `G_micro` and asserts equality with `G_max`.
- `mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl:85-105` — analogous block using `Reduce[... && gphiSq > 0, gphiSq, Reals]` and a Cauchy-saturation `expectZero`.

**Assessment:**
The new sympy assertions are non-tautological. `sp.solve(G_micro_gphi - G_fail, gphi_sq)` produces the root of the equation `rho_star*gphi_sq*Osp^2/(m*cs_star_sq*K_X*N_ss) == G_fail`. The solver re-derives `g_phi^2 = m*cs_star_sq*K_X*N_ss*G_fail/(rho_star*Osp^2)` directly from the gain definition; comparing against the separately-written `g_fail_sq` literal is a real consistency check — if `G_micro` is later altered (e.g. wrong factor of `N_ss`, missing `K_X`), the solver result will diverge from the hand-rearranged form and the residual will be nonzero. Similarly the Cauchy-saturation check (`G_micro|_{Osp^2 = N_ss N_pp} == G_max`) is non-trivial: if `G_max` had been declared with `N_ss` rather than `N_pp` (a plausible bug since both are written in adjacent expressions), the residual would not collapse to zero. The unique-root precondition `assert len(sol_fail) == 1 and len(sol_suff) == 1` guards against the solver returning multiple branches (e.g. via a future change making the equation quadratic). The sympy exec log shows `g_fail^2 from solve(G_micro=G_fail) matches hand-rearranged form = 0`, `g_suff^2 from solve(G_micro=G_suff) matches hand-rearranged form = 0`, and `G_max = G_micro at Cauchy saturation (O_sp^2 = N_ss N_pp) = 0`. No collateral edits beyond the inserts. The original line 102 (`print("\nAll Stage 46 symbolic checks passed.")`) is preserved at the new line 124.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/.../stage063_..._mathematica_audit.wl:85-99` — uses `Reduce[(gMicro /. gPhi^2 -> gphiSq) == gFail && gphiSq > 0, gphiSq, Reals]` rather than `Solve`. The Cases-based root extraction (`Cases[reduceFail, HoldPattern[gphiSq == rhs_] :> rhs, Infinity][[1]]`) replaces the directive's `ToRules` approach due to the documented `ReplaceAll::argt` issue.

**Assessment:**
`Reduce` is structurally different from `sp.solve`: it enumerates branches under explicit domain constraints (`gphiSq > 0`, `Reals`) rather than algebraically isolating the unknown. Since `gphiSq > 0` is imposed, the result is `gphiSq == m csStarSq kX nSS gFail/(rhoStar oSP^2)` — the unique positive root. The `Cases` extraction with `HoldPattern[gphiSq == rhs_]` deviation is acceptable and accomplishes the same goal as `ToRules` (recovering the right-hand-side expression for `gphiSq`); the documented `ReplaceAll::argt` error is a known artifact when `ToRules` returns a `Sequence` against a non-list LHS. Crucially, the Mathematica path now has independent provenance for the threshold expression. The Cauchy-saturation check is left in the same algebraic-substitution form as sympy, which is acceptable because the upstream `gMicro` is independently formed via `FullSimplify` under `$Assumptions`. Exec log shows `g_fail^2 from Reduce[...] matches hand-rearranged form = 0` and `PASS:` lines for both Reduce-derived assertions and the Cauchy check. No `Solve` is used anywhere in the new block (verified via grep).

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
F3 is a documentation-anchor finding only. The F1 inserts include the required comments. In sympy line 117: `# Cauchy-Schwarz: O_sp^2 <= N_ss N_pp, with equality at perfect alignment.` In Mathematica line 101: `(* Cauchy-Schwarz: oSP^2 <= nSS nPP, with equality at perfect alignment (c2 = 1). *)`. Both comments sit directly above the saturating-substitution `expect_zero` / `expectZero` calls.

**Assessment:**
The literal phrase "Cauchy-Schwarz" appears adjacent to the saturation block in both files, and each references the bound in its file's variable convention (`O_sp^2 <= N_ss N_pp` / `oSP^2 <= nSS nPP`). The directive explicitly states that F3 is covered by the F1 comments when applied verbatim, which is the case here. No additional `$Assumptions` change was required (and would have been risky, per the original report).

## Exec log assessment

**SymPy:** exit=0 (inferred — saved output file is fresh and ends with `All Stage 46 symbolic checks passed.`). Notable lines:
- `g_fail^2 from solve(G_micro=G_fail) matches hand-rearranged form = 0`
- `g_suff^2 from solve(G_micro=G_suff) matches hand-rearranged form = 0`
- `G_max = G_micro at Cauchy saturation (O_sp^2 = N_ss N_pp) = 0`

**Mathematica:** exit=0 (inferred — saved output file is fresh, ends with `Stage 063 Mathematica audit passed.`, and contains `PASS:` lines for every assertion). Notable lines:
- `g_fail^2 from Reduce[gMicro==gFail, gphiSq>0] matches hand-rearranged form = 0` + `PASS:`
- `g_suff^2 from Reduce[gMicro==gSuff, gphiSq>0] matches hand-rearranged form = 0` + `PASS:`
- `G_max = G_micro at Cauchy saturation (oSP^2 = nSS nPP) = 0` + `PASS:`

Note: dedicated `stage_063_sympy.log` / `stage_063_mathematica.log` files are not present under `redteam/exec_logs/`, only `stage_063_diff.patch`. The saved `.txt` outputs in `scripts/output/` and `mathematica/output/` serve as the exec evidence and are consistent with successful exit-0 runs (no AssertionError, no `FAIL:` line, terminal `passed.` print reached).

**Output freshness:**
- sympy script mtime: 1779500494; sympy output mtime: 1779500638 (output newer by 144s — fresh).
- mathematica script mtime: 1779500554; mathematica output mtime: 1779500643 (output newer by 89s — fresh).

## Material-change assessment

`material_change`: false.

The edits add new assertions that *confirm* the existing threshold formulas (`g_fail_sq`, `g_suff_sq`, `G_max`) are internally consistent with the gain definition `G_micro` and with the Cauchy-Schwarz alignment bound. No symbolic outputs that downstream units consume are altered — `g_fail_sq`, `g_suff_sq`, `C_fail_sq`, `C_suff_sq`, `G_max`, and the kappa-substituted forms are unchanged. Downstream units consume these expressions as written in the script source, which the new checks merely vouch for; if the new checks had failed, the diagnosis would propagate upstream, not downstream.

## Side observations (non-blocking)

- The directive notes that the placeholder line `reduceFail = Reduce[gMicro /. gPhi^2 -> gphiSq, gphiSq] // Quiet;` may be dropped if the editor allows; Codex correctly omitted it (the final `.wl` only contains the second, constrained `Reduce` assignment).
- The unique-root precondition is asserted in sympy via Python `assert`, which is appropriate (raises before reaching the `expect_zero`). Mathematica relies on `Cases[..., ...][[1]]` which will throw a `Part::partw` if the Cases list is empty; this is acceptable as a "first solution must exist" gate but is less explicit than the sympy precondition. Not blocking.
- The orchestrator did not save dedicated `stage_063_sympy.log` / `stage_063_mathematica.log` files; only the diff was captured. The saved `.txt` outputs are sufficient evidence, but if exec-log capture is supposed to be a separate artifact, the harness may have skipped it for this unit.

## Verdict justification

All three findings are resolved. F1's tautology critique is closed by the new `sp.solve` and Cauchy-saturation `expect_zero` calls on the sympy side and by the analogous `Reduce` / Cauchy-saturation calls on the Mathematica side — these assertions exercise the derivation rather than restating it, and they would fail under plausible bugs (wrong factor in `G_micro`, `N_ss`-vs-`N_pp` swap in `G_max`). F2's transliteration critique is closed by the Mathematica engine now using `Reduce[... && gphiSq > 0, gphiSq, Reals]` — a genuinely different solver path. F3's documentation-anchor requirement is satisfied by the inline `Cauchy-Schwarz` comments in both files. Saved outputs are fresh, contain the new expected zero residuals, and both scripts reach their terminal "passed" prints with no `FAIL:` lines.
