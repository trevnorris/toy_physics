---
unit_id: 082
batch: III.4
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 082

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:103-109` — the two `expect_zero(...)` calls for `Xi_F1(Theta_w) - 136900 Theta_w` and `Xi_F1(Upsilon_w=100 Theta_w) - Xi_F1(Theta_w)` were replaced with plain `print(...)` lines, each suffixed with `"  (display only)"` and preceded by the comment block from the directive.
- `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:94-102` — the matching `expectZero[...]` calls were replaced with `Print[...]` calls carrying the same `(display only)` suffix and equivalent comment.

**Assessment:**
Matches the directive verbatim. The two arithmetic-only assertions no longer enforce equality, so they can no longer mask tautological agreement with a PASS. The accompanying comment explicitly disclaims provenance of `Lambda_ell = 37` / `Upsilon_w = 100 Theta_w` from this stage. The canonical SymPy output (lines 41-42) and Mathematica output (lines 26-27) both display the suffix and show the residual `= 0` as informational text rather than an assertion. No collateral edits to the surrounding `Lambda_ell`/`Xi_F1_*` definitions.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:77-87` — inserted the comment block plus two new `expect_zero` calls: `dR_quad/dzeta_phys + 1` and `dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-)`, using `sp.diff` on `R_quad` and `zeta_req`.
- `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:71-81` — mirror block using `D[...]` with `FullSimplify` and the `$Assumptions` flag, two corresponding `expectZero[...]` calls.

**Assessment:**
The additions match the directive exactly. The new assertions are non-tautological: `R_quad = zeta_req - zeta_phys` makes `d R_quad / d zeta_phys = -1` a real symbolic differentiation result (would fail if `zeta_phys` were silently absorbed into a renamed quantity), and the second check actually compares two independently-evaluated partials (`d R_quad / d Pi_tr` evaluated after substituting `zeta_phys -> zeta_minus`, vs. `d zeta_req / d Pi_tr`) — agreement is non-trivial and exercises that `zeta_minus` is treated as independent of `Pi_tr`. The canonical outputs confirm both new lines: SymPy output 34-35 shows `dR_quad/dzeta_phys + 1 = 0` and `dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-) = 0`; Mathematica output 20-23 shows the corresponding `= 0` and `PASS:` lines.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:33-40` — the literal closed-form assignment `zetaReq = FullSimplify[(PiTr - cMix)/(cMix - epsBlk*(2*cMix - PiTr)), ...]` was replaced with a `Solve[PiTr == cMix*((1 + (1 - 2*epsBlk)*zetaSym)/(1 - epsBlk*zetaSym)), zetaSym]` call. `zetaSym` is allocated via `Unique["zetaSym"]`, the first solution branch is extracted, `ConditionalExpression[x_, _] :> x` strips the defensive wrapper per the project Mathematica idiom note, then `FullSimplify` reduces it. The duplicate `qMap = ...` line is hoisted above the `Solve` call (the original line 34 was effectively repositioned to line 36); no orphan duplicate remains.

**Assessment:**
Matches the directive (Codex took the second variant from the directive's alternative: keep `qMap` defined once before being used in `Solve`). The literal expression `(PiTr - cMix)/(cMix - epsBlk*(2*cMix - PiTr))` no longer appears in the `.wl` script outside any comment (I scanned the file). The Mathematica engine now produces `zetaReq` by inverting `qMap` via `Solve`, so the inverse-map `expectZero` is a genuine cross-engine check rather than a self-reference. Output line 5 prints `zeta_req(Pi_tr,C_mix,eps_blk) = (-cMix + PiTr)/(cMix - 2*cMix*epsBlk + epsBlk*PiTr)`, which is algebraically identical to the SymPy form and to the prior Mathematica output — confirming `Solve` recovered the correct branch and `ConditionalExpression` was successfully stripped. All five downstream assertions (lines 8, 12, 14, 17, 19) and the F2 derivative assertions (lines 21, 23) still PASS.

### F4 — hardcoded_result

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:91-94` — four-line `# TODO(provenance): ...` comment block inserted immediately before `Lambda_ell = sp.Integer(37)`.
- `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:83-87` — equivalent five-line `(* TODO(provenance): ... *)` comment block inserted immediately before `lambdaEll = 37;`.

**Assessment:**
Comment-only edits, matching the directive's wording. No numeric values changed. The TODO markers correctly avoid inventing a stage number, per the directive's guidance.

## Exec log assessment

**SymPy:** exit=n/a. The captured exec log `redteam/exec_logs/stage_082_sympy.log` is absent. However, the canonical SymPy output file `scripts/output/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.txt` has mtime `May 25 00:24` (newer than the script's `May 25 00:22`), confirming it was regenerated post-fix. Notable lines from the regenerated output:

- `inverse map zeta_req(C_mix*Q(zeta)) - zeta = 0`
- `dR_quad/dzeta_phys + 1 = 0`
- `dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-) = 0`
- `Xi_F1(Theta_w) - 136900 Theta_w = 0  (display only)`

All `expect_zero` lines show `= 0`; the F1 demoted lines now carry the `(display only)` suffix.

**Mathematica:** exit=n/a. The captured exec log `redteam/exec_logs/stage_082_mathematica.log` is absent. The canonical Mathematica output `mathematica/output/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.txt` has mtime `May 25 00:24` (newer than the `.wl` mtime `May 25 00:22`), confirming refresh post-fix. Notable lines:

- `PASS: inverse map zeta_req(C_mix*Q(zeta)) - zeta` (line 8) — confirms the Solve-derived `zetaReq` correctly inverts `qMap`.
- `PASS: dR_quad/dzeta_phys + 1` (line 21) and `PASS: dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-)` (line 23) — confirms the F2 additions.
- `Xi_F1(Theta_w) - 136900 Theta_w = 0  (display only)` (line 26) — confirms F1 demotion; no `PASS:` line follows.

The absence of `FAIL:` lines and the `Exit[0]` at the script's end (the script would have exited 1 via `fail[...]` on any failure) imply a clean run.

**Output freshness:** confirmed. Both `.txt` outputs (`May 25 00:24`) are newer than their corresponding scripts (both `May 25 00:22`).

## Material-change assessment

`material_change`: false.

The printed symbolic forms of `zeta_req`, `Q`, `Pi_suff`, `Pi_fail`, `R_quad`, `Xi_F1 from Upsilon_w`, and `Xi_F1 from Theta_w` are unchanged from the pre-fix outputs (still `(-cMix + PiTr)/(cMix - 2*cMix*epsBlk + epsBlk*PiTr)` for `zetaReq`, etc.). F1's demotion is a provenance/assertion-status change (the residual values are still printed as `= 0`). F2 added new printed lines but did not alter any previously-printed value. F3 changed the derivation route in Mathematica but produced the same closed form for `zetaReq`. F4 is comment-only.

## Side observations (non-blocking)

- The script docstring (lines 2-11) and banner (line 32: `STAGE 65 — MASTER QUADRUPOLE RESIDUAL`) still reference "STAGE 65" while the file/audit unit is `082`. The audit report flagged this in passing but did not raise it as a finding. Out of scope for this verification.
- The Mathematica banner reads `STAGE 065` while the SymPy banner reads `STAGE 65`. Minor cosmetic inconsistency; not in any finding.
- The F4 TODOs leave the provenance unfilled (as the directive permitted). A future audit pass may want to chase down the upstream stage that establishes `Lambda_ell = 37` and the `Upsilon_w = 100 Theta_w` convention.

## Verdict justification

All four findings resolved as specified by the directive, with no deviation reported by Codex and none observable in the diff. The refreshed canonical outputs reflect the expected post-fix state: F1's two checks are now display-only with `(display only)` suffixes, F2's two new derivative assertions print and PASS in both engines, F3's Mathematica `zetaReq` is derived via `Solve` (literal closed form gone from the script outside comments) while still producing the algebraically-correct inverse, and F4's TODO comments are in place. No regressions; no new symbolic values affected; downstream stages do not need to re-audit on the basis of changed numerics.

stage 082: verified
