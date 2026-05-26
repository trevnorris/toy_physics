---
unit_id: 013
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T22:15:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: true
---

# Verification — unit 013

## Per-finding outcomes

### F1 — paper_misalignment (paper_missing_script_claim)

**Classification:** resolved

**What changed:**
Per the harness note, F1 was user-resolved as direction (b): trim the scripts.
The diff at `redteam/exec_logs/stage_013_diff.patch` shows Codex removed:
- SymPy: `D2, D4, N0, Ptarget` symbol declarations (line 42), `K1`, `H_even`, `deltaP2`, `deltaP4`, `deltaP2_der`, `deltaP4_der` definitions (former lines 101-128), the `d(delta P2)/dGprime` and `deltaP4`/`G_W` prime assertions, and the entire `qd_matrix`/`sh_matrix` mechanism sieve (former lines 133-158).
- Mathematica: the same `D0, D2, D4, N0` extras from `Clear`/`$Assumptions` (lines 28, 35), and the `K1`, `HEven`, `deltaP2`, `deltaP4`, `deltaP2Der`, `deltaP4Der` definitions plus the M5 `d(delta P2)/dGprime`, `delta P4 Gprime dependence`, and M6 sieve determinant/solve assertions (former lines 119-166).
The closing SymPy print statement was updated to drop the "bottleneck dependencies, G_W transport entry, and mechanism sieve" phrasing (now "Taylor coefficient maps, and Xi_load Pprime dependence"). The Taylor primitive map (`z_0, z_2, z_4, n_0, n_2, n_4`) and the new `Xi_load` paper closed-form check are retained in both engines.

**Assessment:**
The trim is targeted and matches direction (b) from the directive's `## Resolve before fix_loop` block. The script now verifies exactly the four equation-numbered paper deliverables (`z_0, n_0, z_2, z_4`) and the narrative `Xi_load`, with `n_2, n_4` retained as natural chain-rule extensions (which the directive option (b) said "trim... or move to stage 014"; keeping them is conservative and the `n_2, n_4` literals are still independently cross-checked by Mathematica `Series`). No surviving K_1/H_even/deltaP/sieve content. The narrative direction the user picked is fully reflected in the diff.

### F2 — insufficient_verification (Xi closed-form anchor)

**Classification:** resolved

**What changed:**
SymPy `moving_throat_pde_stage013_..._sympy_audit.py:101-107` defines `z0_form`, `n0_form`, `N0_form = P**2/Delta**2`, builds `Xi_paper`, and asserts `assert_zero("Xi matches paper closed form n0/N0 + z0/D0", Xi - Xi_paper)`.
Mathematica `..._mathematica_audit.wl:119-128` mirrors this with `z0Form`, `n0Form`, `N0Form`, `XiPaper`, and `assertZero["M5 Xi matches paper closed form n0/N0 + z0/D0", Xi - XiPaper]`.

**Assessment:**
The insertion matches the directive's required-change block verbatim modulo trivial spacing. The new assertion (a) makes `N_0 = P^2/Delta^2` explicit, (b) tests the full `Xi = n_0/N_0 + z_0/D_0` structure rather than one partial. Both engines use the paper-side polynomial forms (no fresh hardcoded literals). Self-test: dropping `z_0/D_0` from the script-side `Xi` would now produce a nonzero residual `q1/(D0 Delta) - Q d1/(D0 Delta^2)` and trigger `assert_zero`. Non-tautological. The Mathematica output transcript line 11 confirms `OK M5 Xi matches paper closed form n0/N0 + z0/D0 residual = 0`. The pre-existing `dXi/dPprime` check is retained alongside.

### F3 — insufficient_verification (K_1/H_even coefficient guards)

**Classification:** blocked_legitimate

**What changed:**
Codex appended a `## Blocked: F3` block to the directive with the correct reason: "The current stage 013 SymPy and Mathematica scripts no longer define `K1`, `H_even`, `deltaP2`, or `deltaP4`, so the requested coefficient-anchor blocks have no valid insertion point or in-scope symbols."

**Assessment:**
This block is mechanically correct. After the F1 trim (direction b), `K1` and `H_even` are no longer defined anywhere in either stage 013 script, so the F3 patches would reference undefined symbols and the assertions would be vacuous (or fail with `NameError`/`Symbol`-out-of-context). The harness note in the verifier prompt explicitly confirms this block reason is "correct, not a defect" — F3's load-bearing concern (coefficient anchors for `K_1`, `H_even`) belongs to whichever downstream stage owns the gate definitions (stage 014, per the original report's own observation that the literal weights `-1/9, 2/3, -1/27` match stage 014 paper definitions). The verifier hands F3 forward as a follow-on concern for stage 014's audit, not as a stage 013 defect. Per the verifier-prompt classification grid this falls under `blocked_legitimate`; the per-finding classification is treated as resolved-by-trim for the overall verdict rollup because the underlying issue has been migrated out of this unit's scope, not deferred.

## Exec log assessment

**SymPy:** exit=0 (per harness note; scripts run directly under single-seat Mathematica rule; the `redteam/exec_logs/stage_013_sympy.log` is the pre-trim stale stub from May 21 and is not used). Saved output `scripts/output/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.txt`:

```
STEP 11 PROJECTED MAXWELL MOUTH-TAYLOR MASTER AUDIT
Checked one-sided Taylor projection, Taylor coefficient maps, and Xi_load Pprime dependence.
STATUS: PASS
```

The SymPy transcript is bare (no per-assertion residual lines — `assert_zero` only raises on failure), but the "STATUS: PASS" line only prints after all assertions return, and the new `Xi matches paper closed form n0/N0 + z0/D0` assertion is in source.

**Mathematica:** exit=0. Saved output `mathematica/output/..._mathematica_audit.txt` shows 11 `OK ... residual = 0` lines including the three M3 z-coefficient checks, three M4 n-coefficient checks, the two M1/M2 Taylor lemma checks, the M2 second-moment integral, the new `M5 Xi matches paper closed form n0/N0 + z0/D0` line, and the existing `M5 dXi/dPprime` line. Final line: `STAGE 013 MATHEMATICA AUDIT: PASS`.

**Output freshness:** Both `.txt` outputs have mtime 21:54 (May 25), strictly newer than the script mtimes at 21:52 (May 25) — outputs are post-fix.

## Material-change assessment

`material_change`: true.

The trim removed substantive verification content (`K1`, `H_even`, `deltaP2`, `deltaP4`, mechanism sieve) from stage 013. Stage 014 owns `K_1` and `H_even` per paper; the orchestrator should verify that stage 014's audit covers (a) the literal coefficient weights `-1/9, 2/3, -1/27` (the F3 coefficient anchor concern, now homeless), (b) the mechanism-sieve trivial-solution structure, and (c) the `deltaP_2, deltaP_4` formulas if they belong there. Downstream units > 013 should be flagged `upstream_stale: true` so the user can decide narrow vs broad re-audit; in particular stage 014 deserves explicit attention.

The four equation-numbered deliverables (`z_0, n_0, z_2, z_4`) and `Xi_load` are unchanged structurally; their downstream consumers are not affected by the trim.

## Side observations (non-blocking)

- The SymPy `Xi_paper` simplifies symbolically to the same expression as `Xi` (`n_0/N_0 = 2 p1/P - 2 d1/Delta` exactly when `N_0 = P^2/Delta^2`); the new assertion is non-tautological in the sense that flipping a sign or dropping a summand in `Xi` would break it, but it is structurally a definition round-trip rather than a derivation, which is exactly what F2 requested.
- The diff also removed `D2`, `D4`, `N0`, `Ptarget` symbol declarations that were only consumed by the removed `deltaP` block. Clean collateral, no orphan references remain.
- The SymPy assertion `assert_nonzero("Xi should depend on Pprime", sp.diff(Xi, Px))` (line 110) is redundant with the kept `dXi/dPprime` equality check (line 108), but harmless and pre-existing.

## Verdict justification

All three findings are addressed correctly: F1 resolved by Codex trim per user direction (b); F2 resolved by the paper-closed-form `Xi` round-trip insertion in both engines, with Mathematica transcript showing the new `OK` line; F3 legitimately blocked because the target symbols no longer exist after the F1 trim, with the underlying concern correctly redirected to stage 014's scope. Both engines exit 0 and outputs are fresh. No regressions in the diff. Verdict: `verified`. `material_change: true` because downstream stage 014 inherits the F3 coefficient-anchor concern.
