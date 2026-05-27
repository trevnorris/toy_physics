---
unit_id: 086
batch: III.5
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 086

## Per-finding outcomes

### F1 — tautological_check (sympy)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage086_family1_loading_ratio_window_sympy_audit.py:19-23` — new `expect_close(name, value, target, tol)` helper added immediately after `expect_zero`, computing `sp.Abs(sp.N(value - target, 40))`, printing the diff, and raising `AssertionError` when `diff > tol`. Matches the directive's Step 1 verbatim.
- Same file, `:53-56` — the four tautological `expect_zero("rho_X - (1+zeta_X)", ...)` lines (formerly :46-49) have been replaced with four `expect_close("rho_X vs paper", rho_X, sp.Float("<14-digit paper constant>", 30), sp.Float("1e-13", 30))` calls anchoring each rho to the paper-stated `3.46622291347846`, `3.46752913273870`, `3.44257571477179`, `3.46752922945601`. Matches directive Step 2 verbatim.
- Same file, `:26` — banner relabel `"STAGE 69"` → `"STAGE 086"` (collateral fix per stage-086 banner directive; consistent with the unit ID and harmless).

The `rho_X` definitions at lines 43-46, the zeta literals at 38-41, the eps_cap block, and the derivative print block were left untouched, as required.

**Assessment:**
The new checks are non-tautological. Each `expect_close` now compares the *computed* `rho_X = Q(zeta_X; 0)` (numerically evaluated to 30 digits from the upstream zeta literal at 38-41) against a *fixed externally-quoted constant* taken from `paper/stages/stage_086.tex`. If a future edit corrupts any of the four `zeta_*` literals on lines 38-41, the corresponding `rho_X` would shift but the paper-side `sp.Float("3.46...", 30)` constant would not — the diff would exceed `1e-13` and the script would raise. This is a genuine paper-anchor, not a re-derivation of itself.

The 1e-13 tolerance is appropriate: the captured transcript shows the actual diffs are ~5e-17 to 1.2e-16 (i.e., machine-precision noise from the 14-digit input vs. the 30-digit symbolic evaluation), comfortably within tolerance but tight enough that a +1e-12 perturbation on any literal would flip exactly one assertion to FAIL.

### F2 — tautological_check (mathematica)

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage086_family1_loading_ratio_window_mathematica_audit.wl:53-56` — four new `expectApprox` calls inserted immediately after the `zetaMaxNum = …` line (formerly line 51), anchoring each of `zetaSuffChi`, `zetaFailChi`, `zetaSuffJ`, `zetaMaxNum` to the paper-stated 14-digit values `2.46622291347846`, `2.46752913273870`, `2.44257571477179`, `2.46752922945601` with tolerance `10^-14`. Matches the directive insertion-point and content verbatim.
- Same file, `:32` — banner relabel `"STAGE 069"` → `"STAGE 086"` (collateral fix consistent with the unit ID).
- The four ToExpression zeta literals at 48-51, the rho computations at 58-61, and the existing rho `expectApprox` checks at 68-71 are unchanged, as required.

**Assessment:**
The new zeta anchors compare *the literal values typed into the script* (`zetaSuffChi`, …, `zetaMaxNum`) directly against extra-script numeric paper constants `2.46622291347846`, etc. This is precisely the upstream-anchor that was missing. With these in place, the existing `rho_*` `expectApprox` checks at 68-71 become non-tautological: the chain is now (a) zeta-literal pinned to paper by L53-56, (b) rho computed by symbolic substitution, (c) rho pinned to paper by L68-71. A mistyped zeta would break L53-56 first.

The 10^-14 tolerance is tight: transcript shows actual diffs of `0.` (the input precision is exactly 14 digits and the targets are typed as 14 digits, so the difference is exactly zero). A +10^-13 perturbation on any literal would push the diff above 10^-14 and flip exactly one assertion to FAIL.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Q(zeta;0) - (1+zeta) = 0` (symbolic reduction intact)
- `rho_suff^(chi) vs paper diff = 1.214391841494943541093365767218853890273E-16`
- `rho_fail^(chi) vs paper diff = 1.070317295670973629640094647445676925917E-16`
- `rho_suff^(J)   vs paper diff = 8.129994879709555720190268933381680602833E-17`
- `rho_max        vs paper diff = 5.366711433453002476816172462459761803208E-17`

All four new diffs are well below the `1e-13` tolerance (by ~3 orders of magnitude), as expected.

**Mathematica:** exit=0. Notable lines:
- `PASS: dQ exact formula`, `PASS: Q(zeta;0) - (1+zeta)` (structural checks intact)
- `zeta_suff^(chi) vs paper diff = 0.` → `PASS: zeta_suff^(chi) vs paper`
- `zeta_fail^(chi) vs paper diff = 0.` → `PASS: zeta_fail^(chi) vs paper`
- `zeta_suff^(J) vs paper diff = 0.` → `PASS: zeta_suff^(J) vs paper`
- `zeta_max^(F1) vs paper diff = 0.` → `PASS: zeta_max^(F1) vs paper`
- Existing rho `expectApprox` PASS lines retained at 24, 26, 28, 30.

The four new zeta anchors appear before the rho-numeric checks, exactly as the directive required.

**Output freshness:** Both transcripts are newer than their scripts:
- sympy script mtime 1779898627; output mtime 1779899091 (~8 minutes later — fresh).
- mma script mtime 1779898632; output mtime 1779899174 (~9 minutes later — fresh).

## Material-change assessment

`material_change`: **false**.

No derived numeric or symbolic *results* changed. The edits are entirely additive (a new helper + four new assertions per engine) plus cosmetic banner relabels (`"STAGE 69/069"` → `"STAGE 086"`) that do not affect any computation. The `rho_X` values printed at L48-51 (sympy) and L63-66 (mma) are identical to the pre-fix transcripts. Downstream units that depend on Stage 086's outputs see no change.

## Side observations (non-blocking)

- The Codex banner relabel covers both engines (sympy used `"STAGE 69"`, mma used `"STAGE 069"`); the unified `"STAGE 086"` form is now consistent and matches the unit ID. Not part of either finding but flagged in the user message; no objection.
- The Mathematica `expectApprox` targets at lines 68-71 are still the extended-precision `3.46622291347846012143918414949` etc. values. Per F2's intent, these are now legitimate *secondary* anchors because the upstream zeta is independently pinned by L53-56 — not flagged as a regression.
- The notes section 2 / paper card stage_086.tex lines 17-28 referenced by the directive were not read by the verifier (per scope rule: scripts-only), but the directive vouches for them and the chained anchors are internally consistent.

## Verdict justification

Both findings are addressed correctly: F1 swaps four tautological self-identities for four direct numeric anchors of computed `rho_X` against externally-quoted paper constants (1e-13 tolerance, observed diffs ~1e-16), and F2 inserts four upstream-zeta anchors immediately after the zeta literal declarations that pin the literals to the paper-stated 14-digit values (1e-14 tolerance, observed diffs 0). Both engines exit 0 with the expected new transcript lines, and the existing structural assertions (Q-reduction, derivative formulas, eps_blk cap) remain intact. No regressions, no collateral edits beyond the banner relabel called out in the user note, and `material_change` is false — Stage 086's downstream results are unchanged.
