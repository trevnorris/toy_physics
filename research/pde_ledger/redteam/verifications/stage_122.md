---
unit_id: 122
batch: IV.3
verifier_model: claude-opus-4-8
verify_date: 2026-06-01
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 122

Retro-sweep dual-engine retrofit. Stage 122 was audited + verified in batch IV.3
(SymPy-only; its traction-ratio de-tautologization is already in the `.py`). The
sole directive finding here is **F1 = missing_mathematica**: Codex added a brand-new
independent-route Mathematica `.wl`. There is no fresh auditor report — this is
verified against `redteam/directives/stage_122.md` F1 and its claim manifest M1–M5.
Stale IV.3 report content is ignored.

## Per-finding outcomes

### F1 — missing_mathematica (dual-engine retrofit)

**Classification:** resolved

**What changed:**
New file `mathematica/moving_throat_pde_stage122_mouth_source_compensation_test_mathematica_audit.wl`
(133 lines, untracked / new — `git status` shows it as `??`, no prior commit touches it).
The SymPy reference engine
`scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py`
is UNMODIFIED: `git diff HEAD` against it is empty; `git status` lists only the `.wl`
as new. The directive's `## Applied: F1` block claims a single new file with
"deviation: none" — confirmed accurate. (The captured `stage_122_diff.patch` is empty
because the new file is untracked; verified directly from disk + git instead.)

**Assessment:** Correct, independent, and non-tautological on every count.

- **Independent route, not a transliteration.** `r_F1` is reconstructed from the
  geometric relation `Sqrt[12*radius^2/Pi^2 - 1]` with `radius = 37/20`
  (`rFromGeometry`, .wl:52,56) — not a bare numeric — and self-checked via
  `expectZero["r_F1 geometric relation", ...]` (.wl:59). `g_±` is DERIVED by
  `Solve[compensationResidual[rF1, g] == 0, g, Reals]` (.wl:61-62); the two roots are
  pinned to the `±` branches by SIGN of `# - rF1` via `SelectFirst` (.wl:66-73) — not
  hard-typed. The exact surd targets `gMinusExact`/`gPlusExact` (.wl:82-83) appear ONLY
  as the comparison RHS in `expectZero` (.wl:85-86), satisfying the anti-transliteration
  rule. This is a genuinely different decomposition than the SymPy script, which defines
  the roots directly as `rF ∓ sqrt(1+rF**2)/2` (.py:26-27). The defect is derived by
  substituting `gNat -> 1` and `rF1` into the residual polynomial (.wl:92), not re-typed.

- **THE CRITICAL GUARD — M5 traction ratio is genuinely non-tautological.**
  `cStage` (the stage-119 background constant `C`) is carried as a FREE POSITIVE SYMBOL
  (`Clear[..., cStage]` .wl:45; `cStage > 0` in `$Assumptions` .wl:49). Tractions are
  built as `C/g`: `tmNatRaw = cStage/gNat`, `tmMinusRaw = cStage/gMinus`,
  `tmPlusRaw = cStage/gPlus` (.wl:101-103). The ratio is FORMED and Mathematica COMPUTES
  the `C` cancellation: `reduceExact[tmMinusRaw/tmNatRaw]` (.wl:105), reported in the log
  as `(20*gNat*Pi)/(-37*Sqrt[3]+2*Sqrt[4107-100*Pi^2])` (log:32) — i.e. `gNat/gMinus`,
  with `C` gone and `gNat` STILL SYMBOLIC (not pre-set to 1). The residual is held
  symbolic: `tractionResidualMinusBeforeAnsatz = reduceExact[tractionRatioMinusRaw - 1/gMinus]`
  (.wl:107), which the log prints as the NONZERO `(20*(-1 + gNat)*Pi)/(...)` (log:34) =
  `(gNat-1)/gMinus`. That is exactly the required pre-ansatz witness: the residual reaches
  0 ONLY when `gNat -> 1` is substituted LAST (.wl:121-124). The script does NOT define
  `T_ratio := 1/g_±` directly — the `1/gMinus` on .wl:107 is the comparison target for an
  independently-computed ratio, not an X−X self-check. A mis-stated natural normalization
  (`gNat != 1`) would leave the residual nonzero and FAIL. Tautology guard fully satisfied,
  matching the IV.3 SymPy de-tautologization.

- **All M-items closed.** M1 (gminus/gplus exact form ×2, .wl:85-86), M2 (compensation
  quadratic at gminus/gplus ×2, .wl:88-89), M3 (defect closed form, .wl:98), M4 (natural
  off compensation, via `expectNonzero` that double-gates symbolic `res != 0` under
  assumptions AND numeric `Abs[N[res,40]] > 10^-6`, .wl:34-41,99), M5 (traction ratio ±
  identity ×2, .wl:121-128).

- **Cross-engine agreement (numerics match SymPy reference and directive anchors exactly):**
  g_− ≈ 0.75803507894466, g_+ ≈ 2.79795199200529, defect/C_nat ≈ 1.74016524722739,
  T-ratio(−) ≈ 1.31920016339112, T-ratio(+) ≈ 0.35740427386079. Exact symbolic surd forms
  also match the SymPy log.

## Exec log assessment

**SymPy:** exit=0 (reference engine, unchanged from IV.3; log dated 2026-05-29). Notable:
- `g_minus = (-37*sqrt(3) + 2*sqrt(4107 - 100*pi**2))/(20*pi)`
- `numeric defect = 1.7401652472273881285`
- `traction ratio (-) identity = 0`, `traction ratio (+) identity = 0`

**Mathematica:** exit=0. 9 PASS lines (requirement ≥8-9 met): r_F1 geometric relation;
gminus exact form; gplus exact form; compensation quadratic at gminus; compensation
quadratic at gplus; defect closed form; natural off compensation; traction ratio (-)
identity; traction ratio (+) identity. Notable:
- `derived g_minus = (-37*Sqrt[3] + 2*Sqrt[4107 - 100*Pi^2])/(20*Pi)` (log:13) — solved, not typed
- `traction residual (-), before natural ansatz = (20*(-1 + gNat)*Pi)/(-37*Sqrt[3] + 2*Sqrt[4107 - 100*Pi^2])` (log:34) — nonzero pre-ansatz residual, the non-tautology witness
- `traction ratio (-) identity = 0` / `traction ratio (+) identity = 0` (log:38,40)
- `Stage 122 Mathematica audit passed.`, exit_code 0

Comment hygiene: a `grep` for `*)` substrings in the `.wl` finds none.

**Output freshness:** `.wl` mtime 2026-06-01 15:49:26; its committed output
`mathematica/output/...mathematica_audit.txt` mtime 2026-06-01 16:05:46 — output is NEWER
than the script (regenerated post-fix). The SymPy `.txt` (2026-05-29 13:14:39) is older but
its reference engine is unchanged, so no regeneration was required.

## Material-change assessment

`material_change`: **false**. The retrofit adds a second, independent verification engine
only. No derived result changed: the `.wl` re-derives the SAME roots, defect, and ratios as
the unchanged SymPy reference (the IV.3-verified result), and the SymPy `.py` was not touched
(empty `git diff HEAD`). No downstream unit's inputs change; no re-audit warranted.

## Side observations (non-blocking)

- `reduceExact` wraps `FullSimplify[RootReduce[Cancel[Together[...]]], Assumptions -> $Assumptions]`
  plus `stripCE` — robust surd handling consistent with house idioms.
- `expectNonzero` double-gates (symbolic `!= 0` AND a 40-digit numeric floor), a stronger M4
  non-vanishing check than the SymPy `expect_nonzero` (symbolic only). Strengthens, no issue.

## Verdict justification

The single finding (F1, dual-engine retrofit) is fully resolved. Codex added a genuinely
independent Mathematica route that SOLVES the compensation quadratic and selects branches by
sign (not hard-typing roots), derives the defect by substitution, and — critically — computes
the `C` cancellation symbolically while holding `gNat` free, exposing the nonzero pre-ansatz
residual `(gNat-1)/g_±` that collapses to 0 only under the `gNat -> 1` natural-normalization
ansatz applied last. The exec exits 0 with 9 PASS lines; all numerics and exact forms agree
with the unchanged SymPy reference; the reference `.py` is untouched; comment hygiene is clean;
and the output `.txt` was regenerated post-fix. Verdict: verified; material_change false.
