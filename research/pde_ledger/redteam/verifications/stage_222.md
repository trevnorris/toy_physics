---
unit_id: 222
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 222

## Per-finding outcomes

### F1 — paper_misalignment (notes_contradicts_script, user-resolved)

**Classification:** resolved

**What changed:**
`notes/stages/moving_throat_pde_stage222_..._sympy_audit.md:395` — the λ_W=0.2 upper-wall
`R_{Q,*}` table cell was changed from `213.483858657863` to `145.483858657863`. This is the
notes-only edit the user authorized in the `## RESOLVED` block (script is authoritative; no
script change). Codex appended `## Applied: F1` to the directive.

**Assessment:**
Correct and exactly scoped. Acceptance criteria met:
- `213.483858657863` no longer appears anywhere in the notes file (grep: not found).
- `145.483858657863` is present in the 0.2-row R_{Q,*} cell (line 395), confirmed.
- The corrected value is a `+68` integer-part slip with the fractional tail `…858657863`
  unchanged, exactly as the directive described.
- The rest of the scan table (rows 0.4/0.6/0.8/1.0, and the P0/D0/ω_* columns of the 0.2 row)
  is unchanged and matches the SymPy output verbatim (notes 0.2 row P0=0.00594740531769,
  D0=2.82723442158450, ω_*=2.04402272302752 ↔ output line 53).
- No edits to the script, paper.tex, or appendix (only notes/*.md touched per git status for
  this unit's files).

**Cross-engine corroboration:** the NEW Mathematica engine independently computes the 0.2-cell
via `scanRow[1/5]` from `NSolve` roots and the symbolic `rqExpr`, printing
`145.48385865786412765…` (mathematica log line 88) and passing `M6 lambda_W scan row 1`
(residual on the R_Q column ≈ 2.8e-14). So the notes-typo direction (145, not 213) is confirmed
by an independent second engine, not just by the SymPy literal.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New file `mathematica/moving_throat_pde_stage222_..._mathematica_audit.wl` (310 lines, exact
mandated path). It verifies claim manifest M1–M6 and exits 0 with 26 PASS / 0 FAIL.

**Assessment — independence / anti-transliteration (PASS):**
The `.wl` is a genuinely independent re-derivation, not a port of the `.py` choreography:
- **M2 (quartic):** the `.py` hand-expands `F_y` term-by-term (py lines 92–95) and checks
  `cancel(D − F_y/(…))`. The `.wl` instead builds `dBranch = Together[kbExpr − qExpr/deltaExpr]`,
  clears the denominator natively via `Numerator[Together[dBranch]]`, substitutes
  `omega → Sqrt[y]`, and reads degree/coeffs off `Exponent`/`CoefficientList` (wl 138–148).
  It then ALSO builds a separate product-form `fProductY` and cross-checks
  `quarticFromD − fProductY == 0` (wl 157), plus reconstructs `Together[D]` from
  numerator/denominator (wl 153–156). This is the directive-mandated different decomposition
  (native denominator clearing vs. pre-typed hand-expansion), and both the native-numerator and
  the product form independently test degree 4.
- **M3 (residue/linewidth):** non-tautological. `aqqStar = 1/d0p`,
  `gammaStar = gamma5·ω_*^5·N_*/d0p`, and `FullSimplify[aqqStar/gammaStar]` must cancel the free
  pole-derivative symbol `d0p` to reach `27 cs^5/(aScale^5 ω_*^5 N_*)`; residual = 0 (log line 37).
  `d0p` is a free real symbol present in both numerator and denominator — the equality holds only
  because it genuinely cancels, not by self-equality. The survival threshold is derived by
  actually SOLVING `Solve[stagePeak == deltaVreq, rqRatio]` (wl 179) and comparing to the stage
  form, rather than re-asserting a typed literal — an independent re-derivation.
- **M5/M6 (figures):** R_Q is recomputed from primitives — `nFun = Together[pPrimitive^2/deltaExpr^2]`,
  `rqExpr = Together[27 cs^5/(aScale^5 omega^5 nFun)]` (wl 206–208) — and poles come from
  `NSolve[…, WorkingPrecision→50]` (vs SymPy `nroots`), a different root-finder. The expected
  numerics are regression anchors compared against these in-engine derivations; M5/M6 do NOT
  echo the `.py` literals as inputs. The 0.2-cell is re-derived independently (see F1 above).

**Assessment — claim coverage (PASS):** all M1–M6 genuine:
- M1 overlap `2√2/π` via native `Integrate` (residual 0).
- M2 quartic identity + degree 4 (two independent forms).
- M3 cancellation + peak form + threshold (all residual 0).
- M4 static slice Δ0/D0/N0/P0 (residuals ≤ 5e-17).
- M5 wall/internal/coupled roots + 4 R_Q figures (residuals ≤ 2e-13).
- M6 two thresholds + λ_W scan with monotone P0↑ / R_Q↓ + explicit 0.2-cell report.

The `assert_close`/`expectClose` tolerances (1e-12 for static/threshold, 1e-10 for roots/R_Q)
have ample margin over observed residuals (≤ ~2e-13). No fabricated literals: the pole census,
R_Q figures, and scan table are regression anchors over in-engine symbolic/NSolve derivations.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `kappa = 2*sqrt(2)/pi` (line 9); `R_Q,* = 27*c_s**5/(N_star*a**5*omega_star**5)` (line 22).
- 0.2 scan row `(0.2, …, 2.04402272302752, 145.4838586578641)` (line 53) — script value matches
  the corrected notes cell.
- `All Stage 222 symbolic and numerical checks passed.`

**Mathematica:** exit=0, 26 PASS / 0 FAIL. Notable lines:
- `M3 residue/linewidth cancellation residual = 0` (line 37) — d0p cancels.
- `native quartic degree = 4` (line 21); both degree-4 forms PASS.
- `lambda_W=0.2 upper-wall R_Q,* = 145.48385865786412765…` (line 88) — independent corroboration
  of the F1 correction direction.
- `All Stage 222 Mathematica checks passed.`

**Output freshness:** both `.txt` outputs were regenerated post-fix. SymPy output mtime
2026-06-02 16:15:50 (> script mtime 2026-05-11). Mathematica output mtime 2026-06-02 16:15:50
(> `.wl` mtime 2026-06-02 14:30). Both fresh.

## Material-change assessment

`material_change`: false.

F1 is a notes-only typo correction that brings the published cell into line with the
already-verified script value (no derived result changed). F2 adds a second engine that confirms
existing figures; it derives nothing new that downstream units consume. No script numerics, no
paper.tex, no appendix changed. No downstream unit depends on a value that moved.

## Side observations (non-blocking)

- The captured `stage_222_diff.patch` is 0 bytes; the new `.wl` and its output are untracked
  (`??` in git), so a `git diff` against HEAD shows nothing for them. I verified directly from the
  working-tree files and exec logs instead, which is sufficient. Not blocking.
- M2's `quarticFromD` is post-processed with `PowerExpand[… /. omega → Sqrt[y]]`; under the
  positivity assumptions in `$Assumptions` this is sound (varpi,ω_U,ω_W > 0, y the squared
  frequency), and the result is independently cross-checked against the product form, so no
  branch-cut concern. Not blocking.

## Verdict justification

Both findings are resolved. F1's authorized notes-only edit is in place and exactly scoped
(213→145.483858657863, sole cell changed, rest of table intact, no script change), and the new
Mathematica engine independently re-derives 145.4838586578641… for that cell — confirming the
correction direction cross-engine. F2's new `.wl` is a genuine independent re-derivation (native
`Numerator[Together]`/`CoefficientList`/`NSolve`/`Integrate` through a different decomposition than
the SymPy hand-expansion, with a non-tautological d0p cancellation and primitive-derived R_Q),
covers M1–M6 with adequate tolerances, and passes 26/26 with exit 0. Outputs are fresh. No
regressions. `material_change: false`.
