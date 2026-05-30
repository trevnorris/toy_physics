---
unit_id: 175
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 175

This pass verifies the SINGLE re-review finding F1 = R1 (`codex_reviews/stage_175.md`):
the Mathematica `Sigma_N` differential block was a line-by-line transliteration of the
SymPy block (both extracting the first-order log-slope via the same `dlog =
D[Log[.],eps]/.eps->0` primitive), so the differential SLOPE identity was singly-routed
across engines. The stale V.1 F1/F2 were already resolved + PASS in prior work and were
intentionally NOT re-prescribed (per the batch-8 consult); they are not re-evaluated here.

## Per-finding outcomes

### F1 — mathematica_transliteration (R1)

**Classification:** resolved

**What changed:**
Codex applied F1 as the directive's PRIMARY edit (option B SUPPLEMENT), Mathematica-only,
`.py` untouched. Two additions to
`mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl`
(confirmed against the diff patch and the live file):

1. New slope extractor at line 31, immediately after the existing `dlog` helper (lines 26-29):
   `dlogSeries[expr_] := Coefficient[Normal[Series[Log[expr], {eps, 0, 1}]], eps];`
   The load-bearing operation is `Series` + `Coefficient` (Mathematica-native), NOT
   `D[Log[...]]` and NOT a port of SymPy's `sp.diff(sp.log(...))`.

2. Exactly ONE new `expectZero` check at lines 106-108, placed immediately after the
   existing `dlog`-based `Sigma_N - dln(Lambda^2/K)` line (line 100), which is left
   untouched:
   ```
   sigmaNDirectSeries = FullSimplify[2*dlogSeries[exprPoverDeltaPhys] - kappa, Assumptions -> $Assumptions];
   sigmaNShapeSeries  = FullSimplify[dlogSeries[(lambda^2/k) /. subsEps], Assumptions -> $Assumptions];
   expectZero["Sigma_N - dln(Lambda^2/K) [series route]", sigmaNDirectSeries - sigmaNShapeSeries];
   ```
   It reuses the exact live names (`exprPoverDeltaPhys`, `lambda`, `k`, `subsEps`, `kappa`,
   `$Assumptions`); none renamed or redefined. The diff patch contains only these two hunks
   (the helper insertion and the supplemental block + its explanatory comment). No collateral
   edit: `sigmaNDirect`/`sigmaNShape` (lines 98-100), the downstream `sigmaNCommon` (line 120)
   and the `Xi_load` weighted aggregate (lines 124-128) are byte-for-byte unchanged, and the
   `.py` is not in the diff at all. The directive's `## Applied: F1` block records
   `deviation: none`, consistent with the diff.

**Assessment:**
Correct and complete; it addresses the finding on its own merits, not by rubber-stamp.

- Independence (anti-guard 2): the slope is extracted via `Series`/`Coefficient`, a
  structurally distinct route from the mirrored `dlog = D[Log[.],eps]/.eps->0` primitive
  that R1 flagged as a transliteration of the SymPy `sp.diff(sp.log(...))` route. The new
  check is NOT the degenerate `dlogSeries[e] - dlog[e]` on the SAME argument (which would
  only validate two differentiation methods against each other — a method tautology). It
  compares series-route DIRECT (`2*dlogSeries[exprPoverDeltaPhys] - kappa`) against the
  SHAPE target (`dlogSeries[(lambda^2/k)/.subsEps]`), exactly as the directive requires.

- Non-tautological (anti-guard 3): the check can fail. The SHAPE target carries the `1/k`
  factor whose first-order log-slope is `-kappa`; a wrong shape (e.g. `lambda^2*k`) would
  flip that term to `+kappa` and break the identity, and an omitted/wrong `-kappa` on the
  DIRECT side would not cancel. As the directive notes honestly, the `2 dln(P/Delta)` vs
  `2 dln(Lambda)` portion is value-equal — `lambda` is the FullSimplify of `(p/delta)/.subsHat`
  and `exprPoverDeltaPhys` is the same `(p/delta)/.subsHat` then `/.subsEps`. So the
  load-bearing degrees of freedom are the symbolic `-kappa` (= `-delta_K`) term against the
  `/k` in the SHAPE target, PLUS the independent extraction-method coverage. That is exactly
  the intended R1 fix — a structurally independent slope ROUTE, not new physics — and it is
  reported as such (the comment says "the differential identity no longer relies on a
  transliteration of the SymPy dlog route"), not over-claimed.

- `-kappa` (= `-delta_K`) is kept symbolic; no numeric substitution (anti-guard).

- The escape clause (Block + sanctioned MIRROR_POLICY mirror / `## Blocked: F1`) was
  correctly NOT used: the series route reduced to `=== 0` robustly under `FullSimplify`, so
  no waiver applies and the directive's primary independent-route resolution was achieved.

## Exec log assessment

**SymPy:** exit=0. Transcript UNCHANGED (the `.py` was not edited), as required — no
`[series route]` line:
```
Sigma_N - dln(Lambda^2/K) = 0
Common-shape branch Sigma_N + dK = 0
Xi_load (all shapes frozen) + dK = 0
```

**Mathematica:** exit=0. The new line appears at log lines 32-33, immediately after the
existing `dlog` line, both PASS:
```
Sigma_N - dln(Lambda^2/K) = 0
PASS: Sigma_N - dln(Lambda^2/K)
Sigma_N - dln(Lambda^2/K) [series route] = 0
PASS: Sigma_N - dln(Lambda^2/K) [series route]
```
Closing line `Stage 175 Mathematica audit passed.` with `# exit_code: 0`. Every prior
check (B0/Delta/Q/P/Z0/N0, Sigma_B/Sigma_Z, conservative branches, Xi_load) still PASSes —
inserting the `dlogSeries` helper and the one supplemental check introduced no regression.

**Output freshness:** confirmed. The `.wl` was edited at 2026-05-29 23:44:10; both saved
`.txt` outputs were regenerated at 2026-05-29 23:50:24 (newer than the `.wl`). The saved
Mathematica `.txt` contains the `Sigma_N - dln(Lambda^2/K) [series route] = 0` / PASS lines
and `Stage 175 Mathematica audit passed.`; the saved SymPy `.txt` has no `[series route]`
line, consistent with the untouched `.py` (mtime 2026-05-28, well before the edit). Outputs
are fresh and match the captured exec logs.

## Material-change assessment

`material_change`: false.

The edit ADDS a corroborating independent slope check; it does not alter any derived result,
constant, identity target, or printed value that downstream units could consume.
`sigmaNDirect`/`sigmaNShape`/`sigmaNCommon`/`Xi_load` are unchanged, and the SymPy reference
engine is untouched. Nothing moved, so no downstream unit > 175 can depend on a value that
changed; no specific re-audit concern.

## Side observations (non-blocking)

- Line-number drift: the directive/consult cite older line ranges (`dlog` at 26-29, Sigma_N
  block at 95-98, `wl:113`/`wl:118`); the live file now has the `dlog` Sigma_N line at 100
  and the downstream checks at 120-128, due to the V.1 F1/F2 edits and the new 3-line
  supplement. The references still resolve to the correct constructs — purely cosmetic, not
  a finding.
- `dlogSeries` deliberately omits the inner `FullSimplify[..., Assumptions]` wrapping that
  `dlog` applies before differentiating. The directive permitted that wrapping only "if
  needed to land `=== 0`"; it was not needed (the outer `FullSimplify` on
  `sigmaNDirectSeries - sigmaNShapeSeries` reduced the residual to 0), so the simpler form
  matches the directive's prescription verbatim. Not a defect.

## Verdict justification

The single finding F1 (R1 transliteration) is resolved. Codex added the prescribed
Mathematica-native `dlogSeries` (Series+Coefficient) extractor and exactly one supplemental
series-route `expectZero` comparing series-route DIRECT vs the SHAPE target, leaving the
existing `dlog` route, the downstream checks, and the entire SymPy reference engine
untouched — exactly the consult-chosen smaller-blast-radius option B. The new check is
genuinely independent of the SymPy route, is non-tautological (the `/k` SHAPE factor and the
symbolic `-kappa` are load-bearing), and passes `=== 0`; both engines exit 0, the SymPy
transcript is unchanged, and the saved outputs are fresh and contain the new PASS line. No
regression. Verdict: verified, material_change false.
