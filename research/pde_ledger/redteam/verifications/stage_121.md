---
unit_id: 121
batch: IV.3
verifier_model: claude-opus-4-8
verify_date: 2026-06-01
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 121

> Retro-sweep dual-engine retrofit. Stage 121 was already audited + verified in batch IV.3
> (SymPy-only). The only finding to verify is the directive's F1 = `missing_mathematica`:
> Codex added a brand-new independent-route Mathematica `.wl`. Verified against the directive
> F1 and its claim manifest M1–M6 (no fresh auditor report exists for this retrofit). The
> SymPy `.py` is the unchanged reference engine.

## Per-finding outcomes

### F1 — missing_mathematica (dual-engine retrofit)

**Classification:** resolved

**What changed:**
Codex created the new (untracked) file
`mathematica/moving_throat_pde_stage121_geometric_r_selection_mathematica_audit.wl` (122 lines).
No other file was modified — the diff patch (`stage_121_diff.patch`) is 0 bytes because the
new `.wl` is untracked, and `git status --porcelain` shows only the `??` (untracked) `.wl`.
The SymPy reference `scripts/...stage121...sympy_audit.py` is UNCHANGED (last touched
2026-05-27 16:31, the IV.3 batch; its last commit is `b4e02d8` "batch IV.3 close", and it is
not listed as modified). Confirmed: Codex did not modify the reference engine.

**Assessment:** The `.wl` is a genuinely independent route, not a transliteration:

- **Independent derivation of r_geom (M1).** The `.wl` defines the Stage-99 length law
  `stage99TubeLength = (Pi*throatScale/2)*Sqrt[(1+radius^2)/3]`, then `Solve`s
  `tubeLength == stage99TubeLength` for `radius` over `Reals`, and derives the positive root —
  it does NOT retype `Sqrt[12 L²/(π²a²) − 1]`. The surd appears only as the RHS of the M1
  comparison (`derivedRadius^2 - (12*ell^2/(Pi^2*throatScale^2) - 1)`), and the log prints the
  Solve-derived form `Sqrt[12*ell^2 - Pi^2*throatScale^2]/(Pi*throatScale)` — a different
  factored form than the M1 target, so the squared-residual check is a real cross-check, not a
  self-subtraction.
- **Explicit positive-branch selection (guard a).** `positiveRadiusRoots = Select[radiusRoots,
  TrueQ[FullSimplify[# > 0, ...]] &]` followed by `If[Length[...] =!= 1, fail[...]]`. This is an
  explicit positivity filter with a uniqueness guard — NOT an accidental `[0]` index grab (the
  SymPy reference uses `solve(...)[0]` relying on the `positive=True` symbol; the `.wl` route is
  meaningfully different). Branch positivity is backed by `branchAssumptions` declaring the
  radicand positive and the aspect ratio above threshold.
- **Different names / decomposition.** `ell, throatScale, tubeLength, radius, soundSpeed` vs the
  SymPy `L, a, L_W, r, c_s`; the `.wl` keeps `tubeLength` as a distinct symbol used for the M2
  back-substitution and only later substitutes `tubeLength -> ell`.
- **M2 tube-length relation:** back-substitutes `derivedRadius` into the length law and checks it
  returns `ell` exactly → exact 0.
- **M3 nested-radical equality (guard b).** `rF1Derived` is produced by substituting
  `ell -> (37/20)*throatScale` into the Solve-derived `derivedRadius`; the equality
  `rF1Derived - rF1Target` is certified SYMBOLICALLY through `expectZero` (which runs
  `FullSimplify[Together[Expand[...]]]` and requires `=== 0`), not numeric-only. The ≈ decimal is
  reproduced separately via `N[rF1Derived, 20] = 1.77799353547497811851...`.
- **M4 r_c:** `rcF1Derived = FullSimplify[rF1Derived^2]`, log shows `-1 + 4107/(100*Pi^2)`,
  `N = 3.16126101219081227320...`, residual against target exact 0.
- **M5 Omega_W (acknowledged low-information parity):** `halfWavePole = Pi*soundSpeed/(2*tubeLength)`
  evaluated at `tubeLength -> ell` compared to `Pi*soundSpeed/(2*ell)` → exact 0. The directive
  explicitly accepts this as a shallow definitional-parity check; not failed for being shallow.
- **M6 threshold (guard c).** Uses the EXACT symbol `thresholdAspect = Pi/(2*Sqrt[3])` (not a
  decimal); `derivedRadius /. ell -> thresholdAspect*throatScale` FullSimplifies to exact 0.
- **Comment hygiene (guard):** no `*)` premature-close substring in the file (grep confirmed none).

All six manifest items M1–M6 map one-to-one to the six PASS labels in the exec log, and the
`.wl` conclusions agree with the SymPy reference (same r_geom form, r_F1 surd, r_c, threshold,
Omega_W). Cross-engine agreement confirmed.

## Exec log assessment

**SymPy:** exit=n/a. The captured SymPy log (`stage_121_sympy.log`) is empty (1 blank line) —
expected: this is a `.wl`-only retrofit; the SymPy reference was deliberately not re-run/modified.

**Mathematica:** exit=0. All 6 expected PASS lines present, no FAIL:
```
PASS: r_geom closed-form (explicit)
PASS: tube-length relation
PASS: r_F1 symbolic value at L/a = 37/20
PASS: r_c(F1) symbolic value
PASS: Omega_W identification at L_W = L
PASS: r_geom vanishes at existence threshold
...
STAGE 121 MATHEMATICA AUDIT PASSED
# exit_code: 0
```
`grep -c "^PASS:"` = 6 (≥6 required). No parser-skip: every PASS is preceded by its printed
residual `= 0` line, so this is not a passing-rc-with-missing-PASS case.

**Output freshness:** confirmed. The saved
`mathematica/output/...stage121...mathematica_audit.txt` mtime is 2026-06-01 16:05:39, newer than
the `.wl` script mtime 2026-06-01 15:48:26 → regenerated post-fix. The SymPy output `.txt`
(2026-05-27 16:38) matches its unchanged reference engine.

## Material-change assessment

`material_change`: false. The change adds a second (Mathematica) verification engine only. No
derived value, constant, identity target, or paper number was moved: r_geom, r_F1
(≈1.77799353547498), r_c (≈3.16126101219081), the existence threshold (π/(2√3)), and the Omega_W
half-wave relation all match the pre-existing SymPy reference exactly. The reference SymPy script
is byte-for-byte unchanged. Therefore no downstream unit (> 121) depends on any altered result;
no `upstream_stale` propagation is warranted from this retrofit.

## Side observations (non-blocking)

- The new `.wl` is currently untracked (`??`); the orchestrator will need to `git add` it when
  committing the batch. This is expected for a freshly created retrofit file and is not a defect.
- `branchAssumptions` declares both `ell/throatScale` and `tubeLength/throatScale` above threshold
  plus both radicands positive; `$Assumptions` is later reset to `baseAssumptions` before the M6
  threshold substitution (correct — at threshold the strict-positivity branch assumption would be
  inconsistent with r_geom = 0). The scoping is handled correctly.

## Verdict justification

The sole finding (F1, missing_mathematica) is fully resolved. Codex added an independent-route
Mathematica audit that DERIVES the geometric radius from the Stage-99 length law via `Solve` with
an explicit positive-branch `Select`/uniqueness guard, certifies the nested-radical M3 equality
symbolically, uses the exact π/(2√3) threshold symbol for M6, and reaches the same six
conclusions as the unchanged SymPy reference. The exec log exits 0 with all six expected PASS
lines and no FAIL, the output `.txt` was regenerated post-fix, and no derived value moved
(`material_change: false`). Verdict: verified.
