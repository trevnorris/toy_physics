---
unit_id: 166
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T22:55:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 166

This is a recovery-pass verification. The authoritative finding source is the V.2
directive `redteam/directives/stage_166.md` (its single F1, with `## Applied: F1`
block), which encodes the Codex read-only review `redteam/codex_reviews/stage_166.md`
(R1 = `tautological_check`, verdict-table row 13) and the Claude+Codex consult
`_consult_batch7.md` Q5 (DECIDED — REPLACE, do not keep). The prior v2 auditor report
`reports/stage_166.md` carried two findings (F1 general-form asserts; F2 banner +
matrix-inverse block); the directive's "Already applied" section documents those as
confirmed-PASS (review verdict-table rows 1-12) and NOT re-issued. The only live finding
for this pass is the tautological matrix round-trip. I verify against that one finding.
(Note: the older `verifications/stage_166.md` v2 content — which itself flagged this
round-trip tautology in its side observations — is superseded by this file.)

## Per-finding outcomes

### F1 — tautological_check (matrix round-trip vacuous; replace with forward transcription)

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl`.
The diff (`exec_logs/stage_166_diff.patch`) shows the old block at wl:73-76 —
`(* Round-trip ... *)` plus
`expectZero["matrix round-trip", Total[(Mmat . solVec - {dTheta, dKs, dKq, dP})^2]];` —
was DELETED and REPLACED (not supplemented) by the new block now live at wl:73-79:
the four-line comment, then
`fwdLaws = {2*drho, drho + 2*da, dZ + 2*dcs - 2*da, 5*(dcs - da)};` and
`expectZero["matrix forward-transcription", Total[(Mmat . {drho, da, dcs, dZ} - fwdLaws)^2]];`.
This is exactly the directive's "Insert exactly" text. The diff is scoped to these lines
only — no collateral edits. The SymPy script is untouched (mtime 2026-05-28, predates
this pass), as the directive required (SymPy has no matrix route).

**Assessment:** Correct and genuinely non-tautological. Verified against the four crux
specifics:

1. **Uses the symbolic DRIFT vector, not the inverse.** The residual multiplies
   `Mmat . {drho, da, dcs, dZ}` — the same primitive drift symbols `Mmat` is built to
   act on — NOT `Mmat . solVec` and nothing from `Inverse[Mmat]`. The old line computed
   `Mmat . Inverse[Mmat] . v - v ≡ 0` (identity cancellation for any invertible Mmat);
   the new line computes `Mmat . drift - fwdLaws`, which has no such cancellation. A
   wrong `Mmat` coefficient now propagates into the residual ⇒ the assertion FAILS.
   The retained inverse construction (wl:67-68) and the four `matrix drho/da/dcs/dZ`
   target asserts (wl:69-72) still use `Inverse[Mmat]`, but those compare the inverse to
   the paper §2 targets (a meaningful comparison, review rows 9-12 = PASS); only the
   round-trip was the tautology, and it is the only thing replaced.

2. **`fwdLaws` is literally hand-typed and matches the laws/Mmat rows.**
   `fwdLaws = {2*drho, drho + 2*da, dZ + 2*dcs - 2*da, 5*(dcs - da)}` is a literal list,
   not built from `Mmat`/`Solve`/`Inverse`. In column order `(drho, da, dcs, dZ)` it
   equals `Mmat . {drho, da, dcs, dZ}` for the live matrix: rows 1,2 `{2,0,0,0}`,
   `{1,2,0,0}` give `2 drho`, `drho + 2 da`; rows 3,4 `{0,-2,2,1}`,`{0,-5,5,0}` give
   `-2 da + 2 dcs + dZ = dZ + 2 dcs - 2 da` and `-5 da + 5 dcs = 5(dcs - da)`. These are
   the notes §1 boxed branch laws (eq1: 2δlnρ; eq2: δlnρ+2δln a; eq3: δlnZ+2δln c_s−2δln a;
   eq4: 5(δln c_s−δln a)). Residual is therefore identically 0 ⇒ PASS for the live,
   correct matrix; rows 3,4 = `{0,-2,2,1}`,`{0,-5,5,0}` confirmed (consistent with the
   consult). No numeric tolerance introduced — exact symbolic `=== 0` via
   `FullSimplify[Expand[...]]`.

3. **Old round-trip assertion is GONE.** The diff shows it deleted, not duplicated; the
   live file at wl:73-79 contains only `matrix forward-transcription`. No
   `matrix round-trip` / `Mmat . solVec - {dTheta...}` text remains anywhere in the file.
   Replaced, not supplemented.

4. **Exec log shows the new PASS at exit 0.** `exec_logs/stage_166_mathematica.log`
   (exit_code 0) shows, under "Independent matrix-inverse cross-check":
   `matrix forward-transcription = 0` / `PASS: matrix forward-transcription`, plus the
   four prior `matrix drho/da/dcs/dZ = 0` PASS lines and all general-form / forward /
   bundle / frozen asserts still PASSing. No `matrix round-trip` line appears.

**Anti-tautology reasoning (the crux, explicitly checked):** Perturb one `Mmat` entry,
e.g. row 4 `{0,-5,5,0} → {0,-4,5,0}` (a `5 → 4` coefficient error) or a row-3 sign flip.
Then `Mmat . {drho,da,dcs,dZ}` differs from the hand-typed `fwdLaws` in that component;
`Total[(...)^2]` becomes a nonzero sum of squares of a symbolic expression that
`FullSimplify` cannot reduce to 0; `expectZero` calls `fail` → `Exit[1]`. The old
round-trip could not detect this (it stayed 0 for any invertible perturbed matrix). The
new check genuinely has teeth. The directive's anti-tautology and no-fabricated-literal
guards are satisfied.

## Exec log assessment

**SymPy:** exit=0. Untouched this pass. Notable lines:
`drho general = 0`, `da general = 0`, `dcs general = 0`, `dZ general = 0`;
`bundle identity for dcs = 0`, `bundle identity for dZ = 0`; `frozen drho/da/dcs/dZ = 0`;
`rho_w^(chi) = 0.403417022451042341`. The prior-F1 general-inversion asserts still pass.
No matrix route in SymPy by design.

**Mathematica:** exit=0. Notable lines:
`matrix forward-transcription = 0` + `PASS: matrix forward-transcription` (the new check);
`matrix drho/da/dcs/dZ = 0` + PASS (prior matrix target asserts, intact);
`drho/da/dcs/dZ general = 0` + PASS (prior F1 asserts, intact);
`Theta/Ks/Kq/P0 law = 0` + PASS; bundle + frozen blocks all PASS;
`rho_w^(chi) = 0.40341702245104232684...`; final `Stage 166 Mathematica audit passed.`
No `matrix round-trip` line present. PASS count unchanged (one assertion replaced by one),
exactly as the directive predicted.

**Output freshness:** confirmed. `mathematica/output/...txt` mtime 2026-05-29 22:49:21 >
`.wl` mtime 22:46:11; `scripts/output/...txt` mtime 2026-05-29 22:49:21 > `.py` mtime
2026-05-28 15:54:48. Both committed `.txt` outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false. The edit changes only HOW the Mathematica engine self-checks
its matrix transcription (a tautological assertion replaced by a real one). It does not
alter any solved form, printed coefficient, numeric (`ρ_w^(χ)` unchanged), or the matrix
entries themselves — `Mmat` was already correct and is untouched. No derived result that
downstream units depend on changes. This is a how-it's-checked fix, matching consult Q5
("does NOT change Stage 166's published conceptual claim"). Downstream units need not be
re-audited on account of unit 166.

## Side observations (non-blocking)

- The retained matrix block (wl:67-72) still uses `Inverse[Mmat]` to form `solVec` and
  compares it to the paper §2 targets — a meaningful comparison (review rows 9-12 = PASS),
  not a tautology. With the new forward-transcription check now also pinning `Mmat` against
  the hand-typed forward laws, both the forward map and its inverse are independently
  anchored to external (notes §1/§2) targets — a strict improvement. Not a finding.
- The `## Applied: F1` block in the directive accurately describes the edit (files_changed,
  one-sentence summary, deviation: none) — consistent with the diff.

## Verdict justification

The single live finding (R1/F1, tautological matrix round-trip) is fully resolved. The
vacuous `Total[(Mmat . solVec - v)^2]` assertion was deleted and replaced verbatim with
the prescribed hand-typed forward-transcription check
`Total[(Mmat . {drho,da,dcs,dZ} - fwdLaws)^2]`, which is genuinely non-tautological: it
multiplies the symbolic drift vector (not the inverse), `fwdLaws` is a literal list of the
notes §1 boxed laws matching the `Mmat` rows, and a perturbed `Mmat` coefficient would now
make the residual nonzero and fail. Both engines exit 0, the new PASS line appears, the old
line is gone, the prior PASS asserts are intact, and the saved outputs are fresh. No
regressions in the diff, no collateral edits, SymPy untouched as directed. Verdict: verified.
