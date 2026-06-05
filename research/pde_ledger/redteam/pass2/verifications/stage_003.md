---
unit_id: 003
batch: I.1
verifier_model: claude-opus-4-8
verify_date: 2026-06-04
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 003

> Scope note: the original F1 (Mathematica "Lagrangian doubling") and F2 (stale
> output) were WITHDRAWN by the directive as false positives and are NOT
> re-raised here. The two live findings are F3 (notes index garble) and F4
> (`lRed` robustness consolidation).

## Per-finding outcomes

### F3 — paper_misalignment (notes index garble)

**Classification:** resolved

**What changed:**
`notes/stages/moving_throat_pde_stage003_bdg_coupling.md:367-372`. The three
grouped-P2 effective-coefficient index tokens were corrected from the garbled
`d_{237}^{\rm eff}/d_{238}^{\rm eff}/d_{239}^{\rm eff}` to
`d_{2,20}^{\rm eff}/d_{2,21}^{\rm eff}/d_{2,22}^{\rm eff}` in the three formulas
for `\bar d_2`, `a_2`, and `b_2`.

**Assessment:**
Correct and exactly scoped. The full `git diff` of the notes file shows ONLY
the three affected hunks; every changed line is a one-for-one substitution of
the index subscript (`237→2,20`, `238→2,21`, `239→2,22`) with no change to
equation structure, coefficients, or surrounding prose. The corrected indices
match the script channels A ∈ {20,21,22} (`d220/d221/d222` in the SymPy audit)
and the `.tex` card. Nothing else in the notes file changed.

### F4 — robustness consolidation of `lRed`

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:54-60`. `lRed`
is now a single fully-parenthesized assignment `lRed = ( ... );` capturing the
complete Lagrangian once, and the compensating `lRed = lRed + ( ... );` re-add
block (old lines 62-67) was removed. The `Clear[qa0, qL0, xa0, xb0, vqa, vqL,
vxa, vxb, aqa, aqL, axa, axb];` statement is retained.

**Assessment:**
Correct, math no-op, no collateral edits. The consolidated `lRed` contains
every term exactly once: qa/qL inertia (`maa`,`maL`,`mLL`), the K-potential
(`-1/2 kaa qa^2 - kaL qa qL - 1/2 kLL qL^2`), both BdG kinetic+potential terms
(`xa`,`xb` with `wa`,`wb`), and all four wall–matter couplings
(`c1a,c1b,c2a,c2b`). Term-by-term it is identical to the prior net value
(kinetic-only surviving first block + re-add block), so the algebra is
unchanged. This is confirmed independently: the committed output transcript
`mathematica/output/moving_throat_pde_stage003_bdg_mathematica_audit.txt` is
UNMODIFIED in git (only the `.wl` and notes `.md` are modified in the working
tree), i.e. byte-identical to before, and it shows all 19 PASS lines plus
"Stage 003 Mathematica audit passed." The diff is limited to exactly the
intended consolidation; nothing else in the script changed.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script was edited for this unit; no
`stage_003_sympy.log` was captured in `redteam/pass2/exec_logs/` (expected — F3
is a notes-only edit, F4 is Mathematica-only).

**Mathematica:** exit=0 (per orchestrator's independent re-run, already noted in
the directive). No `stage_003_mathematica.log` was captured in
`redteam/pass2/exec_logs/`; per instructions I relied on the committed
transcript instead of executing the script. The committed `.txt` contains 19
`PASS:` lines and the closing "Stage 003 Mathematica audit passed." banner,
consistent with the script's `Exit[0]`.

**Output freshness:** the committed `.txt` is intentionally byte-identical to
the pre-fix transcript (this is the acceptance criterion for the F4 no-op).
`git status` confirms only the `.wl` and notes `.md` are modified; the `.txt`
shows no diff. This is the expected and correct state for a consolidation that
does not change the math.

## Material-change assessment

`material_change`: false.

F3 is a notes documentation fix (no script/result impact). F4 leaves `lRed`
algebraically identical, so every derived quantity (EL equations, D0_eff,
series coefficients, dispersion roots, grouped-P2 projections, harmonic
overlap) is unchanged — verified by the byte-identical transcript. No
downstream unit can be affected.

## Side observations (non-blocking)

- The expected `redteam/pass2/exec_logs/stage_003_sympy.log` and
  `..._mathematica.log` files are absent (the directory contains only the diff
  patches for 003/006/009). This did not block verification because the
  orchestrator's independent re-run result is recorded in the directive and the
  committed transcript corroborates it, but the orchestrator may want to confirm
  log capture is wired up for later units.

## Verdict justification

Both live findings landed exactly as the directive specified: F3 is a tightly
scoped three-token notes correction matching the scripts/card, and F4 is a
clean single-assignment consolidation of `lRed` that removes the fragile re-add
while remaining a verified math no-op (byte-identical committed transcript, 19
PASS, exit 0). No regressions appear in either diff, no collateral edits, no
new findings. Verdict: verified; material_change: false.
