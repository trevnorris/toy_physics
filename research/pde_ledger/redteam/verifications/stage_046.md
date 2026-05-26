---
unit_id: 046
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 046 (batch III.1 v2)

This report supersedes the prior v1 verification (the v1 audit raised two
script-side findings — mathematica_transliteration and insufficient_verification
— which Codex resolved by editing the Mathematica and SymPy scripts; those edits
remain in place and are still correct). The v2 paper-grounded re-audit raised a
single new finding, F1, which is a notes-only paper_misalignment. This report
documents only the v2 resolution.

## Per-finding outcomes

### F1 — paper_misalignment (notes_contradicts_script)

**Classification:** resolved

**What changed:**
Notes-only edits to
`/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md`
at exactly the five coefficient sites the directive enumerated:

- Line 90 (P_R): `230 R delta^3` -> `162 R delta^3`, and `230 R delta xi^2` -> `162 R delta xi^2`.
- Line 132 (P_1): `248 R delta^2 xi` -> `180 R delta^2 xi`.
- Line 133 (P_1): `230 delta^3` -> `162 delta^3`.
- Line 137 (P_2): `237 R^2 xi^4` -> `220 R^2 xi^4`.

`git diff HEAD -- notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md`
shows exactly these five substitutions and nothing else.
`git diff HEAD --` against the SymPy and Mathematica audit scripts
(`scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py`,
`mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl`)
is empty: scripts were not modified in this batch, consistent with the
directive's user-approved direction (a).

A grep across the notes file for `230|248|237` confirms no stale coefficient
survives in any of the disputed sites, and a grep for the corrected
`162|180|220` confirms the new values appear exactly at lines 90, 132, 133,
and 137 — the four lines the Applied block listed.

**Assessment:**
The fix implements the math-correct direction (a) from the directive's
`## Resolve before fix_loop` block, the same direction the user approved in
the batch III.1 Q4 apply session. Both engines independently derive the
corrected coefficients from the shared definition of `F_tr`:

- The Mathematica saved output at
  `mathematica/output/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.txt`
  line 14 (`dF_tr/dR = ...`) contains `162*delta^3*r` and `162*delta*r*xi^2`,
  matching the corrected P_R coefficients in notes line 90.
- The same Mathematica output line 20 (`F_flat - F_tr = ...`) contains
  `180*delta^2*r*xi`, `162*delta^3`, and `220*r^2*xi^4`, matching the
  corrected P_1 and P_2 coefficients in notes lines 132-133 and 137.
- The SymPy script (unchanged) has used the script-side `162/162/180/162/220`
  values since at least the batch III.1 v1 baseline; its
  `expect_zero("dF_tr/dR formula", ...)` and `expect_zero("F_flat - F_tr formula", ...)`
  assertions against those typed polynomials already pass in the standing
  saved output.

No collateral edits in the notes file (the diff is exactly the five intended
replacements). No script change means there is no new assertion that could
have been made tautological, and no engine output that could have drifted.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy exec log was generated for this unit in batch
III.1 v2. This is correct per the Applied block (`Scripts unchanged ... no
script re-run needed`). The standing saved output
`scripts/output/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.txt`
already reflects the script-side coefficients that the notes have now been
corrected to match; the script that produced it is unchanged from when the
auditor cited it as evidence.

**Mathematica:** exit=n/a. Same reasoning. The standing saved output
`mathematica/output/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.txt`
is the independent-engine evidence that the corrected notes coefficients are
what `D[fTr, r]` (line 14) and `Factor[fFlat - fTr]` (line 20) actually
evaluate to.

**Output freshness:** Not applicable. No script was re-run in this batch
because no script was modified. The pre-existing outputs continue to apply
because the scripts they were generated from are byte-identical to the audit
baseline (confirmed by empty `git diff HEAD --` against both audit script
files).

## Material-change assessment

`material_change`: false.

The fix is a notes-only typographical correction. No script, no derived
script output, no paper.tex line, and no downstream stage's symbol
definitions were altered. Stages > 046 consume the boxed monotonicity and
residual-theorem claims from `paper/stages/stage_046.tex`, which were never
quoted in the disputed coefficient form; the notes correction does not
change any consumable result. No `upstream_stale` propagation is warranted.

## Side observations (non-blocking)

- The captured diff at `redteam/exec_logs/stage_046_diff.patch` is from an
  earlier batch III.1 v1 capture window and shows the prior script-side
  edits (the v1 Mathematica transliteration rewrite and v1 SymPy
  insufficient_verification additions), not the v2 notes-only edit. This is
  expected: the v2 flow did not invoke Codex on this unit (notes-only,
  user-resolved), so the orchestrator did not regenerate a `.patch` for the
  notes correction. The authoritative record of the v2 fix is the live
  `git diff HEAD --` on the notes file, which is what this verification
  checks against. Not blocking.
- The directive's `## Applied: F1` block correctly records
  `files_changed: notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:90,132-133,137`,
  `summary: Fixed 5 coefficient typos in notes auxiliary polynomials P_R
  (230->162 twice), P_1 (248->180, 230->162), P_2 (237->220). Scripts
  unchanged (already correct per both engines' independent derivation). Per
  user-approved Q4 (a) in batch III.1 v2.`, and `deviation: none — notes-only
  fix; no script re-run needed.` This matches the on-disk state.
- Pre-existing SymPy script header label `Stage 29` (sympy_audit.py:3 and
  banner) was already noted as a non-blocking labelling inconsistency in the
  v1 verification; it remains pre-existing and outside the scope of this
  v2 finding.

## Verdict justification

The single v2 finding (F1, paper_misalignment notes_contradicts_script) was
resolved by the user-approved direction (a): correct the notes to match the
script-side coefficients that both engines independently derive. All five
prescribed substitutions (P_R `230` -> `162` twice on line 90; P_1
`248` -> `180` on line 132; P_1 `230` -> `162` on line 133; P_2
`237` -> `220` on line 137) are present at the expected lines; no other
notes edits crept in; the audit scripts were not touched (correctly — there
was nothing to fix on the script side); and the standing Mathematica saved
output independently corroborates that the corrected coefficients are what
`D[fTr, r]` and `Factor[fFlat - fTr]` actually evaluate to. `verdict:
verified`, `material_change: false`.
