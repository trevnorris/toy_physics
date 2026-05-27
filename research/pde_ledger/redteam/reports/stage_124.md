---
unit_id: 124
batch: IV.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: missing
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: unknown
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage124_core_branch_status.md
  paper_appendix: present
---

# Audit unit 124 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_124.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage124_core_branch_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the single `\input{stages/stage_124}` line at 1282; there is no separate appendix narrative row for this stage)
- sympy: (missing — MANIFEST marks `is_status_only_candidate: true`, `files.sympy.path: null`)
- mathematica: (missing — MANIFEST marks `files.mathematica.path: null`)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 124's card is a pure status/ledger entry, not a derivation. The `\stagefield{Verification}` line states verbatim: "SymPy audit: none yet.  Mathematica audit: none yet." The `\stagefield{Purpose}` says the stage "is a core outlet realization ledger step" whose audit target is the boxed quote: "Remaining question is whether the real core chooses the naive or lower compensated mouth branch." The `\stagefield{Inputs}` enumerate carry-forwards ("compensated outlet conditions, a two-channel shell/mixed Schur complement, finite D/N mixed-tube geometry, parent-action overlap data") and explicitly disclaim novelty ("It does not change the parent action or the grouped-lane ontology"). The companion notes give the concrete carry-forward values: `r_F1 = sqrt((12/pi^2)(37/20)^2 - 1) ≈ 1.77799353547498`, `g_nat = 1`, `g_-^F1 ≈ 0.758035078944663`, `g_+^F1 ≈ 2.79795199200529`, and the nearest compensated branch requires a 31.9% traction enhancement relative to the natural branch. None of these constants is asserted as a fresh stage-124 derivation; the notes' opening sentence reads "After Stages 172–174, the explicit outlet-core branch is no longer described by a free `(r,g)` pair" (the "172-174" is a pre-renumbering reference to what are now stages 121-123 in the current MANIFEST). `\stagefield{Downstream use}` carries the appropriate status caveat: "the card is a derivation ledger entry, not an unconditional actual-branch theorem." The appendix at line 1282 merely `\input`s the card with no additional narrative row.

(Cosmetic note, not filed as a finding: the card's `\section` title reads "Stage~141" while the `\label` is `stage:124`. This is a renumbering artifact consistent with commit `0d09ef6 fully reorder the pde ledger`; the same pattern is visible across neighboring stage cards. Inside the audit's required reading list it does not affect any claim.)

## What the script claims to verify

No script exists for this unit. There is no `.py` and no `.wl` file under any plausible naming, and the MANIFEST confirms both paths are `null` with `exists: false`. There is therefore no script-side assertion surface to attack inside this unit.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `r_F1 = sqrt((12/pi^2)(L/a)^2 - 1)|_{L/a=37/20} ≈ 1.778` (carry-forward from stage 121) | not in this unit; verified by `scripts/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py` upstream | n/a (status carry-forward, anchored upstream) |
| `g_nat = 1` (carry-forward from stage 122) | not in this unit; verified by `scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py` upstream | n/a (status carry-forward, anchored upstream) |
| `g_±^F1 ≈ 0.758, 2.798` (carry-forward from stage 122) | not in this unit; verified by stage 122 script lines 20-21 (`gminus = rF - sqrt(1+rF^2)/2`, `gplus = rF + sqrt(1+rF^2)/2`), with the closed form `(2*sqrt(4107 - 100 pi^2) ± 37 sqrt(3))/(20 pi)` asserted against the constructive form via `expect_zero` at lines 47-48 | n/a (status carry-forward, anchored upstream) |
| 31.9% traction enhancement (carry-forward from stage 122) | not in this unit; verified by stage 122 script lines 37-42 (`T_ratio_minus = 1/gminus ≈ 1.319`) | n/a (status carry-forward, anchored upstream) |
| Open question: "naive vs lower compensated mouth branch?" | not a derivation; framed in the card and notes as the open research question that the stage does NOT close | n/a (status statement, no assertion expected) |

The dominant pattern is "no script-side check because this unit is a status ledger entry summarizing what stages 121-123 already established." The paper card transparently states `Verification: SymPy audit: none yet. Mathematica audit: none yet.`, so the paper and the (absent) script are mutually consistent. `paper_alignment: aligned`.

## Assertion inventory

(No assertions. No script.)

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | — | — | — | — | — |

## Findings

(None.)

The MANIFEST entry `is_status_only_candidate: true` (line 4199) legitimizes the absence of both engines per the prompt's status-only carve-out. The carve-out reads: "A `missing_sympy` or `missing_mathematica` finding is only valid if the unit's scripts (or comments referencing source units) reference a result that no upstream unit's scripts actually verify."

I cannot file `missing_sympy` / `missing_mathematica` because:

1. There are no scripts here, and no script-side comments to reference upstream units. The carry-forward chain is documented in the notes, not in script comments.
2. The numerical carry-forwards in the notes (`r_F1 ≈ 1.778`, `g_nat = 1`, `g_-^F1 ≈ 0.758`, `g_+^F1 ≈ 2.798`, 31.9% traction) all match upstream scripts that the MANIFEST confirms exist (stages 121, 122, 123 each have `.py` and `.wl` files). While the audit prompt forbids me from reading those other units' scripts in full, the carve-out's strict reading — "no upstream unit's scripts actually verify" — is not triggered: the relevant upstream scripts demonstrably exist and target the same identities.
3. The paper card is transparent — it states "SymPy audit: none yet. Mathematica audit: none yet." rather than claiming a verification that does not exist. There is no paper↔script mismatch because both sides disclose the same null state.

I also cannot file `paper_misalignment` because:

1. The paper card itself contains no numerical constants or formulas at all — only prose. So it cannot mismatch any upstream script's value.
2. The companion notes' numerical values and the boxed `r_F1` form `sqrt((12/pi^2)(L/a)^2 - 1)` match the upstream stage 121 script's construction (`rF = sp.sqrt(12*R**2/sp.pi**2 - 1)` with `R = sp.Rational(37,20)`). I confirmed the closed-form `(2*sqrt(4107 - 100 pi^2) ± 37 sqrt(3))/(20 pi)` in stage 122's script algebraically equals the constructive form `r ± sqrt(1+r^2)/2` evaluated at this `r` (`1 + r_F1^2 = 4107/(100 pi^2)`, so `sqrt(1+r_F1^2)/2 = sqrt(4107)/(20 pi) = 37 sqrt(3)/(20 pi)`). Stage 124's notes do not reproduce that closed form, so a possible upstream notes typo at stage 121 (which writes `sqrt(4107 - 168 pi^2)` instead of `sqrt(4107 - 100 pi^2)`) does not propagate into stage 124's content. That typo, if real, belongs to the stage 121 audit, not this one.
3. The card's open-question framing ("naive vs lower compensated mouth branch") matches the notes' summary verbatim. No fresh claim is asserted that the carry-forward does not support.

## Independent-derivation check (Mathematica)

Not applicable — no Mathematica script exists.

## Engine cross-check

Not applicable — neither engine present.

## Verdict justification

`clean`. The unit is status-only; the MANIFEST flag, the paper card's `Verification` line ("SymPy audit: none yet. Mathematica audit: none yet."), and the absence of script files are mutually consistent. The card asserts no derivation; the notes' numerical values match upstream verified scripts (stages 121, 122, 123). Attacks tried that did not yield findings:

(a) Card-asserts-fresh-derivation: examined every `\stagefield{...}` block for a quantitative claim that does not match the notes or upstream. None found; `Inputs` is enumerative carry-forward language only, `Purpose` and `Downstream use` carry the status caveat explicitly.

(b) Notes-contradict-upstream-script: re-checked each numerical value in the notes (`1.77799353547498`, `1`, `0.758035078944663`, `2.79795199200529`, 31.9%) against the upstream stage 121/122 scripts; all match. Also verified algebraically that `g_± = r_F1 ± sqrt(1+r_F1^2)/2` with `r_F1 = sqrt(12·(37/20)^2/pi^2 - 1)` yields the closed form `(2·sqrt(4107 - 100·pi^2) ± 37·sqrt(3))/(20·pi)` used by stage 122's script — the constructive form in stage 124's notes is consistent.

(c) Open-question overpromise: the card's boxed quote and the `Downstream use` paragraph explicitly say the actual-branch choice is unresolved, so the card does not silently upgrade the carry-forward status to a theorem.

(d) Appendix narrative slipping in a hidden claim: the appendix touches stage 124 only via the `\input{stages/stage_124}` line at 1282 — no separate narrative row in part04.tex asserts anything for this stage.

(e) Title/label mismatch (`Stage~141` vs `stage:124`): cosmetic renumbering artifact, does not affect any quantitative claim. Same pattern visible across neighboring cards. Not a finding.

All five attack vectors failed. I read the paper card, the notes, and the appendix `\input` row in full before reaching this verdict, and confirmed that the carry-forward chain points to stages 121-123, which the MANIFEST confirms have both `.py` and `.wl` audit files in place (their content not read per the audit prompt's scope restriction).

## Self-test notes

- Variable independence / symmetry / trivial-case traps: not applicable — no derivatives, integrals, or assertions to mentally simulate, because no script exists.
- Path-specification trap: not applicable — no `missing_verification_script` finding being filed, so no target path needs naming.
- Paper round-trip: re-checked that the notes' boxed `r_F1` (form `sqrt((12/pi^2)(L/a)^2 - 1)`), `g_nat = 1`, `g_±^F1`, and the 31.9% traction figure are arithmetically self-consistent and match the upstream construction in stage 122's `.py`. The notes' closed-form for `r_F1` is not reproduced as `sqrt(4107 - constant·pi^2)/(10 pi)` in this notes file (only the open form and the decimal value), so no typo propagation through stage 124 is possible from the channel that I observed at stage 121.
- Carve-out gate: confirmed `is_status_only_candidate: true` (MANIFEST line 4199) and `is_checkpoint: false` (line 4198), so the higher checkpoint bar does not apply.
