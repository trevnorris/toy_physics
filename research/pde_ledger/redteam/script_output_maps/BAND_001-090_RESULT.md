# Band 001-090 — RESULT (SCRIPT/OUTPUT numbering pass)

Status: COMPLETE, verified, **uncommitted** (pending user gate). Content-keyed, never
offset-swept. All edits digit-only (strip-digits proof). Both engines exit 0 on every re-run.

## Applied (label-only, orchestrator/agent-applied directly per the pass mandate)
- **172 per-token source edits** (`apply_script_output_band.py`): 18 A2-SELF self-labels
  + 154 A3-CROSS cross-refs (incl. the 2 user-approved author-mislabel fixes: stage036
  `Stage-18`→034 N₋(x), stage089/090 `Stage-075`→080/079 zeta_max ceiling).
- **20 compound trailing-number fixes** (`apply_script_output_compounds.py`): scan caught
  only the first number of `Stage-X/Y` compounds; completed each trailing number
  (e.g. `Stage-003/4`→`Stage-003/021`, `Stage-078/63/64`→`Stage-078/080/081`).
- **6 miss fixes** (`apply_script_output_misses.py`): label tokens the per-line classifier
  flagged on some lines but missed on others / coexisting with a VAR identifier on one line
  (stage040 `Stage18`/`Stage19` expect_zero labels→035/036, `Stage 22`→039; stage045
  `Stage-27`→044). Identifiers (`F_stage18`/`fStage18`) left untouched per user rule.
- **Output refresh** (`rerun_refresh.py`): re-ran 44 `.py` + 19 `.wl` (+2 `.py`/1 `.wl`
  for the misses) → 40 committed `.txt` outputs refreshed; 26 scripts produced
  byte-identical output (comment-only edits — idempotent confirmation).

## Width convention applied (empirically confirmed against canonical files)
`.py` self-labels UNPADDED (`FINAL STAGE-21 LEDGER`, `Stage 30 SymPy audit`); `.py`
docstring-filenames + `.py` cross-refs 3-digit; `.wl` everything 3-digit.

## LEAVE + log (3 items; not stale or out-of-scope-to-fix-as-label-only)
1. **stage021 `.py:3`** docstring old filename `...stage4_maxwell_mixed...` — old STEM ≠
   current stem `reduced_one_port_normal_form`; a number-only edit would mint a non-existent
   filename. Whole-filename rewrite is a content edit, out of scope.
2. **stage047 `.py:121`** neg-control comment "the exact Stage-30 support-loading
   coefficient" — M_supp is *defined* upstream (stage043); non-printed comment; prior pass
   already left it ambiguous. (Owner 043 vs 048 genuinely unsettled.)
3. **stage087 `.py:12-13`** "(post-renumber) … a.k.a former stage 65/69" — intentional
   historical provenance; the author's old numbers may be the true pre-realignment numbers.

## Genuine current-canonical refs LEFT (not stale; cosmetic 2-digit, padding out of scope)
stage021:9 `Stage-3` (=003), stage078:26/32 `Stage-77`/`Stage-75` (=077/075, the live
old-form-vs-new-form two-epoch case), plus all `.py` 2-digit self-banners (the convention).

## Verification
- 101 files changed (61 source: 44 `.py` + 17 `.wl`; 40 outputs `.txt`).
- **strip-digits label-only proof: 101/101 digit-only** (HEAD vs working byte-identical
  after removing all digits → zero equation/value/variable/punctuation/logic bytes changed).
- Residual scan: **0** remaining cited≠canon capitalized 1-2-digit refs (excl. the 3 LEAVEs).
- 0 half-fixed compounds remain. No notes/`.tex` touched. New tooling files clean of the wrong-root path typo.
- 66/66 script re-runs exit 0.

## Pre-existing issue flagged (RESOLVED before the band-1 commit)
A grep for the wrong-root path typo found 6 PRE-EXISTING occurrences in redteam prose
artifacts (REMEDIATION_HANDOFF.md, pass2/batches/batch_VI1.md,
pass2/reports/stage_234|248|251.md, pass2/verifications/stage_216.md) — crept in during
later pass-2 batches after the d3bbd41 purge. Reworded to clear the grep (commit 14603bb);
the repo is clean of the typo string.
