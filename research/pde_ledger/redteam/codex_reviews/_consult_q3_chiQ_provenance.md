# Consult — Q3 reopened: the chi_Q=1 provenance attribution (Claude ↔ Codex)

**Date:** 2026-06-10
**Mode:** Claude (orchestrator) + Codex read-only consult (`codex exec -s read-only`, gpt-5.5; raw: `_consult_q3_chiQ_provenance.raw`; prompt: `redteam/tmp_prompts/consult_q3_chiQ_provenance.md`).
**Why reopened:** in the first escalation consult Q3 was filed as "LEAVE / semantic — flag for the human author." The user clarified there IS no external author (Claude+Codex wrote the ledger) and they haven't reviewed it, so we must DECIDE it ourselves from content.

## Finding (both engines concur — Codex verdict **B = FIX**)
A carry-forward annotation in stage100/stage101 (scripts + stage100 `.wl`) and the stage101 notes attributes **"the chi_Q = 1 identification by DtN comparison" to stage 097 (single_normalization_defect), calling it UPSTREAM.** That is wrong on two counts, verified by content:
1. **stage097** (`single_normalization_defect`) has NO `chi_Q` symbol/solve — it reduces the actual isotropic passive/outgoing branch to a single normalization defect (`N_Q`). It establishes the branch *context*; it does not identify χ_Q.
2. **stage104** (`outgoing_dtn_fingerprint`) proves the exact outgoing l=2 fingerprint; **stage105** (`chiQ_fix_from_outgoing_dtn`, banner "EXACT FIXING OF chi_Q", `.py:68` solves the ω⁵ match) is where **χ_Q=1 is identified by DtN coefficient comparison**. **stage106** already credits 105 ("Stage 105 … fixes chi_Q = 1 from that fingerprint").
3. 105 > 100/101, so the "upstream" wording is **inconsistent**; and the annotation is **self-contradictory** — it carries χ_Q *free* at stage100 (correct: "would falsely constrain") yet says it's "pinned upstream." χ_Q is free precisely because it is fixed DOWNSTREAM at 105.
4. The stage101 notes' "**Stage 80**" = old-numbered 097 (80+17; canonical stage080 is `family1_zeta_thresholds`, unrelated) — so notes and scripts agree on 097, and both are wrong for the *identification by DtN comparison*.

**Math is UNAFFECTED:** χ_Q is correctly free in 100/101; the closure m̂₀²·χ_Q·N_Q=1 is correctly derived; χ_Q=1 is correctly fixed at 105. This is purely a provenance-wording correction.

## Agreed correction (annotation/notes — comment-only, Codex-supplied exact text)
Credit **stage 104** (exact fingerprint) + **stage 105** (χ_Q=1 by DtN coefficient match), DOWNSTREAM; keep stage 097 only as "establishes the single-defect branch context." Sites:
- `stage100.py:11-14` docstring (iii); `stage100.py:25-27` free-symbol comment.
- `stage100.wl:28-31` carry-forward block; `stage100.wl:34` free-parameter comment.
- `stage101.py:11-15` docstring (3).
- `notes/stages/...stage101...md:40` ("Stage 80 … by definition chi_Q=1" paragraph → 097 reduced the branch; 104 proves the fingerprint; 105 fixes χ_Q by the ω⁵ match).
(Exact replacement strings are in the `.raw` transcript.)

## Additional finding — PAPER-CARD OVERCLAIM (Claude-verified; separate, published surface)
`paper/stages/stage_100.tex:13` and `stage_101.tex:13` (identical "Derivation ledger"):
> "The computation isolates the reduced product m̂₀²χ_Q N_Q=1 **and the canonical condition χ_Q=1**."
The second clause **overclaims** — stage100/101's computation keeps χ_Q *free* and isolates only the *product*; the canonical condition χ_Q=1 is isolated at stage105 (from the stage104 fingerprint). Fix = reword to "keeps χ_Q symbolic; the canonical condition χ_Q=1 is fixed downstream at Stage 105 from the Stage 104 fingerprint." This is a genuine `paper_misalignment` (no math impact), in TWO published cards.

## Who-applies / disposition
All edits are content/wording (no equation/value/logic): scripts + notes + cards → **Codex applies, Claude reviews** (per [[feedback-codex-is-fix-applier]] file ownership).

## ✅ APPLIED 2026-06-10 (user chose "all of it"; committed in this commit)
Codex applied the 7 edits across 6 files (`redteam/codex_logs/apply_q3_chiQ_provenance.txt`):
- `stage100.py` ×2 (docstring (iii) + free-symbol comment), `stage100.wl` ×2 (carry-forward block + free-parameter comment), `stage101.py` ×1 (docstring (3)) — provenance: 097 = branch context; 104 = exact fingerprint; 105 = χ_Q=1 fix (downstream).
- `notes/stages/...stage101...md` — "Stage 80 already fixed…" → 097-reduced + 104-proves + 105-fixes; **math displays preserved** (`Yhat_Q^ret`, `sigma_Q^can`).
- `paper/stages/stage_100.tex:13` + `stage_101.tex:13` — overclaim removed: "keeps χ_Q symbolic; the canonical condition χ_Q=1 is fixed downstream at Stage 105 from the Stage 104 outgoing fingerprint" (+ one inert trailing-blank-line trim in stage_101.tex).
- **Orchestrator verification:** full `git diff` review (every change is prose; all LaTeX math + notes displays intact; residual stale-attribution scan empty); re-ran sympy 100/101 + mathematica 100 → **3/3 exit 0, outputs byte-identical**. Now consistent with stage106's citation of 105.
