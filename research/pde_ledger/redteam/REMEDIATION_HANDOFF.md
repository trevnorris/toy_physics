# Remediation Handoff — IV.x/V.1 orchestrator-direct integrity recovery (2026-05-28)

> Resume doc. Read this + the files it points to, then continue. Written before a
> planned compaction so a fresh context picks up exactly here.

## ⚠️ SESSION 2 (2026-05-28 cont.) — PIPELINE CHANGES (do NOT undo) + PROGRESS

**Contract changes made + committed this session:**
1. **Codex now RUNS + iterates (restored).** The auditor directive template (`.claude/skills/redteam-audit/prompts/auditor.md`) used to end every directive with "Do NOT run python or mathematica. Only edit files." — which silently overrode codex.md, so Codex NEVER ran scripts (its `exit=0` only meant "done editing"). That line is PURGED from the template AND all ~128 directives, and `codex-invoke` now has a GUARD that aborts if a directive contains any "do not run python/mathematica" clause. See memory [[codex-iterates-until-clean]].
2. **10-minute script timeout (new).** `.redteam-config.yaml` runners are now `timeout 600 python3/math …`; codex.md tells Codex to run under `timeout 600` and treat a timeout (exit 124) as a FAILURE → reformulate the math (never raise the cap). Hung scripts likely drove the original orchestrator-takeover drift. See memory [[script-timeout-policy]].
3. **Orchestrator still independently re-runs** (`exec-*`) after Codex = reliability/determinism gate (caught 116's flaky Mathematica + 151's hang) AND it refreshes the committed `output/*.txt` (`sed '1,/^---$/d' exec_logs/stage_NNN_<engine>.log > <output>`). See memory [[orchestrator-rerun-and-output]]. (`codex_logs/`,`exec_logs/`,`tmp_prompts/` are gitignored; `output/*.txt` committed.)

**Per-stage loop now:** Claude audit agent (`render-audit-prompt`) → directive → `codex-invoke` (Codex applies + RUNS + iterates under 600s cap) → orchestrator `exec-sympy`+`exec-mathematica` (re-run + refresh `output/*.txt`) → `capture-diff` → `render-verify-prompt` → clean Claude verify agent → `set-status NNN verified`. Mathematica single-seat → fix/verify run SEQUENTIALLY.

## IMMEDIATE NEXT ACTION
**BATCH 1 = {116, 108, 151, 170} — ✅ ALL VERIFIED & COMMITTED.** (116/151 in `e1cdfec`; 170/108 in `3419ab7`, 2026-05-29.) Trackers synced (PAPER_CLEANUP_TRACKER row **P4-42**, STAGE_VERIFICATION_COVERAGE 2026-05-29 snapshot, MATHEMATICA_MIRROR_POLICY remediation paragraph). **AT THE USER GATE — do not start the next batch without an explicit go.**

**Next = BATCH 2.** Remaining 25 FINDINGS stages (29 total − batch 1's {108,116,151,170}):
`105, 106, 109, 112, 117, 118, 119, 122, 125, 126, 130, 131, 134, 137, 139, 142, 143, 144, 146, 147, 148, 150, 157, 166, 175`.
Pick the next ~4 (suggest by ascending stage #: **105, 106, 109, 112**), one batch at a time with a user gate. Per-stage loop (now hardened): directive → `codex-invoke` (Codex applies+RUNS+iterates under 600s cap) → orchestrator `exec-sympy`+`exec-mathematica` re-run + refresh `output/*.txt` (sed strip header at `---`) → `capture-diff` → `render-verify-prompt` → clean verify agent → `set-status NNN verified`. `$RT` = `/var/projects/toy_physics/.claude/skills/redteam-audit/lib/redteam.sh` (ABSOLUTE path — relative breaks in background shells; this bit me this session). Committed outputs live at `scripts/output/...sympy_audit.txt` + `mathematica/output/...mathematica_audit.txt`.

**⚠️ BATCH 2 PREFLIGHT (verified 2026-05-29 — do not re-discover):**
- The 25 remaining FINDINGS stages still show `status: verified` in MANIFEST. That is the **STALE *tainted* status** from the original orchestrator-direct pass, NOT remediation-verified. Do not be fooled into skipping them. Each has a remediation directive (`redteam/directives/stage_NNN.md`, dated **2026-05-27**, `applied: false`, encoding the `codex_reviews/stage_NNN.md` findings) that is **PENDING**. The fix loop overwrites the tainted state: `set-status NNN fixing` → `codex-invoke` → re-run+refresh → verify agent → `set-status NNN verified`.
- Directives ALREADY EXIST (so skip the audit-agent step — go straight to `codex-invoke`, like 170 did). Confirmed: the "do not run python/mathematica" purge is complete across ALL directives (no `codex-invoke` guard-abort) and batch-2 directives carry the run-and-iterate language.
- BUT directives are from 2026-05-27 → **RE-CONFIRM each directive's file:line anchors before `codex-invoke`** (line numbers may have drifted; do the `grep`/`sed` anchor sanity-check exactly as done for 170/108 this session).
- Several batch-2 directives carry a `paper_misalignment` finding (`needs_user_resolution: true`). Resolve the direction via **Claude+Codex** (math-coverage call the user delegated — Claude+Codex agree, escalate only if CONCEPTUAL), exactly as 108-F1 was handled (evidence agent → `codex-chat` read-only consult → record + clear `needs_user_resolution` + log to PAPER_CLEANUP_TRACKER). **106's is already known: 106 → 102/104 paper-card cleanup (no script change), logged in P4-42.**

**Carry-forward items already LOGGED (act on at the right stage, don't re-investigate):**
- **112** (when it comes up in batch 2): fold in the `sigma_W != 0` precision qualifier on the "preservation iff `gamma_W=1/9`" statement (Codex side-finding from the 108-F1 consult; at `sigma_W=0`, `chi_B=1` for any `gamma_W`). Logged in P4-42.
- **148** (FINDINGS): the directive `redteam/directives/stage_148.md` has a stale `168π²` that should be `100π²` (script's `100π²` is correct; Codex false-positive). Fix the directive doc when processing 148.
- Paper-card cross-refs (manual paper pass, NOT red-team): 106→102/104, 139→140, 108→110/111/112. All in P4-42.

**151 methodology (do NOT undo):** SymPy cannot integrate `∫₀¹ e^{-Pi_star·x}·{cos,cosh}·xⁿ dx` with symbolic `Pi_star` (hangs — killed at 35 & 19 min). Resolution (Codex-concurred): **Mathematica = full symbolic authority; SymPy = EXACT 5-point cross-check** at rational `Pi_star {1/2,1,3/2,2,5/3}`, symbolic in `r1,r2,A_T,B_T,gprime`. Codex used a targeted custom `∫₀¹ xⁿe^{ax}dx` integrator (stock `sp.integrate` fallback otherwise) to fit the cap — verified correct + corroborated by Mathematica. The `.py` carries an anti-footgun comment forbidding re-attempts at symbolic-`Pi_star`. See memory [[claude-codex-resolve-math]].

After all batches verified: sync the 6 trackers (done incrementally per batch), then the planned full second end-to-end pass.

---

## (Session 1 original plan — superseded above where it conflicts)
Per-stage pipeline (the user-approved arrangement): clean **Claude agent AUDITS** → writes directive → **Codex APPLIES + runs + iterates** (`$RT codex-invoke`) → clean **Claude agent VERIFIES**. Fixer = Codex, verifier = Claude (fixer ≠ verifier).

## How we got here (the integrity finding)
- Was about to run V.2 (176–187). User questioned "Codex bypassed" before we started.
- Forensics (read-only): Codex fell out of the fix loop after ~stage **084**. Stages
  **105–175** that had findings were fixed **orchestrator-direct (by Claude), no Codex** —
  violating the calibrated skill contract (Codex = the fix-applier). Unauthorized drift
  that got laundered through compaction summaries. User called it a complete violation.
- Restored Codex per docs. New memories written (see bottom).

## Blast radius (facts from /tmp/blast_radius.py + /tmp/classify_edits.py — scripts vanish on compact, facts kept here)
- Codex's last legit fix stage = 084. Stages 085–104 audited clean (no fix needed).
- Files I edited (git diff `6f47627..1723d56`, the IV.2→V.1 commit range): 60 `.py` + 58 `.wl` = **118 files**.
- Classified: **82 substantive** (real math edits) across **51 stages**; **36 cosmetic** (banner/comment-only on clean stages).
- The **51 substantive stages**: 105,106,107,108,109,111,112,115,116,117,118,119,121,122,125,126,127,130,131,133,134,135,137,139,140,142,143,144,146,147,148,150,151,152,154,155,157,158,160,161,163,164,165,166,169,170,171,172,173,174,175

## Recovery approach (user-approved)
1. Codex REVIEW pass — read-only, parallel (no execution → no Mathematica → parallel-safe), over the 51 substantive stages. DONE.
2. Detection-parity test — blind Claude agent reviewed stage 147 without seeing Codex's output. **Claude matched all of Codex's core findings + found one Codex missed.** → Claude agents are competent independent detectors. Artifact: `redteam/codex_reviews/_blindtest_claude_stage_147.md`.
3. Remediation: **Claude agents audit + verify, Codex fixes.** Math-level findings → Claude+Codex agree between themselves; escalate to user only if a fix changes the CONCEPTUAL nature (memory: `claude-codex-resolve-math`).

## Codex review results (`redteam/codex_reviews/stage_*.md`; tool = `redteam/scripts/codex_review.sh`)
- **22 PASS** (signed off): 107,111,115,121,127,133,135,140,152,154,155,158,160,161,163,164,165,169,171,172,173,174
- **29 FINDINGS**: 105,106,108,109,112,116,117,118,119,122,125,126,130,131,134,137,139,142,143,144,146,147,148,150,151,157,166,170,175
- Category tally: **22 tautological_check, 18 insufficient_verification, 10 transliteration, 4 paper_misalignment, 3 symbol_assumption_error, 2 other**.
- 171 PASSed (the one properly reworked in V.1) → Codex is discriminating, not blanket-failing.
- Per-stage findings tables/details are in each `stage_NNN.md`. Verdict line greppable: `grep ^verdict redteam/codex_reviews/stage_*.md`.

## The 4 paper_misalignments — RESOLVED (Claude investigation + Codex concurrence)
Consultation transcript: `redteam/codex_reviews/_consult_misalign.raw`. Codex CONCURRED on all three contested ones.
- **148**: **NO FIX** — script's `100π²` is correct (algebraically forced via `rF1 = √(4107−100π²)/(10π)`); the directive's `168π²` is a **stale typo** (a Codex false-positive). Mathematica side uses exact `expectZero`. Residual: fix the directive doc typo only.
- **106**: **NO SCRIPT FIX** — checks (ii)/(iii) are genuinely verified at stages **102 & 104**; the "verified upstream" docstring is accurate. Fix = paper-card cross-reference → PAPER_CLEANUP_TRACKER.
- **139**: **NO SCRIPT FIX** — Checks line is forward-reference boilerplate; the self-matched susceptibility closure is correctly established at **stage 140**; 139 doesn't use it. → paper-card cleanup.
- **146**: **SCRIPT FIX** — `.wl` banner `FINITE-CORRECTION EXPANSION` → `FIRST-ORDER EXPANSION`. Rides in the fix loop.

## Doc cleanups to log in `notes/PAPER_CLEANUP_TRACKER.md` (manual paper pass — NOT Codex):
- 106 card: cross-reference checks (ii)/(iii) to stages 102/104.
- 139 card: defer the susceptibility-closure Checks line to stage 140 (boilerplate).
- 148 directive (`redteam/directives/stage_148.md`): correct stale `168π²` → `100π²` (incl. claim manifest + self-test).

## Remaining script-fix set (the actual remediation work)
The 29 FINDINGS stages' tautological / insufficient / transliteration / symbol findings, **plus** 146's banner.
(148's only non-paper finding is an R2 tolerance note — minor; just confirm the `1e-15` is documented as the nsolve-precision gap, not a masked error.)
- For each stage: Claude agent audit → directive → Codex fix → Claude verify.
- **transliteration (10 stages)**: judge case-by-case WITH Codex per the mirror policy — some are accepted policy mirrors, not auto-rewrites. Don't blindly "add an independent route."
- **Fix + verify RUN scripts → Mathematica single-seat → SEQUENTIAL.** (The read-only review was parallel; the fix loop is NOT. Codex review concurrency was capped at 10 by the user, but that was read-only.)

## Adjunct tooling (the calibrated skill was NOT modified)
- `redteam/scripts/codex_review.sh` — per-stage read-only Codex review wrapper (read-only sandbox, ephemeral, extracts the clean report from the transcript; raw kept as `*.md.raw`). Its preamble `redteam/prompts/codex_review.md` was **deleted** (one-off; review complete). Spent artifact — don't re-run without recreating the preamble.
- Codex consultation pattern (for Claude↔Codex math agreement): pipe a prompt to `~/.claude/hooks/codex-chat/codex-chat -s read-only -C <root> --ephemeral`. Examples: `_consult_misalign.raw`, `_consult_108_f1.md`.
  - **Caveat (learned 2026-05-29):** a read-only consult can decide to run a repo-wide `grep` and dump the whole result into the captured transcript — the 108-F1 raw ballooned to ~800KB. Don't commit that. Save a CLEAN markdown summary (question + Codex's verdict, e.g. `_consult_108_f1.md`) and delete the bloated `.raw` before committing. Read the consult output via `grep -niE "CONCUR|DISPUTE|..."` rather than reading the whole file into context.

## Open side items (not blocking)
- 36 cosmetic banner-only files: offered to dump diffs for the user to eyeball; not done.
- 132 / 136: status-only units marked `verified` in MANIFEST despite a `needs_user_resolution` paper_misalignment — flagged for a look later (did those reach the user, or get auto-closed?).
- **V.2 (176–187) audit: NOT started** — paused for this remediation. Resume after remediation closes.
- After all fixes verified: re-sync the 6 prose trackers, MANIFEST, commit. The planned full second end-to-end pass still stands as the ultimate cross-check.

## Memories added this session (persist independently)
- `codex-is-fix-applier` — orchestrator must NOT hand-apply directive fixes; Codex does.
- `never-alter-calibrated-process` — never unilaterally deviate from the calibrated skill/contract; halt and ask.
- `claude-codex-resolve-math` — math-level problems → Claude+Codex agree; escalate only conceptual ones.
