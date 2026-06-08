# Pass-2 Batch IV.5 (stages 139–150) — summary

**Part IV.5 — Family-1 actual mouth gains, susceptibility, profile.** 12/12 verified, all
`material_change: false` (NO downstream staling), 0 stop-cold/blocked, all Codex iter-1 exit 0
(no iter-2), 0 Codex deviations. **NO checkpoints in range. NO EM-projected stages.** Status-only by
design (no engines): **141, 145, 149** (consolidation/carry-forward cards; all values trace to
upstream dual-engine stages 137–148 — same class as 103/113/120/124). Value reconciliation: **0
misaligned batch-wide; NO `paper_misalignment` → ZERO paper/notes edits.**

## Disposition
- **7 dual-engine script-side clean:** 140, 142, 143, 144, 147, 148, 150 (`.wl` independence
  orchestrator-confirmed; value-reconciliation 0-misaligned).
- **3 status-only:** 141, 145, 149 (no scripts — by design; carry-forward values consistent,
  upstream-verified in both engines at 137–148).
- **2 script-side findings, both resolved:** 139 (re-author from a numerical mirror), 146 (de-taut
  in both engines).

## 139 — F1 `mathematica_transliteration` → FULL re-author (USER-AUTHORIZED; ORCHESTRATOR OVERTURN)
**The standout.** The audit agent called 139 **CLEAN** (claiming "no meaningful independent route
exists"). The orchestrator's ground-truth `.wl`-vs-`.py` read OVERTURNED that — the `.wl` was a
NUMERICAL mirror of the `.py`: it imported the `Π_*=1.50882951349316` / `S_q=0.658075937605429`
literals + a hardcoded Stage-134 kernel closed form, and carried a self-described **"sanctioned mirror"**
comment (pass-1 / BATCH-5 had ACCEPTED it). The orchestrator surfaced the user-level re-author-vs-accept
call (re-author-vs-accept = USER-LEVEL; surfaced, not reversed unilaterally; same class as IV.1-100 /
IV.3-117). **User AUTHORIZED re-author.** Codex (iter-1) re-authored the `.wl` to independence:
- **Π_*** now derived via a cleared-denominator bracketed `FindRoot` on `g(Π)=g_minus` (correct sign
  `(e^p−1)`, NOT the §131-F3 `(1−e^p)` trap; residual self-check ~1.6e-58).
- **S_q(Π_*)** now derived via an independent source-moment quadrature `∫Σ·K_q` (`NIntegrate`),
  cross-checked against a symbolic `Integrate` route (residual 0).
- The hardcoded Stage-134 kernel closed form is **GONE** (the two literals `Π_*` / `S_q` survive ONLY
  as `expectApprox` cross-check targets); the "sanctioned mirror" comment is removed.
The SymPy `.py` is UNCHANGED (reference engine). All deliverables preserved at 1e-12 (output gained
precision + 4 new independence checks). `material_change: false`. **139 ADDED to the
Independent-Mirror Set.**

## 146 — F1 `tautological_check` / `insufficient_verification` (BOTH engines) → de-taut
**The audit-agent "transliteration" headline was DOWNGRADED by the orchestrator.** 146's core
moment-formula verifications (`py:33-53` / `wl:44-51` — each engine's OWN `Integrate` of Σ·cos and
Σ·K_q vs the closed forms) are independently integrated per engine and were **LEFT INTACT**, so the
`.wl` is NOT a transliteration. The real defect was a shared verification-strength flaw: the
convex-family affine moment-law checks did not exercise the **ε-slope deliverable**:
- g-side subtracted `g_minus` → collapsed to `(1−ε)(g(Π_*)−g_minus)` = root-closure, mislabeled
  "affine law".
- S-side subtracted `Sformula(Π_*)` → the ε-slope cancelled → effectively x−x.

Fixed SYMMETRICALLY in BOTH engines: closed-form intercepts `g_*=gPi(Π_*)` / `S_*=Sformula(Π_*)`;
deformed moment by ONE direct quadrature of the assembled profile `Σ_ε`; own-quadrature `gbar_v` /
`Sbar_v`; non-vacuity slope guards `|gbar_v−g_*|>1e-3` (~0.0936) & `|Sbar_v−S_*|>1e-3` (~0.0969);
residuals `<1e-25` at ε=1/10 and 1/2; honest relabeling. The residual now reduces to
`(1−ε)(quadrature − closed_form at Π_*)` — a falsifiable kernel identity (a wrong gFormula / sFormula /
K_q now FAILS). No deliverable value moved; ZERO paper/notes edits. `material_change: false`. 146's
`.wl` confirmed STILL independent (core integrals per engine).

## Independence (orchestrator ground-truth `.wl`-vs-`.py` read on ALL 9 dual-engine stages)
**Backstop on all 9 (139, 140, 142, 143, 144, 146, 147, 148, 150). The 100/129/134 lesson
reconfirmed: audit agents UNDER-call transliteration; the orchestrator ground-truth read is the
backstop; re-author-vs-accept is USER-LEVEL.**
- 139 = the standout — audit CLEAN, orchestrator OVERTURNED → user-authorized re-author (above). One
  genuine transliteration, now re-authored.
- 146 = audit "transliteration" headline DOWNGRADED by the orchestrator (core moment formulas ARE
  independently integrated per engine; the real defect was a shared verification-strength flaw, fixed
  in both — above).
- 140 / 142 / 143 / 144 / 147 / 148 / 150 confirmed genuinely independent:
  - **144** — Π_* sign `(e^p−1)` correct + residual self-check.
  - **142** — anchors at Stage-131's Π_* via its OWN projection integral.
  - **147** — A_T cross-checked vs in-engine `D[]` autodiff + per-engine quadrature.
  - **148** — cross-engine dT AGREE ~16 figures (the first-pass live bug stays fixed).
  - **143** — global `Reduce` positivity.
  - **150** — per-engine symbolic identity check.
- **0 sanctioned mirrors remain in IV.5** — all 9 dual-engine stages are now genuinely independent;
  1 newly-independent `.wl` (139 re-authored from the numerical mirror).

## Value reconciliation (pass-2 augmentation)
**0 misaligned batch-wide; NO `paper_misalignment` → ZERO paper/notes edits.** No stale `168π²`/`168%`
anywhere (only interior mantissa digits); the canonical Family-1 radius `√(4107−100π²)/(10π)` is used
throughout.

## INFRA
- **4 orchestrator exec runs exit 0** (reliability/determinism gate): 139 sympy byte-identical (`.py`
  unchanged); 139 mma + 146 sympy + 146 mma refreshed from `exec_logs/` via sed.
- `$RT exec-*` writes `exec_logs/` only → orchestrator refreshed the committed `.txt` from the
  authoritative run.
- **Arbiter grep on committed outputs CLEAN of stale self-epoch (NNN−17 = 122–133 band) banners.**
- **Seat policy held:** 139 ∥ 146 = two `.wl`-touching Codex sessions run concurrently at the 2-seat
  cap (flock-safe); orchestrator exec sequential AFTER both Codex done (no 3rd seat).
- Pass-1 `MANIFEST.yaml` untouched (isolation held). 6 trackers synced (PAPER_CLEANUP **P5-14** = ZERO
  paper/notes edits).

## Deferred (NOT this loop)
- **NONE new for IV.5.** The IV.5 cards are numbering-clean (no +17 `\stagefield{Purpose}` self-label
  drift; all titles canonical — UNLIKE IV.3 which had it at 117/119/120/124). The deferred cross-refs
  are all LEGIT upstream citations (NOT self-epoch, already-correct): 148 "Stage 126 closed form"
  (ξ_*), 144 "stage-131 owning value" (Π_*), 142 "Stage 130 §1" (mouth-source law) — no action.

## Independence outcome
1 newly-independent `.wl` (139 re-authored from a numerical mirror → `FindRoot` Π_* + source-moment
quadrature S_q); **0 sanctioned mirrors**; 139 added to the Independent-Mirror Set. All 9 dual-engine
IV.5 stages now genuinely independent. NO checkpoint constant changed (no checkpoints in range). Pass-1
`MANIFEST.yaml` untouched.
