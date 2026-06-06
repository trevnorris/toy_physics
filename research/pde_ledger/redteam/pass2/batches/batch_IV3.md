# Pass-2 Batch IV.3 (stages 115–126) — summary

**Part IV.3 — Core balance, DtN mixed, outlet, positive source.** 12/12 verified, all
`material_change: false` (NO downstream staling), 0 stop-cold/blocked, all Codex iter-1 exit 0
(no iter-2). **NO checkpoints in range.** Status-only by design: **120 & 124** (notes-only, both
engines null, like 103/113). Value reconciliation: **108 deliverable values batch-wide, 1 misaligned**
(the r_F1 appendix surd, resolved — see below) — 115=9, 116=9, 117=16, 118=15, 119=6, 120=6, 121=7,
122=12(1), 123=6, 124=7, 125=6, 126=9.

## Disposition
- **8 script-side clean:** 115, 116, 119, 121, 122, 123, 125, 126 (all dual-engine, `.wl` independence
  confirmed; value-reconciliation 0-misaligned on the script side).
- **2 status-only notes-only:** 120, 124 (no scripts — by design; carry-forward values consistent).
- **2 script-side findings, both resolved:** 117 (re-author), 118 (de-taut).
- **1 published-paper value typo, resolved by Codex:** r_F1 `√(4107−117π²)` → `√(4107−100π²)`.

## 117 — F1 `mathematica_transliteration` → FULL re-author (USER-AUTHORIZED)
Orchestrator ground-truth `.wl`-vs-`.py` read found the `.wl` was a **full line-by-line `Series[]`-based
transliteration of the `.py` across ALL SIX sections** (the audit agent flagged it but scoped the directive
to §5 only; the orchestrator read broadened it). Every §1–§4 class check `Series`-expanded the same rational
generating function and matched the same canonical-even coefficients; §5 re-typed the reduced shell/mixed
core closed forms `rC/rhoC/sigmaC/kappaC/gammaC`. Same defect class as IV.2's 107/110/114 and IV.1's 100.
**User AUTHORIZED re-author** (re-author-vs-accept = USER-LEVEL; surfaced, not reversed unilaterally).
Directive rewritten to a full re-author (requirement+acceptance only; Codex designed the route per
[[feedback-claude-reviews-codex-codes]]). Codex (session 019e9de7):
- **§1–§4** → an **undetermined-coefficient solve** (`jetFromBalance[tag,den,num]`: posit `Y=Σc_k z^k`,
  impose `den·Y−num==0` order-by-order via `Coefficient`/`Solve`, read the canonical-even / odd-norm
  conditions off the derived `c_k`). §3/§4 additionally clear the pole denominator (a different normalization
  than the `.py`'s raw-Λ `−L2/L0` scheme). NO `Series[ratio]` survives.
- **§5** → an explicit **2×2 `coreMatrix={{ks,lam},{lam,−kq dW}}` elimination** (`Solve` the linear system
  `coreMatrix·{sCore,qCore}==coreSource`, form the Schur complement `gs·sCore+gq·qCore`), with `rhoC/sigmaC`
  DERIVED from the elimination + a hard `coreSchurResidual===0` guard tying the named forms back to the raw
  elimination. The re-typed `sigmaC` literal is GONE.
- Only remaining `Series` = the permitted final §5 residual check (`deltaCore − deltaCoreExpected`).
**Committed `.wl` output BYTE-IDENTICAL to HEAD** (method-only; all 12 `expectZero` PASS; deliverables —
scale/argument solutions, the four class matches, core-balance branches, σ_*, kappa0FromTube, the concrete
core collapse, classification rows, survivor set — unchanged). The SymPy `.py` reference engine UNTOUCHED.
`material_change: false`. **117 ADDED to the Independent-Mirror Set.**

## 118 — F1 `tautological_check` → de-taut (routine, Claude+Codex math-coverage)
The "K_q closed form" check was X−X: `K_q` defined directly as `(Zq/mu0)·π²c_s²/(4 L_W²)` then checked vs the
same literal. Tied `K_q`/`kQ` to the independently-computed gradient integral `chi_grad`/`chiGrad`
(`∫(χ')²dz`, asserted `=π²/(4L_W²)` at the "D/N stiffness check") so the closed-form check is now load-bearing
(reduces to `(Zq/mu0)c_s²·(chi_grad−π²/(4L_W²))`, fails if the integral were wrong). Printed `K_q` value
UNCHANGED; both engines exit 0. `material_change: false`. (Sibling closed-form checks `g_q/J_s/g_s/I_q/K_s`
already genuine.)

## Published-paper r_F1 typo — `117π²` → `100π²` (Codex-applied, Claude-reviewed)
Flagged by 122 F1 / 123 F1 as `paper_misalignment` (value_mismatch). The Family-1 geometric radius surd was
published as `r_F1 = √(4107−117π²)/(10π)` in TWO paper files, but `100π²` is arithmetically forced
(`r_geom=√((12/π²)(L/a)²−1)`, `L/a=37/20` ⟹ `12·(37/20)²=4107/100`), is corroborated by the paper's OWN
adjacent numeric `≈1.77799353547498` (the `117π²` form gives ≈1.7295) and downstream `g_-^{F1}≈0.758035…`,
and matches EVERY script and note (stages 121/122/123/126/127/142/143/148). **User authorized "get everywhere
in the documents to match what the scripts derived."** Codex (out-of-band session 019e9dea, separate from the
scripts-only red-team loop per codex.md) changed `117`→`100` in:
- `paper/appendices/stage_appendix_part04.tex:562`
- `paper/parts/part04_geometry_retarded_mouth.tex:576`
Reviewed: only the two surds changed; adjacent numerics/prose untouched. **Pass-2 caught a published-paper
arithmetic typo the first pass missed.** The scripts were always correct → 121/122/123 are script-clean.

## INFRA
- **20 exec exit 0** (10 dual-engine × 2). Orchestrator independent re-run = the reliability/determinism gate:
  **116 (the historically flaky Mathematica stage) ran clean/deterministic**; **117's re-authored `.wl` exit 0**.
- `$RT exec-*` writes `exec_logs/` only → orchestrator sed-refreshed every committed `.txt`. 117 + 118 outputs
  + 7 clean dual-engine outputs byte-identical post-refresh; **116 NORMALIZED** (refresh stripped a stray
  `# exit_code: 0` trailer a prior commit had baked into BOTH committed `.txt` — deliverables identical, NOT
  a math change; same class as IV.2's 108).
- **Arbiter grep on all 10 committed outputs: CLEAN of stale self-epoch (−17 band 098–109) self-banners.**
  Only hit = 116's `gamma0_bare (upstream-carried input, Stage 98)` — a script-embedded numbering CROSS-ref
  (γ₀ provenance), NOT a self-banner → DEFERRED (content-keyed, never offset-sweep).
- **Seat policy held:** 117 + 118 = 2 `.wl`-touching Codex sessions (concurrent, at the 2-seat cap, flock-safe);
  the r_F1 paper fix = a separate out-of-band 0-seat Codex session (no `math -script`); orchestrator exec
  sequential, after all Codex done (no overlap).

## Deferred (NOT this loop — numbering pass / PAPER_CLEANUP; content-keyed, never offset-sweep)
- **DISCOVERY:** `\stagefield{Purpose}` card self-labels drift +17 — 117 "Stage 134", 119 "Stage 136",
  120 "Stage 137", 124 "Stage 141" — a class the numbering reconciliation MISSED (its scan keyed on
  `\section`/`\label`, not `\stagefield{Purpose}`). Parallels the I.2 bare-`Print` miss.
- 122/123 cards say "Mathematica audit: none yet" though the retro-sweep `.wl` now exists (status understatement).
- 121 `.wl` `stage99TubeLength` (−17 pre-renumber label); 116 source "Stage 98" γ₀-provenance cross-ref.
- **Ansatz catalog:** γ₀=(1+r_c)/9 re-confirmed a POSTULATED pure-scale ANSATZ (not derived; 117 §5 / 116 note)
  — already the catalog's first entry.

## Independence outcome
1 newly-independent `.wl` (117 re-authored); **0 sanctioned mirrors**; 117 added to the Independent-Mirror Set.
121/122/123 retro-sweep `.wl` re-confirmed genuinely independent (Reduce/Solve-with-branch-guard, distinct
from the `.py`; all use the correct `100π²`; 123's negative `Xi_v(F1)≈−1.01675633282526` preserved). 119
checked (algebraic-Solve, distinct from 117's Series-port) — independent-enough. NO checkpoint constant changed
(no checkpoints in range). Pass-1 `MANIFEST.yaml` untouched (isolation held).
