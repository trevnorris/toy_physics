---
batch: VII.1
pass: 2
range: 219-230
total_stages: 12
audit_status: complete
closeout_status: complete    # exec re-run + verify + trackers + commit DONE 2026-06-09
script_corrections: 0
findings_count: 5            # ALL card-text-lag paper_misalignment → deferred P4-51 (non-blocking)
needs_user_resolution: 0     # card-text-lag is the STANDING-deferred class, not a new question
checkpoints: [221]           # CLEARED — re-author sufficient (VI.1-218 outcome, not V.3-200)
independence: 12/12 confirmed
material_change: false       # CONFIRMED at close-out (no script edits; all 12 verify material_change:false)
audit_date: 2026-06-09
closeout_date: 2026-06-09
verified: 12/12              # all 12 set-status verified in pass-2 MANIFEST
---

# Red-team batch VII.1 (PASS 2) — Mixed-bundle / resonance / branch-packet / 5PN

## Summary — the FIRST pass-2 batch needing ZERO script corrections

12-stage pass-2 re-audit of VII.1 (219–230) under the v2 paper-grounded auditor
+ the value-reconciliation augmentation. **No `.py` or `.wl` script needs any
correction** — the first pass-2 batch with zero script-side changes. Every prior
pass-2 batch (III.4 … VI.1) carried at least one re-author or de-taut; VII.1 does
not, and the reason is structural (below). All substantive first-pass work was
independently re-confirmed to HOLD. The one checkpoint (221) cleared the higher
bar. The only findings are the standing-deferred card-text-lag class.

**Audit + adjudication are COMPLETE. The mechanical close-out is PENDING (see end).**

### Why this batch is genuinely clean (not an under-call)

VII.1's FIRST pass (2026-06-02) was already executed under the full v2 +
dual-engine regime: it authored the 11 non-checkpoint `.wl` **fresh, from
scratch, as independent routes** (NOT de-transliterated from pre-existing ports),
re-authored the checkpoint 221, and applied the 6 de-taut fixes + 7 notes-typo
corrections in one pass. The earlier bands (IV.x–VI.1) carried `.wl` written under
the older pre-dual-engine / orchestrator-direct regime, so pass-2 kept surfacing
porting residue. V.3 and VI.1 were also retrofit batches, but their lone
insufficiencies (200, 211) were **checkpoint re-authors-FROM-transliteration**;
VII.1's one such re-author (221) was genuinely sufficient, and the other 11 were
never ports to begin with — so there was nothing to catch.

## Independence verification (the durable record — DO NOT re-litigate next session)

Audit agents systematically UNDER-call independence, so the orchestrator did the
ground-truth backstop directly. Result: **12/12 genuinely dual-engine independent,
0 ports, 0 borderline.**

- **Orchestrator ground-truth `.wl`-vs-`.py` read (read in full by the orchestrator):**
  - **221 (checkpoint, scrutinized hardest):** re-author SUFFICIENT (VI.1-218
    outcome, not V.3-200). Decisive: §II derives the load-bearing Stage-220
    derivative identity via native `D[QPi/DeltaPi, portPi]` (wl:96, an **output**)
    then proves it equals the perfect square `(A·G_W+R·G_U)²/Δ_Π²`, whereas the
    `.py` *posits* that form (py:50, an **input**) and only checks the derivative
    matches it. Opposite information flow → genuinely independent. Plus `.wl`-only
    `Residue[]` extraction (wl:78-81) and generic `ComplexExpand` line-shape checks
    (wl:130-131). F1 survival round-trips connect independent expression chains;
    F2 deliverable #9 (linear survival window) is covered in BOTH engines.
  - **220:** independent. Both build the same premise matrix `K_dyn` (unavoidable),
    but each computes `Det`/`Inverse`/`∂_Π` via its OWN native primitives vs the
    posited closed forms (canonical corroboration). `.wl`-only extras the `.py`
    lacks: M6 `CoefficientRules` Laurent-support CLOSURE proof, M7 Π-dependence
    self-test, M9 dissipative-sample guard. NOT the V.3-200 re-typed-algebra signature.
  - **226:** independent. `.wl` recomputes the compiler derivatives, corridor
    matrices, nullspaces, σ-values from the physical primitives natively; pins
    FEWER literals than the `.py` (never re-types matrix/nullspace vectors); M1
    uses `Series` vs `sp.diff`. The `Orthogonalize`-vs-`QR` swap is cosmetic in
    isolation (the subspace projector is unique) but the overall script is a real
    independent recomputation.
- **Calibrated independence-only sweep (focused agent, operation-level + line refs,
  calibrated to the V.3-200 port signature) on the other 9 → ALL independent:**
  - **224** (strong): `.wl` derives the inverse map by `Solve`/`LinearSolve` of the
    forward lane system; `.py` POSITS the inverse coefficients.
  - **227** (strong): `.wl` hand-built log-derivative operator `param·∂/∂param` vs
    `.py` exp-dressing + ε-diff; Λ factorization proved on fresh abstract symbols.
  - **228** (strong): `.wl` EXACT implicit differentiation `dy=−fε/fy` for the
    load-bearing dynamic R_Q/P0 slopes vs `.py` NUMERICAL central finite-difference
    (step 1e-8) — agreement to ~1e-6 then divergence is a signature a port can't fake.
  - **229 / 230** (strong): `.wl` adds `Resolve[ForAll]`/`Reduce` universal-quantifier
    domain proofs the `.py` lacks (the `.py` only spot-checks at sample points).
  - **219 / 222 / 223 / 225**: native `Det`/`Inverse`/`NullSpace`/`Series`/`Integrate`/
    `NSolve` on the shared physical setup vs deliverable closed forms (canonical
    dual-engine corroboration), several with substantive `.wl`-only extras (219's
    structural-family + `PositiveDefinite` proof, 222's numerator-reconstruction
    route, 225's `expectNonZero` negative control).

  Decisive test applied to every stage: *"if a subtle error were introduced into a
  load-bearing deliverable closed form, would the `.wl` catch it via its OWN
  computation?"* — yes for all 12.

## Re-confirmations of first-pass work — ALL HOLD

- **5 numerical typos (the −68/−51 family), each cross-engine-corroborated, NO stale
  survivor in live notes/cards/appendix** (old values survive only in process
  trackers as historical record):
  - 222: λ_W=0.2 upper-wall R_Q = `145.483858657863` (scripts emit it; notes match;
    no surviving `213.x`).
  - 223: λ_W=0.2 wall R_Q = `138.814136942081` / `137.502546600713` (no `206/205`).
  - 227: i=h rigidity det factor = `200+147π²` (DERIVED by `.wl`; no `251+215π²`).
  - 228: δ_1 = `196π²/(98π²−25)` AND reduced-det = `196(200+147π²)` (no `247…`/`251+215π²`).
  - 229: crossover-cubic leading coeff = `121ξ³` (no `189ξ³`).
- **2 renumbers landed** (221: 238→221, 237→220; 225: 240/241/242→223/224/225).
- **6 script-side de-taut/insufficient/hardcoded fixes still hold:** 220 (P_abs
  perfect-square, load-bearing), 221 F1 (survival round-trips) + F2 (deliverable #9),
  223 (compat-surface `sp.solve`), 224 (budgets→ceiling defining-relations),
  225 (`0==0`→`D4=−3D0u2²` + a negative control that genuinely FAILS),
  230 (onset round-trip→`sp.solve`; R_*≈1.229255438463336 / δ_*≈0.723111617875019 emitted+reflected).
- **Card +17 `\stagefield{Purpose}` scan (236–247): CLEAN** — no drift on 219–230 cards.

## Findings (5) — ALL the standing-deferred card-text-lag class

Stages 219, 220, 224, 226, 230 raised `paper_misalignment` for the card
`\stagefield{Verification}` field reading **"Mathematica audit: none yet"** despite
a passing `.wl`. This is a stale STATUS annotation, NOT a value/identity mismatch.
Several "clean"-verdict stages (222/225/227/229) noted the same lag informally —
i.e. it is present batch-wide; the agents merely applied an inconsistent threshold
on whether to formalize it (the documented under-call). **Per the standing user
decision → defer to PAPER_CLEANUP P4-51** (same as V.3's ~9 and VI.1's 16 cards).
Non-blocking; the stages still go `verified` at close-out with the lag deferred.

Per-stage audit verdicts (reports under `redteam/pass2/reports/stage_NNN.md`):
clean = 221(ckpt), 222, 223, 225, 227, 228, 229; findings(card-text-lag only) =
219, 220, 224, 226, 230. Directives written for the 5 (card-text-lag `Resolve`
blocks only; no Codex-applied script change).

## Stale-output (P4-52) + numbering residuals (deferred)

- **Only one stale committed-output banner batch-wide:** 221's MMA `.txt` prints
  "STAGE 204" (predates the re-author; the `.wl` SOURCE correctly prints "STAGE 221").
  The reliability-gate exec re-run + sed-refresh fixes it at close-out. All 11 other
  committed outputs carry canonical STAGE banners.
- 227's `.wl` output band carries stale Stage-242/243 self-labels → existing
  `NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (content-keyed, never offset-sweep).
- Residual project-wide notes-renumber drift (219/222/227 etc.) → deferred post-253
  stem-keyed cleanup, as standing.

## CLOSE-OUT PENDING — next session's first action

Audit + adjudication are done (this file + the 12 reports persist them). The
mechanical close-out has NOT run yet. Steps, in order:
1. **Reliability-gate exec re-run** all 12 (sympy + mma). Refresh committed
   `output/*.txt` via `sed '1,/^---$/d;/^# exit_code:/d'` (this is what refreshes
   221's stale "STAGE 204" mma banner). Avoid MANIFEST races: run scripts directly
   (`python3` / `math -script`) rather than parallel `$RT exec-*`; ≤2 concurrent
   `math -script`.
2. **Verify agents** (clean-context, one per stage) → write
   `redteam/pass2/verifications/stage_NNN.md` → `$RT set-status NNN verified` ×12.
   The 5 card-text-lag findings are paper_misalignment deferred to P4-51 → the
   stages still verify (scripts are clean).
3. **Sync the 6 prose trackers** under `notes/` (P5-20 = ZERO substantive paper/notes
   edits; record the card-text-lag deferral to P4-51).
4. **Commit** (reports + directives + verifications + refreshed outputs + this
   summary) + **doc-sync** the handoff/BATCHES/memories.
5. Then **HALT for the user gate** before VII.2 (sequential-audit-chunks rule).

## CLOSE-OUT RESULT — DONE 2026-06-09

The mechanical close-out ran and confirmed the audit. No surprises.

1. **Reliability-gate exec re-run (24 runs, all DIRECT — `python3`/`math -script`,
   ≤2 `math -script` concurrent, NOT parallel `$RT exec-*`):** all 12 sympy + 12 mma
   exit 0, FAIL=0, deterministic. SymPy 12/12 outputs byte-identical to committed.
   MMA 10/12 byte-identical; the two that changed were **label-only refreshes, no
   numeric change:**
   - **221** committed `.txt` banner `STAGE 204`→`STAGE 221` (the re-author fix the
     auditor flagged) — refreshed.
   - **227** committed `.txt` in-output cross-stage labels `Stage-242 K`→`Stage-225 K`
     and `M4 Stage-243`→`M4 Stage-226` — the `.wl` SOURCE was ALREADY canonical
     (lines 154/243/244 print 225/226); only the committed `.txt` was stale →
     refreshed. **CORRECTION to the "Stale-output" section above: 227's `.wl` source
     does NOT carry stale labels, so there is NO residual from 227 for
     `NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`.**
   Arbiter grep over the committed VII.1 outputs: CLEAN (no surviving `STAGE 204` /
   `Stage-242` / `Stage-243`).
2. **Verify agents (clean-context, one per stage) → all 12 `verified`** (exits 0/0,
   `material_change: false`). The 5 card-text-lag findings (219/220/224/226/230)
   classified as deferred paper_misalignment (P4-51) — non-blocking; stages verify.
   221 checkpoint independence re-confirmed by the agent (the §II opposite-info-flow
   derivation). `$RT set-status NNN verified` ×12 applied; pass-2 MANIFEST shows all
   12 `verified`. Pass-1 `MANIFEST.yaml` untouched.
3. **6 prose trackers synced** (P5-20; the card-text-lag deferral appended to P4-51).
   ZERO substantive paper/notes edits.
4. **Committed** (this file + 12 reports + 5 directives + 12 verifications + refreshed
   221/227 outputs + MANIFEST + 6 trackers) + doc-sync (handoff/BATCHES/memories).

## Cumulative

Pass-2: I.1 … VI.1 verified (218/253) + **VII.1 verified (219–230) = 230/253**.
Remaining pending: VII.2 (231–242), VIII.1 (243–253).
