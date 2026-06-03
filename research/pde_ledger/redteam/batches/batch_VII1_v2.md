---
batch: VII.1
range: 219-230
total_stages: 12
verified: 12
findings_count: 19
findings_resolved: 19
findings_blocked_legitimate: 0
material_change_count: 0
clean_stages: []
status_only: []
dirty_stages: [219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230]
checkpoints: [221]
consult: none
audit_date: 2026-06-02
verify_date: 2026-06-02
status: closed
---

# Red-team batch VII.1 — Mixed-bundle, resonance, branch-packet

## Summary

12-stage audit unit for VII.1 (`Part VII — Mixed-bundle / resonance /
branch-packet`), the **second batch of Part VII**, forward first-pass under the
v2 paper-grounded auditor **WITH the dual-engine rule** in force. **221 is the
only checkpoint** (higher verify bar); no status-only units. All 12 reached
`verified`; **`material_change: false` on all 12**; **0 stop-cold, 0 blocked, 0
needs_rework left open.** Every change is a strengthening / route change /
notes typo correction; no derived value, constant, identity target, or
PUBLISHED paper number on the SCRIPT side moved, so no `upstream_stale`
propagation.

**All 12 codex-invoke runs exited 0 on iteration 1 — no iter-2 reworks this
batch.**

**Cumulative: 218/253 → 230/253 stages red-team verified (90.9%)**; the entire
range 001–230 is now paper-aligned at v2 depth. Stages 231–253 remain `pending`.

### Headline — full dual-engine coverage with zero sanctioned mirrors

VII.1 ran with the dual-engine rule (a Mathematica `.wl` is **REQUIRED wherever
Mathematica CAN independently verify** — the test is "is it possible," not "is
it necessary") in force. **Every stage now has an INDEPENDENT second engine; 0
sanctioned mirrors were accepted in VII.1.**

- **11 stages had NO `.wl`** and got a **NEW independent-route `.wl`**: 219, 220,
  222, 223, 224, 225, 226, 227, 228, 229, 230.
- Checkpoint **221's pre-existing `.wl` was a line-by-line transliteration of the
  `.py`** and was **RE-AUTHORED to a genuinely independent route**: native
  `D[QPi/DeltaPi,portPi]` derivative, `Residue`, `ComplexExpand`, and an
  uncollapsed Breit–Wigner form.

Every new `.wl` was confirmed independent by a clean verify agent — native
primitives via a DIFFERENT decomposition than the SymPy `.py`. Representative
routes: **219** structural family extraction via `Collect`/`CoefficientList`;
**220** Laurent-support `CoefficientRules`; **226** `Orthogonalize` projector
(not QR); **228** `NSolve`+`Series`+implicit-function slopes; **229/230**
`Resolve[ForAll]` universal-quantifier proofs.

All 11 new `.wl` (219, 220, 222–230) plus the re-authored 221 are recorded in
the Independent-Mirror Set in `MATHEMATICA_MIRROR_POLICY.md`.

The labor split was strictly enforced: **Claude reviews** (audit + verify);
**Codex writes ALL script code** (designs and writes the new/re-authored `.wl`);
the directives stated only the requirement + acceptance criteria, never script
code.

## Checkpoint findings (the higher bar, no rubber-stamp)

### 221 — resonance line-shape tradeoff / dispersive no-free-lunch theorem + linear survival window

Cleared the higher checkpoint bar via three substantive closures plus a notes
renumber:

- **`.wl` RE-AUTHORED** — the existing `.wl` was a line-by-line transliteration
  of the `.py`; rebuilt to a genuinely independent route (native
  `D[QPi/DeltaPi,portPi]` derivative, `Residue`, `ComplexExpand`, an uncollapsed
  Breit–Wigner form).
- **F1 `tautological_check`** — the survival round-trips were de-tautologized.
- **F2 `tautological_check` / `insufficient_verification`** — deliverable #9 (the
  linear survival window) had been print-only/tautological; it is now genuinely
  covered in BOTH engines.
- **F4 `paper_misalignment`** — notes renumber (238→221, 237→220).
- **Sanctioned Codex deviation (F3):** Codex used the native
  `D[QPi/DeltaPi,portPi]` instead of the directive's leading-minus form,
  reconciling the Stage-220 identity `∂_Π D_Π = −N` — verified correct (a route
  choice; the identity holds).

**Checkpoint bar MET.**

## The dual-engine-anchored load-bearing constants

The VII.1 load-bearing constants are now anchored in BOTH engines — each
corrected value is independently recomputed by the new `.wl` and corroborates
the verified SymPy script:

- the corrected R_Q figures `145.483858657863` (222, λ_W=0.2),
  `138.814136942081` / `137.502546600713` (223, λ_W=0.2);
- the `i=h` rigidity determinant factor `200+147π²` (227);
- the δ_1 coefficient `196π²/(98π²−25)` and reduced-det `196(200+147π²)` (228);
- the crossover-cubic leading coeff `121ξ³` (229);
- the 230 thresholds `R_*≈1.229255438463336` / `δ_*≈0.723111617875019`.

The seven notes-only typos these corrected were a systematic **+68/+51 additive
family** (see Paper / notes edits below). **The SCRIPTS were always correct**;
the typos lived only in the notes prose. Provenance is logged at the 221
checkpoint in `CHECKPOINT_CONSTANT_PROVENANCE.md`.

## Per-stage findings tally

| Stage | Status | Findings | Notes |
|-------|--------|----------|-------|
| 219 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only) via structural family extraction (`Collect`/`CoefficientList`) |
| 220 | dirty | 1 + `.wl` | New independent-route `.wl` via Laurent-support `CoefficientRules`. F2 `insufficient_verification` → symbolic P_abs perfect-square assert. Informational notes-label drift renumbered by Codex |
| 221 | dirty (ckpt) | 3, `.wl` re-authored | **Checkpoint.** Transliteration `.wl` RE-AUTHORED (native `D[QPi/DeltaPi,portPi]`/`Residue`/`ComplexExpand`/uncollapsed Breit–Wigner); F1 tautological survival round-trips de-tautologized; F2 deliverable #9 (linear survival window) print-only/tautological → genuine dual-engine coverage; F4 notes renumber (238→221, 237→220). One sanctioned Codex deviation (F3): native `D[QPi/DeltaPi,portPi]` vs the directive's leading-minus form, reconciling Stage-220 `∂_Π D_Π = −N` — verified correct. Checkpoint bar MET |
| 222 | dirty | 1 + `.wl` | New independent-route `.wl`. **F2 paper_misalignment (notes-only):** λ_W=0.2 upper-wall R_Q `213.483858657863`→`145.483858657863`; cross-engine corroborated by the new `.wl` |
| 223 | dirty | 2 + `.wl` | New independent-route `.wl`. F2 `tautological_check` (circular compat-surface) → `sp.solve(Eq(K_norm,K_pole),P0_target)`. **paper_misalignment (notes-only):** λ_W=0.2 wall R_Q `206.814136942081`/`205.502546600713`→`138.814136942081`/`137.502546600713`; cross-engine corroborated. Notes-label renumber (240/241/242→223/224/225 + a filename token) as a formal finding |
| 224 | dirty | 1 + `.wl` | New independent-route `.wl`. F2 `hardcoded_result` (budgets-vs-themselves) → defining-relation checks tied to the ceiling. Informational notes-label drift renumbered by Codex |
| 225 | dirty | 1 + `.wl` | New independent-route `.wl`. F2 `tautological_check` (`0==0` one-pole) → tied to the one-pole constraint `D4=−3·D0·u2²` (+ an `expectNonZero` negative control). Notes-label renumber (formal finding) |
| 226 | dirty | + `.wl` | New independent-route `.wl` via an `Orthogonalize` projector (not QR). Informational notes-label drift renumbered by Codex |
| 227 | dirty | 1 + `.wl` | New independent-route `.wl`. **paper_misalignment (notes-only):** i=h rigidity det factor `251+215π²`→`200+147π²`; cross-engine corroborated by the new `.wl` |
| 228 | dirty | 1 + `.wl` | New independent-route `.wl` via `NSolve`+`Series`+implicit-function slopes. **paper_misalignment (notes-only):** δ_1 coeff `247π²/(98π²−25)`→`196π²/(98π²−25)` AND reduced-det `247(251+215π²)`→`196(200+147π²)`; cross-engine corroborated |
| 229 | dirty | 1 + `.wl` | New independent-route `.wl` via a `Resolve[ForAll]` universal-quantifier proof. **paper_misalignment (notes-only):** crossover-cubic leading coeff `189ξ³`→`121ξ³`; cross-engine corroborated |
| 230 | dirty | 1 + `.wl` | New independent-route `.wl` via a `Resolve[ForAll]` universal-quantifier proof. F2 `tautological_check` (onset round-trip) → `sp.solve(Eq(onset,R_star),delta)`; thresholds `R_*≈1.229255438463336`/`δ_*≈0.723111617875019` corroborated. Codex log captured EMPTY (logging anomaly, not a stall — exit 0, all `## Applied` blocks present, independently re-confirmed). Informational notes-label drift renumbered |

**Totals:** 19 findings closed (the missing-`.wl` dual-engine gap on 219/220/222–230
+ 221's transliteration; the script-side strengthenings on 220/223/224/225/230;
the 7 notes-only paper_misalignment typos on 222/223/227/228/229; and the formal
notes-label renumbers on 221/225), 0 blocked, 0 status-only. **11 new
independent-route `.wl` (219, 220, 222–230)** + 1 re-authored checkpoint `.wl`
(221).

## Mathematica mirror policy — full dual-engine, zero sanctioned mirrors

VII.1 ran with the dual-engine rule in force. **0 sanctioned mirrors were
accepted.** All 11 new `.wl` (219, 220, 222–230) are GENUINE INDEPENDENT routes;
the checkpoint 221's pre-existing `.wl` (caught as a transliteration) was
**RE-AUTHORED** to an independent route (native `D[QPi/DeltaPi,portPi]`
derivative, `Residue`, `ComplexExpand`, an uncollapsed Breit–Wigner form). All
12 are added to the Independent-Mirror Set in `MATHEMATICA_MIRROR_POLICY.md`.

## Paper / notes edits (Codex-applied, Claude-reviewed)

This batch APPLIED notes edits (per the file-ownership contract: Codex owns
`paper/*.tex` + `notes/stages/*.md` edits, Claude reviews). **ALL paper_misalignment
items were NOTES-ONLY; the PUBLISHED paper cards/appendices were UNAFFECTED —
they carry abstract forms.** Each correction was cross-engine corroborated by the
new `.wl` independently computing the corrected value. All were
orchestrator-reviewed: correct, isolated, no collateral.

**Numerical typos (a systematic +68/+51 additive family):**

- **222:** notes λ_W=0.2 upper-wall R_Q `213.483858657863` → `145.483858657863`.
- **223:** notes λ_W=0.2 wall R_Q `206.814136942081` / `205.502546600713` →
  `138.814136942081` / `137.502546600713`.
- **227:** notes i=h rigidity det factor `251+215π²` → `200+147π²`.
- **228:** notes δ_1 coeff `247π²/(98π²−25)` → `196π²/(98π²−25)` AND reduced-det
  `247(251+215π²)` → `196(200+147π²)`.
- **229:** notes crossover-cubic leading coeff `189ξ³` → `121ξ³`.

**Notes-label renumbers (labels → canonical):**

- **221:** 238→221, 237→220 (formal finding F4).
- **225:** 240/241/242→223/224/225 + a filename token (formal finding).
- **220/224/226/230:** informational label drift, also renumbered by Codex.

## Out-of-scope residuals — LOGGED, not fixed this batch

Logged to `PAPER_CLEANUP_TRACKER.md` P4-53 so they are not lost. Both are
notes/banner cosmetics with **ZERO effect on verification**:

- **VII.1 notes multi-epoch stage-renumber drift** — the audit-flagged labels
  were renumbered + verified per-stage, but residual stale stage-number
  references remain in the notes of 219, 221 (partial), 222, 223 (its TITLE
  still reads "Stage 240"), 227, 228, 229. These come from MULTIPLE historical
  renumbering epochs with INCONSISTENT offsets within and across files (e.g.
  219's notes self-number "253" while 221/222/225 self-number "238/239/242"), so
  a uniform blind sweep is UNSAFE — needs a careful dedicated pass with
  per-reference content/stem matching (cf. the IV.4/IV.5 computed-offset
  reattribution approach).
- **221's `.wl` top banner** prints a cosmetic stale "STAGE 204" Print string (no
  finding scope, affects no check).

## Iteration-2 reworks

None. All 12 codex-invoke runs exited 0 on iteration 1.

## Orchestrator notes

1. **230's codex log captured EMPTY** — a logging anomaly, NOT a stall: exit 0
   was recorded, all `## Applied` blocks were present, all artifacts were built,
   and the stage was independently re-confirmed + verified.
2. **The VI.1 `.wl` naming bug did NOT recur** — all 12 directive `.wl` targets
   carried the required `_mathematica_audit.wl` suffix; an explicit pre-invoke
   guard worked.

## Consult

None. No math directions required a Claude+Codex resolution; nothing escalated to
the user. (The 7 notes-only paper_misalignment numerical typos were USER-RESOLVED
2026-06-02 — corrected to the already-correct SymPy scripts.)

## Verification

All 12 verification files under `redteam/verifications/stage_219.md` …
`stage_230.md`. Final verdicts:
- `verified` (12): 219–230.
- `needs_rework` → reworked → re-`verified`: none.
- `blocked_unfixable` (0).

Material change: **0** (`material_change: false` on all 12 — new/re-authored
second engine, de-tautologized checks, and notes-only numerical-typo corrections
that align prose to the already-correct script; no SCRIPT-side derived value,
constant, identity target, or PUBLISHED paper number moved).

## Cumulative

Range 001-230 paper-aligned at v2 depth. **230/253 stages red-team verified
(90.9%)** (was 218 after VI.1; VII.1 adds 12 across 219–230). **Second batch of
Part VII**; zero stop-cold, zero material change, full dual-engine coverage with
zero sanctioned mirrors, and the one checkpoint (221) cleared the higher bar.

Next batch (sequential-audit-chunks rule, awaits explicit user authorization):
**VII.2 onward = stages 231–253** (all currently `pending`). The planned full
end-to-end **second pass** remains a later cross-check, only after the first pass
reaches stage 253.
