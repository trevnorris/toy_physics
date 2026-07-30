# The census pilot — what it measured, and what it found about the schema

**Status:** PILOT RESULT. Three stages, both engines, spec→pilot per the recorded sequence.
⛔ This is not the census. It is the check that the schema survives contact before a 43-stage fan-out.

Artifacts: `REPORTED_RESULTS.md` (universe prerequisite) · `pilot/stage016_rows.md` ·
`pilot/stage023_rows.md` · `pilot/stage043_reconciliation.md`.

---

## 1. ⛔ READ THIS BEFORE QUOTING ANY NUMBER BELOW

**`DERIVED = 0` across both row-producing stages.** That is the schema's headline quantity, and it is
the honest output of the rules as written.

⛔ **It does NOT mean "the model derives nothing."** The measured cause is narrower and is stated in §3:
**neither engine cites a source locus anywhere.** Both name their source *stages* and neither gives a
line. Under §3.3 a `C-PEER` leaf requires a cited source locus, so no closure ever reaches a
`C-FIELDEQ` or `C-EXTERNAL` leaf, and every closure bottoms out in free symbols.

⇒ The correct reading is: **no occurrence in these two stages carries evidence that physics entered
from outside the artifact.** Whether the physics *did* enter is a separate question this pass does not
answer — and largely cannot, until the artifacts carry citations.

⚠ Both readings must travel together. Dropping the first makes the ledger look sound; dropping the
second makes it look empty. Neither is measured.

---

## 2. The numbers

| | stage016 | stage023 | stage043 |
|---|---|---|---|
| occurrences | 38 | 136 (66 `.py` / 70 `.wl`) | **0** |
| proposed QIDs | 16 | 65 | — |
| TIER 1 (occ) | `[0, 30]` | `[15, 111]` | — |
| TIER 1 (QID) | `[0, 13]` | `[8, 57]` | — |
| TIER 2 · TIER 3 · `DERIVED` | 0 · 0 · 0 | 0 · 0 · 0 | — |
| blocked rows | **0** | **0** | — |
| conflict set | 0 QIDs | **6 QIDs** | — |

**Combined: 174 occurrences, ≤81 QIDs, TIER 1 = `[15, 141]` occ / `[8, 70]` QID.**
⚠ `81` and `70` are **upper bounds** — cross-stage identity is unadjudicated (§8.3), and every merge
lowers them. ⛔ Do not quote them as counts.

⚠ **The tier-1 spread is dominated by `unadjudicated`, not by evidence** (30 of 38 · 96 of 136). Its
width is set by schema gap G1 below, not by the corpus. ⛔ The lower bound is not a finding yet.

**stage043 contributes zero occurrences.** A manifest row asserts a *category* for a knob, not a
*value* of it, so §7.2 admits none of the 152. Two independent rules reach this. That is what the
stage is; it is not a defect.

---

## 3. ⭐⭐ The mechanism — the most actionable thing the pilot found

**The ledger's scripts name their sources and never cite them.** stage023: `PY:1049-1054`,
`WL:765-770`, comments at `PY:336`, `PY:344`, `WL:196` — all name a stage, none give a locus.

⇒ **This single property drives TIER 3 and `DERIVED` to zero.** 66 of stage023's 96 unadjudicated rows
are `A-REDUCED` with a `C-UNRESOLVED` closure — the reduction is there, the provenance is not.

⭐ **It is cheap to move.** Adding a source locus to each cross-stage import is mechanical, requires no
physics, and would resolve a large fraction of these rows *in either direction*. ⛔ It is not a fix that
flatters: a citation that fails to check out under §9.1 is a finding, not a pass.

---

## 4. Schema gaps, ranked

⭐ **G1 and G2 were each found independently by two legs that did not share context.** That is the
strongest signal available here.

- **G1 — no precedence between an unresolved closure and evidenced debt.** §3.3/§5.7 send a
  `C-UNRESOLVED` closure to `unadjudicated`; §5's projection and §9.0's "axis-C failure does not touch
  axis A" leave an evidenced `A-REDUCIBLE-UNDERIVED` at `tier1-debt`. A row matching both has no rule.
  *(stage016 C1 · stage023 GAP-2, independently.)*
  ⇒ **This sets the tier-1 lower bound.** stage016 alone: `[0,30]` under one reading, `[16,30]` under
  the other. **Fix first.**
- **G2 — the universe rule excludes postulates, which is what tier 1 is FOR.** §7.1 admits only values
  that are numeric or closed-form symbolic. Structural postulates are label-valued — a BC class, a
  time-arrow, an existence claim — so the expression-walk never reaches them, while §3.1's
  `A-IRREDUCIBLE-POSTULATE` exists for exactly them.
  *(stage043 F2, on the 11 `discrete-structural-postulate` IDs · stage023 §F.4, on `pathA29_premise`,
  "the keystone premise" that no expression consumes — independently.)*
  ⇒ **Both stages report `tier1-postulate = 0`, and that zero is an artifact of the rule, not a
  measurement.** ⛔ It must not be reported as a finding until G2 is fixed.
- **G3 — the universe boundary is ambiguous at a verdict token.** `REPORTED_RESULTS` cites `PY:1068`,
  a verdict token whose value data-depends on all seven ladder guards; a strict reading pulls the whole
  dimensional gate in (**+122 occurrences**, roughly doubling the stage). The schema decides neither
  way. *(stage023 GAP-1.)*
- **G4 — §8.4's reconciliation binary is not exhaustive.** A quantity can be **in** the universe yet
  have **zero occurrences**, because no artifact assigns it a value: `REG:b:K0c` carries a dimension and
  a class but no value anywhere. Neither route fits. **Threatens up to 46 of 152.** *(stage043 F1.)*
- **G5 — a named route that needs a framework extension.** §3.1 reads it as debt; §3.1.1's test reads
  it as structural. Moves 7 occurrences between `tier1-debt` and `tier1-structural` on stage023 alone.
- **G6 — missing enum values, six of them**: no axis-B value for a free-symbol declaration; no binding
  site for a Wolfram global; no `should_be_basis` for "should be convention"; no closure hop for a time
  convention that fixes a claimed sign; no axis-A value for an independent variable; §10.2.2 silent on
  QIDs mixing `unadjudicated` with `unclassified-nonfed`.
- **G7 — `C-PEER`'s citation may live in the stage note, not the script**, and the schema does not say
  whose citation counts. Consequential given §3: the notes *do* carry loci where the scripts do not.
- **G8 — §7.3's exclusion vocabulary is factually wrong** for the parallel-track category: it offers
  "retired or excluded", while the register says **PENDING** and the note says **PAUSED**.

Further unresolved: `C-SELF`/`C-MATH` precedence · whether `rank()`/`NullSpace` is a §3.3.1 operator ·
`SELF-REFERENTIAL` missing a round trip through a *different* occurrence · two competing substrate
claims on one occurrence (§10.3 covers only cross-occurrence).

---

## 5. What the pilot validated

- ⭐ **Zero blocked rows across 174 occurrences.** Every occurrence reached a classification, including
  the ones the schema could not decide — those were recorded as contested with both readings, never
  forced. The schema is usable.
- ⭐⭐ **The axis split is demonstrated, not asserted.** **20 of the 21** records behind stage016's
  `21/0` dimensional verdict fall **outside** this universe. The provenance and dimensional-correctness
  populations are near-disjoint *by construction*. ⇒ The standing instruction not to treat `24/0/10`,
  `27/0/7` or `21/0` as an oracle was load-bearing, not caution.
- **The conflict set found real structure**: 6 stage023 QIDs are each computed at one binding site and
  typed as a target at another — exactly the cross-engine asymmetry the census exists to surface.
- **The manifest re-measured clean**: 152 literals, 152 distinct, all ten category counts reproducing
  the script, and `[40,49]` rederived from the disjoint peers (22+18 low, +9 high).
- **Coverage was reported honestly everywhere but one place** — stage023's test-scaffolding/display
  class was excluded by name without an individual count. Recorded as the single uncounted exclusion.

---

## 6. Verdict on fan-out

⛔ **Not ready.** G1 and G2 must be fixed first: G1 because the tier-1 lower bound is currently
undefined, and G2 because it systematically excludes the category tier 1 is expected to be mostly made
of. Fanning out across 43 stages before both are settled multiplies a known defect by 43 — the failure
the spec→pilot→fan-out sequence exists to prevent.

⭐ **Recommended before fan-out, in order:** fix G1 · fix G2 · settle G3's universe boundary · extend
§8.4 for G4 · then re-run this pilot's two row-producing stages and diff against these numbers.

⚠ **A note on what a re-run should NOT be judged by.** These distributions are not a target. If a
schema fix moves them a lot, that is the fix working. ⛔ Reproducing `[15,141]` is not evidence of
anything.
