I'll review this as a document check of the S11c-c2 step record against its sources. First I'll read the governing N6 disposition, tally, physics review, and related records, then form a view of what c2 actually established before opening the step record.I'll start with the N6 disposition, tally, resolution, and physics-review adjudication so I can form a source-grounded view before opening the step record.Next I’ll pull SHARED_PHYSICS, the predecessor records, and the `I_{M→E}` / engine artifacts so the claim set is grounded before I open the step record.I have the N6 and adjudication picture. Next I’ll read SHARED_PHYSICS §5c, the diagnostic for `I_{M→E}`, the predecessor carry-ins, and the SHA/count evidence.I have the source picture. A few SHA and construction details still need checking, then I’ll read the step record against that view.Source view is in place. I’ll read the step record now and check claims against those sources, quoting both sides on anything that would change what it may claim.I have one more construction check on `I_{M→E}` and the SHA parenthetical, then I’ll write the verdict.Source-grounded view (formed before the record): S11c-c2’s **fold wiring + A/C/D1–D6 are per-engine SOUND**; the `8f3a017f` “0 defects” subject is **overstated** and was corrected. **F/G interpretations rest on withdrawn instruments** and numeric re-grounding is **paused**; increment **values** are unaffected. **§5c was mis-specified** (cross-anchoring) and **corrected** (`30d4b72d`) to Eulerian-vs-material-**coordinate** at **fixed** anchoring. Per-engine N6 is **Reading B**: `R_N6 = 18/288` nonzero in **3 of 4** cases **and** `R_cov` no-nonzero — both preserved. Cross-engine work on this box is **N6 only** (no WL self-energy engine, no self-energy comparator). Path B: covariance-channel matched zeros are **`(0)−(0)` vanishing confirmation, not operand agree**; carrier **40** / source **76** / Φ **18** are **UNADJUDICATED DEBT**; leftover **SHAPE un-inspected**. Carries: 3 N6 caveats, face-force + #90 (kinetic **cancels**), 6 §3d, c1 ENERGY, F-wording, F/G paused, census 6-row.

Checked and **not** findings: “0 defects” is disclaimed; no “weak N6” / “c2 has everything it needs” / “just thickness”; `(0)−(0)` vs 40/76/18 split matches the tally; kinetic-vs-face-force distinction is correct; 3 caveats, 6 §3d, ENERGY, F-wording (§5e still “must vanish”), F/G paused, census crosswalk, S11c-d material-debt, no per-substep card; cited git SHAs `16849fc6`, `30d4b72d`, `8f3a017f`, `aa76105a`, `d21c8ff5`, `e11f2f82`, `ae73b884`, `48a0b4e7`, `a094b284`, `0bca95f3`, `2d12f287`, `cfb2494c`, `28f87dec` exist and match their subjects. Artifact listing confirms fold is SymPy-only; only N6 has WL + comparator.

---

### Finding 1 — `I_{M→E}` labeled “MATERIAL-**anchoring** increment” revives the §5c category error

**Record** (`steps/S11c_c2_self_energy_fold.md:132–133`):

> `I_{M→E}` is the **MATERIAL-anchoring increment**, ⛔ NOT a "mapped-to-Eulerian operand." `I_{M→E}^{α,ρ} = extract(close(SLAB_M) − SLAB_M)`

**Sources:**

- `directives/S11c_c2_SHARED_PHYSICS.md:303–315`: N6 is the two **COORDINATE** constructions **at a FIXED anchoring**; “the routes are the **representation** axis, ⛔ never the **anchoring** axis”; `I_{M→E}^{α,ρ} = extract(close(SLAB_M)−SLAB_M)` with **α held fixed**.
- The record’s own `:83–86`: the fold engine’s residual “compared the two **ANCHORINGS** (distinct physics) — the **WRONG object**”; real N6 is “Eulerian-vs-material-**coordinate** within a **FIXED anchoring**.”
- `_measurements/S11c_c2_N6_RESOLVED.md:54–55`: resolve the misleading **“mapped-operand”** label (it is not the fully-mapped-to-common-variables operator) — ⛔ not a rename onto the anchoring axis.
- `scripts/S11c_c2_N6_diagnostic_sympy.py:850–853` and `_measurements/S11c_c2_N6_build_clearance.md:17`: `M` is `build_increment(..., m_coeff, ms, ...)` from the **MATERIAL** face factory at the **same** `(α,ρ)`; **no `T`** on the increment.

The no-`T`/no-pullback / “not mapped-to-Eulerian” / `R_cov`-separate half of the paragraph is **correct**. The name **“MATERIAL-anchoring”** is not: in this ledger **anchoring** = `{LAB_HELD, MATERIAL_ADVECTED}`. A later reader can quote line 132 as if N6 compared two anchorings — the exact mis-specification `30d4b72d` corrected.

**Must change:** replace “MATERIAL-anchoring increment” with **material-coordinate / native material-route increment at fixed anchoring `α`** (or equivalent that cannot be read as the anchoring axis). Keep the no-`T` and `R_cov` sentences.

---

### Finding 2 — `e5cea55b` is attached to the WL **`.out`**; it is the **engine file** sha256

**Record** (`:62–63`):

> the compared `.out` is `ae73b884` (regenerated, sha `e5cea55b`, annex/GIN — `datalad get` after a fresh checkout)

**Sources:**

- `sha256sum mathematica/S11c_c2_N6_mathematica_audit.wl` → `e5cea55b4c06…` (the **`.wl` engine**).
- `sha256sum mathematica/out/S11c_c2_N6_mathematica_audit.out` → `6531a9f11c81…` (the **`.out`**).
- `git cat-file -t e5cea55b` → not a git object.
- `ae73b884` subject: regenerated from certified engine `e11f2f82` **(sha `e5cea55b`)**.
- Disposition `_measurements/S11c_c2_N6_reconcile_disposition.md:32`: `.out` `ae73b884` regenerated from engine `e11f2f82`, sha `e5cea55b`.

`ae73b884` as the compared `.out` commit is right; annex/GIN is right for the `.out`. Putting `sha e5cea55b` in the `.out` parenthetical next to `datalad get` attributes the **engine** content hash to the **transcript**. A verify-the-`.out` reader will not recover `e5cea55b`.

**Must change:** attach `e5cea55b` to `mathematica/S11c_c2_N6_mathematica_audit.wl` / engine `e11f2f82`; keep `ae73b884` as the `.out` commit (and annex/GIN on the `.out` only).

---

**NOT-SOUND** — must-fix: (1) `I_{M→E}` “MATERIAL-anchoring” name; (2) `e5cea55b` attributed to the `.out`.
