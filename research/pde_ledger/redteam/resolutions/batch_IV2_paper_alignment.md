---
batch_id: IV.2
range: 103-114
created_at: 2026-05-27
total_paper_misalignment_findings: 5
clusters: 3
status: awaiting_user_resolution
---

# Batch IV.2 v2 — Paper-Alignment User-Gate Resolutions

Five `paper_misalignment` findings across 4 stages (105, 106, 108, 109, 112) plus a 12-stage banner survey. The 5 findings consolidate into **3 clusters** matching prior-batch patterns; each cluster has one orchestrator-recommended resolution and 1-2 alternatives.

The remaining 11 findings (105 F1, 106 F2/F3/F4, 107 F1, 108 F2/F3, 109 F1/F2, 111 F1, 112 F2) are script-side fixes that do not require user input — they will be applied as soon as the cluster directions below are chosen.

---

## Cluster A — Paper-card `Checks` coverage (106 F1, 109 F3)

**Pattern:** Paper card lists `Checks` items that the stage's scripts do not directly verify. The items are exercised at neighbour stages (upstream for 106; downstream for 109). Mirrors IV.1 Cluster A (091/095/097/099) — docstring carry-forward annotation naming the verifying neighbour stages.

### 106 F1 (script_missing_paper_claim)

- Paper card (`paper/stages/stage_106.tex:21-25`) lists three Checks:
  1. Factorization separability `m̂₀² · χ_Q · N_Q` — partially covered (factored product is displayed and N_Q solved at script level).
  2. Higher odd terms beyond 2.5PN — **upstream at stage 102** (higher_odd_irrelevance, audited IV.1).
  3. l=2 DtN fingerprint z-expansion — **upstream at stage 104** (outgoing_dtn_fingerprint, this batch, verdict clean).
- Notes file already says "From Stage 88 [renumbered to upstream] the canonical compact passive/outgoing grouped-P₂ DtN model gives χ_Q=1" — confirming the carry-forward intent.

### 109 F3 (script_missing_paper_claim)

- Paper card (`paper/stages/stage_109.tex:21-25`) Checks:
  - "Check pure scale and pure argument deformations separately" — pure scale is covered; pure argument is part of the linearization that the script does perform.
  - "Robin and standalone mixed-pole limits before imposing compensation" — **downstream at stages 110 (robin_outlet_model) and 111 (mixed_sidechannel_pole)**.
  - "Compensated branch preserves even coefficients as well as odd normalization" — **downstream at stage 112 (hybrid_robin_mixed_compensation)**.
- Notes file is silent on these neighbour-stage checks (the notes' deliverables are the boxed expansion + closed-form `a₅ = -5b/9 - a₀/27`).

### Resolution options

- **(a) Script-side docstring carry-forward (Recommended).** For both stages: prepend a docstring block naming the neighbour stages that verify the deferred Checks. Mirrors IV.1 Cluster A. No new assertions, no paper edits.
  - 106 docstring add: "Checks (ii) and (iii) of the paper card are verified upstream at stage 102 (higher-odd irrelevance) and stage 104 (outgoing l=2 DtN fingerprint, from which χ_Q=1 is derived); this stage uses χ_Q=1 as a carry-in."
  - 109 docstring add: "Robin and mixed-pole Check items of the paper card are exercised at downstream stages 110 (Robin outlet model) and 111 (standalone mixed-pole); even-coefficient preservation under the full hybrid ansatz is verified at stage 112."

- **(b) Substantive script extension.** For 106: add z-expansion of Λ₂^out and Ŷ_2^out duplicating stage 104's work. For 109: add Robin/mixed-pole parameters to the linearized ansatz, duplicating stages 110/111. Significantly expands script scope; departs from the no-duplication principle the carry-forward chain assumes.

- **(c) Paper-side trim.** Edit `paper/stages/stage_106.tex` and `stage_109.tex` to remove (ii)/(iii) from the Checks list and reference the carry-in stages directly. Affects paper text; deferred (paper changes go to PAPER_CLEANUP_TRACKER for a batch-wide pass).

**Verification grep (for (a)):**
```
grep -E "stage 102|stage 104" scripts/moving_throat_pde_stage106*.py
grep -E "stage 110|stage 111|stage 112" scripts/moving_throat_pde_stage109*.py
```

---

## Cluster B — β-parameterized preservation submanifold (108 F1)

**Pattern:** Notes file boxes a 1-parameter submanifold; scripts verify only the β=1 reduction. Substantive coverage gap (not a carry-forward case — no neighbour stage exercises the general locus). Mirrors IV.1 Cluster B (100 closure derivation) — adds real new asserts to both engines.

### 108 F1 (script_missing_paper_claim)

- Notes (`notes/stages/moving_throat_pde_stage108_robustness_classes.md:84-92`) state the exact preservation condition for the combined scale+argument+additive deformation class:
  ```
  Σ₅ = S(1 - β⁵)/9 - Σ₀/27.
  ```
- SymPy (`scripts/...stage108...sympy_audit.py:50-70`) only verifies the β=1 reduction `Σ₅ = -Σ₀/27` (Class C, additive-only); same for Mathematica.

### Resolution options

- **(a) Substantive script extension (Recommended).** Add to both engines: combine scale+argument+additive ansatz as
  `Λ_add(z) = S · Λ_out(βz) + Σ₀ + Σ₂z² + Σ₄z⁴ + i Σ₅ z⁵`,
  re-solve `(Σ₂, Σ₄)` against the even-moment match (solutions now depend on β), then assert that χ_Q=1 iff `Σ₅ = S(1-β⁵)/9 - Σ₀/27`. The existing β=1 reductions remain as sanity checks. This is the substantive addition; the existing tests are subsumed.

- **(b) Trim notes/paper to match scripts.** Edit `notes/stages/...stage108...md` and the part-04 appendix row to downgrade the preservation submanifold to "β=1 case only," removing the boxed general formula. Paper-side; deferred.

- **(c) Defer general locus to a later audit unit.** Add a `Scope` note in the stage card stating the general locus is not yet script-verified at this unit, and create a follow-up TODO. Documented gap.

**Verification grep (for (a)):**
```
grep -E "beta\\*\\*5|Lambda.*beta\\*z|1 - beta\\*\\*5" scripts/moving_throat_pde_stage108*.py
grep -E "beta\\^5|lambdaOut\\[beta" mathematica/moving_throat_pde_stage108*.wl
```

---

## Cluster C — Stage-label / banner sweep (105 F2, 112 F1, + 8 other stages, paper card titles)

**Pattern:** Scripts banner stale "Stage 87/88/89/90/91/92/93/94/95/97" labels left over from a prior global renumbering. The IV.1 audit found the same pattern and resolved via script-side sweep + paper-side deferred. Two stages have the issue *filed as a finding* (105 F2, 112 F1, 108 F3); the other 8 stages with stale banners are not separately filed (auditor noted but treated as cosmetic), but every IV.2 script is affected. Mirrors IV.1 Cluster C — 24-site sweep with paper card titles deferred to PAPER_CLEANUP_TRACKER.

### Sites mapped (24 total across 10 stages)

| Stage | .py sites | .wl sites |
|---|---|---|
| 104 | banner L18 (`STAGE 87`) | banner L26 (`STAGE 087`) |
| 105 | docstring L3 (`Stage 88`), banner L28 (`STAGE 88`) | banner L26 (`STAGE 088`) |
| 106 | docstring L3 (`Stage 89`), banner L25 (`STAGE 89`) | banner L27 (`STAGE 089`) |
| 107 | banner L25 (`STAGE 90`) | banner L26 (`STAGE 090`) |
| 108 | banner L25 (`STAGE 91`) | banner L26 (`STAGE 091`) |
| 109 | banner L23 (`STAGE 92`) | banner L26 (`STAGE 092`) |
| 110 | print L31 (`'stage93: PASS'`) | banner L26 (`STAGE 093`) |
| 111 | print L36 (`'stage94: PASS'`) | banner L26 (`STAGE 094`) |
| 112 | docstring L3 (`Stage 95`), print L54 (`'stage95: PASS'`) | banner L26 (`STAGE 095`) |
| 114 | banner L18 (`STAGE 97`) | banner L26 (`STAGE 097`) |

### Paper card display titles (deferred)

`paper/stages/stage_NNN.tex` headings use a separate "display Stage" numbering from the internal unit number — e.g., stage 105 card says `\section[Stage~122]{...}`, stage 112 card says `\section[Stage~129]{...}`. These display numbers also drift relative to file labels and were noted in IV.1 Cluster C as deferred to PAPER_CLEANUP_TRACKER.

### Resolution options

- **(a) Script-side banner sweep; paper card titles deferred (Recommended).** Update all 24 script-side sites to read `STAGE 104` through `STAGE 114` (matching the file-path / label-id internal numbering, as IV.1 chose). Paper card display titles `\section[Stage~NNN]{...}` left unchanged and logged in `PAPER_CLEANUP_TRACKER.md` for a later batch-wide paper pass. Mirrors IV.1 exactly.

- **(b) Adopt the paper display numbers (Stage 108-119) script-side.** Re-banner scripts to `STAGE 108` through `STAGE 119` using the public-display numbering. Diverges from IV.1's choice (which used internal unit numbers); breaks the rule that script transcripts grep against the internal stage_NNN.tex labels.

- **(c) Keep stale provenance labels.** Add a `# formerly Stage NN` comment alongside each banner but leave banners as-is. Preserves historical lineage but leaves transcripts misleading.

**Verification grep (for (a)):** all of:
```
grep -lE "STAGE 0?(87|88|89|9[0-7])" scripts/moving_throat_pde_stage10[3-9]*.py scripts/moving_throat_pde_stage11[01234]*.py mathematica/moving_throat_pde_stage10[3-9]*.wl mathematica/moving_throat_pde_stage11[01234]*.wl
# must be empty after sweep
grep -lE "STAGE (10[4-9]|11[01234])" scripts/moving_throat_pde_stage10[3-9]*.py mathematica/moving_throat_pde_stage10[3-9]*.wl
# must enumerate all 10 stages
```
