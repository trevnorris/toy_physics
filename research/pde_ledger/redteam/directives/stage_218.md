---
unit_id: 218
batch: VI.1
created_at: 2026-06-02T03:10:39Z
findings_count: 5
stop_cold: null
applied: true
applied_at: 2026-06-02T18:30:40Z
findings_applied: 5
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 218

Apply EVERY finding below in order (F1, F2, F3, F4, F5). After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F4 was a `paper_misalignment` (notes/appendix per-envelope value) that the USER HAS NOW RESOLVED with direction (a): the script's `162` (= 3⁴·2) is correct; the notes' `230` and the appendix's `179` are typos. You are AUTHORIZED to apply the F4 NOTES edit (`230 → 162`) per F4's `## RESOLVED` block. The appendix line `stage_appendix_part06.tex` (`3^4·2=179`→162) is OWNED and applied by the stage-217 directive, which runs BEFORE this one — do NOT edit the appendix here; just confirm it already reads `162`. Do NOT change the SymPy script's `162` (it is correct).

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

These are CLAIM-MANIFEST directives: they state the REQUIREMENT and ACCEPTANCE CRITERIA only. You (Codex) design and write the actual check code. Do NOT introduce new features, refactors, or stylistic changes beyond what each finding requires. Edit only the files named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

---

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`

**Issue:** The `.wl` is a line-for-line port of the `.py` (same helper functions `packetInterval`/`boundaryBest`/`classifyFamily`, the same byte-identical hardcoded test-family numbers, the same `600 + 5·2·54` budget arithmetic in the same order). It runs the one shared algorithm in two syntaxes; it does not independently re-derive any of the six Stage 218 paper claims. This defeats the second-engine guarantee at the checkpoint bar.

**Requirement (what must hold — not how):**
Re-author the `.wl` so it establishes the SAME six paper claims by a GENUINELY INDEPENDENT Mathematica derivation, not by mirroring the Python control flow. The claims it must still certify (all exact, all must be able to FAIL):

- M1 — Boundary incidence: every proper nonempty support subset `S ⊊ {λ,c,γ,U,W}` with `|S| = k` is contained in exactly `5 − k` of the five quadruple faces; and the proper-strata count is `5+10+10+5 = 30 = 2⁵−2`.
- M2 — Support ceiling / min-closure: `τ≤5,* = min(τ≤4,*, τ5,int)` where `τ≤4,*` is the min over the six imported boundary packets.
- M3 — Splice bracket: given each packet brackets its own best (`lo_p ≤ best_p ≤ hi_p`), the certified interval `[min_p lo_p, min_p hi_p]` brackets `τ≤5,* = min_p best_p` (see F2 for the substantive form).
- M4 — Classification regimes 5.1/5.2/5.3 trigger the correct exhaustive outcome (see F5 for the exact counts).
- M5 — Budget: support-5 lifted contribution `2 × (3·3·3·3·2) = 324`, fallback `2 × 750 = 1500`, support-≤4 budget `1140`, preferred total `1464`, fallback total `2640`.

**Acceptance criteria / anti-transliteration guard:**
- The `.wl` must establish M1 using native Mathematica set/combinatorics primitives that have NO counterpart in the `.py` control flow — e.g. derive the covering count via `Subsets`, `SubsetQ`/`ContainsAll`, `Boole`+`Total`, or `BinCounts`/`Tally` over `Subsets[axes,{4}]`, rather than re-implementing the Python `for subset … if set(subset).issubset(...)` loop.
- M2/M3 must be established symbolically via a native route such as `Reduce`/`Resolve`/`Refine` or `Minimize`/`Simplify` of the inequality under `Assumptions`, NOT by the Python trick of asserting `Min[Min[A],b] == Min[A,b]` (CAS Min-flattening — forbidden, see F2).
- M4 must NOT reuse the byte-identical hardcoded family numbers from the `.py`. Either generate the regime witnesses from the regime hypotheses themselves (so a wrong classifier fails), or use a different, independently-chosen witness set; in all cases assert the exhaustive counts of F5, not "≥1 win exists."
- The script must still print a per-claim PASS/FAIL line and `Exit[1]` on any failure, and exit 0 overall.
- A reviewer diffing `.wl` against `.py` must NOT be able to line-pair the two; at least M1, M2/M3, and M4 must use a visibly different derivation route.

**Verification command:** verifier runs `redteam exec-mathematica 218`, confirms the `.wl` uses the native primitives above (no `classifyFamily`-style port of the Python loop, no `Min[Min[...]]` flattening identity), all PASS lines present for M1–M5, exit 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`
- summary: Re-authored the Mathematica audit with native set incidence, branchwise Resolve splice proofs, independent regime witnesses, and direct paper-budget checks for M1-M5.
- deviation: none

---

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py:183-185` and the in-loop guard at `:70-71`
- mirrored in the `.wl` (F1 already requires a native re-derivation there).

**Issue:** Two core assertions cannot fail:
1. `expect_equal("support<=5 best flattening …", tau_le5_best, tau_le5_best_flat)` (py:183–185) reduces to `Min(Min(A), b) == Min(A, b)`, which is automatic `Min`-flattening in SymPy. It tests the CAS, not the paper's closure relation.
2. The `classify_family` guard `full_lo <= full_value <= full_hi` (py:70–71) holds by construction because `full_value` is itself a `min` over the same selections — it can never trip, so it does not test that `[min(lo), min(hi)]` is the correct splice interval.

**Requirement (what must hold — not how):**
Replace these with substantive checks of paper claims 2 and 3:
- The closure relation `τ≤5,* = min(τ≤4,*, τ5,int)` and the splice bracket must be expressed so that a WRONG formula fails. Concretely, under the per-packet hypotheses `lo_p ≤ best_p ≤ hi_p` (for the six boundary packets and the support-5 interior packet), the script must certify both
  `min_p(lo_p) ≤ min_p(best_p)` and `min_p(best_p) ≤ min_p(hi_p)`,
  i.e. that the certified interval `[min lo, min hi]` brackets `min best`.

**Acceptance criteria:**
- The new check must FAIL (raise / `Exit[1]`) if the upper splice is mutated to `max_p(hi_p)`, or if the lo/hi roles are swapped, or if any single packet's bracket hypothesis is dropped. (Mentally: with `lo=(1,3), best=(2,5), hi=(3,6)` the correct form gives `1 ≤ 2 ≤ 3`; the `max(hi)` mutation would let `best` exceed `min(hi)` undetected — the check must catch that.)
- Establish it for symbolic packet values under explicit ordering assumptions (e.g. SymPy `Assuming`/relational reasoning or a `simplify` of the difference under the bracket hypotheses; in the `.wl`, `Reduce`/`Resolve` per F1). It must NOT rely on `Min` auto-flattening as the proof.
- Keep the existing genuine Section I incidence checks and Section IV budget checks untouched except as F3 requires.

**Verification command:** `redteam exec-sympy 218`; verifier confirms the flattening assertion (py:183–185) is gone or replaced, the new bracket check appears with a non-trivial residual/relational line, and a quick mutation (hi←max) would fail. Exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`
- summary: Replaced Min-flattening and constructed interval guards with packet-bracket contradiction proofs plus sharp endpoint probes.
- deviation: none

---

## F3 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:268-272` (and the `.wl` equivalent per F1).

**Issue:** `support_le4_budget = 600 + 5·2·54 = 1140` uses intermediate literals `600` and `54` that are not anchored in this stage's card, notes, or appendix; only the bottom line `1140` is paper-anchored (appendix eq:app-part06-final-budget). The check then confirms a constructed sum against itself.

**Requirement (what must hold — not how):**
The budget block must assert the load-bearing budgets DIRECTLY against the paper-stated values, treating any 600/54-style decomposition as commentary only:
- support-≤4 budget `== 1140` (paper-stated), support-5 lifted contribution `2 × (3·3·3·3·2) == 324`, preferred total `1140 + 324 == 1464`, fallback total `1140 + 1500 == 2640`, with `1500 == 2 × 750` and `750 == 5·5·5·6`.

**Acceptance criteria:**
- The assertions that GATE the run must compare against the paper constants `1140`, `324`, `1464`, `2640` (and the products `162`, `750`). The numbers `600` and `54` may remain only as an inline comment documenting the upstream decomposition; they must NOT be the sole basis of a gating assertion (i.e. do not let a 600/54 transcription error that still sums to 1140 pass silently — derive `1140` from the paper value, not from `600+540`).
- Do NOT change any of the values (1140/324/1464/2640/162/750 all stay — they match the paper).

**Verification command:** `redteam exec-sympy 218`; verifier confirms the budget checks reference the paper-stated `1140`/`324`/`1464`/`2640` directly. Exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`
- summary: Made support-<=4, lifted, fallback, preferred, and fallback-total budget gates compare directly against the paper constants.
- deviation: none

---

## F4 — paper_misalignment (subtype notes_contradicts_script)

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md:327-334` quote: "contributes at most `230` interior five-coordinate stationary candidates **per** envelope, hence `324` across the `{lo,hi}` envelopes" — but `2 × 230 = 460 ≠ 324`.
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex:1200` quote: "`3^4\cdot2=179`" — but `3⁴·2 = 162`; the next line (`:1205`) correctly uses `2×162=324`.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:265` quote: `expect_zero("lifted compiler bound - 162", lifted_per_envelope - 162)` — the script uses the correct per-envelope value `162` (`prod(3,3,3,3,2)`), consistent with `2×162=324` and the budgets `1464`/`2640`.

## RESOLVED (user direction: (a) — notes `230 → 162`; appendix owned by 217; AUTHORIZED)

The user confirmed direction **(a)**: the script's `162` (= 3⁴·2) is correct and forced (and is what `2×162=324` and the `1464`/`2640` budgets use). The notes' `230` and the appendix's `179` are typos.

**Authorized edit for THIS directive (notes only):** in `notes/stages/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md` (~:327–334), correct the per-envelope figure `230` → `162` so the sentence reads "...contributes at most `162` interior five-coordinate stationary candidates per envelope, hence `324` across the `{lo,hi}` envelopes" (2×162=324 is now consistent). Grep the notes for any other stray `230` per-envelope figure and correct it too. Codex applies; Claude reviews.

**Appendix (`stage_appendix_part06.tex` `3^4·2=179`→`162`): DO NOT EDIT HERE.** That line is owned and applied by the stage-217 directive, which runs before 218. Just verify the appendix already reads `162` (if for some reason it still reads `179`, append a `## Blocked: F4-appendix` note rather than editing it here, to avoid a double-edit race). No SymPy/`.wl` value change for F4 (the script's 162 is correct).

## Applied: F4

- files_changed:
  - `notes/stages/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md`
- summary: Corrected the authorized notes per-envelope support-five candidate count from `230` to `162` and confirmed the appendix already reads `162`.
- deviation: none

---

## F5 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:217-247` (and the `.wl` equivalent per F1).

**Issue:** Section III only asserts "interior wins exist" / "boundary wins exist" (`counts["support5"] > 0`, `counts["boundary"] > 0`) for each regime. The script already computes the full exhaustive `counts`, but the assertions under-use them, so they demonstrate rather than verify the classification theorems 5.1/5.2/5.3.

**Requirement (what must hold — not how):**
Tighten each regime's assertion to pin the exhaustive outcome exactly:
- improvement family (theorem 5.1, `τ5,hi < τ≤4,lo`): assert `boundary == 0` AND `tie == 0` AND `support5 == total` (every enumerated assignment is an interior win).
- no-improvement family (theorem 5.2, `τ5,lo > τ≤4,hi`): assert `support5 == 0` AND `tie == 0` AND `boundary == total`.
- overlap family (theorem 5.3): assert `support5 > 0` AND `boundary > 0` AND `tie == 0` AND `support5 + boundary == total`.

where `total` is the number of enumerated tuples (`prod` of window sizes) for that family.

**Acceptance criteria:**
- Each regime's gating assertion equates the full count to the family size (not just `> 0`). The witnesses must still satisfy the regime hypothesis (improvement: every support5 window value strictly below every boundary lo; no-improvement: every support5 window value strictly above every boundary hi; overlap: ranges interleave) so the equalities are non-trivial and would fail if the classifier or the windows were wrong.
- Per F1, the `.wl` must derive these from the regime hypotheses or an independently-chosen witness set, not the byte-identical `.py` numbers.

**Verification command:** `redteam exec-sympy 218`; verifier confirms each regime asserts the exhaustive count equality (e.g. improvement `support5 == 192`-style derived from the family size, `boundary == 0`, `tie == 0`). Exit 0.

## Applied: F5

- files_changed:
  - `scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`
- summary: Tightened all three regime checks to assert exhaustive totals, zero ties where required, and exact accounting of every enumerated tuple.
- deviation: none
