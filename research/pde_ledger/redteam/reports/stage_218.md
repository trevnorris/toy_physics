---
unit_id: 218
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T03:10:39Z
verdict: findings
stop_cold: null
findings_count: 5
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md"]
  paper_appendix: present
---

# Audit unit 218 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_218.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 67, 1100–1276 read in full)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.txt`

## What the paper claims

Stage 218 is the "Final local mixed-ray closure theorem." Per `\stagefield{Output}` (stage_218.tex:15): *"Final local mixed-ray closure theorem, certified support-\(\le5\) splice interval, preferred total budget \(1464\), fallback budget \(2640\), and finite-candidate status of any remaining ambiguity."* The derivation ledger field enumerates the distinct deliverables (stage_218.tex:13) and the notes/appendix expand them: (1) the **boundary-identification theorem** `∂Δ₅⁺ = 𝒮≤4^loc` with proper-support strata counts `5+10+10+5 = 30 = 2⁵−2` and the incidence rule "a support subset of size k belongs to exactly 5−k quadruple faces" (notes §2, appendix eq:app-part06-boundary-strata-counts); (2) the **support-cardinality ceiling** `1 ≤ #supp ≤ 5` with `τ≤5,* = min(τ≤4,*, τ5,int)` (appendix eq:app-part06-support-ceiling/best-le5-true); (3) the **support-≤5 splice theorem** `τ≤5,lo = min(τ≤4,lo, τ5,lo,int)`, `τ≤5,hi = min(τ≤4,hi, τ5,hi,int)` bracketing `τ≤5,*` (appendix eq:app-part06-final-splice); (4) the **improvement / no-improvement / ambiguous-but-finite classification** (notes §5.1–5.3); (5) the **finite-candidate status** of any residual ambiguity; (6) the **evaluation-budget theorem**: preferred lifted degree pattern `(3,3,3,3,2)`, fallback `750/envelope`, support-≤4 budget `1140`, preferred total `1140+324 = 1464`, fallback `1140+1500 = 2640` (appendix eq:app-part06-final-budget/fallback-budget). The card marks the stage `\StatusExactClosure{}, \StatusNumerical{}, \StatusOpen{}` and explicitly warns the finite sieve is not proof that the PDE realizes the winning branch.

## What the script claims to verify

The SymPy script (and its line-for-line Mathematica twin) verify five blocks: (I) the five quadruple faces are the size-4 facets of the 5-axis simplex, there are 30 proper support strata, and every size-k proper subset is covered by exactly `5−k` faces (genuine incidence combinatorics); (II) symbolic `Min`-nesting identities `Min(Min(...6 packet syms), τ5) == Min(...6 packet syms, τ5)` for best/lo/hi; (III) three hardcoded integer-window "families" fed through a `min`-comparison classifier, asserting the improvement family yields only interior wins, the no-improvement family only boundary wins, the overlap family both, plus an in-loop interval-splice guard; (IV) arithmetic on hardcoded constants: `prod(3,3,3,3,2)=162`, `prod(5,5,5,6)=750`, `600 + 5·2·54 = 1140`, `1140 + 2·162 = 1464`, `1140 + 2·750 = 2640`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) boundary strata counts 30 = 2⁵−2; incidence k → 5−k faces | Section I: `len(proper_faces)`, per-subset `len(covering_faces)==5−len(subset)` | match (genuine combinatorics) |
| (2) support ceiling, `τ≤5,*=min(τ≤4,*,τ5,int)` | Section II symbolic `Min` flattening | partial — only tests CAS Min-nesting, not the min-closure relation as a theorem |
| (3) splice interval brackets `τ≤5,*` | Section II flattening + Section III in-loop guard `full_lo ≤ full_value ≤ full_hi` | mismatch/tautological — guard holds by construction of `min`; does not test the splice *formula* |
| (4) improvement/no-improvement/overlap classification | Section III three hardcoded families | partial — exercises the conditional-theorem regimes on invented data, not the certified ledger |
| (5) finite-candidate status | implicit in finite enumeration | partial — no explicit assertion |
| (6) budget: 162/750/1140/1464/2640 | Section IV products/sums of hardcoded constants | match on bottom line (1464/2640 paper-anchored); intermediate `600`,`54` unanchored |

`paper_alignment: partial` — the bottom-line budgets and the strata counts match exactly; the splice/ceiling theorems are exercised only by tautological Min-flattening and a construction-guaranteed interval guard.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 120,126 | `len(axes)-5==0` | claim 2 (ceiling) | partial (trivially 5−5) |
| A2 | sympy | 122 | quadruple supports == facets | claim 1 | yes |
| A3 | sympy | 128–137 | per-subset `len(covering)==5−k` | claim 1 (incidence) | yes (genuine) |
| A4 | sympy | 183–185 | `Min(Min(A),b)==Min(A,b)` | claim 2/3 | no — CAS Min-flatten, tautological |
| A5 | sympy | 70–71 | `full_lo ≤ min(b,s) ≤ full_hi` | claim 3 (splice) | no — guaranteed by `min` construction |
| A6 | sympy | 220–247 | family win-count classification | claim 4 | partial — invented windows |
| A7 | sympy | 265–266 | `prod(...)-162`, `prod(...)-750` | claim 6 | yes (factorization check) |
| A8 | sympy | 282–284 | `1140`/`1464`/`2640` budget sums | claim 6 | yes on bottom line; 600/54 unanchored |
| B1–B8 | mathematica | 123–298 | identical forms (`expectZero`/`expectTrue`) | same claims | same — line-for-line port |

## Findings

### F1 — mathematica_transliteration

**Severity:** high (checkpoint bar)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl:1-302`
- compare `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:1-287`

**What's wrong:**
The `.wl` is a mechanical line-for-line transliteration of the `.py`, not an independent re-derivation. Three corresponding sections:
- Helper choreography is identical: py `packet_interval`/`boundary_best`/`full_best`/`classify_family` (py:41–80) map one-to-one onto wl `packetInterval`/`boundaryBest`/`classifyFamily` (wl:40–78) with the same `min(min(values))`/`Min[Min/@Values]` structure and the same `counts = {"support5","boundary","tie"}` accumulator loop.
- The hardcoded test families are byte-identical: py `improvement_family = {"support_le3":(10,11),"omit_lambda":(12,13),...}` (py:189–215) vs wl `improvementFamily = <|"support_le3"->{10,11},"omit_lambda"->{12,13},...|>` (wl:198–224) — same invented numbers, same key order.
- The budget block is the same arithmetic in the same order: py `support_le4_budget = 600 + len(quadruple_faces)*2*54` (py:270–272) vs wl `supportLe4Budget = 600 + Length[quadrupleFaces]*2*54` (wl:285–287); same `source_stage` tags (198, 200), same degree patterns.
There is no place where Mathematica derives any of the six paper claims by a route different from SymPy (e.g. via `Subsets`/`Boole`-based incidence proof, `Reduce`/`Resolve`-based symbolic splice inequality, or `Resultant`/`Bezout`-style degree accounting). Both engines compute the same combinatorics on the same data — a single algorithm in two syntaxes. This violates the second-engine policy at the higher checkpoint bar.

**Why this matters:**
The two-engine guarantee is that an independent CAS, asked to establish the *same physical claim* by *its own* algebra, reaches the same answer. A transliteration cannot catch an error in the shared algorithm or shared invented data — it only catches syntax typos. For a closure-theorem checkpoint, that is not real cross-verification.

**Required change:**
See directive F1 — claim-manifest, route left to Codex (anti-transliteration guard naming native Mathematica primitives).

**Verification:** verifier confirms the `.wl` derivation uses native Mathematica primitives absent from the `.py` (e.g. `Subsets`/`Boole`/`Reduce`/`Resolve`), no longer reuses the byte-identical hardcoded family numbers, and still exits 0 with all PASS lines for the six paper claims.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:183-185` and `:70-71`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/...mathematica_audit.wl:192-194` and `:61-63`

**What's wrong:**
Two of the three "splice"/"ceiling" assertions cannot fail for any input:
1. py:183–185 asserts `tau_le5_best == tau_le5_best_flat` where `tau_le5_best = sp.Min(sp.Min(*6syms), tau5)` and the RHS is `sp.Min(*6syms, tau5)`. SymPy `Min` auto-flattens nested `Min`, so the residual is identically `0` regardless of the physics. This tests the CAS's `Min`-flattening, not the paper's closure relation `τ≤5,* = min(τ≤4,*, τ5,int)`. The `.wl` (wl:192–194) does the same with native `Min`.
2. py:70–71 (inside `classify_family`) asserts `full_lo <= full_value <= full_hi` where `full_value = min(boundary_value, support5_value)`, `full_lo = min(all mins)`, `full_hi = min(all maxes)`. Since `full_value` is itself a `min` over the same packet selections, `full_value >= full_lo` and `full_value <= full_hi` hold by construction for every window dictionary — the guard can never trip, so it does not test that `[min(lo), min(hi)]` is the correct certified splice interval (it would still pass if the splice used `max(hi)` for the upper edge).

**Why this matters:**
The splice theorem (paper claim 3) and the ceiling relation (claim 2) are the conceptual heart of the stage, yet they are currently "verified" by identities that hold for the CAS's `Min` semantics alone. A wrong splice formula (e.g. `τ≤5,hi := max(...)`, or pairing lo with hi) would pass unchanged.

**Required change:** see directive F2 — replace with a substantive bracket-monotonicity check stated as a claim-manifest (route to Codex).

**Verification:** the new check fails when the hi-splice is replaced by `max_p hi_p` or when lo/hi roles are swapped; passes for the correct `min`/`min` form. Output gains a non-trivial residual line.

### F3 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:268-272`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/...mathematica_audit.wl:285-287`

**What's wrong:**
The support-≤4 budget is rebuilt as `600 + 5·2·54 = 1140` (py:270–272). The bottom line `1140` is paper-anchored (appendix eq:app-part06-final-budget), but the intermediate constants `support_le3_budget = 600` and `quadruple_eval_per_envelope = 54` have no anchor in this stage's card, notes, or appendix — they are asserted only by an inline comment ("Stage 215 carries the support-<=3 ledger budget 600 ... 54 interior evaluations per envelope"). The check `support_le4_budget - 1140 == 0` therefore passes for any 600/54-style decomposition that happens to sum to 1140; it confirms a number against a constructed sum, not against a derived or paper-cited quantity. Likewise the fallback factorization `(5,5,5,6)` (py:255) is the script's own choice; the paper states only the scalar `750`.

**Why this matters:**
On a checkpoint, every load-bearing constant should trace to a paper/notes value or to an upstream verified result. `600` and `54` are dead intermediate literals; if either were a transcription error that still summed to 1140 the check would not catch it.

**Required change:** see directive F3 — anchor `1140` directly to the paper value and treat the 600/54 split as a comment, OR cite the upstream stage that verifies them. Claim-manifest; route to Codex.

**Verification:** the `support<=4` budget check asserts against the paper-stated `1140` (and `324`, the support-5 contribution `2·162`) directly; intermediate-literal decompositions are commentary only.

### F4 — paper_misalignment (subtype notes_contradicts_script)

**Severity:** low (informational — user-routed)
**Files:**
- notes `...stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md:327-334` quote: *"contributes at most `230` interior five-coordinate stationary candidates **per** envelope, hence `324` across the `{lo,hi}` envelopes"* (2×230 = 460 ≠ 324)
- appendix `stage_appendix_part06.tex:1200` quote: *"3^4·2=179"* (3⁴·2 = 162, not 179); appendix:1205 then uses `2×162=324`
- script `...sympy_audit.py:265` quote: `expect_zero("lifted compiler bound - 162", lifted_per_envelope - 162)`

**What's wrong:**
The script uses `162` per envelope (`prod(3,3,3,3,2)=162`), matching the appendix's *used* value (`2×162=324`, eq:app-part06-five-budget) and the final budgets. But two paper-side documents disagree with that figure: the notes say `230` per envelope (and its own `2×230` does not equal its stated `324`), and the appendix's eq:app-part06-five-bezout literally prints `3^4·2=179` although `3⁴·2 = 162`. So the script is correct and self-consistent; the **paper/notes carry two stale/typo'd per-envelope figures (179, 230)**. Direction of resolution is the user's call — not Codex's.

**Why this matters:**
The script matches the authoritative computed value, but a reader cross-checking against the notes (230) or the bezout line (179) would see three different per-envelope numbers. This is paper-side hygiene that should be reconciled before the checkpoint is signed off.

**Required change:** `## Resolve before fix_loop` — see directive F4. No script edit (script already matches the correct 162/324).

**Verification:** after user resolution, notes line and appendix eq:app-part06-five-bezout read `3⁴·2 = 162` (and notes per-envelope `162`); script unchanged.

### F5 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:187-247`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/...mathematica_audit.wl:196-264`

**What's wrong:**
Section III's three families are invented integer windows (`(10,11)`, `(2,3)`, `(3,7)`, …) chosen so each regime trivially produces the asserted outcome. They exercise the `classify_family` *machinery* (that `min`-comparison sorts a tuple into improvement/boundary/tie buckets) but do not test any paper claim about the actual certified ledger, since the real τ intervals are still symbolic and the PDE branch data are admittedly not yet inserted (notes §8). The classification theorems 5.1–5.3 are conditional ("IF τ5,hi < τ≤4,lo THEN …"); exercising one representative interval-set per regime is acceptable as a *machinery* test, but the assertions currently only confirm "≥1 interior win exists" / "≥1 boundary win exists" — they do not assert the *full* exhaustive count matches the regime, nor that the boundary winner is correctly the `min` over all six boundary packets. As written this is closer to a demonstration than a verification of claim 4.

**Why this matters:**
At the checkpoint bar, claim 4 (the improvement/no-improvement/overlap classification) deserves an assertion that pins the regime exactly (e.g. improvement ⇒ `boundary == 0` AND `support5 == total`) rather than merely "interior wins exist." The script already computes the counts; the assertion just under-uses them.

**Required change:** see directive F5 — tighten the regime assertions to exhaustive equalities (claim-manifest; route to Codex).

**Verification:** improvement family asserts `support5 == total` and `boundary == tie == 0`; no-improvement asserts `boundary == total`; overlap asserts both strictly positive and `support5 + boundary == total` with `tie == 0`. Output shows the exact counts being equated to the family size.

## Independent-derivation check (Mathematica)

Not independent. The `.wl` reproduces the `.py` function-by-function and line-by-line (see F1 for three quoted corresponding sections). Every claim is established by the identical algorithm and the identical hardcoded data in both engines. No native-Mathematica derivation route (`Subsets`/`Boole`/`Reduce`/`Resolve`/`Resultant`) is used. This is a `mathematica_transliteration` finding.

## Engine cross-check

Both engines exit 0 with byte-identical numeric results (strata count 30; improvement 192/0/0; no-improvement 0/192/0; overlap 64/64/0; budgets 1140/1464/2640; bounds 162/750). They "agree" — but trivially, because the `.wl` is a transliteration running the same algorithm on the same data (F1). The agreement carries no independent-verification weight.

## Verdict justification

`verdict: findings`, `stop_cold: null`. The combinatorial incidence proof (Section I, claim 1) is genuine and survives attack — the `5−k` covering-count would fail if any face support were wrong, and the exhaustive 30-strata enumeration is real. The budget bottom lines (1464/2640) match the appendix exactly. But the stage falls short of the checkpoint bar in five ways: (F1) the second engine is a transliteration, not an independent derivation; (F2) the splice and ceiling theorems — the conceptual core — are "verified" only by CAS `Min`-flattening identities and a construction-guaranteed interval guard that cannot fail; (F3) the support-≤4 budget rebuild leans on unanchored intermediate literals (600, 54); (F4) the notes and one appendix equation carry stale/typo'd per-envelope candidate counts (230, 179) that disagree with the correct 162 the script uses — a paper-side discrepancy for the user; (F5) the classification regimes are demonstrated rather than pinned to exhaustive counts. None of these is UNFIXABLE, and none changes a derived constant a downstream unit consumes (the budgets stay 1464/2640; the splice formula stays `min`/`min`), so no `CRITICAL_DOWNSTREAM`. F4 is paper_misalignment and routes to the user before any fix.

## Self-test notes

Checked the directive's prescribed fixes against the standard traps. No `sp.diff`/`D` derivatives are introduced (variable-independence and parity traps N/A — this is a combinatorial/interval stage, no integrals or calculus). Trivial-case pre-check on the proposed substantive splice claim: with per-packet brackets `lo_p ≤ best_p ≤ hi_p`, monotonicity of `min` gives `min(lo) ≤ min(best) ≤ min(hi)`, and a `max`-for-hi or lo/hi-swap mutation breaks the upper inequality — so the claim is non-tautological and the route is left to Codex. Paper round-trip: every prescribed fix re-anchors to the existing paper values (162/324/1140/1464/2640, `min`/`min` splice, `5−k` incidence) and introduces no new constant, so no new paper_misalignment; the 179/230 discrepancy is held for user resolution under F4 rather than silently edited.
