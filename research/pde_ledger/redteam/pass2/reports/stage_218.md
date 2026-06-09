---
unit_id: 218
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md"]
  paper_appendix: present
---

# Audit unit 218 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_218.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (read rows 4-5, 67, 1195-1275 referencing stage 218)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.txt`

## What the paper claims

Stage 218 is the final-closure checkpoint of the Part VI local mixed-ray sieve. Its
`\stagefield{Output}` is: "Final local mixed-ray closure theorem, certified support-`<=5`
splice interval, preferred total budget `1464`, fallback budget `2640`, and finite-candidate
status of any remaining ambiguity." The notes enumerate six exact deliverables: (1) the
**boundary-identification theorem** `\partial\Delta_5^+ = \mathcal S_{\le4}^{\rm loc}` (the five
quadruple faces and their lower-support strata, with proper-strata counts `5+10+10+5 = 30 = 2^5-2`
and the incidence law "a support-`k` subset belongs to exactly `5-k` quadruple faces");
(2) the **support-cardinality ceiling** `1 <= #supp <= 5` with
`\tau_{\le5,*}^{\rm best} = min(\tau_{\le4,*}^{\rm best}, \tau_{5,*}^{\rm best,int})`;
(3) the **support-`<=5` splice theorem** `\tau_{\le5,\min}^{\rm lo/hi} = min(...,...)` with
`\tau_{\le5,\min}^{\rm lo} <= \tau_{\le5,*}^{\rm best} <= \tau_{\le5,\min}^{\rm hi}`;
(4) the **improvement / no-improvement / ambiguous-overlap classification** (regimes 5.1/5.2/5.3);
(5) the statement that residual ambiguity is **finite-candidate**, not continuum; and
(6) the **evaluation-budget theorem**: lifted `3^4·2 = 162` per envelope → `324` across `{lo,hi}`;
imported support-`<=4` budget `1140`; preferred total `1140+324 = 1464`; fallback `750`/envelope →
`1500`, fallback total `1140+1500 = 2640`. This is a checkpoint (`is_checkpoint: True`): both engines
required, alignment must be exact.

## What the script claims to verify

The SymPy script verifies, in four sections: (I) the boundary stratification — it constructs the
five quadruple faces by omit-axis comprehension and asserts each proper nonempty support subset is
covered by exactly `5-len(subset)` quadruple faces, and that there are 30 proper strata; (II) the
splice bracket — it builds `tau_le5_* = Min(tau_le4_*, tau5_*_int)` symbolically and proves
`min(lo) <= min(best) <= min(hi)` by a per-packet propositional-contradiction check
(`simplify_logic(And(lo>best, lo<=best)) is False`) plus numeric endpoint probes; (III) the three
regime classifications by exhaustive enumeration over hand-listed witness windows; (IV) the budget
ledger (`162→324`, `750→1500`, `1140`, totals `1464`/`2640`). The Mathematica script verifies the
same five claims by independent operations (`Subsets`/`ContainsAll`/`Tally` for M1, `Resolve[ForAll]`
real-QE for the M2/M3 splice ceiling, a generator-based witness family `makeBoundaryWindows` for M4,
and the same budget products for M5).

## Paper <-> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) Boundary identification, strata counts 5+10+10+5=30=2^5-2, incidence `5-k` | py §I `proper_faces`+coverage assert; wl M1 `Subsets`+`Boole[ContainsAll]`+`Tally`+`2^5-2-30` | match |
| (2) Support ceiling `tau_le5 = min(tau_le4, tau5_int)` | py §II `tau_le5_best = Min(...)` + closure-operands assert; wl M2 `Resolve[ForAll[..., Min==...]]` | match |
| (3) Splice bracket `lo <= best <= hi` | py §II `certify_splice_bracket` per-packet contradiction + probes; wl M3 `Resolve[ForAll]` branch closures + probes | match |
| (4) Regimes 5.1/5.2/5.3 improvement/no-improvement/overlap | py §III three families + count asserts; wl M4 generator windows + count asserts | match |
| (5) Finite-candidate (not continuum) ambiguity | py §IV canonical-screen count + overlap family is finite enumeration; wl M4 finite Tuples | match (implicit via finiteness of the enumerations) |
| (6) Budget `162/324/750/1500/1140/1464/2640` | py §IV products+asserts; wl M5 products+asserts | match |

`paper_alignment: aligned`. Every deliverable has a corresponding, non-tautological script check in
both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 173-179 | `expect_zero(axes-5)`, `quad count-5`, `quadruple_supports == expected_quadruples` | (1) face set | yes |
| A2 | sympy | 181-190 | per-subset `len(covering_faces) == 5-len(subset)` else raise | (1) incidence law | yes |
| A3 | sympy | 94,106 | `len(packet_names)-7`, `hyp count-14` | (3) structural | partial (scaffolding-count) |
| A4 | sympy | 111-124 | `simplify_logic(And(lo>best, lo<=best)) is false` ×7 lower, ×7 upper | (3) splice bracket | yes |
| A5 | sympy | 128-133 | numeric `min`-endpoint probes equal least lo/hi | (3) splice min | yes |
| A6 | sympy | 239-242 | `tau_le5_closure_operands == (tau_le4_best, tau5_best_int)` | (2) ceiling structure | yes |
| A7 | sympy | 279-326 | regime family count asserts (boundary/tie/total) ×3 | (4) classification | yes |
| A8 | sympy | 343-366 | budget product + sum asserts | (6) budget | yes |
| M1 | math | 61-64 | `AllTrue[#[[3]]==#[[4]]]`, `Tally===5/10/10/5`, `#-30`, `2^5-2-30` | (1) strata+incidence | yes |
| M2 | math | 83-102,135 | `Resolve[ForAll[..., Min[x,y]==smaller], Reals]` | (2) ceiling | yes |
| M3 | math | 104-129,136-140 | `Resolve[ForAll[..., !(lo>best)]]` ×7 + probes | (3) splice bracket | yes |
| M4 | math | 184-213 | generator windows, regime count asserts ×3 | (4) classification | yes |
| M5 | math | 219-239 | budget product + sum asserts | (6) budget | yes |

A3 is a structural scaffolding count (packet-count `==7`, hypothesis-count `==14`); it confirms the
loop dimensions, not the physics, so it is "partial" — but A4/A5 carry the actual splice proof, so no
finding.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py` (mtime 2026-06-03 15:59:11)
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.txt` (mtime 2026-06-02 12:33:13)

**What's wrong:**
The SymPy `.py` was modified (commit `e2a4780`, "numbering reconciliation Phase 1: doc-only
stage-label fixes") AFTER its `.txt` output was captured (commit `a12029e`). The captured output is
therefore stale. The visible symptom: the script's current banner is
`banner("STAGE 218 — FULL SUPPORT-<=5 COMPLETION AND LOCAL MIXED-RAY SEARCH CLOSURE")` (py:136), but
the saved output line 3 still reads `STAGE 201 — FULL SUPPORT-<=5 ...`. The Mathematica output is
fresh (`.wl` mtime 2026-06-02 12:30:19 < `.txt` 2026-06-02 12:33:13).

I checked whether the staleness affects content: it does not. The numbering-reconciliation commit
only updated the banner string `201`→`218`; every numeric/logical line in the stale `.txt`
(`192`/`192`/`64-64` regime counts, all `= 0` residuals, all `True` flags, budget `1140/324/1500/1464/2640`)
matches what the current script would emit. So this is the informational `stale_output` signal, not a
content disagreement.

**Why this matters:**
A stale output banner labels the captured transcript with the wrong stage number, which is exactly the
numbering-drift attractor this project tracks. The orchestrator's independent re-run will refresh it.

**Required change:**
None for Codex to patch in the script. The orchestrator's standard fresh re-run (`exec-sympy 218`)
will regenerate the `.txt` with the corrected `STAGE 218` banner. No script edit is needed.

**Verification:**
After the orchestrator re-runs `python3` on the script, the new `.txt` line 3 reads
`STAGE 218 — FULL SUPPORT-<=5 ...` and all numeric lines are unchanged.

## Independent-derivation check (Mathematica)

The `.wl` was re-authored in pass-1 to an independent route. I attacked each load-bearing object for
shared-operation transliteration; the verdict is INDEPENDENT on the three theorem objects and a
benign shared arithmetic on the budget bookkeeping.

**Object A — boundary completion / strata counts (load-bearing): INDEPENDENT.**
SymPy builds the faces by omit-axis dict comprehension and counts coverage with Python set logic:
```python
covering_faces = [name for name, packet in quadruple_faces.items()
                  if set(subset).issubset(packet["support"])]
expected = len(axes) - len(subset)
```
Mathematica builds the faces with `Subsets[axes, {4}]` and counts incidence with native set predicates:
```wolfram
Total[Boole[ContainsAll[#, support]] & /@ quadrupleFaces]   (* incidence *)
strataTally = Sort[Tally[Length /@ properStrata]];           (* 5/10/10/5 *)
expectZero["M1 2^5 - 2 - 30", (2^Length[axes] - 2) - 30];
```
Different generators (`itertools.combinations` vs `Subsets`), different predicates
(`set.issubset` vs `ContainsAll`/`Boole`), and the `.wl` adds an independent `Tally` cross-check of
the strata histogram and the `2^5-2` identity that the `.py` does not compute by tally. Different
operation extracting the same combinatorial fact → independent.

**Object B — splice bracket / support ceiling (load-bearing, the heart of the stage): INDEPENDENT.**
This is the discriminator the prompt flagged. SymPy proves `lo <= best <= hi` by treating the
relations as boolean atoms and showing a propositional contradiction with finite logic:
```python
branch_counterexample = sp.simplify_logic(sp.And(lo > best, bracket_hypotheses[(name,"lo<=best")]))
expect_true(f"lower splice branch {name} contradicted", branch_counterexample is sp.false)
```
Mathematica proves the SAME bound by symbolic quantifier elimination over the reals:
```wolfram
m3LowerBranches = Table[Resolve[ForAll[{loVars[[i]], bestVars[[i]]},
   Implies[Element[{loVars[[i]], bestVars[[i]]}, Reals] && loVars[[i]] <= bestVars[[i]],
           ! (loVars[[i]] > bestVars[[i]])]], Reals], {i, packetCount}];
```
and proves the ceiling `Min(x,y)==smaller` itself by `Resolve[ForAll[..., Min[x,y]==x], Reals]`
(M2) — something the SymPy side never does symbolically (SymPy only asserts the operand TUPLE is the
right pair, `tau_le5_closure_operands == (tau_le4_best, tau5_best_int)`, py:240). `simplify_logic`
(finite propositional) vs `Resolve[ForAll[...], Reals]` (real QE) is exactly the
"symbolic QE vs finite enumeration" independence signature named in the brief. Different proof method
→ independent. (The two engines' numeric endpoint probes share the literal probe tuples `{1,4,5,6,7,8,9}`
and `{3,6,7,8,9,10,11}`, but those are a tiny corroborating side-check, not the bracket proof, which
is carried by the contradiction/QE machinery above.)

**Object C — regime witnesses (load-bearing classification): INDEPENDENT.**
The two engines exercise the SAME regime logic with DIFFERENT witness data, provably so from the
saved outputs: SymPy's hand-listed families give improvement/no-improvement/overlap totals
`192 / 192 / 64+64=128`; Mathematica's `makeBoundaryWindows[start,width,gap]` generator gives
`256 / 192 / 65+63=128`. The improvement totals (192 vs 256) and the overlap split (64-64 vs 65-63)
differ because the witness windows were independently chosen (SymPy: explicit per-label 2-tuples;
Mathematica: `makeBoundaryWindows[20,2,3]` / `[2,2,3]` / `ConstantArray[{4,8}]`). Same theorem,
genuinely independent witnesses → independent.

**Object D — evaluation budget (bookkeeping, NOT the load-bearing theorem): SHARED ARITHMETIC, benign.**
Both engines compute `Times[3,3,3,3,2]=162`, `2*162=324`, `Times[5,5,5,6]=750`, `2*750=1500`,
`+1140`, totals `1464`/`2640`. This is the same operation (product of the same degree patterns, plus
the same posited `1140` literal) in both engines. There is only one way to take a product, and the
`1140` is explicitly a paper-stated import (py:347-349 comment). A budget tally is not a theorem to
be "independently re-derived"; the patterns `(3,3,3,3,2)` and the products are matched to paper
eqs `app-part06-five-budget`/`app-part06-final-budget`. I do NOT raise `mathematica_transliteration`
on this object: it is a shared confirmation of a posited bookkeeping number, not a shared-operation
port of a derivation. The stage's three actual theorems (A/B/C) are independently derived.

## Engine cross-check

Both engines emit identical bottom-line PASS verdicts on every shared claim, and their numeric
deliverables agree where they should:
- strata counts: py prints 30 proper strata + per-subset `5-k` incidence; wl prints
  `tally {{1,5},{2,10},{3,10},{4,5}}` and `2^5-2-30 = 0`. Agree.
- splice bracket: py all 14 branch-contradictions `False`/contradicted=True; wl all
  `m3Lower/UpperBranches` True. Agree.
- regimes: py `192/192/64-64`, wl `256/192/65-63` — these DIFFER by design (independent witnesses),
  but the regime CONCLUSIONS agree (5.1 → all support5, 5.2 → all boundary, 5.3 → both nonzero, no
  ties in any regime). No `engine_disagreement`: the counts are not claimed to match; the regime
  logic is.
- budget: both `1140 / 324 / 1500 / 1464 / 2640`. Agree.

No `engine_disagreement` finding.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level table (script-emitted result values vs `.tex` card / `.md` notes / part appendix):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| proper-strata counts `5+10+10+5 = 30` | py out L12 `30`; wl out L11/L15 tally; wl M1 | notes L142-147 `5+10+10+5=30=2^5-2`; appendix eq L1215 | MATCH |
| `2^5-2 = 30` identity | wl out L19 `2^5-2-30=0` | notes L145; appendix L1215 | MATCH |
| incidence law (support-`k` in `5-k` faces) | py out L17-46; wl out L13 | notes L148 "size `k` belongs to exactly `5-k`"; appendix L1217 | MATCH |
| boundary id `\partial\Delta_5^+ = S_{<=4}^{loc}` | py §I coverage proof; wl M1 | notes L155 boxed; appendix eq L1220 | MATCH |
| ceiling `tau_le5,* = min(tau_le4,*, tau5,*int)` | py L220 `tau_le5_best`; wl M2 | notes L176-183 boxed; appendix eq L1230-1233 | MATCH |
| splice `tau_le5,lo/hi = min(...)`, bracket `lo<=best<=hi` | py §II; wl M3 | notes L193-223 boxed; appendix eq L1241-1257 | MATCH |
| regimes 5.1/5.2/5.3 outcomes (all-support5 / all-boundary / mixed) | py §III; wl M4 | notes §5.1-5.3 L245-289 | MATCH (qualitative) |
| lifted per-envelope `162` | py L344; wl out L81 `M5 lifted per-envelope - 162 = 0` | notes L328 `162`; appendix eq L1200 `3^4·2=162` | MATCH |
| lifted total `324` | py L363; wl out L83 | notes L332 `\boxed{324}`; appendix eq L1205 `2×162=324` | MATCH |
| fallback per-envelope `750` | py L345; wl out L85 | notes L349 `2×750`; appendix L1207 `bound 750` | MATCH |
| fallback total `1500` | py L364; wl out L87 | notes L349 `=1500`; appendix L1207 `fallback ... budget 1500` | MATCH |
| imported support-`<=4` budget `1140` | py L349; wl out L89 | notes L338 `\boxed{1140}`; appendix eq L1261 | MATCH |
| preferred total `1464` | py L365; wl out L91 | notes L343 `1140+324=1464`; `\stagefield{Output}` L15; appendix eq L1261 | MATCH |
| fallback total `2640` | py L366; wl out L93 | notes L354 `1140+1500=2640`; `\stagefield{Output}` L15; appendix eq L1266 | MATCH |
| fallback degree pattern `(5,5,5,6)` | py L334; wl L218 | not separately stated; product `750` IS stated (appendix L1207) | MATCH (in-script derivation reproduces the paper's 750) |

INTERNAL (scaffolding, no prose expected, no finding): pass/fail flags; per-packet branch
counterexample `False` flags; `splice packet count - 7`, `bracket hypothesis count - 14`,
`boundary packet count - 6` structural counts; the abstract witness window tuples themselves
(`{1,4,5,6,7,8,9}`, `makeBoundaryWindows` ranges); `canonical_screens = ('gradient-optimal','equal-mix')`;
and the dead-metadata fields `source_stage: 198/200` (never read or asserted — see Self-test notes).

reconciliation: complete; 15 deliverable values checked, 0 misaligned.

## Verdict justification

The stage holds up against the paper at the checkpoint bar. Both engines are present and substantive;
all six paper deliverables map to non-tautological checks in both; every emitted deliverable value
reconciles exactly with the card, notes, and appendix (15/15). I attacked the splice bracket for a
buried tautology (the `simplify_logic` contradiction is a real, if simple, mutual-exclusion fact on
two real symbols, not assignment-trivial), the regime asserts for branch coverage (all three regimes
exercised with the correct hypothesis guards), and — hardest — the `.wl` for residual transliteration
after the V.3-200 lesson that a re-author can still be a port. The `.wl` extracts each of the three
load-bearing objects by a genuinely different operation than the `.py` (native `Subsets`/`ContainsAll`/`Tally`
vs `itertools`+`set.issubset`; `Resolve[ForAll[...], Reals]` real-QE vs `simplify_logic` finite
propositional contradiction; generator-windows `256/192/65-63` vs hand-windows `192/192/64-64`). The
only shared-operation object is the budget product, which is a posited bookkeeping tally, not the
stage's theorem, and is correctly anchored to the paper. The single finding is an informational
`stale_output` (SymPy `.txt` banner still says STAGE 201 after the numbering-reconciliation commit
touched the `.py`; numeric content is identical), which the orchestrator's standard re-run clears. No
`paper_misalignment`, no `mathematica_transliteration`, no `engine_disagreement`. Verdict: `findings`
(one low-severity `stale_output`), checkpoint higher bar CLEARED.

## Self-test notes

I checked: (1) Variable independence — N/A here, there are no `sp.diff`/`D` derivatives; the
load-bearing operations are set-combinatorics, propositional/QE bracket proofs, finite enumeration,
and integer products, all of which I traced to confirm they exercise real facts. (2) The
`simplify_logic(And(lo>best, lo<=best))` contradiction is non-tautological (two mutually exclusive
real orderings on the same pair) and the `Resolve[ForAll[..., Reals]]` QE genuinely proves the Min
identity. (3) I re-derived the Mathematica M4 window arithmetic by hand (improvement `2^6·4=256`,
no-improvement `64·3=192`, overlap `64·2=128` split `65/63`) and confirmed it matches the saved
output and the regime hypotheses. (4) I confirmed `source_stage: 198/200` and the "Stage 249" comment
(py:145,152,331,347) are dead/cosmetic numbering-drift labels never read or asserted on — they do not
affect any verified quantity, so no `value_mismatch`; they are the known incomplete-renumber artifact,
not a deliverable. (5) Paper round-trip: the only finding is `stale_output` with no script edit, so it
introduces no new misalignment.
