---
unit_id: 215
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem.md]
  paper_appendix: present
---

# Audit unit 215 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_215.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows: line 61 table, line 236 narrative, lines 1071-1091 budget block, line 1114 five-simplex face reference)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_mathematica_audit.txt`

## What the paper claims

Stage 215 is a combinatorial / interval-arithmetic closure stage (no new constitutive law or optimizer). The card's `\stagefield{Output}` reads: "Support-\(\le4\) finite search ledger, primitive-quadruple winner theorem, global support-cardinality-four improvement/no-improvement theorem, and preferred cumulative budget \(1140\)." The notes enumerate seven deliverables: (1) the exact combinatorial ledger for the five primitive quadruples (\(\binom54=5\), \(\binom53=10\) triples, each quadruple has 4 triple faces, each triple lies in 2 quadruples); (2) the boundary-splice theorem \(\tau_{Q,\min}^{\rm lo,\square}\le\tau_{Q,*}^{\rm best,\square}\le\tau_{Q,\min}^{\rm hi,\square}\) with \(\beta_Q=\min\) over the four triple faces and full = \(\min(\beta_Q,\iota_Q)\); (3) the local quadruple classification (interior-certified \(\tau^{\rm hi,int}<\beta^{\rm lo}\Rightarrow i<b\); boundary-certified \(\tau^{\rm lo,int}>\beta^{\rm hi}\Rightarrow b<i\)); (4) the primitive-quadruple ranking + unique-winner theorem; (5) the global support-4 improvement / no-improvement theorems against the support-\(\le3\) ledger; (6) the up-to-four-coordinate sieve theorem (global \(\min\) of two intervals); (7) the finite budget \(5\times2\times54=540\) interior + imported \(600\) = \(1140\). The card's `\stagefield{Verification}` line states: "SymPy audit: ... Mathematica audit: **none yet**."

## What the script claims to verify

The SymPy script verifies, in four blocks: (I) the exact combinatorics by direct construction of `itertools.combinations`, asserting `#triples=10`, `#quadruples=5`, every quadruple has exactly 4 triple faces, every triple lies in exactly 2 quadruples, every axis in exactly 4 quadruples; (II) the symbolic `Min`-flattening identity `Min(iota, Min(t_ijk,...)) == Min(iota, t_ijk,...)`; (III) the five interval theorems (full-simplex splice, interior/boundary classification ordering, pairwise ranking, unique-winner over five intervals, global splice + improvement/no-improvement) by exhaustive integer enumeration over nested `range` loops; (IV) the budget arithmetic `5*2*54=540`, `600+540=1140`. The Mathematica script (`M1`-`M7`) mirrors the deliverable SET but proves the interval theorems (M3-M6) with `Resolve[ForAll[{reals}, Implies[premise, conclusion]], Reals]` (real quantifier elimination), and reconstructs the budget from factor data `{{10,12},{10,48}}` and `degreePattern={3,3,3,2}`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) combinatorial ledger (5/10/4-faces/2-incidence/4-axis) | py I (lines 38-67); wl M1 (52-82) | match |
| (2) boundary-splice full-simplex interval | py III (97-109); wl M3 (109-120) | match |
| (3a) interior-certified order | py III (113-122); wl M4 (124-134,146) | match |
| (3b) boundary-certified order | py III (123-130); wl M4 (135-147) | match |
| (4) ranking + unique winner | py III (132-159); wl M5 (151-191) | match |
| (5) global improvement / no-improvement | py III (176-195); wl M6 (206-230) | match |
| (6) up-to-four sieve / global splice | py III (161-174); wl M6 (195-205,228) | match |
| (7) budget 540 / 600 / 1140 | py IV (197-208); wl M7 (232-251) | match |
| Min-flattening identity (notes §2 algebra) | py II (78-93); wl M2 (84-105) | match (extra-but-supportive) |
| Card `\stagefield{Verification}`: "Mathematica audit: none yet" | a `.wl` exists and passes (M1-M7) | mismatch (F1) |

`paper_alignment: partial` — all seven physics deliverables map cleanly to substantive, non-tautological checks in both engines; the single discrepancy is the card's verification line still declaring no Mathematica audit, which a pass-1 retrofit invalidated.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 48-49 | `expect_zero(#triples-C(5,3))`, `#quadruples-C(5,4)` | claim 1 | yes |
| A2 | sympy | 54-67 | `raise` if face/incidence/axis counts != 4/2/4 | claim 1 | yes |
| A3 | sympy | 92-93 | `Min(iota,Min(...))==Min(iota,...)` | notes §2 algebra | yes |
| A4 | sympy | 106-108 | `lo<=min(b*,i*)<=hi` over 3136 samples | claim 2 | yes |
| A5 | sympy | 117-130 | interior/boundary ordering over 924/462 samples | claim 3 | yes |
| A6 | sympy | 137-159 | pairwise + unique-winner ordering (896 / 542760) | claim 4 | yes |
| A7 | sympy | 161-195 | global splice + improve/no-improve (3136/924/462) | claims 5,6 | yes |
| A8 | sympy | 207-208 | `quad_eval-540`, `full-1140` | claim 7 | yes |
| M1 | math | 78-82 | `expectZero`/`expectTrue` combinatorics | claim 1 | yes |
| M2 | math | 98-105 | `FullSimplify[nested==flat]` Reals | notes §2 algebra | yes |
| M3 | math | 109-120 | `Resolve[ForAll[...],Reals]` splice | claim 2 | yes |
| M4 | math | 124-147 | `Resolve[ForAll]` interior/boundary order | claim 3 | yes |
| M5 | math | 151-191 | `Resolve[ForAll]` pairwise + 5-interval unique winner | claim 4 | yes |
| M6 | math | 195-230 | `Resolve[ForAll]` global splice + improve/no-improve | claims 5,6 | yes |
| M7 | math | 248-251 | `expectZero` budget 54/600/540/1140 | claim 7 | yes |

Every script-side assertion traces to a specific paper-side deliverable. No orphaned scaffolding. No tautology: A4-A7 / M3-M6 test order relations that genuinely could fail for an incorrect `Min`/interval definition; A1-A2 / M1 count an independently constructed combinatorial object; A8 / M7 reduce a numeric literal to a factorized product (`3*3*3*2`, `10*12+10*48`) rather than asserting a number against itself.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_215.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_mathematica_audit.wl:1-256`

**What's wrong:**
The card's `\stagefield{Verification}` line states:
> "SymPy audit: \StageFile{scripts/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_sympy_audit.py}.  Mathematica audit: **none yet**."

But a Mathematica audit (`..._mathematica_audit.wl`, 256 lines, M1-M7) exists, was added in the pass-1 dual-engine retrofit, and passes (saved output ends "Stage 215 Mathematica audit passed."). The card therefore understates the verification status. (This is a pure prose/status discrepancy, not a math error; routed to the user because Codex may not edit `paper/`.)

**Why this matters:**
The card misrepresents coverage — it claims single-engine when the unit is now dual-engine. Low severity (no math is wrong), but it should be corrected for the card to be an accurate ledger entry, consistent with the rest of the pass-2 retrofit reconciliation.

**Resolution:** see directive `## Resolve before fix_loop`.

**Verification:**
After the user authorizes, the card line 11 should cite the `.wl` audit (e.g. "Mathematica audit: \StageFile{mathematica/...stage215..._mathematica_audit.wl}.") and the `\claimstatus`/Verification wording no longer says "none yet."

## Independent-derivation check (Mathematica)

**Verdict: INDEPENDENT.** The single discriminating operation is how each engine proves the load-bearing interval theorems (claims 2-6, the certified-interval orderings):

1. **SymPy proves by exhaustive discrete enumeration.** E.g. the full-simplex splice (py 97-109):
   ```python
   for beta_lo_v in range(0,6): ... for i_star in range(iota_lo_v, iota_hi_v+1):
       simplex_best = min(b_star, i_star)
       if not (lo_full <= simplex_best <= hi_full): raise AssertionError(...)
   ```
   It samples 3136 ordered integer tuples and checks the inequality pointwise.

2. **Mathematica proves by continuous quantifier elimination.** The corresponding M3 (wl 109-120):
   ```wolfram
   m3Theorem = Resolve[ForAll[{betaLo,betaHi,iotaLo,iotaHi,bBest,iBest},
     Implies[betaLo<=betaHi && iotaLo<=iotaHi && betaLo<=bBest<=betaHi && iotaLo<=iBest<=iotaHi,
       Min[betaLo,iotaLo] <= Min[bBest,iBest] <= Min[betaHi,iotaHi]]], Reals];
   ```
   This is a `Resolve[ForAll[...], Reals]` real-closed-field decision over the *continuum*, not a sampled loop — a strictly stronger and methodologically distinct route. Same physical premise (the same Min-of-interval definitions), but the load-bearing object (the universally-quantified ordering) is *extracted by a different method*: enumerate-and-check (py) vs. decide-over-reals (wl).

3. The same divide holds for M4/M5/M6 vs the py loops at lines 117-159 / 161-195: the unique-winner theorem in particular is sampled over 542760 integer tuples in py but proved by five `Resolve[ForAll]` calls (one per `star`) over all reals in wl (M5, 169-191). The combinatorial block (M1 vs py I) does share the same counting operation, but that is the shared *definitional* premise (the 5-element index set and its subsets), not the load-bearing extraction; sharing it is permitted. The `Min`-flattening identity (M2 vs py II) is a one-line associativity fact, also shared premise.

This is not a port: the discriminator (derive-vs-posit / different-method) is satisfied on every load-bearing object. The "each CAS runs its own simplifier" defense is not invoked here — the methods are genuinely different (discrete sampling vs. real QE).

## Engine cross-check

Both outputs are fresh (see Verdict). Final lines agree: SymPy "All Stage 215 identities and interval theorems verified."; Mathematica "Stage 215 Mathematica audit passed." with every `PASS:` present. Numeric deliverables agree exactly across engines: per-envelope 54, support-\(\le3\) 600, interior 540, full 1140; combinatorics 10/5/faces-4/incidence-2/axis-4. No residual disagreement. `engines_agree: true`.

## Verdict justification

The math holds up under attack. I tried to break it three ways: (a) tautology — A4-A7/M3-M6 are not tautological (they test order relations of independently-defined `Min` expressions over a real or integer domain that could fail for a wrong definition; A8/M7 factorize the literal budgets rather than asserting numbers against themselves); (b) insufficient sampling vs. general claim — the SymPy loops cover the relevant ordered-interval configurations including the strict-separation premises, and the Mathematica side proves the *general* (all-reals) statement via `Resolve`, so the discrete sampling is backstopped by a continuous proof; (c) transliteration — the `.wl` uses real quantifier elimination where the `.py` uses discrete enumeration, a genuinely independent method on every load-bearing object, so it is INDEPENDENT, not a port (unlike sibling 211). I read the paper card, the full notes, and the Part VI appendix rows; all seven physics deliverables reconcile with both engines. The only finding is the low-severity prose discrepancy F1 (card still says "Mathematica audit: none yet" though a passing `.wl` now exists), which requires user resolution and routes to the gate. Verdict: `findings` (1, paper_misalignment), no stop-cold.

## Value Reconciliation (pass-2 augmentation)

Enumerated deliverable values the scripts emit (from `.py`/`.wl` source + saved `.txt`):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `#triples = 10` (\(\binom53\)) | py 45/48, txt 12; wl M1 78, txt 15 | notes §1 line 51 (`\#\mathfrak T_3=\binom53=10`); appendix line 236 narrative | MATCH |
| `#quadruples = 5` (\(\binom54\)) | py 46/49, txt 13; wl M1 79, txt 17 | notes §1 line 60 / §9 line 487; card line 9 (`\binom54=5`); appendix 1071 (`\binom54=5`) | MATCH |
| each quadruple has 4 triple faces | py 52-55, txt 16-20; wl M1 80, txt 19 | notes §1 line 63 ("exactly four codimension-one faces"); card Checks line 20 | MATCH |
| each triple in 2 quadruples | py 57-61, txt 21-30; wl M1 81, txt 21 | notes §1 line 63 ("belongs to exactly two primitive quadruples") | MATCH |
| each axis in 4 quadruples | py 63-67, txt 31-35; wl M1 82, txt 23 | (combinatorial corollary; notes §1) | MATCH (notes-supported) |
| per-envelope candidate bound `54` (`3*3*3*2`) | py 199, (implicit); wl M7 236/243, txt 76 | notes §7 line 419 (`\boxed{54}`); appendix 1069 (`3·3·3·2=54`) | MATCH |
| interior quadruple budget `540` (`5*2*54`) | py 200/207, txt 67/69; wl M7 238/250, txt 78 | notes §7 lines 429-431 (`5×2×54=540`) + §7 line 447 (`540`); appendix 1074 (`5×2×54=540`) | MATCH |
| imported support-\(\le3\) budget `600` (`10*12+10*48`) | py 198/203, txt 66; wl M7 235/237/249, txt 77 | notes §7 line 436 (`\boxed{600}`) | MATCH |
| full support-\(\le4\) budget `1140` (`600+540`) | py 201/208, txt 68/70; wl M7 239/251, txt 79 | card Output line 15 (`1140`); notes §7 line 441 (`600+540=1140`) + §9 line 494 (`1140`); appendix 1079 (`600+540=1140`) | MATCH |
| degree pattern `{3,3,3,2}` | wl M7 234, txt 74 | appendix line 1064 (`(3,3,3,2)`) | MATCH |
| support-\(\le3\) factor data `{{10,12},{10,48}}` | py 198; wl M7 235, txt 75 | (intermediate factorization of the 600 budget; 10 triples × per-envelope) | MATCH (notes §7 line 436 carries 600; factor split is internal) |

Symbolic deliverables (interval theorems) — these are universal order statements, not numeric values; their forms (`\tau_{Q,\min}^{\rm lo,\square}=\min(\beta_Q^{\rm lo},\tau^{\rm lo,int})`, the splice/ranking/improvement inequalities) all appear boxed in notes §§2-7 and are exercised by py III / wl M3-M6 — accounted for as MATCH under the standard cross-check, not re-listed here.

INTERNAL (scaffolding, no finding): loop sample counts (3136, 924, 462, 896, 542760), PASS/FAIL flags, `expect_zero`/`expectZero` residuals (all 0), the per-engine `True` results from `Resolve`, banner strings, `count_*` accumulators.

`reconciliation: complete; 11 deliverable values checked, 0 misaligned.` (The one non-value discrepancy — the card's "Mathematica audit: none yet" verification line — is captured as F1, a status/prose `paper_misalignment`, not a value mismatch.)

## Self-test notes

No `sp.diff`/`D` derivatives in this unit (pure combinatorics + interval order theory), so the variable-independence trap (step 1) is N/A. No unbounded integrals, so the parity trap (step 2) is N/A. Trivial-case pre-check (step 3): the budget reductions reduce to literal 0 (`540-540`, `1140-1140`, `54-54`, `600-600`); the interval-theorem loops/Resolve premises are satisfiable and the conclusions hold for the simplest separated configurations (e.g. interior-certified with `iota_hi < beta_lo` gives `i* < b*`). Paths (step 4): N/A — F1 is a paper_misalignment, no script target. Paper round-trip (step 5): F1 prescribes no script edit (user-gated), so it introduces no new misalignment.
