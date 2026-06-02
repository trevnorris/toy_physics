---
unit_id: 215
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem.md]
  paper_appendix: present
---

# Audit unit 215 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_215.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows for this unit: line 61 status row, lines 1029-1091 four-coordinate boundary-reduction + budget block, line 1114 five-simplex boundary identification, line 1353 stage input)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card's `\stagefield{Output}` reads verbatim: "Support-\(\le4\) finite search ledger, primitive-quadruple winner theorem, global support-cardinality-four improvement/no-improvement theorem, and preferred cumulative budget \(1140\)." The card frames the stage as splicing all primitive-quadruple interiors into the imported support-\(\le3\) boundary ledger, importing Stage 212 boundary packets, Stage 214 quadruple interior packets, and the \(\binom54=5\) primitive-quadruple family. The notes file expand this into seven deliverables: (1) the exact combinatorial ledger for the five primitive quadruples (and the face/incidence bookkeeping \(\#\mathfrak T_3=\binom53=10\), \(\#\mathfrak Q_4=\binom54=5\), four triple faces per quadruple, each triple in two quadruples); (2) the boundary-splice theorem \(\tau_{Q,\min}^{\rm lo,\square}\le\tau_{Q,*}^{\rm best,\square}\le\tau_{Q,\min}^{\rm hi,\square}\) with \(\tau^{\rm lo/hi,\square}=\min(\beta_Q^{\rm lo/hi},\tau_Q^{\rm lo/hi,int})\); (3) the local classification into interior-certified (\(\tau_Q^{\rm hi,int}<\beta_Q^{\rm lo}\)), boundary-certified (\(\tau_Q^{\rm lo,int}>\beta_Q^{\rm hi}\)), and ambiguous classes; (4) the primitive-quadruple ranking theorem (disjoint ordered intervals \(\Rightarrow\) one quadruple certifiably beats another) plus the **unique certified winner** theorem (\(\tau_{Q_\star}^{\rm hi,\square}<\min_{Q\ne Q_\star}\tau_Q^{\rm lo,\square}\)); (5) the global support-cardinality-4 improvement (\(\tau_4^{\rm hi,int}<\tau_{\le3}^{\rm lo}\)) and no-improvement (\(\tau_4^{\rm lo,int}>\tau_{\le3}^{\rm hi}\)) theorems; (6) the up-to-four-coordinate sieve \(\tau_{\le4}^{\rm lo/hi}=\min(\tau_{\le3}^{\rm lo/hi},\tau_4^{\rm lo/hi,int})\) bounding \(\tau_{\le4,*}^{\rm best}\); (7) the finite budget theorem \(54=3\cdot3\cdot3\cdot2\), \(5\times2\times54=540\), \(600+540=1140\) (with the imported \(600=10\times12+10\times48\)). The appendix (lines 1061-1091) confirms the degree pattern \((3,3,3,2)\), the Bézout product 54, and the 540/1140 budget.

## What the script claims to verify

The single SymPy script verifies, in four sections: (I) the combinatorial ledger via `itertools` — \(\#\)triples \(=10\), \(\#\)quadruples \(=5\), and the three incidence counts (4 faces/quadruple, 2 quadruples/triple, 4 quadruples/axis); (II) the Min-flattening identities `Min(iota, Min(a,b,c,d)) == Min(iota,a,b,c,d)` for both lo and hi envelopes, using symbolic `sp.Min`; (III) five exhaustive integer-lattice order theorems — the local full-simplex interval bound, the interior-certified and boundary-certified orderings, the pairwise certified-interval ranking (`U1<L2 => x<y`), the global support-\(\le4\) splice bound, and the global improvement/no-improvement orderings; and (IV) the budget arithmetic using hardcoded literals `support_le3_budget = 600` and `quad_eval_per_envelope = 54`, asserting `5*2*54 - 540 == 0` and `600+540 - 1140 == 0`. The output transcript exits 0 with all checks passing; output mtime (12:49) is newer than script mtime (11:58), so it is fresh.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Combinatorial ledger (10 triples, 5 quadruples, face/incidence counts) | Section I, lines 38-67 | match |
| Boundary-splice interval `min(beta,iota)` bounds best (notes §2) | Section II Min algebra (78-93) + Section III local full-simplex theorem (96-109) | match |
| Local classification: interior-certified / boundary-certified orderings (notes §3.1, §3.2) | Section III (111-130) | match |
| Local classification: ambiguous class is well-defined (notes §3.3) | none (definitional set complement; nothing to assert) | match (n/a) |
| Primitive-quadruple **pairwise** ranking (disjoint ordered intervals) | Section III certified-interval ordering (132-143) | match |
| Primitive-quadruple **unique certified winner** (`tau_Qstar^hi < min over other four`) | not separately tested | partial |
| Global support-4 improvement / no-improvement | Section III (160-179) | match |
| Up-to-four sieve `tau_le4 = min(tau_le3, tau_4int)` bounds best | Section III global splice (145-158) | match |
| Budget: `54 = 3*3*3*2`, `600 = 10*12+10*48`, `540`, `1140` | Section IV (181-192) | partial (54 and 600 are bare hardcoded literals; their paper-stated factor decompositions are not reconstructed) |
| Second independent engine (Mathematica) | absent | missing |

Dominant pattern: most deliverables match, but one is only partially exercised (unique-winner), the budget constants are imported literals rather than reconstructed, and the entire second engine is absent. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 48 | `expect_zero(#triples - C(5,3))` | claim 1 (combinatorics) | yes |
| A2 | sympy | 49 | `expect_zero(#quadruples - C(5,4))` | claim 1 | yes |
| A3 | sympy | 51-55 | loop `face_count == 4` | claim 1 (faces/quadruple) | yes |
| A4 | sympy | 57-61 | loop `incidence == 2` | claim 1 (triple in 2 quads) | yes |
| A5 | sympy | 63-67 | loop `incidence == 4` | claim 1 (axis in 4 quads) | yes |
| A6 | sympy | 92-93 | `Min` flattening lo/hi | claim 2 (splice algebra) | yes |
| A7 | sympy | 97-109 | lattice: `lo_full <= simplex_best <= hi_full` | claim 2 (boundary-splice interval) | yes |
| A8 | sympy | 113-122 | lattice: interior-certified `i*<b*` | claim 3 (interior-certified) | yes |
| A9 | sympy | 123-130 | lattice: boundary-certified `b*<i*` | claim 3 (boundary-certified) | yes |
| A10 | sympy | 132-143 | lattice: `U1<L2 => x<y` | claim 4 (pairwise ranking) | yes |
| A11 | sympy | 145-158 | lattice: global splice bound | claim 6 (sieve) | yes |
| A12 | sympy | 160-171 | lattice: improvement `q4*<s3*` | claim 5 (improvement) | yes |
| A13 | sympy | 172-179 | lattice: no-improvement `s3*<q4*` | claim 5 (no-improvement) | yes |
| A14 | sympy | 191 | `expect_zero(5*2*54 - 540)` | claim 7 (budget) | partial (54 hardcoded) |
| A15 | sympy | 192 | `expect_zero(600+540 - 1140)` | claim 7 (budget) | partial (600 hardcoded) |
| — | sympy | — | unique-winner `min over others` | claim 4 (unique winner) | MISSING |
| — | mathematica | — | (no second engine) | all | MISSING |

All SymPy assertions trace to a paper deliverable and are non-tautological (the lattice loops genuinely exercise the order/min monotonicity — they would fail if any min were structured wrong; the `Min`-flattening is a real symbolic identity; the combinatorial incidences are genuine). No `tautological_check`, no `symbol_assumption_error`. The two gaps are the unverified unique-winner statement and the absent second engine; the budget literals are a weak (anchored) spot rather than a standalone finding.

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Subtype:** missing_mathematica
**Files:**
- `mathematica/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_mathematica_audit.wl` (does not exist)

**What's wrong:**
The manifest records `mathematica.path: null, exists: false` for unit 215, and `is_status_only_candidate: false`, `is_checkpoint: false`. The stage computes content that is fully and independently verifiable: definite combinatorial incidence counts (`Subsets`), the `Min`-flattening identity, six order/interval theorems over the certified-interval lattice, and the budget arithmetic. None of this is a pure status/label carry-forward — every section of the SymPy script asserts a checkable mathematical statement. Therefore Mathematica genuinely CAN independently verify this stage, and the project dual-engine contract requires a `.wl`.

**Why this matters:**
A checkpoint-bar, non-status-only stage with a single engine has no cross-check. The interval theorems in particular are exactly the class of order-monotonicity claims a second engine should confirm by a different route (continuous quantifier elimination via `Reduce` rather than the SymPy integer-lattice enumeration), which also guards against a sampling-coverage blind spot in the `.py`.

**Required change:**
Codex must author a new independent Mathematica audit per the directive's claim manifest (M1-M7) and anti-transliteration guard. See directive F1.

**Verification:**
`redteam exec-mathematica 215` produces the new `.wl` output, the file exists at the target path, and the script exits 0 with all manifest claims (M1-M7) asserted via `Exit[1]`-on-failure guards.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `scripts/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_sympy_audit.py:132-143` (ranking) and `:181-192` (budget)

**What's wrong:**
Two sub-claims of the paper are only partially exercised by the SymPy script:
(a) **Unique certified winner.** The notes §4 state two distinct ranking results — the *pairwise* ordering (verified at lines 132-143) AND the *unique certified winner* theorem: "If for some primitive quadruple \(Q_\star\) one has \(\tau_{Q_\star}^{\rm hi,\square}<\min_{Q\ne Q_\star}\tau_Q^{\rm lo,\square}\), then \(Q_\star\) is the unique certified primitive-quadruple winner." The script tests only the two-quadruple pairwise version; it never tests the explicit "minimum over the other four quadruples" statement that the card's `\stagefield{Output}` calls the "primitive-quadruple winner theorem."
(b) **Budget constants.** Lines 182-183 hardcode `support_le3_budget = 600` and `quad_eval_per_envelope = 54` as bare literals; the assertions at 191-192 only check the products `5*2*54=540` and `600+540=1140`. The paper *derives* both constants — appendix line 1069 gives `3*3*3*2=54` (Bézout from degree pattern (3,3,3,2)) and line 1025 gives `10*12+10*48=600`. The script does not reconstruct either from its stated factor decomposition, so it confirms the products against numbers it asserted itself.

**Why this matters:**
(a) leaves the headline "winner theorem" of the card without a direct script-side check (the pairwise theorem implies it but is not the same statement, and the card names the multi-quadruple winner). (b) is a weak hardcoded-literal spot: if 54 or 600 were transcribed wrong from upstream, the product checks would still pass. Both are minor because the constants ARE anchored in the paper and the pairwise theorem is closely related — hence low severity, not `paper_misalignment`.

**Required change:**
The new Mathematica engine is the cleaner place to close both gaps (manifest M5 covers the explicit unique-winner statement; M7 reconstructs 54 and 600 from their factor decompositions rather than asserting them as literals). Optionally Codex may also add a unique-winner lattice check and the `54=3*3*3*2`, `600=10*12+10*48` reconstructions to the SymPy script — but the load-bearing fix is captured in the M5/M7 claim manifest so the second engine carries it. See directive F2.

**Verification:**
The new `.wl` contains an explicit `min over the other quadruples` unique-winner check (M5) and reconstructs 54 and 600 from `3*3*3*2` and `10*12+10*48` before forming 540/1140 (M7); if Codex also patches the SymPy script, a `min`-over-others unique-winner loop appears near line 143 and the budget literals at 182-183 are replaced by reconstructed products. Scripts exit 0.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot yet be assessed. The directive includes an explicit anti-transliteration guard requiring a different decomposition (continuous `Reduce`/`Resolve` quantifier elimination for the interval theorems instead of the SymPy integer-lattice enumeration; `Subsets` for combinatorics; factor-product reconstruction for the budget constants).

## Engine cross-check

Only one engine present; no cross-check possible. This is precisely the gap F1 addresses.

## Verdict justification

The SymPy script is internally sound: every assertion traces to a paper deliverable, the lattice/min order theorems are non-tautological (they would fail under a wrong min-structure), the `Min`-flattening is a genuine symbolic identity, the combinatorics are exact, and the output is fresh. Attacks tried that failed: I checked for tautology in the lattice loops (they genuinely test min-monotonicity, not a definition echo), for symbol-domain errors (the `tau_*` symbols are declared `real=True`, consistent with closure times; no positivity is needed for the order claims), for missing branches (both lo and hi envelopes are exercised; both improvement and no-improvement directions are covered), and for budget tautology (540/1140 are products of the literals, non-circular given the literals). What does not hold to the checkpoint bar: there is no second engine (F1), and two sub-claims — the explicit unique-winner statement and the factor-decomposition of the budget constants — are only partially exercised (F2). The banner says "STAGE 198" (line 35) and the notes reference Stages 246/248/249, but these are stale-label cosmetics, not verified-claim mismatches; the actual math tested is correct for Stage 215, so this is not a `paper_misalignment`. The banner mislabel on the SymPy side is folded into F2's optional script touch as a no-cost correction. Verdict: `findings`, not stop-cold — no math is unreconcilable and no downstream constant changes (the budget 1140 and the certified-interval structure are confirmed, not altered).

## Self-test notes

I ran the required traps against the M1-M7 claim manifest. (1) Variable independence: no `D[]`/`sp.diff` appears anywhere in this stage — it is order/combinatorics/arithmetic only, so the zero-derivative trap does not apply. (2) Symmetry/parity: no integrals over unbounded domains; the interval theorems are finite quantifier statements, no parity concern. (3) Trivial-case pre-check: M1 lengths reduce to literal 10 and 5; M3/M6 min-monotonicity holds at the simplest profile (all bounds equal -> all comparisons are equalities, consistent with the non-strict `<=`); M4/M5 strict orderings hold when the certifying gap is positive; M7 reduces to `54`, `600`, `540`, `1140` literally. (4) Path: target `.wl` lives under `mathematica/` with the exact sibling naming convention (`..._mathematica_audit.wl`), stated in full in the directive. (5) Paper round-trip: M7 reconstructs 54 from `3*3*3*2` and 600 from `10*12+10*48` exactly as the appendix states (lines 1069, 1025), introducing no new constant; no new `paper_misalignment` created.
