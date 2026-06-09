---
unit_id: 212
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.md]
  paper_appendix: present
---

# Audit unit 212 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_212.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows: L55 summary, L988-1027 up-to-three-coordinate splice, L991 pair/triple counts, L1025 budget)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.txt`

## What the paper claims

Stage 212 is the global ranking complement of Stage 211: it splices each Stage 211 primitive-triple interior interval to the imported Stage 209 pairwise-boundary intervals, producing one certified closed-simplex interval per primitive triple. `\stagefield{Output}`: "Support-\(\le3\) finite search ledger, primitive-triple ranking theorem, global three-coordinate improvement/no-improvement theorem, and finite evaluation budget \(600\)." Distinct deliverables (from notes Section "main outputs" + Sections 2-7): (1) combinatorial ledger — \(\binom52=10\) pairs, \(\binom53=10\) triples, each pair in exactly 3 triples; (2) exact boundary-splice theorem \(\tau_{T,\min}^{\rm lo,\triangle}\le\tau_{T,*}^{\rm best,\triangle}\le\tau_{T,\min}^{\rm hi,\triangle}\) via nested-Min flattening; (3) local triple classification (interior/boundary/ambiguous order theorems); (4) primitive-triple interval-ranking partial order; (5) global three-coordinate improvement / no-improvement theorems vs. the pairwise ledger; (6) support-\(\le3\) splice sandwich; (7) finite evaluation budget \(120+480=600\). The card states "Mathematica audit: none yet" (L11).

## What the script claims to verify

The SymPy script verifies, in four blocks: (I) the combinatorial ledger via `itertools.combinations` plus exhaustive incidence counts (pair→3, axis→6) and `binomial` cross-checks; (II) nested-Min flattening `Min(iota, Min(a,b,c)) == Min(iota,a,b,c)` for the lo/hi boundary splice; (III) the five interval/order theorems (local sandwich, interior/boundary-certified order, triple ranking, global splice sandwich, improvement/no-improvement) by **exhaustive enumeration over finite integer grids** (3136/924/462/896/3136/924/462 ordered samples); (IV) the budget arithmetic `10*2*6=120`, `10*2*24=480`, sum `600`. The Mathematica script mirrors the deliverables but proves block-III over the FULL reals via `Resolve[ForAll[..., Implies[...], Reals]]`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| 10 pairs / 10 triples / pair-in-3 | py I (L52-59), wl M1 (L61-64) | match |
| boundary-splice nested-Min | py II (L89-90), wl M2 (L68-79) | match |
| local full-simplex sandwich theorem | py III.1 (L98-111), wl M3a (L83-96) | match |
| interior/boundary-certified order | py III.2 (L114-133), wl M3b (L98-126) | match |
| primitive-triple ranking partial order | py III.3 (L136-147), wl M3c (L128-141) | match |
| support-\(\le3\) global splice sandwich | py III.4 (L150-163), wl M3d (L143-156) | match |
| improvement / no-improvement theorems | py III.5 (L166-185), wl M3d (L158-186) | match |
| budget 120 / 480 / 600 | py IV (L192-204), wl M4 (L190-200) | match |
| card L11 "Mathematica audit: none yet" | a passing `.wl` is present | **mismatch** (card stale vs. retrofit) |

`paper_alignment: partial` — all mathematical deliverables match; the lone discrepancy is the card's stale "none yet" verification-status line, which is a paper-side prose status, not a math defect.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 52-53 | `expect_zero(len-binomial)` | (1) counts | yes |
| A2 | sympy | 55-65 | incidence raise-on-!=3/!=6 | (1) incidence | yes |
| A3 | sympy | 89-90 | `expect_true` Min flatten lo/hi | (2) splice | yes |
| A4 | sympy | 98-185 | finite-grid enumeration, raise on violation | (2)-(6) interval/order | yes (finite-sample) |
| A5 | sympy | 202-204 | `expect_zero(budget-120/-480/-600)` | (7) budget | yes |
| A6 | math | 61-64 | `expectExact`/`expectTrue` counts+incidence | (1) | yes |
| A7 | math | 68-79 | `===` Min flatten + fail-guard | (2) | yes |
| A8 | math | 83-186 | `Resolve[ForAll[...,Reals]]===True` x7 | (2)-(6) over Reals | yes (full-domain) |
| A9 | math | 198-200 | `expectExact` budget 120/480/600 | (7) | yes |

All rows are non-tautological: counts are computed by one route (`combinations`/`Subsets`) and checked against an independent route (`binomial`/`Binomial`); the interval theorems are quantified inequalities that can genuinely fail for a mis-stated bound; the budget literals are checked against the product decomposition `len(pairs)*2*6` etc. (not against themselves).

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (stale verification-status line)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_212.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.wl:1-206`

**What's wrong:**
The card's `\stagefield{Verification}` line states (L11): "SymPy audit: ... . Mathematica audit: none yet." But a full Mathematica audit exists and PASSES — `mathematica/moving_throat_pde_stage212_..._mathematica_audit.wl` (206 lines, M1-M4) with a fresh passing transcript (`mathematica/output/..._mathematica_audit.txt`, "Stage 212 Mathematica audit passed."). This `.wl` was added in the pass-1 dual-engine retrofit; the card's status line was not updated to reflect it.

**Why this matters:**
The card understates the unit's verification coverage and would mislead a reader (or a downstream coverage tracker) into thinking the stage is single-engine. No math is wrong; it is a documentation-status mismatch.

**Required change:**
(paper-side, user-routed — see `## Resolve before fix_loop` in the directive). Codex must NOT auto-edit the card. Direction is the user's: update L11 to cite the present `.wl`.

**Verification:**
Card L11 should name `mathematica/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.wl` instead of "none yet."

## Independent-derivation check (Mathematica)

**Verdict: INDEPENDENT.** The load-bearing objects are the interval/order theorems (py Section III, wl M3). The discriminating operation differs by method:

- SymPy extracts them by **finite-grid enumeration**: nested integer loops `for beta_lo_v in range(0,6): ...` testing the inequality at 3136/924/462/896 discrete sample points (py L98-185). This SAMPLES a finite integer subset; it does not decide the universal statement.
- Mathematica extracts them by **symbolic quantifier elimination over the reals**:
  ```
  m3a = Resolve[ForAll[{betaLo,betaHi,iotaLo,iotaHi,bBest,iBest},
        Implies[... , Min[betaLo,iotaLo] <= Min[bBest,iBest] <= Min[betaHi,iotaHi]]], Reals];
  ```
  (wl L83-96, repeated for M3b/M3c/M3d at L98-186). `Resolve[ForAll[..., Reals]]` is a real-closed-field decision procedure — it proves the inequality for ALL reals, a strictly different (and stronger) operation than enumerating an integer grid.

This is the derive-vs-posit / different-method discriminator: same physical premise and same interval definitions (allowed), but the METHOD extracting the load-bearing universally-quantified result genuinely DIFFERS (sample-the-grid vs. decide-over-Reals). The "each CAS runs its own simplifier" defense is not invoked here — the engines use categorically different proof strategies. The trivially-shared blocks (M1 combinatorics via `combinations`/`Subsets`; M2 Min auto-flattening; M4 arithmetic `10*12`, `10*48`) are non-load-bearing scaffolding/definitions and share-of-premise is allowed; no monomial-Jacobian or Series-port pattern is present. Not a port.

## Engine cross-check

Both transcripts agree at the claimed level. SymPy: "All Stage 212 identities and interval theorems verified." (counts 10/10, incidence 3/6, Min-flatten True, all interval grids pass, budget 120/480/600). Mathematica: every M1-M4 check "PASS", "Stage 212 Mathematica audit passed." (counts 10/10, incidence {3..}/{6..}, M2/M3 all `True`, budget 120/480/600). No residual, sign, or factor disagreement. `engines_agree: true`.

## Verdict justification

The math holds. I attacked: (a) tautology — counts are computed one way and checked against an independent route, interval theorems are falsifiable quantified inequalities, budget literals are checked against a product decomposition, none are `x==x`; (b) finite-sample weakness of the SymPy block-III — it only enumerates integer grids, which alone would be `insufficient_verification` for a real-domain claim, BUT the Mathematica `Resolve[...,Reals]` fully closes the universal statement, so the unit as a whole proves the claim; (c) transliteration — the `.wl` uses a categorically different proof method (quantifier elimination over Reals) for the load-bearing theorems, so it is INDEPENDENT; (d) note-typo regression (188→120 and the 246/243/245/247→212/209/211/213 renumber) — re-confirmed HELD: no stray `188`, no stale `243/245/246/247`, all stage references are the renumbered 209/211/212/213, budget `120` correctly present. I read the paper card, notes, and appendix; the script's verified claims match the paper's deliverables exactly. The single finding is a paper-side stale verification-status line ("Mathematica audit: none yet" while a passing `.wl` exists), which routes to the user; hence `verdict: findings`, `paper_alignment: partial`, no stop-cold.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 9 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| #pairs = 10 | py L49/out L12; wl M1/out L14 | notes L50 `\binom52=10`; appendix L991; card L9 (triples side) | MATCH |
| #triples = 10 | py L50/out L13; wl M1/out L16 | card L9 `\binom53=10`; notes L59; appendix L991 | MATCH |
| pair incidence = 3 | py L57/out L16-25; wl out L12 | notes L62 ("exactly three primitive triples") | MATCH |
| nested-Min splice form (lo/hi) | py II/out L40-42; wl M2/out L26-28 | notes L135-147 / appendix L1009-1014 | MATCH (symbolic deliverable) |
| pairwise budget = 120 | py L198/out L60; wl L194/out L52 | notes L429 `10\times12=120`; appendix L1025 | MATCH |
| triple interior budget = 480 | py L199/out L61; wl L195/out L54 | notes L451/L467 `10\times48=480`; appendix L1025 | MATCH |
| full budget = 600 | py L200/out L62; wl L196/out L54 | card L15 Output "budget 600"; notes L460/L515; appendix L1025 | MATCH |
| per-envelope pair count = 12 (=2x6) | py L192-193/out implicit | notes L424-426 ("12 ... six lower + six upper"); appendix L1025 | MATCH |
| per-envelope triple count = 24/48 | py L194-195/out implicit | notes L444-447 ("24 ... per envelope, hence 48"); appendix L982-985, L1025 | MATCH |

INTERNAL (scaffolding, no prose expected): axis incidence = 6 (combinatorial corollary, not a stated deliverable); the literal multipliers `2` (lower+upper envelopes) and `6`/`24` per-envelope constants that exist only to drive the 120/480 products (the products themselves are reconciled above); pass/fail flags; sample counts 3136/924/462/896 (verification-loop metadata).

Every emitted deliverable value reconciles to the card, notes, or appendix. The note-typo corrections (188→120; renumber to 209/211/212/213) HOLD.

## Self-test notes

Checked: (1) variable independence — no `diff`/`D` derivatives in this unit (pure order/interval logic), trap N/A. (2) symmetry/parity — no unbounded integrals, trap N/A. (3) trivial-case — the interval theorems hold at boundary samples (e.g. degenerate `betaLo=betaHi`), confirmed by the grids passing and by `Resolve[...,Reals]===True`. (4) path specs — no missing-script directive (both engines present). (5) paper round-trip — the only finding is paper-side and user-routed; no script edit prescribed, so no risk of introducing a new misalignment. The lone finding is a stale card status line, routed to the user.
