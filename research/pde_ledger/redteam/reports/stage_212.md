---
unit_id: 212
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.md]
  paper_appendix: present
---

# Audit unit 212 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_212.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (read table row line 55; splice + budget subsection lines 988–1027)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card's `\stagefield{Output}` is verbatim: "Support-\(\le3\) finite search ledger, primitive-triple ranking theorem, global three-coordinate improvement/no-improvement theorem, and finite evaluation budget \(600\)." The notes expand this into six exact deliverables: (1) the combinatorial ledger for five primitive free directions (\(\#\mathfrak P_2=\binom52=10\), \(\#\mathfrak T_3=\binom53=10\), each pair in exactly three triples), (2) the boundary-splice theorem \(\tau_{T,\min}^{\rm lo,\triangle}=\min(\beta_T^{\rm lo},\tau_{ijk,\min}^{\rm lo,int})\) with \(\beta_T=\min\) over the three edge cones, (3) the local triple classification into interior-certified / boundary-certified / ambiguous, (4) the primitive-triple ranking theorem (\(\tau_{T_a}^{\rm hi}<\tau_{T_b}^{\rm lo}\Rightarrow\) every \(T_a\) winner beats every \(T_b\) winner), (5) the global improvement/no-improvement theorems against the pairwise ledger, and (6) the finite evaluation budget \(188+480=600\) (notes §7) — where the budget arises as 10 pairs × 12 evals + 10 triples × 48 evals. The appendix (line 1026) states this budget as \(10\times12+10\times48=600\).

## What the script claims to verify

The SymPy script verifies, in four sections: (I) the combinatorial ledger — it generates pairs/triples via `itertools.combinations` and checks `#pairs - binomial(5,2)==0`, `#triples - binomial(5,3)==0`, pair-incidence \(=3\), axis-incidence \(=6\); (II) the boundary-splice algebra in `sp.Min` form, asserting the nested min flattens to a flat min over all four arguments; (III) five exhaustive integer-sample order/interval theorems (local full-simplex sandwich, interior-certified order, boundary-certified order, two-triple ranking, global splice, global improvement/no-improvement); and (IV) the budget \(10\times2\times6=120\), \(10\times2\times24=480\), total \(600\), asserting each against the literal target. All checks pass (output exit 0). The script's banner mislabels the stage as "STAGE 195" although the closing line correctly says "All Stage 212 identities... verified."

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (1) Combinatorial ledger 10/10, pair→3 triples | lines 52–59 | match |
| (1a) axis→6 triples (notes "each pair belongs to three", axis count is script-extra but consistent) | lines 61–65 | match (extra detail, harmless) |
| (2) Boundary-splice nested min | lines 75–90 | partial (structural flatten only; see F-notes) |
| (3) Classification (int/bdry/amb) | lines 114–133 | match (order-logic) |
| (4) Primitive-triple ranking | lines 136–147 | match |
| (5) Global improvement / no-improvement | lines 166–185 | match |
| (5b) Global up-to-3 splice sandwich | lines 150–163 | match |
| (6) Budget 120 + 480 = 600 | lines 192–204 | match to .tex/appendix; MISMATCH to notes literal "188" |
| Second engine (Mathematica) | — | missing (F1) |

Dominant pattern is `match`; one paper-side artifact (notes literal `188`) disagrees with the script's correct `120`, and the second engine is absent. Set `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 52 | `expect_zero(#pairs - binomial(5,2))` | claim 1 | yes |
| A2 | sympy | 53 | `expect_zero(#triples - binomial(5,3))` | claim 1 | yes |
| A3 | sympy | 58–59 | `if incidence != 3: raise` | claim 1 (pair→3 triples) | yes |
| A4 | sympy | 63–65 | `if incidence != 6: raise` | claim 1 (axis→6, extra) | yes |
| A5 | sympy | 89 | `full_lo == Min(iota_lo,tij,tik,tjk)` | claim 2 | partial (structural) |
| A6 | sympy | 90 | `full_hi == Min(iota_hi,...)` | claim 2 | partial (structural) |
| A7 | sympy | 108–109 | local full-simplex sandwich raise | claim 2/5b | yes |
| A8 | sympy | 123–124 | interior-certified order raise | claim 3 | yes |
| A9 | sympy | 129–130 | boundary-certified order raise | claim 3 | yes |
| A10 | sympy | 145 | two-triple ranking raise | claim 4 | yes |
| A11 | sympy | 160–161 | global splice sandwich raise | claim 5b | yes |
| A12 | sympy | 176 | improvement order raise | claim 5 | yes |
| A13 | sympy | 182 | no-improvement order raise | claim 5 | yes |
| A14 | sympy | 202 | `expect_zero(pair_eval_total - 120)` | claim 6 | yes |
| A15 | sympy | 203 | `expect_zero(triple_eval_total - 480)` | claim 6 | yes |
| A16 | sympy | 204 | `expect_zero(full_eval_total - 600)` | claim 6 | yes |

A5/A6 are "partial": SymPy auto-flattens nested `Min`, so once `full_lo = Min(iota_lo, Min(...))` is constructed, equality with the flat `Min` is essentially structural. It still confirms the splice definition collapses to a single min, so it is a weak but legitimate sanity check, not a hard tautology that hides a physics error.

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Subtype:** missing_mathematica
**Files:**
- `(missing)` — expected `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.wl`

**What's wrong:**
This is a non-status-only, non-checkpoint stage (`is_status_only_candidate: False`, checkpoint `False`). The card's `\stagefield{Verification}` states "Mathematica audit: none yet." Per the project dual-engine contract and the rendered prompt (line 118: "both scripts are required, and missing scripts are findings"), a second engine is required wherever Mathematica CAN independently verify the math. Every claim this stage proves is independently verifiable in Mathematica: the combinatorics via `Subsets`/`Binomial`/`Count`, the nested-min splice via Mathematica's native `Min` flattening, the four interval/order theorems via symbolic quantifier elimination (`Reduce`/`Resolve`), and the budget by direct arithmetic. There is no genuine impossibility, so the gap is a finding rather than a legitimate single-engine carve-out.

**Why this matters:**
The interval/order and combinatorial theorems are the load-bearing content of this stage; with only one engine they rest on a single implementation. A second independent derivation (using a different decomposition — symbolic quantifier elimination instead of finite integer enumeration) is what the project's anti-transliteration second-engine policy exists to provide.

**Required change:**
See directive F1: create the `.wl` against the claim manifest M1–M4, independently derived, not a port of the `.py`.

**Verification:**
`redteam exec-mathematica 212` runs the new `.wl`, it exits 0, and the new checks for M1–M4 appear and pass.

### F2 — paper_misalignment

**Severity:** low
**Subtype:** notes_contradicts_script
**Files:**
- notes `/var/projects/toy_physics/.../moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.md` (actual path `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.md`), §7.1
- sympy `.../moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.py:202`

**What's wrong:**
The notes §7.1 render the pairwise budget as `10\times 12 = 188`, then state `188+480=600` in §7.3. The arithmetic is internally inconsistent: \(10\times12=120\), not \(188\); and \(188+480=668\), not \(600\). The script (line 202) and the appendix (line 1026, `10\times12+10\times48=600`) and the stage card (`finite evaluation budget 600`) all use the correct \(120\). So only the notes literal `188` is wrong; it is a transcription typo on the paper side that the script does NOT follow.

Notes quote (§7.1):
> So the full pairwise ledger costs at most \(\boxed{10\times 12 = 188}\) exact candidate evaluations.

Script quote (line 202):
> `expect_zero("pairwise budget - 120", pair_eval_total - 120)`  with `pair_eval_total = len(pairs) * 2 * pair_eval_per_envelope = 10*2*6 = 120`

**Why this matters:**
The script is correct, so verification still passes, but the notes text states a number that does not reconcile with its own total or with the card/appendix. Left alone, a reader of the notes sees a contradictory budget table. Codex must not touch notes/, and the script is already correct, so this is routed to the user.

**Required change:**
None for Codex. See `## Resolve before fix_loop` in the directive — the user decides whether to correct the notes literal `188` → `120` (paper-side prose edit, outside red-team scope).

**Verification:**
No script change. If the user authorizes the notes edit, §7.1 reads `10\times12=120` and §7.3 reads `120+480=600`.

### F3 — paper_misalignment

**Severity:** low
**Subtype:** notes_contradicts_script
**Files:**
- script banner `.../moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.py:35`
- (context) notes title says "Stage 212" but body repeatedly says "Stage 246" (e.g. §intro "after Stage 246", §8/§9)

**What's wrong:**
The script's banner (line 35) prints `STAGE 195 — FULL PRIMITIVE-TRIPLE RANKING...`, while the file is the stage-212 audit and its closing line (206) correctly says "All Stage 212 identities and interval theorems verified." Separately, the notes body refers to this stage as "Stage 246" throughout (e.g. "So Stage 246 is the global ranking complement of Stage 245", §8/§9), although the notes title and the stage card name it Stage 212. The banner mislabel (`195`) is a copy-from-template residue in the SCRIPT and is safely correctable in-place; the notes "Stage 246/245/243" naming is a paper-side prose inconsistency outside red-team scope.

Script quote (line 35):
> `banner("STAGE 195 — FULL PRIMITIVE-TRIPLE RANKING AND THE UP-TO-THREE-COORDINATE SIEVE")`

**Why this matters:**
The wrong banner number breaks traceability (a reader of the saved output sees "STAGE 195" at the top of a stage-212 transcript). The notes "Stage 246/245/243" references are a naming-convention mismatch that may confuse cross-stage citation but does not affect any verified identity.

**Required change:**
Script side only (directive F2): change the banner string `STAGE 195` → `STAGE 212` at line 35. The notes naming is flagged in the directive `## Resolve before fix_loop` for user awareness; Codex must not edit notes.

**Verification:**
After Codex edits line 35, the saved output's header banner reads `STAGE 212 — FULL PRIMITIVE-TRIPLE RANKING...` and the script still exits 0.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot be assessed yet. The directive's claim manifest explicitly requires the new `.wl` to use a DIFFERENT decomposition than the `.py`: symbolic quantifier elimination (`Reduce`/`Resolve` over `Reals`) for the interval/order theorems instead of the `.py`'s finite integer enumeration, and `Subsets` for combinatorics instead of `itertools.combinations`. A line-by-line port of the integer-loop logic would be rejected as `mathematica_transliteration` at verification.

## Engine cross-check

Only one engine present; n/a. The SymPy output is fresh (output mtime 2026-05-11T12:49:24 > script mtime 2026-05-11T11:58:53) and exits 0 with all checks passing, so no `stale_output` finding.

## Verdict justification

The SymPy script substantively and non-tautologically verifies every deliverable the stage card and notes claim: the 10/10 combinatorics with correct incidence, the boundary-splice nested-min, the four interval/order theorems via exhaustive ordered-integer enumeration, and the 120+480=600 budget that matches the card and appendix. Attacks tried and failed: I checked whether the order-logic loops were vacuous (they are not — the inner sample loops execute on hundreds–thousands of samples, confirmed in the output counts 3136/924/462/896/...), whether the budget literal was a self-confirming hardcode (it is derived from `len(pairs)*2*6`, not pre-baked), and whether the Min-flatten check hides a sign/branch error (it does not — it is a weak structural check but the substantive splice logic is carried by section III). The verdict is `findings` (not clean) for three reasons: the required second engine is absent (F1, the headline issue); the notes render a self-inconsistent pairwise budget literal `188` that disagrees with the correct `120` the script uses (F2, paper-side, user-resolution); and the script banner mislabels the stage as 195 (F3, fixable script-side). None propagate downstream (no derived constant or sign changes), so no `CRITICAL_DOWNSTREAM`; the math is fully reconcilable, so no `UNFIXABLE`.

## Self-test notes

I ran the required traps against the claim manifest M1–M4 for the new `.wl`. Variable-independence: M1–M4 contain no `D[]`/derivative, so the identically-zero-derivative trap does not apply. Symmetry/parity: no unbounded integrals, so the parity trap does not apply. Trivial-case pre-check: `Length[Subsets[{a,b,c,d,e},{2}]]=10` and `={3}]=10`, `Min[i,Min[a,b,c]]==Min[i,a,b,c]` holds under Mathematica's native Min-flattening, and `10*12+10*48=600` reduces to a literal True — all manifest checks pass their trivial substitutions, and the budget round-trip matches the card/appendix (120, 480, 600), so the prescribed fixes introduce no new paper_misalignment.
