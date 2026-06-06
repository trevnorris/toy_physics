---
unit_id: 101
batch: IV.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage101_natural_source_map_reduction.md]
  paper_appendix: present
---

# Audit unit 101 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_101.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage101_natural_source_map_reduction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (row at line 1236: `\input{stages/stage_101}`; appendix is an index of `\input` lines, no per-stage prose row beyond this)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.txt`

## What the paper claims

Stage 101 is a retarded 2.5PN factorization ledger step. The card's `Derivation ledger` states: "The computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\)," and the block quote states: "On the natural point-particle source map, \(\widehat m_0\to1\), so the last obstruction is purely \(\chi_Q\)." The notes make the deliverables explicit: (a) from the exact Stage-83 factorization `mhat_0^2 chi_Q N_Q = 1`, taking `mhat_0 -> 1` gives `N_Q = 1/chi_Q`; (b) the canonical compact outgoing branch has `chi_Q = 1`, hence `N_Q = 1`; (c) defining the defect `Delta_Q := chi_Q - 1`, the exact gap is `N_Q = 1/(1+Delta_Q)`, i.e. `N_Q - 1 = -Delta_Q/(1+Delta_Q) = -Delta_Q + O(Delta_Q^2)`. The card's `Checks` checklist adds two more items: (2) higher odd terms begin beyond the 2.5PN coefficient, and (3) the outgoing l=2 DtN fingerprint vs the normalized `z = omega a/c_s` expansion — but the notes (authoritative intent for this stage) develop ONLY the factorization/defect content; (2) and (3) are carry-forward checklist items, not results this stage derives.

## What the script claims to verify

The SymPy docstring states the stage "owns Check (1) — the factorization `mhat_0^2 chi_Q N_Q = 1` keeping source, conservative, and outgoing factors separate," and explicitly carries Check (2) forward to stage 102 and Check (3) to stage 097, citing notes lines 41-51 for the upstream chi_Q=1 attribution. The assertions: (A1) substituting the candidate `NQ = 1/chiQ` at `mhat0=1` into the INPUT factorization `mhat0^2*chiQ*NQ - 1` zeroes the residual; (A2) the candidate `{mhat0=1, chiQ=1, NQ=1}` zeroes it; (A3) the candidate `{mhat0=1, chiQ=1+DeltaQ, NQ=1/(1+DeltaQ)}` zeroes it; (A4) the independently computed order-2 truncation of `1/(1+DeltaQ)-1` equals `-DeltaQ+DeltaQ^2`. The Mathematica script mirrors these with native `Solve`/`FullSimplify`/`Series` and an `expectZero` harness anchored to the same input factorization.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Factorization `mhat_0^2 chi_Q N_Q = 1`; `mhat_0->1` ⇒ `N_Q = 1/chi_Q` | sympy L26-34 / wl L33-43 (candidate `1/chiQ` zeroes input factorization) | match |
| Canonical branch `chi_Q = 1` ⇒ `N_Q = 1` | sympy L35-36 / wl L44-45 | match |
| Exact defect `N_Q = 1/(1+Delta_Q)` | sympy L46-47 / wl L53-54 | match |
| Linearization `N_Q - 1 = -Delta_Q + O(Delta_Q^2)` | sympy L38-51 / wl L47-52 | match |
| Check (2) higher odd terms beyond 2.5PN | (carried forward to stage 102, docstring L7-10) | match (carry-forward, by design) |
| Check (3) outgoing l=2 DtN fingerprint | (carried forward to stage 097, docstring L11-15) | match (carry-forward, by design) |

`paper_alignment: aligned`. The four core deliverables are each exercised by a non-tautological assertion in both engines. Checks (2)/(3) are checklist items the card forwards to other stages; the notes (authoritative for this stage's intent) do not claim 101 derives them, and the docstring documents the forward anchors with explicit provenance, so this is not a `script_missing_paper_claim`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 33 | `simplify((mhat0^2*chiQ*NQ - 1).subs({mhat0:1, NQ:1/chiQ})) == 0` | natural branch `N_Q = 1/chi_Q` | yes |
| A2 | sympy | 35 | `simplify((mhat0^2*chiQ*NQ - 1).subs({mhat0:1,chiQ:1,NQ:1})) == 0` | canonical `chi_Q=1 ⇒ N_Q=1` | yes |
| A3 | sympy | 46 | `simplify((mhat0^2*chiQ*NQ - 1).subs({mhat0:1,chiQ:1+DeltaQ,NQ:1/(1+DeltaQ)})) == 0` | exact defect `N_Q=1/(1+Delta_Q)` | yes |
| A4 | sympy | 50 | `expand(series_delta - (-DeltaQ + DeltaQ**2)) == 0` | linearization `N_Q-1 = -Delta_Q+O(Delta_Q^2)` | yes |
| A5 | mathematica | 42-43 | `expectZero[(mHat0^2*chiQ*nQ - 1) /. {mHat0->1, nQ->1/chiQ}]` | natural branch `N_Q = 1/chi_Q` | yes |
| A6 | mathematica | 44-45 | `expectZero[(...) /. {mHat0->1, chiQ->1, nQ->1}]` | canonical `chi_Q=1 ⇒ N_Q=1` | yes |
| A7 | mathematica | 51-52 | `expectZero[Expand[seriesDelta - (-deltaQ + deltaQ^2)]]` | linearization | yes |
| A8 | mathematica | 53-54 | `expectZero[(...) /. {mHat0->1, chiQ->1+deltaQ, nQ->1/(1+deltaQ)}]` | exact defect | yes |

All eight load-bearing checks substitute a CANDIDATE value into the INDEPENDENT input equation `mhat0^2*chiQ*nQ - 1` (or, for A4/A7, compare an independently computed series against the paper's stated polynomial). None re-tests the solved form against its own definition. Each would fail if the proposed reduction were wrong (e.g. a candidate `NQ = chiQ` at `mhat0=1` would leave residual `chiQ^2 - 1 != 0`), so the checks are non-tautological.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` derives the result natively: `nQSol = First[Solve[mHat0^2*chiQ*nQ == 1, nQ]]` (L33) then `FullSimplify` (L34), `Series[exprDelta, {deltaQ,0,2}]` (L48), and an `expectZero` harness (L20-24) that runs `FullSimplify[Together[Expand[...]]]` under `$Assumptions`. This is the Mathematica-idiomatic route, not a transliteration of Python control flow. Because the underlying claim is a one-line algebraic identity (`mhat0^2 chi_Q N_Q = 1` solved for `N_Q`), structural similarity between the two engines is mathematically unavoidable and is not a `mathematica_transliteration` defect — the comment at `.wl` L40-41 documents the same F1 de-tautologization (anchor to the input factorization, not to the solved form's own definition) was applied independently. No transliteration finding.

## Engine cross-check

Both engines produce identical results:

| quantity | sympy output | mathematica output |
|---|---|---|
| `N_Q` from factorization | `1/(chiQ*mhat0**2)` (L1) | `1/(chiQ*mHat0^2)` (L5) |
| `mhat0->1` | `1/chiQ` (L2) | `chiQ^(-1)` (L6) |
| `chiQ=1` | `1` (L3) | `1` (L7) |
| `N_Q - 1` | `-DeltaQ/(DeltaQ + 1)` (L4) | `-1 + (1 + deltaQ)^(-1)` (L12, same value) |
| small-`DeltaQ` series | `DeltaQ**2 - DeltaQ` (L5) | `-deltaQ + deltaQ^2` (L13, same) |
| all residual checks | `0` / AUDIT PASSED (L6-8) | `0` / PASS×4 / audit passed (L8-19) |

Engines agree at the level claimed. `engines_agree: true`.

## Verdict justification

`clean`. I read the card, notes, and appendix index before opening the scripts and built the model: this stage's owned deliverable is the factorization reduction `mhat_0^2 chi_Q N_Q = 1 ⇒ N_Q = 1/chi_Q` (natural branch), `N_Q = 1` (canonical), and the defect linearization `N_Q - 1 = -Delta_Q + O(Delta_Q^2)`. Every one of these is exercised by a non-tautological assertion in BOTH engines. Attacks tried and failed: (i) tautology — A1-A3/A5-A6/A8 substitute a candidate into the INDEPENDENT input equation, not the solved form's own definition, so a wrong candidate would leave a nonzero residual; (ii) series trap — A4/A7 compute the truncation independently and compare to the paper's stated `-Delta_Q + Delta_Q^2`, which I confirmed matches the order-2 expansion of `-Delta_Q/(1+Delta_Q)`, capturing the load-bearing linear coefficient `-1`; (iii) symbol-domain — `mHat0>0, chiQ>0` and Reals are exactly the physical setup (normalization factors are positive), and positivity is not used to hide a branch (the identities are rational and hold for all nonzero chiQ); (iv) missing-deliverable — card checks (2)/(3) are carry-forward checklist items explicitly forwarded to stages 102/097 with provenance, and the notes (authoritative intent) do not assign their derivation to this stage. Outputs are fresh (both .txt mtimes ~3h after their scripts). Both engines present, agree, and pass.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level table:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `N_Q = 1/(chi_Q mhat_0^2)` (exact factorized form) | py L27 / out L1; wl L36 / out L5 | notes L8 (`mhat_0^2 chi_Q N_Q = 1`), card L13 | MATCH (solved form of card's pinned product) |
| `N_Q = 1/chi_Q` (natural branch `mhat_0->1`) | py L28 / out L2; wl L37 / out L6 | notes L29 | MATCH |
| `N_Q = 1` (canonical `chi_Q=1`) | py L29 / out L3; wl L38 / out L7 | notes L55, card L13/L16 | MATCH |
| `N_Q - 1 = -Delta_Q/(1+Delta_Q)` | py L40 / out L4; wl L49 / out L12 | notes L73 | MATCH |
| small-`Delta_Q` series `-Delta_Q + Delta_Q^2` (order-2 trunc.) | py L41 / out L5; wl L50 / out L13 | notes L77 (`N_Q - 1 = -Delta_Q + O(Delta_Q^2)`) | MATCH (linear term matches; quadratic is the truncation's next term) |
| `Delta_Q := chi_Q - 1` (defect definition) | implicit in py L46/L23, wl L53 | notes L66 | MATCH |

INTERNAL scaffolding (no prose expected): residual `= 0` outputs for the four `expectZero`/`assert` checks, PASS flags, the "AUDIT PASSED" / "audit passed" banners, the `sol_NQ` solve handle.

reconciliation: complete; 6 deliverable values checked, 0 misaligned

## Self-test notes

Checked: (1) variable independence — no `sp.diff`/`D` in this script, only `Solve`/`series`, so no identically-zero-derivative trap. (2) parity/integrals — none present (pure algebra). (3) trivial-case — substituting the candidate `NQ=1/chiQ` into `mhat0^2*chiQ*NQ-1` at `mhat0=1` gives `chiQ/chiQ - 1 = 0` correctly, and a deliberately wrong candidate (`NQ=chiQ`) would give `chiQ^2-1 != 0`, confirming the checks can fail; the order-2 series of `-Delta_Q/(1+Delta_Q)` is `-Delta_Q + Delta_Q^2`, matching A4/A7. (4) paths — n/a, no missing-script directive. (5) paper round-trip — all six deliverable values reconcile to notes/card; no fix prescribed, so no new misalignment introduced. No directive written (zero findings).
