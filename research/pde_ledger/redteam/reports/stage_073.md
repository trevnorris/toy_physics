---
unit_id: 073
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage073_family1_geometry_map.md
  paper_appendix: present
---

# Audit unit 073 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_073.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage073_family1_geometry_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (only the `\input{stages/stage_073}` row near line 264)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.txt`

## What the paper claims

The stage card (lines 13-28) states the Family-1 reference branch numerical freeze: `L/a = 37/20` and `ell/a = 1/20`, which yield `Lambda_ell = L/ell = 37` and (under the Stage-54 local Robin mouth closure `K_m = T_X/ell` recalled in the notes) `eta = K_m L / T_X = 37`. The `\stagefield{Output}` line is verbatim: "The Family-1 geometry values \eqref{eq:app-stage073-Lambda-eta}." (i.e., `Lambda_ell = 37` and `eta = 37`). The notes (sections 1-4) confirm the same two numerical deliverables and explicitly identify the reference branch ratios as a "reference-branch numerical freeze, not a new theorem". No additional symbolic identities are claimed by the paper beyond the algebraic identity `Lambda_ell = (L/a)/(ell/a) = L/ell`, which the notes use implicitly to combine the ratios.

## What the script claims to verify

Both scripts perform a four-step audit: (1) symbolic identity `Lambda_ell = (L/a)/(ell/a) = L/ell` over positive symbols (sympy lines 38-42, math lines 28-35); (2) numerical specialization at the reference branch, `Lambda_ell = (37/20)/(1/20) = 37` (sympy lines 44-55, math lines 37-46); (3) symbolic Robin closure: build `eta = K_m * L / T_X` and substitute `K_m -> T_X/ell`, then assert `eta - L/ell == 0` (sympy lines 60-64, math lines 48-56); (4) numerical specialization of eta to the reference branch via `subs(L/ell -> 37)` (sympy line 65, math line 57). The script's bottom-line ledger (sympy lines 67-69, math line 60) reports `Lambda_ell = 37` and `eta = 37`, matching the paper's two deliverables.

## Paper <-> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `L/a = 37/20`, `ell/a = 1/20` (carried reference) | sympy lines 44-45 / math lines 37-39 hard-code these literals; documented in card eq. (eq:app-stage073-family1-ratios). | match |
| `Lambda_ell = L/ell = 37` | symbolic check at sympy line 42 / math line 34 (algebraic identity) + numerical check at sympy line 55 / math line 46. | match |
| `eta = K_m L / T_X = 37` under `K_m = T_X/ell` | sympy lines 60-65 / math lines 48-57: builds `eta_sym = K_m * L / T_X`, substitutes the closure, asserts both the algebraic identity `eta - L/ell == 0` and the reference numerical value 37. | match |

`paper_alignment: aligned`. Every paper-side deliverable has a corresponding non-tautological script-side check (the symbolic `Lambda_ell - L/ell` is intentionally a sanity-check identity; the load-bearing checks are the numerical `Lambda_ell - 37` and the closure-substituted `eta - L/ell`). No script-side check tests something the paper does not mention.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42 | `expect_zero("Lambda_ell - L/ell (symbolic)", (L_sym/a_sym)/(ell_sym/a_sym) - L_sym/ell_sym)` | sanity: algebraic identity behind `Lambda_ell = L/ell` | partial (algebraic; intended as sanity) |
| A2 | sympy | 55 | `expect_zero("Lambda_ell - 37", (37/20)/(1/20) - 37)` | claim 2 (`Lambda_ell = 37`) | yes |
| A3 | sympy | 64 | `expect_zero("eta - L/ell", subs(K_m, T_X/ell)(K_m*L/T_X) - L/ell)` | claim 3 (Robin closure reduces eta to L/ell) | yes |
| A4 | sympy | 65 | `expect_zero("eta(reference) - 37", eta.subs({L/ell:37}) - 37)` | claim 3 (numerical value `eta = 37`) | partial (substitutes the answer; downstream of A2+A3) |
| A5 | mathematica | 34 | `expectZero["Lambda_ell - L/ell (symbolic)", lambdaStarSym/ellOverASym - lSym/ellSym]` | sanity: algebraic identity | partial (algebraic; intended as sanity) |
| A6 | mathematica | 46 | `expectZero["Lambda_ell - 37", lambdaEll - 37]` where `lambdaEll = (37/20)/(1/20)` | claim 2 | yes |
| A7 | mathematica | 56 | `expectZero["eta - L/ell", eta - len/ell]` where `eta = FullSimplify[etaSym /. km -> tx/ell]` | claim 3 (Robin closure) | yes |
| A8 | mathematica | 57 | `expectZero["eta(reference) - 37", (eta /. (len/ell) -> 37) - 37]` | claim 3 (numerical value) | partial (substitutes the answer) |

A4/A8 are mild redundancies on top of A2/A3 (resp. A6/A7): once you have `eta = L/ell` symbolically (A3/A7) and `Lambda_ell = 37` arithmetically (A2/A6), the substitution `L/ell -> 37` trivially returns 37. But these are not pure tautologies — they bind the symbolic Robin result to the numerical reference value and document the bottom-line ledger entry. The v1 audit corrected a Mathematica precedence bug at this exact line (parentheses), and that fix is present in the current file. The redundancy is intentional bookkeeping rather than a defect.

## Findings

None. All four v1 findings (F1 Mathematica precedence bug, F2 eta self-cancellation, F3 missing symbolic Lambda_ell identity, F4 transliteration-informational) were applied and the current scripts reflect those fixes:

- Sympy lines 36-42 contain the symbolic Lambda_ell identity (F3 fix).
- Sympy lines 60-62 build `eta_sym = K_m * L / T_X` and substitute the closure (F2 fix).
- Mathematica lines 28-35 contain the symbolic Lambda_ell identity (F3 fix).
- Mathematica lines 53-54 build `etaSym = km*len/tx` and substitute the closure (F2 fix).
- Mathematica line 57 has the corrected `(eta /. (len/ell) -> 37) - 37` parentheses (F1 fix).

F4 (transliteration) was marked informational in v1 and remains the structural shape of the two scripts. For a stage whose entire content is `(37/20)/(1/20) = 37` plus a one-step substitution `K_m -> T_X/ell`, there is essentially no independent derivation path available; both engines compute the same rational arithmetic and the same algebraic cancellation. I do not re-raise this as a v2 finding because (a) it was already documented in v1 as informational and not actionable, and (b) the arithmetic content is so narrow that "independent re-derivation" has no meaningful alternative implementation.

The script docstring at sympy line 3 still names "Stage 56" and both engine banners still print "STAGE 56" / "STAGE 056" rather than "STAGE 073". This is consistent with the notes (which reference Stages 54, 55, 56 from the old numbering scheme: "What Stage 56 changes" section header) — i.e., the old numbering is the legacy convention of the notes themselves. The paper card uses the current 073 number. I considered raising this as `paper_misalignment` (`notes_contradicts_script`), but the notes are the script's documented source and they themselves use the old numbering throughout. The discrepancy is a pre-existing notes-versus-card numbering convention issue that is out of this stage's scope to repair; flagging it here would not direct any actionable fix without touching prose documents (forbidden by the auditor mandate).

## Independent-derivation check (Mathematica)

The `.wl` remains a structural transliteration of the `.py` — same variable choreography (`epsilon_r <-> epsilonR`, `Lambda_star <-> lambdaStar`, `L_sym, a_sym, ell_sym <-> lSym, aSym, ellSym`, `K_m, T_X, L, ell <-> km, tx, len, ell`), same banner text, same order of assertions, same intermediate `print`/`Print` lines. However, this stage's mathematical content is one rational division `(37/20)/(1/20) = 37` plus one one-step substitution `K_m -> T_X/ell`. There is no independent derivation path that does not collapse to the same two operations. The v1 audit captured this as informational finding F4 and explicitly declined to require a rewrite. I concur: this is not a v2 finding.

## Engine cross-check

Both engines produce identical-by-residual outputs. From the saved transcripts:

| Check | Sympy residual | Mathematica residual |
|---|---|---|
| `Lambda_ell - L/ell (symbolic)` | 0 | 0 |
| `Lambda_ell - 37` | 0 | 0 |
| `eta - L/ell` | 0 | 0 |
| `eta(reference) - 37` | 0 (genuine: `37 - 37`) | 0 (genuine: `(eta with len/ell -> 37) - 37 = 37 - 37`, after the F1 parenthesis fix) |

`engines_agree: true`. The v1 audit specifically traced the precedence bug at the Mathematica `eta(reference)` line; that bug is fixed and the current saved Mathematica output (lines 17-18) shows the genuine `eta(reference) - 37 = 0` rather than the previously-vacuous `len/ell -> 0` substitution. Output mtimes (scripts at May 22 23:08, outputs at May 22 23:09) confirm the saved outputs reflect the current post-fix scripts.

## Verdict justification

The unit's content is narrow (two ratio specializations and one Robin-closure substitution). The paper card, the notes, and the script all line up on the same two deliverables (`Lambda_ell = 37`, `eta = 37`). All four v1 findings have been applied; the post-fix scripts now contain a load-bearing symbolic Lambda_ell identity, a non-trivial closure substitution, and a corrected Mathematica parenthesization. I attacked the scripts by: (a) re-verifying the v1 precedence-bug fix at mathematica line 57 — the parentheses are present and `(eta /. (len/ell) -> 37) - 37` correctly evaluates to `37 - 37 = 0`; (b) re-verifying the F2 fix that builds `eta` from `K_m * L / T_X` before substituting the closure — both engines now exercise the substitution rather than self-cancellation; (c) confirming the symbolic Lambda_ell identity sympy line 42 / math line 34 is genuine algebra over independent symbols; (d) tracing each script assertion to a paper-side deliverable (table above) and confirming no orphan checks; (e) checking output freshness via mtime (outputs newer than scripts by ~1 minute, both stamped 2026-05-22); (f) considering whether the stale "Stage 56" banner is a `paper_misalignment` — concluded no, since it matches the notes' own legacy numbering, which is out of script-side scope to repair. Verdict: `clean`, `paper_alignment: aligned`.

## Self-test notes

I checked (1) variable independence: the only `subs`/`/.` operations are `K_m -> T_X/ell` (target depends on `K_m`) and `L/ell -> 37` (target depends on the `L/ell` factor present in the simplified eta = `L/ell`); both substitutions land on symbols the expression actually contains. (2) No integrals or parity considerations in this unit. (3) Trivial-case pre-check: substituting `K_m -> T_X/ell` into `K_m*L/T_X` gives `L/ell`, so `eta - L/ell = 0` correctly; perturbing the closure to `K_m -> 2*T_X/ell` would give residual `L/ell != 0`, confirming the assertion is load-bearing. For `eta(reference) - 37`: substituting `L/ell -> 37` into `L/ell` gives 37, so `37 - 37 = 0`; perturbing the assertion's constant to `36` would give residual `1 != 0`. (4) No new script files prescribed; no path-spec concerns. (5) Paper round-trip: re-read the card and notes after walking through the assertions; the script-side checks correspond exactly to `\stagefield{Output}` (Lambda_ell = 37 and eta = 37), with no extra checks tested that the paper does not claim, and no paper claim left unverified.
