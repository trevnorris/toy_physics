---
unit_id: 090
batch: III.5
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage090_updated_reduced_status.md]
  paper_appendix: present
---

# Audit unit 090 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_090.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage090_updated_reduced_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 158, 324, 341, 359 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.txt`

## What the paper claims

Stage 090 is a status-consolidation checkpoint recording the clean reduced-theorem
status after Stages 088–089. Its `\stagefield{Output}` is: *"The support/source side
is finished under the minimal module. The next theorem gate is derivation of the
module itself."* The card's derivation enumerates six settled facts inside the reduced
Family-1 branch; the load-bearing quantitative one is item iv: the natural
contact-plus-pole reading gives `\(\rho_\alpha=4/3\)`, and the downstream-use line
refers to the `\(3/4+1/4\)` conservative quadrupole module. The notes add the full
locked triple — `rho_alpha = 4/3`, `zeta_req = 1/3`, `Pe_req = 0` — and state the
branch succeeds with large margin, lying strictly inside the exact Family-1 success
region with zero required transport bias. The appendix rows (158, 341, 359) carry the
same triple verbatim and mark the verdict StatusExactClosure for the support/source
side / StatusOpen for deriving the module. So the deliverables are: the locked triple
{rho_alpha=4/3, zeta_req=1/3, Pe_req=0}, the contact-plus-pole coefficients {3/4, 1/4},
and the branch-ordering claim (this triple sits strictly inside the carried Family-1
success/ceiling/zero-bias thresholds).

## What the script claims to verify

Both engines (a) re-derive `rho_alpha = 1/c_contact` and `zeta_req = c_pole/c_contact`
from the carried contact/pole coefficients `c_contact=3/4`, `c_pole=1/4`, and assert
they equal 4/3 and 1/3; (b) assert the module-normalization relation
`zeta_req = rho_alpha - 1` (equivalent to `c_pole = 1 - c_contact`); (c) declare three
carried Family-1 threshold decimals (`rho_suff^(chi)=3.46622291347846`,
`zeta_max^(F1)=2.46752922945601`, `A_F1=1.00005192880220`) and assert three strict
branch-ordering inequalities (`rho_alpha < rho_suff^(chi)`, `zeta_req < zeta_max^(F1)`,
`zeta_req < A_F1`); and (d) assert the carry-forward proxy `Pe_req = 0`. The docstring
explicitly frames this as a checkpoint-consistency audit, not a fresh derivation, and
names the upstream source for each carried threshold.

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| `rho_alpha = 4/3` (card iv, notes, appendix) | `rho_alpha = 1/c_contact`; `expect_zero(rho_alpha - 4/3)` (py:55,63 / wl:40,54) | match |
| `zeta_req = 1/3` (notes 18, appendix 341/359) | `zeta_req = c_pole/c_contact`; `expect_zero(zeta_req - 1/3)` (py:56,64 / wl:41,55) | match |
| `Pe_req = 0` (notes 19, appendix 341/359) | `expect_zero(Pe_req)` carry-forward proxy (py:104-105 / wl:64-65) | match (proxy, see Verdict) |
| `3/4 + 1/4` module coefficients (card 25, appendix) | `c_contact=3/4`, `c_pole=1/4` (py:52-53 / wl:38-39) | match |
| branch sits strictly inside success region / below ceilings (notes 5,6,21; appendix 341) | three `<` inequalities vs carried thresholds (py:88-99 / wl:57-59) | match |
| module normalization (contact+pole sum to 1) — implicit in `3/4+1/4` | `expect_zero(zeta_req - (rho_alpha - 1))` (py:65 / wl:56) | match (extra-but-correct: the one substantive own-derivation) |

`paper_alignment: aligned` — every stated deliverable has a faithful script-side
counterpart; no `mismatch`, no `missing`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 63 | `expect_zero(rho_alpha - 4/3)` | rho_alpha=4/3 | partial (arithmetic of 1/(3/4)) |
| A2 | sympy | 64 | `expect_zero(zeta_req - 1/3)` | zeta_req=1/3 | partial (arithmetic of (1/4)/(3/4)) |
| A3 | sympy | 65 | `expect_zero(zeta_req - (rho_alpha - 1))` | module normalization c_pole=1-c_contact | yes (non-tautological) |
| A4 | sympy | 88-91 | `expect_true(rho_alpha < rho_suff_chi)` | branch inside success region | yes |
| A5 | sympy | 92-95 | `expect_true(zeta_req < zeta_max_f1)` | below hard ceiling | yes |
| A6 | sympy | 96-99 | `expect_true(zeta_req < A_F1)` | below zero-bias baseline | yes |
| A7 | sympy | 104-105 | `expect_zero(Pe_req)` | Pe_req=0 | no (carry-forward proxy, labeled) |
| B1 | mathematica | 54 | `expectZero[rhoAlpha - 4/3]` | rho_alpha=4/3 | partial |
| B2 | mathematica | 55 | `expectZero[zetaReq - 1/3]` | zeta_req=1/3 | partial |
| B3 | mathematica | 56 | `expectZero[zetaReq - (rhoAlpha - 1)]` | module normalization | yes (non-tautological) |
| B4 | mathematica | 57 | `expectTrue[rhoAlpha < rhoSuffChi]` | branch inside success region | yes |
| B5 | mathematica | 58 | `expectTrue[zetaReq < zetaMaxF1]` | below hard ceiling | yes |
| B6 | mathematica | 59 | `expectTrue[zetaReq < aF1]` | below zero-bias baseline | yes |
| B7 | mathematica | 64-65 | `expectZero[peReq]` | Pe_req=0 | no (carry-forward proxy, labeled) |

A3/B3 is the load-bearing substantive own-assertion: `zeta_req = rho_alpha - 1`
reduces to `c_pole/c_contact = (1-c_contact)/c_contact`, i.e. `c_pole = 1 - c_contact`.
This CAN fail (it would fail for any contact/pole split that does not sum to 1), so it
genuinely verifies the contact-plus-pole module normalization the paper's `3/4+1/4`
asserts. A1/A2 are arithmetic on the carried coefficients (low independent value but
correctly anchor rho_alpha and zeta_req to the paper's quoted values). A4-A6 are robust
branch-ordering checks with large margins (≈2.13, ≈2.13, ≈0.67). A7 is an honestly
labeled carry-forward proxy, not a re-derivation (see Verdict).

## Findings

None. No script-side defect rises to a finding under any of the ten categories.

Adversarial attempts that FAILED to produce a finding:

- **tautological_check on A3/B3:** Tested whether `zeta_req - (rho_alpha - 1)` is
  algebraically guaranteed. It is not: it encodes `c_pole = 1 - c_contact`, a real
  normalization fact that would fail for a non-summing coefficient split. Non-tautological.
- **hardcoded_result on the three threshold decimals:** The carried decimals
  (`rho_suff_chi`, `zeta_max_f1`, `A_F1`) are NOT re-derived in stage 090 and would
  not change the pass/fail of the three inequalities (huge margins). However (i) the
  docstring (py:14-19) and the provenance tracker
  (`notes/CHECKPOINT_CONSTANT_PROVENANCE.md` §"Stage 090", lines 440-445) explicitly
  declare each as carry-forward with a named source; (ii) the values reconcile exactly
  with what the immediately-upstream sibling checkpoint stage 089 actually DERIVES
  in-script (089 output: `A_F1 = 1.000051928802195…`, `zeta_max = 2.467529229456012…`,
  `rho_suff = 3.466222913478464…`). They fall into the "carried forward with source
  anchor" bucket, which the checkpoint provenance rule permits. Not a finding.
- **insufficient_verification (checkpoint carries no substantive own-assertion):** The
  checkpoint does carry one substantive non-tautological own-derivation (A3/B3) plus
  three robust inequality checks. Not pure restatement. Not a finding.
- **value_mismatch (carried decimals vs paper / vs 089):** All carried decimals match
  089's derived values and the appendix-quoted triple. No mismatch.
- **symbol_assumption_error:** No symbol domain is load-bearing here; the only symbol-ish
  objects are exact Rationals and high-precision Floats used in numeric comparisons. None.
- **stale_output:** Outputs are FRESH — sympy .txt (2026-05-27 10:24:55) > .py
  (10:22:25); mathematica .txt (10:26:37) > .wl (10:22:21). No finding.
- **stale self-label (numbering):** The two 090 scripts' OWN self-labels are correctly
  "Stage 090" (py:3,47,50,103; wl:32,63,68 banner "STAGE 090"). No off-by-17 self-label
  in the scripts. (Cross-references to other stages — "Stage 69"/"63/64"/"62"/"075" as
  provenance sources, py:17-19,68-70,101; wl:61,65 — are DEFERRED per the numbering
  policy; noted, not findings. One stale self-label "Stage 73" lives in
  `notes/CHECKPOINT_CONSTANT_PROVENANCE.md:451`, but that is a process tracker, explicitly
  out of audit scope.)

## Independent-derivation check (Mathematica)

The `.wl` is NOT a mere transliteration. It is a short status-consolidation script, so
both engines necessarily share the same tiny choreography (declare 3/4, 1/4; form
1/cContact and cPole/cContact; three inequalities; Pe_req proxy). The Mathematica
script independently derives `rhoAlpha = 1/cContact` and `zetaReq = cPole/cContact`
(wl:40-41) rather than echoing a SymPy-computed literal — the comment at wl:34-37
explicitly notes "the Mathematica engine derives both here rather than hardcoding the
answer," which the source confirms. The normalization relation A3/B3 is computed by
each engine's own simplifier (`sp.simplify`/`sp.expand` vs `FullSimplify[Together[...]]`).
For a consolidation stage of this size, this is the appropriate degree of independence;
there is no shared intermediate-variable choreography that would constitute a
transliteration. Not a `mathematica_transliteration` finding.

## Engine cross-check

The two engines agree on every check:

| quantity | sympy output | mathematica output |
|---|---|---|
| rho_alpha | `4/3` (exact Rational) | `1.33333333333333…` (N[…,25]) |
| zeta_req | `1/3` (exact Rational) | `0.33333333333333…` (N[…,25]) |
| `rho_alpha - 4/3` | `0` PASS | `0` PASS |
| `zeta_req - 1/3` | `0` PASS | `0` PASS |
| `zeta_req - (rho_alpha - 1)` | `0` PASS | `0` PASS |
| inequality 1/2/3 | `True/True/True` | `True/True/True` PASS |
| Pe_req | `0` PASS | `0` PASS |

SymPy keeps the rationals exact; Mathematica reports the 25-digit numeric form of the
same rationals. No residual, sign, or factor disagreement.

## Verdict justification

`clean`. I read the paper card, the stage notes, and the appendix rows before opening
the scripts, then attacked the scripts against that model. The paper's deliverables —
the locked triple {rho_alpha=4/3, zeta_req=1/3, Pe_req=0}, the 3/4+1/4 coefficients,
and the branch-ordering claim — are each exercised by a faithful script-side check in
both engines, with exact value agreement and fresh outputs.

This is a CHECKPOINT, so it gets the higher bar: both engines present, exact paper
alignment, and at least one substantive non-tautological assertion of its own. It
**clears the higher bar.** The substantive own-assertion is `zeta_req = rho_alpha - 1`
(A3/B3), which verifies the contact-plus-pole normalization `c_pole = 1 - c_contact`
and is not algebraically guaranteed by construction; combined with the three robust
inequality checks, the consolidation is backed by real checks rather than pure
restatement. The carried threshold decimals are honestly labeled carry-forwards (not
smuggled literals) and reconcile exactly with the in-script derivations of the
immediately-upstream sibling checkpoint 089. The one genuinely tautological line,
`expect_zero(Pe_req)` on `Integer(0)` (A7/B7), is explicitly and accurately captioned
as a "carry-forward proxy for the locked triple value Pe_req = 0" rather than a
re-derivation; for a status-consolidation checkpoint whose Pe_req=0 is decided upstream
by the Stage-075 transport map, a labeled carry-forward proxy is the honest
representation and does not, on its own, drop the stage below the bar given that A3-A6
carry the substantive load. No `paper_misalignment`, no `UNFIXABLE`, no
`CRITICAL_DOWNSTREAM`.

## Value Reconciliation (pass-2 augmentation)

Enumerated every RESULT/deliverable value the two scripts emit (from source + committed
outputs; nothing was run). Deliverable-level table:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `c_contact = 3/4` | py:52 / wl:38; sympy out L5, math out L5 | `stage_090.tex:25` ("3/4+1/4"); appendix:341,359 | MATCH |
| `c_pole = 1/4` | py:53 / wl:39; sympy out L6, math out L6 | `stage_090.tex:25` ("3/4+1/4"); appendix:341,359 | MATCH |
| `rho_alpha = 4/3` | py:55,63 / wl:40,54; sympy out L7, math out L7 | `stage_090.tex:18`; notes:17,33; appendix:154,341,359 | MATCH |
| `zeta_req = 1/3` | py:56,64 / wl:41,55; sympy out L8, math out L8 | notes:18; appendix:341,359 (omitted from terse .tex card — present in notes ⇒ MATCH per guard) | MATCH |
| `Pe_req = 0` | py:104 / wl:64; sympy out L23, math out L24 | notes:19; appendix:341,359 | MATCH |
| `zeta_req = rho_alpha - 1` (normalization identity) | py:65 / wl:56; sympy out L11, math out L16 | implicit in `3/4+1/4` (`stage_090.tex:25`); notes contact-plus-pole framing | MATCH (relation, not a pinned number) |

INTERNAL / carried-input scaffolding (not stage-090 deliverables; account-only, no finding):
`rho_suff^(chi) = 3.46622291347846`, `zeta_max^(F1) = 2.46752922945601`,
`A_F1 = 1.00005192880220` (carried Family-1 thresholds, each declared with a named
upstream source in the docstring and provenance tracker; all three match the values
sibling stage 089 derives in-script); the three printed margins
(`2.1328895801…`, `2.1341958961…`, `0.6667185954…`); the three boolean inequality
flags (`True`); and the FINAL LEDGER prose block.

reconciliation: complete; 6 deliverable values checked, 0 misaligned.

## Self-test notes

I checked: (1) variable independence — no `sp.diff`/`D[]` appears, so the
zero-derivative trap is not in play; (2) the tautology trap on the two arithmetic
asserts (A1/A2, low-value but correct anchors) and confirmed A3/B3 is genuinely
non-tautological (`c_pole = 1 - c_contact` can fail); (3) trivial-case on A7/B7 —
`expect_zero(Integer(0))` is tautological in isolation but explicitly labeled a
carry-forward proxy, which I judged acceptable for a consolidation checkpoint that
carries substantive load in A3-A6; (4) carry-forward value integrity — the three
threshold decimals match stage 089's in-script derived values to printed precision, so
no `value_mismatch`. No directive is written (zero findings).
