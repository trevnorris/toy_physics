---
unit_id: 069
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26T17:40:48Z
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
    - notes/stages/moving_throat_pde_stage069_final_reduced_verdict.md
  paper_appendix: present
---

# Audit unit 069 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_069.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage069_final_reduced_verdict.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 116, 256, 318 reference this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.txt`

## What the paper claims

The Stage 069 card's `\stagefield{Output}` is "the final pre-Family-1 support/source
verdict \eqref{eq:app-stage069-three-zone}" — the boxed three-zone classification:
`W_wall < Pe_req/Delta_inf` ⇒ universal fail; `W_wall > Pe_req/Delta_0` ⇒ universal
matched success; `Pe_req/Delta_inf ≤ W_wall ≤ Pe_req/Delta_0` ⇒ profile-sensitive band.
The card's body text adds that the sech–Gaussian benchmark (Stages 067–068) shows
realistic smooth mismatch need not significantly enlarge the middle band, but the actual
PDE branch still has to supply its own profiles. The notes file (Section 3) elaborates
that the only *truly* profile-sensitive sub-bands are
`(Pe_req/Delta_inf, P_res * Pe_req/Delta_inf)` on the failure side and
`(Pe_req/Delta_0, P_res * Pe_req/Delta_0)` on the success side, each of relative width
`P_res − 1 ≈ 0.56%`. The appendix row 116 summarizes it as "Universal fail, matched
success, and profile-sensitive middle band." This is a consolidation stage: the inputs
come from Stage 066 (the matched window) and Stage 068 (`P_res`), and the deliverable is
the unified three-zone statement plus the explicit side-band widths.

## What the script claims to verify

Both scripts (SymPy and Mathematica) assemble the three-zone verdict from upstream
symbols `Pe_req`, `Delta_0`, `Delta_gap = Delta_inf − Delta_0`, and the resonance penalty
parameter (`Pres_gap` in SymPy; `Cres2Prim` in Mathematica). They verify: (1) the matched
window endpoints `[Pe_req/Delta_inf, Pe_req/Delta_0]` arise from a generating function
`W_match(Delta_eff) = Pe_req/Delta_eff` evaluated at the endpoints, and that this
function is monotonically decreasing in `Delta_eff` (genuine derivative check); (2) the
resonance penalty `P_res` recovered as a band-edge ratio `Wfail_res/Wfail_match` matches
`(1 + Pres_gap)`, and the failure- and success-side ratios agree with each other; (3)
the two resonance-shifted thresholds equal `Pe_req/(C_res^2 Delta_inf)` and
`Pe_req/(C_res^2 Delta_0)`; (4) the matched window width equals
`Pe_req(Delta_inf − Delta_0)/(Delta_0 Delta_inf)` (partial-fraction identity); (5) the
two side-band widths equal `Pe_req(1 − C_res^2)/(C_res^2 Delta)` and the relative-width
identity `P_res − 1 = (1 − C_res^2)/C_res^2`; (6) all expected positivities. The
docstring is explicit that this is a consolidation script and the upstream derivations
(Stages 066, 068) are not re-checked here.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check(s) | Status |
|---|---|---|
| Universal-fail edge `W_wall < Pe_req/Delta_inf` (paper eq. three-zone) | `Wfail_match = Pe_req/Deltainf`; `matched fail edge from W_match(Delta_inf) = 0` (sympy L99–102 / math L103–106) | match |
| Universal-success edge `W_wall > Pe_req/Delta_0` | `Wsuff_match = Pe_req/Delta_0`; `matched success edge from W_match(Delta_0) = 0` (sympy L103–106 / math L107–110) | match |
| Three-zone partition structure | Combined effect of `Wsuff_match − Wfail_match > 0` (sympy L119 / math L132) plus the monotonicity check on `W_match(Delta_eff)` (sympy L107–110 / math L111–114) | match |
| `P_res ≈ 1.005612...` carry-forward (Stage 068) | `P_res` carried symbolically as `1 + Pres_gap` (sympy) or `1/Cres2Prim` (math); cross-ratio assertions anchor it to two independent ratios; numeric value not asserted here (it lives upstream in Stage 068) | match (consolidation-appropriate) |
| Profile-sensitive failure side-band of width `P_res − 1` (notes Section 3) | `delta_fail = Pe_req(1−Cres2)/(Cres2 Deltainf)` (sympy L151–154 / math L140–143); relative-width identity (sympy L159–162 / math L148–151) | match |
| Profile-sensitive success side-band | `delta_succ = Pe_req(1−Cres2)/(Cres2 Delta_0)` (sympy L155–158 / math L144–147) | match |
| Resonance threshold form `Pe_req/(C_res^2 Delta)` | `Wfail_res − Pe_req/(Cres2 Deltainf) = 0` and success analog (sympy L133–140 / math L88–95) | match |
| Side-band interpolation wedge (paper does not explicitly require this) | `W_failure_band = Wfail_match + u_fail * delta_fail` with `u_fail ∈ (0,1)`; positivity wedge (sympy L167–173 / math L155–161) | extra (harmless scaffolding consistent with side-band geometry) |

Front-matter `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 99–102 | `expect_zero("matched fail edge from W_match(Delta_inf)")` | matched-fail edge | yes (generating-function form ↔ endpoint identity) |
| A2 | sympy | 103–106 | `expect_zero("matched success edge from W_match(Delta_0)")` | matched-success edge | yes |
| A3 | sympy | 107–110 | `expect_positive("W_match decreasing in Delta_eff")` | matched-threshold monotonicity | yes (genuine derivative check; value is `1`) |
| A4 | sympy | 113 | `expect_zero("P_res − 1/C_res^2")` | linkage between P_res and C_res^2 | partial (true by construction of `Cres2 = 1/Pres`) |
| A5 | sympy | 114–117 | `expect_zero("matched window width")` | window-width partial-fraction identity | yes (substantive algebraic check) |
| A6 | sympy | 118 | `expect_positive("Delta_inf − Delta_0")` | gap positivity | yes (symbol-assumption echo) |
| A7 | sympy | 119 | `expect_positive("matched success threshold − matched fail threshold")` | window-width positivity | yes |
| A8 | sympy | 123–127 | `expect_zero("P_res from band-edge ratio matches (1 + Pres_gap)")` | resonance ratio identity | partial (Wfail_res is defined as Pres·Wfail_match, so ratio equals Pres by construction; still useful as linkage) |
| A9 | sympy | 128–132 | `expect_zero("P_res from success-band ratio agrees with failure-band ratio")` | failure ratio ↔ success ratio agreement | yes (independent ratios computed from two threshold pairs) |
| A10 | sympy | 133–140 | `expect_zero` resonance edges minus `Pe_req/(C_res^2 Delta)` | resonance threshold form | yes |
| A11 | sympy | 141–142 | `expect_positive("1 − C_res^2")`, `expect_positive("P_res − 1")` | resonance reduces coherence | yes (uses Pres_gap > 0 assumption) |
| A12 | sympy | 151–158 | `expect_zero` side-band widths minus `Pe_req(1−Cres2)/(Cres2 Delta)` | side-band width identity | yes (substantive algebra) |
| A13 | sympy | 159–162 | `expect_zero("P_res − 1 − (1−Cres2)/Cres2")` | relative-width identity | partial (algebraic identity true once `Cres2 = 1/Pres`) |
| A14 | sympy | 163–164 | `expect_positive` failure- and success-side widths | side-band positivity | yes |
| A15 | sympy | 167–173 | `expect_positive` wedge — failure/success-band point strictly between matched and resonance edges | side-band interior is non-empty | yes (monotone wedge for u ∈ (0,1)) |
| B1 | math | 84 | `expectPositive["1 − C_res^2"]` | C_res^2 < 1 | yes (uses Cres2Prim < 1 assumption) |
| B2 | math | 85 | `expectPositive["P_res − 1"]` | Pres > 1 | yes |
| B3 | math | 86 | `expectZero["P_res − 1/C_res^2"]` | P_res ↔ C_res^2 linkage | partial (true by `Pres = 1/Cres2` definition) |
| B4 | math | 87 | `expectPositive["Delta_inf − Delta_0"]` | gap positivity | yes |
| B5 | math | 88–95 | `expectZero` resonance edges − `PeReq/(Cres2 Delta)` | mirror of A10 | yes |
| B6 | math | 97 | `expectZero["Pres-PresGap consistency via Solve"]` | Solve round-trip on `Pres == 1 + presGapFree` | partial (tautological by construction, but routed via Solve rather than substitution to distinguish from SymPy path) |
| B7 | math | 103–110 | matched fail/success edges via `WMatchGen` | mirror of A1/A2 | yes |
| B8 | math | 111–114 | `expectPositive["W_match decreasing in Delta_eff"]` | mirror of A3 | yes |
| B9 | math | 117–126 | band-edge ratio cross-check | mirror of A8/A9 | yes |
| B10 | math | 128–131 | `expectZero["matched window width"]` | mirror of A5 | yes |
| B11 | math | 132 | `expectPositive["matched success threshold − matched fail threshold"]` | mirror of A7 | yes |
| B12 | math | 140–147 | side-band width identities | mirror of A12 | yes |
| B13 | math | 148–151 | `expectZero["P_res − 1 − (1−Cres2)/Cres2"]` | mirror of A13 | partial |
| B14 | math | 152–153 | side-band width positivities | mirror of A14 | yes |
| B15 | math | 155–161 | side-band wedge positivities | mirror of A15 | yes |

The "partial" rows correspond to algebraic identities that follow once `Cres2 = 1/Pres`
is defined (or vice versa, in Mathematica's reverse parameterization). None are
tautological in the strict "x = expr; assert x == expr" sense — each tests an identity
that combines at least two definitions. The substantive consolidation checks
(matched-window generating function, monotonicity in `Delta_eff`, the
partial-fraction matched-window width, the band-edge ratio cross-check between failure
and success sides, and the side-band width identities) are anchored.

## Findings

None.

## Independent-derivation check (Mathematica)

The Mathematica script is **not** a transliteration of the SymPy script. The two engines
use genuinely different primitive parameterizations of the resonance penalty:

- **SymPy** (L67–72): primitive `Pres_gap`, derive `Pres = 1 + Pres_gap`, then
  `Cres2 = 1/Pres`. Resonance threshold `Wfail_res = Pres * Wfail_match` (multiplication).
- **Mathematica** (L35–47): primitive `Cres2Prim` with assumption `0 < Cres2Prim < 1`,
  derive `Pres = 1/Cres2`, then extract `PresGap` via
  `Solve[Pres == 1 + presGapFree, presGapFree]`. Resonance threshold
  `WfailResViaCres2 = WfailMatch / Cres2` (division). The body comment at lines 43–45
  states this explicitly: "Parameterize via Cres2 as primitive; derive Pres = 1/Cres2
  and verify it equals (1 + PresGap) by Solve. This routes the resonance-penalty
  identity through a different algebraic operation than the SymPy script."

Three corresponding sections justify the independence:

1. SymPy L72: `Cres2 = sp.simplify(1 / Pres)` (Cres2 derived from Pres). Mathematica
   L46–47: `Cres2 = Cres2Prim; Pres = FullSimplify[1/Cres2, ...]` (Pres derived from
   Cres2; the primitive direction is reversed).
2. SymPy L78–79: `Wfail_res = sp.simplify(Pres * Wfail_match)` (multiplication by Pres).
   Mathematica L66–67: `WfailResViaCres2 = FullSimplify[WfailMatch / Cres2, ...]`
   (division by Cres2). Algebraically equivalent under `Pres = 1/Cres2`, but routed
   through a different intermediate quantity.
3. SymPy L123–127 builds `Pres_from_ratio = Wfail_res/Wfail_match` and `expect_zero`s
   the difference from `(1 + Pres_gap)` via direct simplification. Mathematica L48–53
   uses `Solve[Pres == 1 + presGapFree, presGapFree]` to extract `PresGapFromSolve`
   and stores it as `PresGap`; then L97 asserts a "Pres-PresGap consistency via Solve"
   check. The same identity is verified through a different solver call.

The two engines therefore verify that the `(Pres, Cres2, Pres_gap)` triple is mutually
consistent regardless of which member is taken as primitive. This is the right shape of
engine independence for a checkpoint consolidation stage.

## Engine cross-check

Both saved outputs (sympy `.txt` and mathematica `.txt`) show every assertion passing.
Selected residuals/forms agree under the substitution `Pres_gap = (1 − Cres2Prim)/Cres2Prim`:

| Quantity | SymPy form | Mathematica form |
|---|---|---|
| Matched fail threshold | `Pe_req/(Delta_0 + Delta_gap)` | `PeReq/(Delta0 + DeltaGap)` |
| Matched success threshold | `Pe_req/Delta_0` | `PeReq/Delta0` |
| Resonance fail threshold | `Pe_req*(Pres_gap + 1)/(Delta_0 + Delta_gap)` | `PeReq/(Cres2Prim*(Delta0 + DeltaGap))` |
| Resonance success threshold | `Pe_req*(Pres_gap + 1)/Delta_0` | `PeReq/(Cres2Prim*Delta0)` |
| Failure-side width | `Pe_req*Pres_gap/(Delta_0 + Delta_gap)` | `(PeReq − Cres2Prim*PeReq)/(Cres2Prim*(Delta0 + DeltaGap))` |
| Success-side width | `Pe_req*Pres_gap/Delta_0` | `(PeReq − Cres2Prim*PeReq)/(Cres2Prim*Delta0)` |

Under `Pres_gap = (1 − Cres2Prim)/Cres2Prim` (equivalently `Cres2Prim = 1/(1 + Pres_gap)`),
all rows are identical: `1 + Pres_gap = 1/Cres2Prim`, `Pres_gap = (1 − Cres2Prim)/Cres2Prim`.

Output mtimes (sympy `1779501932`, math `1779501936`) are both newer than their script
mtimes (sympy `1779501776`, math `1779501844`). No staleness. No engine disagreement.

## Verdict justification

This is a checkpoint consolidation stage and both engines do the consolidation correctly.
The matched window endpoints `[Pe_req/Delta_inf, Pe_req/Delta_0]` are checked via a
parameterized generating function `W_match(Delta_eff) = Pe_req/Delta_eff` with
endpoint-equality and an actual `-d/dDelta_eff` derivative used to confirm monotonicity.
The matched-window width is verified as the partial-fraction identity
`Pe_req(Delta_inf − Delta_0)/(Delta_0 Delta_inf)`, a real algebraic check. The resonance
penalty `P_res` is anchored from two independent ratios (failure-side band-edge ratio
and success-side band-edge ratio) plus an internal consistency cross-check that the two
agree. The side-band widths and the relative-width identity
`P_res − 1 = (1 − C_res^2)/C_res^2` are derived and checked. The two engines use
genuinely different primitive parameterizations (SymPy: `Pres_gap` primitive, Cres2 derived;
Mathematica: `Cres2Prim` primitive, Pres derived via `Solve`), and the algebra is routed
through different operations (substitution vs. `Solve` extraction; multiplication vs.
division). Both saved outputs are fresh and every assertion passes.

I tried to break the script along these lines and could not find a defect:

- **Tautologies**: looked for `x = expr; assert x == expr` patterns. The matched-window
  generating-function check is structurally close but is justified — it tests that the
  endpoint values reproduce the generating function, which is the actual paper form of
  the matched-window theorem (Stage 066's `\eqref{eq:app-stage066-Wwindow}`). The Solve
  round-trip in Mathematica (B6) is a minor sanity check, not a load-bearing assertion.
- **Hardcoded `P_res ≈ 1.005612`**: not present. The numeric value lives upstream
  (Stage 068's `\eqref{eq:app-stage068-Pres}`), and the consolidation script correctly
  carries the penalty symbolically. The docstring (sympy L24–37) makes this provenance
  explicit and notes that the assertions are conditional on Stages 066 and 068.
- **Missing branches**: the script covers both the failure-side and the success-side
  resonance bands plus the side-band interpolation wedge. The paper's three-zone
  partition is the union of the verified pieces.
- **Symbol-assumption errors**: all symbols are positive real; the Mathematica side adds
  `0 < Cres2Prim < 1` consistent with `C_res^2 < 1` from Stage 067. No contradictory
  declarations.
- **Paper misalignment**: the paper card's `Output` matches the script's three-zone
  structure; the notes' refined side-band statement matches the explicit width
  assertions; the appendix row matches the verdict summary.
- **Transliteration**: the two engines pick different primitives and use different
  algebraic routings, so the line-by-line port concern from the previous audit no
  longer applies.

The checkpoint bar (both engines required, paper alignment exact, substantive assertions,
constants traced to prior stages) is met. The previous audit's three findings
(F1 tautological, F2 hardcoded provenance, F3 transliteration) have all been remedied
in the current state of the scripts as recorded in the prior directive. Verdict: clean.

## Self-test notes

(1) **Variable independence**: For `sp.diff(W_match_generator, Delta_eff)` (sympy L109)
and `D[WMatchGen[DeltaEff], DeltaEff]` (math L113), the function genuinely depends on
the differentiation variable, so the derivative is `−Pe_req/Delta_eff^2`, not zero. The
combined expression `−d/dDelta_eff * Delta_eff^2 / Pe_req` evaluates to `1` (positive),
confirmed by both saved outputs. No identically-zero-derivative trap.
(2) **Symmetry/parity**: No integrals over symmetric domains appear; only algebraic
identities and a single derivative. No parity trap.
(3) **Trivial-case pre-check**: For `Wsuff_match − Wfail_match −
Pe_req(Delta_inf − Delta_0)/(Delta_0 Delta_inf)` (sympy L114–117), substituting
`Pe_req=1, Delta_0=1, Delta_gap=1` gives `1/1 − 1/2 − 1·1/(1·2) = 1/2 − 1/2 = 0`,
matching the assertion. For the side-band width identity with the same substitution
and `Pres_gap = 1` (so `Cres2 = 1/2`), `delta_fail = 1·1/2 = 1/2` and
`Pe_req(1−Cres2)/(Cres2 Deltainf) = 1·(1/2)/((1/2)·2) = 1/2` — match.
(4) **Paper round-trip**: no fixes prescribed (clean verdict), so no risk of introducing
a new `paper_misalignment`.
