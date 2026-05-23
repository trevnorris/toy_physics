---
unit_id: 069
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 069 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.txt`

## What the script claims to verify

The script docstring asserts four substantive claims: (1) the matched-branch Stage-49 window is exactly `[Pe_req/Delta_inf, Pe_req/Delta_0]` carried forward from Stage 066; (2) the Stage-51 resonance family shifts both thresholds by the exact penalty factor `P_res = 1/C_res^2` carried forward from Stage 068; (3) the profile-sensitive failure-side and success-side bands have the exact widths claimed in the paper; (4) the reduced three-zone verdict is the correct matched-branch theorem, with the resonance refinement living in two narrow side-bands of widths `P_res - 1 = (1 - C_res^2)/C_res^2`. Because this is a CHECKPOINT stage that consolidates results from earlier stages, the assertions must demonstrate that the *consolidation* (combining stage-49 and stage-51 results into the three-zone verdict) is mathematically faithful, not merely re-define the symbols and check they obey their own definitions.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 74 | `Pres - 1/Cres2 == 0` | no (Cres2 := 1/Pres by line 59) |
| A2 | sympy | 75-78 | `Wsuff_match - Wfail_match - Pe_req*(Deltainf-Delta0)/(Delta0*Deltainf) == 0` | no (algebraic identity from `Wsuff_match := Pe_req/Delta0`, `Wfail_match := Pe_req/Deltainf`) |
| A3 | sympy | 79 | `Deltainf - Delta0 > 0` | partial (positivity of `Delta_gap`; trivial from assumption) |
| A4 | sympy | 80 | `Wsuff_match - Wfail_match > 0` | partial (follows trivially from `Deltainf > Delta0` and positive Pe_req) |
| A5 | sympy | 81-84 | `Wfail_res - Pe_req/(Cres2*Deltainf) == 0` | no (Wfail_res := Pres*Pe_req/Deltainf and Pres = 1/Cres2 by definition) |
| A6 | sympy | 85-88 | `Wsuff_res - Pe_req/(Cres2*Delta0) == 0` | no (same as A5 with Delta0) |
| A7 | sympy | 89 | `1 - Cres2 > 0` | partial (Cres2 = 1/(1+Pres_gap) < 1 trivially given Pres_gap > 0) |
| A8 | sympy | 90 | `Pres - 1 > 0` | partial (Pres = 1+Pres_gap, Pres_gap > 0; trivial) |
| A9 | sympy | 99-102 | `delta_fail - Pe_req*(1-Cres2)/(Cres2*Deltainf) == 0` | no (`delta_fail := Wfail_res - Wfail_match = (Pres-1)*Pe_req/Deltainf`; using `Pres - 1 = (1-Cres2)/Cres2` recovers the asserted form by definition) |
| A10 | sympy | 103-106 | `delta_succ - Pe_req*(1-Cres2)/(Cres2*Delta0) == 0` | no (same structure as A9) |
| A11 | sympy | 107-110 | `(Pres - 1) - (1 - Cres2)/Cres2 == 0` | no (pure algebra: `(1 - 1/Pres)/(1/Pres) = Pres - 1`) |
| A12 | sympy | 111-112 | `delta_fail > 0`, `delta_succ > 0` | partial (follows from `Pres > 1` and positive Pe_req) |
| A13 | sympy | 118-121 | four convex-combination positivity checks for `W_failure_band`, `W_success_band` | partial (interior of an interval; given `0 < u_fail < 1`, `0 < u_succ < 1` from `v_fail, v_succ > 0`, this is convex-combination trivia) |
| B1 | math | 57 | `Pres - 1/Cres2 == 0` | no (identical structure to A1) |
| B2 | math | 58-61 | matched window width | no (identical structure to A2) |
| B3 | math | 62 | `Deltainf - Delta0 > 0` | partial |
| B4 | math | 63 | matched success threshold - matched fail threshold > 0 | partial |
| B5 | math | 64-67 | resonance fail threshold form | no (identical structure to A5) |
| B6 | math | 68-71 | resonance success threshold form | no |
| B7 | math | 72 | `1 - Cres2 > 0` | partial |
| B8 | math | 73 | `Pres - 1 > 0` | partial |
| B9 | math | 81-84 | failure-side width form | no |
| B10 | math | 85-88 | success-side width form | no |
| B11 | math | 89-92 | `(Pres-1) - (1-Cres2)/Cres2 == 0` | no |
| B12 | math | 93-94 | failure-side width > 0, success-side width > 0 | partial |
| B13 | math | 99-102 | convex-combination positivity | partial |

Every row marked "no" feeds `tautological_check`; rows marked "partial" are positivity facts that follow trivially from the declared symbol assumptions (`positive=True`) and so do not exercise any physical claim of the docstring.

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py:54-90` (definitions and the first nine assertions)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py:93-110` (delta-band definitions and assertions)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl:35-94` (mirror block)

**What's wrong:**
The script defines its symbolic objects so that every "expect_zero" assertion reduces to an algebraic identity that cannot fail. Concretely:

1. **Lines 57-61** define `Deltainf := Delta_0 + Delta_gap`, `Pres := 1 + Pres_gap`, and crucially `Cres2 := 1/Pres`. Then line 74's assertion `expect_zero("P_res - 1/C_res^2", Pres - 1/Cres2)` reduces to `Pres - 1/(1/Pres) = Pres - Pres = 0` by direct substitution of the definition. **The check is `A - A == 0`.**

2. **Lines 63-66** define `Wfail_match := Pe_req/Deltainf`, `Wsuff_match := Pe_req/Delta_0`, `Wfail_res := Pres*Wfail_match`, `Wsuff_res := Pres*Wsuff_match`. Then:
   - Line 75-78 asserts `Wsuff_match - Wfail_match == Pe_req*(Deltainf-Delta0)/(Delta0*Deltainf)`. Substituting the definitions gives `Pe_req/Delta0 - Pe_req/Deltainf = Pe_req*(Deltainf - Delta0)/(Delta0*Deltainf)` — a pure arithmetic identity for combining fractions, independent of any physics.
   - Line 81-84 asserts `Wfail_res == Pe_req/(Cres2*Deltainf)`. But `Wfail_res := Pres*Pe_req/Deltainf` and `1/Cres2 = Pres` by line 59, so this is `Pres*Pe_req/Deltainf == Pres*Pe_req/Deltainf`. Trivially true.
   - Lines 85-88 are the same identity with `Delta0` in place of `Deltainf`.

3. **Lines 93-94** define `delta_fail := Wfail_res - Wfail_match = (Pres-1)*Pe_req/Deltainf` and `delta_succ := Wsuff_res - Wsuff_match = (Pres-1)*Pe_req/Delta0`. Then lines 99-106 assert these equal `Pe_req*(1-Cres2)/(Cres2*Deltainf)` and `Pe_req*(1-Cres2)/(Cres2*Delta0)`. Using `Cres2 = 1/Pres`, `(1 - 1/Pres)/(1/Pres) = Pres - 1`, so the RHS reduces to `(Pres-1)*Pe_req/Deltainf` and `(Pres-1)*Pe_req/Delta0` — exactly the LHS by definition. **Lines 107-110 then directly assert `(Pres - 1) - (1 - Cres2)/Cres2 == 0`, which is the literal algebraic identity used to make 99-106 trivial.**

4. **Lines 89, 90, 79, 80, 111, 112** are positivity assertions on quantities that are products / sums / quotients of symbols declared `positive=True`. Each is immediate from the symbol assumption list. For example, `expect_positive("Delta_inf - Delta_0", Deltainf - Delta0)` evaluates to `Delta_gap`, which is positive by declaration on line 54.

5. **Lines 115-121** introduce convex combinations `W_failure_band = Wfail_match + u_fail * delta_fail` with `u_fail = v_fail/(1+v_fail) ∈ (0,1)` and assert (a) `W_failure_band - Wfail_match > 0` and (b) `Wfail_res - W_failure_band > 0`. With `u_fail ∈ (0,1)` and `delta_fail > 0`, (a) is `u_fail * delta_fail > 0` and (b) is `(1 - u_fail) * delta_fail > 0` — convex-combination trivia, not a physics claim.

The Mathematica script (lines 35-102) reproduces every one of these definitions and assertions in the same order with the same algebraic structure, so the same tautology critique applies block-for-block (`expectZero["P_res - 1/C_res^2", Pres - 1/Cres2]` on wl:57 is identical to A1, etc.).

**Why this matters:**
This unit is flagged as a *checkpoint*: it is supposed to consolidate Stage 049's matched-branch window with Stage 051/068's resonance penalty into the final three-zone verdict. A faithful checkpoint must verify that the consolidation step itself is mathematically correct — at minimum, that the matched-window edges and the penalty factor used here equal the values derived in their source stages. The script does the opposite: it introduces every relevant quantity (`Pe_req`, `Delta_0`, `Delta_inf`, `P_res`, `C_res^2`, `Wfail_match`, `Wsuff_match`, `Wfail_res`, `Wsuff_res`) as a free symbol or a definition built from free symbols, and then checks they satisfy identities that are built into those definitions. If Stage 066 derived `Wfail_match = 2*Pe_req/Deltainf` (off by a factor of 2) and Stage 068 derived `P_res = 1 + 2*Pres_gap` (off by a factor of 2), this script would still PASS every assertion, because it never reaches back to the source derivations. The docstring's promise ("the matched-branch Stage-49 window is exactly `[Pe_req / Delta_inf, Pe_req / Delta_0]`") is not exercised; only the post-Stage-49 ledger bookkeeping is.

**Required change:**
Replace the tautological identity block with at least two substantive consolidation checks. The minimum non-trivial set:

(a) **Anchor `Wfail_match` and `Wsuff_match` to a derivation, not a definition.** Stage 49's matched-window theorem is "on the matched branch, `J_eff = 0` at the wall iff `W_wall = Pe_req / Delta_phi_eff`, with the failure edge corresponding to the largest effective drop `Delta_inf` and the success edge to the smallest `Delta_0`." Encode this as an actual extremization. In SymPy, introduce a symbol `Delta_eff` ranging over `[Delta_0, Delta_inf]`, define `W_match(Delta_eff) = Pe_req/Delta_eff`, and `expect_zero` that the failure edge `W_match(Delta_inf)` is the minimum (i.e., `sp.diff(W_match, Delta_eff)` is negative and the boundary values give `[Pe_req/Delta_inf, Pe_req/Delta_0]`):

```python
Delta_eff = sp.symbols('Delta_eff', positive=True)
W_match = Pe_req / Delta_eff
expect_zero("W_match at Delta_eff=Delta_inf equals fail edge",
            W_match.subs(Delta_eff, Deltainf) - Wfail_match)
expect_zero("W_match at Delta_eff=Delta_0 equals success edge",
            W_match.subs(Delta_eff, Delta0) - Wsuff_match)
expect_positive("W_match decreasing in Delta_eff",
                -sp.diff(W_match, Delta_eff))
```

These checks tie `Wfail_match`/`Wsuff_match` to a generating function with a real monotonicity statement, so a wrong sign on `Delta_eff` or a wrong exponent on `Pe_req` would surface.

(b) **Anchor `P_res` to its Stage-68 definition, not declare it as a free symbol.** Stage 068's resonance penalty is `P_res = (some function of the resonance parameters)`. At a minimum, derive `P_res` from a residue / amplification computation in this script (or import the symbolic expression from the same form used in Stage 068) and `expect_zero` that the derived form simplifies to `1/Cres2`. If the Stage-68 form is unavailable, the script must at least show **two independent derivations** of `P_res`: e.g., (i) as the ratio of resonance-band edge to matched edge, and (ii) as `1 + Pres_gap` from a generic gap expansion, and assert they agree. The current `Pres = 1 + Pres_gap` followed by `Cres2 = 1/Pres` makes both definitions identical by fiat.

(c) **Remove or relabel the trivial-positivity assertions.** `Delta_inf - Delta_0 > 0` (line 79) and `P_res - 1 > 0` (line 90) follow from `Delta_gap > 0` and `Pres_gap > 0` declared on line 54. These add no information beyond the assumption list. Either drop them, or change them to conditional positivity checks that exercise an inequality not built into the assumptions (e.g., comparing the failure band's width to the success band's width: `delta_fail < delta_succ` iff `Delta_inf > Delta_0`, which is a non-tautological consequence of the assumption set).

(d) **Mirror these changes in the Mathematica script** with an independent derivation route (see F3).

**Verification:**
After the edit, the .txt outputs should show non-trivial residual lines such as `W_match at Delta_eff=Delta_inf equals fail edge = 0` AND a derivation of `P_res` printed before the line `P_res - 1/C_res^2 = 0` (e.g., `P_res derived from resonance ratio = (1 + Pres_gap)` or analogous). The number of `expect_zero` calls should not decrease without compensating derivation lines. Both scripts must still exit 0.

### F2 — hardcoded_result

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py:54-59`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl:35-43`

**What's wrong:**
The provenance comment on lines 20-23 states:

> The matched window `[Pe_req/Delta_inf, Pe_req/Delta_0]` is the Stage 066 theorem surface carried forward verbatim. `P_res = 1/C_res^2` is the Stage 068 resonance penalty; this script keeps those symbols unchanged and only assembles the final verdict from them.

But the script does not "carry the symbols forward" from any source — it introduces them as free `sp.symbols(...)` on lines 54-55, all `positive=True`. There is no import, no reference to a Stage-066/068 output value, no cross-engine cache, no symbolic substitution check. The expressions `Pe_req`, `Delta_0`, `Delta_gap`, `Pres_gap` are abstract placeholders with arbitrary positive values, and the matched-window form `Pe_req/Deltainf` is *defined here* (line 63) rather than derived elsewhere and checked here. The provenance comment is therefore unsupported by the script: nothing in the script ties these expressions to the upstream stages they claim to carry. The Mathematica script does the same — `Clear[...]` on line 35 followed by `$Assumptions = Element[...]` on line 36 introduces the symbols fresh.

If the upstream stages (066 and 068) compute concrete symbolic forms for `Wfail_match` and `P_res` (e.g., `Wfail_match = Pe_req/(Delta0 + Delta_gap)` from a specific BVP-derived effective gap, or `P_res = 1 + (some expansion coefficient)`), the checkpoint should *import* or *re-derive* those forms and assert they match. As written, the script is a self-contained algebra exercise with no upstream linkage.

**Why this matters:**
A checkpoint that "assembles the final verdict" from upstream results must, at minimum, name and reproduce those upstream results. The current script is structurally equivalent to writing a paper that says "let `X = Y` from Stage A, `Z = 1/X` from Stage B, then `XZ = 1`" — the conclusion is correct, but it carries no information about whether Stages A and B were correct. The audit's higher bar for checkpoints (per the prompt) is precisely to catch this pattern: a checkpoint must verify the consolidation, not perform fresh definitions.

**Required change:**
Add explicit upstream-anchor assertions at the start of the script (after line 66, before line 74). At minimum:

(a) Define each upstream quantity in symbolic terms that match its source stage's derivation, not as a free symbol. For `Delta_inf`, instead of `Deltainf = Delta_0 + Delta_gap`, use the Stage-066-form gap if known (e.g., `Delta_inf = Delta_0 + ∫ K(ξ) dξ` or whatever symbolic expression Stage 066 produced). If the exact upstream form is genuinely opaque from this unit's vantage, the minimum acceptable substitute is a **symbolic carry-forward assertion**: `expect_zero("matched fail edge matches Stage 066", Wfail_match - <Stage 066 form imported as a symbolic expression>)`, where `<Stage 066 form>` is constructed from primitive symbols and any reduction Stage 066 performs.

(b) For `P_res`, do the same against Stage 068. If Stage 068 derives `P_res` from a resonance integral or amplification factor, write that integrand/factor symbolically here and `expect_zero` its evaluation equals `1 + Pres_gap` (or whatever the stage-68 closed form is). Currently `Pres = 1 + Pres_gap` is the definition AND the result.

(c) If neither (a) nor (b) is possible without reading upstream scripts, the directive should add a comment block at the top of the script identifying the exact lines of the Stage-066 and Stage-068 scripts that produce `Wfail_match` and `P_res`, AND add a sentence in the docstring that this unit's assertions are *conditional* on those upstream derivations being correct (i.e., the script does not independently verify them). This is a weaker remediation but at least makes the limitation explicit.

The fully-correct remediation is option (a)+(b); the documentation-only remediation (c) is a fallback if it can be demonstrated that the upstream symbolic forms cannot be expressed in this unit's symbol set without expanding scope.

**Verification:**
After the edit, the .txt outputs must show either (i) new lines anchoring `Wfail_match` and `P_res` to upstream-stage symbolic forms (e.g., `Wfail_match matches Stage 066 form = 0`), or (ii) explicit docstring text marking the assertions as conditional on Stages 066 and 068. The script must still exit 0.

### F3 — mathematica_transliteration

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl:35-102` (entire body)

**What's wrong:**
The Mathematica script is a line-by-line transliteration of the SymPy script with no independent derivation. Every variable name, every definition, every assertion appears in the same order with the same algebraic structure. Three corresponding excerpts:

1. **Identical free-symbol declaration and the same `Cres2 := 1/Pres` definition.**
   - SymPy line 54-59: `Pe_req, Delta0, Delta_gap = sp.symbols("Pe_req Delta_0 Delta_gap", positive=True, real=True); Pres_gap, v_fail, v_succ = sp.symbols("Pres_gap v_fail v_succ", positive=True, real=True); Deltainf = sp.simplify(Delta0 + Delta_gap); Pres = sp.simplify(1 + Pres_gap); Cres2 = sp.simplify(1 / Pres)`.
   - Mathematica line 35-43: `Clear[PeReq, Delta0, DeltaGap, PresGap, Deltainf, Pres, Cres2, uFail, uSucc, vFail, vSucc]; $Assumptions = Element[{PeReq, Delta0, DeltaGap, PresGap, vFail, vSucc}, Reals] && PeReq > 0 && Delta0 > 0 && DeltaGap > 0 && PresGap > 0 && vFail > 0 && vSucc > 0; Deltainf = FullSimplify[Delta0 + DeltaGap, ...]; Pres = FullSimplify[1 + PresGap, ...]; Cres2 = FullSimplify[1/Pres, ...]`.

   Same symbol names (PEP-8 → camelCase rewrite), same definitions in the same order, same `Cres2 := 1/Pres` rather than independent derivation.

2. **Identical threshold definitions in the same order.**
   - SymPy line 63-66: `Wfail_match = sp.simplify(Pe_req / Deltainf); Wsuff_match = sp.simplify(Pe_req / Delta0); Wfail_res = sp.simplify(Pres * Wfail_match); Wsuff_res = sp.simplify(Pres * Wsuff_match)`.
   - Mathematica line 47-50: `WfailMatch = FullSimplify[PeReq/Deltainf, ...]; WsuffMatch = FullSimplify[PeReq/Delta0, ...]; WfailRes = FullSimplify[Pres WfailMatch, ...]; WsuffRes = FullSimplify[Pres WsuffMatch, ...]`.

   Same four definitions, same order, same algebraic dependence (Wfail_res = Pres * Wfail_match).

3. **Identical assertion sequence.**
   - SymPy lines 74-90 issue (in order): `Pres - 1/Cres2 == 0`, matched window width identity, `Deltainf - Delta0 > 0`, `Wsuff_match - Wfail_match > 0`, resonance fail threshold identity, resonance success threshold identity, `1 - Cres2 > 0`, `Pres - 1 > 0`.
   - Mathematica lines 57-73 issue **the same 8 assertions in the same order with the same names**: `expectZero["P_res - 1/C_res^2", ...]`, matched window width, `Delta_inf - Delta_0`, matched success minus fail, resonance fail form, resonance success form, `1 - C_res^2`, `P_res - 1`.

The same correspondence holds for lines 93-112 (sympy) ↔ lines 75-102 (wl): delta_fail / delta_succ definitions, three identity assertions, two positivity assertions, four convex-combination positivity assertions in identical order with identical names.

**Why this matters:**
The two-engine policy requires independent derivation paths so that a mistake or simplifier artifact in one engine cannot survive in the other by mere transcription. Here, both scripts share the same fundamental design flaw (F1: tautological definitions; F2: no upstream anchor), and the Mathematica engine cannot catch it because it makes the same definitions in the same order. If `Cres2 := 1/Pres` is the wrong relationship (e.g., it should be `Cres2 := 1/Pres^2`), both scripts will agree on the wrong answer and pass. A truly independent Mathematica derivation would, for instance, start from a different parameterization (e.g., directly use `Cres2` as the free parameter and derive `Pres = 1/Cres2`, then verify it equals `1 + PresGap` from a Stage-068-derived expression), or use `Reduce` / `Solve` to extract the resonance penalty from a different algebraic relation. Neither happens.

**Required change:**
Restructure the Mathematica script to derive the consolidation differently. Concrete options (any **one** is sufficient, but option (i) is the cleanest):

(i) **Reverse the parameterization.** Take `Cres2` as the primitive symbol (positive, less than 1) and derive `Pres = 1/Cres2`. Then on the resonance side, instead of `WfailRes = Pres * WfailMatch`, write `WfailRes = WfailMatch / Cres2` and verify it equals `(1 + PresGap) * WfailMatch` for `PresGap = 1/Cres2 - 1`. This makes the SymPy script's `Cres2 := 1/Pres` definition a *theorem* in Mathematica rather than a definition, so a wrong identity would surface.

(ii) **Use `Solve` to extract the resonance penalty.** Instead of defining `Pres = 1 + PresGap` directly, write the resonance condition symbolically (e.g., from the Stage-068 amplification expression or band-edge ratio `WfailRes / WfailMatch`) and use `Solve` to extract `Pres` as a function of `PresGap`. Verify the extracted form equals `1 + PresGap`. This routes the assertion through a different algebraic operation than SymPy's direct substitution chain.

(iii) **At minimum, change at least 4 of the 17 mirrored assertions in (3) above** so they probe different algebraic relationships (e.g., assert `WfailRes * Cres2 == WfailMatch` rather than `WfailRes == WfailMatch / Cres2`; assert `delta_fail / delta_succ == Delta0/Deltainf` rather than the two separate width-form identities). The new assertions must not share the same trivial-by-construction structure as SymPy's.

In all three options, the directive must rename at least two intermediate variables to avoid the line-by-line correspondence pattern (e.g., rename `WfailRes` to something derived from a primitive resonance ratio, and reorder so the assertion ordering differs from SymPy's).

**Verification:**
After the edit, the Mathematica .txt should show a visibly different derivation chain from the SymPy output: at least one `Solve[..., Pres]` or `Reduce` result printed, OR a different intermediate variable for the resonance penalty, OR a different ordering of the assertions. Both scripts must still exit 0 and their final residuals (on assertions that overlap in claim) must still agree.

## Independent-derivation check (Mathematica)

The Mathematica script is **not** an independent re-derivation. Its body (wl:35-102) mirrors the SymPy script (py:54-121) section-for-section with identical variable names (modulo camelCase), identical definitions, and identical assertion ordering. The three excerpts in F3 demonstrate the mapping. Because both engines share the tautological definition `Cres2 := 1/Pres` (sympy:59, wl:43), the Mathematica engine provides no cross-check on the consolidation step — a wrong relationship between `Pres` and `Cres2` would pass in both.

## Engine cross-check

Both engines produce zero residuals on every assertion (sympy .txt lines 17-35; wl .txt lines 17-52). The symbolic outputs agree as expected. However, because both scripts assert the same identities in the same order on the same trivially-related definitions, this agreement does not constitute independent verification — it just confirms two engines agree on basic algebra. See F3.

## Verdict justification

The verdict is `findings`: three high-severity findings affecting a checkpoint stage. What holds up: the algebra itself is correct — every assertion reduces to zero, and the positivity facts genuinely hold given the symbol assumptions. The SymPy and Mathematica engines agree at the symbolic level. What does not hold up: as a checkpoint of Stage 049 / Stage 051 / Stage 068, the script verifies none of those upstream results. Every threshold (`Wfail_match`, `Wsuff_match`), every penalty (`P_res`, `C_res^2`), and every width (`delta_fail`, `delta_succ`) is defined locally as a free-symbol composition or a one-line algebraic formula, and the subsequent `expect_zero` checks confirm those definitions obey the algebraic identities built into them. The checkpoint stage requires a higher bar; under that bar, F1 (tautological definitions throughout), F2 (no upstream anchor for the claimed carry-forward), and F3 (Mathematica transliteration) are all genuine. Attacks the script **did** survive: the positivity declarations are internally consistent (no contradictory assumptions); the convex-combination band points (W_failure_band, W_success_band) really do lie strictly between the matched and resonance edges given `u_fail, u_succ ∈ (0,1)`; and the algebraic identity `Pres - 1 = (1 - Cres2)/Cres2` is correctly stated.

## Self-test notes

(1) **Variable independence**: For the directive's proposed `sp.diff(W_match, Delta_eff)` (F1, item (a)), `W_match = Pe_req/Delta_eff` actually depends on `Delta_eff`, so the derivative is `-Pe_req/Delta_eff^2` (non-zero, negative — matches the proposed `expect_positive(-sp.diff(...))`). No identically-zero-derivative trap.
(2) **Symmetry/parity**: No integrals over symmetric domains are introduced by the directive — only algebraic identities and a single derivative. No parity trap.
(3) **Trivial-case pre-check**: For F1's `W_match.subs(Delta_eff, Deltainf) - Wfail_match`, both sides equal `Pe_req/Deltainf` by definition; the assertion is meant to *anchor* `Wfail_match` to the generating function, so the residual is correctly zero — but the new assertion structure (deriving from a parameterized form) catches errors that the current `expect_zero(Pres - 1/Cres2, ...)` does not, because the generating function `W_match(Delta_eff)` is genuinely a function of `Delta_eff` and an erroneous reuse of `Delta_0` in place of `Deltainf` would surface. For F2's symbolic carry-forward assertion, substituting concrete values (e.g., `Pe_req=1, Delta_0=1, Delta_gap=1` giving `Wfail_match = 1/2`) gives a literal nonzero value, so the proposed comparisons are non-trivial.
(4) **Path specifications**: This audit raises no `missing_verification_script` findings; both target files exist at the paths confirmed in the prompt (sympy in `scripts/`, mathematica in `mathematica/`).
