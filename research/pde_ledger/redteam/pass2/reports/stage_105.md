---
unit_id: 105
batch: IV.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md]
  paper_appendix: present
---

# Audit unit 105 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_105.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows/eqs: `eq:app-part04-retarded-one-pole`, `eq:app-part04-sigmaQcan`, `eq:app-part04-Lambda-out-dtn`, `eq:app-part04-Yout-dtn`, `eq:app-part04-chiQ-equals-one`, `eq:app-part04-chiQ-general`; result-anchor MTDC-T8.2 line 1176; card checks 1244)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.txt`

## What the paper claims

Stage 105 (checkpoint, anchor MTDC-T8) fixes the last reduced 2.5PN scalar `chi_Q` of the canonical retarded grouped-`P2` one-pole module. The card states the result as `\widehat m_0^2 chi_Q N_Q = 1` together with the canonical condition `\chi_Q=1`, and quotes: *"Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\)."* The notes/appendix give the mechanism: the retarded module `\widehat Y_Q^{ret}(\omega)=3/4+(1/4)/(1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{can}\omega^5)` expands to `1+a^2\omega^2/(9c_s^2)+4a^4\omega^4/(81c_s^4)+i\chi_Q a^5\omega^5/(27c_s^5)`, with `\Omega_Q=3c_s/(2a)`, `\sigma_Q^{can}=9/(8\Omega_Q^5)=4a^5/(27c_s^5)`. Crucially, the **outgoing target** coefficient `i z^5/27` (i.e. `i a^5\omega^5/(27c_s^5)`) is asserted by the paper to come from the exact spherical Hankel DtN operator `\Lambda_2^{out}(z)=z\,d/dz\ln h_2^{(1)}(z)=-3+z^2/3+z^4/9+i z^5/9`, normalized as `\widehat Y_2^{out}=-3/\Lambda_2^{out}=1+z^2/9+4z^4/81+i z^5/27`. Matching the `O(\omega^5)` coefficients fixes `\chi_Q=1`. The card's third explicit Check requires: *"Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion."* A second deliverable is the general deformed result `\chi_Q=\xi_Q` (the only remaining 2.5PN obstruction is the deviation of the actual branch from the canonical `\xi_Q=1`). Deliverables: (D1) `\sigma_Q^{can}=4a^5/(27c_s^5)`; (D2) the retarded series coefficients (z²,z⁴,z⁵); (D3) `\chi_Q=1` from the outgoing DtN fingerprint match; (D4) `\chi_Q=\xi_Q` general.

## What the script claims to verify

Both engines: (1) confirm `\sigma_Q^{can}=(9/8)/\Omega_Q^5=4a^5/(27c_s^5)`; (2) series-expand the retarded module to `O(\omega^5)` and assert its z²,z⁴ real coefficients and its `i\chi_Q a^5/(27c_s^5)` imaginary z⁵ coefficient; (3) **solve / `Reduce` the equation `[retarded i\omega^5 coeff] = a^5/(27c_s^5)` for `\chi_Q`** and assert the solution is 1; (4) build the deformed operator `\Lambda_2^{def}=-3+z^2/3+z^4/9+i\xi_Q z^5/9`, normalize via `-3/\Lambda` (SymPy `series`) or polynomial inversion of `\Lambda\cdot Y=-3` (Mathematica), and assert its z⁵ coefficient is `i\xi_Q/27` (i.e. `\chi_Q=\xi_Q`). The bottom-line assertions are the `\chi_Q-1=0` check (step 3) and the deformed coefficient checks (step 4).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 `\sigma_Q^{can}=4a^5/(27c_s^5)` | `expect_zero("sigma_Q^can - 4 a^5/(27 c_s^5)", ...)` (py:37; wl:37) | match |
| D2 retarded series coeffs (z²,z⁴,z⁵) | `expect_zero` on coeff(2/4/5) (py:45–47; wl:48–50) | match |
| D3 `\chi_Q=1` **from outgoing DtN fingerprint** | solve/Reduce retarded ω⁵ coeff `= a^5/(27c_s^5)` → `\chi_Q-1=0` (py:49–55; wl:52–59) | **partial** — pins `\chi_Q=1` but matches against a HARDCODED literal `a^5/(27c_s^5)`, NOT the Hankel DtN value; the DtN fingerprint identity (`z d/dz ln h_2^{(1)}`) is never exercised (card Check #3 not performed) |
| D4 `\chi_Q=\xi_Q` general | deformed branch z⁵ coeff `= i\xi_Q/27` (py:64–66; wl:75–78) | match (and Mathematica route is genuinely independent) |

`paper_alignment: aligned` — every emitted value matches the card/notes/appendix; no `paper_misalignment`. The defect is in *anchoring strength*, not in value correctness.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 37 | `simplify(sigma_can - 4a^5/(27c_s^5)) == 0` | D1 | yes (derives `(9/8)/Ω⁵` then compares) |
| A2 | sympy | 45 | `coeff(ω,2) - a²/(9c_s²) == 0` | D2 | yes |
| A3 | sympy | 46 | `coeff(ω,4) - 4a⁴/(81c_s⁴) == 0` | D2 | yes |
| A4 | sympy | 47 | `coeff(ω,5)/i - χ_Q a⁵/(27c_s⁵) == 0` | D2 | yes (real exercise of the series) |
| A5 | sympy | 49–55 | `solve(coeff(ω,5)/i == a⁵/(27c_s⁵), χ_Q)[0]==1` | D3 | **no** — RHS hardcoded; reduces to `χ_Q·K=K` |
| A6 | sympy | 64 | deformed `coeff(z,2)-1/9==0` | D4 | yes |
| A7 | sympy | 65 | deformed `coeff(z,4)-4/81==0` | D4 | yes |
| A8 | sympy | 66 | deformed `coeff(z,5)/i - ξ_Q/27==0` | D4 | yes |
| B1 | math | 37 | `expectZero[sigmaQ - 4a⁵/(27c_s⁵)]` | D1 | yes |
| B2 | math | 48–50 | real ω²,ω⁴ + imag ω⁵ coeff checks | D2 | yes |
| B3 | math | 52–59 | `Reduce[imagPart5 - a⁵/(27c_s⁵)==0,χ_Q]`; `expectZero[χ_Q-1]` | D3 | **no** — same hardcoded-RHS structure as A5 |
| B4 | math | 75–78 | deformed b0,b2,b4,b5 via polynomial inverse | D4 | yes (independent route) |

A5/B3 are the load-bearing `\chi_Q=1` pin and are the rows that fail the "Anchored?" test.

## Findings

### F1 — insufficient_verification

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:49-55`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:52-59`

**What's wrong:**

105 owns the canonical pin `\chi_Q=1`. The paper grounds this pin in a genuine identity: the *outgoing target* coefficient `i a^5\omega^5/(27c_s^5)` is the normalized exact-Hankel DtN fingerprint `\widehat Y_2^{out}=-3/\Lambda_2^{out}`, `\Lambda_2^{out}(z)=z\,d/dz\ln h_2^{(1)}(z)` (appendix `eq:app-part04-Lambda-out-dtn`/`eq:app-part04-Yout-dtn`), and the card's Check #3 explicitly is "Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion."

Neither script computes that fingerprint. The pin is produced by:

SymPy (py:49–52):
```python
chi_sol = sp.solve(sp.Eq(Yret_series.coeff(omega, 5)/sp.I, a**5/(27*c_s**5)), chi_Q)[0]
```
Mathematica (wl:52–53):
```
chiBranch = FullSimplify[Reduce[imagPart5 - aThroat^5/(27*cSound^5) == 0, chiQ, Reals], ...]
```

The left side, `Yret_series.coeff(omega,5)/i` (= `imagPart5`), is **by construction** exactly `chi_Q * a^5/(27c_s^5)`: it was generated from `sigma_can`, which line 33/36 *defined* as `(9/8)/Omega_Q^5 = 4a^5/(27c_s^5)`. The right side is the same literal `a^5/(27c_s^5)`, typed in directly. So the solved equation is identically `chi_Q · K = K` with `K = a^5/(27c_s^5)`, whose root is `chi_Q = 1` for **any** nonzero `K` — the equation cannot detect whether `K` is the true outgoing DtN coefficient. The DtN operator `z d/dz ln h_2^{(1)}(z)` is never evaluated, never series-expanded, never normalized; the "fingerprint match" the card promises is replaced by self-matching a hardcoded copy of the answer. Re-derived independently: yes, `chi_Q=1` is the *correct* value (the retarded ω⁵ coefficient is `i·chi_Q·a^5/(27c_s^5)` and the genuine normalized Hankel value is `i·a^5/(27c_s^5)`), but the script's assertion does not *test* that the two agree — it asserts agreement by typing the same number on both sides.

Note the asymmetry: the deformed half (A6–A8 / B4) *does* derive its target coefficients from the operator (`-3/\Lambda` in SymPy; polynomial inversion of `\Lambda\cdot Y=-3` in Mathematica), so `\chi_Q=\xi_Q` is genuinely exercised. Only the canonical `\chi_Q=1` pin — the checkpoint's headline result — is left self-referential.

**Why this matters:**

This is a checkpoint, and 105 is the unit that establishes `\chi_Q=1` as the DtN l=2 fingerprint pin carried forward by stages 097/100/106 and the 2.5PN/4PN bridge (107–113). A reader auditing the chain would expect the pin to be forced by the Hankel DtN identity, per Check #3. As written, a transcription error or sign error in the *true* DtN coefficient could not be caught here: the script would still "pass" because both sides carry the same literal. The pin is correct, but it is asserted, not verified.

**Required change:**

Make the canonical match exercise the actual outgoing DtN fingerprint rather than a hardcoded literal, in BOTH engines. Construct the exact `l=2` outgoing DtN operator from the spherical Hankel function and take its series, then normalize and read off the imaginary z⁵ coefficient as the *target*, instead of typing `a^5/(27c_s^5)`:

- Define the outgoing radial solution `h2 = SphericalHankelH1[2, z]` (Mathematica) / `sympy.functions.special.bessel.hankel1`-route or the closed form `h_2^{(1)}(z) = -((3/z^2 - 1) + i*3/z)/z * exp(i z)` (SymPy), then
- `Lam_out = z * D[Log[h2], z]` and `Lam_out_series = Series[Lam_out, {z,0,5}]` (Mathematica) / `sp.series(z*sp.diff(sp.log(h2), z), z, 0, 6)` (SymPy);
- assert `Lam_out_series == -3 + z^2/3 + z^4/9 + i z^5/9` (this is the genuine fingerprint exercise — Check #3);
- `Y2_out = sp.series(-3/Lam_out, z, 0, 6)`; read its imaginary z⁵ coefficient `T := Y2_out.coeff(z,5)/i` (should equal `1/27`), translate to ω via `z = a ω / c_s` to get `T_omega = a^5/(27 c_s^5)` **as a derived quantity**, and
- solve `Yret_series.coeff(omega,5)/i == T_omega` for `chi_Q`, asserting `chi_Q - 1 == 0`.

Then the `chi_Q=1` pin is forced by the genuine Hankel DtN coefficient. (If 104 `outgoing_dtn_fingerprint` already verifies the Hankel `→ -3+z^2/3+z^4/9+i z^5/9` expansion and 105 is intended to consume that as a carry-forward, the alternative acceptable fix is to make the match symbolic against a single named DtN-coefficient symbol `T_dtn`, derive `chi_Q = T_dtn / (a^5/(27c_s^5)) · (a^5/(27c_s^5))`-style so the literal appears once, and add an in-script `expect_zero` tying `T_dtn` to the 104-verified value — but the self-test prefers deriving the Hankel coefficient locally so the pin is not self-referential.)

**Verification:**

After the fix, the script must contain a new derived `Lam_out`/`Y2_out` (or named DtN-coefficient symbol) whose imaginary z⁵ coefficient is computed, and the `chi_Q` solve must reference that derived quantity rather than a typed `a^5/(27c_s^5)` literal on the RHS. New output line should show the Hankel-derived `Lambda_2^out` expansion equalling `-3+z^2/3+z^4/9+i z^5/9` and `chi_Q - 1 = 0` still passing. Both engines exit 0.

### F2 — mathematica_transliteration

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:33-59`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:30-55`

**What's wrong:**

The **retarded half** of the `.wl` is a close port of the `.py`. Corresponding lines:

| step | SymPy | Mathematica |
|---|---|---|
| pole/sigma | `Omega_Q=3c_s/(2a)`, `sigma_can=(9/8)/Omega_Q**5` (py:32–33) | `omegaPole=3cSound/(2aThroat)`, `sigmaQ=(9/8)/omegaPole^5` (wl:33–34) |
| series black-box | `sp.series(Yret, omega, 0, 6)` (py:40) | `Series[retFactored, {omega,0,5}]` (wl:42) |
| coeff extraction + same target | `coeff(ω,5)/i` vs hardcoded `a^5/(27c_s^5)` (py:47,50) | `imagPart5` vs hardcoded `aThroat^5/(27cSound^5)` (wl:50,53) |
| assertion order | sigma → ω²,ω⁴,ω⁵ → χ_Q→1 (py:37,45–55) | sigma → ω²,ω⁴,ω⁵ → χ_Q→1 (wl:37,48–59) |

Same variable choreography, same series black-box, the SAME hardcoded RHS literal, identical assertion ordering. The only divergences in this half are presentational: Mathematica pre-combines `3/4+1/4/(...)` into a single fraction (wl:40–41) and uses `ComplexExpand`/`Re`/`Im` to split the projection. That is line-by-line porting, not an independent re-derivation, for the retarded section.

**Mitigation (why low, and why I do NOT escalate to a full transliteration verdict):** the **deformed half is genuinely independent.** SymPy uses `series(-3/Lam_def)` (py:60); Mathematica instead sets a trial polynomial `b0+...+b5 z^5`, forms `outgoingOp*trial`, expands, and solves the linear system `Lambda·Y+3=0` coefficient-by-coefficient via `Solve` (wl:67–73). That is a structurally different algebraic route (polynomial inversion vs direct series of a reciprocal) and is exactly the kind of second-engine independence the policy wants. So the `.wl` is **partial**: independent on D4, transliterated on D1–D3.

**Why this matters:**

For a checkpoint owning the canonical pin, the retarded-half port means both engines share the same self-referential weakness in F1 with no cross-check value — Mathematica cannot catch a transcription error in the SymPy retarded route because it re-types the same literal and the same module. Fixing F1 (deriving the DtN coefficient from `h_2^{(1)}` in each engine) is the natural opportunity to also de-transliterate: have one engine reach the outgoing coefficient by `Series[z D[Log[SphericalHankelH1[2,z]],z]]` and the other by the recurrence/explicit-Hankel-rational form, so the retarded-half match no longer mirrors.

**Required change:**

Subsumed by F1's fix: when the canonical match is re-grounded in the Hankel DtN, route the two engines through different constructions of `\Lambda_2^{out}` (e.g. Mathematica via `SphericalHankelH1[2,z]` + `Series`; SymPy via the explicit rational closed form of `h_2^{(1)}` or vice-versa) so the retarded/canonical half is no longer a line-by-line port. No separate edit beyond F1 is required if F1 is implemented with deliberately distinct DtN constructions per engine; this finding documents the requirement.

**Verification:**

After F1's fix, confirm the `.wl` retarded/canonical section does not mirror the `.py` step-for-step on the same hardcoded literal — at minimum the two engines construct `\Lambda_2^{out}` (or the target z⁵ coefficient) by visibly different operations, and the RHS literal `a^5/(27c_s^5)` no longer appears typed on both sides.

## Independent-derivation check (Mathematica)

Partial. **Deformed branch (independent):** SymPy `Y_def = series(-3/Lam_def, z, 0, 6)` (py:60) vs Mathematica polynomial inversion — trial `b0+b1 z+...+b5 z^5`, `prodTrunc=Series[outgoingOp*trialBranch]`, `Solve[Thread[(prodCoeffs+{3,0,..})==0]]` (wl:66–73). Genuinely different algebra reaching the same `1+z^2/9+4z^4/81+i\xi_Q z^5/27`. **Retarded/canonical branch (transliteration):** see F2 — same pole/sigma construction, same `series`/`Series` black-box on the same module, same hardcoded RHS, same assertion order; Mathematica differs only in pre-combining the fraction and `ComplexExpand`-splitting. Net `.wl` independence verdict: **partial.**

## Engine cross-check

Outputs agree exactly. SymPy retarded series (txt:9–14): `1 + a²ω²/(9c_s²) + 4a⁴ω⁴/(81c_s⁴) + i·a⁵χ_Q·ω⁵/(27c_s⁵)`. Mathematica (txt:9–10): real projection `1 + aThroat²ω²/(9cSound²) + 4aThroat⁴ω⁴/(81cSound⁴)`, imaginary projection `aThroat⁵·chiQ·ω⁵/(27cSound⁵)`. Both yield `chi_Q = 1` (sympy txt:19; math `chiQ == 1` txt:17–18). Deformed branch identical: `1 + z²/9 + 4z⁴/81 + i·ξ_Q·z⁵/27` (sympy txt:21–25; math txt:21). All `expect_zero`/`expectZero` checks read 0/PASS. No `engine_disagreement`.

## Verdict justification

The math is correct — I re-derived both `\chi_Q=1` and `\chi_Q=\xi_Q` by hand from the script construction and confirmed the retarded ω⁵ coefficient is `i\chi_Q a^5/(27c_s^5)`, the polynomial-inverse deformed coefficients are `b0=1,b2=1/9,b4=4/81,b5=i\xi_Q/27`, and `\sigma_Q^{can}=(9/8)(2a/3c_s)^5=4a^5/(27c_s^5)`. Every emitted value reconciles with the card, notes, and appendix (paper_alignment aligned). The verdict is **findings**, not clean, because this is a checkpoint owning the canonical `\chi_Q=1` pin and that pin is **not forced by a genuine identity**: the load-bearing match (A5/B3) solves `[retarded ω⁵ coeff] = a^5/(27c_s^5)` where the LHS is by construction `\chi_Q·a^5/(27c_s^5)` and the RHS is a hardcoded copy of the same literal, so the equation reduces to `\chi_Q·K=K → \chi_Q=1` for any K and never evaluates the outgoing DtN fingerprint `z\,d/dz\ln h_2^{(1)}(z)` that the card's Check #3 demands (F1, insufficient_verification, medium). The retarded half of the `.wl` additionally re-types the `.py` choreography and the same literal (F2, mathematica_transliteration, low; the deformed half is genuinely independent, so the `.wl` is partial). Both findings are script-side and fixable by deriving the Hankel DtN coefficient locally (and routing the two engines through distinct constructions). No `paper_misalignment`, no stop-cold.

**Checkpoint higher-bar verdict: not-cleared.** chi_Q=1 derivation result: the value is *correct* (independently re-derived), but the in-script assertion does not *exercise* the DtN-fingerprint identity — it self-matches a hardcoded copy of the target coefficient, so the canonical pin is asserted rather than forced. The higher bar (substantive, non-tautological, fingerprint actually exercised) is not met until F1 is fixed.

## Self-test notes

Checked: (1) Variable independence — no spurious `diff`; the retarded series genuinely depends on ω and χ_Q, deformed on z and ξ_Q. (2) The `solve`/`Reduce` for χ_Q is linear with a unique real root; `positive=True`/`chiQ>0` assumptions don't change it (no symbol_assumption_error). (3) Trivial-case: substituting the literal `K=a^5/(27c_s^5)` into A5/B3 confirms the equation is `χ_Q·K=K`, root 1 independent of K — this is precisely the insufficiency in F1. (4) Re-derived σ, the retarded series, and the polynomial-inverse coefficients by hand; all match outputs. (5) Paper round-trip: F1's prescribed fix derives the SAME coefficient `a^5/(27c_s^5)` from `h_2^{(1)}`, introducing no new constant and no new paper_misalignment.

## Value Reconciliation (pass-2 augmentation)

Outputs are fresh (sympy .txt 2026-05-29 10:40 > .py 2026-05-27 15:12; math .txt 2026-05-29 10:40 > .wl 2026-05-29 10:36); reconciliation based on script source + committed outputs (no execution).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `\Omega_Q = 3c_s/(2a)` | py:32 / wl:33; sympy txt:5, math txt:5 | notes:29,36; appendix `eq:app-part04-retarded-one-pole` (uses Ω_Q); card body (Inputs) | MATCH |
| `\sigma_Q^{can} = 4a^5/(27c_s^5)` | py:33,37 / wl:34,37; sympy txt:6, math txt:6 | notes:41 (boxed); appendix:275 `eq:app-part04-sigmaQcan` | MATCH |
| retarded `Y_Q^{ret}` z² coeff `a^2/(9c_s^2)` | py:45 / wl:48; sympy txt:10–14, math txt:9 | notes:50; appendix `eq:app-part04-Yout-dtn` (=z²/9 in z-vars) | MATCH |
| retarded z⁴ coeff `4a^4/(81c_s^4)` | py:46 / wl:49; sympy txt:11, math txt:9 | notes:51; appendix:321 | MATCH |
| retarded i z⁵ coeff `χ_Q a^5/(27c_s^5)` | py:47 / wl:50; sympy txt:11, math txt:10 | notes:53; appendix:321 (=i z⁵/27) | MATCH |
| `\chi_Q = 1` (canonical pin) | py:53–55 / wl:58–59; sympy txt:19, math txt:17–18 | notes:69 (boxed); card line 16, appendix:326 `eq:app-part04-chiQ-equals-one` | MATCH (value); see F1 re anchoring |
| deformed branch `Y_2^{def}=1+z²/9+4z⁴/81+iξ_Q z⁵/27` | py:60–66 / wl:73–78; sympy txt:21–25, math txt:21 | notes:87–92 | MATCH |
| `\chi_Q = \xi_Q` (general) | py:70 / wl (implicit via b5=iξ_Q/27); sympy txt:32 | notes:95 (boxed); appendix `eq:app-part04-chiQ-general` (general form) | MATCH |

INTERNAL scaffolding (no finding): `omegaPole^2`/`retFactored` intermediate fraction (wl:40–41); `imagProjection`/`realPart`/`retWindow` intermediate projections; trial polynomial `b0..b5`, `prodCoeffs`, `coeffSys`, `branchRules` (wl:67–73); `chi_sol`/`chiBranch`/`chiWitness` solve handles; all `expect_zero`/`expectZero` residual prints and PASS flags.

reconciliation: complete; 8 deliverable values checked, 0 misaligned. (The F1/F2 findings are anchoring/independence defects, not value mismatches — every emitted value agrees with the card and notes.)
