---
unit_id: 152
batch: IV.6
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: insufficient
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage152_family1_actual_correction.md]
  paper_appendix: present
---

# Audit unit 152 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_152.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage152_family1_actual_correction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only inclusion line at L1338; no separate prose row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage152_family1_actual_correction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage152_family1_actual_correction_mathematica_audit.txt`

## What the paper claims

The stage card (purpose + boxed quote) states: "Full mouth profile broadens the source and shifts the point to (Pi_corr, T_m,corr)." It declares the stage tests three checks: that first-order profile deformations are centered before covariance formulas are used; that the rigidity kernel (not branch ambiguity) controls non-exponential corrections; and that the one-step nonlinear correction remains within the reduced mouth-layer regime. The notes file supplies the load-bearing numerical deliverables (boxed equations): `Cov_*(c,R_*) ≈ 0.0648069687666328`, `Cov_*(K_q,R_*) ≈ 0.0388718368650403`, `delta g_act ≈ -0.0648069687666328`, `delta S_act ≈ -0.0388718368650403`, `delta Pi_act ≈ 0.907084414842908`, `delta T_m,act ≈ 0.271653979462338`, `Pi_corr ≈ 2.41591392833607`, `T_m,corr ≈ 1.17313803363654`, the one-step Picard iterate `Sigma_1(x) = e^{-Pi_* x - R_*(x)} / Z_1`, its moments `g_1 ≈ 0.684423574065325`, `S_1 ≈ 0.616333130570251`, retuned `Pi_1 ≈ 2.53914847609768`, `T_m,1 ≈ 1.21036942084359`, and effective interpolation `lambda_eff^(Pi) ≈ 0.380487632771110`, `lambda_eff^(T) ≈ 0.378939241176339` (i.e. `lambda_eff ≈ 0.38`, broadening fraction ≈ 0.62).

## What the script claims to verify

The SymPy script (mpmath-numerical) constructs `R_*(x) = Sigma_m_*·(4 T_s(x) - T_q(x)) - Pi_* x` from hardcoded canonical inputs and integrates Sigma_*(t)-weighted covariances against `c(t)=cos(pi x/2)` and `K_q(t)=cosh((pi/2)(1-x))/cosh(pi/2)`. It then computes `delta g_act`, `delta S_act`, `delta Pi_act`, `delta T_m,act`, `Pi_corr`, `T_corr`, the one-step iterate `sigma_1` and its moments `g_1`, `S_1`, the retuned `Pi_1`, `T_1`, and the effective interpolation parameters `lambda_eff^(Pi)`, `lambda_eff^(T)`. Apart from forward-difference prints of `R'(0)` and printed values of all of these, the only `expect_close` assertion in the script is a single linearized-covariance consistency identity (line 118-120) that checks `⟨c·(1-(R-⟨R⟩))⟩ - g_* ≈ -Cov(c,R)`. The Mathematica script independently derives all upstream constants from `gFormula[p]`, `sFormula[p]`, and the `gMinus` Family-1 root, computes the same downstream quantities, then asserts five `expectApprox` calls against the notes-quoted targets (delta_Pi, delta_T, lambda_eff^Pi, lambda_eff^T, and the covariance consistency identity).

## Paper ↔ script cross-check

| Paper-side deliverable | SymPy check | Mathematica check | Status |
|---|---|---|---|
| `Cov_*(c,R_*) ≈ 0.0648...` | printed only | printed only | partial (Mathematica `delta Pi_act scale` indirectly anchors via `deltaPi = covCR/gPrimeStar`) |
| `Cov_*(K_q,R_*) ≈ 0.0389...` | printed only | printed only | partial (anchored indirectly through `delta Tm_act scale`) |
| `delta g_act ≈ -0.0648...` | printed only | printed only | partial (anchored via covariance identity assertion) |
| `delta S_act ≈ -0.0389...` | printed only | printed only | partial (anchored indirectly) |
| `delta Pi_act ≈ 0.907084...` | printed only | `expectApprox` vs `0.907084414842908` (line 126) | match (Mathematica only) |
| `delta T_m,act ≈ 0.271654...` | printed only | `expectApprox` vs `0.271653979462338` (line 127) | match (Mathematica only) |
| `Pi_corr ≈ 2.4159...` | printed only | printed only | partial |
| `T_m,corr ≈ 1.1731...` | printed only | printed only | partial |
| `Sigma_1` definition | implemented (lines 93-98) | implemented (line 98) | match |
| `g_1 ≈ 0.6844...`, `S_1 ≈ 0.6163...` | printed only | printed only | partial |
| `Pi_1 ≈ 2.5391...`, `T_m,1 ≈ 1.2104...` | printed only | printed only | partial |
| `lambda_eff^(Pi) ≈ 0.3805...` | printed only | `expectApprox` vs `0.380487632771110` (line 128) | match (Mathematica only) |
| `lambda_eff^(T) ≈ 0.3789...` | printed only | `expectApprox` vs `0.378939241176339` (line 129) | match (Mathematica only) |

`paper_alignment: aligned` — when the assertions are present they match the paper. Numbers in printed transcripts match notes to all displayed digits. The defect is coverage-on-the-sympy-side, not paper↔script disagreement.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 118-120 | `expect_close(linearized cov identity, delta_g, tol=1e-9)` | delta g_act identity `delta_g = -Cov(c,R)` | partial (re-derives Cov(c,R) from same R_star/Sigma_star and confirms the linearized formula = covariance moment; consistency check, not a paper-value anchor) |
| A2 | mathematica | 125 | `expectApprox[linearized cov, deltaG, 10^-9]` | mirrors A1 | partial (same nature as A1) |
| A3 | mathematica | 126 | `expectApprox[deltaPi, 0.907084414842908, 10^-6]` | delta Pi_act ≈ 0.9070844... | yes |
| A4 | mathematica | 127 | `expectApprox[deltaT, 0.271653979462338, 10^-6]` | delta T_m,act ≈ 0.2716540... | yes |
| A5 | mathematica | 128 | `expectApprox[lambdaEffPi, 0.380487632771110, 10^-6]` | lambda_eff^(Pi) ≈ 0.3805... | yes |
| A6 | mathematica | 129 | `expectApprox[lambdaEffT, 0.378939241176339, 10^-6]` | lambda_eff^(T) ≈ 0.3789... | yes |

The SymPy script's one assertion (A1) is a self-consistent rearrangement of the covariance moment definition; it does not anchor against any paper-quoted numeric target. All numerical scale values (delta Pi_act, delta T, Pi_corr, T_corr, g_1, S_1, Pi_1, T_1, lambda_eff) are produced only as `print` statements.

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py:80-115`

**What's wrong:**
The SymPy script computes every load-bearing numerical deliverable that the notes file boxes (delta Pi_act, delta T_m,act, Pi_corr, T_corr, g_1, S_1, Pi_1, T_1, lambda_eff^(Pi), lambda_eff^(T)) but only `print`s them. The script's only `expect_close` assertion (lines 118-120) is a self-consistency identity that algebraically rewrites `Cov(c,R) = ⟨c·R⟩ - ⟨c⟩⟨R⟩`; that assertion would pass whether or not the upstream constants `Pi_star`, `Sigma_m_star`, `gprime_star`, `AT`, `BT` produced the right physical values. The Mathematica counterpart correctly asserts five anchored scale checks (`expectApprox` at lines 126-129 against the notes-quoted constants `0.907084414842908`, `0.271653979462338`, `0.380487632771110`, `0.378939241176339`). Quoting the SymPy script:

```
deltaPi = Cov_cR / gprime_star
deltaT = AT * delta_g + BT * delta_S
...
print("delta Pi_act      =", deltaPi)
print("delta Tm_act      =", deltaT)
...
print("lambda_eff^(Pi)   =", lam_Pi)
print("lambda_eff^(T)    =", lam_T)
```

Notes (`moving_throat_pde_stage152_family1_actual_correction.md`, boxed):
```
\delta\Pi_{\rm act} \approx 0.907084414842908,
\delta\widehat T_{m,{\rm act}} \approx 0.271653979462338.
```
and
```
\lambda_{\rm eff}^{(\Pi)}\approx 0.380487632771110,
\lambda_{\rm eff}^{(T)}\approx 0.378939241176339.
```

No `expect_close` in the SymPy script binds these printed values to the paper-quoted targets, so any drift in inputs (`Pi_star`, `Sigma_m_star`, `AT`, `BT`, or the hardcoded baselines `1.69941496131430`, `2.08240814741023`) could silently change the prints while the script still exits 0.

**Why this matters:**
The single-engine policy depends on each engine producing a substantive PASS/FAIL signal on the paper-side claims. As written, the SymPy script always exits 0 regardless of the deliverable values; a regression that altered the canonical inputs would not be caught by the SymPy audit even though Mathematica would catch it. This violates the explicit charter that scripts independently verify the paper's claims.

**Required change:**
Add `expect_close` assertions in the SymPy script anchoring the computed deliverables to the notes-quoted targets, with the same tolerance pattern the Mathematica script uses (`1e-6` for scale checks). Specifically, after line 115 in `scripts/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py`, before the existing `expect_close` on the covariance consistency identity, add:

```python
expect_close("delta Pi_act scale", float(deltaPi), 0.907084414842908, tol=1e-6)
expect_close("delta Tm_act scale", float(deltaT), 0.271653979462338, tol=1e-6)
expect_close("lambda_eff^(Pi) scale", float(lam_Pi), 0.380487632771110, tol=1e-6)
expect_close("lambda_eff^(T) scale", float(lam_T), 0.378939241176339, tol=1e-6)
```

These are exactly the four anchored checks the Mathematica script already runs (`mathematica/moving_throat_pde_stage152_family1_actual_correction_mathematica_audit.wl:126-129`). Do not modify the rest of the script; do not change the existing covariance-consistency assertion.

**Verification:**
After Codex applies, the SymPy script must contain four new `expect_close` calls referencing the constants `0.907084414842908`, `0.271653979462338`, `0.380487632771110`, `0.378939241176339`. Re-running `redteam exec-sympy 152` should still exit 0, and the saved output transcript should add lines reporting nonzero `diff` values for each of the four checks below the tolerance.

## Independent-derivation check (Mathematica)

The Mathematica script is NOT a transliteration of the SymPy script. SymPy hardcodes `Pi_star`, `Sigma_m_star`, `g_star`, `S_star`, `gprime_star`, `AT`, `BT`, `Tstar` as mpmath literals and then evaluates downstream integrals. Mathematica re-derives those same canonical inputs from first principles:

```
gFormula[p_] := 2*p*(2*p*Exp[p] + Pi)/((4*p^2 + Pi^2)*(Exp[p] - 1));
sFormula[p_] := p*(kap*Tanh[kap] + p*(Exp[-p]*Sech[kap] - 1))/((1 - Exp[-p])*(kap^2 - p^2));
rF1 = Sqrt[(12*(37/20)^2)/Pi^2 - 1];
gMinus = rF1 - Sqrt[1 + rF1^2]/2;
pStar = p /. FindRoot[gFormula[p] == gMinus, {p, 1.5}, ...];
...
aT = N[-(9/(40*tStar))*(1/(gPrimeStar*(1 - sStar/4)) + pStar*sPrimeStar/(4*gPrimeStar*(1 - sStar/4)^2)), 40];
bT = N[(9/(40*tStar))*pStar/(4*(1 - sStar/4)^2), 40];
```

The lambda_eff baseline is also derived independently: SymPy uses hardcoded `1.69941496131430` and `2.08240814741023` from "Stage 148 affine laws", while Mathematica computes `gU = 2/Pi`, `sU = 2 Tanh[Pi/2]/Pi`, `gD = Pi/4`, `sD = (1+Sinh[Pi/2])/(2 Cosh[Pi/2])` and then `(deltaPi - deltaPiU)/(deltaPiD - deltaPiU)`. The two engines therefore exercise structurally independent paths to the same numerical answers, which is what the second-engine policy requires.

## Engine cross-check

Both engine outputs agree to >=14 decimal places on every computed deliverable:

| Quantity | SymPy | Mathematica |
|---|---|---|
| Cov_*(c,R_*) | 0.0648069687666328 | 0.0648069687666327 |
| Cov_*(K_q,R_*) | 0.0388718368650403 | 0.0388718368650402 |
| delta Pi_act | 0.907084414842908 | 0.9070844148429054 |
| delta T_m,act | 0.271653979462338 | 0.2716539794623376 |
| Pi_corr | 2.41591392833607 | 2.4159139283360610 |
| T_corr | 1.17313803363654 | 1.1731380336365416 |
| g_1 | 0.684423574065325 | 0.684423574065324 |
| S_1 | 0.616333130570251 | 0.616333130570250 |
| Pi_1 | 2.53914847609768 | 2.539148476097673 |
| T_1 | 1.21036942084359 | 1.210369420843591 |
| lambda_eff^(Pi) | 0.380487632771111 | 0.380487632771111 |
| lambda_eff^(T) | 0.378939241176339 | 0.378939241176337 |

Engines agree. No `engine_disagreement` finding.

## Verdict justification

The paper-side claims are entirely numerical: thirteen boxed scale values across the notes file. The Mathematica script verifies the four most load-bearing of them (delta Pi_act, delta T_m,act, lambda_eff^Pi, lambda_eff^T) via anchored `expectApprox` calls against the exact notes-quoted constants, computes the remaining ones consistently, and reaches them through an independent derivation path (not a port of the SymPy script). The SymPy script computes the same quantities to the same precision but only `print`s them — the sole assertion is a covariance moment identity that would pass regardless of whether the canonical inputs were correct. That is `insufficient_verification` on the SymPy side. Both engines produce identical numerical results matching the paper's boxed targets to all displayed digits, so paper alignment is good and no `paper_misalignment`, `tautological_check`, `mathematica_transliteration`, `engine_disagreement`, or `hardcoded_result` finding applies. Output transcripts are both newer than their generating scripts, so no `stale_output`. Attacks tried that failed: (i) verifying the SymPy `lam_Pi` baseline literals `1.69941496131430` and `2.08240814741023` reproduce Mathematica's `(deltaPiU, deltaPiD)` chain — they do, given Mathematica's `gU=2/Pi`, `gD=Pi/4`, and `gprime_star`; (ii) checking R'(0) ≈ 0 tangent condition — confirmed by the forward-difference outputs (≈ -1.31e-6 at h=1e-6 and ≈ -1.31e-5 at h=1e-5, scaling linearly in h means R'(0)=0 and R''(0) ≈ -2.62); (iii) checking that `sigma_1` definition matches notes — does; (iv) checking parity / symmetry of the integrands on [0,1] — no parity issue since domain is asymmetric and weights are explicit. Single script-side finding: SymPy missing anchored scale assertions.

## Self-test notes

Checked traps: (1) Variable independence — no `sp.diff` or `D[expr,var]` issues; both scripts are numerical. (2) Symmetry/parity — integration domain is [0,1] (one-sided), so no odd/even cancellation traps. (3) Trivial-case pre-check — the four proposed `expect_close` calls compare floats from running SymPy code to literals from the notes; running mentally, `deltaPi` is computed as `Cov_cR / gprime_star` where `Cov_cR ≈ 0.0648069687666328` and `gprime_star = 0.0714453558083195`, giving `≈ 0.907084414842908` which matches the literal to `1e-6`. Same logic for `deltaT`, `lam_Pi`, `lam_T`. (4) Path specifications — Required change targets `scripts/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py`, correct location. (5) Paper round-trip — the new asserts use exactly the constants quoted in the notes; no new paper_misalignment introduced.
