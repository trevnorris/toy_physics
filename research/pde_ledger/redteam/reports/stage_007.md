---
unit_id: 007
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 007 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.txt`
- mathematica output: `(missing)`

## What the script claims to verify

The script contrasts two derivations of the effective Maxwell coupling in a zero-mode reduction with a Gaussian transverse profile `Z(w) = exp(-w^2/lambda^2)`. The reduction-first answer is `mu0_eff^(red) = mu0/Z_int` (with `Z_int = sqrt(pi)*lambda`). The projection-first answer, under the assumptions `F^{mu nu}(x,w) = f^{mu nu}(x)`, `J^nu = j^nu(x) S(w)`, and `int W dw = 1`, becomes `mu0_eff^(proj) = mu0 * I_WS / I_WZ` with `I_WZ = int W Z dw`, `I_WS = int W S dw`. The script then evaluates concrete cases: smooth Gaussian `W`, `S`; a w-dependent field mutation (to exhibit how non-zero-mode content leaks); the matched-observer case `W = Z/Z_int` with a delta-localized source (claiming projection produces `mu0/Z2_int`, a factor `sqrt(2)` larger than reduction); and a regulated sharp-observer / sharp-source case (eps -> 0) showing the `I_WS` integral diverges, so exact delta/delta is ill-defined.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 79-82 | `I_WZ_smooth - lam/sqrt(lam**2 + sigma**2) == 0` | yes |
| A2 | sympy | 83-86 | `I_WS_smooth - 1/(sqrt(pi)*sqrt(sigma**2+tau**2)) == 0` | yes |
| A3 | sympy | 87-90 | `int W*(Z*df - mu0*S*j) dw - (I_WZ*df - mu0*I_WS*j) == 0` | partial (linearity of `w`-integral; cannot fail as written) |
| A4 | sympy | 95-98 | field mutation delta `- eta*lam^3*sigma^2/(2*(lam^2+sigma^2)^{3/2}) == 0` | yes |
| A5 | sympy | 99 | field mutation delta `!= 0` | yes |
| A6 | sympy | 104-108 | source mutation delta `- eta*mu0*x*sigma^2*tau^2/(2*sqrt(pi)*(sigma^2+tau^2)^{3/2}) == 0` | yes |
| A7 | sympy | 109 | source mutation delta `!= 0` | yes |
| A8 | sympy | 117 | `Z_int - sqrt(pi)*lam == 0` | yes |
| A9 | sympy | 118 | `Z2_int - sqrt(2*pi)*lam/2 == 0` | yes |
| A10 | sympy | 133 | `I_WZ_match - Z2_int/Z_int == 0` | partial (algebraic restatement of W_match := Z/Z_int) |
| A11 | sympy | 134 | `I_WZ_match - sqrt(2)/2 == 0` | yes |
| A12 | sympy | 135 | `mu0_proj_match/mu0_red - sqrt(2) == 0` | yes |
| A13 | sympy | 136 | `mu0_proj_match/mu0_red - 1 != 0` | partial (numeric sanity guard) |
| A14 | sympy | 155 | `lim_{eps->0+} I_WZ_eps - 1 == 0` | yes |
| A15 | sympy | 156 | `I_WS_eps - sqrt(2)/(2*sqrt(pi)*eps) == 0` | yes |
| A16 | sympy | 157-160 | `lim_{eps->0+} I_WZ_eps - sqrt(2)/2 != 0` | partial (numeric sanity guard) |

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Subtype:** `missing_mathematica`
**Files:**
- `(missing)` — no `.wl` script exists for unit 007

**What's wrong:**
The unit ships only a SymPy verification script (`moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py`). No Mathematica script accompanies it. A search of `/var/projects/toy_physics/research/pde_ledger/scripts/` for any `stage007*` file returns the sympy file only. Per the audit prompt, unit 007 is marked `is_status_only_candidate: False`, so the second-engine policy applies: both engines must independently derive the unit's results.

The SymPy script itself is internally substantive: it verifies (i) closed-form Gaussian overlap integrals `I_WZ`, `I_WS` for smooth `W`, `S`, `Z`, (ii) explicit Gaussian-moment computations for `w`-dependent field and source mutations against analytic targets, (iii) `Z_int = sqrt(pi)*lambda` and `Z2_int = sqrt(2*pi)*lambda/2`, (iv) `mu0_proj_match / mu0_red = sqrt(2)`, and (v) the regulator analysis `I_WZ(eps) = lambda/sqrt(eps^2 + lambda^2)` with `I_WS(eps) = sqrt(2)/(2*sqrt(pi)*eps) -> infinity`. All of these need an independent Mathematica derivation to satisfy the second-engine cross-check.

**Why this matters:**
Without a second engine, the entire numerical/symbolic story (especially the `sqrt(2)` ratio between projection-first and reduction-first couplings, which is the central qualitative claim of the unit) rests on a single CAS. Any unnoticed SymPy quirk in the Gaussian moment evaluations or in `sp.DiracDelta` handling at line 124 (`I_WS_match = W_match.subs(w, 0)` is a hand-installed delta-sampling, not an integration — only Mathematica can adversarially recompute this independently) would propagate undetected.

**Required change:**
Create `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl` that independently verifies the claims listed in the Claim manifest in the directive. The script must not be a line-by-line transliteration of the Python (no shared intermediate variable names, no identical algebraic choreography). It should derive the Gaussian integrals via `Integrate[..., {w, -Infinity, Infinity}]` natively in Mathematica and form its own `mu0_eff` ratios.

**Verification:**
After Codex creates the file, `redteam exec-mathematica 007` must succeed with exit code 0 and produce a saved output at `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.txt` containing numeric/symbolic matches to the SymPy results: `Z_int = sqrt(pi) lambda`, `Z2_int = sqrt(pi/2) lambda`, `I_WZ_match = 1/sqrt(2)`, `mu0_proj_match/mu0_red = sqrt(2)`, `I_WZ_eps = lambda/sqrt(eps^2 + lambda^2)`, and the field/source mutation Gaussian-moment closed forms used at lines 95-108.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` script exists. See F1.

## Engine cross-check

Not applicable — only the SymPy engine is present. See F1.

## Verdict justification

The SymPy script is internally correct and non-trivially exercises its claims. The Gaussian overlap integrals at lines 79-86 match by direct hand-evaluation (`int exp(-w^2/sigma^2)*exp(-w^2/lambda^2)/(sqrt(pi)*sigma) dw = lambda/sqrt(sigma^2+lambda^2)`, etc.), the second-moment integrals at lines 95-108 reproduce `sigma^2*lambda^3 / (2*(sigma^2+lambda^2)^{3/2})` correctly, the regulator computations at lines 145-160 are consistent (`I_WZ_eps = lambda/sqrt(eps^2+lambda^2) -> 1` as `eps -> 0`; `I_WS_eps = sqrt(2)/(2*sqrt(pi)*eps)` diverges), and the central qualitative claim — that the matched-observer projection coupling exceeds the reduction-first coupling by a factor `sqrt(2)` — is anchored to `Z2_int/Z_int = sqrt(2)/2` plus `W_match(0) = 1/(sqrt(pi)*lambda)` via the delta-sampling identity. The few mildly tautological assertions (A3, A10, A13, A16) are sanity guards layered on top of the substantive checks; they do not rise to standalone `tautological_check` findings because the meat of each section is verified elsewhere by non-tautological assertions (A1, A2, A4, A6, A11, A12, A14, A15). Saved output (mtime 2026-05-11 12:38) is newer than the script (mtime 2026-05-04 12:00) so freshness is fine. The single finding is the missing Mathematica engine, blocking the second-engine policy.
