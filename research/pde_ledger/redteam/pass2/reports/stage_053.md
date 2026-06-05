---
unit_id: 053
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage053_overlap_boost_window.md]
  paper_appendix: present
---

# Audit unit 053 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_053.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage053_overlap_boost_window.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 84)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage053_overlap_boost_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage053_overlap_boost_mathematica_audit.txt`

## What the paper claims

Stage 053 computes the maximum overlap boost the first exponential lowest-lane source family can supply, and the universal overlap ceiling. With the lowest D/N mode `chi_0(s)=sqrt(2/L) sin(pi s/2L)` and the normalized exponential source `sigma_alpha(s)=alpha e^{alpha s/L}/(e^alpha-1)` (total source strength `=L`), the `\stagefield{Output}` is: the boxed overlap-boost formula `Omega_exp(alpha)=pi*alpha(2*alpha*e^alpha+pi)/[(4*alpha^2+pi^2)(e^alpha-1)]` (eq. app-stage053-Omega-exp) and the boxed ceiling `Omega_0^2 <= pi^2/4` (eq. app-stage053-overlap-ceiling). Two further boxed claims: the endpoint values `Omega_exp(0)=1` and `lim_{alpha->inf}Omega_exp = pi/2`. The notes additionally state the building blocks `I_W=2 sqrt(2L)/pi`, `max chi_0 = sqrt(2/L)`, `Omega_max = pi/2`, the small-asymmetry linear coefficient `(4-pi)/(2pi)`, and the rescue criterion: pure overlap rescue is possible only if `zeta_req <= pi^2/4`. The appendix row (line 84) summarizes: "Exponential source family and ceiling `Omega_0^2<=pi^2/4`."

## What the script claims to verify

Both engines re-derive the geometric quantities by direct symbolic integration: `I_W = integral chi_0 = 2 sqrt(2L)/pi`, `chi_0_max = sqrt(2/L)`, then `Omega_max = L*chi_0_max/I_W` and `A_I,max = Omega_max^2`, asserting `Omega_max - pi/2 == 0` and `A_I,max - pi^2/4 == 0`. They then build `sigma_alpha`, assert it has unit-throat total strength (`integral - L == 0`), integrate `I_alpha = integral sigma_alpha*chi_0`, form `Omega_alpha = I_alpha/I_W`, and assert it equals the paper's boxed closed form (`Omega_alpha - Omega_alpha_simpler == 0`). They take the two limits (asserting `Omega(0)-1==0`, `Omega(inf)-pi/2==0`), expand the small-alpha series and assert the linear coefficient equals `(4-pi)/(2pi)`. The rescue criterion `A_I,max - zeta_req` is only printed (a definitional restatement of the already-asserted `A_I,max=pi^2/4`).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Boxed `Omega_exp(alpha)` closed form | sympy L63 / wl L65 assert `Omega_alpha - Omega_alpha_simpler == 0` (integration vs closed form) | match |
| Boxed ceiling `Omega_0^2 <= pi^2/4` | sympy L42-43 / wl L48-49 assert `Omega_max=pi/2`, `A_I,max=pi^2/4` | match |
| Boxed endpoints `Omega_exp(0)=1`, `lim=pi/2` | sympy L69-70 / wl L76-77 assert both | match |
| `sigma_alpha` unit total strength `=L` | sympy L57 / wl L64 assert `Sigma_total - L == 0` | match |
| `I_W = 2 sqrt(2L)/pi`, `max chi_0 = sqrt(2/L)` (notes) | computed + printed; feed asserted `Omega_max` | match |
| small-alpha linear coeff `(4-pi)/(2pi)` (notes) | sympy L76 / wl L78 assert | match |
| rescue criterion `zeta_req <= pi^2/4` (notes) | printed only (sympy L82-83); content = asserted `A_I,max=pi^2/4` | match (definitional) |

`paper_alignment: aligned` — every load-bearing deliverable has a faithful, non-tautological script-side check, and all values reconcile.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42 | `Omega_max - pi/2 == 0` | ceiling (via Omega_max) | yes |
| A2 | sympy | 43 | `A_I,max - pi^2/4 == 0` | ceiling | yes |
| A3 | sympy | 57 | `Sigma_total - L == 0` | equal-strength normalization | yes |
| A4 | sympy | 63 | `Omega_alpha - closed_form == 0` | boxed Omega_exp formula | yes |
| A5 | sympy | 69 | `Omega0 - 1 == 0` | endpoint Omega(0)=1 | yes |
| A6 | sympy | 70 | `Omegainf - pi/2 == 0` | endpoint lim=pi/2 | yes |
| A7 | sympy | 76 | `linear_coeff - (4-pi)/(2pi) == 0` | small-alpha expansion (notes) | yes |
| A8 | mathematica | 48 | `omegaMax - Pi/2 == 0` | ceiling | yes |
| A9 | mathematica | 49 | `aIMax - Pi^2/4 == 0` | ceiling | yes |
| A10 | mathematica | 64 | `sigmaTotal - ell == 0` | normalization | yes |
| A11 | mathematica | 65 | `omegaAlpha - omegaAlphaSimple == 0` | boxed formula | yes |
| A12 | mathematica | 76 | `omega0 - 1 == 0` | endpoint | yes |
| A13 | mathematica | 77 | `omegaInf - Pi/2 == 0` | endpoint | yes |
| A14 | mathematica | 78 | `linearCoeff - (4-Pi)/(2Pi) == 0` | small-alpha expansion | yes |

No assertion is tautological: every quantity entering an assertion is produced by an independent symbolic integration / limit / series expansion, then compared to the paper's closed form. A4/A11 in particular are the strong checks — the boxed formula is compared against a freshly-integrated `I_alpha/I_W`, so a wrong boxed coefficient would fail. The `zeta_req` criterion is print-only (sympy L82-83, wl L82-83) but its mathematical content (`A_I,max=pi^2/4`) is already asserted by A2/A9, so it needs no separate assertion.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage053_overlap_boost_sympy_audit.txt` (mtime 2026-05-22 17:37)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage053_overlap_boost_mathematica_audit.txt` (mtime 2026-05-22 17:37)
- vs scripts `..._sympy_audit.py` and `..._mathematica_audit.wl` (mtime 2026-06-03 15:59)

**What's wrong:**
Both committed output transcripts predate the current scripts by ~12 days, and the content disagrees in the banner self-label. The `.wl` script now prints `STAGE 053 — EXACT OVERLAP-BOOST WINDOW` (wl:32) but the saved Mathematica output still shows `STAGE 036 — EXACT OVERLAP-BOOST WINDOW` (mathematica output:3). The SymPy script banner currently reads `STAGE 53` (py:25) while the saved SymPy output reads `STAGE 36` (sympy output:3). Confirms the transcripts are stale.

Additionally there are stale numbering self-labels in the script sources themselves (a known low-severity numbering-label class, non-blocking — orchestrator resolves scope):
- SymPy docstring `py:3` — `"""Stage 36 SymPy audit: ..."""` (should be 053).
- SymPy banner `py:25` — `banner("STAGE 53 — ...")` is unpadded ("53" not "053"); the `.wl` uses the correct `STAGE 053`.

All numeric/symbolic results in the stale transcripts otherwise match what the current scripts compute (verified line-by-line against the source); only the banner label and zero-padding differ.

**Why this matters:**
The transcripts are advertised in the stage card `\stagefield{Verification}` as the audit record; a stale banner/label undermines traceability even though the math is unchanged. Refreshing the outputs is cheap and the verifier re-runs anyway.

**Required change:**
Re-run both scripts to regenerate the committed `.txt` transcripts (the orchestrator's independent re-run handles this). Optionally fix the two SymPy self-labels: `py:3` docstring `Stage 36 -> Stage 053`, and `py:25` banner `STAGE 53 -> STAGE 053` to match the `.wl`.

**Verification:**
After refresh, the SymPy output banner reads `STAGE 053` (matching the fixed `py:25`), the Mathematica output banner reads `STAGE 053` (matching wl:32), and all PASS lines remain present.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration: it performs its own symbolic integration in the Wolfram engine (`Integrate[chi0,{s,0,ell}]` wl:38, `Integrate[sigmaAlpha chi0,{s,0,ell}]` wl:53) rather than echoing SymPy's algebra. It uses `ell` (vs SymPy `L`), `FullSimplify`/`Together` (vs `simplify`/`expand`), and `Series[...,{alpha,0,2}]`/`Coefficient` for the small-alpha coefficient. The parallel high-level structure (same seven physical checks in the same order) is inherent to the stage's short linear derivation, not a copied variable choreography — each engine reaches the boxed forms by its own integral evaluation. Acceptable second engine.

## Engine cross-check

Final symbolic forms agree exactly:
- `Omega_alpha`: sympy `pi*alpha*(2*alpha*exp(alpha)+pi)/((4*alpha**2+pi**2)*(exp(alpha)-1))` (sympy output:31) == wl `(alpha*Pi*(2*alpha*E^alpha+Pi))/((-1+E^alpha)*(4*alpha^2+Pi^2))` (mathematica output:17).
- `I_alpha`, `I_W`, `Omega_max=pi/2`, `A_I,max=pi^2/4`, endpoints `1` and `pi/2`, linear coeff `(4-pi)/(2pi)` — all identical across both transcripts; every assertion is `0`/PASS. The Mathematica `Limit::alimv` warnings (mathematica output:23-25) are benign (limit-variable assumptions ignored); the limits still evaluate to `1` and `Pi/2`. Engines agree.

## Verdict justification

verdict: findings — the only finding is a low-severity `stale_output` (transcripts predate the scripts; banner self-label drifted 036/53→053, plus two stale SymPy self-labels). No script-side math defect: all seven physical deliverables (boxed `Omega_exp` formula, ceiling `pi^2/4`, both endpoints, unit normalization, small-alpha linear coefficient) are asserted non-tautologically by comparing fresh symbolic integrals/limits/series against the paper's closed forms, and the two engines agree. I attacked the closed-form check for tautology (A4/A11) — it is not: `Omega_alpha` is the integrated `I_alpha/I_W`, compared to the hand-written `Omega_alpha_simpler`, so a wrong boxed coefficient would surface. I checked the `positive=True`/`ell>0`,`alpha>0` assumptions — they are justified by the physical setup (`L>0` throat length, `alpha>0` bottom-bias) and do not over-constrain the integrals. The print-only `zeta_req` criterion is backed by the already-asserted `A_I,max=pi^2/4`. Paper, notes, and appendix all read and consistent.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `I_W = 2 sqrt(2L)/pi` | py:31 / sympy out:11; wl:38 / wl out:6 | notes:48 `I_W = 2 sqrt(2L)/pi` | MATCH |
| `max chi_0 = sqrt(2/L)` | py:32 / sympy out:12; wl out:7 | notes:74,78 `max chi_0 = sqrt(2/L)` | MATCH |
| `Omega_max = pi/2` | py:33,42 / sympy out:13; wl out:8 | tex:33 `lim=pi/2`; notes:82,86 `Omega_0 <= pi/2` | MATCH |
| `A_I,max = pi^2/4` | py:34,43 / sympy out:14; wl out:9 | tex:39 `Omega_0^2<=pi^2/4`; notes:90 `A_I<=pi^2/4`; appendix:84 | MATCH |
| `sigma_alpha total = L` | py:48,57 / sympy out:29; wl out:18 | tex:18 `int sigma=L`; notes:110 | MATCH |
| `I_alpha = 2 sqrt(2L) alpha(2 alpha e^a+pi)/[(4a^2+pi^2)(e^a-1)]` | py:49 / sympy out:30; wl out:16 | notes:115-116 (identical) | MATCH |
| `Omega_exp(alpha) = pi alpha(2 alpha e^a+pi)/[(4a^2+pi^2)(e^a-1)]` | py:56,59 / sympy out:31; wl out:17 | tex:23-25 (boxed); notes:122-123,183 | MATCH |
| `Omega_exp(0) = 1` | py:65 / sympy out:34; wl out:26 | tex:31; notes:127 | MATCH |
| `lim_{alpha->inf} Omega_exp = pi/2` | py:66 / sympy out:35; wl out:27 | tex:33; notes:129 | MATCH |
| small-alpha linear coeff `(4-pi)/(2pi) = 2/pi-1/2` | py:74 / sympy out:39; wl out:29 | notes:135-138 `(2/pi-1/2)` | MATCH |
| rescue criterion `zeta_req <= pi^2/4` | py:82 / sympy out:45; wl:83 / wl out:36 | tex:41; notes:156,188 | MATCH |

INTERNAL (scaffolding, no prose expected): `chi0(s)` echo print, residual/`expect_zero` "= 0" lines, PASS flags, the full small-alpha series quadratic term (only the linear coeff is a stated deliverable; the `alpha^2` term is intermediate).

reconciliation: complete; 11 values checked, 0 misaligned

## Self-test notes

Variable-independence: no `diff`/`D` on a wrong variable here — all checks are integrals/limits/series in `s`/`alpha`, each genuinely dependent. Parity/symmetry: the integrals are over the bounded `[0,L]` with strictly positive integrands, so no symmetric-domain vanishing trap. Trivial-case: at `alpha->0` the integrated `Omega_alpha` correctly limits to `1` (matches A5/A12) and the closed form's `alpha->0` agrees; the boxed-formula check A4/A11 compares an independently-integrated ratio against the literal, so it cannot pass trivially. No directive-prescribed new check is needed (the sole finding is informational stale_output).
