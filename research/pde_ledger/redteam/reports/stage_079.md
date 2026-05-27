---
unit_id: 079
batch: III.4
auditor_model: claude-opus-4-7[1m]
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage079_family1_quadrupole_pe_map.md
  paper_appendix: present
---

# Audit unit 079 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_079.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage079_family1_quadrupole_pe_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 136)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage079_family1_quadrupole_pe_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage079_family1_quadrupole_pe_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.txt`

## What the paper claims

The paper card converts a required support ratio into the transport bias required on the explicit Family-1 branch. `\stagefield{Output}{Family--1 demand map \eqref{eq:app-stage079-zetaF1} and ceiling \eqref{eq:app-stage079-ceiling}.}` The two boxed equations are `zeta_F1(Pe) = A_F1 Omega_Pe^2` and `zeta_max^(F1) = A_F1 pi^2/4 ≈ 2.46752922945601`, where `A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2)` with `kappa_F1 = 12321/5`, `y_F1 tan(y_F1) = 37` on `0 < y_F1 < pi/2`, giving `y_F1 ≈ 1.52948248371470` and `A_F1 ≈ 1.00005192880220`. The notes add (Section 2) the explicit closed form `Omega_Pe = pi Pe (2 Pe exp(Pe) + pi)/[(4 Pe^2 + pi^2)(exp(Pe)-1)]` with continuous extension `Omega_0 = 1` and (Section 4) the small-Pe expansion `Omega_Pe = 1 + ((4-pi)/(2 pi)) Pe + O(Pe^2)`, hence `zeta_F1(Pe) = A_F1 [1 + ((4-pi)/pi) Pe + O(Pe^2)]`. The part appendix row at line 136 summarizes: "Demand-to-transport map ... `zeta_F1(Pe) = A_F1 Omega_Pe^2` and hard ceiling."

## What the script claims to verify

Both engines instantiate the same Family-1 inputs (`kappa_F1 = 12321/5`, `eta_F1 = 37`), solve the Robin root `y_F1 tan(y_F1) = 37` near 1.53, form `A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2)`, define the closed form for `Omega(Pe)`, and verify: (i) `Omega(0+) = 1`, (ii) `Omega(inf) = pi/2`, (iii) `zeta_F1(0+) = A_F1`, (iv) `zeta_F1(inf) = A_F1 pi^2/4`, (v) the small-Pe linearization `zeta_F1(Pe) = A_F1 (1 + ((4-pi)/pi) Pe) + O(Pe^2)`. The Mathematica file additionally exercises an independent slope check `Omega'(0+) = (4-Pi)/(2 Pi)` via symbolic differentiation `D[omega, pe]`, a code path distinct from SymPy's `series`. SymPy carries the Robin residual at runtime (`~1.4e-14`); Mathematica carries it at WorkingPrecision 80 (`~2.3e-55`).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `kappa_F1 = 12321/5` (Inputs) | sympy:33 / wl:37 — both assign `Rational(12321,5)` / `12321/5` | match |
| `eta_F1 = 37` (Inputs) | sympy:34 / wl:38 — both assign `Integer(37)` / `37` | match |
| `y_F1 tan(y_F1) = 37`, `0 < y_F1 < pi/2` | sympy:36 `nsolve(...,1.53,...)` / wl:39 `FindRoot[..., {y,1.53}, WorkingPrecision->80]`; Robin residual asserted | match |
| `y_F1 ≈ 1.52948248371470` | sympy output: `1.52948248371469963...`; wl output: `1.5294824837146996499...` | match |
| `A_F1 = (kappa+pi^2/4)/(kappa+y_F1^2)` | sympy:37 / wl:40 — identical formula | match |
| `A_F1 ≈ 1.00005192880220` | sympy output: `1.00005192880219520...`; wl output: `1.00005192880219532865...` | match (to ~15 digits) |
| Demand map `zeta_F1(Pe) = A_F1 Omega_Pe^2` | sympy:55-63 builds and limit-checks `A_F1*Omega^2`; wl:67-70 does `Limit[aF1*omega^2,...]` on product | match |
| Ceiling `zeta_max^(F1) = A_F1 pi^2/4 ≈ 2.46752922945601` | sympy:63 `expect_small("... - A_F1*pi^2/4", zeta_inf - A_F1*pi^2/4)`; wl:70 `zetaInfSym = Limit[aF1*omega^2, pe -> Infinity]` checked vs `aF1*Pi^2/4` | match |
| Small-Pe expansion `zeta_F1 = A_F1(1+((4-pi)/pi)Pe)+O(Pe^2)` (notes Section 4) | sympy:66-69 / wl:72-76 — both check series coefficients vs `aF1*(1+((4-Pi)/Pi)*pe)` | match |
| Notes Section 4 implicit `Omega'(0+) = (4-pi)/(2 pi)` | wl:77-79 — independent `Limit[D[omega,pe], pe -> 0]` vs `(4-Pi)/(2*Pi)` | match (Mathematica only; SymPy covers same content through its `series` call) |

Every paper-side deliverable has a matching script-side check; no rows are `partial`, `missing`, `mismatch`, or `extra`. Front-matter `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 52 | `expect_small("Omega(0+) - 1", Omega0 - 1)` (sp.limit→1) | underpins `zeta_F1(0+) = A_F1` | yes |
| A2 | sympy | 53 | `expect_small("Omega(inf) - pi/2", Omega_inf - pi/2)` | underpins ceiling `A_F1 pi^2/4` | yes |
| A3 | sympy | 62 | `expect_small("zeta_F1(0+) - A_F1", zeta0 - A_F1)` (limit of `A_F1*Omega^2`) | demand map at Pe=0 | yes |
| A4 | sympy | 63 | `expect_small("zeta_F1(inf) - A_F1 pi^2/4", zeta_inf - A_F1*pi^2/4)` | ceiling (Output, boxed eq) | yes |
| A5 | sympy | 69 | `expect_small("small-Pe coefficient check", series - A_F1*(1+((4-pi)/pi)*Pe))` | notes Sec. 4 linearization | yes |
| A6 | mathematica | 46 | `expectApprox["Robin residual", yF1*Tan[yF1] - etaF1, 0, 1e-28]` | Inputs: `y_F1 tan(y_F1) = 37` | yes |
| A7 | mathematica | 58 | `expectZero["Omega(0+) - 1", omega0 - 1]` (symbolic Limit) | underpins `zeta_F1(0+)` | yes |
| A8 | mathematica | 59 | `expectZero["Omega(inf) - Pi/2", omegaInf - Pi/2]` | underpins ceiling | yes |
| A9 | mathematica | 69 | `expectApprox["zeta_F1(0+) - A_F1", N[zeta0Sym - aF1, 50], 0, 1e-40]` where `zeta0Sym = Limit[aF1*omega^2, pe->0]` | demand map at Pe=0 | yes (post-fix; was `no` in v1) |
| A10 | mathematica | 70 | `expectApprox["zeta_F1(inf) - A_F1 Pi^2/4", N[zetaInfSym - aF1*Pi^2/4, 50], 0, 1e-40]` where `zetaInfSym = Limit[aF1*omega^2, pe->Infinity]` | ceiling (Output) | yes (post-fix; was `no` in v1) |
| A11 | mathematica | 75 | `expectApprox["small-Pe constant coefficient check", Coefficient[seriesDiff,pe,0], 0, 1e-28]` | notes Sec. 4 | yes |
| A12 | mathematica | 76 | `expectApprox["small-Pe linear coefficient check", Coefficient[seriesDiff,pe,1], 0, 1e-28]` | notes Sec. 4 | yes |
| A13 | mathematica | 79 | `expectApprox["Omega'(0+) - (4-Pi)/(2 Pi)", N[omegaPrime0 - (4 - Pi)/(2*Pi), 50], 0, 1e-40]` where `omegaPrime0 = Limit[D[omega,pe], pe->0]` | notes Sec. 4 (implicit `Omega'(0+)`) | yes (post-fix; new in v2) |

## Findings

(None.)

The v1 audit raised three findings (F1 `hardcoded_result`, F2 `tautological_check`, F3 `mathematica_transliteration`). The v1 directive landed all three:

- F1/F2: lines 67-68 of the Mathematica script no longer compare `zeta0`/`zetaInf` against hardcoded SymPy decimals; instead `zeta0Sym = FullSimplify[Limit[aF1*omega^2, pe -> 0, ...]]` and `zetaInfSym = FullSimplify[Limit[aF1*omega^2, pe -> Infinity]]` are formed, and `expectApprox` at tolerance `10^-40` compares them to `aF1` and `aF1*Pi^2/4`. The symbolic limit is now taken on the *product*, not on `omega` alone followed by multiplication. (Mathematica output: `zeta_F1(0+) - A_F1 diff = 0``49.7`, `zeta_F1(inf) - A_F1 Pi^2/4 diff = 0``49.3` — both reduce to numerical zero at 50-digit precision.)
- F3: lines 77-79 add the independent slope check `omegaPrime0 = FullSimplify[Limit[D[omega, pe], pe -> 0, ...]]` and `expectApprox["Omega'(0+) - (4-Pi)/(2 Pi)", N[omegaPrime0 - (4-Pi)/(2*Pi), 50], 0, 1e-40]`. Mathematica's output reports `Omega'(0+) = -1/2 + 2/Pi`, equal to `(4-Pi)/(2*Pi)`, and the residual is zero to 99-digit precision. This exercises symbolic differentiation (`D`), a different code path from SymPy's `sp.series`.

No new findings appear in v2. The paper-alignment gate is satisfied: every boxed paper equation, every `\stagefield{Inputs}` constant, every notes-derived expansion has a matching script-side check, and the numerical values in both transcripts match the paper's quoted constants to the precision the paper states (`A_F1 ≈ 1.00005192880220` matches to 14 digits; `zeta_max^(F1) ≈ 2.46752922945601` matches to 14 digits). The appendix row at part03 line 136 ("`zeta_F1(Pe)=A_F1 Omega_Pe^2` and hard ceiling") matches the script's bottom-line ledger output.

## Independent-derivation check (Mathematica)

The Mathematica script remains, by construction, a parallel computation of the same closed form for `Omega(Pe)` (neither script re-derives `Omega(Pe)` from the underlying PDE in this stage — both take the Stage-39 carry-forward as input). However, after the F3 fix the Mathematica script does at least one substantive computation through an independent code path:

- SymPy verifies the small-Pe linearization via `sp.series(zeta_F1, Pe, 0, 2)` (Taylor series).
- Mathematica verifies the same content via `Limit[D[omega, pe], pe -> 0, Direction -> "FromAbove"]` (symbolic differentiation followed by limit) and compares to `(4-Pi)/(2*Pi)`. This is a different symbolic-engine code path than `Series`; if Mathematica's `D` or `Limit` disagreed with SymPy's `series` on the closed form, the cross-check would catch it.

The closed-form-for-`Omega(Pe)` itself is taken as a Stage-39 carry-forward in both scripts; that is a legitimate carry-forward (not a transliteration) provided Stage 39's audit is sound. No `mathematica_transliteration` finding in v2.

## Engine cross-check

Both engines pass all checks. Side-by-side:

- Robin residual: SymPy `-1.42e-14` (within `1e-12` tolerance); Mathematica `~2.3e-55` (within `1e-28` tolerance). Both PASS.
- `Omega(0+) = 1`: SymPy exact zero; Mathematica exact zero. Both PASS.
- `Omega(inf) = Pi/2`: SymPy exact zero; Mathematica exact zero. Both PASS.
- `zeta_F1(0+) - A_F1`: SymPy `-3.3e-16`; Mathematica `0` at 49.7-digit precision. Both PASS.
- `zeta_F1(inf) - A_F1 Pi^2/4`: SymPy `-9.4e-16`; Mathematica `0` at 49.3-digit precision. Both PASS.
- Small-Pe series: SymPy `0`; Mathematica constant and linear coefficients both `0` at 49+-digit precision. Both PASS.
- `Omega'(0+) - (4-Pi)/(2 Pi)`: Mathematica only — reduces symbolically to `0` at 99-digit precision (Mathematica reports `Omega'(0+) = -1/2 + 2/Pi`). SymPy's series A5 lands on the same identity via a different route. Cross-engine consistency confirmed.

`engines_agree: true` in front-matter.

## Verdict justification

Clean. The v1 directive's three findings landed correctly and the v2 paper-alignment audit confirms every paper-side deliverable maps to a script-side check that is non-tautological and exercises the right derivation path. Attacks tried that failed:

- Re-derived the small-Pe slope by hand: numerator `pi Pe(2 Pe exp(Pe) + pi) ≈ pi^2 Pe + 2 pi Pe^2 + ...`; denominator `(4 Pe^2 + pi^2)(exp(Pe)-1) ≈ pi^2 Pe + pi^2 Pe^2/2 + 4 Pe^3 + ...`; ratio `≈ 1 + ((4-pi)/(2 pi)) Pe + O(Pe^2)`. Matches `Omega'(0+) = (4-Pi)/(2 Pi) ≈ 0.13662`. Mathematica's `D[omega, pe]` slope check is correctly anchored.
- Re-derived the `Pe → infinity` limit: numerator ~`2 pi Pe^2 exp(Pe)`; denominator ~`4 Pe^2 exp(Pe)`; ratio `→ pi/2`. Correct.
- Verified the Robin root: `y tan(y) = 37` with `y = 1.5294824837...` gives `1.5294824837... * tan(1.5294824837...) ≈ 37.0000...`. Correct.
- Checked numeric agreement of `A_F1`: SymPy 30-digit value `1.00005192880219520364707725466` vs Mathematica 50-digit value `1.00005192880219532865933408371...`. Agreement to ~15 digits is consistent with SymPy's mpmath default precision (the `tol=1e-34` in `sp.nsolve` is a stopping criterion, not a precision setting; the Robin residual of `1.4e-14` confirms only ~14-digit precision was actually achieved). Within all stage assertion tolerances. Not a finding.
- Checked the paper's quoted value `A_F1 ≈ 1.00005192880220`: matches Mathematica to all 15 quoted digits.
- Checked the paper's quoted `zeta_max^(F1) ≈ 2.46752922945601`: matches both engines to all 15 quoted digits.
- The `Limit::alimv` warnings in the Mathematica output ("Assumptions that involve the limit variable are ignored") are benign — they refer to `pe > 0` being dropped when `pe` is the limit variable, which doesn't invalidate the result.
- Output files are fresher than scripts (`stale_output` not triggered): script `.wl` mtime 2026-05-22 23:34:30; output `.txt` mtime 2026-05-22 23:35:35; SymPy `.py` mtime 2026-04-01 (older but the output mtime 2026-05-22 23:35:29 is fresh — script has not changed since first authored).

`paper_alignment: aligned`. No `paper_misalignment` of any subtype. `stop_cold: null`. `findings_count: 0`. `verdict: clean`. No directive written.

## Self-test notes

- Variable independence: this audit produced no new findings, so no `D[..., var]` or `sp.diff(..., var)` proposals to validate.
- Symmetry/parity: no new unbounded integrals introduced; not applicable.
- Trivial-case pre-check: for the existing F3 fix already in place, `Limit[D[omega, pe], pe -> 0]` evaluates to `(4-Pi)/(2 Pi)` per Mathematica's own output line — non-trivially nonzero, no `assert_zero` is trivially true. The existing `expectApprox` against `0` tests a real residual.
- Paper round-trip: every script-side check traces to a paper-side claim (Inputs constants, boxed equations, Output line, or notes Section 4). No script-side check exercises material the paper does not mention; no paper-side claim is left unchecked. The constants `12321/5`, `37`, and the numeric tail of `A_F1` / `zeta_max^(F1)` all match across paper and scripts.
