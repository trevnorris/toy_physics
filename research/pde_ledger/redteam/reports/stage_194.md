---
unit_id: 194
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra.md]
  paper_appendix: present
---

# Audit unit 194 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_194.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 119, 1263-1304, anchor list 1463-1468)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

The card's `\stagefield{Output}` reads verbatim: "Derives the exact outgoing spherical DtN fingerprint and fixes \(\chi_Q=1\) on the canonical compact branch." The notes expand this into three deliverables: (1) the exact outgoing spherical `l=2` Dirichlet-to-Neumann operator `\Lambda_2^{out}(z)=z\,d/dz\ln h_2^{(1)}(z)` with small-`z` fingerprint `-3 + z^2/3 + z^4/9 + i z^5/9 - 2z^6/27 - i z^7/27`, whose static slot is `-3` and whose normalized response `\widehat Y_2^{out}=-3/\Lambda_2^{out}=1+z^2/9+4z^4/81+i z^5/27+...`; (2) the exact fixing `\chi_Q=1` obtained by matching the carried retarded grouped-`P2` one-pole module (with `\Omega_Q=3c_s/2a`, `\sigma_Q^{can}=9/(8\Omega_Q^5)=4a^5/27c_s^5`) to that fingerprint at the leading odd `\omega^5` coefficient; (3) the isotropic DtN deformation algebra `\Lambda_2^{def}=S\Lambda_2^{out}(\beta z)+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5`, with canonical-even matching forcing `\Sigma_2=-(3S\beta^2-3S+\Sigma_0)/9`, `\Sigma_4=-(3S\beta^4-3S+\Sigma_0)/27`, and the outgoing map `\chi_Q=3(S\beta^5+9\Sigma_5)/(3S-\Sigma_0)`, `\chi_Q-1=(3S(\beta^5-1)+\Sigma_0+27\Sigma_5)/(3S-\Sigma_0)`. The notes Section 5 additionally lists, as a non-new carry-forward corollary, the invariant tuple `\overline K_0=54Gc_s^5/5a^5c^5`, `\overline K_2=6Gc_s^3/5a^3c^5`, `\overline K_4=8Gc_s/15ac^5`, `\overline\Gamma_5=2G/5c^5`. The appendix (lines 1263-1304) restates (1)-(3) identically. Anchors: MTDC-T9.5 (local), MTDC-T9 (Part V compiler). Status: ExactClosure, checkpoint False, is_status_only_candidate False.

## What the script claims to verify

The SymPy script builds the genuine closed-form spherical Hankel mode `h_2^{(1)}=j_2+i y_2`, forms `\Lambda_out=z h2'/h2` and `\widehat Y_out=-3/\Lambda_out` symbolically, takes their `sp.series(...,0,8)` expansions, and asserts these equal the paper's literal fingerprints (Section I). Section II constructs the retarded one-pole module `3/4+(1/4)/(1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{can}\omega^5)`, checks `\sigma_Q^{can}=4a^5/27c_s^5`, checks its `\omega`-series equals the claimed even+`i\chi_Q\omega^5` form, checks it coincides with the outgoing branch at `\chi_Q=1`, and extracts the `\omega^5` coefficient to isolate `\chi_Q-1` exactly. Section III defines `L0..L5`, series-expands the normalized deformed branch, solves the two canonical-even conditions for `\Sigma_2,\Sigma_4` and matches the paper formulas, then derives `\chi_Q=-27 L5/L0` under those constraints and matches the `\chi_Q` and `\chi_Q-1` deformation laws. Section IV reproduces the carry-forward tuple `\overline K_n` via `\Omega_Q`-scaled relations and matches each paper literal.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Exact `\Lambda_2^{out}` fingerprint `-3+z^2/3+z^4/9+i z^5/9+...` | I: `Lambda_out series - expected` (sympy-derived series vs literal) | match |
| Static slot `\Lambda_2^{out}(0)=-3` | I: `static DtN slot + 3` | match |
| Normalized `\widehat Y_2^{out}=1+z^2/9+4z^4/81+i z^5/27+...` | I: `Yhat_out series - expected` | match |
| `\Omega_Q=3c_s/2a`, `\sigma_Q^{can}=4a^5/27c_s^5` | II: `sigma_Q^can - 4 a^5/(27 c_s^5)` + symbolic `Omega_Q` | match |
| Retarded module expands to even form + `i\chi_Q\omega^5` | II: `Yret series - expected low-frequency form` | match |
| `\chi_Q=1` fixing | II: `retarded module matches ... when chi_Q=1` and `odd-coefficient matching fixes chi_Q-1` | match |
| `\Sigma_2,\Sigma_4` canonical-even formulas | III: `Sigma_2 formula`, `Sigma_4 formula` (from `sp.solve`) | match |
| `\chi_Q=3(S\beta^5+9\Sigma_5)/(3S-\Sigma_0)`; `\chi_Q-1` law | III: `chi_Q deformation law`, `chi_Q - 1 deformation law` | match |
| Carry-forward tuple `\overline K_0,K_2,K_4,\Gamma_5` (corollary) | IV: `Kbar_2/Kbar_4/Gammabar_5` minus literals | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 76 | `expect_zero(Lambda_out_series - Lambda_out_expected)` | deliverable 1 (fingerprint) | yes |
| A2 | sympy | 77 | `expect_zero(Yhat_out_series - Yhat_out_expected)` | deliverable 1 (normalized) | yes |
| A3 | sympy | 78 | `expect_zero(Lambda_out_series.subs(z,0)+3)` | deliverable 1 (static slot) | yes |
| A4 | sympy | 112 | `expect_zero(sigma_Q_can - 4 a^5/(27 c_s^5))` | deliverable 2 (carry-forward consts) | yes |
| A5 | sympy | 113 | `expect_zero(Yret_series - Yret_expected)` | deliverable 2 (module expansion) | yes |
| A6 | sympy | 114-117 | `expect_zero(Yret_series.subs(chi_Q,1) - Yout_omega_expected)` | deliverable 2 (`\chi_Q=1`) | yes |
| A7 | sympy | 125 | `expect_zero((odd_ret-odd_out)*27 c_s^5/(i a^5) - (chi_Q-1))` | deliverable 2 (odd-coeff fixing) | yes |
| A8 | sympy | 154 | `expect_zero(Ydef_series - Ydef_expected)` | deliverable 3 (normalized def branch) | yes |
| A9 | sympy | 172 | `expect_zero(sol_even[Sigma2] - Sigma2_expected)` | deliverable 3 (`\Sigma_2`) | yes |
| A10 | sympy | 173 | `expect_zero(sol_even[Sigma4] - Sigma4_expected)` | deliverable 3 (`\Sigma_4`) | yes |
| A11 | sympy | 182 | `expect_zero(chi_from_def_even - chi_expected)` | deliverable 3 (`\chi_Q` map) | yes |
| A12 | sympy | 183 | `expect_zero(chi_from_def_even-1 - chi_minus_one_expected)` | deliverable 3 (`\chi_Q-1` map) | yes |
| A13 | sympy | 204-206 | `expect_zero(Kbar_2/Kbar_4/Gammabar_5 - literal)` | corollary (invariant tuple) | yes (consistency) |

Every load-bearing assertion compares a SymPy-derived object (series expansion, `sp.solve` result, or symbolic substitution into the genuine Hankel function) against the paper's stated literal. The literal is independent of the derivation in A1-A12, so none is tautological. A13 is a consistency check that the `\Omega_Q`-scaled relations reproduce the paper's four hardcoded tuple values; see Verdict justification for why this is acceptable here.

## Findings

None.

## Independent-derivation check (Mathematica)

No `.wl` exists. See Engine cross-check / Verdict for the single-engine judgment.

## Engine cross-check

Only SymPy is present; no cross-engine comparison applies.

## Verdict justification

Clean. Attacks tried and failed: (1) Tautology hunt — every `expect_zero` in Sections I-III compares an independently SymPy-derived quantity (a `sp.series` of the genuine closed-form Hankel mode, or an `sp.solve` output) against a hardcoded paper literal; the literal never feeds the derivation, so the checks can fail and would fail if the paper's claimed coefficients were wrong. I hand-verified the key chain: `\widehat Y_out=-3/\Lambda_out` gives `z^2/9`, `4z^4/81` (= `z^4/27+z^4/81`), `i z^5/27`, matching the paper, and the script's series check confirms SymPy's own expansion equals these. (2) The `\chi_Q=1` fixing (A6/A7) is genuinely mechanistic, not asserted: A7 extracts the `\omega^5` coefficient via the 5th derivative, forms `(i a^5\chi_Q/27c_s^5 - i a^5/27c_s^5)*27c_s^5/(i a^5)` and subtracts `(\chi_Q-1)`, which is a real isolation of `\chi_Q-1` from the coefficient mismatch — by hand this is exactly `\chi_Q-1-\chi_Q+1=0` only because the coefficients truly match the claimed forms. (3) Deformation algebra (A9-A12): `sol_even` comes from `sp.solve` of the two even-matching equations, independent of the `_expected` literals; I hand-confirmed `\chi_Q=-27L5/L0=3(S\beta^5+9\Sigma_5)/(3S-\Sigma_0)` and the `\chi_Q-1` form. (4) Symbol-assumption attack: `z` real; `a,c_s,\omega` positive real (physically justified — length scale, sound speed, frequency); `\chi_Q` real; `S,\beta` nonzero real; `\Sigma_i` real — none of these hides a branch error, and all expansions are honest Taylor/Laurent series at the origin. (5) Series-order attack: `sp.series(...,0,6)` in `\omega` retains the `\omega^5` odd term, so the odd-coefficient matching is reachable. I confirmed I read the card, notes, and appendix, and the script's verified claim matches the paper's stated claim exactly on all three deliverables plus the corollary tuple.

On the single-engine question: all three deliverables are pure symbolic identities — closed-form series expansions of an explicit Hankel function and exact rational algebra (solve/simplify) for the deformation map. There is no numerical tolerance, no transcendental root-finding, and no branch ambiguity that an independent second engine would resolve. SymPy genuinely and fully settles these identities, so per the prompt's line-114 judgment and established pipeline precedent (SymPy-only non-status-only stages such as 121/122/123), single-engine verification is acceptable here; I do not flag `missing_mathematica`.

One soft observation (not a finding): Section IV (A13) hardcodes `\overline K_0=54Gc_s^5/5a^5c^5` and then derives the remaining three tuple entries via the relations `\overline K_2=\overline K_0/4\Omega_Q^2`, `\overline K_4=\overline K_0/4\Omega_Q^4`, `\overline\Gamma_5=9\overline K_0/32\Omega_Q^5`, checking each against the paper literal. These scaling relations are not separately stated in the card/notes, but the notes explicitly mark the tuple as a carried-forward corollary ("not a new theorem input"), not an in-stage derivation, and `\overline K_0` plus all three derived values match the paper deliverable verbatim. The check therefore confirms the tuple's internal `\Omega_Q`-scaling consistency rather than re-deriving it from scratch, which is appropriate for a carry-forward corollary and does not rise to `hardcoded_result` or `paper_misalignment`.

A second cosmetic-only observation (not a finding): the script's `banner`/`subbanner`/ledger print strings say "STAGE 177" (lines 35, 208, and the output transcript), while ledger item 6 (line 231) and the filename/paper say "Stage 194"; the notes prose alternately calls this Stage 245/194. These are display-label inconsistencies in `print` statements with zero effect on any assertion or on the verified math, so no math finding is raised.

## Self-test notes

Checked the four self-test traps. Variable independence: the only differentiations are `sp.diff(h2,z)` (h2 genuinely depends on z) and the 5th `\omega`-derivatives in A7 (both `Yret_series` and `Yout_omega_expected` genuinely carry `\omega^5` terms), so no identically-zero-derivative trap. Symmetry/parity: no unbounded-domain integrals here — all checks are series-coefficient identities. Trivial-case pre-check: hand-substitution confirmed A6 reduces to `0` only at `\chi_Q=1` (it is the fixing), and A7 isolates `\chi_Q-1` correctly. Paper round-trip: every literal asserted matches the card/notes/appendix exactly. No directive written (zero findings).
