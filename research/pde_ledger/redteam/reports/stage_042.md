---
unit_id: 042
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 042 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.txt`

## What the script claims to verify

Under a rank-2 "support completion" with two loading directions, the script claims that the selected lower-mode eigenvector ratio `e1/e0` has a closed-form rank-2 expression (verified by matching two independent row-relations against the closed form), that the outgoing and source overlaps reduce exactly to the named expressions `Z_expected` and `S_expected`, and that the generalized normalization function `F_(q,r,t) = Z*S/(1-xi)` matches its closed form. It further claims that if the support tracks the mixed direction (`r=q`), `F` collapses exactly to the Stage-23 two-vector function; that the source-tied split-U specialization (`r=t=sqrt(lam0)`, `q=sqrt(lam0) R_U` with `lam0=2/9`) yields a specific closed form `F_src(xi,delta;m,R_U)`; that `F_src(R_U=1) = F_flat`; and that the first-order Taylor coefficients of `n_src` and `F_src/F_flat` in `R_U` about `R_U=1` match named closed forms `H_n_expected` and `H_F_expected`, and the full first-order series expansions agree with the linearized forms.

## Assertion inventory

| #  | Script       | Line | Form                                                          | Anchored to claim? |
|----|--------------|------|---------------------------------------------------------------|--------------------|
| A1 | sympy        | 72   | `simplify(ratio_row1 - ratio_expected) == 0`                  | yes                |
| A2 | sympy        | 73   | `simplify(ratio_row2 - ratio_expected) == 0`                  | yes                |
| A3 | sympy        | 89   | `simplify(Z_overlap - Z_expected) == 0`                       | yes                |
| A4 | sympy        | 90   | `simplify(S_overlap - S_expected) == 0`                       | yes                |
| A5 | sympy        | 101  | `simplify(F_general - F_expected) == 0`                       | yes                |
| A6 | sympy        | 112  | `simplify(F_track - F_stage23) == 0`                          | yes                |
| A7 | sympy        | 130  | `simplify(F_src_direct - F_src) == 0`                         | yes                |
| A8 | sympy        | 138  | `simplify(F_src(R_U=1) - F_flat) == 0`                        | yes                |
| A9 | sympy        | 152  | `simplify(H_n_src - H_n_expected) == 0`                       | yes                |
| A10| sympy        | 165  | `simplify(H_F_src - H_F_expected) == 0`                       | yes                |
| A11| sympy        | 171-174 | `simplify(series(n_src, R_U=1+eps, ord 2) - n_linear) == 0`| yes                |
| A12| sympy        | 175-178 | `simplify(series(F_ratio, R_U=1+eps, ord 2) - F_linear) == 0`| yes              |
| B1 | mathematica  | 56-57   | `expectZero` row1, row2 versus ratioExpected                 | yes                |
| B2 | mathematica  | 80-82   | `expectZero` Z_overlap, S_overlap, F_general                 | yes                |
| B3 | mathematica  | 94      | `expectZero` tracking collapse                                | yes                |
| B4 | mathematica  | 112     | `expectZero` source-tied specialization                       | yes                |
| B5 | mathematica  | 120     | `expectZero` F_src(rU=1) - F_flat                             | yes                |
| B6 | mathematica  | 146-147 | `expectZero` H_n, H_F                                          | yes                |
| B7 | mathematica  | 148-155 | `expectZero` linear series expansions of nSrc and fRatio      | yes                |

All rows verified non-tautological after attempting to break them (see verdict justification).

## Findings

None.

## Independent-derivation check (Mathematica)

The Mathematica script is closely parallel to the SymPy script in structure: same section headings (1-6), same intermediate variable names (`nReq <-> n_req`, `ratioRow1 <-> ratio_row1`, `dQR <-> D_qr`, `a1`, `b1`, `dSrc`, `fSrc`, `fRatio`, `hNSrc`, `hFSrc`), same comparison targets, and same evaluation order. This is parallel structure rather than first-principles re-derivation, but the script's claims are themselves a set of named closed-form algebraic identities (e.g., "this rational function equals that rational function under simplification"). Independent re-derivation of an algebraic identity is not meaningfully different across CAS engines — both engines build the same two closed forms and verify their residual under their respective simplification routines (`simplify+expand` for SymPy at line 38 of the sympy script; `FullSimplify[Together[Expand[...]], Assumptions -> $Assumptions]` for Mathematica at line 28 of the wl). The cross-check is real (different simplification kernels, different normal forms — Mathematica returns `H_F^(src)` in partial-fraction form while SymPy returns it as a single rational, see line 62 of the wl output vs line 132-136 of the py output, and both equal the named closed form when normalized). This is not flagged as `mathematica_transliteration`.

## Engine cross-check

Both engines report PASS on every assertion (sympy exit 0, output line 155; Mathematica exit 0, output line 74). Residuals all explicitly print as 0. Final closed forms agree when normalized:

- `e1/e0`: SymPy prints `-(m*q - m*r + r*xi)/(-delta + m*q^2 - m*q*r - xi)`; Mathematica prints `(m*(q-r) + r*xi)/(delta + m*q*(-q+r) + xi)` — identical under sign-flip of numerator and denominator.
- `H_n^(src)`: both print `-4*m*xi / (9*delta + 11*xi)`.
- `H_F^(src)`: SymPy prints the combined-fraction form `4*(18*d^2*m + 9*d^2*x + 22*d*m*x + 18*d*x^2 + 11*x^3) / ((9*d + 11*x)*(9*d^2 + 18*d*x + 11*x^2))`; Mathematica prints the partial-fraction form `(4*(xi/(delta + 11*xi/9) + (2*delta*m)/(2*xi^2/9 + (delta+xi)^2))) / 9`. Manual normalization (factor 4/9 against the two denominators) confirms equality.

Engines agree: true.

## Verdict justification

The audit holds up. Attacks tried and failed:

1. **Tautology probe on 25.1**: `n_req` is defined as a specific rational function; `ratio_row1` and `ratio_row2` are then independent rational combinations involving `n_req`. The claim is that both reduce to `ratio_expected = (m*(q-r)+r*xi)/(delta+xi-m*q*(q-r))`. Hand-checking row 1: `xi - m - n_req` with `n_req`'s numerator subtracted yields, after a common denominator, `(m*q + n_req*r) * [(m*(q-r) + r*xi)/(delta+xi-m*q*(q-r))]` — the algebraic structure is non-trivial and the assertion can fail if `n_req` is wrong.
2. **Tautology probe on 25.2**: hand-derived `(1 + q*ratio_expected)^2 / (1 + ratio_expected^2)` and confirmed it equals `(delta + (1+q*r)*xi)^2 / D_qr` only because of the specific choice of `D_qr = (delta+xi-m*q*(q-r))^2 + (m*(q-r)+r*xi)^2`. Not tautological.
3. **R_U=1 limit (25.5)**: hand-checked that `F_src` at `R_U=1` evaluates to `(delta+(1+lam0)*xi)^4 / ((1-xi)*((delta+xi)^2 + lam0*xi^2)^2) = F_flat`. The numerator factor `a1^2 b1^2` becomes `(delta+(1+lam0)*xi)^4` and `dSrc` becomes `(delta+xi)^2 + lam0*xi^2`. Real identity.
4. **First-order coefficients (25.6)**: hand-computed `H_n_expected = -2*lam0*m*xi/(delta+(1+lam0)*xi) = -4*m*xi/(9*delta+11*xi)` with `lam0=2/9`; matches output line 127-129. Hand-derived `H_F_expected` and confirmed it factors to the form Mathematica reports. Real identity.
5. **Series expansion (25.6 endings)**: SymPy's `sp.series(..., eps, 0, 2).removeO()` returns terms up to order 1 in eps; Mathematica's `Series[..., {eps, 0, 1}]` returns terms up to order 1. Both convention-correct.
6. **Symbol-assumption probe**: `delta>0, xi>0` declared in both engines; `m, q, r, t` real (unrestricted sign) in both. No assumption silently kills branches. `R_U > 0` declared in both. No tautology via overstrong assumptions.
7. **Freshness**: sympy output (2026-05-11 12:42) > sympy script (2026-04-01 12:39); mathematica output (2026-05-11 12:49) > mathematica script (2026-05-11 11:56). Both fresh.
8. **Status-only carve-out**: `is_status_only_candidate: False` — both engines required and both present.

Verdict: clean.

## Self-test notes

- **Variable independence**: For `sp.diff(n_src, R_U)` (line 148) and `sp.diff(F_ratio, R_U)` (line 155), checked that both `n_src` and `F_src` explicitly contain `R_U` in their definitions (line 143 and 117-121 respectively). Derivatives are non-trivial; the post-substitution `R_U->1` evaluation produces a non-zero expression, matching the printed output `-4*m*xi/(9*delta+11*xi)`.
- **Parity / branch cuts**: no integrals or unbounded sums in this unit; n/a. R_U=1 substitution does not pinch the denominator of any expression (checked dSrc at R_U=1: `(delta+xi)^2 + lam0*xi^2 > 0`).
- **Trivial-case substitution**: at `R_U=1` (so eps=0), `n_src` reduces to `xi(delta+xi)/(delta+(1+lam0)xi) - m = G_flat - m`, which is exactly the zeroth-order term in `n_linear`. The linear-expansion assertion thus reduces non-trivially in the first-order term to the H_n_expected check, which is independently asserted in A9.
- No directive written (findings_count = 0).
