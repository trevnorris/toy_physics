---
unit_id: 048
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

# Audit unit 048 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage048_support_compensation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.txt`

## What the script claims to verify

The scripts verify, for a "support compensation theorem" governing a tracking-branch
mixed-load model, four claims: (1) the closed-form for the critical load
`M_crit = lim_{xi -> 1^-} G_tr(xi) = 9(1+delta)/(2R^2+9*delta+9)` and the explicit
factored forms of `dG_tr/dxi` and the gap `M_crit - G_tr`; (2) the analytic
inverse map between the support-enhancement factor
`S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta*eps)` and `zeta`, namely
`zeta_req = (S_req-1)/(1+eps(S_req-2))` (with a parallel formula for
`zeta_crit`), plus the pole `zeta -> 1/eps^-` of `S` diverging to +oo on the
physical branch `0 < eps < 1`; (3) two explicit feasibility margins,
`1/eps - zeta_req` and `zeta_crit - zeta_req`, against their factored
closed-forms; and (4) the implicit-differentiation relation
`dxi/dzeta = M_mix * dS/dzeta / (dG_tr/dxi)` against an explicit closed-form.
Mathematica adds two independent limit checks not present in SymPy: a softening
coefficient `lim_{xi -> 1^-} (1-xi) F_tr` and a pole coefficient
`lim_{zeta -> 1/eps^-} (1/eps - zeta) S = (1-eps)/eps^2`.

## Assertion inventory

| #  | Script      | Line | Form | Anchored to claim? |
|----|-------------|------|------|--------------------|
| A1 | sympy       | 52-56  | `simplify(dG_dxi - 9*(...)/(2R^2 xi+9 delta+9 xi)^2) == 0` | yes |
| A2 | sympy       | 57     | `simplify(F_tr.subs(xi,0) - 1) == 0` | yes |
| A3 | sympy       | 61-66  | `simplify(M_crit - G_tr - 9(1-xi)(...)/((...)(...))) == 0` | yes |
| A4 | sympy       | 73     | `simplify(S.subs(zeta,0) - 1) == 0` | yes |
| A5 | sympy       | 74     | `diff(S,zeta) - (1-eps)/(1-zeta eps)^2 == 0` | yes |
| A6 | sympy       | 81-82  | `assert limit_phys == oo` for physical-branch S divergence | yes |
| A7 | sympy       | 87     | `simplify(S.subs(zeta,zeta_req) - S_req) == 0` | yes |
| A8 | sympy       | 88     | `simplify(S.subs(zeta,zeta_crit) - S_crit) == 0` | yes |
| A9 | sympy       | 99-102 | `simplify((1/eps - zeta_req) - (1-eps)/(eps(1+eps(S_req-2)))) == 0` | yes |
| A10| sympy       | 103-106 | `simplify((zeta_crit - zeta_req) - (S_crit-S_req)(1-eps)/((...)(...))) == 0` | yes |
| A11| sympy       | 112-117 | `simplify(dxi_dzeta - M_mix(1-eps)(...)^2/((1-zeta eps)^2 * 9 (...))) == 0` | yes |
| B1 | mathematica | 49-52 | `expectZero(dG_tr/dxi formula)` | yes |
| B2 | mathematica | 53    | `expectZero(F_tr(xi=0)-1)` | yes |
| B3 | mathematica | 64    | `expectZero((1-xi) F_tr softening coefficient)` | yes |
| B4 | mathematica | 68-72 | `expectZero(M_crit - G_tr formula)` | yes |
| B5 | mathematica | 76    | `expectZero(S(zeta=0)-1)` | yes |
| B6 | mathematica | 77    | `expectZero(dS/dzeta - (1-eps)/(1-zeta eps)^2)` | yes |
| B7 | mathematica | 84    | `expectZero(pole coefficient formula = (1-eps)/eps^2)` | yes |
| B8 | mathematica | 86-88 | `Solve[sEnhance==sReq,zeta]` and zetaReq compared back through S | yes |
| B9 | mathematica | 93    | `expectZero(inverse map S(zeta_req)-S_req)` | yes |
| B10| mathematica | 94    | `expectZero(inverse map S(zeta_crit)-S_crit)` | yes |
| B11| mathematica | 100   | `expectZero(pole margin formula)` | yes |
| B12| mathematica | 101-104 | `expectZero(branch margin formula)` | yes |
| B13| mathematica | 108-112 | `expectZero(dxi_phys/dzeta formula)` | yes |

Every row is non-tautological: the LHS in each assertion is *computed*
(via `diff`, `limit`, `factor`, `simplify`, `Solve`) and then matched against
a hardcoded closed-form RHS. If the algebra were wrong, the residual would
not simplify to 0.

## Findings

None.

## Independent-derivation check (Mathematica)

Mathematica is NOT a transliteration of SymPy.

- Mathematica uses `Solve[sEnhance == sReq, zeta, Reals]` (line 86) to *derive*
  `zeta_req` from the inverse-map equation, then compares against the same
  closed form SymPy hardcodes. SymPy never solves; it asserts the inverse.
  This is genuinely independent algebra.
- Mathematica adds a softening coefficient check via
  `Limit[(1-xi)*fTr, xi -> 1, "FromBelow"]` (lines 55-64), compared against an
  expanded explicit limit value. SymPy has no analogue — it only checks
  `limit(F_tr, xi -> 1^-) = oo` qualitatively.
- Mathematica adds a pole coefficient check via
  `Limit[(1/eps - zeta) * sEnhance, zeta -> 1/eps, "FromBelow"]` against
  `(1-eps)/eps^2` (lines 79-84). SymPy has no analogue.

The two engines share the same defining expressions for `G_tr`, `F_tr`, and
`S` (these are physical inputs, not derivations), but the verifications they
perform diverge in non-trivial ways. The Mathematica side does strictly more
independent work than SymPy.

## Engine cross-check

Where both engines test the same identity, both report residual `0` (PASS).
Specifically:

| Identity                          | SymPy residual | Mathematica residual |
|-----------------------------------|----------------|----------------------|
| dG_tr/dxi formula                 | 0              | 0 (PASS)             |
| F_tr(xi=0) - 1                    | 0              | 0 (PASS)             |
| M_crit - G_tr formula             | 0              | 0 (PASS)             |
| S(zeta=0) - 1                     | 0              | 0 (PASS)             |
| dS/dzeta - (1-eps)/(1-zeta eps)^2 | 0              | 0 (PASS)             |
| Inverse map S(zeta_req) - S_req   | 0              | 0 (PASS)             |
| Inverse map S(zeta_crit) - S_crit | 0              | 0 (PASS)             |
| pole margin formula               | 0              | 0 (PASS)             |
| branch margin formula             | 0              | 0 (PASS)             |
| dxi_phys/dzeta formula            | 0              | 0 (PASS)             |

No disagreement.

Both engines also report identical explicit forms for `M_crit`
(`9(1+delta)/(9+9*delta+2R^2)`), `dG_tr/dxi`, `S(zeta;eps)`,
`zeta_req`, `zeta_crit`, and `dxi_phys/dzeta` (up to trivial rearrangement
like `(xi-1)` vs `-(1-xi)`).

## Verdict justification

Verdict: clean.

Attacks tried and failed:

1. *Tautology probe.* Each `expect_zero` LHS is reached by a non-trivial chain
   (`sp.diff` then `sp.factor`, or `sp.simplify(sp.limit(...))`); the RHS is a
   distinct hardcoded closed-form. Replacing the RHS with anything but the
   true algebraic value would fail `simplify(... - ...) == 0`. No row of the
   inventory is algebraically guaranteed by construction.
2. *Hardcoded-result probe.* The RHS formulas (e.g. for `dG_tr/dxi`, gap,
   margins, `dxi/dzeta`) are explicit closed-forms. They are not numeric
   constants pulled from upstream; they are symbolic shapes the assertions
   are claiming. Each can be derived by hand (I rederived `dG_tr/dxi`'s
   numerator from the quotient rule and matched it to
   `9(2R^2 xi^2 + 9 delta^2 + 18 delta xi + 9 xi^2)`).
3. *Symbol-domain probe.* SymPy declares all symbols `positive=True, real=True`.
   The `S` pole-divergence assertion needs `0 < eps < 1`; this is routed via
   the `nu` substitution (`eps_phys = 1/(1+nu)` with `nu > 0`) before taking
   the limit and checking it equals `+oo` — so the domain constraint is
   enforced where it matters. Mathematica imposes `0 < eps < 1` globally in
   `$Assumptions`. No domain abuse.
4. *Branch-coverage probe.* The scripts only consider the tracking branch
   `G_tr`. The docstring/comments restrict scope to "tracking branch" — they
   never claim coverage of other branches. So this is in-scope only for the
   tracking branch, and the tracking-branch checks fully exercise the
   stated claim.
5. *Transliteration probe.* Mathematica performs an independent `Solve` for
   `zeta_req` and two additional limit-coefficient checks that SymPy lacks.
   Same physical inputs, divergent verification strategies.
6. *Freshness probe.* SymPy output mtime (May 11 12:48) equals SymPy script
   mtime; Mathematica output mtime (May 11 12:51) is later than script mtime
   (May 11 11:56). Captured `EXIT_CODE: 0` and PASS lines. Outputs are fresh
   and consistent with the present scripts.

The scripts hold up. No findings.

## Self-test notes

- *Variable independence*: `D[gTr, xi]` and `sp.diff(G_tr, xi)` are taken on
  expressions that explicitly contain `xi` (numerator `9 xi (delta+xi)`,
  denominator `9 delta + (9+2R^2) xi`); derivative is non-trivially nonzero,
  as the script's output confirms. Likewise `D[sEnhance, zeta]` /
  `sp.diff(S, zeta)` on `S(zeta;eps)` which depends on `zeta`. No
  identically-zero-derivative trap.
- *Parity / symmetry*: No integrals over unbounded domains appear; the entire
  unit is algebraic/limit-based. The pole and softening limits use a one-sided
  `Direction -> FromBelow` / `dir="-"` to handle the simple pole at `zeta = 1/eps`
  and at `xi = 1` correctly. The signs of the residual factorizations
  (`-(xi-1)` vs `(1-xi)`) are consistent across both engines and identities.
- *Trivial-case substitution*: At `zeta = 0`, `S = 1` (checked, A4/B5). At
  `xi = 0`, `F_tr = 1` (checked, A2/B2). The hardcoded `dG_tr/dxi` formula at
  `xi = 0` reduces to `9 * 9 delta^2 / (9 delta)^2 = 1`, matching the numerator's
  algebraic form at xi=0. Self-consistent.
- *Path / engine presence*: Both engines present at the expected `scripts/`
  and `mathematica/` paths; no missing-script directive needed.
