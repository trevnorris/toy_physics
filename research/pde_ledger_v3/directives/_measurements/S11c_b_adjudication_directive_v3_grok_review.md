# Grok leg measurements — S11c-b adjudication directive v3 (decision-list, round 3)

Artifact: `directives/S11c_b_adjudication_build_directive.md` (v3).
Sources read before judging: comparator, v1 reconcile, PY/WL engines, SHARED_PHYSICS, `~/.s11_build/S11c_b_reconcile_run.out`.

## Bridge A θ² coefficient alignment (no Mathematica)

From sources:
- WL:472-479: `uniformCoefficients` entry `bRho anchoredWidth` with `uniformFactors` `1/2` on `thetaWave^2`
  ⇒ WL θ² energy coefficient = `(1/2)*bRho*anchoredWidth`
- WL:1625-1629: homogeneous block coefficient `bRho WZero` with factor `1/2`
  ⇒ `(1/2)*bRho*WZero`
- PY:1135: catalogue `(btheta**2, B_rho_3 * W_bg / (2 * W0))`
  ⇒ PY θ² energy coefficient = `B_rho_3*W_bg/(2*W0)`
- spec:102: `B_ρ⁽³⁾ ≡ B_ρ W₀` inside `½ B_ρ⁽³⁾ θ²`

Equating WL and PY with `anchoredWidth = W_bg`:
```
(1/2)*bRho*W_bg = B_rho_3*W_bg/(2*W0)  ⇒  B_rho_3 = bRho*W_0  ⇒  bRho ↦ B_rho_3/W_0
```
Uniform limit `W_bg = WZero = W_0`:
```
WL (1/2)*bRho*W_0 = PY B_rho_3/2  ⇒  same identity
```

## PY selected-index recomputation (literal)

```
$ cd research/pde_ledger_v3/scripts && python3 - <<'PY'
import S11c_b_brane_operator_sympy_audit as A
cands = A.enumerate_new_candidates(tuple(A.bg))
exprs = tuple(e for _, e in cands)
sigs = A.basis_euler_signatures(exprs, A.basis_fields)
selected, omitted = A.quotient_independent_indices(exprs, sigs)
print(tuple(i+1 for i in selected))
print(tuple(i+1 for i in omitted))
PY
(1, 4, 5, 6, 7, 9, 10, 13)
(2, 3, 8, 11, 12, 14, 15)
```

## Boolean / TextAtom vs sp.Expr (literal)

```
Equivalent is Expr? False
Not is Expr? False
Eq is Expr? False
TextAtom is Expr? False
Symbol is Expr? True
```

## Jet spelling of the mandated rename (literal)

```
canon_jet_name('theta_d1') -> 'theta_d1'
canon_jet_name('grad_theta_1') -> 'grad_theta_1'   # no-op; trailing _1 is not _dN
naive _dN decode: theta_d1 -> (theta, (1,)); grad_theta_1 -> (grad_theta_1, ())
```
PY:205-207 defines `grad_theta_{i}` (not `theta_d{i}`).

## DivGrad rename absence (literal)

WL names `gammaWidthDivGradTheta/Ew`, `gammaModulusDivGradTheta/Ew` — none present in `H.WL_TO_PY_RENAME`.
Mapped gammas target only `{01,04,05,06,09,13}` per source.

## Run residual kinds (v1 reconcile; no Bridge A)

`SLAB_OPERATOR_TERM_ORIGINS` by OBJECT: ADVECTIVE algebraic; BULK_ENERGY STRUCTURE; KINETIC STRUCTURE;
FACE_FLUX/ACCUMULATION/FACE/FLUX TextAtom (unjoined).
`ADMISSIBILITY_RESIDUAL`: STRUCTURAL_ATOM or STRUCTURE only (no bare algebraic Expr).
`MU_THETA_OPERATOR` residual contains both `bRho` and `B_rho_3`, plus `gamma_s11cb_w_bg_07` and unmapped `gammaWidthDivGrad*`.
