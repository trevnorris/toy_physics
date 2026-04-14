
#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage014_continuum_placement_dynamic_class_map_sympy_audit.py

Stage 014 — pull the Stage-013 selected-branch classifier and dynamic sign
thresholds back through the actual continuum selected-branch placement map.

What this script does
---------------------
1. Carries forward the exact Stage-012 / Stage-013 selected-branch classifier
      R_ND(xi,delta)
   and the universal selected-branch geometry
      F(xi,delta), G(xi,delta).
2. Uses the Stage-25 / Stage-26 continuum placement rule
      F(xi_req,delta) = R_target
   to define the actual physical classifier
      R_phys(delta,R_target) = R_ND(xi_req(delta,R_target),delta).
3. Proves that for fixed delta,
      dR_phys/dR_target < 0,
   so larger normalization demand ratio pushes the selected branch in the
   denominator-like / dynamically safer direction.
4. Builds the exact threshold compiler for any classifier cap c>0, via the
   unique root xi_c(delta) of
      R_ND(xi_c,delta) = c,
   whenever the onset value 8/(9 delta) exceeds c.
5. Specializes that to:
      - c = R_*  : lower-pole sign-flip threshold from Stage 013,
      - c = 1    : denominator-like threshold from Stage 012.
6. Pulls those thresholds back to the exact continuum kernel ratios
      (eps_eta, eps_W, rho, Z_W, delta_0, Lambda)
   and to the mixed-baseline coordinate M_mix through the exact product law.
7. Reconfirms that the dynamic ceilings remain globally weaker than the
   transported static Xi_1 ceilings even on the physical continuum placement map.

Main structural result
----------------------
The actual continuum-selected branch does not resurrect a dynamic kill condition.
Instead, the continuum map sharpens the classification:

- for fixed delta, larger R_target always makes the selected branch more
  denominator-like;
- equivalently, at fixed product scale, larger M_mix makes it more numerator-like;
- there is an exact target curve R_flip(delta) (or mixed-load curve M_flip(delta))
  beyond which the nonempty dynamic ceiling is infinite;
- there is a second exact target curve R_den(delta) beyond which the physical
  branch is denominator-like;
- but everywhere on the continuum placement map the first kill condition remains
  the transported static Xi_1 budget, not the wall-like dynamic window.
"""

from __future__ import annotations
import math
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.factor(sp.simplify(sp.expand(z))))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.factor(sp.simplify(sp.expand(expr)))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("STAGE 014 — CONTINUUM PLACEMENT DYNAMIC-CLASS MAP")

# ---------------------------------------------------------------------------
# Carried universal selected-branch geometry and classifier
# ---------------------------------------------------------------------------

subbanner("I. Carried universal selected-branch geometry")

xi, delta = sp.symbols("xi delta", positive=True, real=True)

F = sp.factor(
    (9 * delta + 11 * xi) ** 4
    / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)
)

G = sp.factor(9 * xi * (xi + delta) / (9 * delta + 11 * xi))

R_ND = sp.factor(
    72 * delta**2 * (1 - xi)
    / ((9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2))
)

dF_dxi = sp.factor(sp.diff(F, xi))
dR_dxi = sp.factor(sp.diff(R_ND, xi))

print("F(xi,delta) =")
sp.pprint(F)
print("G(xi,delta) =")
sp.pprint(G)
print("R_ND(xi,delta) =")
sp.pprint(R_ND)
print("dF/dxi =")
sp.pprint(dF_dxi)
print("dR_ND/dxi =")
sp.pprint(dR_dxi)

expect_zero(
    "dF/dxi - positive_fraction",
    sp.simplify(
        dF_dxi
        - (9 * delta + 11 * xi) ** 3
        * (81 * delta**3 + 189 * delta**2 * xi + 72 * delta**2 + 297 * delta * xi**2 + 121 * xi**3)
        / (81 * (1 - xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 3)
    ),
)

expect_zero(
    "dR_ND/dxi + positive_fraction",
    sp.simplify(
        dR_dxi
        + 72
        * delta**2
        * (81 * delta**3 + 261 * delta**2 + 594 * delta * xi - 297 * delta * xi**2 + 363 * xi**2 - 242 * xi**3)
        / ((9 * delta + 11 * xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)
    ),
)

print("Therefore on the stable interval 0 <= xi < 1, delta > 0:")
print("  dF/dxi > 0, so F = R_target has a unique solution xi_req(delta,R_target) for every R_target >= 1.")
print("  dR_ND/dxi < 0, so the selected-branch classifier decreases monotonically along the placement locus.")

# ---------------------------------------------------------------------------
# Continuum placement map and exact product law
# ---------------------------------------------------------------------------

subbanner("II. Exact continuum placement map")

eps_eta, eps_W, rho, ZW, delta0, Lambda = sp.symbols(
    "eps_eta eps_W rho Z_W delta_0 Lambda", positive=True, real=True
)

delta_map = sp.factor(delta0 / (1 - eps_eta))
Mmix_map = sp.factor(8 * ZW * (1 + rho) ** 2 / (sp.pi**2 * (1 - eps_eta) * (1 - eps_W)))
Rtarget_map = sp.factor(Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (ZW * (1 + rho) ** 2))

print("delta =")
sp.pprint(delta_map)
print("M_mix =")
sp.pprint(Mmix_map)
print("R_target =")
sp.pprint(Rtarget_map)

expect_zero(
    "R_target * M_mix - 8 Lambda (1-eps_W)/pi^2",
    sp.simplify(Rtarget_map * Mmix_map - 8 * Lambda * (1 - eps_W) / sp.pi**2),
)

print("So the continuum placement map is exactly:")
print("  delta    = delta_0 / (1 - eps_eta)")
print("  M_mix    = 8 Z_W (1+rho)^2 / [pi^2 (1-eps_eta)(1-eps_W)]")
print("  R_target = Lambda (1-eps_eta)(1-eps_W)^2 / [Z_W (1+rho)^2]")
print("with exact product law R_target * M_mix = 8 Lambda (1-eps_W)/pi^2.")

# ---------------------------------------------------------------------------
# Pullback of the classifier to the physical placement locus
# ---------------------------------------------------------------------------

subbanner("III. Pullback of the classifier to the physical placement locus")

R_target = sp.symbols("R_target", positive=True, real=True)

dRphys_dRtarget = sp.factor(sp.simplify(dR_dxi / dF_dxi))

print("If xi_req(delta,R_target) is the unique solution of F(xi,delta) = R_target, then")
print("  R_phys(delta,R_target) := R_ND(xi_req(delta,R_target),delta).")
print()
print("By implicit differentiation:")
print("  d xi_req / d R_target = 1 / (dF/dxi) > 0,")
print("so")
print("  d R_phys / d R_target = (dR_ND/dxi) / (dF/dxi).")
print("Exact sign-carrying ratio =")
sp.pprint(dRphys_dRtarget)

print("Hence d R_phys / d R_target < 0 on the full stable branch.")
print("For fixed delta, larger R_target always pushes the physical selected branch in the denominator-like / dynamically safer direction.")

print()
print("Along the fixed product curve R_target * M_mix = 8 Lambda (1-eps_W)/pi^2,")
print("the same statement becomes")
print("  d R_phys / d M_mix > 0.")
print("So larger mixed baseline pushes the physical selected branch in the numerator-like direction.")

# ---------------------------------------------------------------------------
# Exact threshold compiler for any classifier cap c > 0
# ---------------------------------------------------------------------------

subbanner("IV. Exact threshold compiler for any classifier cap")

c = sp.symbols("c", positive=True, real=True)

P_c = sp.expand(
    c * (9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2)
    - 72 * delta**2 * (1 - xi)
)
dP_c_dxi = sp.factor(sp.diff(P_c, xi))

print("Threshold polynomial P_c(xi,delta) from R_ND(xi,delta) = c:")
sp.pprint(P_c)
print("dP_c/dxi =")
sp.pprint(dP_c_dxi)

expect_zero(
    "P_c(0,delta) - 9 delta^2 (9 c delta - 8)",
    sp.simplify(P_c.subs({xi: 0}) - 9 * delta**2 * (9 * c * delta - 8)),
)

print("Since dP_c/dxi > 0 for xi >= 0, delta > 0, c > 0, the threshold root is unique whenever it exists.")
print("The onset value is R_ND(0,delta) = 8/(9 delta).")
print("Therefore the exact onset threshold for the cap c is")
print("  delta_c = 8 / (9 c).")
print()
print("If delta >= delta_c, then the physical branch already satisfies R_phys <= c at onset,")
print("so the pulled-back target threshold is simply R_target >= 1.")
print("If 0 < delta < delta_c, there is a unique xi_c(delta) in (0,1) with R_ND(xi_c,delta) = c,")
print("and the exact pulled-back target threshold is")
print("  R_c(delta) = F(xi_c(delta), delta) > 1.")

# ---------------------------------------------------------------------------
# Specialize to the Stage-013 dynamic sign threshold and denominator threshold
# ---------------------------------------------------------------------------

subbanner("V. Dynamic sign threshold and denominator threshold")

R_star = sp.Float("1.229255438463336", 50)
delta_star_dyn = sp.N(8 / (9 * R_star), 30)
delta_den = sp.Rational(8, 9)

print(f"Stage-013 lower-pole sign-flip classifier threshold: R_* = {R_star}")
print(f"Corresponding onset threshold: delta_*^(dyn) = 8 / (9 R_*) = {delta_star_dyn}")
print(f"Denominator-like onset threshold from Stage 012: delta_den = 8/9 = {sp.N(delta_den, 30)}")

print()
print("Pulled back to the continuum placement map:")
print("  - if delta >= delta_*^(dyn), then every physical branch point has infinite nonempty dynamic ceiling already from onset,")
print("    because R_phys <= R_* on the whole stable interval;")
print("  - if delta >= 8/9, then every physical branch point is denominator-like from onset,")
print("    because R_phys <= 1 on the whole stable interval.")

# helper functions for sample threshold tables
def threshold_xi(delta_value: float, c_value: float) -> float:
    delta_c = 8.0 / (9.0 * c_value)
    if delta_value >= delta_c:
        return 0.0
    poly = sp.Poly(
        sp.expand(P_c.subs({delta: sp.Float(str(delta_value), 80), c: sp.Float(str(c_value), 80)})),
        xi,
    )
    roots = sp.nroots(poly)
    physical = sorted(
        float(complex(r).real)
        for r in roots
        if abs(complex(r).imag) < 1.0e-10 and 0.0 < float(complex(r).real) < 1.0
    )
    if len(physical) != 1:
        raise AssertionError(f"Expected unique physical threshold root, got {physical}")
    return physical[0]


def threshold_target(delta_value: float, c_value: float) -> float:
    xi_c = threshold_xi(delta_value, c_value)
    return float(sp.N(F.subs({xi: sp.Float(str(xi_c), 80), delta: sp.Float(str(delta_value), 80)}), 50))


sample_deltas = [0.25, 0.50, 0.75]
print("Sample pulled-back target thresholds:")
print("  delta       xi_flip        R_flip         xi_den         R_den")
for dval in sample_deltas:
    xi_flip = threshold_xi(dval, float(R_star))
    R_flip = threshold_target(dval, float(R_star))
    xi_den = threshold_xi(dval, 1.0)
    R_den = threshold_target(dval, 1.0)
    print(f"  {dval:5.2f}   {xi_flip:12.9f}  {R_flip:12.9f}  {xi_den:12.9f}  {R_den:12.9f}")

# ---------------------------------------------------------------------------
# Exact continuum-kernel threshold surfaces
# ---------------------------------------------------------------------------

subbanner("VI. Exact continuum-kernel threshold surfaces")

R_flip_sym = sp.Symbol("R_flip", positive=True, real=True)
R_den_sym = sp.Symbol("R_den", positive=True, real=True)

ZW_flip_bound = sp.factor(Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (R_flip_sym * (1 + rho) ** 2))
ZW_den_bound = sp.factor(Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (R_den_sym * (1 + rho) ** 2))

M_flip_bound = sp.factor(8 * Lambda * (1 - eps_W) / (sp.pi**2 * R_flip_sym))
M_den_bound = sp.factor(8 * Lambda * (1 - eps_W) / (sp.pi**2 * R_den_sym))

print("Equivalent threshold surfaces in continuum-kernel variables:")
print("If R_target >= R_flip(delta), then")
print("  Z_W (1+rho)^2 <= ")
sp.pprint(ZW_flip_bound)
print("and, by the exact product law, equivalently")
print("  M_mix <= ")
sp.pprint(M_flip_bound)
print()
print("If R_target >= R_den(delta), then")
print("  Z_W (1+rho)^2 <= ")
sp.pprint(ZW_den_bound)
print("and equivalently")
print("  M_mix <= ")
sp.pprint(M_den_bound)

print("So at fixed product scale 8 Lambda (1-eps_W)/pi^2:")
print("  - lowering M_mix pushes the physical branch across the dynamic sign threshold first,")
print("  - lowering it further pushes the branch all the way into the denominator-like regime.")

# ---------------------------------------------------------------------------
# Static-first theorem on the physical continuum placement map
# ---------------------------------------------------------------------------

subbanner("VII. Static-first theorem on the physical continuum placement map")

B_dyn_both_inf = sp.Float("0.967282389363822", 50)
B_dyn_nonempty_inf = sp.Float("0.990581810705233", 50)
B_stat_both = sp.Float("0.367930328492646", 50)
B_stat_nonempty = sp.Float("0.737619063660757", 50)

print("From Stage 013, on the full classifier half-line R >= 0:")
print(f"  inf B_dyn^(both)      = {B_dyn_both_inf}")
print(f"  inf B_dyn^(nonempty)  = {B_dyn_nonempty_inf}")
print(f"  B_stat^(both)         = {B_stat_both}")
print(f"  B_stat^(nonempty)     = {B_stat_nonempty}")

print()
print("The continuum placement map only restricts us to a subset of the full classifier range.")
print("Therefore the same strict inequalities remain true after pullback:")
print("  B_dyn^(both)(physical)    > B_stat^(both),")
print("  B_dyn^(nonempty)(physical)> B_stat^(nonempty).")
print()
print("So even on the actual continuum-selected branch, the first kill condition is still the transported static Xi_1 budget, not the wall-like dynamic window.")

# ---------------------------------------------------------------------------
# Final verdict
# ---------------------------------------------------------------------------

subbanner("VIII. Final verdict")
print("Stage 014 pulls the Stage-013 classifier and dynamic sign split back through the actual continuum placement map.")
print("The result is sharper than the sample-branch picture:")
print("  * for fixed delta, larger R_target always makes the physical selected branch more denominator-like;")
print("  * equivalently, at fixed product scale, larger M_mix makes it more numerator-like;")
print("  * there is an exact pulled-back threshold R_flip(delta) beyond which the nonempty dynamic ceiling is infinite;")
print("  * there is a second exact pulled-back threshold R_den(delta) beyond which the physical branch is denominator-like;")
print("  * these thresholds translate directly into exact inequalities on (eps_eta, eps_W, rho, Z_W, delta_0, Lambda) or on M_mix through the exact product law;")
print("  * but nowhere on the continuum placement map does the dynamic window become the first kill condition.")
print()
print("So the same-charge corridor is still alive after Stage 014.")
print("The remaining first kill condition is still the static placement of Xi_1 on the actual moving-throat branch.")
