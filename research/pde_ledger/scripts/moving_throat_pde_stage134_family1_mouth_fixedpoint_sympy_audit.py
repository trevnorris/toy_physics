
#!/usr/bin/env python3
"""
moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py

SymPy audit for the Family-1 shell + first mixed D/N tube reduction.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

Pi = sp.symbols("Pi", positive=True, real=True)
M_s, M_q = sp.symbols("M_s M_q", real=True)

def S(Pi, kappa):
    return sp.simplify(
        Pi * (kappa * sp.tanh(kappa) + Pi * (sp.exp(-Pi) / sp.cosh(kappa) - 1))
        / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2))
    )

banner("STAGE 134 — FAMILY-1 FIXED-POINT REDUCTION")
kk = sp.symbols("kk", positive=True, real=True)
S_shell = sp.simplify(sp.limit(S(Pi, kk), kk, 0))
S_q = sp.simplify(S(Pi, sp.pi/2))
Pi_eq = sp.simplify(M_s + M_q * S_q)

print("S_shell =", S_shell)
print("S_q(Pi) =")
sp.pprint(S_q)
print("Fixed-point law Pi =")
sp.pprint(Pi_eq)

Pi_star = sp.Float("1.50882951349316")
S_star = sp.N(S_q.subs(Pi, Pi_star), 30)
print("S_q(Pi_star) =", S_star)
print("Canonical gain line: M_s = Pi_star - S_q(Pi_star) M_q")
sp.pprint(sp.N(Pi_star - S_star*M_q, 30))

# --- Substantive checks ---

# Check 1: the kappa -> 0 limit of S equals 1 exactly (static shell channel).
residual_shell = sp.simplify(S_shell - 1)
print("S_shell - 1 =", residual_shell)
assert residual_shell == 0, f"static shell limit failed: residual={residual_shell}"
print("OK: S_shell = 1")

# Check 2: S_q evaluated at three independent numeric Pi values matches
# high-precision targets computed independently via mpmath (see
# redteam/resolutions/batch_IV4_paper_alignment.md M3 for the verification).
# Targets are NOT computed by re-typing S(Pi, pi/2) here; they are literals
# from a separate mpmath run.
S_q_at_half = sp.N(S_q.subs(Pi, sp.Rational(1, 2)), 30)
S_q_at_one  = sp.N(S_q.subs(Pi, 1), 30)
S_q_at_two  = sp.N(S_q.subs(Pi, 2), 30)
target_half = sp.Float("0.608336415687717065435990381419", 30)
target_one  = sp.Float("0.633127670034487546375729566676", 30)
target_two  = sp.Float("0.681366857005321783286541952613", 30)
for name, got, want in [
    ("S_q(1/2)", S_q_at_half, target_half),
    ("S_q(1)",   S_q_at_one,  target_one),
    ("S_q(2)",   S_q_at_two,  target_two),
]:
    diff = abs(sp.N(got - want, 30))
    print(f"{name} = {got}  (target {want}, diff {diff})")
    assert diff < sp.Float("1e-12"), f"{name} mismatch: diff={diff}"
print("OK: S_q matches independent numeric targets at Pi=1/2, 1, 2")

# Check 3: S_q(Pi_*) ≈ 0.658075937605428 (carried from notes; Pi_* = 1.50882951349316
# is numerically located in stage 131, not re-derived here).
S_star_target = sp.Float("0.658075937605428", 30)
diff_star = abs(sp.N(S_star - S_star_target, 30))
assert diff_star < sp.Float("1e-12"), f"S_q(Pi_star) mismatch: diff={diff_star}"
print("OK: S_q(Pi_star) matches notes value 0.658075937605428")

# Note (no in-stage gain-line assertion): the canonical gain line
#   M_s = Pi_* - S_q(Pi_*)*M_q
# is printed symbolically above for the transcript only. The intercept is the
# imported literal Pi_* (owned by stage 131; see stage 131 note) and the slope
# is -S_q(Pi_*), already validated against the notes value in Check 3 above.
# Re-asserting intercept == Pi_* and slope == -S_q(Pi_*) here would be an X-X
# tautology (it compares constants already inserted into the script), so it is
# intentionally omitted. The substantive deliverable that uses this gain line --
# outlet consistency of the gain pair (M_s, M_q) -- is verified downstream at
# Stage 135 (outlet-consistent mouth closure); susceptibility closure at Stage 140.
