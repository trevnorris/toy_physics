from sympy import symbols, pi, sqrt, solve, diff, simplify, besselj, N

# Define symbols
a, L, Pvac, k, x01, rho0, Gamma, Cq = symbols('a L Pvac k x01 rho0 Gamma Cq', real=True, positive=True)

# --- Part 1: Stability Check ---
# Mode Energy: E ~ sqrt((x01/a)^2 + (pi/L)^2)
modeEnergy = k * sqrt((x01/a)**2 + (pi/L)**2)

# Vacuum Work: W = Pvac * Volume
volume = pi * a**2 * L
enthalpy = modeEnergy + Pvac * volume

# Derivatives
dPda = diff(enthalpy, a)
dPdl = diff(enthalpy, L)

# Solve dH/da = 0 and dH/dL = 0 for Pvac
sol_Pvac_a = solve(dPda, Pvac)[0]
sol_Pvac_l = solve(dPdl, Pvac)[0]

# Equate Pvacs to find geometric constraint
constraint = sol_Pvac_a - sol_Pvac_l
# Solve for L in terms of a
sol_L = solve(constraint, L)

print("Possible solutions for L:", sol_L)

# Assume the real positive solution is the one we want.
# Based on the Mathematica output, L/a is a constant.
# Let's extract L/a
L_val = sol_L[0] # Picking the first one, usually the relevant one for simple constraints
ratio_expr = simplify(L_val / a)
print("Symbolic L/a ratio:", ratio_expr)

# Numerical value
# x01 is approx 2.4048
x01_val = 2.4048
ratio_num = ratio_expr.subs(x01, x01_val).evalf()
print("Numerical L/a ratio:", ratio_num)


# --- Part 2: Hierarchy Check ---
# Mass ~ Volume * rho0
# Charge ~ Area * Circulation (Gamma)

# We assume stable geometry, so L = ratio * a
stable_L = ratio_expr * a
mass_stable = rho0 * pi * a**2 * stable_L
# mass_stable simplifies to rho0 * pi * ratio * a^3

charge_stable = rho0 * (pi * a**2) * Gamma

# Force Ratio ~ Q^2 / M^2
force_ratio = (charge_stable**2) / (mass_stable**2)
force_ratio_simplified = simplify(force_ratio)

print("\nForce Ratio Expression:", force_ratio_simplified)

# Check scaling with 'a'
# Q^2 ~ a^4
# M^2 ~ a^6
# Ratio ~ 1/a^2
print("Power of 'a' in Force Ratio:", diff(simplify(force_ratio), a) * a / force_ratio) # This gives the exponent if it's a power law

"""
Possible solutions for L: [sqrt(2)*pi*a/x01]
Symbolic L/a ratio: sqrt(2)*pi/x01
Numerical L/a ratio: 1.84750621180903

Force Ratio Expression: Gamma**2*x01**2/(2*pi**2*a**2)
Power of 'a' in Force Ratio: -2
"""
