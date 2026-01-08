from sympy import symbols, Function, diff, simplify
from sympy.vector import CoordSys3D, curl, divergence, gradient

# Initialize Coordinate System
N = CoordSys3D('N')
x, y, z = N.x, N.y, N.z
t = symbols('t')
c = symbols('c')  # sound speed / light speed effective

# Define Potentials A (Vector) and phi (Scalar) as functions of space and time
# We use generic functions for components
Ax = Function('Ax')(x, y, z, t)
Ay = Function('Ay')(x, y, z, t)
Az = Function('Az')(x, y, z, t)
phi = Function('phi')(x, y, z, t)

# Construct Vector Potential A
A = Ax * N.i + Ay * N.j + Az * N.k

# --- Step 1: Define EM Fields from Fluid Potentials ---
# B = Curl A
B = curl(A)

# E = -Grad phi - dA/dt
E = -gradient(phi) - diff(A, t)

# --- Step 2: Calculate LHS of Ampere-Maxwell Law ---
# LHS = Curl B - (1/c^2) * dE/dt
# We want to see what this simplifies to.
LHS = curl(B) - (1/c**2) * diff(E, t)

# --- Step 3: Apply Vector Identity and Gauge Condition ---
# Vector Identity: curl(curl(A)) = grad(div(A)) - laplacian(A)
# We can compute curl(curl(A)) directly or rely on sympy to handle it.
# Let's check the divergence of A and E.

# We impose the Lorenz Gauge condition for the effective theory:
# Div A + (1/c^2) d(phi)/dt = 0  =>  Div A = -(1/c^2) d(phi)/dt
div_A = divergence(A)
lorenz_gauge_condition = -(1/c**2) * diff(phi, t)

# Substitute Divergence of A with the Gauge term in the expression?
# Sympy might not simplify vector calculus identities automatically 
# without explicit component expansion, so let's look at the components.

# Calculate the x-component of the LHS
LHS_x = LHS.dot(N.i)

# We want to show LHS_x corresponds to -Box(Ax) = - (Laplacian Ax - 1/c^2 d^2Ax/dt^2)
# Box operator on a scalar f: grad^2 f - 1/c^2 d^2f/dt^2
# But wait, the standard identity is:
# Curl(Curl A) = Grad(Div A) - Laplacian A
# So Curl B = Grad(Div A) - Laplacian A
# LHS = Grad(Div A) - Laplacian A - 1/c^2 d/dt (-Grad phi - dA/dt)
#     = Grad(Div A) - Laplacian A + Grad(1/c^2 dphi/dt) + 1/c^2 d^2A/dt^2
#     = Grad( Div A + 1/c^2 dphi/dt ) - (Laplacian A - 1/c^2 d^2A/dt^2)

# If Lorenz Gauge holds (Div A + 1/c^2 dphi/dt = 0), the first term vanishes.
# Then LHS = - (Laplacian A - 1/c^2 d^2A/dt^2) = - Box A.
# If the wave equation holds (Box A = -mu0 J), then LHS = mu0 J.
# This recovers Ampere's Law: Curl B - 1/c^2 dE/dt = mu0 J.

# Let's verify this algebra with Sympy to be sure.
# We will construct the term "Remainder" = LHS - (-Laplacian A + 1/c^2 d^2A/dt^2)
# And see if it equals Grad( Div A + 1/c^2 dphi/dt )

# Laplacian of a Vector in Cartesian coordinates is just Laplacian of components
laplacian_A = (diff(Ax, x, 2) + diff(Ax, y, 2) + diff(Ax, z, 2)) * N.i + \
              (diff(Ay, x, 2) + diff(Ay, y, 2) + diff(Ay, z, 2)) * N.j + \
              (diff(Az, x, 2) + diff(Az, y, 2) + diff(Az, z, 2)) * N.k

# Box A (D'Alembertian)
box_A = laplacian_A - (1/c**2) * diff(A, t, 2)

# The "Source" J is defined by the wave equation Box A = -J (ignoring mu0 for now)
# So we expect LHS to be equal to -Box A + Gauge_Terms

difference = LHS - (-box_A)

# We expect this difference to be Grad( Div A + 1/c^2 dphi/dt )
gauge_term = divergence(A) + (1/c**2) * diff(phi, t)
grad_gauge = gradient(gauge_term)

# Check if difference - grad_gauge is zero
check = simplify((difference - grad_gauge).dot(N.i)) # Check x component

print("Verification of Ampere-Maxwell derivation:")
print(f"Is the identity (Curl B - 1/c^2 dE/dt) == (-Box A + Grad(Gauge)) valid? {check == 0}")

if check == 0:
    print("\nCONCLUSION:")
    print("If the fluid satisfies the Wave Equation (Box A = -mu0 J) and the Lorenz Gauge condition,")
    print("then the Ampere-Maxwell Law is automatically satisfied.")
    print("Therefore, we do NOT need a new hydrodynamic derivation for the displacement current;")
    print("it is a mathematical consequence of the acoustic wave equation derived in Papers I-III.")
else:
    print("Verification failed. Manual inspection needed.")

"""
Verification of Ampere-Maxwell derivation:
Is the identity (Curl B - 1/c^2 dE/dt) == (-Box A + Grad(Gauge)) valid? True

CONCLUSION:
If the fluid satisfies the Wave Equation (Box A = -mu0 J) and the Lorenz Gauge condition,
then the Ampere-Maxwell Law is automatically satisfied.
Therefore, we do NOT need a new hydrodynamic derivation for the displacement current;
it is a mathematical consequence of the acoustic wave equation derived in Papers I-III.
"""
