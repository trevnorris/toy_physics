import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ==========================================
# 1. PHYSICAL CONSTANTS (Dimensionless)
# ==========================================
G = 1.0
M = 1.0
c0 = 1.0  # Speed of light at infinity
mu = G * M

# ==========================================
# 2. THE SUPERFLUID ENVIRONMENT GENERATOR
# ==========================================

def getEnvironment(rVec, mode):
    """
    Returns local sound speed (cs) and flow velocity (vFlow)
    based on the active Physics Mode.
    """
    r = np.linalg.norm(rVec)

    # --- DEFAULT VALUES ---
    cs = c0
    vFlow = np.zeros(2)

    # --- PHYSICS COMPONENT 1: BULK TENSION (Space/Refraction) ---
    # Vacuum is Stiff (n=3).
    # Creates a refractive index n ~ 1 + GM/(r c^2)
    # This contributes exactly 2GM/bc^2 bending (Half of GR).
    if mode == 'refractionOnly' or mode == 'combined' or mode == 'splitModel':
        cs = c0 * (1.0 - mu / (r * c0**2))

    # --- PHYSICS COMPONENT 2: INFLOW DRAG (Time/Ether Wind) ---
    # The fluid flows into the sink.
    if mode == 'flowEscape' or mode == 'combined':
        # CASE A: Full Escape Velocity (v^2 = 2GM/r)
        # This usually corresponds to full g_00 time dilation.
        # In this Hamiltonian formulation, it might be too strong.
        vMag = np.sqrt(2.0 * mu / r)
        vFlow = -vMag * (rVec / r)

    elif mode == 'flowOrbital' or mode == 'splitModel':
        # CASE B: Orbital Velocity (v^2 = GM/r)
        # We hypothesize the bulk flow is slower than the brane pinch.
        # This should reduce drag bending by half.
        vMag = np.sqrt(1.0 * mu / r)
        vFlow = -vMag * (rVec / r)

    return cs, vFlow

# ==========================================
# 3. HAMILTONIAN RAY EQUATIONS
# ==========================================

def rayEquations(t, state, mode):
    """
    Equations of motion for a phonon in moving fluid.
    Derived from Hamiltonian H = cs*|k| + v . k
    """
    pos = state[:2]
    k = state[2:] # Momentum vector

    kMag = np.linalg.norm(k)
    kHat = k / kMag
    rVal = np.linalg.norm(pos)

    # Get local properties
    cs, vFlow = getEnvironment(pos, mode)

    # 1. GROUP VELOCITY (dx/dt)
    dxdt = cs * kHat + vFlow

    # 2. MOMENTUM EVOLUTION (dk/dt = -Grad H)
    # Numerical Gradient for robustness
    eps = 1e-5 * rVal

    # Gradient in X
    posXp = pos + np.array([eps, 0])
    posXm = pos - np.array([eps, 0])
    csXp, vXp = getEnvironment(posXp, mode)
    csXm, vXm = getEnvironment(posXm, mode)
    HXp = csXp * kMag + np.dot(vXp, k)
    HXm = csXm * kMag + np.dot(vXm, k)
    dHdx = (HXp - HXm) / (2 * eps)

    # Gradient in Y
    posYp = pos + np.array([0, eps])
    posYm = pos - np.array([0, eps])
    csYp, vYp = getEnvironment(posYp, mode)
    csYm, vYm = getEnvironment(posYm, mode)
    HYp = csYp * kMag + np.dot(vYp, k)
    HYm = csYm * kMag + np.dot(vYm, k)
    dHdy = (HYp - HYm) / (2 * eps)

    dkdt = -np.array([dHdx, dHdy])

    return np.concatenate([dxdt, dkdt])

# ==========================================
# 4. SIMULATION RUNNER
# ==========================================

def runSimulation(impactB, mode):
    # Start far away
    startX = -1000.0 * impactB
    startY = impactB
    startK = np.array([1.0, 0.0]) # Moving +x

    initState = np.array([startX, startY, startK[0], startK[1]])
    tSpan = [0, 2500.0 * impactB]

    # Pass 'mode' to the solver using lambda
    sol = solve_ivp(lambda t, y: rayEquations(t, y, mode),
                    tSpan, initState, rtol=1e-9, atol=1e-9)

    # Get final angle
    kFinal = sol.y[2:, -1]
    thetaFinal = np.arctan2(kFinal[1], kFinal[0])

    # Calculate Ratio to GR (4GM/bc^2)
    thetaGR = 4.0 * mu / (impactB * c0**2)
    ratio = np.abs(thetaFinal) / thetaGR

    return ratio

# ==========================================
# 5. EXECUTE ISOLATION TESTS
# ==========================================

bTest = 1000.0 # Use large b for cleanest 1PN approximation

print(f"{'TEST CASE':<25} | {'PHYSICS INPUT':<30} | {'PREDICTED RATIO':<15} | {'ACTUAL RATIO':<15}")
print("-" * 95)

# --- TEST 1: REFRACTION ONLY ---
# Checks the Stiff Vacuum (n=3) math.
# Should give 2GM (Ratio 0.5).
ratio1 = runSimulation(bTest, 'refractionOnly')
print(f"{'1. Refraction Only':<25} | {'n=3 Vacuum, No Flow':<30} | {'0.50':<15} | {ratio1:<15.4f}")

# --- TEST 2: FLOW ONLY (ESCAPE) ---
# Checks the standard 'Ether Drag' assumption (v = sqrt(2GM/r)).
# Previous runs suggested this was too strong (Ratio ~1.0).
ratio2 = runSimulation(bTest, 'flowEscape')
print(f"{'2. Flow (Escape Vel)':<25} | {'v = sqrt(2GM/r)':<30} | {'1.00           ':<15} | {ratio2:<15.4f}")

# --- TEST 3: THE SPLIT MODEL ---
# Checks the proposed fix:
# Space = Stiff Vacuum (n=3) -> 0.5
# Time  = Orbital Flow (v = sqrt(GM/r)) -> 0.5 (Hypothesized)
# Total -> 1.0
ratio3 = runSimulation(bTest, 'splitModel')
print(f"{'3. Split Model':<25} | {'n=3 Vacuum + v=sqrt(GM/r)':<30} | {'1.00 (Target)':<15} | {ratio3:<15.4f}")

print("-" * 95)

# ==========================================
# 6. FINAL VERDICT
# ==========================================
if abs(ratio3 - 1.0) < 0.05:
    print("\n[SUCCESS] The 'Split Model' matches General Relativity!")
    print("Conclusion: 1PN Lensing is resolved by a Stiff Vacuum (n=3) and Orbital Bulk Flow.")
else:
    print("\n[INCONCLUSIVE] Tuning required. Check Flow contribution.")
