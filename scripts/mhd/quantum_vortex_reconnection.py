import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

"""
3D GPE SIMULATION: QUIET START & DEALIASING
-------------------------------------------
Fixes applied:
1. Imaginary Time Relaxation (damps initial shock).
2. Phase Pinning (prevents vortex annihilation during relaxation).
3. 2/3 Rule Dealiasing (prevents spectral leakage).
4. Energy Audit (Helmholtz Decomposition).
"""

# --- GPU SETUP ---
try:
    import cupy as cp
    xp = cp
    use_gpu = True
    print("[System] CuPy detected. GPU Mode Engaged.")
except ImportError:
    import numpy as cp
    xp = np
    use_gpu = False
    print("[System] CuPy not found. CPU Mode (Slower).")

# --- CONFIGURATION ---
N = 128              # Grid Size 
L = 16.0            # Domain Size
DT = 0.05           # Time step
STEPS_RELAX = 100   # Imaginary time steps (The "Quiet Start")
STEPS_REAL = 400    # Real time steps (The Experiment)
PLOT_EVERY = 20     

# --- GRID ---
x = xp.linspace(-L/2, L/2, N)
y = xp.linspace(-L/2, L/2, N)
z = xp.linspace(-L/2, L/2, N)
X, Y, Z = xp.meshgrid(x, y, z, indexing='ij')

dx = x[1] - x[0]
k = xp.fft.fftfreq(N, d=dx/(2*xp.pi)) * 2*xp.pi
KX, KY, KZ = xp.meshgrid(k, k, k, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2

# --- DEALIASING MASK (2/3 Rule) ---
k_max = k.max()
# Zero out top 1/3 of frequencies to stop spectral leakage
dealias_mask = (xp.abs(KX) < 2/3 * k_max) & \
               (xp.abs(KY) < 2/3 * k_max) & \
               (xp.abs(KZ) < 2/3 * k_max)

# --- INITIAL CONDITION: ORTHOGONAL TUBES ---
def vortex_tube_phase(X_grid, Y_grid, x0, y0):
    # Just the phase winding (for pinning)
    return xp.arctan2(Y_grid - y0, X_grid - x0)

def vortex_tube_density(X_grid, Y_grid, x0, y0):
    # Initial guess for density core
    r = xp.sqrt((X_grid - x0)**2 + (Y_grid - y0)**2)
    return (r**2) / (r**2 + 1.0) # Simple Padé

# Tube 1 (Along X)
phi1 = vortex_tube_phase(Y, Z, 0.0, -1.2)
rho1 = vortex_tube_density(Y, Z, 0.0, -1.2)

# Tube 2 (Along Y)
phi2 = vortex_tube_phase(X, Z, 0.0, 1.2)
rho2 = vortex_tube_density(X, Z, 0.0, 1.2)

# Combined State
psi_phase = phi1 + phi2
psi_rho = rho1 * rho2
psi = xp.sqrt(psi_rho) * xp.exp(1j * psi_phase)

# Save the analytic phase for pinning
initial_phase = xp.angle(psi)

# --- HELMHOLTZ DECOMPOSITION (ENERGY METRIC) ---
def get_compressible_energy(psi_curr):
    # Mass current j = Im( conj(psi) * grad_psi )
    psi_k = xp.fft.fftn(psi_curr)
    grad_x = xp.fft.ifftn(1j * KX * psi_k)
    grad_y = xp.fft.ifftn(1j * KY * psi_k)
    grad_z = xp.fft.ifftn(1j * KZ * psi_k)
    
    jx = xp.imag(xp.conj(psi_curr) * grad_x)
    jy = xp.imag(xp.conj(psi_curr) * grad_y)
    jz = xp.imag(xp.conj(psi_curr) * grad_z)
    
    # Project onto k (Compressible part)
    jk_x = xp.fft.fftn(jx)
    jk_y = xp.fft.fftn(jy)
    jk_z = xp.fft.fftn(jz)
    
    k_dot_j = KX*jk_x + KY*jk_y + KZ*jk_z
    
    # Avoid div/0
    K2_safe = K2.copy()
    K2_safe[0,0,0] = 1.0
    
    jc_x = KX * (k_dot_j / K2_safe)
    jc_y = KY * (k_dot_j / K2_safe)
    jc_z = KZ * (k_dot_j / K2_safe)
    
    # Zero DC
    jc_x[0,0,0] = 0; jc_y[0,0,0] = 0; jc_z[0,0,0] = 0
    
    # Energy = sum |j_c|^2
    E_c = 0.5 * xp.sum(xp.abs(jc_x)**2 + xp.abs(jc_y)**2 + xp.abs(jc_z)**2)
    return float(E_c) / (N**3) # Normalize by grid size

# --- VISUALIZATION SETUP ---
history_t = []
history_Ec = []

plt.ion()
fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(211, projection='3d')
ax2 = fig.add_subplot(212)

print("--- PHASE 1: QUIET START (Relaxation) ---")
print("Damping initial phonon shock while pinning vortices...")

# --- PHASE 1: IMAGINARY TIME RELAXATION ---
for step in range(STEPS_RELAX):
    # Imaginary time propagator: exp(-H * dt)
    # 1. Potential Step (Real space)
    rho = xp.abs(psi)**2
    # Decay operator for relaxation: exp(-V * dt)
    psi *= xp.exp(-rho * DT)
    
    # 2. Pinning (Force Phase back to analytic)
    # This keeps the vortex lines where we want them, but lets density profile smooth out
    # psi = |psi_relaxed| * exp(i * initial_phase)
    current_amp = xp.abs(psi)
    psi = current_amp * xp.exp(1j * initial_phase)
    
    # 3. Kinetic Step (Fourier space)
    # Decay operator: exp(-0.5 * k^2 * dt)
    psi_k = xp.fft.fftn(psi)
    psi_k *= xp.exp(-0.5 * K2 * DT)
    psi_k *= dealias_mask # Apply 2/3 Rule
    psi = xp.fft.ifftn(psi_k)
    
    # Renormalize to background density = 1
    # (Simple trick for relaxation stability)
    norm_factor = 1.0 / xp.max(xp.abs(psi))
    psi *= norm_factor

print("Relaxation Complete. System is 'Quiet'.")
print("--- PHASE 2: UNITARY EVOLUTION (The Experiment) ---")

# --- PHASE 2: REAL TIME EVOLUTION ---
for step in range(STEPS_REAL):
    
    # 1. Potential (Unitary: exp(-i * V * dt))
    rho = xp.abs(psi)**2
    psi *= xp.exp(-1j * rho * DT)
    
    # 2. Kinetic (Unitary: exp(-i * K * dt))
    psi_k = xp.fft.fftn(psi)
    psi_k *= xp.exp(-0.5j * K2 * DT)
    psi_k *= dealias_mask # Apply 2/3 Rule
    psi = xp.fft.ifftn(psi_k)
    
    if step % PLOT_EVERY == 0:
        # Measure Energy
        Ec = get_compressible_energy(psi)
        history_t.append(step * DT)
        history_Ec.append(Ec)
        
        # Plot
        if use_gpu:
            d_plot = rho.get()
            X_c, Y_c, Z_c = X.get(), Y.get(), Z.get()
        else:
            d_plot = rho
            X_c, Y_c, Z_c = X, Y, Z
            
        ax1.clear()
        mask = d_plot < 0.1
        stride = 2
        ax1.scatter(X_c[mask][::stride], Y_c[mask][::stride], Z_c[mask][::stride], c='k', s=1, alpha=0.5)
        ax1.set_xlim(-L/2, L/2); ax1.set_ylim(-L/2, L/2); ax1.set_zlim(-L/2, L/2)
        ax1.set_title(f"Step {step}: Topology")
        
        ax2.clear()
        ax2.plot(history_t, history_Ec, 'r-', linewidth=2)
        ax2.set_xlabel("Time")
        ax2.set_ylabel("Compressible Energy (Sound)")
        ax2.set_title("Acoustic Emission (Look for the Jump)")
        ax2.grid(True)
        
        plt.pause(0.001)

plt.ioff()
plt.show()
