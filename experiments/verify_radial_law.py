import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

from superfluid_lib.core import UnitScaler, Grid3D
from superfluid_lib.solvers import PoissonSolverFFT

try:
    import cupy as cp
    to_cpu = cp.asnumpy
except ImportError:
    import numpy as cp
    to_cpu = lambda x: x

def run_experiment():
    print("=== EXPERIMENT 1.1: Radial Law Verification (Near-Field) ===")

    # 1. Setup Simulation Environment
    # N=128, L=100
    N = 128
    L_sim = 100.0
    grid = Grid3D(N, L_sim)
    scaler = UnitScaler(L0=1.0, T0=1.0, M0=1.0)

    # 2. Create Mass Source
    # TUNED: Sigma = 2.5 * dx.
    # This is the limit of "safe" smoothness. It allows us to measure
    # the field much closer to the center before PBC artifacts take over.
    sigma = 2.5 * grid.dx
    r_sq = grid.X**2 + grid.Y**2 + grid.Z**2
    rho = cp.exp(-r_sq / (2 * sigma**2))

    # Normalize
    current_mass_integral = cp.sum(rho) * (grid.dx**3)
    desired_mass = 1000.0
    rho *= (desired_mass / current_mass_integral)

    print(f"initialized: Grid=128^3, Mass=1000.0, Sigma={sigma:.2f}")

    # 3. Solve Potential
    solver = PoissonSolverFFT(grid)
    phi = solver.solve(rho, scaler.G_sim)

    # 4. Compute Force via Exact Spectral Derivative
    print("computing: Spectral Gradient...")
    KX, KY, KZ, K_sq = grid.get_k_space()
    phi_k = cp.fft.fftn(phi)

    # F = -Grad(Phi) -> -ik * Phi_k
    fx = cp.fft.ifftn(-1j * KX * phi_k).real
    fy = cp.fft.ifftn(-1j * KY * phi_k).real
    fz = cp.fft.ifftn(-1j * KZ * phi_k).real

    f_mag = cp.sqrt(fx**2 + fy**2 + fz**2)

    # 5. Data Analysis
    r_flat = to_cpu(cp.sqrt(r_sq).flatten())
    f_flat = to_cpu(f_mag.flatten())

    # FILTER: The "Newtonian Sweet Spot"
    # Min: 4.5*sigma (Just outside the Gaussian core)
    # Max: L/6 (Far enough from periodic boundaries)
    # For L=100, this is roughly r in [8.7, 16.6]
    r_min = 4.5 * sigma
    r_max = L_sim / 6.0

    mask = (r_flat > r_min) & (r_flat < r_max)
    r_data = r_flat[mask]
    f_data = f_flat[mask]

    print(f"analyzing: {len(r_data)} grid points in zone [{r_min:.1f}, {r_max:.1f}]")

    # 6. Regression
    log_r = np.log10(r_data)
    log_f = np.log10(f_data)
    slope, intercept = np.polyfit(log_r, log_f, 1)

    print(f"RESULT: Slope = {slope:.5f} (Target: -2.00000)")

    # 7. Visualization
    plt.figure(figsize=(10, 6))
    plt.scatter(r_data[::10], f_data[::10], s=1, alpha=0.2, color='blue', label='Sim Data')
    fit_y = (10**intercept) * (r_data**slope)
    plt.plot(r_data, fit_y, 'r--', linewidth=2, label=f'Fit: Slope={slope:.4f}')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Distance r')
    plt.ylabel('Force |F|')
    plt.title(f'Phase 1.1: Radial Law Check\nSlope: {slope:.5f}')
    plt.legend()
    plt.grid(True, which='both', alpha=0.5)

    output_file = os.path.join(os.path.dirname(__file__), 'radial_law_nearfield.png')
    plt.savefig(output_file)
    print(f"saved: {output_file}")

if __name__ == "__main__":
    run_experiment()
