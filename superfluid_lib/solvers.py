import numpy as np
import warnings

try:
    import cupy as cp
except ImportError:
    import numpy as cp

class PoissonSolverFFT:
    """
    Spectral solver for the static Poisson equation:
        ∇^2 Φ_P = 4π G_sim ρ
    Uses FFT to invert the Laplacian (Φ_P is the Newtonian reference potential).
    """
    def __init__(self, grid):
        """
        Args:
            grid: An instance of core.Grid3D.
        """
        self.grid = grid

        # Precompute -1/k^2 for the inverse Laplacian
        KX, KY, KZ, K_sq = grid.get_k_space()

        # Avoid division by zero at k=0 (the 'DC' component).
        # We set the inverse to 0 at k=0, which effectively subtracts
        # the infinite background potential (Jeans swindle).
        with np.errstate(divide='ignore', invalid='ignore'):
            self.inv_k_sq = -1.0 / K_sq

        # Fix the singularity at index [0,0,0]
        self.inv_k_sq[0, 0, 0] = 0.0

    def solve(self, rho, G_sim):
        """
        Returns the potential Phi given density rho.
        Phi_k = -4 * pi * G * rho_k / k^2
        """
        # 1. Forward FFT
        rho_k = cp.fft.fftn(rho)

        # 2. Spectral inversion
        # Factor 4*pi*G comes from the Poisson equation definition in info.md
        phi_k = (4 * cp.pi * G_sim) * rho_k * self.inv_k_sq

        # 3. Inverse FFT (return real part to strip numerical noise)
        phi = cp.fft.ifftn(phi_k).real
        return phi


class WaveEquationSolver:
    """
    Evolves the scalar gravitational potential Φ:
        d^2Φ/dt^2 = c_s^2 (∇^2 Φ - 4π G_sim ρ) - γ dΦ/dt
    In the static limit this solution equals the Poisson potential Φ_P.
    """
    def __init__(self, grid, c_s, gamma=0.0):
        self.grid = grid
        self.c_s = c_s      # Sound speed (Simulation Units)
        self.gamma = gamma  # Damping factor

        # Precompute Laplacian operator in k-space (-k^2)
        KX, KY, KZ, K_sq = grid.get_k_space()
        self.laplacian_k = -K_sq
        self._warned_cfl = False

        if self.c_s == 0:
            self.cfl_dt_limit = np.inf
        else:
            self.cfl_dt_limit = grid.dx / (self.c_s * np.sqrt(3.0))

    def max_stable_dt(self, safety=1.0):
        """
        Returns the CFL-stable timestep (optionally scaled by safety).
        """
        return safety * self.cfl_dt_limit

    def step(self, phi_curr, phi_prev, rho, dt, G_sim):
        """
        Performs one Verlet time-step.
        Returns: phi_next
        """
        if dt > self.cfl_dt_limit and not self._warned_cfl:
            warnings.warn(
                f"WaveEquationSolver dt={dt:.3g} exceeds CFL limit {self.cfl_dt_limit:.3g} "
                f"(dx={self.grid.dx:.3g}, c_s={self.c_s:.3g}). Results may be unstable.",
                RuntimeWarning
            )
            self._warned_cfl = True

        # 1. Compute spatial Laplacian of current field via FFT
        phi_k = cp.fft.fftn(phi_curr)
        lap_phi_k = phi_k * self.laplacian_k
        lap_phi = cp.fft.ifftn(lap_phi_k).real

        # 2. Compute Source Term (Force drive)
        # The wave tries to relax toward the Poisson solution, so source is 4*pi*G*rho
        source = 4 * cp.pi * G_sim * rho

        # 3. Wave Acceleration
        # a = c_s^2 * (Laplacian - Source) - Damping
        # Note: We damp velocity estimated as (phi_curr - phi_prev) / dt
        vel_est = (phi_curr - phi_prev) / dt
        accel = (self.c_s**2) * (lap_phi - source) - (self.gamma * vel_est)

        # 4. Verlet Integration
        # phi_next = 2*phi_curr - phi_prev + a * dt^2
        phi_next = 2 * phi_curr - phi_prev + accel * (dt**2)

        return phi_next
