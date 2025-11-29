import numpy as np

try:
    import cupy as cp
except ImportError:
    import numpy as cp

from .core import HAS_GPU

class DyonEnsemble:
    """
    Data structure holding the state of N particles (Dyons).
    """
    def __init__(self, N):
        self.N = N
        self.pos = cp.zeros((N, 3), dtype=np.float64)
        self.vel = cp.zeros((N, 3), dtype=np.float64)
        self.acc = cp.zeros((N, 3), dtype=np.float64)
        self.mass = cp.ones(N, dtype=np.float64)      # Sink Q
        self.circ = cp.zeros((N, 3), dtype=np.float64) # Vector Gamma (Spin axis)

    def add_particle(self, idx, pos, vel, mass, circ_vec=(0,0,0)):
        if idx >= self.N: raise ValueError("Index out of bounds")
        self.pos[idx] = cp.array(pos)
        self.vel[idx] = cp.array(vel)
        self.mass[idx] = mass
        self.circ[idx] = cp.array(circ_vec)


class ForceAccumulator:
    """
    Computes accelerations from Fields and Fluid Dynamics.
    """
    def __init__(
        self,
        grid,
        poisson_solver,
        wave_solver=None,
        beta=0.0,
        enable_beta_inertia=False,
    ):
        """
        Parameters
        ----------
        grid : Grid3D
        poisson_solver : PoissonSolverFFT
        wave_solver : WaveEquationSolver or None
        beta : float
            Dimensionless coefficient for position-dependent inertia σ.
            beta = 0.0  -> no inertia correction (pure scalar sector).
            beta = 1.5  -> full 1PN toy model matching GR precession.
        enable_beta_inertia : bool
            If True, apply the β-based inertia correction to particle
            accelerations using the toy model Lagrangian.
        """
        self.grid = grid
        self.poisson = poisson_solver
        self.wave = wave_solver

        self.beta = beta
        self.enable_beta_inertia = enable_beta_inertia

        # Extract c_s from wave solver if available; needed for σ = -β Φ_P / c_s^2
        self.c_s = getattr(wave_solver, "c_s", None) if wave_solver is not None else None

    def deposit_mass(self, ensemble):
        """
        Cloud-in-Cell (trilinear) mass deposition with periodic wrapping.
        Returns a density field (mass / volume).
        """
        N = self.grid.N
        dx = self.grid.dx
        L = self.grid.L

        # Convert to grid coordinates in index space
        grid_coords = (ensemble.pos + (L / 2.0)) / dx
        base_idx = cp.floor(grid_coords).astype(cp.int64)
        frac = grid_coords - base_idx

        fx, fy, fz = frac[:, 0], frac[:, 1], frac[:, 2]
        wx0, wx1 = 1.0 - fx, fx
        wy0, wy1 = 1.0 - fy, fy
        wz0, wz1 = 1.0 - fz, fz

        rho = cp.zeros((N, N, N), dtype=cp.float64)

        offsets = [
            (0, 0, 0, wx0 * wy0 * wz0),
            (0, 0, 1, wx0 * wy0 * wz1),
            (0, 1, 0, wx0 * wy1 * wz0),
            (0, 1, 1, wx0 * wy1 * wz1),
            (1, 0, 0, wx1 * wy0 * wz0),
            (1, 0, 1, wx1 * wy0 * wz1),
            (1, 1, 0, wx1 * wy1 * wz0),
            (1, 1, 1, wx1 * wy1 * wz1),
        ]

        for ox, oy, oz, weight in offsets:
            idx_x = (base_idx[:, 0] + ox) % N
            idx_y = (base_idx[:, 1] + oy) % N
            idx_z = (base_idx[:, 2] + oz) % N

            flat_idx = idx_x * (N * N) + idx_y * N + idx_z
            cp.add.at(rho.ravel(), flat_idx, ensemble.mass * weight)

        rho /= dx**3
        return rho

    def compute_gravity_field(self, ensemble, rho_field, G_sim, phi_wave=None, use_wave=False):
        """
        Compute gravitational acceleration from a single scalar potential Φ.

        Args:
            ensemble: Particle ensemble to update.
            rho_field: Mass density on the grid.
            G_sim: Gravitational constant in simulation units.
            phi_wave: Dynamical scalar potential Φ evolved by WaveEquationSolver.
            use_wave: If True and phi_wave is provided, use it as Φ; otherwise use Φ_P.

        Returns:
            phi_total: Potential used for dynamics.
            phi_P: Static Poisson solution for diagnostics.
        """
        phi_P = self.poisson.solve(rho_field, G_sim)

        if use_wave and (phi_wave is not None):
            phi_total = phi_wave
        else:
            phi_total = phi_P

        if HAS_GPU:
            # Keep gradients on device to avoid host/device copies
            grad_x, grad_y, grad_z = cp.gradient(phi_total, self.grid.dx)
            gx = -grad_x
            gy = -grad_y
            gz = -grad_z
        else:
            grad_x, grad_y, grad_z = np.gradient(phi_total, self.grid.dx)
            gx = -grad_x
            gy = -grad_y
            gz = -grad_z

        ax = self._sample_field(gx, ensemble.pos)
        ay = self._sample_field(gy, ensemble.pos)
        az = self._sample_field(gz, ensemble.pos)

        ensemble.acc[:, 0] += ax
        ensemble.acc[:, 1] += ay
        ensemble.acc[:, 2] += az

        if (
            self.enable_beta_inertia
            and self.beta != 0.0
            and (self.c_s is not None)
        ):
            self._apply_beta_inertia_correction(
                ensemble,
                phi_P=phi_P,
                ax_scalar=ax,
                ay_scalar=ay,
                az_scalar=az,
            )

        return phi_total, phi_P

    def _apply_beta_inertia_correction(self, ensemble, phi_P, ax_scalar, ay_scalar, az_scalar):
        """
        Apply the β-based position-dependent inertia correction implied by
        the toy-model Lagrangian:

            L = 1/2 m (1 + σ(x)) v^2 - m Φ_eff(x),

        with σ(x) = -β Φ_P(x) / c_s^2 and the acceleration law

            a_full = 1/(1+σ) * [ g_eff - (v·∇σ) v + 1/2 v^2 ∇σ ],

        where g_eff is the already-applied scalar acceleration.
        This computes a_corr = a_full - g_eff and adds it to ensemble.acc.
        """
        beta = float(self.beta)
        c_s = float(self.c_s)
        if c_s == 0.0 or beta == 0.0:
            return

        dx = self.grid.dx

        # --- 1. Gradient of Φ_P on the grid ---
        # We only ever use ∇Φ_P to build ∇σ at particle positions.
        grad_phiP_x, grad_phiP_y, grad_phiP_z = cp.gradient(phi_P, dx)

        # --- 2. Sample Φ_P and ∇Φ_P at particle positions ---
        # σ_p = -β Φ_P / c_s^2 at particles
        phiP_p = self._sample_field(phi_P, ensemble.pos)
        sigma_p = -beta * phiP_p / (c_s ** 2)

        gpx_p = self._sample_field(grad_phiP_x, ensemble.pos)
        gpy_p = self._sample_field(grad_phiP_y, ensemble.pos)
        gpz_p = self._sample_field(grad_phiP_z, ensemble.pos)

        # We no longer need the grid gradients; they can be GC'd later.
        # (Optional: for CuPy you can also free the memory pool here.)

        # ∇σ at particles: ∇σ = -β / c_s^2 ∇Φ_P
        factor = -beta / (c_s ** 2)
        gsx = factor * gpx_p
        gsy = factor * gpy_p
        gsz = factor * gpz_p

        # --- 3. Velocities and base scalar acceleration (g_eff) ---
        vx = ensemble.vel[:, 0]
        vy = ensemble.vel[:, 1]
        vz = ensemble.vel[:, 2]

        gx = ax_scalar
        gy = ay_scalar
        gz = az_scalar

        v2 = vx * vx + vy * vy + vz * vz
        vdot_gradsigma = vx * gsx + vy * gsy + vz * gsz

        one_plus_sigma = 1.0 + sigma_p
        inv_one_plus_sigma = 1.0 / one_plus_sigma

        a_full_x = inv_one_plus_sigma * (gx - vdot_gradsigma * vx + 0.5 * v2 * gsx)
        a_full_y = inv_one_plus_sigma * (gy - vdot_gradsigma * vy + 0.5 * v2 * gsy)
        a_full_z = inv_one_plus_sigma * (gz - vdot_gradsigma * vz + 0.5 * v2 * gsz)

        a_corr_x = a_full_x - gx
        a_corr_y = a_full_y - gy
        a_corr_z = a_full_z - gz

        ensemble.acc[:, 0] += a_corr_x
        ensemble.acc[:, 1] += a_corr_y
        ensemble.acc[:, 2] += a_corr_z

    def compute_magnus_force(self, ensemble, fluid_vel_field, rho_fluid=1.0):
        """
        Computes Hydrodynamic Lorentz Force.
        F_magnus = rho * (v_rel x Gamma)
        v_rel = v_particle - v_fluid_local
        """
        # Unpack fluid velocity field (Expect tuple of 3D arrays: vx, vy, vz)
        vfx_grid, vfy_grid, vfz_grid = fluid_vel_field

        # Sample fluid velocity at particle positions
        vfx = self._sample_field(vfx_grid, ensemble.pos)
        vfy = self._sample_field(vfy_grid, ensemble.pos)
        vfz = self._sample_field(vfz_grid, ensemble.pos)

        # Relative Velocity (Particle - Fluid)
        # If particle moves WITH fluid, force is zero.
        # If particle moves AGAINST fluid, force is max.
        # Note: Standard Magnus is v_rel x Gamma.
        vr_x = ensemble.vel[:, 0] - vfx
        vr_y = ensemble.vel[:, 1] - vfy
        vr_z = ensemble.vel[:, 2] - vfz

        # Cross Product: v_rel x Gamma
        # Gamma is ensemble.circ (Nx3)
        # Fx = vy*Gz - vz*Gy
        Fx = rho_fluid * (vr_y * ensemble.circ[:, 2] - vr_z * ensemble.circ[:, 1])
        Fy = rho_fluid * (vr_z * ensemble.circ[:, 0] - vr_x * ensemble.circ[:, 2])
        Fz = rho_fluid * (vr_x * ensemble.circ[:, 1] - vr_y * ensemble.circ[:, 0])

        # Apply to Acceleration (a += F/m)
        ensemble.acc[:, 0] += Fx / ensemble.mass
        ensemble.acc[:, 1] += Fy / ensemble.mass
        ensemble.acc[:, 2] += Fz / ensemble.mass

    def _sample_field(self, field_grid, positions):
        """Trilinear interpolation helper (Cloud-in-Cell sampling)."""
        grid_coords = (positions + (self.grid.L / 2.0)) / self.grid.dx
        coords_transposed = grid_coords.T

        if HAS_GPU:
            import cupyx.scipy.ndimage as ndimage
        else:
            import scipy.ndimage as ndimage

        return ndimage.map_coordinates(field_grid, coords_transposed, order=1, mode='wrap')

    def compute_halo_force(self, ensemble, v0_halo, r_core=1.0):
        """
        Applies the inward pressure gradient from a rotating superfluid halo.
        F = - rho * (v0^2 / r)
        Assumes rho=1 for simplicity.
        """
        # Calculate r for all particles
        pos = ensemble.pos
        r_sq = pos[:, 0]**2 + pos[:, 1]**2 + pos[:, 2]**2
        r = cp.sqrt(r_sq) + 1e-12

        # Smooth the core to avoid singularity at r=0
        # F_mag = v0^2 / (r + r_core)
        f_mag = (v0_halo**2) / (r + r_core)

        # Direction is Inward (-pos / r)
        fx = -f_mag * (pos[:, 0] / r)
        fy = -f_mag * (pos[:, 1] / r)
        fz = -f_mag * (pos[:, 2] / r)

        # Apply Accel
        ensemble.acc[:, 0] += fx / ensemble.mass
        ensemble.acc[:, 1] += fy / ensemble.mass
        ensemble.acc[:, 2] += fz / ensemble.mass


class SymplecticStepper:
    def half_kick(self, ensemble, dt):
        ensemble.vel += 0.5 * ensemble.acc * dt

    def drift(self, ensemble, dt):
        ensemble.pos += ensemble.vel * dt

    def clear_acc(self, ensemble):
        ensemble.acc.fill(0.0)
