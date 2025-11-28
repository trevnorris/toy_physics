import numpy as np

# Try to import cupy for GPU acceleration; fallback to numpy if not found
try:
    import cupy as cp
    HAS_GPU = True
    print(" [System] CuPy detected. Running in GPU mode.")
except ImportError:
    import numpy as cp
    HAS_GPU = False
    print(" [System] CuPy not found. Running in CPU mode (NumPy).")

class UnitScaler:
    """
    Manages the translation between Physical Units (SI) and Simulation Units.

    Based on the scaling relations defined in the info.md Section 7.1:
    X = x_tilde * L0
    t = t_tilde * T0
    m = m_tilde * M0
    """

    G_PHYS_SI = 6.67430e-11  # m^3 kg^-1 s^-2

    def __init__(self, L0, T0, M0):
        """
        Initialize the scaler with characteristic scales.

        Args:
            L0 (float): Characteristic length in meters (e.g., 1 AU or 1 kpc).
            T0 (float): Characteristic time in seconds (e.g., 1 year).
            M0 (float): Characteristic mass in kg (e.g., 1 Solar Mass).
        """
        self.L0 = L0
        self.T0 = T0
        self.M0 = M0

        # Calculate Dimensionless G (G_sim)
        # Derived from: G_phys * M0 / L0^2 * (T0^2 / L0) = G_sim
        # Simplified: G_sim = (G_phys * M0 * T0^2) / L0^3
        self.G_sim = (self.G_PHYS_SI * self.M0 * (self.T0**2)) / (self.L0**3)

        print(f" [Units] Initialized Scaler.")
        print(f"         L0={L0:.2e}m, T0={T0:.2e}s, M0={M0:.2e}kg")
        print(f"         derived G_sim = {self.G_sim:.6f} (dimensionless)")

    # --- Converters: Physical -> Sim ---
    def to_sim_len(self, val_m): return val_m / self.L0
    def to_sim_time(self, val_s): return val_s / self.T0
    def to_sim_mass(self, val_kg): return val_kg / self.M0
    def to_sim_vel(self, val_ms): return val_ms * (self.T0 / self.L0)

    # --- Converters: Sim -> Physical ---
    def to_phys_len(self, val_sim): return val_sim * self.L0
    def to_phys_time(self, val_sim): return val_sim * self.T0
    def to_phys_mass(self, val_sim): return val_sim * self.M0
    def to_phys_vel(self, val_sim): return val_sim * (self.L0 / self.T0)


class Grid3D:
    """
    A 3D Cartesian grid wrapper using CuPy/NumPy.
    Handles mesh generation and frequency-space (k-space) setup for FFT solvers.
    """

    def __init__(self, N, L_sim):
        """
        Args:
            N (int): Resolution per side (total grid points = N^3).
            L_sim (float): Side length of the box in Simulation Units.
        """
        self.N = N
        self.L = L_sim
        self.dx = L_sim / N

        # Spatial coordinates (centered at 0)
        # Range: [-L/2, L/2)
        lin = cp.linspace(-L_sim/2, L_sim/2, N, endpoint=False)

        # Create Meshgrid (Index ordering: 'ij' usually preferred for matrix math)
        self.X, self.Y, self.Z = cp.meshgrid(lin, lin, lin, indexing='ij')

        # Pre-allocate complex arrays for wave functions to save memory alloc time later
        self.shape = (N, N, N)

        print(f" [Grid3D] Initialized {N}^3 grid. Box Size={L_sim:.2f}. dx={self.dx:.4f}")

    def get_k_space(self):
        """
        Generates wave vectors (kx, ky, kz) and k^2 for spectral solvers.
        """
        # Frequency components associated with the grid
        k_freq = cp.fft.fftfreq(self.N, d=self.dx) * 2 * cp.pi

        # Meshgrid for k-space
        KX, KY, KZ = cp.meshgrid(k_freq, k_freq, k_freq, indexing='ij')

        # k^2 magnitude (Squared Laplacian operator in Fourier space)
        # Added epsilon to avoid division by zero at k=0
        K_sq = KX**2 + KY**2 + KZ**2

        return KX, KY, KZ, K_sq


class Grid4D:
    """
    STUB: Placeholder for future 4D bulk operations.
    Keeps the API dimension-agnostic.
    """
    def __init__(self, N, L_sim):
        print(" [Grid4D] Warning: 4D Grid is currently a stub.")
        self.ndim = 4
        self.N = N
        self.L = L_sim
