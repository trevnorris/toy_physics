import math

class HutchisonExperimentConfig:
    def __init__(self, sample_length_mm=10.0, sample_density=2700.0, target_mode="levitation"):
        """
        Initialize the Experiment Calculator with Sample Specs.

        Args:
            sample_length_mm (float): Length of the object (e.g., 10mm).
            sample_density (float): Density in kg/m^3 (Al=2700, Cu=8960).
            target_mode (str): "levitation" (high power) or "softening" (low power).
        """
        # --- Model A Constants (from your papers) ---
        self.c_star = 578.0       # Speed of sound in vacuum enthalpy (m/s)
        self.rho_star = 1.2       # Vacuum density coupling (kg/m^3)
        self.chi = 1.0            # Geometric mode factor
        self.g = 9.81

        # --- Sample Specs ---
        self.L = sample_length_mm / 1000.0  # Convert to meters
        self.rho_m = sample_density
        self.mode = target_mode

        # --- Cavity Specs (Assumptions) ---
        self.Q_cavity = 500.0     # Quality factor of the aluminum test chamber
        self.cavity_harmonic = 5  # Which harmonic the box uses (to keep box size manageable)

    def calculate_resonance(self):
        """Step 1: Calculate the Target Frequency (f0)"""
        # f0 = chi * c_star / (2 * L)
        self.f0 = (self.chi * self.c_star) / (2.0 * self.L)
        self.wavelength = self.c_star / self.f0
        return self.f0

    def calculate_thresholds(self):
        """Step 2: Calculate Required Pressure (pA) for the effect"""
        # Levitation Threshold: Force > Gravity
        # p_lev = 2 * c_star * sqrt(rho_m * rho_star * g / k) where k = pi/L
        k = math.pi / self.L
        p_lev = 2.0 * self.c_star * math.sqrt(self.rho_m * self.rho_star * self.g / k)

        # Softening Threshold: Empirical estimate from model plots (~30% of p_lev)
        p_soft = p_lev * 0.35

        if self.mode == "levitation":
            self.target_pA = p_lev * 1.10 # Add 10% safety margin
        else:
            self.target_pA = p_soft

        return self.target_pA

    def design_cavity(self):
        """Step 3: Design the Chamber Dimensions"""
        # The box must be an integer multiple of half-wavelengths to resonate
        # L_box = n * (lambda / 2)
        half_wave = self.wavelength / 2.0
        self.box_length = self.cavity_harmonic * half_wave

        # Optimal Sample Placement: Node location (L/4 of the fundamental wave)
        # The fundamental wave in the sample is length L.
        # The node for max force gradient is at L/4 from the center of the standing wave.
        self.sample_pos = self.wavelength / 8.0 # Lambda/8 is the max gradient point
        return self.box_length

    def calculate_power_budget(self):
        """Step 4: Calculate Watts required to maintain target Pressure"""
        # 1. Energy Density in the vacuum field: u = p^2 / (2 * rho * c^2)
        u_stored = (self.target_pA**2) / (2.0 * self.rho_star * self.c_star**2)

        # 2. Total Energy in the box volume (assume cubic box)
        volume = self.box_length**3
        total_energy = u_stored * volume

        # 3. Power Loss = Energy * Omega / Q
        omega = 2.0 * math.pi * self.f0
        self.power_watts = (total_energy * omega) / self.Q_cavity
        return self.power_watts

    def generate_report(self):
        self.calculate_resonance()
        self.calculate_thresholds()
        self.design_cavity()
        self.calculate_power_budget()

        print(f"--- HUTCHISON EXPERIMENT SPEC SHEET ---")
        print(f"Target Object: {self.L*1000:.1f}mm {'Aluminum' if self.rho_m==2700 else 'Custom'}")
        print(f"Goal: {self.mode.upper()} Effect")
        print(f"\n1. FREQUENCY LOCK")
        print(f"   Target Frequency (f0): {self.f0/1000:.3f} kHz")
        print(f"   Wavelength (lambda):   {self.wavelength*1000:.1f} mm")

        print(f"\n2. APPARATUS DIMENSIONS")
        print(f"   Cavity Size (Internal): {self.box_length*1000:.1f} mm (Cube)")
        print(f"   Sample Placement (z):   {self.sample_pos*1000:.1f} mm from wall (Max Gradient Node)")
        print(f"   WARNING: Do not place sample in center ({self.box_length*500:.1f} mm). Lift is ZERO there.")

        print(f"\n3. POWER REQUIREMENTS")
        print(f"   Target Vacuum Pressure: {self.target_pA/1000:.2f} kPa")
        print(f"   Min. Amplifier Power:   {self.power_watts:.1f} Watts (RMS)")
        print(f"   Recommended Source:     {math.ceil(self.power_watts * 2)} Watt Amplifier (for headroom)")

        print(f"\n4. TUNING PROTOCOL")
        print(f"   - Step A: Warm up at {self.f0/1000:.3f} kHz at 10% power.")
        print(f"   - Step B: Rapid Sweep +/- 500 Hz to capture phase.")
        print(f"   - Step C: Overshoot Power to {math.ceil(self.power_watts * 1.5)}W, then drop to {math.ceil(self.power_watts)}W.")

# --- RUN THE CALCULATOR ---
# Example 1: Standard 10mm Aluminum Sample for Levitation
experiment = HutchisonExperimentConfig(sample_length_mm=10.0, sample_density=2700, target_mode="levitation")
experiment.generate_report()

print("\n" + "="*40 + "\n")

experiment = HutchisonExperimentConfig(sample_length_mm=10.0, sample_density=2700, target_mode="softening")
experiment.generate_report()

# Example 2: Larger 20mm Sample (easier frequency, harder power)
# experiment_large = HutchisonExperimentConfig(sample_length_mm=20.0, sample_density=2700, target_mode="levitation")
# experiment_large.generate_report()
