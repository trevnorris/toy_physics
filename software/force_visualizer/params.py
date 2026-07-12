"""Single parameter and provenance block for all four visualized sectors.

Cross-sector status is taken from ``pathA_40_cone_lock.md`` (``delta_r=2``:
``lambda_gamma`` and ``c_E`` are independent calibration inputs) and
``pathA_41_ng5_second_medium_drift.md`` (``rho_B0``, ``chi_c`` and ``C_hu``
are active free-unreduced quantities).
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from enum import Enum
from math import sqrt
from typing import Dict, Iterable


class Provenance(str, Enum):
    """Required status vocabulary for numerical values shown by the scenes."""

    DERIVED_FORM = "DERIVED-FORM"
    CALIBRATED_MAGNITUDE = "CALIBRATED-MAGNITUDE"
    FREE_UNREDUCED = "FREE-UNREDUCED"


@dataclass(frozen=True)
class ParameterInfo:
    """Human-readable provenance attached to a value in :class:`ModelParameters`."""

    status: Provenance
    source: str


@dataclass(frozen=True)
class ModelParameters:
    """Internally consistent hand-picked values used by every scene.

    The values are dimensionless visualization units.  Their *forms* follow
    the cited research, while dimensional magnitudes remain calibrated or
    unreduced exactly as recorded below.
    """

    # Shared speed / EM family.
    c_s: float = 10.0
    lambda_gamma: float = 1.0
    rho_br: float = 1.0
    mu_R: float = 100.0
    c_E: float = 10.0
    rho_B0: float = 1.0
    chi_c: float = 4.0
    C_hu: float = 0.10

    # Effective magnitudes.
    G: float = 1.0
    Q_E: float = 2.0
    aT: float = 1.0
    # CALIBRATED: freely chosen sim-deferred amplitude.  Its resulting
    # scalar/transverse ratio is a display calibration, not a model prediction.
    aL: float = 0.022360679774997897
    N_u: float = 1.0

    # Throat/partner geometry used by charge and scalar departures.
    ell: float = 0.80
    mouth_half_width_b: float = 0.40
    M_h: float = 1.0
    yukawa_fraction: float = 0.20

    # Characterized gravity-return and light-longitudinal display inputs.
    epsilon0: float = 0.20
    epsilon1: float = 0.20
    radiation_length_a: float = 0.50
    J_phase: float = 0.50
    kappa_phase: float = 1.0
    longitudinal_display_fraction: float = 0.22

    # The 2.5PN paper calls the native normalization genuinely blocked.
    # One means "show the standard Burke--Thorne benchmark" only.
    radiation_reaction_benchmark_scale: float = 1.0

    @property
    def c_gamma(self) -> float:
        """Light speed, satisfying both calibrated cone locks."""

        return self.lambda_gamma * self.c_s

    @property
    def B_eff(self) -> float:
        """Derived density modulus ``rho_B0**2 / chi_c``."""

        return self.rho_B0**2 / self.chi_c

    @property
    def N0(self) -> float:
        """Derived electric wall zero-mode norm ``8/(3 ell)``."""

        return 8.0 / (3.0 * self.ell)

    @property
    def yukawa_mass(self) -> float:
        """Derived wall-shape partner gap ``sqrt(3)/ell``."""

        return sqrt(3.0) / self.ell

    def validate(self) -> None:
        """Reject parameter sets that violate the shared derived identities."""

        positive = (
            self.c_s,
            self.lambda_gamma,
            self.rho_br,
            self.mu_R,
            self.c_E,
            self.rho_B0,
            self.chi_c,
            self.G,
            self.Q_E,
            self.ell,
            self.mouth_half_width_b,
            self.M_h,
            self.kappa_phase,
        )
        if any(value <= 0.0 for value in positive):
            raise ValueError("positive model parameters must be strictly positive")
        if not 0.0 <= self.yukawa_fraction:
            raise ValueError("yukawa_fraction must be non-negative")
        if not 0.0 <= self.epsilon0 or not 0.0 <= self.epsilon1:
            raise ValueError("return partition strengths must be non-negative")
        if abs(self.c_gamma**2 - self.mu_R / self.rho_br) > 1e-12:
            raise ValueError("shared EM block must satisfy c_gamma^2=mu_R/rho_br")
        if self.C_hu**2 >= self.B_eff * self.M_h * self.c_E**2:
            raise ValueError("scalar block violates B_eff*K_h-C_hu^2>0")


PARAMETER_INFO: Dict[str, ParameterInfo] = {
    "c_s": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "ledger_stage005: derived EOS form; background-state magnitude chosen"),
    "lambda_gamma": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_40: independent cone-lock calibration, delta_r=2"),
    "c_gamma": ParameterInfo(Provenance.DERIVED_FORM, "pathA_36 transverse: c_gamma^2=mu_R/rho_br; magnitude follows calibrated block"),
    "rho_br": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_41: Route-A sim-deferred brane inertia"),
    "mu_R": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_41: Route-A sim-deferred shear modulus"),
    "c_E": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_40: independent electric-cone lock; default set equal to c_gamma"),
    "rho_B0": ParameterInfo(Provenance.FREE_UNREDUCED, "pathA_41 active_irreducible compression parameter"),
    "chi_c": ParameterInfo(Provenance.FREE_UNREDUCED, "pathA_41 active_irreducible susceptibility"),
    "C_hu": ParameterInfo(Provenance.FREE_UNREDUCED, "pathA_41 active_irreducible embedding overlap"),
    "B_eff": ParameterInfo(Provenance.DERIVED_FORM, "pathA_36: B_eff=rho_B0^2/chi_c"),
    "G": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_29 fixes 1/r^2 form only; effective Newton constant is free"),
    "Q_E": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_38 calibrated electric source anchor"),
    "aT": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_39 sim-deferred moving-source transverse amplitude"),
    "aL": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_39 freely chosen sim-deferred moving-source longitudinal amplitude; channel ratio is not predicted"),
    "N_u": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_39 keeps f_u normalization symbolic; display value chosen"),
    "ell": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_41 calibrated throat geometry input"),
    "mouth_half_width_b": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_41 calibrated throat geometry input"),
    "M_h": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_41 calibrated h-sector geometry input"),
    "N0": ParameterInfo(Provenance.DERIVED_FORM, "pathA_38 goldstone.N0_norm=8/(3 ell)"),
    "yukawa_mass": ParameterInfo(Provenance.DERIVED_FORM, "pathA_38 wall-shape gap m=sqrt(3)/ell"),
    "yukawa_fraction": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_38 derives Yukawa form/gap but not a universal main-source residue"),
    "epsilon0": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_29 return continuity partition; residual epsilon0/(1+epsilon0)"),
    "epsilon1": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_29 dipole return continuity partition"),
    "radiation_length_a": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_29 low-frequency residual length scale"),
    "J_phase": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_36 Josephson magnitude not fixed by frozen definitions"),
    "kappa_phase": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "pathA_36 conventional phase stiffness magnitude"),
    "longitudinal_display_fraction": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "visual amplitude only; pathA_36 fixes DOF and pole form, not scene amplitude"),
    "radiation_reaction_benchmark_scale": ParameterInfo(Provenance.CALIBRATED_MAGNITUDE, "4d_2_5pn: GENUINE_BLOCKED native normalization; 1 selects Burke-Thorne benchmark"),
}


DEFAULT_PARAMS = ModelParameters()
DEFAULT_PARAMS.validate()


def value_of(params: ModelParameters, name: str) -> float:
    """Return a stored or derived parameter by its public provenance name."""

    return float(getattr(params, name))


def labeled_value(params: ModelParameters, name: str, precision: int = 4) -> str:
    """Format a number with its mandatory status label."""

    info = PARAMETER_INFO[name]
    return f"[{info.status.value}] {name}={value_of(params, name):.{precision}g}"


def parameter_names() -> Iterable[str]:
    """Names physically stored in the shared dataclass (for audit tests)."""

    return (field.name for field in fields(ModelParameters))
