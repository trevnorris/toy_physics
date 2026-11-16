"""Measure the brane reaction to translation for the isolated target-rate preparations."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from lateral_mouth_acceleration import Translation
from lateral_reservoir_field import LateralSlab
from lateral_source_rate import check, verified_packet


def extend_archived_field(slab, layers, normal):
    """Recover the full interior with the archived mouth derivative held fixed."""
    field = np.broadcast_to(layers[0], slab.eigenvalues.shape).copy()
    field[-4:] = layers
    for _ in range(500):
        potential = (abs(field)**8-1)/4
        real = slab.elliptic(2*(field.real-1)-potential*field.real, normal.real, 2)
        imag = slab.elliptic(-potential*field.imag, normal.imag, 0)
        updated = 1 + real + 1j*imag
        delta = np.max(abs(updated-field))
        if not np.isfinite(delta):
            raise ArithmeticError("Nonfinite interior extension")
        if delta < 2e-13:
            mismatch = np.max(abs(updated[-4:]-layers))
            if mismatch > 2e-11:
                raise ArithmeticError(f"Interior misses the archived mouth layers: {mismatch}")
            return updated, mismatch
        field = updated
    raise ArithmeticError("Interior extension did not converge")


class CalibratedTranslation(Translation):
    """Translate the common rate-0.25 intake shape across a rate-specific base field."""

    def __init__(self, slab, base, current, rate, density_ref):
        self.slab = slab
        self.f = base
        self.face = base[-1]
        self.amp = abs(self.face)
        self.n = abs(base)**2
        self.j = current
        self.rate = rate
        x = slab.x[:, None]
        support = np.maximum(1-x*x-slab.r**2, 0)
        self.g = support**3
        self.speed = current/self.n[-1]
        self.normal_factor = -self.g*(1-1/self.amp) + 1j*self.speed
        self.potential = (self.n**4-1)/4

        interp = CubicSpline(slab.x, density_ref, axis=0)
        shape = self.g/density_ref
        normalization = np.sum(shape*self.n[-1]*slab.dv)
        dg = 6*x*support**2
        dshape = (dg + shape*interp(slab.x, 1))/density_ref
        dnorm = np.sum(dshape*self.n[-1]*slab.dv)
        dspeed = (rate*dshape-self.speed*dnorm)/normalization
        self.displaced_normal = (-dg*(1-1/self.amp)+1j*dspeed)*self.face

        def trace_normal(offset):
            g = np.maximum(1-(x-offset)**2-slab.r**2, 0)**3
            shape_at_offset = g/interp(slab.x-offset)
            speed = rate*shape_at_offset/np.sum(shape_at_offset*self.n[-1]*slab.dv)
            return (-g*(1-1/self.amp)+1j*speed)*self.face

        eps = 1e-4*slab.dx
        numerical = (trace_normal(eps)-trace_normal(-eps))/(2*eps)
        slope_error = np.max(abs(numerical-self.displaced_normal))
        trace_error = np.max(abs(trace_normal(0)-self.normal_factor*self.face))
        if not np.isfinite([slope_error, trace_error]).all() or slope_error > 3e-8 or trace_error > 2e-11:
            raise ArithmeticError(f"Translated intake law mismatch: {slope_error=} {trace_error=}")
        print(f"translation_derivative_error={slope_error:.6g} base_trace_error={trace_error:.6g}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("calibration", type=Path)
    parser.add_argument("isolated", type=Path)
    parser.add_argument("rate", type=float, choices=(0.125, 0.5))
    args = parser.parse_args()
    archive = json.loads(Path("docs/_measurements/lateral_target_rate.json").read_text())
    hashes = {Path(k): v for run in archive["runs"] for group in ("reference_fields", "fields")
              for k, v in run[group].items()}
    with verified_packet(args.calibration, hashes) as packet:
        x, r, dz = packet["x"], packet["r"], float(packet["dz"])
        slab = LateralSlab(round(r[-1]+r[0]), round(1/(2*r[0])), round(1/dz))
        if (int(packet["d"]) != 0 or not np.array_equal(x, slab.x)
                or not np.array_equal(r, slab.r) or dz != slab.dz):
            raise ValueError("Wrong isolated calibration geometry")
        density_ref = abs(packet["field"][-1])**2
    with verified_packet(args.isolated, hashes) as packet:
        if (not np.array_equal(packet["centers"], [0]) or not np.array_equal(packet["rates"], [args.rate])
                or not np.array_equal(packet["x"], x) or not np.array_equal(packet["r"], r)
                or float(packet["dz"]) != dz):
            raise ValueError("Isolated target preparation mismatch")
        layers, currents = packet["layers"], packet["currents"]
    face = layers[-1]
    g = np.maximum(1-x[:, None]**2-r**2, 0)**3
    normal = (-g*(1-1/abs(face))+1j*currents[0]/abs(face)**2)*face
    base, layer_error = extend_archived_field(slab, layers, normal)
    experiment = CalibratedTranslation(slab, base, currents[0], args.rate, density_ref)
    residual, bath, rates = check(slab, base, experiment.normal_factor*base[-1], currents, np.array([args.rate]))
    p, static_residual = experiment.static_tangent()
    print(f"extent={r[-1]+r[0]:g} dx={slab.dx:.8g} dz={dz:.8g} rate={rates[0]:.12g} "
          f"bath={bath:.12g} layer_error={layer_error:.6g} base_residual={residual:.6g} "
          f"static_residual={static_residual:.6g}", flush=True)
    for omega in (0., 0.005, 0.01):
        h, residual, steps = experiment.dynamic_tangents(p, omega)
        carry, stress, direct = experiment.load(h)
        reaction = -(carry+stress)
        print(f"omega={omega:g} carried_per_A={carry:.12g} stress_per_A={stress:.12g} "
              f"compensating_per_A={reaction:.12g} direct_compensating_per_A={-direct:.12g} "
              f"per_rate={reaction/args.rate:.12g} per_rate_squared={reaction/args.rate**2:.12g} "
              f"residual={residual:.6g} steps={steps}", flush=True)


if __name__ == "__main__":
    main()
