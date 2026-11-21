"""Full-depth momentum of a trial inertial order front at the archived intake."""
import json
import math

from scipy.integrate import quad

from mouth_bend_depth_participation import FOLDER, spectral_readout
from mouth_carrier_bend_budget import settings
from mouth_projected_bend_spectrum import fields


def full_bend():
    point, _, _, flow, model = settings()
    previous = json.loads((FOLDER / "mouth_bend_coefficients.json").read_text())["runs"]
    reference = fields(previous["bend_projection_wave015"]["stdout_lines"][0])
    radial_mass = float(reference["breathing"])
    radius = model["phase_A"]["q_star"][0]
    half_length = float(point["L"]) / 2
    tension = model["parameters"]["sigma"]
    # Three-area of both straight sidewalls in the four-dimensional preparation.
    wall_area = 2 * half_length * 4 * math.pi * radius**2
    mass_per_area = radial_mass / wall_area
    intrinsic_speed = math.sqrt(tension / mass_per_area)
    # Divide out the peak speed. Integrating the normal's x component twice
    # (in the imposed velocity and in its projection onto x) gives 1/3.
    def shape(u):
        return math.sin(math.pi * u)**2

    x_momentum = radial_mass / 3 * quad(shape, 0, 1)[0]
    inertia = radial_mass / 3 * quad(lambda u: shape(u)**2, 0, 1)[0]
    print(f"preparation intake={float(flow['reference_inward_current']):.12g} "
          f"radius={radius:.12g} H={half_length:.12g} M_radial={radial_mass:.12g} "
          f"tension={tension:.12g}")
    print(f"homogeneous mass_per_area={mass_per_area:.12g} c_order={intrinsic_speed:.12g} "
          f"both_Px_over_bdot={x_momentum:.12g} both_weighted_P_over_bdot={inertia:.12g}")
    assert math.isclose(inertia, radial_mass / 8, rel_tol=1e-12)
    assert math.isclose(x_momentum, radial_mass / 6, rel_tol=1e-12)

    for grid in (48, 64):
        fold = fields(previous[f"bend_projection_fold{grid}"]["stdout_lines"][0])
        critical_mass = float(fold["wall_mass_at_fold"])
        target = 3 * critical_mass / radial_mass
        # For a fixed allowed contrast, the best possible normalized profile
        # takes its larger value exactly where the local weight exceeds target.
        edge = math.asin(target**0.25) / math.pi
        center_width = 1 - 2 * edge
        def weight(u):
            return shape(u)**2
        center_integral = quad(weight, edge, 1 - edge)[0]
        excess = center_integral - target * center_width
        mean = 3 / 8
        contrast = 1 + (target - mean) / excess
        normalization = 1 + (contrast - 1) * center_width
        direct = (mean + (contrast - 1) * center_integral) / normalization
        print(f"fold cells={grid} critical_mass={critical_mass:.12g} "
              f"target_mean_f2={target:.12g} center_lo={edge:.12g} "
              f"center_hi={1-edge:.12g} min_mass_contrast={contrast:.12g} "
              f"min_speed_ratio={math.sqrt(contrast):.12g} attained_mass={radial_mass*direct/3:.12g}")
        assert math.isclose(direct, target, rel_tol=1e-12)
        if grid == 64:
            for omega in (.01, .02):
                speed = 1e-5 * omega
                print(f"drive omega={omega:g} bdot_peak={speed:.12g} "
                      f"both_homogeneous_Px={speed*x_momentum:.12g} "
                      f"both_homogeneous_weighted_P={speed*inertia:.12g} "
                      f"both_weighted_P_at_fold={speed*critical_mass:.12g}")
    spectral_readout(inertia, "homogeneous_front")


if __name__ == "__main__":
    full_bend()
