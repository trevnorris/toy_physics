"""Check the archived collar bounds and an independent age-resolved patch response."""
import cmath
import json
import math
import re
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
ARCHIVE = ROOT / "research/pde_ledger_v3/_measurements/mouth_collar_charge_transmission.json"


def rounded_field(line, name):
    return Fraction(re.search(rf"{name}=([0-9.]+)", line).group(1))


def angular_checks():
    # A resident with independent residence and alignment clocks registers if the
    # latter rings first. The expected remaining residence is still tau.
    alpha = Fraction(1, 256)
    dwell = (1, 128)
    delayed = tuple(2 * tau * alpha / (alpha + Fraction(1, tau)) for tau in dwell)
    immediate = tuple(Fraction(2 * tau) for tau in dwell)
    runs = json.loads(ARCHIVE.read_text())["runs"].values()
    rows = [row for run in runs for row in run["stdout_lines"]
            if "wrong_cut_lower_per_G=" in row]
    weak_margins, strong_margins, cases = [], [], set()
    for line in rows:
        case = line.split(" degrees=", 1)[0]
        cases.add(case)
        floors = re.search(r"wrong_cut_lower_per_G=\[([^]]+)\]", line).group(1)
        floors = [Fraction(entry.strip()) for entry in floors.split(",")]
        carried = (7 * rounded_field(line, "rim_field_per_cut_upper")
                   + rounded_field(line, "radial_vacuum_upper")
                   + delayed[0] * rounded_field(line, "carried_attachment_per_feed_upper"))
        weak_cost = min(rounded_field(line, "charged_prism_lower"),
                        rounded_field(line, "tie_lower"), delayed[1] * floors[-1])
        screened = rounded_field(line, "preferred_trial_upper")
        strong_cost = min(immediate[1] * floors[-1], immediate[0] * floors[0])
        weak_margins.append(weak_cost - carried)
        strong_margins.append(strong_cost - screened)
    print(f"delayed_feed_rim={delayed[0]},{delayed[1]} "
          f"immediate_feed_rim={immediate[0]},{immediate[1]}")
    print(f"archived_cases={len(cases)} delayed_margin_K_ge={min(weak_margins)} "
          f"immediate_margin_K_ge={min(strong_margins)}")
    assert len(cases) == len(rows) == 8 and min(weak_margins + strong_margins) > 0


def spring_geometry():
    k, lever, visibility = 1.9, 1.4, .23
    x, angle = .3 - .4j, 1.65
    # A registered population with mean unit arm direction of modulus visibility.
    def energy(angle):
        total = 0.0
        for phase, weight in ((0, (1 + visibility) / 2), (math.pi, (1 - visibility) / 2)):
            arm_slip = lever * (cmath.exp(1j * angle) - cmath.exp(1j * phase))
            total += weight * k / 4 * (abs(x + arm_slip)**2 + abs(x - arm_slip)**2)
        return total

    reference = k * lever**2 * (1 - visibility)
    measured = energy(angle) - reference
    predicted = k * abs(x)**2 / 2 + k * lever**2 * visibility * (1 - math.cos(angle))
    print(f"arm_energy_minus_constant={measured:.10g} separated_energy={predicted:.10g} "
          f"error={abs(measured - predicted):.3g}")
    assert abs(measured - predicted) < 1e-13


def simpson(function, right=48., intervals=4096):
    step = right / intervals
    odd = sum(function(i * step) for i in range(1, intervals, 2))
    even = sum(function(i * step) for i in range(2, intervals, 2))
    return step / 3 * (function(0.) + function(right) + 4 * odd + 2 * even)


def cohort(age, s, omega, center, incoming):
    # Zero extension and velocity incoming(t-age) at birth. Each exponential
    # below integrates a harmonic drive in the single-parcel oscillator.
    def memory(z):
        return (1 - cmath.exp(-z * age)) / z

    plus, minus = memory(s - 1j * omega), memory(s + 1j * omega)
    old_velocity = incoming * cmath.exp(-s * age)
    extension = center * (plus + minus) / 2 - old_velocity * math.sin(omega * age) / omega
    velocity = center * omega * (plus - minus) / (2j) + old_velocity * math.cos(omega * age)
    return extension, velocity


def translation_checks():
    # Rescale time by this patch's dwell and mass by a parcel mass. Unit
    # arrival rate gives one mean resident. The exit and resident age laws
    # coincide only because the departure hazard is constant.
    for omega in (.5, 2.):
        k = omega**2
        for s, center, incoming in ((0, 1, .25), (.37j, 1, 0), (.37j, 0, 1),
                                    (.37j, 1, 1)):
            def age_mean(index):
                return simpson(lambda a: math.exp(-a) * cohort(a, s, omega, center, incoming)[index])

            force = k * age_mean(0)
            momentum = age_mean(1)
            impedance = k * ((s + 1) * center - incoming) / ((s + 1)**2 + k)
            balance = s * momentum + momentum - incoming
            print(f"patch w_tau={omega:g} s_tau={s} center={center:g} incoming={incoming:g} "
                  f"cohort_force={force:.9g} moment_force={impedance:.9g} "
                  f"momentum_residual={abs(force - balance):.3g}")
            assert max(abs(force - impedance), abs(force - balance)) < 2e-8
            if s == 0:
                def exit_energy(age):
                    e, v = cohort(age, s, omega, center, incoming)
                    return math.exp(-age) * (k * abs(e)**2 + abs(v)**2) / 2

                net_export = simpson(exit_energy) - abs(incoming)**2 / 2
                print(f"dc_center_work={force.real * center:.9g} net_exit_energy={net_export:.9g}")
                assert abs(net_export - force.real * center) < 2e-8
        # At fixed incoming velocity a slow center acceleration adds an
        # O(s) force; the signed coefficient is different from the DC drag.
        h = 1e-5
        slow = simpson(lambda a: math.exp(-a) *
                       cohort(a, 1j * h, omega, 1, 0)[0]).imag * k / h
        derivative = k * (k - 1) / (k + 1)**2
        print(f"w_tau={omega:g} center_slope_per_resident_mass={slow:.9g} "
              f"moment_slope={derivative:.9g}")
        assert abs(slow - derivative) < 2e-8


if __name__ == "__main__":
    angular_checks()
    spring_geometry()
    translation_checks()
