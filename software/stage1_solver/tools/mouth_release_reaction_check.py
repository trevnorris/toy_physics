"""Independent parcel-energy and reaction check for steady local release tests."""
import math
from dataclasses import dataclass

from numpy.polynomial.legendre import leggauss


PERIOD = 2 * math.pi
DECISION_RATE = 1 / 128
REGISTRATION_RATE = 1 / 256


def phase_density(age, total_rate):
    """Sum decision density, optionally weighted by surviving swivel, over cycles."""
    return DECISION_RATE * math.exp(-total_rate * age) / -math.expm1(-total_rate * PERIOD)


@dataclass
class Parcel:
    velocity: float
    extension: float = 0.0

    def advance(self, anchor, duration):
        # In the frame of the constant-speed anchor, the oscillator rotates
        # the vector (v - anchor, extension) without changing its length.
        relative = self.velocity - anchor
        cosine = math.cos(duration)
        sine = math.sin(duration)
        return Parcel(anchor + relative * cosine + self.extension * sine,
                      self.extension * cosine - relative * sine)


def check_case(label, policy, lever, incoming, rim, receiver, nodes):
    reaction_rim = reaction_bulk = net_export = lost = net_kinetic = 0.
    count = tail = max_residual = 0.
    for age, weight in nodes:
        total = phase_density(age, DECISION_RATE)
        swivel = phase_density(age, DECISION_RATE + REGISTRATION_RATE)
        delayed = {"registered": total - swivel, "all": total, "none": 0.}[policy]
        captured = Parcel(incoming).advance(rim, age)
        first = captured.velocity - incoming
        for continue_motion, density in ((True, delayed), (False, total - delayed)):
            if continue_motion:
                anchor = (1 - lever) * rim + lever * receiver
                exit_state = captured.advance(anchor, PERIOD - age)
                tail += weight * density * (PERIOD - age)
            else:
                exit_state = captured
            second = exit_state.velocity - captured.velocity
            rim_impulse = first + (1 - lever) * second
            receiver_impulse = lever * second
            heat = exit_state.extension ** 2 / 2
            kinetic = (exit_state.velocity ** 2 - incoming ** 2) / 2
            work = rim * rim_impulse + receiver * receiver_impulse
            max_residual = max(max_residual, abs(work - kinetic - heat))
            factor = weight * density
            count += factor
            reaction_rim += factor * rim_impulse
            reaction_bulk += factor * receiver_impulse
            net_export += factor * (exit_state.velocity - incoming)
            lost += factor * heat
            net_kinetic += factor * kinetic

    print(f"{label} policy={policy} h={lever:g} (U,V,W)=({incoming:g},{rim:g},{receiver:g}) "
          f"rim={reaction_rim:.12g} receiver={reaction_bulk:.12g} net_export={net_export:.12g} "
          f"net_kinetic={net_kinetic:.12g} release_energy={lost:.12g} "
          f"mean_extra_residence={tail:.12g} parcel_energy_residual={max_residual:.3g}")
    assert abs(count - 1.) < 1e-12
    assert abs(reaction_rim + reaction_bulk - net_export) < 1e-12
    assert max_residual < 1e-12
    return reaction_rim


def main():
    positions, weights = leggauss(32)
    nodes = [(PERIOD * (1 + x) / 2, PERIOD * w / 2) for x, w in zip(positions, weights)]
    a = 1 / (1 + DECISION_RATE ** 2)
    survival = DECISION_RATE / (DECISION_RATE + REGISTRATION_RATE)
    unregistered = survival / (1 + (DECISION_RATE + REGISTRATION_RATE) ** 2)
    b = a - unregistered
    print(f"immediate_gain={a:.12g} unregistered_gain={unregistered:.12g} "
          f"registered_gain={b:.12g} registered_fraction={1 - survival:.12g}")
    for policy in ("none", "registered", "all"):
        returned = {"none": 0., "registered": b, "all": a}[policy]
        for h in (0., .5, 1.):
            c = a - returned + h * returned
            g = a - returned + h ** 2 * returned
            held = check_case("held", policy, h, -1., 0., 0., nodes)
            tow = check_case("tow", policy, h, 0., 1., 0., nodes)
            assert abs(held - c) < 1e-12 and abs(tow - g) < 1e-12
            assert g + 1e-12 >= c ** 2 / a
    h = .5
    for u, v, w in ((1., 0., 1.), (0., 1., 1.), (.5, .5, .5), (.3, -.1, .2)):
        rim = check_case("independent_drives", "registered", h, u, v, w, nodes)
        expected = ((unregistered + h**2 * b) * (v - w)
                    + (unregistered + h * b) * (w - u))
        assert abs(rim - expected) < 1e-12


if __name__ == "__main__":
    main()
