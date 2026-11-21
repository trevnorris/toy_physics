"""Steady transverse momentum and work on a draining patch during spring release.

All impulse and energy outputs are per unit mass throughput. Clocks and spring
units match the previous registration test. A reaction assigned to a second
support has a separately recorded conjugate velocity and work.
"""
import math
from typing import NamedTuple

import numpy as np
from numpy.polynomial.legendre import leggauss


TAU = 128.0
REGISTER_TIME = 256.0
OMEGA = 1.0
PERIOD = 2 * math.pi
NODES, WEIGHTS = leggauss(80)
PHASE = math.pi * (1 + NODES)


def decision_measure(age_rate):
    """Fold decision ages onto one oscillator period, retaining the exit hazard."""
    phase_rate = age_rate / OMEGA
    return (math.pi * WEIGHTS * np.exp(-phase_rate * PHASE)
            / (TAU * OMEGA * -math.expm1(-PERIOD * phase_rate)))


ALL = decision_measure(1 / TAU)
UNREGISTERED = decision_measure(1 / TAU + 1 / REGISTER_TIME)
REGISTERED = ALL - UNREGISTERED


def advance(velocity, stretch_speed, anchor_speed, phase):
    """Propagate one spring with continuous stretch. stretch_speed = omega * stretch."""
    relative = velocity - anchor_speed
    cosine, sine = np.cos(phase), np.sin(phase)
    return (anchor_speed + relative * cosine + stretch_speed * sine,
            stretch_speed * cosine - relative * sine)


class Reading(NamedTuple):
    rim: float
    support: float
    export: float
    kinetic_change: float
    release_heat: float
    extra_dwell: float


def readout(unload_weights, share, incoming, rim, support, phase_error=0.0):
    """Unloading springs can remain linked to the rim or share a moving support."""
    pre_velocity, pre_stretch = advance(incoming, 0.0, rim, PHASE)
    anchor = (1 - share) * rim + share * support
    tail_phase = PERIOD - PHASE + phase_error
    post_velocity, post_stretch = advance(pre_velocity, pre_stretch, anchor, tail_phase)
    direct_weights = ALL - unload_weights
    initial_impulse = ALL @ (pre_velocity - incoming)
    tail_impulse = unload_weights @ (post_velocity - pre_velocity)
    rim_impulse = initial_impulse + (1 - share) * tail_impulse
    support_impulse = share * tail_impulse
    export = (direct_weights @ (pre_velocity - incoming)
              + unload_weights @ (post_velocity - incoming))
    kinetic_change = (direct_weights @ (pre_velocity**2 - incoming**2)
                      + unload_weights @ (post_velocity**2 - incoming**2)) / 2
    release_heat = (direct_weights @ pre_stretch**2 + unload_weights @ post_stretch**2) / 2
    return Reading(rim_impulse, support_impulse, export, kinetic_change, release_heat,
                   unload_weights @ tail_phase / OMEGA)


def main():
    strength = .01
    logged = []
    capture = 1 - np.cos(PHASE)
    full, locked = ALL @ capture, REGISTERED @ capture
    print(f"registered_departure_fraction={sum(REGISTERED):.12g} "
          f"initial_impulse_coefficient={full:.12g} registered_part={locked:.12g} "
          f"unregistered_part={full-locked:.12g}")
    controls = ((1., 0., 0.), (0., 1., 0.), (0., 0., 1.), (1., 1., 1.),
                (.375, -.75, .5))
    for population, weights in (("direct", np.zeros_like(ALL)),
                                ("registered", REGISTERED), ("all", ALL)):
        for share in ((0.,) if population == "direct" else (0., .5, 1.)):
            held, tow, supported = [readout(weights, share, *(strength*x for x in drive))
                                    for drive in controls[:3]]
            c, g = -held.rim / strength, tow.rim / strength
            print(f"release={population} share={share:g} held_gain={c:.12g} "
                  f"tow_drag={g:.12g} support_drive_rim={supported.rim/strength:.12g} "
                  f"tow_support={tow.support/strength:.12g} tow_net_export={tow.export/strength:.12g} "
                  f"tow_exit_kinetic={tow.kinetic_change/strength**2:.12g} "
                  f"tow_release_heat={tow.release_heat/strength**2:.12g} "
                  f"extra_dwell={tow.extra_dwell:.12g}")
            logs = []
            for drive in controls:
                for scale in (strength, -strength, strength/2, -strength/2):
                    u, v, w = (scale*x for x in drive)
                    data = readout(weights, share, u, v, w)
                    logs.extend([abs(data.rim + data.support - data.export),
                                 abs(v*data.rim + w*data.support - data.kinetic_change - data.release_heat),
                                 abs(data.rim - (c*(w-u) + g*(v-w))),
                                 max(0., -data.release_heat)])
            amount = weights @ capture
            logs.extend([abs(c - (full - amount + share*amount)),
                         abs(g - (full - amount + share**2*amount)),
                         max(0., c**2 - full*g)])
            logged.extend(logs)
    # A late release probes the precision required by exact phase recovery.
    late = math.pi / 6
    for weights, label in ((REGISTERED, "registered"), (ALL, "all")):
        held = readout(weights, 0., strength, 0., 0., late)
        tow = readout(weights, 0., 0., strength, 0., late)
        print(f"late_phase=pi/6 release={label} rim_only "
              f"held_gain={-held.rim/strength:.12g} tow_drag={tow.rim/strength:.12g}")
        logged.extend([abs(tow.rim - tow.export),
                       abs(strength*tow.rim - tow.kinetic_change - tow.release_heat)])
    # With a stationary receiving support, only incoming sideways kinetic power
    # and actuator work are available here. Inspect the mixed drive explicitly.
    u, v = strength, 2 * strength
    for weights, label in ((REGISTERED, "registered"), (ALL, "all")):
        data = readout(weights, .5, u, v, 0.)
        available = u**2/2 + v*data.rim
        spent = u**2/2 + data.kinetic_change + data.release_heat
        print(f"mixed release={label} share=0.5 U={u:g} V={v:g} W=0 "
              f"available_per_U2={available/u**2:.12g} outgoing_plus_heat_per_U2={spent/u**2:.12g}")
        logged.extend([abs(available-spent), max(0., -available)])
    hypothetical = u**2/2 - u*v
    print(f"mixed proposed_held_gain=1 proposed_tow_drag=0 U={u:g} V={v:g} W=0 "
          f"available_per_U2={hypothetical/u**2:.12g}")
    print(f"largest_balance_or_control_residual={max(logged):.3g}")
    assert abs(np.sum(ALL) - 1) < 1e-13 and np.min(REGISTERED) > 0
    assert max(logged) < 1e-12


if __name__ == "__main__":
    main()
