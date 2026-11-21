"""Intrinsic wall-front timing controls and a finite-spacing bound for bend mass."""
import json
import math

import numpy as np
from scipy.integrate import quad

from mouth_bend_depth_participation import FOLDER, spectral_readout
from mouth_carrier_bend_budget import settings
from mouth_projected_bend_spectrum import fields


def bend_weight(u):
    return math.sin(math.pi*u)**4


def transit_bounds(slowness, n, half_length, reference_speed, radial):
    """Use only the segment transit times and the known total wall inventory."""
    edges = np.linspace(0, 1, n+1)
    times = []
    coarse_moment = spread = mean_square = 0.
    for a, b in zip(edges[:-1], edges[1:]):
        width = b-a
        # slowness is dimensionless, reference_speed / local_speed.
        time = half_length/reference_speed * quad(slowness, a, b)[0]
        times.append(time)
        mean = (reference_speed/half_length) * time/width
        mean_weight = quad(bend_weight, a, b)[0]/width
        mean_square += width * mean**2
        coarse_moment += width * mean**2 * mean_weight
        spread += mean**2 * quad(lambda u: (bend_weight(u)-mean_weight)**2, a, b)[0]
    residual = 1-mean_square
    assert residual >= -2e-13
    residual = max(residual, 0.)
    mixed = 2*math.sqrt(residual*spread)
    limits = max(0., coarse_moment-mixed), min(1., coarse_moment+mixed+residual)
    return times, tuple(radial*x/3 for x in limits)


def compare():
    previous = json.loads((FOLDER/"mouth_bend_full_material.json").read_text())["runs"]
    lines = previous["bend_full_material01"]["stdout_lines"]
    prep = fields(next(s for s in lines if s.startswith("preparation ")))
    radial = float(prep["M_radial"])
    half = float(prep["H"])
    speed = float(fields(next(s for s in lines if s.startswith("homogeneous ")))["c_order"])
    print(f"same_intake={prep['intake']} H={half:.12g} radial_mass={radial:.12g} "
          f"speed_for_uniform_wall={speed:.12g}", flush=True)
    _, _, _, _, archive = settings()
    assert float(prep["tension"]) == archive["parameters"]["sigma"]
    transit = []
    for label, shift in (("center", 0.), ("shifted", .15)):
        # Identical material contrast and total inventory. Shift only its axial
        # placement; both halves carry the same profile. Storage is normal only.
        avg = 1 + 8*math.comb(8, 4)/2**8
        def slowness(u):
            return math.sqrt((1 + 8*math.sin(math.pi*(u-shift))**8)/avg)
        norm = quad(lambda u: slowness(u)**2, 0, 1)[0]
        mass = radial/3 * quad(lambda u: bend_weight(u)*slowness(u)**2, 0, 1)[0]
        assert abs(norm-1)<1e-12
        print(f"profile={label} depth_of_min_c={.5+shift:g} c_min={speed*math.sqrt(avg)/3:.12g} "
              f"c_max={speed*math.sqrt(avg):.12g} bend_mass={mass:.12g}", flush=True)
        for n in (1, 4, 8, 16):
            times, (low, high) = transit_bounds(slowness, n, half, speed, radial)
            print(f"times profile={label} stations_intervals={n} total={sum(times):.12g} "
                  f"M_low={low:.12g} M_high={high:.12g} "
                  f"segment_times={','.join(f'{value:.11g}' for value in times)}", flush=True)
            assert low-1e-10 <= mass <= high+1e-10
            if n == 1:
                transit.append(times[0])
        spectral_readout(mass, label)
    assert abs(transit[0]-transit[1])<1e-10


if __name__ == "__main__":
    compare()
