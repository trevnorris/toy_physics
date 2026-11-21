"""Eliminate the unknown contact birth rate between two measured carrier paths."""
import json

import numpy as np

from mouth_contact_age import PRIOR, emit, trace_path


def residence_bounds(first, second):
    # Every interval in this union belongs to one linear segment of each
    # reconstructed inverse path. This avoids smoothing out a concentration of
    # births in a narrow band of order.
    def inverse(path):
        t, _, _, q = path
        knots, times = q[::-1].copy(), t[::-1].copy()
        if max(abs(knots[[0, -1]]-[.25, .75])) > 1e-8:
            raise ValueError("Mismatched endpoints")
        knots[0], knots[-1] = .25, .75
        slopes = -np.diff(times)/np.diff(knots)
        if not np.all(slopes > 0):
            raise ValueError("Order must fall continuously on each path")
        return knots, slopes
    q1, s1 = inverse(first)
    q2, s2 = inverse(second)
    edges = np.union1d(q1, q2)
    mids = edges[:-1]+np.diff(edges)/2
    s1 = s1[np.searchsorted(q1[1:-1], mids)]
    s2 = s2[np.searchsorted(q2[1:-1], mids)]
    times = np.diff(edges)
    ratios = s2/s1
    # A linear renewal law within each segment integrates exactly at midpoint.
    exponent = (times*(1-mids)) @ s2 / ((times*(1-mids)) @ s1)
    # Extreme weights supported on one bin attain the bounds on these paths.
    i, j = np.argmin(ratios), np.argmax(ratios)
    return dict(exponent_min=ratios[i], exponent_max=ratios[j], linear_rate_exponent=exponent,
                min_bin=edges[i:i+2].tolist(), max_bin=edges[j:j+2].tolist(),
                pieces=len(mids), example_F_near=.9,
                deep_interval=[.9**ratios[j], .9**ratios[i]], deep_for_linear_rate=.9**exponent)


def run():
    archive = json.loads(PRIOR.read_text())["runs"]["mouth_mixed_order_ripple02"]
    patches = [json.loads(line) for line in archive["stdout_lines"]]
    cases = {(p["mesh"],p["profile"],p["state"]):p for p in patches if p["test"]=="archive_patch"}
    for mesh in (34,48):
        for profile in ("cubic","delayed"):
            for resolution in (1000,2000,4000):
                paths = [trace_path(cases[mesh,profile,state],resolution) for state in ("near","deep")]
                result = residence_bounds(*paths)
                emit(profile=profile, mesh=mesh, steps_per_radius=resolution, **result)
                assert result["exponent_min"] <= result["linear_rate_exponent"] <= result["exponent_max"]
                assert result["deep_interval"][0] <= result["deep_for_linear_rate"] <= result["deep_interval"][1]


if __name__ == "__main__":
    run()
