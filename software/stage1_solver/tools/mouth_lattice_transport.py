"""Kinematic controls for a crystalline marker in an ideal spherical feeding flow."""
import hashlib
import json
from math import exp, log, pi
from pathlib import Path

from scipy.integrate import quad, solve_ivp
from scipy.special import exp1


ARCHIVE = (Path(__file__).resolve().parents[3] / "research/pde_ledger_v3"
           / "_measurements/mouth_live_waist_order.json")


def intake(archive, radius):
    run = archive["runs"][f"r{radius}_broad_72"]
    transcript = run["literal_stdout"]
    if (run["worker_exit"] != 0 or not run["limits_verified"] or run["swap"] != 0
            or hashlib.sha256(transcript.encode()).hexdigest() != run["stdout_sha256"]):
        raise ValueError("Archived current missing its validated output")
    rows = [json.loads(line) for line in transcript.splitlines()]
    if rows[0]["radius"] != radius or rows[0]["point"] != "broad":
        raise ValueError("Unexpected current preparation")
    state, = [row for row in rows if row["kind"] == "coupled" and row["action"] == 1]
    if abs(state["current"] - state["terminal_current"]) > 1e-8 or state["current"] <= 0:
        raise ArithmeticError("Archived carrier intake mismatch")
    return state["current"]


def shell_trajectory(start_ratio):
    # Dimensionless time A*t/r_final^3; stretches start at unity.
    def rhs(time, state):
        radius, radial, hoop = state
        return -1 / radius**2, 2 * radial / radius**3, -hoop / radius**3

    duration = (start_ratio**3 - 1) / 3
    result = solve_ivp(rhs, (0, duration), (start_ratio, 1, 1),
                       method="DOP853", rtol=2e-12, atol=2e-14)
    if not result.success:
        raise RuntimeError(result.message)
    return duration, result.y[:, -1]


def main():
    archive = json.loads(ARCHIVE.read_text())
    for mouth_radius in (1, 2):
        current = intake(archive, mouth_radius)
        # Extrapolated 3D spherical stream, unit integrated feed density.
        strength = current / (4 * pi)
        read_radius = 4 * mouth_radius
        for ratio in (1.25, 2., 4.):
            duration, (end, radial, hoop) = shell_trajectory(ratio)
            print(json.dumps(dict(entrance=mouth_radius, D=current,
                                  shell_radius=read_radius, start_radius=ratio*read_radius,
                                  time_in_ideal_flow=duration*read_radius**3/strength,
                                  radial_stretch=radial, hoop_stretch=hoop,
                                  volume_ratio=radial*hoop**2, log_hoop=log(hoop))))
            error = max(abs(end-1), abs(radial/ratio**2-1), abs(hoop*ratio-1),
                        abs(radial*hoop**2-1))
            if error > 1e-8:
                raise ArithmeticError("Lagrangian shell integration failed exact-map check")

    # Candidate relaxation control. x=r^3/(3*A*tau); no physical tau is assigned.
    for x in (.25, 1., 4.):
        elastic_hoop = -exp(x)*exp1(x)/3
        integral, uncertainty = quad(lambda y: -exp(-x*(y**3-1))/y, 1, float("inf"))
        memory_ratio = (1 + 1/x)**(1/3)
        print(json.dumps(dict(x=x, steady_log_hoop=elastic_hoop,
                              independent_integral=integral,
                              one_relaxation_time_upstream_radius_ratio=memory_ratio)))
        if abs(integral-elastic_hoop) > 1e-10 or uncertainty > 1e-8:
            raise ArithmeticError("Upstream memory integral failed closed-form check")


if __name__ == "__main__":
    main()
