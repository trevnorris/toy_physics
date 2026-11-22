"""Finite-reservoir regression checks for the rigid moving-mouth force."""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np


TOOLS = Path(__file__).resolve().parents[1] / "tools"


def test_external_push_work_matches_open_reservoir_energy():
    output = subprocess.check_output([
        sys.executable, str(TOOLS / "mouth_rigid_push_energy.py"),
        "--mesh", "12", "--omega", "0.01", "--rate", "0.25",
    ], text=True)
    rows = [json.loads(line) for line in output.splitlines()]
    push = next(r for r in rows if r["test"] == "unit_velocity")
    boundary = next(r for r in rows if r["test"] == "boundary_outward_power_per_speed_squared")
    np.testing.assert_allclose(push["actuator_mean_power_per_speed_squared"],
                               boundary["remote_total"], rtol=1e-7, atol=1e-12)


def test_incident_pressure_time_term_at_both_reservoirs(monkeypatch):
    monkeypatch.syspath_prepend(str(TOOLS))
    from mouth_carrier_free_center import momentum_force_rows, motion_source
    from mouth_empty_bend_stiffness import quadrature
    from mouth_flared_receiver import EXIT, FlaredDomain

    dom = FlaredDomain(10, "cubic", 10, 9)
    f = np.ones(len(dom.points), complex)
    static, temporal = momentum_force_rows(dom, f, motion_source(dom, f))
    # A uniform trial phase slope has h=i*r on both cylindrical reservoir
    # sides. Each side contributes -4π R³/3 per unit depth from the temporal
    # pressure. The roof normal has no x component.
    radius = dom.points[:, 0]
    h = np.r_[np.zeros_like(radius), radius]
    rq, _, radial_gradient, _ = quadrature(dom)
    stored = dom.weight @ (radial_gradient @ radius + 2*(dom.N @ radius)/rq) / 3
    port = -4*np.pi*dom.far**3*(1 + dom.roof - EXIT)/3
    np.testing.assert_allclose(static @ h, 0., atol=1e-11)
    np.testing.assert_allclose(temporal @ h, stored + port, rtol=1e-10, atol=1e-9)
