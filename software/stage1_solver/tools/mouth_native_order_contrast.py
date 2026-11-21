"""Read the initial response to a compact order offset on archived fixed carriers.

The Q2 weighted residual and mass matrix are those of the passive Allen–Cahn
solve. Measure the response where the imposed offset is initially uniform, on
one carrier streamline per preparation. Time is the original model clock.
"""
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.sparse.linalg import splu

from mouth_free_contact import Carrier
from mouth_free_order_diffuse import NU, Order, emit, read


MEASUREMENTS = Path("research/pde_ledger_v3/_measurements")
LEVELS = (.1, .25, .5, .75)
AMPLITUDE = 1e-5


def locate_probes(domain, psi, chi, r_exit, w_exit):
    """Trace one parcel forward/backward from its archived quarter-order point."""
    def radial_slope(w, radius):
        value, dr, dw = read(domain, psi, radius[0], w)
        ur, uw = (dr/value).imag, (dw/value).imag
        if uw <= 0 or not 0 < radius[0] < domain.radius(w)[0]:
            raise ValueError("Parcel left the inward conduit flow")
        return [ur/uw]

    def reach(target, end):
        def crossing(w, radius):
            return read(domain, chi, radius[0], w)[0]-target
        crossing.terminal = True
        path = solve_ivp(radial_slope, (w_exit, end), [r_exit], rtol=2e-9,
                         atol=2e-11, max_step=.002, events=crossing, dense_output=True)
        if not path.success or len(path.t_events[0]) != 1:
            raise ValueError(f"Streamline failed to reach order {target}")
        return path

    upstream = reach(max(LEVELS), w_exit-.18)
    downstream = reach(min(LEVELS), w_exit+.18)
    result = {}
    for level in LEVELS:
        path = upstream if level > .25 else downstream
        if level == .25:
            w = w_exit
        elif level in (min(LEVELS), max(LEVELS)):
            w = path.t_events[0][0]
        else:
            w = brentq(lambda z: read(domain, chi, path.sol(z)[0], z)[0]-level,
                       min(path.t[0], path.t[-1]), max(path.t[0], path.t[-1]))
        result[level] = float(path.sol(w)[0]), float(w)
    return result


def flat_window(offset, flat, cutoff):
    t = np.clip((abs(offset)-flat)/(cutoff-flat), 0., 1.)
    return 1.-t**3*(10.-15.*t+6.*t*t)


def rate_field(domain, mass_factor, force):
    rate = np.zeros(len(domain.points))
    rate[domain.free] = mass_factor.solve(-force)
    return rate


def measure(patch, domain, psi, model, mass_factor, hashes):
    order_path = Path(patch["order_file"])
    with order_path.open("rb") as stream:
        digest = file_digest(stream, "sha256").hexdigest()
    if digest != patch["order_sha256"] or digest != hashes[str(order_path)]:
        raise ValueError(f"Changed order input {order_path}")
    with np.load(order_path) as saved:
        chi = saved["chi"]
    probes = locate_probes(domain, psi, chi, patch["r"], patch["w"])
    rc, wc = probes[.5]
    bump = (flat_window(domain.points[:, 0]-rc, .13, .235)
            * flat_window(domain.points[:, 1]-wc, .17, .30))
    if np.any(bump[~domain.free]):
        raise ValueError("Pulse intersects a distant imposed order value")
    residual, jac = model.measure(chi)
    free_bump = bump[domain.free]
    jac_bump = jac @ free_bump
    tangent = rate_field(domain, mass_factor, jac_bump)
    plus = model.measure(chi+AMPLITUDE*bump)[0]
    minus = model.measure(chi-AMPLITUDE*bump)[0]
    centered = rate_field(domain, mass_factor, (plus-minus)/(2*AMPLITUDE))
    weighted_norm_rate = -(free_bump @ jac_bump)/(free_bump @ (model.mass @ free_bump))
    rows = []
    for level, (r, w) in probes.items():
        q, qr, qw = read(domain, chi, r, w)
        value, dr, dw = read(domain, psi, r, w)
        ur, uw = (dr/value).imag, (dw/value).imag
        drift = ur*qr+uw*qw
        contrast, hr, hw = read(domain, bump, r, w)
        advection = ur*hr+uw*hw
        # The weighted semidiscrete mass solves return Eulerian derivatives.
        # On the material parcel the sampled contrast also advects.
        fe = (read(domain, tangent, r, w)[0]+advection)/contrast
        fd = (read(domain, centered, r, w)[0]+advection)/contrast
        local = -model.scale*(1.-6.*q*(1.-q))
        rows.append(dict(level=level, r=r, w=w, chi=q, parcel_order_drift=drift,
                         offset=contrast, offset_gradient=[hr, hw],
                         continuum_material_tangent=local, fe_material_tangent=fe,
                         centered_material_tangent=fd, fixed_rate_max_tangent=drift/q))
    emit(test="compact_order_offset", profile=patch["profile"], state=patch["state"],
         mesh=patch["mesh"], intake=patch["intake"], order=str(order_path),
         order_sha256=digest, carrier=patch["carrier_file"], carrier_sha256=patch["carrier_sha256"],
         nu=NU, width=model.delta, amplitude=AMPLITUDE,
         maximum_stationary_weak_residual=np.max(abs(residual)),
         weighted_initial_norm_log_slope=weighted_norm_rate, probes=rows)
    for row in rows:
        if (abs(row["level"]-row["chi"]) > 1e-9 or abs(row["offset"]-1.) > 1e-12
                or np.linalg.norm(row["offset_gradient"]) > 1e-11
                or abs(row["centered_material_tangent"]-row["fe_material_tangent"]) > 1e-5
                or row["parcel_order_drift"] >= 0.):
            raise ArithmeticError("Probe location, flat offset, or residual difference check failed")


def main():
    mixed = json.loads((MEASUREMENTS/"mouth_mixed_order_ripple.json").read_text())
    depth = json.loads((MEASUREMENTS/"mouth_free_order_depth.json").read_text())
    patches = [json.loads(s) for s in mixed["runs"]["mouth_mixed_order_ripple02"]["stdout_lines"]]
    patches = [p for p in patches if p["test"] == "archive_patch"]
    for source in sorted({p["carrier_file"] for p in patches}):
        group = [p for p in patches if p["carrier_file"] == source]
        if any(p["carrier_sha256"] != depth["input_carriers_sha256"][source] for p in group):
            raise ValueError("Carrier digests disagree across archives")
        carrier = Carrier(Path(source).parent)
        key = f"mouth_free_order_{group[0]['profile']}_n{group[0]['mesh']}_d04"
        setup = json.loads(depth["runs"][key]["stdout_lines"][0])
        if setup["source"] != source or setup["nu"] != NU:
            raise ValueError("Order operator and carrier archive disagree")
        model = Order(carrier.domain, carrier.field, setup["delta"])
        mass_factor = splu(model.mass.tocsc())
        for patch in group:
            measure(patch, carrier.domain, carrier.field, model, mass_factor,
                    depth["diffuse_orders_sha256"])


if __name__ == "__main__":
    main()
