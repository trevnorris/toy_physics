"""Two-sided held-shape work using matched material coordinates and the new symmetric transverse mode."""
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_live_waist_fields import Carrier
from mouth_live_waist_order import Order
from mouth_sealed_shape_relax import Meridian

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "research/pde_ledger_v3/_measurements"


def main():
    wave_manifest = json.loads((DATA/"mouth_unified_transverse.json").read_text())
    order_manifest = json.loads((DATA/"mouth_live_waist_order.json").read_text())
    previous = json.loads((DATA/"mouth_coupled_waist_balance.json").read_text())
    # The wave manifest separately pins the two hashes and the source snapshots.
    signed = ROOT/"_scratch/s11c/unified_transverse01/signed/stdout"
    guard = json.loads((signed.parent/"outcome.json").read_text())
    digest = hashlib.sha256(signed.read_bytes()).hexdigest()
    if (guard["exitCode"] != 0 or not guard["limitsVerified"]
            or digest != wave_manifest["runs"]["signed"]["transcript_sha256"]):
        raise ValueError("The signed-wave input has not been verified")
    series = []
    for line in signed.read_text().splitlines():
        wave = json.loads(line)
        case = wave["case"]
        prep = json.loads(order_manifest["runs"][case]["literal_stdout"].splitlines()[0])
        if prep["radius"] not in (1, 2) or prep["point"] != "broad" or prep["J"] != 1.:
            raise ValueError("Unexpected held state")
        path = Path("_scratch/s11c/live_waist_dynamics01/fields") / (case+".npz")
        fingerprint = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        if not (fingerprint == order_manifest["saved_field_sha256"][str(path)]
                == wave["frozen_material_order_sha256"]):
            raise ValueError("Order and wave use different material preparations")
        with np.load(ROOT/path) as source:
            chi = source["order"]
        radius, n, m = prep["radius"], prep["mesh"], prep["field_mesh"]
        old, = [r for r in previous["readings"] if
                (r["radius"], r["state"], r["order_mesh"]) == (radius, "broad", n)]
        x0 = np.array(prep["coordinates"])
        for axis, key, index, factor in (("waist", "radial", 2, 1/radius),
                                         ("log_depth", "axial", 0, 1.)):
            for pair in wave["signed_wave_forces"][axis]["symmetric"]:
                h = pair["offset"]
                legs = []
                for sign in (-1., 1.):
                    x = x0.copy()
                    x[index] += factor*sign*h
                    shape = Meridian(radius, x)
                    eo, _ = Order(shape, n).energy_gradient(chi)
                    ec, _, flow = Carrier(shape, n, 2*n, m).response(chi, BIAS)
                    ext, = [r["exterior_order"] for r in old["energies"][key]
                            if r["displacement"] == h*sign]
                    freq = pair["high_frequency" if sign > 0 else "low_frequency"]
                    legs.append(dict(order_inside=eo, order_outside=ext, carrier_two_sided=ec,
                                     stored_wave=freq, current=flow["current"],
                                     terminal=flow["terminal_current"], depth=shape.length))
                lo, hi = legs
                gap = 2*h if key=="radial" else hi["depth"]-lo["depth"]
                forces = {f:(lo[f]-hi[f])/gap for f in
                          ("order_inside", "order_outside", "carrier_two_sided", "stored_wave")}
                curl_wave, = [p["force"] for p in wave["signed_wave_forces"][axis]["curl"] if p["offset"]==h]
                old_total = old["forces"][key]["center_half_step" if h==.001 else "center_full_step"]["total"]
                curl_control = sum(forces.values()) - forces["stored_wave"] + curl_wave
                print(json.dumps(dict(case=case, axis=key, offset=h, inward_currents=[lo["current"],hi["current"]],
                                     legs=legs, forces=forces, held_zero_crystal_prestress_total=sum(forces.values()),
                                     old_stationary_total=old_total, frozen_curl_control=curl_control,
                                     input_wave_sha256=digest, input_order_sha256=fingerprint)), flush=True)
                if abs(curl_control/old_total-1.) > .015 or abs(forces["stored_wave"]-pair["force"])>1e-8:
                    raise ArithmeticError("Archived curl control or physical wavelength denominator mismatch")
                series.append(sum(forces.values()))
    if len(series)!=16:
        raise ValueError("Incomplete signed wave and flow pairs")


if __name__ == "__main__":
    main()
