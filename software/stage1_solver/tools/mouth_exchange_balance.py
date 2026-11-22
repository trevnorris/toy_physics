"""Combine the measured radius-one coefficients without choosing a wave occupation."""
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RECORDS = ROOT / "research/pde_ledger_v3/_measurements"


def main():
    inputs = {}
    hashes = {}
    for name in ("mouth_conversion_shape", "mouth_exchange_shear", "mouth_exchange_carrier"):
        path = RECORDS / (name + ".json")
        data = path.read_bytes()
        inputs[name] = json.loads(data)
        hashes[path.name] = hashlib.sha256(data).hexdigest()

    archive = inputs["mouth_conversion_shape"]["primary"]["exchange"]
    literal = archive["literal_stdout"]
    if hashlib.sha256(literal.encode()).hexdigest() != archive["stdout_sha256"]:
        raise ValueError("Reference measurement digest mismatch")
    shape, = (r for r in json.loads(literal) if
              (r["waist"], r["exit_to_waist"], r["cosines"]) == (1, 1.5, 10))
    wave = inputs["mouth_exchange_shear"]
    primary = wave["summary"]["n24_p3_f0"]
    wave_literal = wave["runs"]["n24_p3_f0"]["literal_stdout"]
    if hashlib.sha256(wave_literal.encode()).hexdigest() != wave["runs"]["n24_p3_f0"]["stdout_sha256"]:
        raise ValueError("Wave measurement digest mismatch")
    center_mode, = (r for r in map(json.loads, wave_literal.splitlines())
                    if r["test"] == "mode" and r["step"] == 0)

    carrier = inputs["mouth_exchange_carrier"]
    for record in (wave["reduction"], carrier["reduction"]):
        if record["guard_exit"] or record["worker_exit"] or record["swap_bytes"]:
            raise ValueError("Invalid sector reduction")
    full = carrier["finite"]["mesh64_drive1"]
    pair, = (r for r in full["signed"] if r["half_width"] == .005)
    center, = (r for r in full["states"] if r["step"] == 0)
    if abs(pair["virtual_stiffness"] - pair["energy_stiffness"]) > .00001:
        raise ValueError("Energy and stress stiffness differ")
    threshold = -shape["area_stiffness"] / primary["second_per_action"]
    if abs(threshold - primary["wave_only_required_action_over_tension"]) > 1e-9:
        raise ValueError("Wave and geometry benchmarks disagree")

    joint = {
        "inputs_sha256": hashes,
        "coordinate": "nominal neck displacement divided by reference radius; fixed volume and axial length",
        "units": "each column multiplies sigma R^3, J c_T/R or the carrier energy unit, respectively",
        "first_energy_derivative": {"surface": 0, "wave": primary["first_per_action"],
                                    "carrier": pair["virtual_slope"]},
        "second_energy_derivative": {"surface": shape["area_stiffness"],
                                     "wave": primary["second_per_action"],
                                     "carrier": pair["energy_stiffness"]},
        "fixed_end_first_derivative": {"wave": wave["summary"]["interior_n24"]["frequency_first"],
                                       "carrier": carrier["finite"]["held_ends"]["slope"]},
        "wave_only_action_threshold_Jc_over_sigmaR4": threshold,
        "threshold_stored_wave_energy_over_lateral_surface_energy":
            center_mode["frequency"] * threshold / shape["area"],
        "flow_fraction_through_side": center["current_by_cut"]["side_layer"] / center["current"],
        "energy_loss_fraction_at_widest_neck_over_400pi":
            max(primary["retention"], key=lambda row: row["step"])["loss_fraction_at_400pi"],
    }
    print(json.dumps(joint, indent=2))


if __name__ == "__main__":
    main()
