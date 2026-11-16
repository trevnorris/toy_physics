"""Read the frozen one-throat carrier fields and the saved Stage-1 current and drive diagnostics."""

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np


def frozen_current(path):
    packet = json.loads(path.read_bytes())
    mesh, fields = packet["grid"], packet["fields"]
    cfg = packet["config"]["branch"]
    real, imag, aw = (np.asarray(fields[key]) for key in ("psi_R0", "psi_I0", "A_w0"))
    bc = packet["convergence_evidence"]["boundaries"]["matter_w_lower"]
    if not packet["solver"]["converged"] or mesh["w_min"] != 0 or cfg["r_mouth"] <= 0:
        raise ValueError(f"Unexpected one-throat preparation in {path}")
    if (bc["kind"], bc["alpha"], bc["beta"], bc["gamma"]) != ("dirichlet", 1, 0, 0):
        raise ValueError(f"Cannot infer the face current from these boundary data in {path}")
    if any(v.shape != (mesh["nr"], mesh["nw"]) or not np.isfinite(v).all()
           for v in (real, imag, aw)):
        raise ValueError(f"Invalid saved carrier fields in {path}")

    # Cell-centered operator matches coupled_branch._matter_number_current:
    # centered interior differences, second-order one-sided end differences.
    dw = mesh["dw"]
    def axial_current(re, im, dre, dim, gauge):
        return (cfg["hbar"] * (re * dim - im * dre)
                - cfg["gauge_charge"] * gauge * (re**2 + im**2)) / cfg["particle_mass"]

    current = axial_current(real, imag, np.gradient(real, dw, axis=1, edge_order=2),
                            np.gradient(imag, dw, axis=1, edge_order=2), aw)
    # Homogeneous Dirichlet gives ghost = -interior for each carrier component.
    # The face trace is zero and the gauge term vanishes for any finite face A_w.
    # Positive w points from the lower face into this sleeve.
    face, dface = [], []
    for component in (real, imag):
        interior = component[:, 0]
        ghost = -interior
        face.append((interior + ghost) / 2)
        dface.append((interior - ghost) / dw)
    j_face = axial_current(*face, *dface, 0.0)
    # Integrate only the aperture under the saved 4*pi*r^2 dr face measure.
    edges = np.minimum(mesh["r_faces"], cfg["r_mouth"])
    weights = (4 * np.pi / 3) * np.diff(edges**3)
    intake = float(weights @ j_face)
    print(f"{path.name} radius={cfg['r_mouth']:g} max_abs_jw_cells={np.abs(current).max():.12g} "
          f"max_abs_jw_face={np.abs(j_face).max():.12g} signed_inward_number_rate={intake:.12g}")


def saved_diagnostics(root):
    counts = Counter()
    drives = Counter()
    manifests = sorted(root.rglob("manifest.json"))
    for path in manifests:
        data = json.loads(path.read_bytes())
        result = data.get("results", {})
        if not isinstance(result, dict):
            continue
        surrogate = result.get("surrogate_values", {})
        if not isinstance(surrogate, dict):
            surrogate = {}
        for section in (result, surrogate):
            if "spatial_current_l2" in section:
                value = section["spatial_current_l2"]
                if not np.isfinite(value) or value < 0:
                    raise ValueError(f"Invalid current norm in {path}")
                counts[data["benchmark"], value > 0] += 1
        if "omega" in result and "wall_lower_endpoint_eta_abs" in surrogate:
            drives[data["benchmark"], data["config"]["p2_tangent"]["spherical_l"]] += 1
    print(f"Stage-1 run manifests scanned={len(manifests)} current-norm entries={sum(counts.values())}")
    for (benchmark, positive), count in sorted(counts.items()):
        print(f"  {benchmark} current_norm_{'positive' if positive else 'zero'}={count}")
    for (benchmark, ell), count in sorted(drives.items()):
        print(f"  {benchmark} frequency-domain wall responses with ell={ell}: {count}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("frozen_fields", type=Path, nargs="+")
    parser.add_argument("--runs", type=Path, required=True)
    args = parser.parse_args()
    for path in args.frozen_fields:
        frozen_current(path)
    saved_diagnostics(args.runs)
