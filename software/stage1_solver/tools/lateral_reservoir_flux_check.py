"""Read the off-diagonal carrier momentum flux directly from the archived complex field."""

import json
from hashlib import sha256
from pathlib import Path

import numpy as np


def load_from_field(path):
    with np.load(path, allow_pickle=False) as packet:
        key = "field" if int(packet["d"]) == 0 else "layers"
        psi = packet[key][-3:]
        x, r, dz = packet["x"], packet["r"], float(packet["dz"])
    dx = x[1] - x[0]
    # Differentiate entirely from the brane side of the cut, without the
    # prescribed mouth current, density slope or original force reader.
    dx_psi = np.gradient(psi[-1], dx, axis=0, edge_order=2)
    dw_psi = np.gradient(psi, dz, axis=0, edge_order=2)[-1]
    dw_n = np.gradient(abs(psi)**2, dz, axis=0, edge_order=2)[-1]
    flux = (dx_psi.conj()*dw_psi).real - np.gradient(dw_n, dx, axis=0, edge_order=2)/4
    mask = x[:, None]**2 + r**2 < 1
    return -np.sum((flux*2*np.pi*r*dx**2)[mask])


def main():
    archive = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())
    fields = {Path(key): digest for run in archive["runs"]
              for key, digest in run.get("fields", {}).items()}
    folder = Path("_scratch/s11c/lateral_bath_e48_c4_01/fields")
    values = {}
    for d in (0, 3, 4, 6, 8):
        path = folder/f"d{d}.npz"
        if sha256(path.read_bytes()).hexdigest() != fields[path]:
            raise ValueError(f"Archive mismatch: {path}")
        values[d] = load_from_field(path)
    coefficient = 0.25**2/(4*np.pi)  # Comparison only, after integrating every field
    print(f"alone={values[0]:.12g} ideal_d2F={coefficient:.12g}")
    for d in (3, 4, 6, 8):
        values[d] -= values[0]
        if not np.isfinite(values[d]):
            raise ArithmeticError("Nonfinite load")
        print(f"d={d} field_flux={values[d]:.12g} scaled={d*d*values[d]:.12g} "
              f"ideal={coefficient/d**2:.12g}")
    print(f"F6/F3={values[6]/values[3]:.12g} F8/F4={values[8]/values[4]:.12g}")


if __name__ == "__main__":
    main()
