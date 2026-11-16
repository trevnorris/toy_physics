"""Separate regular brane momentum near the mouth from a possible surface term.

The cosine field p gives the first momentum variation per imposed velocity
amplitude, or stored momentum rate per acceleration at the matched crossings.
Integrate over the entire periodic brane plane so both sides of the rim are
included. The last plane uses its trapezoidal volume weight exactly once.
"""

from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from mouth_acceleration_response import load_preparation
from mouth_carrier_momentum import integrate
from mouth_feed_weak_check import dx4


def main():
    manifest = json.loads(Path("docs/_measurements/mouth_acceleration_response.json").read_text())
    fields = {Path(name): digest for run in manifest["runs"]
              for name, digest in run.get("fields", {}).items()}
    budget = { (int(row["n"]),row["law"],float(row["omega"])): float(row["stored_rate_per_A"])
               for run in manifest["runs"] for line in run["stdout"] if "stored_rate_per_A=" in line
               for row in [dict(word.split("=", 1) for word in line.split())] }
    seen = set()
    for path, digest in sorted(fields.items()):
        if sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Archive checksum error: {path}")
        with np.load(path) as packet:
            case = (int(packet["level"]), str(packet["law"]), float(packet["omega"]))
            p = packet["p"]
        n, law, omega = case
        patch, base, _ = load_preparation(n, law)
        def dx(field):
            return dx4(field.swapaxes(0, 1), patch.dx).swapaxes(0, 1)
        response = (p.conj()*dx(base) + base.conj()*dx(p)).imag
        whole = integrate(response, patch.dx, patch.dz)
        if case in seen or abs(whole - budget[case]) > 1e-14:
            raise ArithmeticError("Repeated case or incompatible carrier inventory")
        seen.add(case)
        planes = patch.dx**3 * response.sum(axis=(1,2,3))
        unsigned = patch.dx**3 * abs(response).sum(axis=(1,2,3))
        end_weight = patch.dz * planes[-1] / 2
        print(f"n={n} law={law} omega={omega:g} whole={whole:.12g} "
              f"end_plane_per_depth={planes[-1]:.12g} end_half_cell={end_weight:.12g}")
        for cells in (1,2,4):
            # A pillbox of positive thickness has half weights at both of
            # its ends. Its sum with the complementary slab restores the
            # original inventory including the shared cut.
            layer = integrate(response[-cells-1:], patch.dx, patch.dz)
            complement = integrate(response[:-cells], patch.dx, patch.dz)
            absolute_layer = integrate(abs(response[-cells-1:]), patch.dx, patch.dz)
            thickness = cells*patch.dz
            if (abs(layer+complement-whole)>1e-13 or abs(layer)>absolute_layer+1e-13
                    or absolute_layer>thickness*max(unsigned[-cells-1:])+1e-13):
                raise ArithmeticError("Layer inventory or absolute-momentum bound failed")
            print(f"  thickness={thickness:.12g} layer={layer:.12g} "
                  f"absolute_layer={absolute_layer:.12g} abs_per_depth={absolute_layer/thickness:.12g}")
    if len(seen) != 8 or seen != set(budget):
        raise ArithmeticError("Incomplete acceleration data")


if __name__ == "__main__":
    main()
