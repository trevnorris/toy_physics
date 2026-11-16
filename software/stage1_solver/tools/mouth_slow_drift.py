"""Read the velocity-odd momentum flux through an isolated translating mouth."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from lateral_mouth_acceleration import Translation
from lateral_paired_acceleration import target_load
from lateral_source_rate import verified_packet
from mouth_feed_weak_check import dx4


def read_drift(packet, radius):
    x, r = packet["x"], packet["r"]
    dx = 2*r[0]
    slab = SimpleNamespace(dx=dx, dz=float(packet["dz"]), dv=2*np.pi*r*dx**2)
    f, u = packet["field"][-4:], packet["u"][-4:]
    face, du = f[-1], u[-1]
    n, amp = abs(f)**2, abs(face)
    g = np.maximum(1-(x[:, None]/radius)**2-(r[None, :]/radius)**2, 0)**3
    # The stored normal derivative contains the boundary speed in its phase part.
    speed = (packet["normal"]/face).imag
    if (float(packet["radius"]) != radius or x[2] >= -radius or x[-3] <= radius
            or r[-1] <= radius or np.any(speed[g == 0] != 0)):
        raise ValueError("Local packet does not contain the complete specified intake")
    j = speed*n[-1]
    rate = np.sum(slab.dv*j)
    fields = SimpleNamespace(slab=slab, face=face, n=n, f=f, amp=amp, g=g,
                             speed=speed, j=j, rate=rate)
    carry, stress, interior = Translation.load(fields, u)
    dj = speed*2*(face.conj()*du).real
    vx = (face.conj()*dx4(face, dx)).imag/n[-1]
    d_rate = np.sum(slab.dv*dj)
    # The archived tangent holds normal speed fixed. Undo the fixed-rate
    # normalization in Translation.load, even though symmetry nulls this term.
    carry += d_rate*np.sum(slab.dv*j*vx)/rate
    redistributed = np.sum(slab.dv*dj*vx)
    lab = carry+stress

    def odd_reading(step):
        plus = target_load(face+step*du, g, speed, slab)
        minus = target_load(face-step*du, g, speed, slab)
        return (plus-minus)/(2*step)

    return dict(rate=rate, rate_slope=d_rate, rest=target_load(face, g, speed, slab),
                redistribution=redistributed, velocity=carry-redistributed, carry=carry,
                stress=stress, lab=lab, interior=interior, mouth_frame=lab-rate,
                fd_1em4=odd_reading(1e-4), fd_5em5=odd_reading(5e-5))


def main():
    folder = Path("docs/_measurements")
    packets = json.loads((folder/"lateral_radius_momentum.json").read_text())
    held = json.loads((folder/"lateral_mouth_radius.json").read_text())
    for row in packets["runs"]:
        grid, radius = row["grid"], row["radius"]
        path, = map(Path, row["fields"])
        with verified_packet(path, {path: row["fields"][str(path)]}) as packet:
            values = read_drift(packet, radius)
        print(f"grid={grid} radius={radius:g} " +
              " ".join(f"{k}={v:.12g}" for k, v in values.items()), flush=True)
        v = values
        if (not np.isfinite(list(v.values())).all() or v["rate"] <= 0
                or abs(v["rest"]) > 1e-11 or abs(v["rate_slope"]) > 1e-10 or
                not np.allclose([v["fd_1em4"], v["fd_5em5"]], v["lab"],
                                rtol=1e-6, atol=2e-13)):
            raise ArithmeticError("Steady or velocity-odd load check failed")
        if radius == 1:
            comparison, = [h for h in held["runs"] if h["grid"] == grid and h["radius"] == 1]
            for line in comparison["stdout"]:
                if line.startswith(("d=4 ", "d=8 ")):
                    data = dict(word.split("=") for word in line.split())
                    magnitude = float(data["load"])
                    print(f"grid={grid} radius=1 d={data['d']} held={magnitude:.12g} "
                          f"lab_at_V0.01_over_held={0.01*v['lab']/magnitude:.12g}",
                          flush=True)


if __name__ == "__main__":
    main()
