"""Compare velocity-odd mixed flux from the mouth law and feeding-side layers."""

import json
from pathlib import Path

import numpy as np

from lateral_source_rate import verified_packet
from mouth_feed_weak_check import dx4


def read_flux(packet):
    f, u = packet["field"][-4:], packet["u"][-4:]
    face, du = f[-1], u[-1]
    x, r, dz = packet["x"], packet["r"], float(packet["dz"])
    dx, radius = 2*r[0], float(packet["radius"])
    weight = 2*np.pi*r*dx**2
    g = np.maximum(1-(x[:, None]/radius)**2-(r[None, :]/radius)**2, 0)**3
    mouth = g > 0
    # Include all cells touched by D_x acting on a supported mouth quantity.
    halo = np.logical_or.reduce([np.roll(mouth, k, axis=0) for k in range(-2, 3)])
    dn = 2*(f.conj()*u).real
    normal = packet["normal"]
    dnormal = du*(normal/face)-g*face*dn[-1]/(2*abs(face)**3)
    density_normal_change = 2*(du.conj()*normal+face.conj()*dnormal).real

    def sums(normal_f, normal_u, normal_dn):
        product = (dx4(du, dx).conj()*normal_f + dx4(face, dx).conj()*normal_u).real
        divergence = dx4(normal_dn, dx)/4
        flux = (product-divergence)*weight
        return flux[mouth].sum(), flux[halo].sum(), np.sum(product*weight*halo)

    def from_layers(q):
        return (11*q[-1]-18*q[-2]+9*q[-3]-2*q[-4])/(6*dz)

    law_cut, law_halo, law_link = sums(normal, dnormal, density_normal_change)
    layer_cut, layer_halo, layer_product = sums(from_layers(f), from_layers(u), from_layers(dn))
    return dict(law_cut=law_cut, law_halo=law_halo, law_link=law_link, layer_cut=layer_cut,
                layer_halo=layer_halo, layer_product_halo=layer_product)


def main():
    archive = json.loads(Path("docs/_measurements/lateral_radius_momentum.json").read_text())
    for row in archive["runs"]:
        path, = map(Path, row["fields"])
        with verified_packet(path, {path: row["fields"][str(path)]}) as packet:
            values = read_flux(packet)
        print(f"grid={row['grid']} radius={row['radius']:g} " +
              " ".join(f"{k}={v:.12g}" for k, v in values.items()), flush=True)
        if (not np.isfinite(list(values.values())).all() or
                abs(values["law_halo"]-values["law_link"]) > 5e-15):
            raise ArithmeticError("Full supported normal-density derivative did not telescope")


if __name__ == "__main__":
    main()
