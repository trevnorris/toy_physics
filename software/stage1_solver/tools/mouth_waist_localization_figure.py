"""Plot the relaxed order and measured axial wave-energy distribution in a held throat."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mouth_live_waist_order import state
from mouth_live_waist_wave import OccupiedWave
from mouth_sealed_shape_relax import Meridian


ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "_scratch/s11c/live_waist_dynamics01"
OUTPUT = ROOT / "software/stage1_solver/reports/figures/mouth_live_waist_localization.png"


def main():
    case = BASE / "r2_tight_72"
    outcome = json.loads((case/"outcome.json").read_text())
    if not outcome["limitsVerified"] or outcome["exitCode"] or outcome["childOutcome"]["exitCode"]:
        raise ValueError("Unverified mode source")
    raw = [json.loads(line) for line in (case/"stdout").read_text().splitlines()]
    prep = raw[0]
    center, = [r for r in raw if r["kind"]=="coupled" and r["action"]==1]
    path = BASE / "fields/r2_tight_72.npz"
    with np.load(path) as data:
        chi, vector = data["order"], data["vector"]
        if abs(float(data["frequency"])-center["wave_frequency"])>1e-9:
            raise ValueError("Wrong eigenvector or state")
    shape = Meridian(2, state(2, "tight")["coordinates"])
    nz, nr = np.array(chi.shape)-1
    mesh = prep["field_mesh"]
    wave = OccupiedWave(shape, nr, nz, mesh, 4.2)
    weight = wave.project @ chi.T.ravel()
    v = np.r_[vector, 0j][wave.ids]
    mass = weight * np.einsum("ei,eij,ej->e", v.conj(), wave.mc, v).real
    column = mass.reshape(mesh, -1).sum(axis=0)
    if not np.all(column >= 0):
        raise ValueError("Wave kinetic energy left its positive domain")
    fraction = np.clip(np.cumsum(column)/column.sum(), 0, 1)
    boundaries = np.linspace(0, shape.length, len(column)+1)
    w = np.linspace(0, shape.length, nz+1)
    r = shape.radius(w)[0]
    where = w[np.argmin(r)]
    signed_x = np.linspace(-1, 1, 2*nr+1)
    signed_chi = np.c_[chi[:, :0:-1], chi]
    fig, (upper, lower) = plt.subplots(2, 1, figsize=(8.1, 5.5), sharex=True,
                                      gridspec_kw={"height_ratios":[2.0,1.]})
    axw = np.broadcast_to(w[:, None], signed_chi.shape)
    colors = upper.pcolormesh(axw.T, (r[:, None]*signed_x).T, signed_chi.T,
                              shading="gouraud", cmap="viridis", vmin=.5, vmax=1.)
    upper.plot(w, r, color="#253e55", lw=1)
    upper.plot(w, -r, color="#253e55", lw=1)
    upper.set(ylabel="Radial position", title="Order inside the held radius-2 throat")
    fig.colorbar(colors, ax=upper, label="Ordered fraction", pad=.015)
    lower.plot(boundaries, np.r_[0., fraction]*100, color="#04724d", lw=2.3)
    lower.set(xlim=(0,shape.length), ylim=(-2,103), ylabel="Cumulative cavity\nkinetic energy (%)",
              xlabel="Depth from the entrance toward the bulk (reference radii)")
    lower.grid(axis="y", alpha=.23)
    for ax in (upper,lower):
        ax.axvline(where, ls="--", lw=1.4, color="#b34a1d")
        ax.spines[["top","right"]].set_visible(False)
    lower.annotate("Narrowest section", xy=(where,55), xytext=(5,0), textcoords="offset points",
                   color="#944018", fontsize=9)
    # Reserve the same physical width below the color bar for exact x registration.
    fig.subplots_adjust(left=.14, bottom=.12, right=.89, hspace=.12, top=.93)
    pos = upper.get_position()
    lo = lower.get_position()
    lower.set_position([pos.x0, lo.y0, pos.width, lo.height])
    fig.savefig(OUTPUT, dpi=155)
    print(json.dumps(dict(case=case.name, radius=2, axis_order=center["axis_order"],
                          stored_frequency=center["wave_frequency"],
                          axial_edges=boundaries.tolist(), cumulative_fraction=fraction.tolist(),
                          snapshot_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                          figure_sha256=hashlib.sha256(OUTPUT.read_bytes()).hexdigest())))


if __name__ == "__main__":
    main()
