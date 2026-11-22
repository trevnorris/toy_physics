"""Plot two actual meridians from each guarded radius-one and radius-two descent."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mouth_sealed_shape_relax import Meridian


ROOT = Path(__file__).resolve().parents[3]
RECORD = ROOT / "research/pde_ledger_v3/_measurements/mouth_sealed_shape_relax.json"
FIGURE = ROOT / "software/stage1_solver/reports/figures/mouth_sealed_shape_relax.png"


def main():
    training = json.loads(RECORD.read_text())["training"]
    figure, panels = plt.subplots(1, 2, figsize=(9, 4.25))
    for panel, radius in zip(panels, (1, 2)):
        states = training[str(radius)]["states"]
        admitted = [s for s in states if s["maximum_abs_principal_delta"] <= .1]
        for state, label, color, style in (
            (states[0], "Starting profile", "#576773", "--"),
            (admitted[-1], "After relaxation within curvature range", "#126975", "-"),
        ):
            profile = Meridian(radius, state["coordinates"])
            t = np.linspace(0, 1, 501)
            width, _, _ = profile.profile(t)
            if abs(width[0]-radius) > 1e-9 or abs(min(width)-state["minimum_radius"]) > .002:
                raise ValueError("Plotted meridian differs from the recorded state")
            w = t*profile.length
            panel.plot(w, width, color=color, ls=style, lw=2, label=label)
            panel.plot(w, -width, color=color, ls=style, lw=2)
            panel.plot([w[-1], w[-1]], [-width[-1], width[-1]], color=color, ls=style, lw=1)
        panel.axhline(0, color="#aab4b6", lw=.6)
        panel.set(xlabel="Depth into the bulk (reference radii)", title=f"Entrance radius {radius}")
        panel.set_aspect("equal", adjustable="box")
        panel.spines[["right", "top"]].set_visible(False)
    panels[0].set_ylabel("Signed radial position (reference radii)")
    handles, labels = panels[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    figure.tight_layout(rect=(0, .10, 1, 1))
    figure.savefig(FIGURE, dpi=160)
    print(json.dumps({"measurement_sha256": hashlib.sha256(RECORD.read_bytes()).hexdigest(),
                      "figure_sha256": hashlib.sha256(FIGURE.read_bytes()).hexdigest()}))


if __name__ == "__main__":
    main()
