"""Compare the matched tilt energies and resolve the former barrier in a close-up."""
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    directory = Path("docs/_measurements")
    paths = json.loads((directory / "relaxed_neck_paths.json").read_text())
    reference = json.loads((directory / "flare_elastic_barrier.json").read_text())
    pattern = r"^beta=(\S+) delta_area=(\S+) delta_bend=(\S+) delta_frequency=(\S+)$"
    records = reference["runs"]["sweep"]["stdout_claim_lines"]
    rows = [list(map(float, match.groups())) for line in records if (match := re.match(pattern, line))]
    beta, area, bend, frequency = np.asarray(rows).T
    action = paths["free"]["action"]
    fixed = area + 5*bend + action*frequency
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True,
                             gridspec_kw={"width_ratios":[1.1,1]})
    for ax in axes:
        ax.axhline(0, color="#92979b", linewidth=0.8)
        ax.plot(beta, fixed, color="#b54c31", label="Fixed shape", linewidth=2.1)
        for name, color, label in [("clamped", "#227b80", "Shape free, angle fixed"),
                                   ("free", "#40478e", "Shape and angle free")]:
            states = paths[name]["states"]
            xs = [s["beta"] for s in states]
            ys = [s["energy"]-states[0]["energy"] for s in states]
            ax.plot(xs, ys, color=color, label=label, linewidth=2.1)
        ax.spines[["top","right"]].set_visible(False)
        ax.set_xlabel("Tilt parameter β")
    axes[0].set(xlim=(0, 3), ylim=(-52, 4), title="Path to loss of the count",
                ylabel="Energy change  ΔE / (σR³)")
    axes[0].axvspan(2, 3, color="#e2ebeb", alpha=0.6, zorder=-1)
    axes[0].legend(loc="lower left", frameon=False, fontsize=9)
    axes[1].set(xlim=(0, 1.4), ylim=(-3.6, 0.95), title="Near the former barrier")
    axes[1].annotate("+0.419", (1.101,0.419), xytext=(0,7), textcoords="offset points",
                     ha="center", color="#82331e")
    output = Path("software/em_charge_attribute/reports/relaxed_neck.png")
    fig.savefig(output, dpi=160)
    plt.close(fig)
    for name, path in paths.items():
        states = path["states"]
        initial = np.array(states[0]["parts"])
        for s in states:
            if s["beta"] in (1.1, 3.):
                wall, flex, wave = (np.array(s["parts"])-initial)*[1.,path["stiffness"],path["action"]]
                print(f"{name} beta={s['beta']:g} tension={wall:+.8f} "
                      f"bending={flex:+.8f} wave={wave:+.8f} net={wall+flex+wave:+.8f}")


if __name__ == "__main__":
    main()
