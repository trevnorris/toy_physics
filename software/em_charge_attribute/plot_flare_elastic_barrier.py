"""Resolve the intermediate scale-balanced barrier and plot its energy profile."""
import argparse
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

from flare_elastic_barrier import wall, archived_frequencies
from wave_at_level_patch import Collar


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sweep_stdout", type=Path)
    parser.add_argument("figure", type=Path)
    parser.add_argument("--peak-stdout", type=Path, help="render from the completed fine peak")
    args = parser.parse_args()
    pattern = r"^beta=([\d.]+) delta_area=(\S+) delta_bend=(\S+) delta_frequency=(\S+)$"
    rows = np.array([[float(v) for v in match] for match in re.findall(
        pattern, args.sweep_stdout.read_text(), re.MULTILINE)])
    rows = rows[rows[:, 0] <= 3]
    tilt, areas, bends, frequency_changes = rows.T
    if tilt[0] != 0 or tilt[-1] != 3 or not np.all(np.diff(tilt) > 0):
        raise ValueError("Expected the ordered sweep from zero to complete count loss")
    stiffness = 5.0
    baseline = wall(0, 64).sum(axis=0)

    def balanced_action(omega):
        return (3*baseline[0]+stiffness*baseline[1])/omega

    action = balanced_action(archived_frequencies(24)[0.])
    profile = areas + stiffness*bends + action*frequency_changes
    peak = profile.argmax()
    print(f"coarse peak tilt={tilt[peak]:.6f} energy={profile[peak]:+.9f}", flush=True)
    if not 0 < peak < len(tilt)-1:
        raise ArithmeticError("No interior peak to refine")
    if args.peak_stdout:
        fields = dict(item.split("=") for item in args.peak_stdout.read_text().splitlines()[-1].split())
        if float(fields["k"]) != stiffness:
            raise ValueError("The saved peak uses a different stiffness")
        at_peak, height, at_loss = (float(fields[key]) for key in
                                    ("fine_peak_tilt", "fine_peak_energy", "energy_at_loss"))
    else:
        at_peak, height, at_loss = refine(tilt, peak, baseline, stiffness)
    draw(tilt, profile, at_peak, height, at_loss, args.figure)


def refine(tilt, peak, baseline, stiffness):
    solver = Collar(36, 14)
    frequencies = archived_frequencies(36)
    action = (3*baseline[0]+stiffness*baseline[1])/frequencies[0.]

    def energy(beta):
        if beta not in frequencies:
            frequencies[beta] = solver.solve(beta)[0]
        da, db = wall(beta, 64).sum(axis=0)-baseline
        return da+stiffness*db+action*(frequencies[beta]-frequencies[0.])

    peak_fit = minimize_scalar(lambda x: -energy(x), bounds=(tilt[peak-1], tilt[peak+1]),
                               method="bounded", options={"xatol": 2e-5})
    if not peak_fit.success:
        raise ArithmeticError(peak_fit.message)
    at_peak, height = peak_fit.x, -peak_fit.fun
    at_start, at_loss = energy(2), energy(3)
    print(f"k={stiffness:g} j={action:.9f} fine_peak_tilt={at_peak:.7f} "
          f"fine_peak_energy={height:+.9f} energy_at_onset={at_start:+.9f} "
          f"energy_at_loss={at_loss:+.9f}", flush=True)
    if not height > max(energy(at_peak-0.02),energy(at_peak+0.02),0) or at_loss >= 0:
        raise ArithmeticError("Fine energies do not resolve the barrier example")
    return at_peak, height, at_loss


def draw(tilt, profile, at_peak, height, at_loss, figure):
    fig, ax = plt.subplots(figsize=(7.1, 3.9), constrained_layout=True)
    ax.axvspan(2, 3, color="#e5f0f1")
    ax.axhline(0, lw=0.9, color="#677381")
    ax.plot(tilt, profile, color="#155765", lw=2.5)
    ax.scatter([at_peak, 3], [height, at_loss], color="#155765", zorder=4)
    ax.annotate(f"Barrier  +{height:.2f}", (at_peak,height), xytext=(0,11),
                textcoords="offset points", ha="center")
    ax.annotate(f"Count reaches zero\n{at_loss:.2f}", (3,at_loss), xytext=(-52,0),
                textcoords="offset points", ha="right", va="center",
                arrowprops={"arrowstyle":"-", "color":"#677381", "lw":0.8})
    ax.text(2.04, 0.91, "Count starts to fall", transform=ax.get_xaxis_transform(),
            fontsize=9, color="#335563")
    ax.set(xlabel="Tilt parameter β", ylabel="Energy change  ΔE / (σR³)",
           xlim=(-0.04,3.04), ylim=(at_loss-0.8,height+1.2))
    ax.spines[["right","top"]].set_visible(False)
    fig.savefig(figure, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
