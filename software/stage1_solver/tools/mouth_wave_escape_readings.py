"""Collect end-point fractions and plot the saved release inside one flared throat."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


def emit(**data):
    print(json.dumps(data, allow_nan=False), flush=True)


def collect(folder):
    if json.loads((folder / "outcome.json").read_text())["exitCode"] != 0:
        raise ValueError(f"Unsuccessful run at {folder}")
    histories = {}
    for line in (folder / "stdout").read_text().splitlines():
        trial = json.loads(line)
        if trial.get("test") != "torsional_pulse":
            continue
        filename = Path(trial["history"])
        if hashlib.sha256(filename.read_bytes()).hexdigest() != trial["history_sha256"]:
            raise ValueError(f"Changed history {filename}")
        with filename.open() as stream:
            rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
        last = rows[-1]
        if last != trial["readings"][-1] or not 11.99 < last["t"] <= 12:
            raise ValueError(f"Final sample mismatch at {filename}")
        returned = last["feed"] + last["feed_absorbed"]
        layers = last["upper_layer"] + last["lower_layer"]
        total = returned + last["tail"] + layers + last["ordered"] + last["receiver"]
        if abs(1 - total) > 1e-8 or trial["maximum_abs_relative_balance_residual"] > 1e-8:
            raise ValueError(f"Unbalanced history {filename}")
        key = trial["state"], trial["pulse"]
        if key in histories:
            raise ValueError(f"Duplicate history {key} in {folder}")
        histories[key] = trial, rows
        emit(folder=str(folder), state=key[0], pulse=key[1], t=last["t"],
             feed_including_absorbed=returned, absorbed_at_remote_feed=last["feed_absorbed"],
             throat_strong=last["ordered"], throat_transition=layers, throat_weak=last["tail"],
             receiving_chamber=last["receiver"], weak_at_release=rows[0]["tail"],
             initial_energy=trial["initial_physical_energy"], balance_error_max=trial["maximum_abs_relative_balance_residual"])
    if len(histories) != 4:
        raise ValueError(f"Incomplete set of releases at {folder}")
    for pulse in ("entrance_velocity", "cavity_displacement"):
        energies = [histories[state, pulse][0]["initial_physical_energy"] for state in ("near", "deep")]
        emit(folder=str(folder), pulse=pulse, deep_over_near_initial_energy=energies[1] / energies[0])
    return histories


def plot(histories, destination):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), sharey=True, layout="constrained")
    for ax, state, title in zip(axes, ("near", "deep"), ("Order ends near the waist", "Order extends into the flare"), strict=True):
        _, rows = histories[state, "cavity_displacement"]
        times = [r["t"] for r in rows]
        series = [
            ("Returned to feeding brane", "#19687a", "-", [r["feed"] + r["feed_absorbed"] for r in rows]),
            ("Throat, shear order below 10%", "#a4551f", "--", [r["tail"] for r in rows]),
            ("Throat, shear order at least 10%", "#42396d", "-.", [r["ordered"] + r["upper_layer"] + r["lower_layer"] for r in rows]),
        ]
        for label, color, style, fraction in series:
            ax.plot(times, fraction, color=color, ls=style, lw=2.25, label=label)
        ax.set(xlim=(0, 12), xticks=(0, 3, 6, 9, 12), ylim=(-.035, 1.035), xlabel=r"Time ($R/c_T$)", title=title)
        ax.grid(axis="y", color="#d9dfe0", lw=.7)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].yaxis.set_major_formatter(PercentFormatter(1))
    axes[0].set_ylabel("Fraction of initial release energy")
    axes[1].legend(loc="center right", fontsize=8, frameon=True, facecolor="white", framealpha=.95)
    fig.suptitle("The same displacement released inside the throat", fontsize=14, weight="semibold")
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=185)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", type=Path, nargs="+", help="Plot the first folder and summarize every folder")
    parser.add_argument("--figure", type=Path, required=True)
    args = parser.parse_args()
    for i, folder in enumerate(args.folders):
        data = collect(folder)
        if i == 0:
            plot(data, args.figure)
