"""Validate the signed strip preparations and display their through-thickness stress."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "_scratch/s11c/ordered_strip_bend01"
RECORD = ROOT / "research/pde_ledger_v3/_measurements/mouth_ordered_strip_bend.json"
FIGURE = ROOT / "software/stage1_solver/reports/figures/mouth_ordered_strip_stress.png"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked(name):
    directory = BASE / name
    read = lambda f: json.loads((directory / f).read_text())
    out, child = read("outcome.json"), read("child-outcome.json")
    inv, lim, dur = read("invocation.json"), read("limit-validation.json"), read("duration-validation.json")
    if (out["exitCode"] or child["exitCode"] or child["guardReason"] or not out["limitsVerified"]
            or not lim["verified"] or lim["actual"]["memory.swap.max"] != "0"
            or not dur["verified"] or inv["wallDeadlineSeconds"] is not None):
        raise ValueError(f"Unverified experiment: {directory}")
    path = ROOT / inv["command"][inv["command"].index("--output") + 1]
    return dict(command=inv["command"], start_utc=inv["startedUtc"], verified_exit=out["exitCode"],
                memory_cap_bytes=inv["memoryMax"], swap_cap_bytes=inv["swapMax"],
                stdout=(directory / "stdout").read_text().splitlines(), data_sha256=digest(path),
                results=json.loads(path.read_bytes()))


def main():
    runs = {case: checked(case) for case in ("main", "resolution")}
    primary = runs["main"]["results"]
    fig, axes = plt.subplots(2, 1, figsize=(7.7, 5.8), sharex=True, layout="constrained")
    profile_checks = []
    for panel, key, title in zip(axes, ("native", "registered_trial"),
                                  ("Current strip equations", "With registered crystalline layers (trial)")):
        data = primary[key]
        plus, minus = data["training"][1], data["training"][0]
        moment = 0.
        probes = []
        for pp, nn in zip(plus["stress_profiles"], minus["stress_profiles"]):
            z = np.array(pp["depth"])
            if pp["depth"] != nn["depth"]:
                raise ValueError("Mismatched depth at opposite bends")
            order = (np.array(pp["order_axial_stress"]) - nn["order_axial_stress"]) / 2.
            crystal = (np.array(pp["crystal_axial_stress"]) - nn["crystal_axial_stress"]) / 2.
            total = order + crystal
            moment += simpson(z * total, x=z)
            first = pp["region"] == 0
            panel.plot(z, total, color="#19446a", lw=2., label="Total" if first else None)
            if key == "registered_trial":
                panel.plot(z, crystal, color="#bd501e", ls="--", lw=1.7,
                           label="Crystalline contribution" if first else None)
            if pp["region"] == 1:
                for location in (-.25, .25):
                    probes.append(dict(depth=location, odd_order=float(np.interp(location, z, order)),
                                       odd_crystal=float(np.interp(location, z, crystal)),
                                       odd_total=float(np.interp(location, z, total))))
        stress_integral = (plus["holding_moment"] - minus["holding_moment"]) / 2.
        profile_checks.append(dict(law=key, interior_samples=probes, plotted_moment=moment,
                                   full_quadrature_moment=stress_integral))
        panel.set_title(title, fontsize=12, loc="left")
        panel.set_ylabel("Signed axial stress")
        panel.axhline(0., color=".55", lw=.7, zorder=0)
        for p in (-.5, .5):
            panel.axvline(p, color=".55", lw=.95, ls=":", zorder=0)
        if key == "registered_trial":
            panel.legend(loc="upper left", frameon=False)
        panel.spines[["top", "right"]].set_visible(False)
    axes[1].set_xlabel("Depth across the strip, in units of its thickness")
    axes[1].set_xlim(-.8, .8)
    fig.suptitle("Stress that reverses under opposite bends", fontsize=15)
    fig.savefig(FIGURE, dpi=160)
    package = dict(input_checkpoint="25b20c94", runs=runs, profile_checks=profile_checks,
                   figure_path=str(FIGURE.relative_to(ROOT)), figure_sha256=digest(FIGURE),
                   sources={str(p.relative_to(ROOT)):digest(p) for p in
                            (Path(__file__).resolve(), Path(__file__).with_name("mouth_ordered_strip_bend.py"),
                             Path(__file__).with_name("mouth_curved_order_patch.py"))})
    RECORD.write_text(json.dumps(package, indent=2, allow_nan=False) + "\n")
    print(json.dumps(profile_checks))
    for row in profile_checks:
        if abs(row["plotted_moment"] - row["full_quadrature_moment"]) > 5.e-6:
            raise ValueError("Figure integral disagrees with measured moment")


if __name__ == "__main__":
    main()
