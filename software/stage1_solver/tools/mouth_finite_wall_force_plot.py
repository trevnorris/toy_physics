"""Plot finite-box corrected force and port powers from the preserved audit."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
RECORD = ROOT / "research/pde_ledger_v3/_measurements/mouth_finite_wall_force_audit.json"
FIGURE = ROOT / "software/stage1_solver/reports/figures/mouth_finite_wall_force.png"


def rows(record, run, test):
    data = map(json.loads, record["runs"][run]["literal_stdout"].splitlines())
    return [r for r in data if r.get("test") == test]


def main():
    record = json.loads(RECORD.read_text())
    push = rows(record, "corrected", "corrected_push")
    series = rows(record, "scaling", "finite_box_scaling")
    slowest = min(series, key=lambda row: row["omega"])
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.8, 4.3), layout="constrained")
    colors = ("#005a70", "#bf650d", "#70428f")
    sizes = ((9, 10), (13.5, 15), (18, 20))
    box_labels = []
    for (radius, depth), color in zip(sizes, colors):
        sample = sorted((r for r in push if r["mesh"] == 34 and r["far"] == radius),
                        key=lambda row: row["omega"])
        w = [r["omega"] for r in sample]
        if len(w) != 6 or len(set(w)) != 6:
            raise ValueError(f"Incomplete frequency sweep for radius {radius}")
        q = np.array([r["signed_push"] for r in sample])
        # Windows contain 2P_out/|V|². In the exp(+iωt) convention its
        # negative can be compared directly to ω Im(1/χ_a).
        check = -np.array([r["independent_windows"] for r in sample])
        np.testing.assert_allclose(q, [-np.real(complex(*r["impedance"])) for r in sample],
                                   rtol=1e-9, atol=1e-15)
        line, = a.plot(w, 1e6*q, color=color, marker="o", markersize=6,
                       markerfacecolor="white", label=f"{radius:g} / {depth:g}")
        a.plot(w, 1e6*check[:, 0], color=color, marker="x", linestyle="none", markersize=5)
        a.plot(w, 1e6*check[:, 1], color=color, marker="+", linestyle="none", markersize=5)
        box_labels.append(line)
    keys = [Line2D([], [], color="#303030", marker=m, markerfacecolor="white",
                   linestyle="none", label=name) for m, name in
            (("o", "Push"), ("x", "Window 1"), ("+", "Window 2"))]
    box_legend = a.legend(handles=box_labels, title="Radius / depth", loc="center left",
                          frameon=False, fontsize=9)
    a.add_artist(box_legend)
    a.legend(handles=keys, loc="center right", frameon=False, fontsize=9)
    a.axhline(0, color="#666666", lw=.7, ls="--")
    a.set(xscale="log", xlabel=r"Angular frequency $\omega$",
          ylabel=r"$\omega\,\mathrm{Im}(1/\chi_a)$ or $-2\overline{P}_{\rm window}/|V|^2$  [$10^{-6}$]",
          title="Corrected push and two enclosing windows")
    a.xaxis.set_major_locator(FixedLocator([.000225, .001, .005, .02]))
    a.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    a.xaxis.set_minor_formatter(NullFormatter())

    radii = np.array([s[0] for s in sizes])
    for key, power, name, marker, color in (
            ("feed_resistance", 3, "Brane (outward)", "o", "#216e69"),
            ("bulk_resistance", 4, "Bulk (inward)", "s", "#b34727")):
        values = 1e6*np.abs(slowest[key]["measured"])
        guide = values[0]*(radii/radii[0])**(-power)
        b.plot(radii, values, color=color, marker=marker, ms=7, ls="none", label=name)
        b.plot(radii, guide, color=color, ls="--", lw=1.2,
               label=rf"$R^{{-{power}}}$ through first point")
    b.set(xscale="log", yscale="log", xlabel=r"Reservoir radius $R$",
          ylabel=r"$2|\overline{P}_{\rm port}|/|V|^2$  [$10^{-6}$]",
          title=f"Port magnitudes at ω = {slowest['omega']:g}")
    b.xaxis.set_major_locator(FixedLocator(radii))
    b.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    b.xaxis.set_minor_formatter(NullFormatter())
    b.yaxis.set_major_locator(FixedLocator([1, 2, 4, 8, 16]))
    b.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    b.yaxis.set_minor_formatter(NullFormatter())
    b.legend(loc="upper right", fontsize=8, frameon=False)
    for axis in (a, b):
        axis.grid(alpha=.2)
        axis.set_axisbelow(True)
        axis.title.set_fontsize(11)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
