"""Reduce archived 0PN push, stream and tow outputs. No field solver is called."""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MEASUREMENTS = ROOT / "research/pde_ledger_v3/_measurements"


def archived_rows(record):
    if "literal_stdout" in record:
        body = record["literal_stdout"]
        exits = record["guard_exit"], record["child_exit"]
    else:
        body = "\n".join(record["stdout_lines"]) + "\n"
        exits = record["exit_code"], record["worker_exit"]
    digest = hashlib.sha256(body.encode()).hexdigest()
    if (digest != record["stdout_sha256"] or exits != (0, 0)
            or record["limits_verified"] is not True or record["swap_cap_bytes"] != 0):
        raise ValueError(f"Unverified guarded input for {record['command']}")
    return [json.loads(line) for line in body.splitlines()]


def rigid_response(rows):
    drives = {}
    for r in rows:
        if r["test"] == "free_response":
            w = r["omega"]
            a = complex(*r["acceleration"])
            if abs(a - 1j * w * complex(*r["velocity"])) > 1e-12:
                raise ValueError("Velocity and acceleration disagree")
            pair = drives.setdefault(w, {})
            if r["drive"] in pair:
                raise ValueError("Duplicate drive")
            pair[r["drive"]] = a
    result = []
    for w, pair in sorted(drives.items()):
        push, stream = pair["push"], pair["incident"]
        load = stream / push
        result.append(dict(omega=w, a_per_F=[push.real, push.imag],
                           a_per_U=[stream.real, stream.imag],
                           actuator_resistance=(1j * w / push).real,
                           minus_H_per_U=[load.real, load.imag],
                           stream_quadrature_per_omega=load.imag / w))
    return result


def make_figure(data, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.8), layout="constrained")
    variables = ("actuator_resistance", "a_per_U", "stream_quadrature_per_omega")
    labels = (r"$\mathrm{Re}(F/V)$", r"$\mathrm{Im}(a/U)$",
              r"$\mathrm{Im}(-H/U)\,/\,\omega$")
    for ax, key, ylabel in zip(axes, variables, labels):
        for case, color, marker in zip(data, ("#086f78", "#a23e27", "#464a65"), ("o", "+", "x")):
            points = case["responses"]
            values = [p[key][1] if key == "a_per_U" else p[key] for p in points]
            ax.plot([p["omega"] for p in points], values, color=color, marker=marker,
                    markersize=8, markerfacecolor="none", linestyle="none",
                    label=f"mesh {case['mesh']}, intake {case['intake_nominal']:g}")
        ax.set(xlabel=r"Angular frequency $\omega$", ylabel=ylabel, xticks=(.01, .02, .04))
        fmt = ScalarFormatter(useOffset=False, useMathText=True)
        fmt.set_powerlimits((-3, 3))
        ax.yaxis.set_major_formatter(fmt)
        ax.grid(alpha=.23)
    axes[0].legend(frameon=False, fontsize=8)
    for ax, title in zip(axes, ("Push in resting baths", "Free response in stream", "Reconstructed held load")):
        ax.set_title(title, fontsize=11)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()
    src = json.loads((MEASUREMENTS / "mouth_carrier_0pn.json").read_text())["runs"]
    occ = json.loads((MEASUREMENTS / "mouth_occupied_slow_push.json").read_text())["guarded_runs"]
    hashes = {}

    def read(archive, name):
        rows = archived_rows(archive[name])
        hashes[name] = archive[name]["stdout_sha256"]
        return rows

    rigid = []
    steady = []
    for name in ("carrier_free_c26_01", "carrier_free_c34_01", "carrier_free_rest_01"):
        rows = read(src, name)
        prep = rows[0]
        rigid.append(dict(run=name, mesh=prep["mesh"], far=prep["far"], roof=prep["roof"],
                          intake_nominal=0.25 if prep["rate"] else 0, intake=prep["rate"],
                          responses=rigid_response(rows)))
    for suffix in ("c26_01", "c34_01", "c48_01", "ext14_01", "ext20_01",
                   "restmatched26_02", "still_01"):
        name = "carrier_inertia_" + suffix
        rows = read(src, name)
        prep = rows[0]
        tow = next(r for r in rows if r.get("drive") == "tow")
        steady.append(dict(run=name, mesh=prep["mesh"], far=prep["far"], roof=prep["roof"],
                           intake=prep["rate"], actuator_wall=-tow["direct_wall"],
                           actuator_outer_cut=-tow["outer_cut"],
                           actuator_interior=[-f for f in tow["interior"]]))

    occupied = []
    for suffix in ("base26_01", "base34_01", "large13p5_01", "large18_01"):
        name = "occupied_push_" + suffix
        rows = read(occ, name)
        if rows[0]["test"] == "preparation":
            n = rows[0]["mesh"]
        else:
            bath = Path(occ[name]["command"][-1]).name
            n = read(occ, bath)[0]["n"]
        setup = next(r for r in rows if r["test"] == "setup")
        points = []
        for r in rows:
            if r["test"] == "push":
                a = r["acceleration_per_force"]
                chi = complex(a["real"][0], a["imag"][0])
                points.append(dict(tension=r["tension"], omega=r["omega"],
                                   a_per_F=[chi.real, chi.imag],
                                   actuator_resistance=(1j * r["omega"] / chi).real))
        occupied.append(dict(run=name, mesh=n, far=setup["far"], roof=setup["roof"],
                             wave_action=setup["wave_action"], responses=points))
    if [len(c["responses"]) for c in rigid] != [3, 2, 2] or any(
            len(c["responses"]) != 12 for c in occupied):
        raise ValueError("Incomplete archived drive coverage")
    if args.figure:
        make_figure(rigid, args.figure)
    print(json.dumps(dict(convention="exp(+i omega t), external actuator force",
                          rigid=rigid, steady_tow=steady, occupied=occupied,
                          input_stdout_sha256=hashes), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
