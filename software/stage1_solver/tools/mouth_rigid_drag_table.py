"""Reduce the finite-box coherent-carrier push runs and plot separate stream controls.

Run at the repository root. --collect assembles the durable record from original
guarded outputs. Given this record, the table and plot can be regenerated without
scratch fields or any solver. Frequencies are angular and all quantities are in
the nondimensional carrier units.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


SCRATCH = Path("_scratch/s11c")
SLOW = SCRATCH / "rigid_energy_slow01"
TOOLS = Path("software/stage1_solver/tools")
MEASURE = Path("research/pde_ledger_v3/_measurements")
SOURCE = TOOLS / "mouth_rigid_push_energy.py"
RECORD = MEASURE / "mouth_rigid_drag_table.json"
PLOT = Path("software/stage1_solver/reports/figures/mouth_rigid_drag_table.png")
PREPARATIONS = {
    "n26": (26, 9, 10, "mouth_flared_cubic_n26"),
    "n34": (34, 9, 10, "mouth_flared_cubic_n34"),
    "n48": (48, 9, 10, "mouth_flared_cubic_n48"),
    "r13p5": (34, 13.5, 15, "occupied_push_bath13p5_n34_01"),
    "r18": (34, 18, 20, "occupied_push_bath18_n34_01"),
}


def load(path):
    return json.loads(path.read_text())


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def close(x, y):
    return math.isclose(x, y, rel_tol=1e-10, abs_tol=1e-17)


def collect_case(folder, prep, frequency, expected_hash):
    mesh, far, roof, source_dir = PREPARATIONS[prep]
    invocation = load(folder / "invocation.json")
    outcome = load(folder / "outcome.json")
    child = load(folder / "child-outcome.json")
    limits = load(folder / "limit-validation.json")
    duration = load(folder / "duration-validation.json")
    command = invocation["command"]
    body = (folder / "stdout").read_text()
    rows = [json.loads(line) for line in body.splitlines()]
    require(outcome["exitCode"] == child["exitCode"] == 0 and child["guardReason"] is None
            and child == outcome["childOutcome"] and outcome["limitsVerified"]
            and duration["verified"] and duration["actual"]["RuntimeMaxUSec"] == "infinity"
            and limits["verified"] and limits["actual"]["memory.swap.max"] == "0"
            and int(limits["actual"]["memory.max"]) == invocation["memoryMax"]
            and invocation["unit"] == outcome["unit"] and invocation["swapMax"] == 0,
            f"Invalid guard at {folder}")
    source = SCRATCH / source_dir / "field.npz"
    require(command[1] == str(SOURCE) and Path(command[3]) == source.parent
            and sha(source) == rows[0]["sha256"] == expected_hash
            and rows[0]["far"] == far and rows[0]["roof"] == roof, f"Preparation at {folder}")
    saved = Path(command[command.index("--save-response") + 1])
    push, boundary, a, b = rows[1:5]
    require(push["test"] == "unit_velocity" and push["dynamic_residual"] < 2e-8
            and push["omega"] == frequency == float(command[command.index("--omega") + 1])
            and a["test"] == b["test"] == "enclosing_outward_power_per_speed_squared"
            and a["window"] != b["window"], f"Invalid readout at {folder}")
    with np.load(saved, allow_pickle=False) as data:
        require(str(data["background_sha256"]) == expected_hash
                and float(data["omega"]) == frequency
                and complex(data["impedance"]) == complex(*push["actuator_impedance"])
                and data["harmonic_real"].shape == data["harmonic_imag"].shape == (push["nodes"],),
                f"Response field mismatch at {saved}")
    return dict(preparation=prep, mesh=mesh, far=far, roof=roof, omega=frequency,
                command=command, source_field_sha256=expected_hash, response_field_sha256=sha(saved),
                guard_exit=outcome["exitCode"], worker_exit=child["exitCode"], limits_verified=True,
                memory_cap_bytes=invocation["memoryMax"], swap_cap_bytes=invocation["swapMax"],
                stdout_sha256=sha(folder / "stdout"), literal_stdout=body)


def collect():
    design = load(SLOW / "design.json")
    done = load(SLOW / "DONE.json")
    # This table describes the original direct-wall runs, whose pinned source
    # remains with the saved batch even as the active force calculation evolves.
    require(design["source_sha256"] == sha(SLOW / "instrument.py")
            == "18847c31a68c3b5e7850dffe72cf276c3c915e4b9d8e192d9bef08804376f47f",
            "Pinned instrument differs")
    inputs = [(SLOW / (k + "_" + w[1:]), k, float(w)) for k, w in design["cases"]]
    require(len(inputs) == len(done) == 17 and all(r["exit"] == 0 for r in done)
            and {r["name"] for r in done} == {p.name for p, _, _ in inputs}, "Incomplete slow batch")
    for k in PREPARATIONS:
        suffix = "001" if k == "n34" else "control01"
        inputs.append((SCRATCH / f"rigid_push_energy_{k}_{suffix}", k, .01))
    field_hashes = load(MEASURE / "mouth_order_and_flare.json")["output_field_sha256"]
    occupied = load(MEASURE / "mouth_occupied_slow_push.json")["guarded_runs"]
    expected = {}
    for k, (mesh, far, roof, name) in PREPARATIONS.items():
        if name in occupied:
            rec = occupied[name]
            require(rec["stdout_sha256"] == hashlib.sha256(rec["literal_stdout"].encode()).hexdigest()
                    and rec["guard_exit"] == rec["child_exit"] == 0 and rec["limits_verified"]
                    and rec["swap_cap_bytes"] == 0, f"Nested preparation {name}")
            state = json.loads(rec["literal_stdout"])
            require(state["n"] == mesh and (state["far"], state["roof"]) == (far, roof)
                    and abs(state["bulk"]-.25) < 3e-9, f"Nested state {name}")
            expected[k] = state["field_sha256"]
        else:
            expected[k] = field_hashes[str(SCRATCH / name / "field.npz")]
    runs = {str(p):collect_case(p, k, w, expected[k]) for p,k,w in inputs}
    stream_archive = load(MEASURE / "mouth_drag_response_audit.json")
    raw = (json.dumps(stream_archive["stdout"], indent=2, allow_nan=False) + "\n").encode()
    require(hashlib.sha256(raw).hexdigest() == stream_archive["stdout_sha256"]
            and stream_archive["worker_exit"] == stream_archive["guard_exit"] == 0
            and stream_archive["verified"] and stream_archive["swap_cap_bytes"] == 0,
            "Stream archive verification")
    return dict(convention="exp(+i omega t), positive applied force on mouth",
                per_unit="outward cycle-mean excess power divided by |V|²; twice this has F/V units",
                source=str(SOURCE), source_sha256=design["source_sha256"], runs=runs,
                stream_archive="mouth_drag_response_audit.json", stream_stdout_sha256=stream_archive["stdout_sha256"],
                stream=stream_archive["stdout"]["rigid"])


def reduce(record):
    rows = []
    for folder, run in record["runs"].items():
        require(hashlib.sha256(run["literal_stdout"].encode()).hexdigest() == run["stdout_sha256"], folder)
        prep, push, boundary, *other = map(json.loads, run["literal_stdout"].splitlines())
        w = push["omega"]
        z = complex(*push["actuator_impedance"])
        chi = 1j*w/z
        primary = w * (1/chi).imag
        fw = boundary["remote"]["feed"]["total"]
        bw = sum(port["total"] for key,port in boundary["remote"].items() if key != "feed")
        require(close(-primary, z.real) and close(-primary/2,push["actuator_mean_power_per_speed_squared"])
                and close(fw+bw, boundary["remote_total"]), folder)
        near, far = other[:2]
        for cut in (near,far):
            require(close(cut["parts"]["feed"]["total"]+cut["parts"]["receiver"]["total"],cut["total"]), folder)
        signed = {r["applied_force"]:r for r in other[2:]}
        amp = max(signed)
        require(set(signed) == {0,amp,-amp,amp/2,-amp/2}, folder)
        for force,r in signed.items():
            velocity = complex(*r["velocity"])
            require(close(abs(velocity-force/z),0) and close(r["actuator"],.5*force*velocity.real)
                    and close(r["remote"],abs(velocity)**2*(fw+bw)), folder)
        rows.append(dict(mesh=run["mesh"], far=prep["far"], roof=prep["roof"], omega=w,
                         chi_a=[chi.real,chi.imag], omega_im_inverse_chi=primary,
                         feed_resistance=2*fw, bulk_resistance=2*bw,
                         net_port_resistance=2*(fw+bw), window1=2*near["total"], window2=2*far["total"]))
    return sorted(rows,key=lambda r:(r["mesh"],r["far"],r["omega"]))


def draw(rows, stream):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

    fig, axes = plt.subplots(1,3,figsize=(12.5,4.4),layout="constrained")
    colors = ("#216788","#ac5224","#267662")
    for i, (far,roof) in enumerate(((9,10),(13.5,15),(18,20))):
        sub = [r for r in rows if r["mesh"] == 34 and r["far"] == far]
        x = [r["omega"] for r in sub]
        style = dict(color=colors[i],marker=("o","s","D")[i],ms=5,lw=1.2)
        axes[0].plot(x,[1e6*r["omega_im_inverse_chi"] for r in sub],label=f"{far:g} / {roof:g}",**style)
        axes[1].plot(x,[1e6*r["net_port_resistance"] for r in sub],**style)
    for u in stream:
        r = u["responses"]
        axes[2].plot([k["omega"] for k in r],[1e3*k["a_per_U"][1] for k in r],
                     color={18:"#575963",26:"#aa5729",34:"#216788"}[u["mesh"]],
                     marker={18:"x",26:"o",34:"+"}[u["mesh"]],markersize=7,markerfacecolor="none",
                     linestyle="none",label=f"mesh {u['mesh']}, intake {u['intake_nominal']:g}")
    axes[0].set_ylabel(r"$\omega\,\mathrm{Im}(1/\chi_a)\quad[10^{-6}]$")
    axes[1].set_ylabel(r"$2\,\overline{P}_{\rm out}/|V|^2\quad[10^{-6}]$")
    axes[2].set_ylabel(r"$\mathrm{Im}(a/U_s)\quad[10^{-3}]$")
    axes[0].legend(title="Bath radius / depth",fontsize=8,frameon=False)
    axes[1].axhline(0,lw=.8,ls="--",color=".45",zorder=0)
    axes[2].legend(fontsize=8,frameon=False)
    for ax,title in zip(axes,("Push susceptibility • mesh 34","Power into both reservoirs","Separate stream drive • bath 9 / 10")):
        ax.set_title(title,fontsize=11,pad=10)
        ax.set_xlabel(r"Angular frequency $\omega$")
        ax.grid(alpha=.18)
    for ax in axes[:2]:
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(FixedLocator((.000225,.001,.005,.02)))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:f"{v:g}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
    axes[2].set_xticks((.01,.02,.04))
    fig.savefig(PLOT,dpi=170)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--collect",action="store_true",help="Build the record from guarded source jobs")
    args = p.parse_args()
    record = collect() if args.collect else load(RECORD)
    table = reduce(record)
    draw(table,record["stream"])
    if args.collect:
        record.update(reducer=str(Path(__file__).relative_to(Path.cwd())), reducer_sha256=sha(Path(__file__)),
                      table=table, figure_sha256=sha(PLOT))
        RECORD.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    else:
        require(record["table"] == table,"Saved table differs from the source records")
    print("mesh bath omega             push_q       feed       bulk        net    window1    window2  (all readouts ×10⁻⁶)")
    for r in table:
        values = [r[k]*1e6 for k in ("omega_im_inverse_chi","feed_resistance","bulk_resistance",
                                      "net_port_resistance","window1","window2")]
        print(f"{r['mesh']:>4} {r['far']:>4g} {r['omega']:>8g} " + " ".join(f"{a:+10.5f}" for a in values))


if __name__ == "__main__":
    main()
