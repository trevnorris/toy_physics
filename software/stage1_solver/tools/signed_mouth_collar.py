"""Run the short wall/wave segment with the archived radius-one rim drives."""

import argparse
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from stage1_solver.mouth_collar import fit_terminal_load, solve_segment


ROOT = Path(__file__).resolve().parents[3]
RUNS = ROOT / "_scratch/s11c"
DRIVES = ("negative", "negative_half", "zero", "positive_half", "positive")
BASE_FIELD = RUNS / "mouth_feed_depth01/fields/n32_restoring_slope_u0.npz"


def load_drives():
    provenance = (ROOT/"docs/_measurements/mouth_prescribed_acceleration.txt").read_text()
    def verified(path):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if f"{digest}  {path.relative_to(ROOT)}" not in provenance:
            raise ValueError(f"Archive provenance mismatch for {path}")
        return digest
    baseline = verified(BASE_FIELD)
    drives = {}
    for name in DRIVES:
        directory = RUNS/f"mouth_prescribed_accel_{name}02"
        verified(directory/"outcome.json")
        verified(directory/"stdout")
        row = json.loads((directory/"stdout").read_text())
        outcome = json.loads((directory/"outcome.json").read_text())
        if not outcome["limitsVerified"] or outcome["exitCode"] or outcome["childOutcome"]["exitCode"]:
            raise ValueError("Invalid source drive")
        sample = row["32"]
        times = np.array([p["t"] for p in sample])
        a = row["acceleration"]
        np.testing.assert_allclose([p["x"] for p in sample], a*times**2/2, atol=1e-22,rtol=1e-12)
        np.testing.assert_allclose([p["inward_rate"] for p in sample],0.25,atol=1e-13,rtol=0)
        drives[name] = a, times
    return baseline, drives


def reference_calibration():
    """Independent velocity and traction traces of two outgoing elastic modes."""
    angle = 0.3
    rot = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    speed = np.array([1.2,0.8])
    modulus = np.array([1.08,0.4])
    def readings(omega, phase):
        t = np.linspace(0,4,37)[:,None]
        def field(z):
            return np.sin((t-z/speed)*omega+phase)
        h = 1e-5
        v = omega*np.cos(t*omega+phase)
        stress = -modulus*(field(h)-field(-h))/(2*h)
        return v@rot.T, stress@rot.T
    train = readings(np.array([1.3,2.4]),0.25)
    held = readings(np.array([2.1,0.8]),1.1)
    return fit_terminal_load(*train,*held,1e-7)


def demo_profiles(z):
    # Known standing-mode amplitude on the benchmark patch, transverse norm I=1.
    amplitude = np.cos(np.pi*z/(2*1.85))
    weight = amplitude**2/2
    return (np.array([np.ones_like(z),2*weight]),
            np.array([np.ones_like(z),1.125*weight]), 2*amplitude**2)


def record(path, readings, amplitude, tag):
    np.savez_compressed(path, acceleration=amplitude, scope=tag, **readings)
    summary = {k:readings[k][:,-1].tolist() for k in ("stored","cut_flux","cut_impulse")}
    summary.update(wall_wave_cut_positions=readings["displacement"][:,
                   np.searchsorted(readings["z"],readings["cut_depth"]),-1].tolist(),
                   maximum_balance_error=float(abs(readings["balance_residual"]).max()))
    print(json.dumps(dict(case=path.stem,t=readings["times"][-1],**summary),allow_nan=False),flush=True)


def run_set(directory, tag, drives, z, cut, profiles, load, extra_times=(), body_force=None):
    results = {}
    for name,(a,t) in drives.items():
        times = np.union1d(t,extra_times)
        forcing = None if body_force is None else lambda time: body_force(name,time)
        row = solve_segment(z,cut,*profiles,load,a,times,distributed_force=forcing)
        record(directory/f"{tag}_{name}.npz",row,a,tag)
        scale = max(abs(row["stored"]).max(), abs(row["cut_impulse"]).max(),
                    abs(row["rim_impulse"]).max(),abs(row["source_impulse"]).max(),1e-28)
        if abs(row["balance_residual"]).max() > 2e-9*scale:
            raise ArithmeticError("Impulse balance failed")
        np.testing.assert_allclose(row["displacement"][0,0], a*times**2/2,rtol=1e-9,atol=1e-20)
        results[name] = row
    # Unforced, reversed and held-out half-strength readings in the linear fixture.
    if body_force is None:
        for name,multiple in (("zero",0),("negative",-1),("positive_half",0.5),("negative_half",-0.5)):
            for key in ("stored","cut_impulse","cut_flux","displacement"):
                primary = results["positive"][key]
                allowed = 3e-7*abs(primary).max()+5e-15
                if abs(results[name][key]-multiple*primary).max() > allowed:
                    raise ArithmeticError(f"{name} scaling failed for {key}")
    return results["positive"]


def demo(directory, drives):
    load, errors = reference_calibration()
    print(json.dumps(dict(calibration="outgoing-mode benchmark",load=load.tolist(),
                          fit_and_holdout_relative_error=errors)), flush=True)
    cut_z, end = 1.85/4, 1.85*5/16
    for window in ("archive","continued"):
        full = []
        for n in (64,128):
            if window == "archive":
                # Resolve the short wavefront from the saved 0.001 drive history.
                near = np.r_[np.linspace(0,0.002,n+1)[:-1],np.linspace(0.002,cut_z,n+1)]
                extra = ()
            else:
                near = np.linspace(0,cut_z,n+1)
                extra = np.linspace(0.005,1,200)
            z = np.r_[near,np.linspace(cut_z,end,n//4+1)[1:]]
            cut = len(near)-1
            row = run_set(directory,f"benchmark_{window}_{n}",drives,z,cut,
                          demo_profiles(z),load,extra)
            full.append(row)
        # Report absolute changes too, particularly for the causally unreached cut.
        for key in ("stored","cut_flux","cut_impulse"):
            coarse,fine = (r[key][:,-1] for r in full)
            delta = abs(fine-coarse)
            print(json.dumps(dict(window=window,refinement=key,absolute_change=delta.tolist())),flush=True)
            if window == "continued":
                np.testing.assert_allclose(coarse,fine,rtol=0.025,atol=5e-10)
    # Independent continuum control. Matched wall load, no coupling, no wave drive.
    mu,tension,coupling = demo_profiles(z)
    row = solve_segment(z,cut,mu,tension,0*coupling,np.diag([1.,0.]),0.02,np.array([1.]))
    travel = max(0,1-cut_z)
    actual = np.array([row["stored"][0,-1],row["cut_flux"][0,-1],row["cut_impulse"][0,-1]])
    expected = 0.02*np.array([(1-travel**2)/2,travel,travel**2/2])
    print(json.dumps(dict(outgoing_wall=actual.tolist(),continuum=expected.tolist(),
                          undriven_wave=float(abs(row["displacement"][1]).max()))),flush=True)
    np.testing.assert_allclose(actual,expected,rtol=0.002,atol=1e-12)
    if np.any(row["displacement"][1]):
        raise ArithmeticError("Unforced independent wave moved")


def supplied(directory, drives, data, baseline):
    if data["baseline_sha256"] != baseline or data["inward_rate"] != 0.25:
        raise ValueError("Material packet belongs to a different drain")
    load, errors = fit_terminal_load(**data["terminal_readings"])
    print(json.dumps(dict(calibration="supplied terminal readings",load=load.tolist(),
                          fit_and_holdout_relative_error=errors)),flush=True)
    z = np.asarray(data["z"],dtype=float)
    amplitude = np.asarray(data["standing_amplitude"],dtype=float)
    if amplitude.shape != z.shape or not np.isfinite(amplitude).all():
        raise ValueError("Supply the standing profile at every depth")
    norm = np.asarray(data["transverse_gradient_norm"],dtype=float)
    if norm.shape not in ((),z.shape) or not np.isfinite(norm).all() or np.any(norm < 0):
        raise ValueError("Invalid transverse gradient normalization")
    weight = norm*amplitude**2/2
    mass = np.stack([np.broadcast_to(data["wall_line_inertia"],z.shape),
                     weight*np.asarray(data["shear_inertia"])])
    stiffness = np.stack([np.broadcast_to(data["wall_tension"],z.shape),
                          weight*np.asarray(data["shear_modulus"])])
    coupling = np.broadcast_to(data["coupling_per_length"],z.shape)
    # Explicit null declares zero incremental distributed traction in this model.
    histories = data["distributed_force"]
    force = None
    if histories is not None:
        traces = {}
        for name,(a,t) in drives.items():
            samples = np.asarray(histories[name]["times"],dtype=float)
            values = np.asarray(histories[name]["values"],dtype=float)
            if (samples.ndim != 1 or len(samples)<2 or np.any(np.diff(samples)<=0)
                    or samples[0]>0 or samples[-1]<t[-1] or values.shape!=(len(samples),2,len(z))
                    or not np.isfinite(samples).all() or not np.isfinite(values).all()):
                raise ValueError("Force histories must cover each drive on the same grid")
            traces[name] = samples,values
        def force(name,time):
            samples,values = traces[name]
            j = np.searchsorted(samples,time,side="right").clip(1,len(samples)-1)
            u = (time-samples[j-1])/(samples[j]-samples[j-1])
            return (1-u)*values[j-1]+u*values[j]
    run_set(directory,"supplied",drives,z,data["cut_index"],(mass,stiffness,coupling),load,
            body_force=force)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--demo",action="store_true")
    group.add_argument("--material",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    baseline,drives = load_drives()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=args.output.parent) as temp:
        work = Path(temp)/"results"
        work.mkdir()
        if args.demo:
            demo(work,drives)
        else:
            supplied(work,drives,json.loads(args.material.read_text()),baseline)
        work.rename(args.output)


if __name__ == "__main__":
    main()
