#!/usr/bin/env python3
"""
stepE_sweep_convergence.py

Automated, resumable sweep runner for:
  stepE_euler_throat_longitudinal.py

Designed to run unattended on a cloud GPU instance. It:
- builds a parameter grid (preset or explicit lists),
- runs StepE once per grid point (each run in its own subprocess to avoid GPU memory-pool carryover),
- parses StepE's JSON diagnostics in real-time,
- prints compact live progress lines while each run is evolving,
- writes one JSON summary row per run to a JSONL file (append-only),
- optionally saves full raw stdout logs per run.

This is the “Next Step 1 — Convergence & robustness sweep” runner described in poisson_01.md.

Example (GPU, small preset):
  python stepE_sweep_convergence.py --preset small --backend cupy --dtype float32 \
    --steps 10000 --diag_every 400 --out stepE_sweep_small.jsonl

Example (custom lists):
  python stepE_sweep_convergence.py --backend cupy --dtype float32 \
    --N_list 256,384,512 --L_list 150,200,300 --sigma_list 1,2,3 \
    --gamma_drag_list 0,0.01,0.03 --ramp_time_list 0,2 \
    --steps 12000 --diag_every 600 --out stepE_sweep.jsonl

Resume behavior:
- By default, if a run_id already exists in the output JSONL, that run is skipped.
- Use --rerun to force re-running everything.

Notes:
- For big N (>=512), float32 is recommended unless you have a very large GPU.
- You can reduce per-run stdout volume by keeping --sample_points 0 (default here).
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple


# ----------------------------
# Small utilities
# ----------------------------

def _parse_csv_list(s: Optional[str], *, cast) -> Optional[List[Any]]:
    """Parse 'a,b,c' into a list with casting; returns None if s is None."""
    if s is None:
        return None
    s = s.strip()
    if not s:
        return None
    out: List[Any] = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        out.append(cast(part))
    return out


def _stable_json_dumps(x: Any) -> str:
    return json.dumps(x, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def make_run_id(params: Dict[str, Any]) -> str:
    """Stable run id from params (excluding non-physics keys)."""
    payload = _stable_json_dumps(params).encode("utf-8")
    return hashlib.sha1(payload).hexdigest()[:12]


def load_existing_run_ids(path: Path) -> Set[str]:
    """Return set of run_id already present in a JSONL file."""
    if not path.exists():
        return set()
    done: Set[str] = set()
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except Exception:
                continue
            rid = obj.get("run_id")
            if isinstance(rid, str):
                done.add(rid)
    return done


def safe_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def fmt(x: Optional[float], *, width: int = 8, prec: int = 3) -> str:
    if x is None:
        return " " * width
    try:
        return f"{x:{width}.{prec}g}"
    except Exception:
        return str(x)[:width].rjust(width)


# ----------------------------
# Presets
# ----------------------------

PRESETS: Dict[str, Dict[str, List[Any]]] = {
    # Very small "sanity sweep" — should finish quickly and tell you if the pipeline works.
    "tiny": {
        "N": [256],
        "L": [200.0],
        "sigma": [2.0],
        "gamma_drag": [0.03],
        "ramp_time": [2.0],
    },
    # Default: small but informative (shows effect of drag and resolution).
    "small": {
        "N": [256, 512],
        "L": [200.0],
        "sigma": [2.0],
        "gamma_drag": [0.0, 0.03],
        "ramp_time": [2.0],
    },
    # Medium: closer to the grid suggested in poisson_01.md but still finite.
    "medium": {
        "N": [256, 384, 512],
        "L": [150.0, 200.0, 300.0],
        "sigma": [1.0, 2.0, 3.0],
        "gamma_drag": [0.0, 0.01, 0.03],
        "ramp_time": [0.0, 2.0, 5.0],
    },
    # Full: matches the full grid suggested in poisson_01.md (can be huge).
    "full": {
        "N": [256, 384, 512, 768],
        "L": [150.0, 200.0, 300.0],
        "sigma": [1.0, 2.0, 3.0, 4.0],
        "gamma_drag": [0.0, 0.01, 0.03, 0.1],
        "ramp_time": [0.0, 1.0, 2.0, 5.0],
    },
}


# ----------------------------
# Pass criteria (from poisson_01.md suggestion)
# ----------------------------

@dataclass(frozen=True)
class PassCriteria:
    slope_tol: float = 0.05
    amp_tol: float = 0.05
    frac_long_min: float = 0.999
    lag_max: float = 0.05


DEFAULT_CRIT = PassCriteria()


def evaluate_pass(diag: Dict[str, Any], *, mdot_target: float, crit: PassCriteria) -> Dict[str, Any]:
    """Return per-metric pass booleans + derived ratios."""
    fits = diag.get("fits", {}) if isinstance(diag.get("fits"), dict) else {}
    fftm = diag.get("fft_metrics", {}) if isinstance(diag.get("fft_metrics"), dict) else {}
    mono = diag.get("monopole_fit", {}) if isinstance(diag.get("monopole_fit"), dict) else {}

    dpsi_clean_slope = safe_float(fits.get("dpsi_clean_slope"))
    g_slope = safe_float(fits.get("g_slope"))
    lag = safe_float(fftm.get("lag_rms_over_psiP_rms"))
    frac_long = safe_float(fftm.get("frac_long"))

    A_fit = safe_float(mono.get("A_fit"))
    A_expected = safe_float(mono.get("A_expected"))
    mdot_from_flux = safe_float(mono.get("Mdot_from_flux"))

    A_rel = None
    if A_fit is not None and A_expected not in (None, 0.0):
        A_rel = A_fit / A_expected

    mdot_rel = None
    if mdot_from_flux is not None and mdot_target not in (0.0, None):
        mdot_rel = mdot_from_flux / mdot_target

    # Evaluate
    pass_dpsi = None
    if dpsi_clean_slope is not None:
        pass_dpsi = abs(dpsi_clean_slope + 1.0) <= crit.slope_tol

    pass_g = None
    if g_slope is not None:
        pass_g = abs(g_slope + 2.0) <= crit.slope_tol

    pass_A = None
    if A_rel is not None:
        pass_A = abs(A_rel - 1.0) <= crit.amp_tol

    pass_mdot = None
    if mdot_rel is not None:
        pass_mdot = abs(mdot_rel - 1.0) <= crit.amp_tol

    pass_frac = None
    if frac_long is not None:
        pass_frac = frac_long >= crit.frac_long_min

    pass_lag = None
    if lag is not None:
        pass_lag = lag <= crit.lag_max

    # All-pass if every metric that exists passes (missing values don't fail)
    flags = [pass_dpsi, pass_g, pass_A, pass_mdot, pass_frac, pass_lag]
    known = [f for f in flags if f is not None]
    all_pass = bool(known) and all(known)

    return {
        "derived": {
            "dpsi_clean_slope": dpsi_clean_slope,
            "g_slope": g_slope,
            "lag_rms_over_psiP_rms": lag,
            "frac_long": frac_long,
            "A_fit": A_fit,
            "A_expected": A_expected,
            "A_rel": A_rel,
            "Mdot_from_flux": mdot_from_flux,
            "Mdot_rel": mdot_rel,
        },
        "pass": {
            "dpsi_clean_slope": pass_dpsi,
            "g_slope": pass_g,
            "A_rel": pass_A,
            "Mdot_rel": pass_mdot,
            "frac_long": pass_frac,
            "lag": pass_lag,
            "all": all_pass,
        },
    }


# ----------------------------
# Runner
# ----------------------------

def iter_grid(grid: Dict[str, Sequence[Any]]) -> Iterable[Dict[str, Any]]:
    keys = list(grid.keys())
    for vals in itertools.product(*(grid[k] for k in keys)):
        yield {k: vals[i] for i, k in enumerate(keys)}


def run_stepE_one(
    *,
    stepE_path: Path,
    params: Dict[str, Any],
    fixed_args: Dict[str, Any],
    crit: PassCriteria,
    out_jsonl: Path,
    log_dir: Optional[Path],
    skip_if_done: bool,
    done_ids: Set[str],
    stream_child_json: bool,
    stop_on_pass: bool,
    min_step_for_stop: int,
) -> Tuple[bool, Dict[str, Any]]:
    """Run one configuration. Returns (completed_ok, summary_row)."""

    # Split: physics-like params used in run_id + everything else
    run_params = {**params, **fixed_args}

    # Make a run_id from the physics knobs + a few numerics that change the outcome
    rid_inputs = {
        # physics-ish
        "N": run_params["N"],
        "L": run_params["L"],
        "sigma": run_params["sigma"],
        "gamma_drag": run_params["gamma_drag"],
        "ramp_time": run_params["ramp_time"],
        "Mdot": run_params["Mdot"],
        "rho0": run_params["rho0"],
        "cs": run_params["cs"],
        "refill_momentum": run_params["refill_momentum"],
        # numerics that change outcome
        "dtype": run_params["dtype"],
        "cfl": run_params["cfl"],
        "dt": run_params["dt"],
        "dt_max": run_params["dt_max"],
        "rho_floor": run_params["rho_floor"],
        "steps": run_params["steps"],
        "diag_every": run_params["diag_every"],
        "fit_rmin": run_params["fit_rmin"],
        "fit_rmax": run_params["fit_rmax"],
        "nbins": run_params["nbins"],
        "remove_mean_momentum": run_params["remove_mean_momentum"],
    }
    run_id = make_run_id(rid_inputs)

    if skip_if_done and run_id in done_ids:
        row = {
            "event": "sweep_skip",
            "run_id": run_id,
            "params": rid_inputs,
            "reason": "already_in_output_jsonl",
            "ts": time.time(),
        }
        print(json.dumps(row), flush=True)
        return True, row

    # Build command
    cmd = [
        sys.executable,
        "-u",
        str(stepE_path),
        "--backend", str(run_params["backend"]),
        "--dtype", str(run_params["dtype"]),
        "--N", str(int(run_params["N"])),
        "--L", str(float(run_params["L"])),
        "--rho0", str(float(run_params["rho0"])),
        "--cs", str(float(run_params["cs"])),
        "--sigma", str(float(run_params["sigma"])),
        "--Mdot", str(float(run_params["Mdot"])),
        "--ramp_time", str(float(run_params["ramp_time"])),
        "--refill_momentum", str(run_params["refill_momentum"]),
        "--gamma_drag", str(float(run_params["gamma_drag"])),
        "--rho_floor", str(float(run_params["rho_floor"])),
        "--steps", str(int(run_params["steps"])),
        "--diag_every", str(int(run_params["diag_every"])),
        "--cfl", str(float(run_params["cfl"])),
        "--nbins", str(int(run_params["nbins"])),
        "--fit_rmin", str(float(run_params["fit_rmin"])),
        "--fit_rmax", str(float(run_params["fit_rmax"])),
        "--sample_points", str(int(run_params["sample_points"])),
    ]
    if run_params["dt"] is not None:
        cmd += ["--dt", str(float(run_params["dt"]))]
    if run_params["dt_max"] is not None:
        cmd += ["--dt_max", str(float(run_params["dt_max"]))]
    if run_params["remove_mean_momentum"]:
        cmd += ["--remove_mean_momentum"]

    # Logging
    log_fp = None
    if log_dir is not None:
        log_dir.mkdir(parents=True, exist_ok=True)
        log_fp = (log_dir / f"{run_id}.log").open("w", encoding="utf-8")

    start_ts = time.time()
    start = {
        "event": "run_start",
        "run_id": run_id,
        "stepE": str(stepE_path),
        "cmd": cmd,
        "params": rid_inputs,
        "ts": start_ts,
    }
    print(json.dumps(start), flush=True)

    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env,
    )

    last_diag: Optional[Dict[str, Any]] = None
    last_diag_pass: Optional[Dict[str, Any]] = None
    last_diag_step: Optional[int] = None

    try:
        assert proc.stdout is not None
        for raw_line in proc.stdout:
            line = raw_line.rstrip("\n")
            if log_fp is not None:
                log_fp.write(line + "\n")
                log_fp.flush()

            # Try parse JSON (StepE prints JSON lines)
            obj = None
            try:
                obj = json.loads(line)
            except Exception:
                obj = None

            if stream_child_json:
                # Echo the child's JSON (or raw) line
                print(line, flush=True)

            if isinstance(obj, dict) and obj.get("event") == "diag":
                last_diag = obj
                try:
                    last_diag_step = int(obj.get("step"))
                except Exception:
                    last_diag_step = None

                # Evaluate pass metrics for live display
                evald = evaluate_pass(obj, mdot_target=float(run_params["Mdot"]), crit=crit)
                last_diag_pass = evald

                # Compact live line
                d = evald["derived"]
                p = evald["pass"]
                msg = (
                    f"[{run_id}] step={obj.get('step')} "
                    f"lag={fmt(d.get('lag_rms_over_psiP_rms'), prec=3)} "
                    f"fracL={fmt(d.get('frac_long'), prec=4)} "
                    f"dpsi={fmt(d.get('dpsi_clean_slope'), prec=4)} "
                    f"g={fmt(d.get('g_slope'), prec=4)} "
                    f"Arel={fmt(d.get('A_rel'), prec=4)} "
                    f"Mrel={fmt(d.get('Mdot_rel'), prec=4)} "
                    f"PASS={p.get('all')}"
                )
                # Always print compact progress even if not streaming child JSON
                if not stream_child_json:
                    print(msg, flush=True)

                # Optional early stop
                if stop_on_pass and p.get("all") and (last_diag_step is not None) and (last_diag_step >= min_step_for_stop):
                    # terminate gently; we already have the diag we want
                    proc.terminate()
                    break

        rc = proc.wait()
    except KeyboardInterrupt:
        proc.terminate()
        rc = proc.wait()
        raise
    finally:
        if log_fp is not None:
            log_fp.close()

    end_ts = time.time()

    summary: Dict[str, Any] = {
        "event": "run_done",
        "run_id": run_id,
        "params": rid_inputs,
        "returncode": int(rc),
        "wall_s": float(end_ts - start_ts),
        "ts": end_ts,
        "last_diag": last_diag,
        "last_diag_eval": last_diag_pass,
    }

    # Write one JSONL row per run (append-only, so we can resume)
    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with out_jsonl.open("a", encoding="utf-8") as f:
        f.write(json.dumps(summary) + "\n")

    print(json.dumps(summary), flush=True)

    ok = (rc == 0) and (last_diag is not None)
    return ok, summary


def main() -> None:
    ap = argparse.ArgumentParser()

    # Which StepE to run
    ap.add_argument("--stepE", type=str, default=None, help="Path to stepE_euler_throat_longitudinal.py (defaults to sibling file).")

    # Sweep definition
    ap.add_argument("--preset", choices=sorted(PRESETS.keys()), default="small", help="Which preset grid to run.")
    ap.add_argument("--N_list", type=str, default=None, help="Comma list overriding preset N, e.g. '256,512'.")
    ap.add_argument("--L_list", type=str, default=None, help="Comma list overriding preset L, e.g. '150,200,300'.")
    ap.add_argument("--sigma_list", type=str, default=None, help="Comma list overriding preset sigma, e.g. '1,2,3,4'.")
    ap.add_argument("--gamma_drag_list", type=str, default=None, help="Comma list overriding preset gamma_drag, e.g. '0,0.03,0.1'.")
    ap.add_argument("--ramp_time_list", type=str, default=None, help="Comma list overriding preset ramp_time, e.g. '0,2,5'.")

    # Fixed physics params
    ap.add_argument("--Mdot", type=float, default=20.0)
    ap.add_argument("--rho0", type=float, default=1.0)
    ap.add_argument("--cs", type=float, default=7.0710678118654755)
    ap.add_argument("--refill_momentum", choices=["comoving", "at_rest"], default="at_rest")

    # Numerics / backend
    ap.add_argument("--backend", choices=["auto", "numpy", "cupy"], default="auto")
    ap.add_argument("--dtype", choices=["float32", "float64"], default="float32")
    ap.add_argument("--cfl", type=float, default=0.25)
    ap.add_argument("--dt", type=float, default=None)
    ap.add_argument("--dt_max", type=float, default=None)
    ap.add_argument("--rho_floor", type=float, default=1e-6)

    ap.add_argument("--steps", type=int, default=10000)
    ap.add_argument("--diag_every", type=int, default=400)

    # Diagnostics parameters (keep consistent across sweep to compare)
    ap.add_argument("--nbins", type=int, default=180)
    ap.add_argument("--fit_rmin", type=float, default=8.0)
    ap.add_argument("--fit_rmax", type=float, default=50.0)
    ap.add_argument("--sample_points", type=int, default=0, help="Forwarded to StepE. Keep at 0 for compact logs.")

    ap.add_argument("--remove_mean_momentum", action="store_true")

    # Output
    ap.add_argument("--out", type=str, default="stepE_sweep.jsonl", help="JSONL file to append one summary row per run.")
    ap.add_argument("--log_dir", type=str, default="stepE_sweep_logs", help="Directory to store per-run raw stdout logs.")
    ap.add_argument("--no_logs", action="store_true", help="Disable per-run log files.")
    ap.add_argument("--rerun", action="store_true", help="Do not skip already-completed run_ids in --out.")

    # Printing & early stop
    ap.add_argument("--stream_child_json", action="store_true", help="Print the full JSON lines from StepE (very verbose).")
    ap.add_argument("--stop_on_pass", action="store_true", help="Terminate a run early once it passes all criteria at a diag.")
    ap.add_argument("--min_step_for_stop", type=int, default=3000)

    # Pass criteria overrides
    ap.add_argument("--slope_tol", type=float, default=DEFAULT_CRIT.slope_tol)
    ap.add_argument("--amp_tol", type=float, default=DEFAULT_CRIT.amp_tol)
    ap.add_argument("--frac_long_min", type=float, default=DEFAULT_CRIT.frac_long_min)
    ap.add_argument("--lag_max", type=float, default=DEFAULT_CRIT.lag_max)

    # Optional limit for debugging
    ap.add_argument("--max_runs", type=int, default=None, help="If set, run at most this many configurations.")

    args = ap.parse_args()

    crit = PassCriteria(
        slope_tol=float(args.slope_tol),
        amp_tol=float(args.amp_tol),
        frac_long_min=float(args.frac_long_min),
        lag_max=float(args.lag_max),
    )

    here = Path(__file__).resolve().parent
    stepE_path = Path(args.stepE).resolve() if args.stepE is not None else (here / "stepE_euler_throat_longitudinal.py")
    if not stepE_path.exists():
        raise FileNotFoundError(f"Could not find StepE script at: {stepE_path}")

    out_path = Path(args.out).resolve() if Path(args.out).is_absolute() else (here / args.out)
    log_dir = None if args.no_logs else (Path(args.log_dir).resolve() if Path(args.log_dir).is_absolute() else (here / args.log_dir))

    # Build grid (preset overridden by explicit lists)
    preset = PRESETS[args.preset]
    grid = {
        "N": _parse_csv_list(args.N_list, cast=int) or preset["N"],
        "L": _parse_csv_list(args.L_list, cast=float) or preset["L"],
        "sigma": _parse_csv_list(args.sigma_list, cast=float) or preset["sigma"],
        "gamma_drag": _parse_csv_list(args.gamma_drag_list, cast=float) or preset["gamma_drag"],
        "ramp_time": _parse_csv_list(args.ramp_time_list, cast=float) or preset["ramp_time"],
    }

    fixed_args = {
        "backend": args.backend,
        "dtype": args.dtype,
        "Mdot": float(args.Mdot),
        "rho0": float(args.rho0),
        "cs": float(args.cs),
        "refill_momentum": args.refill_momentum,
        "cfl": float(args.cfl),
        "dt": args.dt if args.dt is None else float(args.dt),
        "dt_max": args.dt_max if args.dt_max is None else float(args.dt_max),
        "rho_floor": float(args.rho_floor),
        "steps": int(args.steps),
        "diag_every": int(args.diag_every),
        "nbins": int(args.nbins),
        "fit_rmin": float(args.fit_rmin),
        "fit_rmax": float(args.fit_rmax),
        "sample_points": int(args.sample_points),
        "remove_mean_momentum": bool(args.remove_mean_momentum),
    }

    # Resume / skip logic
    done_ids = load_existing_run_ids(out_path)
    skip_if_done = not args.rerun

    sweep_init = {
        "event": "sweep_init",
        "ts": time.time(),
        "stepE": str(stepE_path),
        "out": str(out_path),
        "log_dir": None if log_dir is None else str(log_dir),
        "preset": args.preset,
        "grid": grid,
        "fixed_args": fixed_args,
        "pass_criteria": crit.__dict__,
        "already_done_run_ids": len(done_ids),
        "skip_if_done": skip_if_done,
        "notes": "This sweep runs StepE once per grid point and appends one summary JSON row per run.",
    }
    print(json.dumps(sweep_init), flush=True)

    # Enumerate configs
    configs = list(iter_grid(grid))
    if args.max_runs is not None:
        configs = configs[: max(0, int(args.max_runs))]

    total = len(configs)
    ok_count = 0
    for i, cfg in enumerate(configs, start=1):
        banner = {"event": "sweep_progress", "i": i, "total": total, "cfg": cfg, "ts": time.time()}
        print(json.dumps(banner), flush=True)

        ok, _row = run_stepE_one(
            stepE_path=stepE_path,
            params=cfg,
            fixed_args=fixed_args,
            crit=crit,
            out_jsonl=out_path,
            log_dir=log_dir,
            skip_if_done=skip_if_done,
            done_ids=done_ids,
            stream_child_json=bool(args.stream_child_json),
            stop_on_pass=bool(args.stop_on_pass),
            min_step_for_stop=int(args.min_step_for_stop),
        )
        if ok:
            ok_count += 1

    sweep_done = {
        "event": "sweep_done",
        "ts": time.time(),
        "total": total,
        "completed_ok": ok_count,
        "out": str(out_path),
        "log_dir": None if log_dir is None else str(log_dir),
    }
    print(json.dumps(sweep_done), flush=True)


if __name__ == "__main__":
    main()
