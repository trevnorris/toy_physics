#!/usr/bin/env python3
"""
Analyze slabScalingOutputs produced by the Mathematica WL sweep.

Expected files:
  - slabScalingOutputs/sweep_summary.csv (optional but preferred)
  - slabScalingOutputs/*-rows.csv (per-case raw points)

It will:
  1) Parse parameters from filenames like:
       aQ0p5-k1p-EOSIncompressible-DP1pem7-FM3D-Hg1pe6-rows.csv
  2) Load & concatenate all rows into one DataFrame
  3) Compute per-case summaries (ratio stats, simple slopes) from rows
  4) Optionally merge with sweep_summary.csv if present
  5) Export combined datasets + a few plots for regime transitions

Usage:
  python analyze_slab_outputs.py --dir slabScalingOutputs --out analysis_out --plots
"""

from __future__ import annotations

import argparse
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# Filename parsing utilities
# -----------------------------

def parse_mathematica_token_number(tok: str) -> float:
    """
    Parse tokens like:
      0p5   -> 0.5
      1p25  -> 1.25
      1pe6  -> 1e6
      1pem7 -> 1e-7
    Also accepts normal floats/ints.

    Returns float.
    """
    tok = tok.strip()
    # Common cleanups
    tok = tok.replace('"', '').replace("'", "")
    # Direct numeric?
    try:
        return float(tok)
    except ValueError:
        pass

    # Match like 1pem7 or 1pe6 or 0p5
    # Interpret:
    #   p -> decimal point
    #   pe -> *10^
    #   pem -> *10^-   (m meaning minus)
    # Examples:
    #   1pem7 -> 1e-7
    #   2pe10 -> 2e10
    #   0p5 -> 0.5
    s = tok
    # scientific exponent
    if "pem" in s:
        base, exp = s.split("pem", 1)
        base = base.replace("p", ".")
        return float(base) * (10.0 ** (-int(exp)))
    if "pe" in s:
        base, exp = s.split("pe", 1)
        base = base.replace("p", ".")
        return float(base) * (10.0 ** (int(exp)))

    # decimal only
    if "p" in s:
        return float(s.replace("p", "."))
    raise ValueError(f"Unrecognized numeric token: {tok!r}")


@dataclass(frozen=True)
class CaseParams:
    alphaQ: float
    kappa: float
    EOS: str
    DeltaPcrit: float
    FlowModel: str
    Hgeom: float

    def as_dict(self) -> Dict[str, object]:
        return {
            "alphaQ": self.alphaQ,
            "kappa": self.kappa,
            "EOS": self.EOS,
            "DeltaPcrit": self.DeltaPcrit,
            "FlowModel": self.FlowModel,
            "Hgeom": self.Hgeom,
        }


CASE_RE = re.compile(
    r"""
    aQ(?P<aQ>[^-]+)
    -k(?P<k>[^-]+)
    -EOS(?P<EOS>[^-]+)
    -DP(?P<DP>[^-]+)
    -FM(?P<FM>[^-]+)
    -Hg(?P<Hg>[^-]+)
    """,
    re.VERBOSE
)

def parse_case_params_from_filename(name: str) -> Optional[CaseParams]:
    """
    Extract CaseParams from a filename stem containing the case id.
    Returns None if it doesn't match.
    """
    stem = Path(name).stem
    # remove trailing "-rows" if present
    stem = re.sub(r"-rows$", "", stem)

    m = CASE_RE.search(stem)
    if not m:
        return None

    aQ = parse_mathematica_token_number(m.group("aQ"))
    kappa = parse_mathematica_token_number(m.group("k"))
    eos = m.group("EOS")
    dp = parse_mathematica_token_number(m.group("DP"))
    fm = m.group("FM")
    hg = parse_mathematica_token_number(m.group("Hg"))

    return CaseParams(alphaQ=aQ, kappa=kappa, EOS=eos, DeltaPcrit=dp, FlowModel=fm, Hgeom=hg)


# -----------------------------
# Math helpers for summaries
# -----------------------------

def safe_log10(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.full_like(x, np.nan, dtype=float)
    mask = x > 0
    out[mask] = np.log10(x[mask])
    return out

def linfit_slope_intercept(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
    """
    Return (slope, intercept, r2) for y ~ intercept + slope*x.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]; y = y[mask]
    if len(x) < 2:
        return (np.nan, np.nan, np.nan)
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    yhat = slope*x + intercept
    ss_res = np.sum((y - yhat)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
    return (float(slope), float(intercept), float(r2))


# -----------------------------
# Loading
# -----------------------------

def load_all_rows(outdir: Path) -> pd.DataFrame:
    files = sorted(outdir.glob("*-rows.csv"))
    if not files:
        raise FileNotFoundError(f"No *-rows.csv found in {outdir}")

    all_parts = []
    skipped = 0

    for f in files:
        params = parse_case_params_from_filename(f.name)
        if params is None:
            skipped += 1
            continue

        df = pd.read_csv(f)

        # normalize columns (case-insensitive)
        rename_map = {}
        for col in df.columns:
            c = col.strip().lower()
            if c == "m":
                rename_map[col] = "M"
            elif c == "q":
                rename_map[col] = "Q"
            elif c in ("h_eff", "heff"):
                rename_map[col] = "Heff"
            elif c in ("v_inf", "vinf", "v∞"):
                rename_map[col] = "vinf"
            elif c in ("heffoverhgeom", "heff_over_hgeom", "heff/hgeom", "ratio"):
                rename_map[col] = "HeffOverHgeom"
            elif c == "regime":
                rename_map[col] = "Regime"

        df = df.rename(columns=rename_map)

        # Required numeric columns
        required = {"M", "Q", "Heff", "vinf", "HeffOverHgeom"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{f.name} missing columns: {missing}. Found: {list(df.columns)}")

        # If Regime missing, compute it (default thresholds; tweak if you want)
        if "Regime" not in df.columns:
            rlow, rhigh = 0.3, 3.0
            ratio = pd.to_numeric(df["HeffOverHgeom"], errors="coerce")
            df["Regime"] = np.select(
                [
                    ratio.isna(),
                    ratio < rlow,
                    ratio > rhigh
                ],
                [
                    "no-root",
                    "3D-dominated",
                    "2D-dominated"
                ],
                default="transition"
            )

        # attach params
        for k, v in params.as_dict().items():
            df[k] = v
        df["caseFile"] = f.name

        all_parts.append(df)

    out = pd.concat(all_parts, ignore_index=True)
    if skipped:
        print(f"[warn] skipped {skipped} *-rows.csv files that didn't match the case-id pattern")
    return out


def load_sweep_summary_if_present(outdir: Path) -> Optional[pd.DataFrame]:
    f = outdir / "sweep_summary.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f)
    return df


# -----------------------------
# Per-case summaries from rows
# -----------------------------

def summarize_cases_from_rows(rows: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["alphaQ", "kappa", "EOS", "DeltaPcrit", "FlowModel", "Hgeom"]

    def one_case(g: pd.DataFrame) -> pd.Series:
        g = g.copy()
        g = g.sort_values("M")

        # logs for fits
        logM = safe_log10(g["M"].to_numpy())
        logH = safe_log10(g["Heff"].to_numpy())
        logv = safe_log10(g["vinf"].to_numpy())

        bH, aH, r2H = linfit_slope_intercept(logM, logH)
        bTF, aTF, r2TF = linfit_slope_intercept(logv, logM)

        ratio = pd.to_numeric(g["HeffOverHgeom"], errors="coerce").to_numpy(dtype=float)
        M = pd.to_numeric(g["M"], errors="coerce").to_numpy(dtype=float)

        good = np.isfinite(ratio) & np.isfinite(M) & (M > 0)
        ratio = ratio[good]
        M = M[good]

        rlow, rhigh = 0.3, 3.0
        mask_low = (ratio < rlow)
        mask_high = (ratio > rhigh)

        # low-regime fit: logM ~ a + b logv
        bTF_low, aTF_low, r2TF_low = linfit_slope_intercept(
            safe_log10(g["vinf"].to_numpy()[good][mask_low]),
            safe_log10(g["M"].to_numpy()[good][mask_low])
        )

        bTF_high, aTF_high, r2TF_high = linfit_slope_intercept(
            safe_log10(g["vinf"].to_numpy()[good][mask_high]),
            safe_log10(g["M"].to_numpy()[good][mask_high])
        )

        ratio_min = np.nanmin(ratio) if ratio.size else np.nan
        ratio_med = np.nanmedian(ratio) if ratio.size else np.nan
        ratio_max = np.nanmax(ratio) if ratio.size else np.nan

        # closest point to ratio=1 (always useful)
        if ratio.size:
            idx = int(np.nanargmin(np.abs(np.log10(ratio) - 0.0)))
            Mclosest = float(M[idx])
            rclosest = float(ratio[idx])
        else:
            Mclosest = np.nan
            rclosest = np.nan

        # detect transition mass where ratio crosses 1 in log space
        crosses1 = bool(np.isfinite(ratio_min) and np.isfinite(ratio_max) and (ratio_min < 1.0 < ratio_max))
        Mtrans = np.nan
        status = "NoData"

        if ratio.size == 0:
            status = "NoData"
        elif crosses1:
            status = "CrossesInRange"
            order = np.argsort(M)
            M2 = M[order]
            lr = safe_log10(ratio[order])
            lm = safe_log10(M2)
            for i in range(len(M2) - 1):
                y1, y2 = lr[i], lr[i + 1]
                if np.isfinite(y1) and np.isfinite(y2) and (y1 <= 0 <= y2):
                    x1, x2 = lm[i], lm[i + 1]
                    if y2 != y1 and np.isfinite(x1) and np.isfinite(x2):
                        xt = x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
                        Mtrans = float(10 ** xt)
                    break
        else:
            # no crossing in sampled range: indicate where it lies
            if ratio_max < 1.0:
                status = "AboveRange"   # would cross at higher M than sampled (if ever)
            elif ratio_min > 1.0:
                status = "BelowRange"   # would cross at lower M than sampled (if ever)
            else:
                status = "NoCross"

        return pd.Series({
            "nPoints": int(len(g)),
            "ratioMin": ratio_min,
            "ratioMed": ratio_med,
            "ratioMax": ratio_max,
            "crosses1": crosses1,
            "Mtrans": Mtrans,
            "MtransStatus": status,
            "MclosestTo1": Mclosest,
            "ratioClosestTo1": rclosest,
            "bH_all": bH,
            "aH_all": aH,
            "R2_H_all": r2H,
            "bTF_all": bTF,
            "aTF_all": aTF,
            "R2_TF_all": r2TF,
            "bTF_low": bTF_low, "R2_TF_low": r2TF_low,
            "bTF_high": bTF_high, "R2_TF_high": r2TF_high,
            "nLow": int(np.sum(mask_low)), "nHigh": int(np.sum(mask_high)),
        })

    summary = rows.groupby(group_cols, dropna=False).apply(one_case).reset_index()
    return summary


# -----------------------------
# Plotting
# -----------------------------

def plot_transition_examples(rows: pd.DataFrame, outdir: Path, max_plots: int = 12) -> None:
    """
    Make a small set of example plots that clearly show:
      - ratio vs M (where it crosses 1)
      - vinf vs M
    Picks cases where ratioMax is large (i.e. likely to cross).
    """
    summary = summarize_cases_from_rows(rows)
    summary = summary.sort_values("ratioMax", ascending=False).head(max_plots)

    for _, s in summary.iterrows():
        mask = (
            (rows["alphaQ"] == s["alphaQ"]) &
            (rows["kappa"] == s["kappa"]) &
            (rows["EOS"] == s["EOS"]) &
            (rows["DeltaPcrit"] == s["DeltaPcrit"]) &
            (rows["FlowModel"] == s["FlowModel"]) &
            (rows["Hgeom"] == s["Hgeom"])
        )
        g = rows.loc[mask].copy().sort_values("M")

        fig = plt.figure(figsize=(10, 4))

        ax1 = fig.add_subplot(1, 2, 1)
        ax1.loglog(g["M"], g["HeffOverHgeom"])
        ax1.axhline(1.0, linestyle="--")
        ax1.set_xlabel("M")
        ax1.set_ylabel("Heff/Hgeom")
        ax1.set_title("Crossover ratio")

        ax2 = fig.add_subplot(1, 2, 2)
        ax2.loglog(g["M"], g["vinf"])
        ax2.set_xlabel("M")
        ax2.set_ylabel("vinf")
        ax2.set_title("Plateau speed")

        title = (
            f"aQ={s['alphaQ']}, k={s['kappa']}, EOS={s['EOS']}, DP={s['DeltaPcrit']}, "
            f"FM={s['FlowModel']}, Hg={s['Hgeom']}"
        )
        fig.suptitle(title, fontsize=10)
        fig.tight_layout()

        fname = (
            f"example-aQ{s['alphaQ']}-k{s['kappa']}-EOS{s['EOS']}-DP{s['DeltaPcrit']}"
            f"-FM{s['FlowModel']}-Hg{s['Hgeom']}.png"
        ).replace("/", "_")
        fig.savefig(outdir / fname, dpi=200)
        plt.close(fig)

def plot_only_crossing_cases(rows: pd.DataFrame, outdir: Path, alphaQ: float = 1.0) -> None:
    summ = summarize_cases_from_rows(rows)
    summ = summ[(summ["crosses1"] == True) & (summ["alphaQ"] == alphaQ)].copy()
    if summ.empty:
        print(f"[info] no crossing cases found for alphaQ={alphaQ}")
        return

    # Limit to a manageable set (most interesting = widest ratio span)
    summ["ratioSpan"] = summ["ratioMax"] / summ["ratioMin"]
    summ = summ.sort_values("ratioSpan", ascending=False).head(20)

    fig1, ax1 = plt.subplots(figsize=(9, 6))
    fig2, ax2 = plt.subplots(figsize=(9, 6))

    for _, s in summ.iterrows():
        mask = (
            (rows["alphaQ"] == s["alphaQ"]) &
            (rows["kappa"] == s["kappa"]) &
            (rows["EOS"] == s["EOS"]) &
            (rows["DeltaPcrit"] == s["DeltaPcrit"]) &
            (rows["FlowModel"] == s["FlowModel"]) &
            (rows["Hgeom"] == s["Hgeom"])
        )
        g = rows.loc[mask].copy().sort_values("M")

        label = f"k={s['kappa']}, EOS={s['EOS']}, DP={s['DeltaPcrit']:.1e}, FM={s['FlowModel']}, Hg={s['Hgeom']:.1e}"

        ax1.loglog(g["M"], g["HeffOverHgeom"], marker="o", linewidth=1, label=label)
        ax2.loglog(g["M"], g["vinf"], marker="o", linewidth=1, label=label)

        # mark transition mass
        if np.isfinite(s["Mtrans"]):
            ax1.axvline(s["Mtrans"], linestyle="--", linewidth=0.8)
            ax2.axvline(s["Mtrans"], linestyle="--", linewidth=0.8)

    ax1.axhline(1.0, linestyle="--", linewidth=1)
    ax1.set_xlabel("M")
    ax1.set_ylabel("Heff/Hgeom")
    ax1.set_title(f"Crossover cases (alphaQ={alphaQ}): ratio crosses 1")
    ax1.legend(fontsize=7, ncol=1, loc="best")
    fig1.tight_layout()
    fig1.savefig(outdir / f"crossing_ratio_alphaQ_{alphaQ}.png", dpi=200)
    plt.close(fig1)

    ax2.set_xlabel("M")
    ax2.set_ylabel("vinf")
    ax2.set_title(f"Crossover cases (alphaQ={alphaQ}): vinf vs M with transition markers")
    ax2.legend(fontsize=7, ncol=1, loc="best")
    fig2.tight_layout()
    fig2.savefig(outdir / f"crossing_vinf_alphaQ_{alphaQ}.png", dpi=200)
    plt.close(fig2)



def plot_heatmap_like_scatter(summary: pd.DataFrame, outdir: Path) -> None:
    """
    A simple “heatmap-ish” view without seaborn:
    - scatter of Hgeom vs DeltaPcrit, colored by bTF_all, one panel per alphaQ
    (matplotlib default colormap)
    """
    alphas = sorted(summary["alphaQ"].unique())
    for aQ in alphas:
        s = summary[summary["alphaQ"] == aQ].copy()
        if s.empty:
            continue

        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(1, 1, 1)

        x = s["Hgeom"].to_numpy(dtype=float)
        y = s["DeltaPcrit"].to_numpy(dtype=float)
        c = s["bTF_all"].to_numpy(dtype=float)

        sc = ax.scatter(x, y, c=c)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Hgeom")
        ax.set_ylabel("DeltaPcrit")
        ax.set_title(f"bTF_all across levers (alphaQ={aQ})")
        fig.colorbar(sc, ax=ax, label="bTF_all (logM vs logv slope)")
        fig.tight_layout()

        fig.savefig(outdir / f"scatter_bTF_alphaQ_{aQ}.png", dpi=200)
        plt.close(fig)


# -----------------------------
# Main
# -----------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="slabScalingOutputs", help="directory with CSV outputs")
    ap.add_argument("--out", default="analysis_out", help="output directory for merged data + plots")
    ap.add_argument("--plots", action="store_true", help="generate plots")
    ap.add_argument("--alphaQ_focus", type=float, default=1.0, help="alphaQ value to focus crossover plots on")
    ap.add_argument("--max_cross_plots", type=int, default=20, help="max number of crossing cases to plot")
    args = ap.parse_args()

    in_dir = Path(args.dir).expanduser().resolve()
    out_dir = Path(args.out).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[info] reading rows from: {in_dir}")
    rows = load_all_rows(in_dir)

    # Save combined rows
    rows_csv = out_dir / "all_rows.csv"
    rows.to_csv(rows_csv, index=False)
    print(f"[ok] wrote {rows_csv}")

    # Parquet is much faster/smaller if you have pyarrow installed
    try:
        rows_parq = out_dir / "all_rows.parquet"
        rows.to_parquet(rows_parq, index=False)
        print(f"[ok] wrote {rows_parq}")
    except Exception as e:
        print(f"[warn] parquet export skipped: {e}")

    # Summaries from rows
    case_summary = summarize_cases_from_rows(rows)
    summ_csv = out_dir / "case_summary_from_rows.csv"
    case_summary.to_csv(summ_csv, index=False)
    print(f"[ok] wrote {summ_csv}")

    # Merge with sweep_summary.csv if present
    sweep = load_sweep_summary_if_present(in_dir)
    if sweep is not None:
        key_cols = ["alphaQ", "kappa", "EOS", "DeltaPcrit", "FlowModel", "Hgeom"]
        if all(c in sweep.columns for c in key_cols):
            merged = case_summary.merge(sweep, on=key_cols, how="left", suffixes=("_from_rows", "_from_wl"))
            merged_csv = out_dir / "merged_case_summary.csv"
            merged.to_csv(merged_csv, index=False)
            print(f"[ok] wrote {merged_csv}")
        else:
            print("[warn] sweep_summary.csv present but missing expected key columns; skipping merge")

    # Better console summary: focus on TRUE crossover cases (crosses ratio=1) for alphaQ focus & FlowModel=Crossover
    focus = args.alphaQ_focus
    cross = case_summary[
        (case_summary["alphaQ"] == focus) &
        (case_summary["FlowModel"].astype(str).str.lower() == "crossover") &
        (case_summary["crosses1"] == True)
    ].copy()

    if not cross.empty:
        cross["ratioSpan"] = cross["ratioMax"] / cross["ratioMin"]
        cross["closeness"] = np.abs(np.log10(cross["ratioClosestTo1"]))  # smaller is better
        cross = cross.sort_values(["closeness", "ratioSpan"], ascending=[True, False])

        print(f"\nTop crossover cases (alphaQ={focus}, FlowModel=Crossover, crosses in-range):")
        print(cross.head(12)[[
            "alphaQ", "kappa", "EOS", "DeltaPcrit", "FlowModel", "Hgeom",
            "ratioMin", "ratioMed", "ratioMax",
            "Mtrans", "MtransStatus",
            "bTF_all", "R2_TF_all"
        ]].to_string(index=False))
    else:
        print(f"\n[info] No in-range crossover cases found for alphaQ={focus} with FlowModel=Crossover.")

    # Still useful: show extremes, but label them correctly
    print("\nTop 10 cases by ratioMax (these are often 'deep 2D', not necessarily crossings):")
    print(case_summary.sort_values("ratioMax", ascending=False).head(10)[
        ["alphaQ", "kappa", "EOS", "DeltaPcrit", "FlowModel", "Hgeom",
         "ratioMin", "ratioMed", "ratioMax", "crosses1", "MtransStatus", "bTF_all", "R2_TF_all"]
    ].to_string(index=False))

    # Plots
    if args.plots:
        plots_dir = out_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        # Primary paper-facing plot: only true crossing cases, alphaQ focus
        plot_only_crossing_cases(rows, plots_dir, alphaQ=focus)

        # Secondary diagnostics (optional, can be noisy)
        plot_transition_examples(rows, plots_dir, max_plots=12)
        plot_heatmap_like_scatter(case_summary, plots_dir)

        # If you want: plot only true crossings for alphaQ focus with a smaller cap
        # (plot_only_crossing_cases already caps internally at 20 by default; adjust in that function if desired)

        print(f"[ok] plots written to {plots_dir}")


if __name__ == "__main__":
    main()

