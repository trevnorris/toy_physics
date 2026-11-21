"""Test projected bend properties with the measured complex one-arm carrier load."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, minimize_scalar
from scipy.signal import find_peaks

from mouth_carrier_bend_budget import settings

MEASUREMENTS = Path("research/pde_ledger_v3/_measurements")


def fields(line):
    return dict(item.split("=", 1) for item in line.split() if "=" in item)


def spectrum(lines, ws, wd, coupling, gap, wall_mass, stiffness, wave=None):
    # Carrier is accounted once, at the actual two sidebands, including its
    # reactive part. The selected dipole is treated explicitly.
    kernels = []
    for tone in ("difference","sum"):
        rows = sorted((fields(s) for s in lines if s.startswith("tone="+tone+" ")),
                      key=lambda r:float(r["offset"]))
        offset, omega, re, power = (np.array([float(r[k]) for r in rows])
                                     for k in ("offset","omega","induced_real","power_per_b2"))
        assert np.all(np.diff(offset)>0) and np.all(power>0)
        kernels.append((PchipInterpolator(offset,re,extrapolate=False),
                        PchipInterpolator(offset,2*power/omega,extrapolate=False)))
    def absorb(x):
        z = wd+x
        diagonal = []
        for (real,imag),nu in zip(kernels,(z-ws,z+ws)):
            diagonal.append(stiffness-wall_mass*nu**2-real(x)-1j*imag(x))
        off = 0
        if wave is not None:
            direct,central,lower,upper_re,upper_im,loss = wave
            spectator = central(x)
            diagonal[0] -= spectator+lower(x)
            diagonal[1] -= spectator+upper_re(x)+1j*loss*upper_im(x)
            off = direct/2-spectator
        low,high = diagonal
        # (1,1) D^-1 (1,1)^T for the two coupled bend sidebands.
        self_energy = coupling**2/4*(low+high-2*off)/(low*high-off**2)
        return (1/(wd**2-z**2-self_energy)).imag
    grid = np.linspace(offset[0],offset[-1],6001)
    y = absorb(grid)
    maxima = find_peaks(y)[0]
    i = maxima[np.argmax(y[maxima])]
    peak = minimize_scalar(lambda x:-absorb(x),bounds=(grid[i-1],grid[i+1]),
                           method="bounded",options={"xatol":1e-14}).x
    half = absorb(peak)/2
    left_ids = np.flatnonzero((grid<peak)&(y<half))
    right_ids = np.flatnonzero((grid>peak)&(y<half))
    assert left_ids.size and right_ids.size, "Half-height outside measured band"
    j,k = left_ids[-1],right_ids[0]
    lo = brentq(lambda x:absorb(x)-half,grid[j],peak)
    hi = brentq(lambda x:absorb(x)-half,peak,grid[k])
    width = hi-lo
    # Historical ideal selective readout, equal areas and a 1:3 width ratio.
    x = np.linspace(-gap-8*width,8*width,40001)
    a = width/2
    count = len(find_peaks(a/(x*x+a*a)+3*a/((x+gap)**2+9*a*a))[0])
    return peak,width,(hi-peak)/(peak-lo),count


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--projection",type=Path,required=True)
    ap.add_argument("--extension",type=Path,required=True)
    ap.add_argument("--wave-response",type=Path,required=True)
    ap.add_argument("--threshold-only",action="store_true")
    ap.add_argument("--cells",type=int,choices=(48,64),required=True)
    args = ap.parse_args()
    point,energy,bend,_,archive = settings()
    prior = json.loads((MEASUREMENTS/"mouth_open_reservoirs.json").read_text())["runs"]
    key = {48:"open_reservoir_lines48",64:"final_open_lines64"}[args.cells]
    lines = prior[key]["stdout_lines"]+args.extension.read_text().splitlines()
    proj = args.projection.read_text().splitlines()
    geo, = (fields(s) for s in proj if s.startswith("geometry "))
    wave, = (fields(s) for s in proj if s.startswith("support "))
    ws,wd = float(point["omega_s"]),float(point["omega_d"])
    g = float(bend["overlap"])*np.sqrt(2*archive["parameters"]["I_wave"]/ws)
    gap = abs(float(energy["U_over_action"]))
    static, = (float(fields(s)["static_induced"]) for s in lines if s.startswith("spectrum "))
    direct_K = float(wave["direct_K"])
    tension = float(geo["tensile_both"])
    wave_rows = [fields(s) for s in args.wave_response.read_text().splitlines() if s.startswith("wave ")]
    offsets = [float(r["offset"]) for r in wave_rows]
    dispersion = [PchipInterpolator(offsets,[float(r[k]) for r in wave_rows],extrapolate=False)
                  for k in ("central","lower","upper_real","upper_imag")]
    if args.threshold_only:
        # The fold is the double stationary point of the specified unit-area
        # Lorentzian pair. Its location is independent of both wall masses.
        fold = .503226486865
        interval = (float(geo["uniform_wall_mass"]),float(geo["central_wall_mass"]))
        def clearance(m):
            result = spectrum(lines,ws,wd,g,gap,m,tension+direct_K,(direct_K,*dispersion,1))
            return result[1]/gap-fold
        sampled = [clearance(m) for m in np.linspace(*interval,25)]
        assert sampled[0]>0>sampled[-1] and np.all(np.diff(sampled)<0)
        threshold = brentq(clearance,*interval)
        print(f"threshold cells={args.cells} bracket={interval[0]:.12g},{interval[1]:.12g} "
              f"fold_FWHM_over_gap={fold:.12g} wall_mass_at_fold={threshold:.12g} "
              f"fold_FWHM={fold*gap:.12g} branch_samples={len(sampled)}",flush=True)
        return
    cases = [("archived_trial",12.,6.+static,None,"archived")]
    for label,mass in (("uniform_sheet",geo["uniform_wall_mass"]),("central_sheet",geo["central_wall_mass"])):
        for loss in (0,1):
            cases.append((label,float(mass),tension+direct_K,(direct_K,*dispersion,loss),
                          "carrier_and_scalar" if loss else "carrier_only"))
    for label,mass,k,wave_matrix,channel in cases:
        shift,width,asymmetry,count = spectrum(lines,ws,wd,g,gap,mass,k,wave_matrix)
        if label=="archived_trial":
            original, = (float(fields(s)["FWHM"]) for s in lines if s.startswith("spectrum "))
            assert abs(width-original)<1e-8
        print(f"peaks cells={args.cells} sheet={label} loss={channel} M_wall={mass:.12g} K_diagonal={k:.12g} "
              f"shift={shift:.12g} FWHM={width:.12g} gap={gap:.12g} width_over_gap={width/gap:.12g} "
              f"right_over_left={asymmetry:.9g} selective_count={count}",flush=True)


if __name__=="__main__":
    main()
