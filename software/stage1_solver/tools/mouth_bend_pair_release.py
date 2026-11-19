"""Phase-locked two-tone preparation of the trapped dipole, starting at rest."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_dipole_bow_decay import generator
from mouth_driven_bend import integrate, saved_row


def run(ws, wd, g, damping, force):
    period = 2*np.pi/ws
    drive_phase = np.exp(1j*wd*period)
    def matrix(t):
        return generator(t, ws, wd, g, 12., 6., damping)
    def homogeneous(t, state):
        return (matrix(t) @ state.reshape(4,4)).ravel()
    propagation = integrate(homogeneous, period, np.eye(4).ravel()).y[:,-1].reshape(4,4)
    def forcing(t, state):
        derivative = matrix(t) @ state
        derivative[3] += 2*force*np.cos(ws*t)*np.exp(1j*wd*t)/12.
        return derivative

    # Accumulate the measured one-period offset across physical drive cycles.
    # Start at rest. Both tones have the same advance at the end of a cycle.
    offset = integrate(forcing, period, np.zeros(4,complex)).y[:,-1]
    drive_map = np.zeros((5,5), complex)
    drive_map[:4,:4] = propagation
    drive_map[:4,4] = offset
    drive_map[4,4] = drive_phase
    augmented = np.array([0,0,0,0,1.],complex)
    driven_periods = 1
    cycle = drive_map.copy()
    while True:
        state = cycle @ augmented
        start = state[:4]/state[4]
        cycle_residual = np.linalg.norm(propagation@start+offset-drive_phase*start)
        if cycle_residual < 1e-14*force:
            break
        cycle = cycle @ cycle
        driven_periods *= 2
        if not np.isfinite(cycle).all():
            raise ArithmeticError("Nonfinite driven period map")

    # Pick the origin of the applied drive phase so the dipole reaches its
    # positive turning point at release. Both members of the pair use it.
    predicted = 2*force/g
    time = period*np.arange(256)/256
    steady = integrate(forcing, period, start, dense=True).sol(time)
    wave_error = max(abs(steady[0]-predicted*np.exp(1j*wd*time)))/predicted
    wall_residual = max(abs(steady[2]))/force
    error = abs(start[1]-1j*wd*predicted)/predicted
    print(f"pair C={damping:g} force_each={force:g} cycles_from_rest={driven_periods} "
          f"max_wall_per_force={wall_residual:.4g} q_release={start[0].real:.12g} "
          f"q_pred={predicted:.12g} relative_wave_error={wave_error:.3g} velocity_error={error:.3g}")
    assert max(wall_residual, wave_error, error)<2e-8
    # The real no-drive trajectory begins at the measured release state; its
    # wall signal alone identifies the slowly decaying dipole component.
    x = start.real.copy()
    trace = []
    for index in range(160):
        trace.append(x[2]/force)
        if index in (0, 1, 2, 4, 16, 64, 128):
            print(f"pair_free C={damping:g} support_cycles={index} time={period*index:.12g} "
                  f"wall={x[2]:.12g} wave={x[0]:.12g}")
        x = propagation @ x
    history = np.lib.stride_tricks.sliding_window_view(np.asarray(trace),5)
    coeff, _, rank, _ = np.linalg.lstsq(history[:,:4],history[:,4],rcond=None)
    roots = np.roots(np.r_[1.,-coeff[::-1]])
    folded = np.angle(roots)/period
    frequency = folded + ws*np.rint((wd-folded)/ws)
    pole = np.argmin(abs(frequency-wd))
    rate = -2*np.log(abs(roots[pole]))/period
    print(f"pair_free_pole C={damping:g} rank={rank} ring_frequency={frequency[pole]:.12g} "
          f"energy_rate={rate:.12g} max_sampled_wall={max(abs(np.asarray(trace)))*force:.12g}")
    assert rank==4 and rate>0 and max(abs(np.asarray(trace)))>.001


def main(manifest):
    record = json.loads(manifest.read_text())
    levels = saved_row(record,"smooth_dipole_l015_01","point n=1 ")
    wall = saved_row(record,"mouth_wave_linewidth_budget01","bow h=0.015 ")
    ws, wd = float(levels["omega_s"]), float(levels["omega_d"])
    g = float(wall["overlap"])*np.sqrt(2/ws)  # archived support action is one
    for damping in (.1,6.5):
        run(ws,wd,g,damping,.001)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spectrum_manifest",type=Path)
    main(parser.parse_args().spectrum_manifest)
