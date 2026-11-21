"""Try to read retained shear from unforced fluctuations on the two saved paths.

The archived state fixes the mean order, not the stochastic physics. In the
benchmark the dynamics supply counts and a separate passive contrast. Only
these recorded signals enter the inversions. A separate driven parcel supplies
the force comparison.
"""
import json
from pathlib import Path

import numpy as np

from mouth_contact_age import FREQUENCIES, PRIOR, emit, readings, trace_path


SITES = 48
PARCELS = 65536
BINS = 16
RESET_RATES = (0., .8)


def prepare(fields):
    path_time, _, _, path_order = trace_path(fields, 2000)
    # The environment is shared across the independent material parcels.
    # Follow the measured mean with piecewise-constant continuous-time on/off
    # generators. The added gross activity here is an illustrative material input.
    times = np.linspace(0., path_time[-1], 16*BINS+1)
    q = np.interp(times, path_time, path_order)
    dt = np.diff(times)
    mid = (q[1:]+q[:-1])/2
    strength = .6+1.8*mid**2
    log_m = np.log(q[1:]/q[:-1])-strength*dt
    lost_slope = -np.expm1(log_m)
    p01 = q[1:] * (-np.expm1(-strength*dt))
    p10 = lost_slope-p01
    if min(p01.min(), p10.min()) < 0:
        raise ValueError("Invalid microscopic rates")
    death_integral = -log_m*p10/lost_slope
    observed_times = times[::16]
    mean = q[::16]
    slopes = np.exp(log_m.reshape(BINS, 16).sum(axis=1))
    birth = mean[1:]-slopes*mean[:-1]
    loss = 1-slopes-birth
    survival = np.exp(-death_integral.reshape(BINS, 16).sum(axis=1))
    return observed_times, mean, birth, loss, survival


def sample_colors(positive, previous, current, keep, rng):
    old_positive = np.zeros_like(previous)
    live = keep > 0
    old_positive[live] = rng.hypergeometric(positive[live], (previous-positive)[live], keep[live])
    return old_positive + rng.binomial(current-keep, .5)


def observe(times, mean, birth, loss, survival, seed):
    rng = np.random.default_rng(seed)
    counts = np.empty((len(times), PARCELS), dtype=np.int16)
    counts[0] = rng.binomial(SITES, mean[0], PARCELS)
    initial_positive = rng.binomial(counts[0], .5)
    colors = [np.empty_like(counts) for _ in RESET_RATES]
    positive = [initial_positive.copy() for _ in RESET_RATES]
    for color in colors:
        color[0] = 2*initial_positive-counts[0]
    for j, dt in enumerate(np.diff(times)):
        k = counts[j].astype(int)
        same_occupancy = rng.binomial(k, 1-loss[j])
        counts[j+1] = same_occupancy+rng.binomial(SITES-k, birth[j])
        for h, reset in enumerate(RESET_RATES):
            # The endpoint occupation can be 1 after an intervening break.
            # A reset restarts the shear reference without changing occupancy.
            p_keep = survival[j]*np.exp(-reset*dt)/(1-loss[j])
            keep = rng.binomial(same_occupancy, p_keep)
            positive[h] = sample_colors(positive[h], k, counts[j+1], keep, rng)
            colors[h][j+1] = 2*positive[h]-counts[j+1]
    return counts, colors


def interval(mean, slopes):
    birth = mean[1:]-slopes*mean[:-1]
    loss = 1-slopes-birth
    if min(slopes.min(), birth.min(), loss.min()) < 0 or slopes.max() > 1:
        raise ValueError("Count data inconsistent with the assumed independent-site transitions")
    return mean[0]/mean[-1]*np.array([np.prod(slopes/(1-birth)), np.prod(1-loss)])


def from_counts(counts):
    mean = counts.mean(axis=1)/SITES
    residual = counts-mean[:, None]*SITES
    slopes = (residual[1:]*residual[:-1]).mean(axis=1)/(residual[:-1]**2).mean(axis=1)
    return interval(mean, slopes)


def from_contrast(color, times):
    # Per-contact contrast is centered in the ensemble by the labeling
    # procedure; fixed scatterer backgrounds must be removed in a real reading.
    final = color[-1].astype(float)
    correlation = (color*final).mean(axis=1)/np.mean(final**2)
    drive = np.exp(-1j*times[:, None]*FREQUENCIES)
    increments = np.vstack([drive[0], np.diff(drive, axis=0)])
    return correlation @ increments / drive[-1]


def driven_force(times, mean, birth, loss, survival, reset, seed):
    # Direct census of extensions in a different, driven parcel. Each cohort
    # shares its shear at birth; strain stays constant between the weak steps.
    rng = np.random.default_rng(seed)
    n = 2_000_000
    cohorts = np.zeros(len(times), dtype=int)
    cohorts[0] = rng.binomial(n, mean[0])
    for j, dt in enumerate(np.diff(times)):
        live = cohorts[:j+1].copy()
        k = live.sum()
        same = rng.binomial(live, 1-loss[j])
        keep = rng.binomial(same, survival[j]*np.exp(-reset*dt)/(1-loss[j]))
        cohorts[:j+1] = keep
        cohorts[j+1] = (same-keep).sum()+rng.binomial(n-k, birth[j])
    end = np.exp(-1j*times[-1]*FREQUENCIES)
    references = np.vstack([np.zeros_like(end), np.exp(-1j*times[:-1, None]*FREQUENCIES)])
    return (cohorts @ (end-references))/(cohorts.sum()*end)


def run_case(fields, seed):
    times, mean, birth, loss, survival = prepare(fields)
    counts, colors = observe(times, mean, birth, loss, survival, seed)
    slope = 1-birth-loss
    emit(test="occupancy", state=fields["state"], intake=fields["intake"], bins=BINS,
         parcel_sites=SITES, parcels=PARCELS, traversal=times[-1], archived_order=[mean[0], mean[-1]],
         sampled_order=(counts[[0,-1]].mean(axis=1)/SITES).tolist(),
         one_bin_exact=interval(mean[[0,-1]], np.array([slope.prod()])).tolist(),
         one_bin_sampled=from_counts(counts[[0,-1]]).tolist(),
         all_bins_exact=interval(mean, slope).tolist(), all_bins_sampled=from_counts(counts).tolist(),
         sample_count_file_hash_input=fields["order_sha256"])
    for i, (reset, color) in enumerate(zip(RESET_RATES, colors, strict=True)):
        pred = from_contrast(color, times)
        measured = driven_force(times, mean, birth, loss, survival, reset, seed+100+i)
        batches = np.stack([from_contrast(batch, times) for batch in np.hsplit(color, 16)])
        se = np.sqrt(np.var(batches.real, axis=0, ddof=1)+np.var(batches.imag, axis=0, ddof=1))/4
        exact = mean[0]/mean[-1]*survival.prod()*np.exp(-reset*times[-1])
        emit(test="contrast", state=fields["state"], reset_per_time=reset, exact_retained=exact,
             omega=FREQUENCIES.tolist(), predicted_normalized_force=readings(pred),
             independently_driven_force=readings(measured), absolute_error=abs(pred-measured).tolist(),
             contrast_complex_standard_error=se.tolist())
        if np.any(abs(pred-measured) > 6*se+.004) or abs(exact-measured[0]) > .004:
            raise ArithmeticError("Passive contrast and direct force are inconsistent")


if __name__ == "__main__":
    archive = json.loads(Path(PRIOR).read_text())["runs"]["mouth_mixed_order_ripple02"]
    cases = [json.loads(s) for s in archive["stdout_lines"]]
    for i, state in enumerate(("near", "deep")):
        case = next(p for p in cases if p["test"]=="archive_patch" and
                    (p["state"], p["mesh"], p["profile"])==(state, 48, "cubic"))
        run_case(case, 4178+i)
