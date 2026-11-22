"""Measure one-transit shear retention and neighbor survival in archived sample 87."""
from math import sqrt
from statistics import covariance, mean, variance

from scipy.stats import t

from mouth_bulk_shear_precision import emit, frames, manifest


def retained_fraction(prompt, late):
    """Invert the paired t interval for mean(late - r * prompt) containing zero."""
    n = len(prompt)
    factor = t.ppf(.975, n - 1) ** 2 / n
    p, l = mean(prompt), mean(late)
    a = p * p - factor * variance(prompt)
    b = p * l - factor * covariance(prompt, late)
    c = l * l - factor * variance(late)
    if a <= 0:
        raise ValueError("Prompt mean is unresolved; confidence set need not be bounded")
    root = sqrt(b * b - a * c)
    return [(b - root) / a, (b + root) / a]


def neighbor_memory(rows):
    n = rows[0]["initial"]
    lost = [row["ever_missing"] for row in rows]
    assert n > 0 and lost[0] == 0 and lost == sorted(lost)
    assert all(row["initial"] == n and 0 <= row["original_missing"] <= m <= n
               for row, m in zip(rows, lost, strict=True))
    # Discrete first detection, censored at the last observation. For a contact
    # first missing at frame j, exactly j preceding frames count toward its time.
    dt = rows[-1]["elapsed"] / (len(rows) - 1)
    assert all(abs(row["elapsed"] - j * dt) < 1.e-12 for j, row in enumerate(rows))
    return dict(initial_pairs=n, first_detected_missing_by_end=lost[-1],
                missing_at_end=rows[-1]["original_missing"],
                seen_missing_then_present_at_end=lost[-1] - rows[-1]["original_missing"],
                average_time_to_first_saved_loss_censored_at_transit=
                dt * sum(1 - m / n for m in lost[:-1]))


def read_rung(runs, rung, temp):
    observations = {amplitude: [] for amplitude in (.004, .002)}
    for history in range(3):
        prefix = f"r{rung}h{history}_"
        base, undriven = frames(runs[prefix + "undriven"])
        assert (base["phase"], base["rung"], base["history"], base["temperature"],
                base["sites"], base["applied_pressures"], base["shears"]) == (
                    "undriven", rung, history, temp, 128, [0.], [0.])
        emit("neighbor_memory", temperature=temp, history=history,
             transit=base["transit"], **neighbor_memory(undriven))
        for label, amp in (("shear", .004), ("half_shear", .002)):
            setup, driven = frames(runs[prefix + label])
            for key in ("rung", "history", "temperature", "axis", "transit", "refine",
                        "source_sha256", "dt", "pressure_time", "volume_mobility"):
                assert setup[key] == base[key]
            assert setup["phase"] == "shear" and setup["shears"] == [amp, -amp]
            assert setup["applied_pressures"] == [0., 0.]
            slopes = []
            for row, zero in zip(driven, undriven, strict=True):
                assert row["elapsed"] == zero["elapsed"]
                plus, minus = row["readings"]
                slope = (plus["shear"] - minus["shear"]) / (2 * amp)
                assert abs(slope - row["shear_slope"]) < 1.e-9
                slopes.append(slope)
            observations[amp].append((slopes[0], mean(slopes[-12:])))
    for amp, pairs in observations.items():
        prompt, late = zip(*pairs)
        bounds = retained_fraction(prompt, late)
        emit("shear_retention", temperature=temp, strain=amp, histories=len(pairs),
             prompt=prompt, final_quarter=late, ratio_of_means=mean(late)/mean(prompt),
             paired_fieller_95=bounds)


if __name__ == "__main__":
    records = manifest("mouth_bulk_temperature_steps")["runs"]
    for rung, temperature in ((0, .008), (4, 1.024)):
        read_rung(records, rung, temperature)
