# Parent Throat Action — Monopole JSONL Fast-Screen

## Purpose

The current `single_throat_monopole.py` solver emits JSON diagnostics, not full
field snapshots.

That means the exact runtime monitor cannot be applied to its stdout directly.
But we can still get a useful **partial falsification screen** immediately by
reading the live diagnostic stream and checking whether the solver is even close
to the expected \(1/r\) pressure and \(1/r^2\) acceleration behavior.

That screen lives in
[single_throat_monopole_jsonl_fastscreen.py](/home/trevnorris/Downloads/em_projected/single_throat_monopole_jsonl_fastscreen.py).

---

## 1. What it uses

The solver already prints diagnostic events with:

- `mach_max`
- `fits.dP_slope`
- `fits.geff_slope`
- `fits.dP_npts`
- `fits.geff_npts`

The fast-screen reads the last diagnostic event and applies fail-fast checks.

---

## 2. Current thresholds

Default early-kill thresholds are:

- `dP_slope` target `-1 +/- 0.35`
- `geff_slope` target `-2 +/- 0.35`
- `mach_max <= 0.6`
- at least `8` fit points in each slope estimate

Verdict rules:

- `FAIL`: slope targets missed or flow is too compressible
- `INCOMPLETE`: not enough fit support, no diagnostics, or malformed JSONL
  lines were encountered without a harder failure
- `PASS`: all checks satisfied

This is deliberately weaker than the full snapshot-based runtime monitor.

---

## 3. Self-test

[step_37_parent_throat_action_monopole_jsonl_fastscreen_sympy.py](/home/trevnorris/Downloads/em_projected/step_37_parent_throat_action_monopole_jsonl_fastscreen_sympy.py)
checks six cases:

1. a good log -> `PASS`
2. a bad-slope / high-Mach log -> `FAIL`
3. a weakly sampled log -> `INCOMPLETE`
4. a threshold-boundary log -> `PASS`
5. a just-outside-threshold log -> `FAIL`
6. a malformed JSONL log -> `INCOMPLETE`

So the JSONL screen is doing real discrimination, not just parsing lines.

---

## 4. CLI

```bash
python single_throat_monopole_jsonl_fastscreen.py monopole.log --output-json monopole_verdict.json
```

This is the shortest path to falsification for the existing monopole solver as
it is currently written.

It should be treated as a front-end sieve, not a replacement for the full
snapshot-based monitor and fail-fast classifier once field dumps are available.
