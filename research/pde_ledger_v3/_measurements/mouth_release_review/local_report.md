# Independent local review of r1

**Wall derivation and executable — scoped clear for the normalized, fixed-coefficient scalar E5 wall component.**

**0PN requirements and two-mouth release record — scoped clear for its conditional predictions and identified missing physical inputs.**

The two-unknown traction balance gives `v=(T+γV)/(K+γ)` and `f_wall=γ(T−KV)/(K+γ)`. Opposite wall and carrier forces give slip power `−γ(v−V)²`. Independent signed variations in incident speed, wall speed and both speeds reproduce the two printed slopes and their sum. Equality of the opposite slopes requires `T_U+T_V=K`. Retaining the archived homogeneous equation along the whole isolated-translation branch gives `T_V=0` there.

The revised executable rejected all **18 one-sided controls** with a nonmatching operand in the intended check. These cover an inconsistent traction operator or Rayleigh input, an independently changed exterior operator or archived trace, either force sign, omitted wall work, wrong domain or normal condition, faulty forced trace, loss of forced response, interchanged slopes, an incorrect sum, and zeroing either reported drift observable. The independent result also matches every reported algebraic operand of the frozen baseline.

The corrected aperture requirement includes both carried and stress currents, the rim and the assigned storage/exchange. Equal mass rates cancel the frame shift but can leave a carried jump `ṁ(v⁻−v⁺)`. The single neutral carrier with a shared tangential phase and equal pointwise current cancels it on the specified stationary branch. The independently inspected sources supply brane-side loads, particular wall/coupling contributions and selected benchmark outlets. None selects the total native mouth force.

The corrected release preserves independent initial data and history and permits unsteady pressure to respond. Independent probes at the same prepared state give both columns of the onset response. The nonlinear coupled witness gives a zero first-order residual for `α(0)+χH` and a nonzero quadratic remainder. Absolute release uses both total pair supports. Isolated controls and contributions from other force sectors must be accounted for when extracting the gravity comparison. The two signed source gradients give the candidate's near-rest Newton normalizations under its stated far-zone assumptions.

**Limits.** Physical radius-selected intake, exterior matching, net incident and translation slopes, and the native pair's onset response remain undetermined. The check supplies neither a physical acceleration nor an orbital drag measurement.

**Evidence and reproduction.** `sources.txt` holds exact source commands and readings. `witness.py` accepts a frozen tree or an exact checkout and a disposable directory. It checks all input digests, the recorded stdout, and each mutant’s literal output. The launch below used 1 GiB, zero swap, no deadlines, and exited 0. For reproduction, use a fresh guard-log path. `witness.stdout.txt` is the full 29,671-byte stdout, SHA256 `fe549ee091455ad3e7585284e4a34b4216c6c186cea3a002512a3fea2c5bd96c`.

```sh
python3 scripts/s11c_guarded_run.py --pool mouth-release-local-review-r1 --memory-gib 1 --log-directory _scratch/s11c/mouth_release_review01/local1/guard -- python3 _scratch/s11c/mouth_release_review01/local1/witness.py _scratch/s11c/mouth_release_review01/r1/frozen _scratch/s11c/mouth_release_review01/local1/copies
```
