# S16 build run record

Contract: the user's supplied sections of `S16_decision_list.md` at `0464b80a`.
Deliverable: independently constructed SymPy and Wolfram engines, final-code
transcripts from serial guarded runs, and this run record.

Files (relative to `research/pde_ledger_v3/`):

- `scripts/S16_moving_region.py`; `scripts/out/S16_moving_region.out`.
- `mathematica/S16_moving_region.wl`; `mathematica/out/S16_moving_region.out`.
- `directives/_measurements/S16_build_run_record.md`.

Final commands, executed from `/var/projects/toy_physics`:

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s16/guard/wl_author_03 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S16_moving_region.wl
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s16/guard/py_author_02 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S16_moving_region.py
```

Each final guard receipt reports `exitCode: 0`, `limitsVerified: true`,
`childOutcome.exitCode: 0`, `guardReason: null`, and `stderrBytes: 0`.
No admission refusal or guard kill occurred. Limits were 50 GiB, zero swap,
32 tasks and no deadline. Runtime below is the guard's child wall time;
peak is the maximum recorded cgroup `memory.peak`, including the supervisor,
not an isolated kernel RSS. Guard logs remain in the command's scratch directory.
The existing `s16/guard` symlink resolves into `_scratch/s11c/s16`.

| Engine | Child runtime (s) | Guard wall time (s) | Peak cgroup bytes |
|---|---:|---:|---:|
| Wolfram | 5.670735234976746 | 5.740474468970206 | 223051776 |
| SymPy | 14.319565332960337 | 14.38509699900169 | 71172096 |

SHA256 (the final source and the transcript produced by it):

- Wolfram source: `de23aa97f71d61d4264740efeaf6dd4ef71f76c3154d39dcb377d471aef4c76e`.
- Wolfram transcript: `b33858a1d35147cc3d2ecd307b693a9d1d6e17f85f15341af7c61e0f1764b698`.
- SymPy source: `12c291b99cd6547366baf6bdb5adc4d5856c0f75f0c4c548e648fa824ee78590`.
- SymPy transcript: `893063aa36a70516546160c389c9b7b2737e139418e066c3800c3af325827529`.

The Wolfram source was written first and uses weak differentiation of a shell
indicator; it imports no SymPy code or output. SymPy uses moving-limit Leibniz
differentiation and product-rule integration by parts. Outputs were not compared.

Each transcript prints the count, dipole, second-moment tensor, momentum,
material and momentum surface terms, material face sources, integrated force,
relative and absolute centre kinematics, and received total momentum separately
from momentum rate and surface motion. Derivatives of the rate terms also print
local, centre-motion and radius-motion contributions. Named rate terms have
printed definitions and derivatives; ungoverned field rates remain derivatives.
Normalized objects carry `N != 0`; unnormalized objects do not. Dimensional
operands and residuals are printed, with additive zero left dimension-unspecified.

Unproduced in both engines: separately identified stress and advected-momentum
constituents of `Pi`. The supplied total flux has no decomposition law. The
combined received-momentum object is printed using that total flux and `f`.
No other requested integrated object is omitted. No special field was substituted
to close this gap. No comparison, review, or physical acceptance is claimed.

Commits: none; the build remains uncommitted for the orchestrator's review under
`CLAUDE.md` G3. The pre-existing unrelated skill edit was left untouched.
