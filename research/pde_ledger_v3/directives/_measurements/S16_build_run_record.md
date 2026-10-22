# S16 build run record: B1–B3 repair

Contract: the supplied sections at `0464b80a`, with B1–B3 as stated in
`S16_build_review_r0_disposition.md`. Reviewed baseline: `b60e38fb`.

Changed files, relative to `research/pde_ledger_v3/`: the two engines
`scripts/S16_moving_region.py` and `mathematica/S16_moving_region.wl`, their
`out/S16_moving_region.out` transcripts, and this record.

Final commands, run serially from `/var/projects/toy_physics`:

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s16/guard/repair_r1_wl_01 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S16_moving_region.wl
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s16/guard/repair_r1_py_02 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S16_moving_region.py
```

Each final guard receipt reports `exitCode: 0`, `limitsVerified: true`,
`childOutcome.exitCode: 0`, `guardReason: null`, and `stderrBytes: 0`.
No admission refusal or guard kill occurred. Each used 50 GiB, zero swap,
32 tasks and no deadline. Child runtime and guard wall time are the receipt's
separate measurements. Peak is maximum sampled cgroup `memory.peak`, including
the supervisor, not isolated CAS RSS. Logs remain under the command's scratch
path; `s16/guard` resolves into `_scratch/s11c/s16`.

| Engine | Child runtime (s) | Guard wall time (s) | Peak cgroup bytes |
|---|---:|---:|---:|
| Wolfram | 6.719703507027589 | 6.786985001002904 | 223035392 |
| SymPy | 163.41480918001616 | 163.47749282000586 | 77463552 |

Final transcript and source SHA256:

- `mathematica/out/S16_moving_region.out`: `72d352df5f631ae5d71772ec992306fd5e1eb0c83843a69f0d65fae526699567`.
  Source: `b106e096dacf0e55e4e1c52ff8bde9a13daebd344c00fe696e2b7387845a599d`.
- `scripts/out/S16_moving_region.out`: `da72cd98811b03a85d766743e0abad0df80b979f7b2d26cc0dfe230f8f5d1435`.
  Source: `8031d3e8e9c90e83837b76a7478bef4f4907aa3149ab6bc888d2eb7321ab6a96`.

Wolfram constructs the balances by differentiating the shell indicator and
integrating spatial flux derivatives by parts. It differentiates each motion
source of the window before multiplying by the integrand. Integrand derivatives
use independent geometric-jet coordinates; window derivatives retain the whole
integrand as a coefficient, including its velocity dependence. It imports nothing
from SymPy.

SymPy solves each supplied balance and integrates by parts by extracting the
spatial flux-derivative coefficients. Rest surface terms contract the resulting
currents with the outward normal. Motion terms multiply the pulled integrand by
the centre or radius normal speed. The Leibniz derivative supplies endpoint
orientation; the generic pullback derivative is a diagnostic, not the source of
the motion terms. Integral derivatives use normal-speed transport and
geometric-jet partial derivatives, with a separate partial-time coordinate.

Both engines derive and print the moving-reference contribution apart from the
integrated material flux. Geometric differentiation is separate from local field
evolution. Each final transcript prints direct derivative operands, sums of
origin pieces, and residual integrands for 44 rate-term derivatives; every such
residual integrand is `0`. No expected-value guard is used. Outputs were not
compared, and no review verdict is claimed.

Each transcript prints the shared normalized domain once as `domain.normalized`,
adding nonzero count to its geometric domain. Individual normalized-object lines
do not repeat that condition. `domain.unnormalized` omits the count restriction.

Unproduced in Wolfram: separately identified stress and advected-momentum
constituents of the supplied total `Pi`; no decomposition law is supplied.
Unproduced in SymPy: the same separately identified constituents, for that reason.
Each prints combined received momentum. No other requested integrated object is
omitted, and no special field closes this gap.

Commits: no new commit; the repair awaits orchestrator review. The unrelated
pre-existing skill edit was left untouched.
