# Selected mixed term survives the closed face and reference maps

2026-10-01. **The saved bare mixed height–slope term produces a nonzero change
in the selected upper-face pressure response, and that change survives the
reference-pressure conversion.** It is not absorbed into the existing iterated
first-shape contribution. This result is for prescribed normal velocity at the
saved point, not transverse light excitation, a complete slab correction or loss.

The guarded worker completed in **8.037 seconds**, with **154,464,256 bytes
(147.3 MiB)** peak sampled cgroup memory. Inspection verified 413 metadata checks,
74 unchanged source pins, 76 exact snapshots, both journal operations and all
55 literal scalar entries of the saved zero residuals. Scientific stderr is
empty and stdout is byte-identical to checks.json. All prior files remain intact.
See [completion evidence](S11c_upstream_mixed_trace_completion.json).

## Actual comparison

The point remains strict rest bulk, LAB_HELD/RHO4_CONSTANT, omega=3, c_s0=10,
rho_m=1/10, W0=1, L_W=10, conserved edge momenta (1/5,1/10), profile input 0
and output 1/10. Physical permeability and memory are retained. Eta and sigma_W
remain independent, with the common wave amplitude factored out.

Let D denote the whole inherited bare coefficient at this selected transfer,
including its conserved-edge delta factors before reduction. The old native
three-leg impedance has a zero [0,2] direct slot. The diagnostic inserts
eta*sigma_W*D once and otherwise holds the original flat/first-shape entries
fixed. Its actual saved difference has only that corner entry:

`delta P_face = eta*sigma_W*D * F`, per unit prescribed V, with

`F = (60 + 91 i) / [105 + 241 i + sqrt(3) (30 + 250 i)]`.

Both external response denominators are finite and nonzero. The factor has no
free middle momentum, depth, shape or integration variable. The saved exact
factor norm is positive. Thus it multiplies the inherited nonzero whole tanh
convolution without another middle integration. The existing first-shape
composition is present once in each compared response and cancels from their
difference. Both original boundary-closure residuals and the resolvent-difference
identity are saved as exact zero matrices.

The native upper-face trace is `p + eta*w1*dp/dw/2` at W0=1. Its saved value
coefficient is 1 and its profile-independent jet coefficient is 0. The computed
mixed reference-pressure difference is identical to the physical-face difference;
the lab-normal jet difference is `i*sqrt(3)/10` times it. The trace reconstruction,
height reconstruction, trace-difference and normal-jet identities all pass.
These are dependent outputs of one conversion, not independent physics tests.

All factor/integrand labels explicitly say **per unit prescribed normal velocity**.
The full native source, including its chemical term, is saved; the chemical-channel
response is not reported. No incoming slab eigenmode is constructed.

The addressed direct-slot omission gives exact zero; sign reversal changes the
sign. Dropping the input-side response factor gives a nonzero saved movement.
Both outgoing middle-branch domain certificates have positive real parts away
from their endpoints. The inherited endpoint/convergence evidence is reused,
not recalculated. The full/reduced kernel units join with exactly two factored
edge deltas; no delta(0), extra 2pi or second convolution is introduced.

## Consequence and limits

The earlier bare-term result now extends through these **selected one-face
closure and reference-coordinate maps**. Those maps do not remove the addition.
Deciding its effect on the closed slab still requires the actual pressure/jet
consumer contraction, both faces where applicable, and the relevant source or
mode. Neither transverse excitation nor a change to the saved benchmark follows
from this diagnostic alone. No producer, production operator or defect solve
was changed or rerun. The next bounded decision is that consumer-level check.

A future implementation must not pass an already integrated D unchanged into
c2 `kernel_apply(second=...)`, whose existing second-slot path integrates the
middle momentum again. This run provides no production routing implementation.
Pure eta²/sigma² terms, full second-order shape, calibrated speed ratio, drain,
loss/current magnitude, Green/FORM/A11/A12 and radiating-witness claims stay open.

## Review and execution provenance

Both literal build reports remain **NEEDS REVISION**. Both supported the selected
mathematics and identified a Python namespace mismatch before computation.
Reviewed baseline and reports are preserved at `8305e075`. Local namespace,
evidence-ordering and per-unit-V label repairs are in `9d49ec43`; eleven stdlib
tests pass, including the old NameError and corrected native-fragment execution
using stand-ins. All scientific assignment ASTs remain identical. No fresh
independent CLEAR or result review is claimed. See the review disposition and
repair record; author disposition is not independent clearance.

Launch command:

`python3 research/pde_ledger_v3/_measurements/S11c_upstream_mixed_trace_launch.py`

The exact command, worker/helper/review hashes and user authority are in the gate.
Runtime root: `_scratch/s11c/s11c-mixed-trace-20261001/diagnostic-01`.
The worker's literal output is `upstream_mixed_trace.stdout`; complete operands
and returns are in `complete/` (39 files, 270,090 bytes).

Actual containment was 4 GiB native/cgroup in the 16 GiB pool, zero swap and memory
events, one CPU/thread, 32 tasks, 4 GiB host reserve and desktop-managed priority.
RuntimeMaxUSec=infinity, Restart=no; no deadline, fallback or retry. The local hook
was armed before launch. Completion inspection used only JSON, source, hashes
and opaque bytes. The diagnostic is stopped complete; no further job launched.
