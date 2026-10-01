# Mixed height–slope term: focused diagnostic result

2026-10-01. **The selected bare upper-face mixed term survives the actual tanh
profile and common outgoing bulk dispersion.** The saved coefficient and sign
ingredients establish a nonzero source action at this diagnostic point. This
does not establish excitation by a transverse light mode, a missing term in
the complete closed response, or a leakage value.

The guarded worker completed in **5.247 seconds**, with **77,709,312 bytes
(74.1 MiB)** peak sampled cgroup memory. All 32 saved scalar identity residuals
are exactly zero; native slope omission and sheet-reversal controls respond.
The [completion record](S11c_upstream_mixed_tanh_completion.json) contains the
actual operands, receipt hashes, resources and 115 successful metadata checks.

## Actual selected result

The inputs are strict rest bulk, LAB_HELD/RHO4_CONSTANT, omega=3, c_s0=10,
rho_m=1/10, W0=1, L_W=10 and conserved edge momenta (1/5,1/10). The prescribed
normal-velocity source has profile-direction input momentum zero; output
momentum is Q=1/10. It is not an on-shell transverse slab mode.

The original boundary expansion gives the pressure coefficient per height and
slope amplitude

`-i H omega rho_m q_i / (q_h q_o)`

at zero profile input momentum. The arbitrary-input coefficient, four solved
amplitudes and four original boundary residuals are saved as well. The actual
native flat and first-height/first-slope expressions agree with this derivation.

With the native half-height and half-slope factors routed from `shape_source`,
the coefficient per independent eta*sigma_W is the action

`-(sqrt(3)/2) Integral[ A(t) A(Q-t) / q(t), t over the real line ]`,

where `A(t)=5t/[2 sinh(5 pi t)]`, `A(0)=1/(2 pi)` and
`q(t)^2=1/25-t^2`. The branch is positive real for |t|<1/5 and positive
imaginary outside. The two edge delta functions are factored out with the
native Fourier convention; there is no extra convolution factor of 2 pi.

The recorded evenness identity, positive-transfer flag and removable value
give A(t)>0 on the real line. The prefactor is negative. Hence the propagating
interior contributes a strictly negative real part; the exterior is purely
imaginary and cannot cancel it. Both branch numerators are finite/nonzero,
the endpoint-square identities give integrable inverse-square-root behavior,
and the recorded tail coefficient is 5 with exponential decay. The sign
argument applies to the unsymmetrized convolution; its symmetrized form is
integral-equivalent, not a separate pointwise sign assertion. No integral
magnitude was evaluated.

The tanh height transform's constant delta and principal-value term remain
explicit. Multiplication by the height transfer removes the delta and yields
the regular A(t)/i factor; the slope transform is regular at zero transfer.
This prevents a discarded contact term from being mistaken for cancellation.
Both transfer assignments are included once by the full convolution.

At the control momentum t=1/20, the saved integrand is
`-sqrt(5)/(32 sinh(pi/4)^2)`. Omitting the addressed native tangential-normal
term gives zero; reversing that term or the physical sheet flips the sign.
The two unused edge-tilt components are exactly zero. Reduced kernel units
are M L^-2 T^-1, one length above impedance after the edge deltas are removed.

## Review status and limits

The literal reports remain **Claude NEEDS REVISION / Grok CLEAR FOR THIS
FOCUSED MIXED-GRADE INSTRUMENT**. Both supported the selected mathematics.
The reviewed baseline is `6d645f57`; the local persistence/native-routing
repairs are `c13e23d0`. No fresh independent build/result CLEAR is claimed.
The gate records `independentBuildClearance=false` and the standing user
tooling-fix authority. Actual runtime joins now pass. The unchanged gamma and
hyperbolic reductions also returned exact zero residuals, so the predicted
simplifier failure did not occur.

This is the existing rectangular eta*sigma grade, not a complete physical
second-order shape expansion after substituting sigma_W=eta W0/L_W. The native
three-leg bare slot is literally zero, but its closed permeable resolvent,
reference-pressure conversion and possible downstream accounting were not
computed here. The result therefore does **not** by itself prove the entire
closed operator is incomplete or quantify a correction to the saved benchmark.
The prior selected uniform results, historical benchmark and review debts are
unchanged. No drain or calibrated-light inference follows.

The next bounded question is to trace this term through the existing closed
face response and reference-trace map, checking whether it survives or is
already accounted for, before deciding what needs rebuilding. No such repair,
producer run or defect sweep was launched. This diagnostic is stopped complete.

## Reproduction and containment

Launch command:

`python3 research/pde_ledger_v3/_measurements/S11c_upstream_mixed_tanh_launch.py`

The exact guarded/supervised worker command and hashes are in
`S11c_upstream_mixed_tanh_gate.json`. Runtime root:
`_scratch/s11c/s11c-mixed-tanh-20261001/diagnostic-01`.
Literal worker output is `upstream_mixed_tanh.stdout`, byte-identical to
`complete/checks.json`; scientific stderr is empty. All three journal stages
completed, all 18 source posthashes and 20 snapshots match. The 52 result files
total 99,447 bytes; all inputs, returns and side evidence remain preserved.

Actual containment was 4 GiB cgroup/native inside the 16 GiB pool, zero swap,
one assigned CPU/thread, 32 tasks and a 4 GiB host reserve. All sampled memory
events were zero. RuntimeMaxUSec was infinity and Restart was no; no deadline,
fallback or retry occurred. Completion inspection used JSON/source/hash reads
only, without scientific restoration. Scratch remains ignored and uncommitted.
