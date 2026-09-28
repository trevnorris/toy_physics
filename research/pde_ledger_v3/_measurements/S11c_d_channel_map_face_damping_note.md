# Damping source inspection — stop before repair

2026-09-28. The user approved a bounded repair in principle, but explicitly
required the damping type and size question first, then a go/no-go on the
classification rule. No worker/plan repair, review submission, scientific
payload restoration or new job was performed.

**Type.** The saved inputs activate the physical affinity/permeability channel:
`Lambda_A_0 = 1/100`, `tau_A = 1/10`; `Lambda_V_0 = Lambda_X_0 = 0`.
The supplied law is `J = Lambda_A(omega)*affinity + Lambda_V(omega)*V`, with
`Lambda_I(omega) = Lambda_I_0/(1 - i*omega*tau_I)`.
These are physical interface memory/transfer coefficients, not a numerical
outgoing-wave epsilon. See `S11c_c1_SHARED_PHYSICS.md:146-158` and the native
implementation at `S11c_d_mixing_scattering_sympy_audit.py:3060-3070`.

The separate spatial Fourier/Abel regulators occur in the broader reduction
and finite-profile quadrature. The accepted uniform-source checkpoint records
`regulator: false` for REFERENCE, LEFT and RIGHT; its producer checks this before
acceptance. The raw frequency pencils derive from those saved uniform symbols.
Thus removing an Abel regulator does not remove the physical closure response
in these end pencils. Outgoing bulk radiation is a separate physical effect.

**Size that the inputs actually establish.** Frequencies are in inverse
`T_ref`; the saved memory time is `0.1 T_ref`.

| omega | omega*tau_A |
|---:|---:|
| 0.1 | 0.01 |
| 1 | 0.1 |
| 4 | 0.4 |

The inverse memory time is `10/T_ref`, above this frequency window. This is
**not** a modal damping rate. In the inherited L/T/M schema `Lambda_A_0` has
exponents `(-5,1,1)`, whereas omega has `(0,-1,0)`; dividing the coefficient
0.01 by omega is not a dimensionless damping/frequency comparison. Neither
omega*tau_A nor the kernel coefficient alone establishes weak/strong modal
attenuation. The actual modal decay/frequency ratio, or spatial attenuation
relative to propagation, remains unmeasured across this window and depends on
the coupled end pencil and momentum.

**Implications for the proposed rule, not an adopted rule.** Setting tau to
zero gives the instantaneous kernel `Lambda_A_0`; it does not turn off the
permeability channel. Turning that channel off changes a physical coupling,
not a numerical regularizer. Even then outgoing acoustic radiation can remain,
so a root staying complex need not be material absorption or a closed channel.
The inherited specification explicitly separates bulk radiation from closure
dissipation and warns that Im(k) alone does not identify the mechanism. Strong
attenuation likewise cannot be relabeled entirely as absorption without that
separation. No damping homotopy or new channel definition is selected here.

The accompanying [source inspection](S11c_d_channel_map_face_damping_source_inspection.json)
records exact source hashes, literal equations, metadata and elementary rational
time-scale arithmetic. It contains no new attenuation result. Stop for the
user's decision; the existing build remains NEEDS REVISION and its science
budget is unused.
