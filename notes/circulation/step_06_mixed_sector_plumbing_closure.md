# Step 6 — Mixed-sector plumbing closure and theorem status

## Purpose

The previous steps derived the magnetism-like force law after adding a fixed-current/current-loop closure. This step states the corresponding condition in the actual throat ontology: the mixed sector must supply a plumbing law that maps exterior fluxoid/circulation into an effective localized current or dipole moment.

The relevant microscopic channels are
\[
A_w,
\qquad
F_{\mu w},
\qquad
J^w,
\qquad
E_w=F_{w0},
\qquad
C_a=F_{aw}.
\]
These channels are suppressed only in the strict far-field zero-mode brane limit. They are not erased from the microscopic model.

## Minimal reduced mixed-sector transfer

A one-port reduced Maxwell/mixed block has a static transfer factor of the form
\[
N_0
=
\frac{(\Omega_U^2 g_W+R g_U)^2}
{(\Omega_U^2\Omega_W^2-R^2)^2}.
\]
This is nonnegative when the denominator is nonzero. It can amplify or suppress the outgoing/current-like response, but it does not by itself choose the sign of the circulation source.

The sign still comes from the effective current mapping.

## Minimal plumbing ansatz

Write the effective current/dipole closure as
\[
I_A^{\rm eff}
=\Lambda_A\,\sigma_A\,n_A,
\qquad
\Lambda_A\in\mathbb R
\]
where the sign of \(\Lambda_A\) is part of the unresolved plumbing law. Then the leading coaxial force is
\[
F_d
= -\frac{K}{d^4}\,
\Lambda_1\Lambda_2 N_0\sigma_1\sigma_2 n_1n_2,
\qquad K>0.
\]

The transfer factor is a magnitude factor in this reduced closure. If \(N_0=0\), this channel carries no current-like force at this order. If \(N_0>0\), it cannot flip the sign.

## Facing-mouth result

For two facing mouths,
\[
\sigma_1\sigma_2=-1.
\]
Therefore
\[
F_d
= +\frac{K}{d^4}\Lambda_1\Lambda_2 N_0 n_1n_2.
\]
If \(\Lambda_1\Lambda_2N_0>0\), attraction requires \(F_d<0\), so
\[
\boxed{n_1n_2<0.}
\]
If \(\Lambda_1\Lambda_2N_0<0\), the sign rule reverses. If \(\Lambda_1\Lambda_2N_0=0\), this reduced channel gives no leading force.

Thus, under the additional identical passive mixed/current condition \(\Lambda_1\Lambda_2N_0>0\):
\[
\boxed{\text{opposite local swirl labels attract for facing mouths.}}
\]

## What is still open

The remaining mathematical gap is not the sign algebra. The sign algebra is now clear.

The remaining gap is the PDE-side plumbing law:
\[
\boxed{
\text{derive }\Lambda_A\text{ and its sign from }A_w/F_{\mu w}/J^w\text{ transport on the actual moving-throat branch.}
}
\]

Without that law, the honest theorem is still a no-universal-sign statement. With an identical passive current-like closure that makes \(\Lambda_1\Lambda_2N_0>0\), the magnetism-like facing-mouth rule follows.

## Final conditional force law

For two coaxial facing mouths in the finite-mouth far-field regime,
\[
\boxed{
F_d
=+\frac{3\mu_0\pi R_1^2R_2^2}{2d^4}
\Lambda_1\Lambda_2 N_0 n_1n_2
\left[1-\frac{5}{2}\frac{R_1^2+R_2^2}{d^2}+\cdots\right]
}
\]
where \(F_d<0\) means attraction.

So, if \(\Lambda_1\Lambda_2N_0>0\), opposite local swirl labels attract.

## SymPy audit

The script verifies the square-over-square transfer-factor structure, includes \(N_0\) in the current-like force magnitude, records an explicit square nonnegativity certificate for \(N_0\), and leaves \(\Lambda_A\) real rather than positive. It does not derive the PDE-side plumbing law for \(\Lambda_A\).

Run:

```bash
python step_06_mixed_sector_plumbing_closure_sympy.py
```
