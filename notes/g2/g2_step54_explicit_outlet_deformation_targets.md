# Step 54 — If the electron sliver exists, it lives in a very small outgoing branch deformation

## Goal

Step 53 fixed the canonical compact outgoing branch itself:
\[
\chi_Q=1,
\qquad
N_Q=1
\]
on the natural point-particle source-map branch. So the only way to keep the
electron-point quartic sliver alive is to deform that outgoing branch.

The later moving-throat notes already give the exact isotropic DtN deformation
algebra for that purpose. Once the canonical even `l=2` fingerprint is preserved,
the full isotropic outgoing deformation reduces to one scalar law:
\[
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
\]
They also classify the robust harmless classes and the genuinely dangerous ones:
pure scale is invisible, and pure scale+argument deformation collapses back to
\(\beta=1\) if the even moments stay canonical. Only a genuine isotropic
throat-core self-energy can move \(\chi_Q\). fileciteturn19file2turn20file18

---

## Step 54A — Exact isotropic deformation law

Write the deformed `l=2` DtN branch as
\[
\Lambda_2^{\rm def}(z)
=
S\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5.
\]

Demanding that the lower even moments stay canonical fixes
\[
\Sigma_2=-\frac{3S\beta^2-3S+\Sigma_0}{9},
\qquad
\Sigma_4=-\frac{3S\beta^4-3S+\Sigma_0}{27},
\]
and leaves
\[
\boxed{
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
}
\]

So the last isotropic outgoing branch data are:

- argument deformation `\beta`,
- static additive core shift `\Sigma_0`,
- odd throat-core outlet `\Sigma_5`.

Overall scale `S` is not itself dangerous; it only enters through the ratios above. fileciteturn19file2turn20file18

---

## Step 54B — The electron target as an exact branch condition

From the carried g-2 chain,
\[
\chi_e=\frac{1}{1+\delta},
\qquad
\delta:=\Lambda_1 f,
\]
with the electron-point value
\[
\delta \approx 3.24737004004\times10^{-4},
\qquad
\chi_e\approx 0.999675368415848.
\]

So the electron sliver is a very small deformation away from the canonical
branch, not a large reorganization of the outgoing sector.

---

## Step 54C — Exact explicit realizations of the electron target

### Pure additive / Robin-like isotropic core

Set
\[
\beta=1,
\qquad
\Sigma_5=0.
\]
Then
\[
\chi_Q=\frac{3S}{3S-\Sigma_0},
\]
so the electron target requires
\[
\boxed{
\Sigma_0=-3S\,\delta.
}
\]

In the simplest normalized case `S=1`, this is just a tiny negative static core shift.

Equivalently, for the pure Robin outlet
\[
\chi_Q^{\rm R}=\frac{3}{3-\rho_R},
\]
the electron target requires
\[
\boxed{
\rho_R=-3\delta.
}
\]

### Compensated Robin–mixed outlet

The compensated hybrid class gives
\[
\chi_Q^{\rm hyb}
=
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
\]
Solving for the electron target yields
\[
\boxed{
\gamma_W=
\frac{\sigma_W+\delta}{9\sigma_W(1+\delta)}.
}
\]

So on the compensated hybrid outlet the electron sliver can be carried by a very
small odd mixed-channel detuning away from the canonical value
\[
\gamma_W=\frac19.
\]
fileciteturn21file7

---

## Step 54D — Linearized tangent law

Expanding around the canonical branch with
\[
S=1+\varepsilon s,
\qquad
\beta=1+\varepsilon b,
\qquad
\Sigma_0=\varepsilon a_0,
\qquad
\Sigma_5=\varepsilon a_5,
\]
gives the exact linearized branch-selection rule
\[
\boxed{
\chi_Q
=
1+
\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)
+O(\varepsilon^2).
}
\]

So to match the leading electron sliver
\[
\chi_e = 1-\Lambda_1 f + O(f^2),
\]
one needs
\[
\boxed{
5b+\frac{a_0}{3}+9a_5=-\Lambda_1.
}
\]

That is the cleanest first-order branch-selection target for the anomaly.

---

## What Step 54 changes

This step sharpens the status from “one remaining branch datum” to:

- the canonical outgoing branch is exact and robust;
- the electron-point sliver, if kept, is a tiny isotropic deformation away from it;
- and that deformation can already be expressed in exact outlet variables.

So the model is no longer missing a broad unknown mechanism.
It is missing only the actual branch-selection law that decides whether the
isotropic outgoing branch stays canonical or carries one of these small explicit
deformations.
