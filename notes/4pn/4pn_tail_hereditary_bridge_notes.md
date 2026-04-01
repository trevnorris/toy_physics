# 4PN tail / hereditary bridge — result ledger

## What this step accomplished

This step freezes the **nonlocal 4PN bridge coefficient**.

The local 4PN sector was already assembled exactly in the previous referee step. What remained open was the separate hereditary/tail side. The new audit shows that the 4PN tail coefficient is **not** a new independent datum. It is fixed by the *same* quadrupole normalization that the 2.5PN program had already isolated as the last remaining universal gap.

Equivalently, once the orbital/worldtube STF quadrupole is canonically normalized, the conservative 4PN tail coefficient follows automatically.

---

## GR reference structure

The standard conservative 4PN dynamics splits into

\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}}.
\]

On the action side, the conservative tail term may be written in the time-symmetric form

\[
S_{\mathrm{tail}}
=
\frac{G^2 M}{5c^8}
\operatorname{Pf}_{2s/c}
\iint \frac{dt\,dt'}{|t-t'|}
I_{ij}^{(3)}(t)
I_{ij}^{(3)}(t').
\]

On the Hamiltonian side the corresponding nonlocal term is

\[
H_{\mathrm{tail}}(t)
=
-
\frac{G^2 M}{5c^8}
I_{ij}^{(3)}(t)
\operatorname{Pf}_{2s/c}
\int_{-\infty}^{+\infty}
\frac{dv}{|v|}
I_{ij}^{(3)}(t+v).
\]

So the universal GR tail coefficient is

\[
C_{\mathrm{tail}}^{\mathrm{GR}}=
\frac{G^2 M}{5c^8}.
\]

---

## Canonical bridge to the 2.5PN Burke–Thorne coefficient

The canonical Burke–Thorne local odd quadrupole coefficient is

\[
\gamma_{\mathrm{GR}}=
\frac{2G}{5c^5}.
\]

The new audit verifies the exact identity

\[
\boxed{
C_{\mathrm{tail}}^{\mathrm{GR}}
=
\frac{GM}{2c^3}\,\gamma_{\mathrm{GR}}.
}
\]

This is the cleanest structural statement of the 4PN bridge.

It says that the conservative hereditary coefficient is exactly the 2.5PN quadrupole coefficient multiplied by the standard monopole-scattering factor

\[
\frac{GM}{2c^3}.
\]

So if the 2.5PN quadrupole coefficient is correct, then the 4PN tail coefficient is fixed automatically.

---

## Harmonic-frequency structure of the tail kernel

Using a regulated finite-part proxy for the kernel

\[
2\int_0^{\infty} d\tau\,e^{-\varepsilon\tau}
\frac{\cos(\omega\tau)-1}{\tau}
=
-\ln\!\left(1+\frac{\omega^2}{\varepsilon^2}\right),
\]

the nonlocal tail action reduces on a monochromatic mode to the expected logarithmic form

\[
K_{\mathrm{tail}}(\omega)
\sim
-
\frac{G^2 M}{5c^8}
\omega^6
\ln|\omega|.
\]

Relative to the local Burke–Thorne kernel

\[
K_{\mathrm{BT}}(\omega)\sim \gamma_{\mathrm{GR}}\omega^5,
\]

the ratio is

\[
\frac{K_{\mathrm{tail}}}{K_{\mathrm{BT}}}
\sim
\frac{GM\,\omega}{c^3}\ln|\omega|,
\]

which is the expected 1.5PN tail promotion above the leading 2.5PN quadrupole reaction channel.

---

## Toy-model quadrupole branch version

The 2.5PN quadrupole program already reduced the remaining universal normalization problem to the canonical invariant quantity

\[
\gamma_{\mathrm{quad}}^{\mathrm{eff}}
=
\mathcal N_Q\frac{a^5}{27c_s^5}
=
\overline\Gamma_5
=
9\frac{\overline K_2^{5/2}}{\overline K_0^{3/2}}.
\]

The new bridge result then gives the toy-model 4PN tail coefficient as

\[
\boxed{
C_{\mathrm{tail}}^{\mathrm{toy}}
=
\frac{GM}{2c^3}
\gamma_{\mathrm{quad}}^{\mathrm{eff}}.
}
\]

Equivalently,

\[
C_{\mathrm{tail}}^{\mathrm{toy}}
=
\frac{GM}{2c^3}\,\overline\Gamma_5
=
\frac{9GM}{2c^3}
\frac{\overline K_2^{5/2}}{\overline K_0^{3/2}}.
\]

So the entire hereditary 4PN bridge is controlled by the same normalized quadrupole branch data as the 2.5PN odd channel.

---

## Exact matching target inherited from the 2.5PN normalization problem

Imposing Burke–Thorne normalization gives

\[
\gamma_{\mathrm{quad}}^{\mathrm{eff}}=
\frac{2G}{5c^5},
\]

which is equivalent to

\[
\mathcal N_Q^{\mathrm{target}}
=
\frac{54Gc_s^5}{5a^5c^5}.
\]

Substituting that into the 4PN bridge coefficient yields exactly

\[
C_{\mathrm{tail}}^{\mathrm{toy}}
=
\frac{GM}{2c^3}
\frac{2G}{5c^5}
=
\frac{G^2 M}{5c^8}
=
C_{\mathrm{tail}}^{\mathrm{GR}}.
\]

So the GR 4PN hereditary coefficient follows automatically once the 2.5PN quadrupole normalization gap is closed.

---

## Main theorem statement from this step

Let the canonically normalized orbital/worldtube STF quadrupole branch satisfy the isotropic passive/outgoing low-frequency relations already isolated in the 2.5PN audit. Then:

1. the unique compatible conservative hereditary 4PN coefficient is
   \[
   C_{\mathrm{tail}}=
   \frac{GM}{2c^3}\,\gamma_{\mathrm{quad}}^{\mathrm{eff}},
   \]
2. therefore the 4PN tail sector introduces **no new independent quadrupole normalization datum**,
3. and full GR matching of the hereditary coefficient is equivalent to the already known 2.5PN normalization target
   \[
   \gamma_{\mathrm{quad}}^{\mathrm{eff}}=
   \frac{2G}{5c^5}.
   \]

So the remaining full-4PN gap is still the same narrow passive/outgoing quadrupole normalization problem already identified by the 2.5PN program.

---

## What remains open after this step

This step does **not** close the full conservative 4PN theorem unconditionally.

What remains open is:

1. the final derivation of the passive/outgoing quadrupole normalization from the completed moving-throat PDE,
2. the explicit insertion of that normalized tail module into a full end-to-end 4PN referee/master script.

But the bottleneck is now much sharper:

- the **local** 4PN sector is already closed,
- the **tail** coefficient is now reduced to the same invariant normalization problem as 2.5PN,
- so the next full-4PN theorem step is an assembly step, not a new open-ended search.
