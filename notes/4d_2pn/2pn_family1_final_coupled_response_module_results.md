
# 2PN Family-1 final coupled response module

## What this step closes

This step does two things at once.

First, it goes **beyond the separated sidewall-plus-endcap approximation** and evaluates the full
coupled Family-1 Thomas–Fermi wall profile directly in the `(x,y)` throat coordinates.

Second, it packages that coupled wall completion together with the already-passed exact
support/source law and odd/even port structure into a **single low-frequency Family-1 throat-response module**
that reproduces the conservative added 2PN cross sector exactly at zero frequency.

So this is the closest thing yet to the actual finish line.

---

## 1. Exact static support/source law remains unchanged

The local static even-wall law is still

\[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\qquad
t(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4,
\]

with exact source profile

\[
S(\mu)=10-\frac{63}{2}z_{\rm base}(\mu)
      =\frac{7}{16}+\frac{45}{16}\mu^2.
\]

That part does not move.

---

## 2. Final throat module reproduces the exact added 2PN cross block

Using the already-fixed odd dipole wake plus the canonical \(P_0 \oplus P_2\) even ports, the final low-frequency module is

\[
L^{\rm add}_{\rm odd}
=
\frac12(v_A^2+v_B^2)L^{\rm wake}_{1PN}
-\frac{15}{4}(U_A+U_B)(\mathbf v_A\!\cdot\!\mathbf v_B),
\]

\[
L^{\rm add}_{\rm even}
=
\Pi_A\!\cdot\!\Pi_B
+
U_A\,J\!\cdot\!\Pi_B
+
U_B\,J\!\cdot\!\Pi_A
+
\bigl(J\!\cdot\!J-\Delta_{\rm geom}\bigr)U_AU_B,
\]

with

\[
J=\left(\frac{4}{\sqrt5},\frac54\right),
\qquad
\Delta_{\rm geom}=\frac{281}{80}.
\]

The symbolic residual against the exact solved added 2PN comparable-mass target vanishes identically.

So the constructive zero-frequency 2PN cross sector is now closed.

---

## 3. Direct coupled Family-1 full-profile wall completion

On the balanced thin-layer-consistent reference branch

\[
\alpha_r=10,\qquad \varepsilon_r=0.05,\qquad p_r=2,
\]
\[
\chi_{\rm cap}=4,\qquad \varepsilon_z=0.05,\qquad
\alpha_{\rm cap}=4\varepsilon_z\chi_{\rm cap}=0.8,\qquad p_z=2,
\]

use the full coupled TF profile

\[
\widetilde\rho(x,y)=
\Bigl[
1-y^2
-\alpha_r\,S\!\left(\frac{x-1}{\varepsilon_r}\right)^{p_r}
-\alpha_{\rm cap}\,S\!\left(\frac{y-1}{2\varepsilon_z}\right)^{p_z}
\Bigr]_+^{1/4},
\]

on \(0 \le x \le 1\), \(0 \le y \le 1\), with symmetry factor \(2\) in \(y\).

The exact coupled integrals are

\[
I_2 = 2\int_0^1\!\!\int_0^1 x^2\widetilde\rho\,dx\,dy,
\qquad
I_4 = 2\int_0^1\!\!\int_0^1 x^4\widetilde\rho\,dx\,dy,
\]
\[
I_w = 2\int_0^1\!\!\int_0^1 \frac{y^2}{4}x^2\widetilde\rho\,dx\,dy.
\]

Normalizing by the sharp-wall filled-to-endcap TF baseline gives

\[
R_{\rm mass}^{\rm full}\approx 0.886313972989725,
\]
\[
\widehat M_{aa}^{\rm full}\approx 0.563114968953987,
\qquad
\widehat M_{LL}^{\rm full}\approx 0.065829228119349.
\]

This is the first direct **coupled** sidewall+endcap check, rather than a separated leading-order composition.

---

## 4. Separated-order approximation is already very good

Compared to the carried-forward separated-order wall completion,

\[
R_{\rm mass}^{\rm sep}\approx 0.8846236634,
\qquad
\widehat M_{aa}^{\rm sep}\approx 0.5623810783,
\qquad
\widehat M_{LL}^{\rm sep}\approx 0.0671965962,
\]

the direct coupled full-profile values differ by only

\[
\frac{R_{\rm mass}^{\rm full}-R_{\rm mass}^{\rm sep}}
     {R_{\rm mass}^{\rm sep}}
\approx +1.91\times 10^{-3},
\]
\[
\frac{\widehat M_{aa}^{\rm full}-\widehat M_{aa}^{\rm sep}}
     {\widehat M_{aa}^{\rm sep}}
\approx +1.30\times 10^{-3},
\]
\[
\frac{\widehat M_{LL}^{\rm full}-\widehat M_{LL}^{\rm sep}}
     {\widehat M_{LL}^{\rm sep}}
\approx -2.03\times 10^{-2}.
\]

So the earlier separated reduction was already excellent for the total mass scale and radial inertia, and still very respectable for the axial inertia.

That is an important closure check.

---

## 5. Final coupled monopole breathing response

Feeding the exact coupled wall data into the already-fixed geometry Hessian at the EM worked point gives

\[
\lambda_-\approx 6.405572392138922,
\qquad
\lambda_+\approx 254.444968136936126,
\]

with positive residues

\[
R_-\approx 0.002552474771738,
\qquad
R_+\approx 0.386733239513976,
\qquad
R_-+R_+=\frac{109}{280}.
\]

So the final monopole channel is

\[
K_{00}(s)=
-\frac{757}{2520}
+\frac{R_-}{1-s/\lambda_-}
+\frac{R_+}{1-s/\lambda_+}.
\]

At zero frequency,
\[
K_{00}(0)=\frac{4}{45},
\]
exactly as required.

Relative to the sharp-wall TF baseline, the physical monopole pole ratios are

\[
\frac{\Omega_-^2}{\Omega_{-,{\rm TF}}^2}
\approx 1.219665554412172,
\qquad
\frac{\Omega_+^2}{\Omega_{+,{\rm TF}}^2}
\approx 1.206678094538845.
\]

The one-pole Padé reduction is still excellent:

\[
\lambda_{\rm eff}\approx 202.923516367519028,
\qquad
\max {\rm rel.err.}\approx 6.89\times 10^{-5}
\quad\text{on } 0\le s\le 0.1\lambda_-.
\]

So the old single breathing auxiliary remains a very controlled low-frequency reduction of the exact two-pole geometry channel.

---

## 6. Final state of play

At conservative 2PN order, the constructive throat program is now essentially closed in the following sense:

- **exact local static support/source wall law**: closed,
- **exact zero-frequency odd + \(P_0 \oplus P_2\) port reconstruction**: closed,
- **exact algebraic match to the solved added comparable-mass 2PN cross block**: closed,
- **static monopole closure \(109/280\) from geometry Hessian**: closed,
- **bulk TF inertia scale**: closed,
- **radial sidewall and endcap dynamic corrections**: closed,
- **direct coupled full-profile sidewall+endcap check**: now done.

The remaining open quantities are narrower than before:

1. the **non-monopole dynamic pole scales** (\(\ell=1,2\) channels), which are true inner-throat PDE observables rather than conservative 2PN necessities;
2. a fully analytic derivation of the local wall profiles from the microscopic soft-wall traction PDE, rather than from reduced matching.

Those no longer block the conservative 2PN derivation itself.

So the practical finish-line statement is:

\[
\boxed{
\text{the conservative 2PN throat-response module is now effectively in hand.}
}
\]
