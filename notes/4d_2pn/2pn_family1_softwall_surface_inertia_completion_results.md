
# 2PN Family-1 soft-wall surface inertia completion

## What this step adds

The previous geometry-breathing chain already fixed the **bulk** part of the monopole dynamics from the \(n=5\) Thomas–Fermi interior and reduced the remaining open item to the genuine **soft-wall / surface inertia** completion.

This step derives the leading radial soft-wall correction from the actual Family-1 wall profile
\[
V_{\rm wall}(r;a)=V_0\,S\!\left(\frac{r-a}{d_r}\right)^p,
\qquad
S(x)=\frac{1+\tanh x}{2},
\]
using the same frozen \(n=5\) hydrostatic closure inside the throat.

The result is that the remaining dynamic completion is **not** a new phenomenological term. It is a derived thin-wall correction controlled by the actual wall parameters
\[
\alpha_0 \equiv \frac{V_0}{\mu_c},
\qquad
p,
\qquad
\varepsilon_r \equiv \frac{d_r}{a_0}.
\]

---

## 1. Local \(n=5\) soft-wall profile

At fixed axial slice, the radial hydrostatic/TF wall profile is
\[
f_{\alpha,p}(\xi)=\left(1-\alpha\,S(\xi)^p\right)_+^{1/4},
\qquad
\xi=\frac{r-a}{d_r}.
\]

The turning point is
\[
\xi_*(\alpha,p)=\operatorname{arctanh}\!\bigl(2\alpha^{-1/p}-1\bigr).
\]

So
\[
\xi_*(0)=0 \iff \alpha=2^p.
\]

This is useful physically:

- \(\alpha_0>2^p\): support already turns off **inside** the nominal radius \(a_0\),
- \(\alpha_0<2^p\): support still spills a bit **outside** \(a_0\).

---

## 2. Universal thin-wall moment expansion

Define defect moments against the sharp-wall step:
\[
m_k(\alpha,p)=\int_{-\infty}^{\infty}\xi^k\Bigl[f_{\alpha,p}(\xi)-\Theta(-\xi)\Bigr]\,d\xi.
\]

For the 3-ball cross-section, the relevant radial moments are
\[
J_2=\frac13+\varepsilon_r m_0+2\varepsilon_r^2 m_1+\varepsilon_r^3 m_2+O(\varepsilon_r^4),
\]
\[
J_4=\frac15+\varepsilon_r m_0+4\varepsilon_r^2 m_1+6\varepsilon_r^3 m_2+O(\varepsilon_r^4).
\]

Therefore the cross-sectional mass factor is
\[
\boxed{
R_{\rm mass}=3J_2
=1+3\varepsilon_r m_0+6\varepsilon_r^2 m_1+3\varepsilon_r^3 m_2+O(\varepsilon_r^4),
}
\]
and the radial inertia coefficient becomes
\[
\boxed{
M_{aa}=\frac{J_4}{J_2}
=
\frac35+\frac65\varepsilon_r m_0
+\left(\frac{42}{5}m_1-\frac{18}{5}m_0^2\right)\varepsilon_r^2
+\left(\frac{54}{5}m_0^3-\frac{162}{5}m_0m_1+\frac{81}{5}m_2\right)\varepsilon_r^3
+O(\varepsilon_r^4).
}
\]

So the leading dynamic correction is especially simple:

> the same wall moment \(m_0\) controls both the mass-scale correction and the radial inertia shift at \(O(\varepsilon_r)\).

---

## 3. Axially averaged moments on the filled-to-endcap \(n=5\) branch

On the carried-forward filled-to-endcap \(n=5\) TF branch,
\[
\mu(y)=\mu_c(1-y^2),
\qquad
y=\frac{2w}{L}.
\]

So the local wall ratio becomes
\[
\alpha(y)=\frac{\alpha_0}{1-y^2}.
\]

The relevant averaged wall moments are
\[
\bar m_k(\alpha_0,p)=
\frac{1}{2c_0}
\int_{-1}^{1}(1-y^2)^{1/4}
\,m_k\!\left(\frac{\alpha_0}{1-y^2},p\right)\,dy,
\]
with
\[
c_0=\frac12\int_{-1}^{1}(1-y^2)^{1/4}\,dy.
\]

Then the full radial soft-wall completion of the monopole inertia sector is
\[
\boxed{
\rho_{\rm eff}^{\rm TF+wall}
=
\rho_{\rm eff}^{\rm TF}(5)
\Bigl[
1+3\varepsilon_r\bar m_0
+6\varepsilon_r^2\bar m_1
+3\varepsilon_r^3\bar m_2
+O(\varepsilon_r^4)
\Bigr],
}
\]
\[
\boxed{
\widehat M_{aa}^{\rm TF+wall}
=
\frac35+\frac65\varepsilon_r\bar m_0
+\left(\frac{42}{5}\bar m_1-\frac{18}{5}\bar m_0^2\right)\varepsilon_r^2
+\left(\frac{54}{5}\bar m_0^3-\frac{162}{5}\bar m_0\bar m_1+\frac{81}{5}\bar m_2\right)\varepsilon_r^3
+O(\varepsilon_r^4),
}
\]
while
\[
\widehat M_{LL}=\frac1{14}
\]
stays unchanged at this radial-order completion.

---

## 4. Dynamic monopole response with the wall correction

The carried-forward static geometry Hessian and volume-work coupling are unchanged, so the static monopole closure remains
\[
\Delta K_{00}=\frac{109}{280}.
\]

Only the dynamic inertia changes. The corrected geometry response is therefore
\[
Y_{\rm geom}^{\rm TF+wall}(s)=
\frac{1}{\Sigma_*}\,
\bar g^T\!\left(\widehat H-s\,\widehat M^{\rm TF+wall}\right)^{-1}\bar g,
\]
with
\[
s=\omega^2\,
\frac{\rho_{\rm eff}^{\rm TF+wall}V_0a_0^2}{\Sigma_*}.
\]

So the physical monopole poles are
\[
\Omega_\pm^2(\alpha_0,p,\varepsilon_r)
=
\frac{\Sigma_*}{\rho_{\rm eff}^{\rm TF+wall}V_0a_0^2}
\lambda_\pm\!\left(\widehat M^{\rm TF+wall}\right).
\]

At the EM worked point already fixed in the earlier scripts,
\[
\Lambda_{\rm EM}\approx 1.8474865771,
\qquad
\Sigma_*\approx 0.20761432918,
\]
the sharp-wall TF baseline remains
\[
\lambda_-\approx 5.92556258,\qquad
\lambda_+\approx 237.91117494,
\]
with residues
\[
R_-\approx 0.00262800,\qquad
R_+\approx 0.38665771.
\]

The leading physical pole shifts are
\[
\frac{\Omega_-^2}{\Omega_{-,{\rm sharp}}^2}
=
1-3.40828621\,\varepsilon_r\bar m_0+O(\varepsilon_r^2),
\]
\[
\frac{\Omega_+^2}{\Omega_{+,{\rm sharp}}^2}
=
1-4.59171379\,\varepsilon_r\bar m_0+O(\varepsilon_r^2).
\]

So on the steep-wall branch, where \(\bar m_0<0\), the soft wall **raises** the physical poles.

---

## 5. Representative steep-wall branch

Take the first genuinely steep Family-1 reference branch
\[
p=2,\qquad \alpha_0=10,\qquad \varepsilon_r=0.05.
\]

Then
\[
2^p=4<\alpha_0,
\qquad
\xi_*(0)\approx -0.38558107,
\]
so support already turns off inside the nominal radius at the throat center.

The averaged wall moments are
\[
\bar m_0\approx -0.65067123,\qquad
\bar m_1\approx 0.25044370,\qquad
\bar m_2\approx -0.15585783.
\]

Therefore
\[
\boxed{
R_{\rm mass}\approx 0.90609752,
\qquad
\widehat M_{aa}^{\rm TF+wall}\approx 0.56238115.
}
\]

The corrected dimensionless poles are
\[
\boxed{
\lambda_-\approx 6.00228906,
\qquad
\lambda_+\approx 250.58092901,
}
\]
with residues
\[
\boxed{
R_-\approx 0.00289847,
\qquad
R_+\approx 0.38638724.
}
\]

So the actual physical pole-squared ratios relative to the sharp-wall TF baseline are
\[
\boxed{
\frac{\Omega_-^2}{\Omega_{-,{\rm sharp}}^2}\approx 1.11792424,
\qquad
\frac{\Omega_+^2}{\Omega_{+,{\rm sharp}}^2}\approx 1.16240703.
}
\]

The one-pole Padé reduction remains excellent:
\[
\lambda_{\rm eff}\approx 192.25314580,
\qquad
\max{\rm\ rel.\ err.}\approx 7.84\times 10^{-5}
\]
on the natural low-frequency band \(0\le s\le 0.1\lambda_-\).

So the old breathing auxiliary still survives, but now with a derived soft-wall correction.

---

## 6. Contrast branches

Two useful contrasts from the prototype:

### Stiffer wall, same \(p\)
\[
p=2,\quad \alpha_0=20,\quad \varepsilon_r=0.05
\]
gives
\[
R_{\rm mass}\approx 0.87606211,\qquad
\widehat M_{aa}\approx 0.54983763,
\]
and
\[
\frac{\Omega_-^2}{\Omega_{-,{\rm sharp}}^2}\approx 1.16125672,\qquad
\frac{\Omega_+^2}{\Omega_{+,{\rm sharp}}^2}\approx 1.22438713.
\]

### Sharper power, same central wall ratio
\[
p=4,\quad \alpha_0=10,\quad \varepsilon_r=0.05
\]
gives
\[
R_{\rm mass}\approx 0.99001853,\qquad
\widehat M_{aa}\approx 0.59627288,
\]
and
\[
\frac{\Omega_-^2}{\Omega_{-,{\rm sharp}}^2}\approx 1.01136447,\qquad
\frac{\Omega_+^2}{\Omega_{+,{\rm sharp}}^2}\approx 1.01510709.
\]

So the dynamic wall correction is modest and controllable in the steep Family-1 regime.

---

## Bottom line

This step closes the dynamic wall piece much more tightly:

1. the bulk TF core still fixes the dominant inertia scale,
2. the radial Family-1 soft wall gives a derived thin-wall correction to both the overall mass scale and \(M_{aa}\),
3. the static \(109/280\) monopole closure is untouched,
4. and the dynamic monopole channel remains a positive two-pole Stieltjes response with an excellent one-pole reduction.

So the remaining dynamic “surface inertia” is no longer a placeholder. It is already a controlled function of
\[
\left(\frac{V_0}{\mu_c},\,p,\,\frac{d_r}{a_0}\right).
\]

The next natural tightening is now very specific:

- add the **endcap soft-wall layer** to refine \(M_{LL}\),
- and/or tie this same wall profile directly to the earlier **tangential traction/support** law so the static and dynamic wall sectors come from one boundary-layer model.
