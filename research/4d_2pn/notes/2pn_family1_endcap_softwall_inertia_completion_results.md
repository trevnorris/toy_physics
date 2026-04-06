# 2PN Family-1 endcap soft-wall inertia completion

## What this step adds

The previous chain had already fixed:

- the bulk Thomas–Fermi inertia scale on the filled-to-endcap `n=5` branch,
- the radial Family-1 sidewall correction to the monopole breathing channel,
- and the reduced geometry-Hessian / two-pole breathing closure.

The remaining missing dynamic wall ingredient was the **endcap soft-wall layer**.

This step derives that correction and folds it into the carried-forward monopole response.
The punchline is simple:

1. the endcap layer is **parametrically weaker** than the sidewall,
2. its first nontrivial effect scales as `eps_z^(5/4)` on the frozen `n=5` branch,
3. and once it is included, the current program already has a near-final **full wall-completed monopole breathing branch**, at least to separated leading order.

---

## 1. Why the endcap scaling is different

On the filled-to-endcap TF branch,

axial support already vanishes at the cap:

\[
\rho_0(u) \propto (1-u^2)^{1/4},
\qquad u=\frac{2w}{L}.
\]

So near the right endcap, with
\[
u = 1 + \varepsilon_z x,
\qquad \varepsilon_z = \frac{2 d_z}{L},
\]
we have
\[
1-u^2 = \varepsilon_z\bigl(-2x-\varepsilon_z x^2\bigr).
\]

That means a genuine thin-cap layer cannot have a wall amplitude of order `mu_c`.
To stay localized on an `O(d_z)` layer, the cap must scale as
\[
\frac{V_{\rm cap}}{\mu_c} = 2\varepsilon_z\,\alpha_z\,S(x)^p,
\qquad S(x)=\frac{1+\tanh x}{2}.
\]

This is the key distinction from the sidewall.
Because the TF profile already dies at the endcap, the first wall correction is suppressed by one extra quarter-power.

---

## 2. Local reduced cap profile and universal defect moment

With the scaling above, the local reduced cap profile is
\[
g_{\alpha,p}(x)=\bigl(-x-\alpha S(x)^p\bigr)_+^{1/4}.
\]

Relative to the sharp-cap baseline `(-x)_+^(1/4)`, define the universal defect moments
\[
\nu_k(\alpha,p)=\int_{-\infty}^{\infty}x^k\Bigl[g_{\alpha,p}(x)-(-x)_+^{1/4}\Bigr]dx.
\]

Then the filled-to-endcap TF integrals obey
\[
c_0^{\rm cap}=c_0+2^{1/4}\varepsilon_z^{5/4}\nu_0+O(\varepsilon_z^{9/4}),
\]
\[
c_2^{\rm cap}=c_2+2^{1/4}\varepsilon_z^{5/4}\nu_0+O(\varepsilon_z^{9/4}),
\]
with the carried-forward sharp-cap baseline
\[
c_0=\frac{\sqrt\pi\,\Gamma(5/4)}{2\Gamma(7/4)},
\qquad
c_2=\frac{2}{7}c_0.
\]

So the effective bulk inertia scale and the axial breathing metric become
\[
\boxed{
\frac{\rho_{\rm eff}^{\rm TF+cap}}{\rho_{\rm eff}^{\rm TF}}
=1+A_{\rm cap}\,\nu_0\,\varepsilon_z^{5/4}+\cdots,
}
\]
\[
\boxed{
\widehat M_{LL}^{\rm TF+cap}
=\frac{1}{14}+B_{\rm cap}\,\nu_0\,\varepsilon_z^{5/4}+\cdots,
}
\]
with exact coefficients
\[
A_{\rm cap}=\frac{2^{1/4}}{c_0}
=\frac{6\,2^{1/4}\Gamma(3/4)}{\sqrt\pi\,\Gamma(1/4)}
\approx 1.3606190066912236,
\]
\[
B_{\rm cap}=\frac{5\,2^{1/4}}{28c_0}
\approx 0.24296767976628994.
\]

This is the clean main scaling result:

> on the frozen `n=5` branch, the endcap wall correction is `O(eps_z^(5/4))`, not `O(eps_z)`.

---

## 3. Representative steep-cap branch and direct full-profile check

Take the representative cap branch
\[
\varepsilon_z=0.05,
\qquad \alpha_z=1,
\qquad p_z=2.
\]

The local turning point solving
\[
x+\alpha S(x)^p=0
\]
is
\[
x_*\approx -0.1720646550263600,
\]
and the universal defect moment is
\[
\nu_0\approx -0.1297171210945550.
\]

Using the asymptotic formulas above gives
\[
c_0^{\rm asym}\approx 0.8703719198752305,
\qquad
c_2^{\rm asym}\approx 0.2460725021866305,
\]
\[
\widehat M_{LL}^{\rm asym}\approx 0.0706802737334132.
\]

Direct numerical integration of the full soft-cap profile gives
\[
c_0^{\rm full}\approx 0.8703584556098921,
\qquad
c_2^{\rm full}\approx 0.2461237735473917,
\]
\[
\widehat M_{LL}^{\rm full}\approx 0.0706960942244548.
\]

So already at `eps_z = 0.05` the leading asymptotic is very good:

- `c0` relative error `~ 1.55 × 10^-5`,
- `M_LL` relative error `~ 2.24 × 10^-4`.

This makes the cap reduction numerically trustworthy in the same regime used by the earlier sidewall work.

---

## 4. Dynamic monopole response with the endcap correction

On the carried-forward EM worked point,
\[
\Lambda_{\rm EM}=\frac{\sqrt2\pi}{x_{01}}\approx 1.8474865771,
\qquad
\rho=\frac{1}{10},
\qquad
\beta=12,
\qquad
\Sigma_*\approx 0.20761432918,
\]
the sharp-wall TF baseline had
\[
\lambda_-\approx 5.92556258,
\qquad
\lambda_+\approx 237.91117494,
\]
with residues summing to the carried-forward monopole target
\[
R_-+R_+=\frac{109}{280}.
\]

The new cap branch gives
\[
R_{\rm cap}=\frac{c_0^{\rm cap}}{c_0}\approx 0.99581161464,
\qquad
\widehat M_{LL}^{\rm cap}\approx 0.07069609422.
\]

Feeding that into the same reduced geometry response yields
\[
\lambda_-\approx 5.97431790,
\qquad
\lambda_+\approx 238.41448955,
\]
with residues
\[
R_-\approx 0.00258517,
\qquad
R_+\approx 0.38670054,
\qquad
R_-+R_+=\frac{109}{280}.
\]

The physical pole-squared ratios relative to the sharp-wall TF baseline are
\[
\boxed{
\frac{\Omega_-^2}{\Omega_{-,\rm sharp}^2}\approx 1.01246857,
\qquad
\frac{\Omega_+^2}{\Omega_{+,\rm sharp}^2}\approx 1.00633046.
}
\]

So the endcap correction is real, but modest — exactly what the `eps_z^(5/4)` scaling suggested.

The one-pole Padé reduction remains excellent:
\[
\lambda_{\rm eff}\approx 189.46282891,
\qquad
\max \text{ rel. err.}\approx 6.98\times 10^{-5}
\]
on the natural low-frequency band `0 ≤ s ≤ 0.1 lambda_-`.

---

## 5. Leading separated-order full-wall composite branch

Now combine:

- the carried-forward radial sidewall result from the previous step,
  \[
  R_{\rm side}\approx 0.90609752477,
  \qquad
  \widehat M_{aa}\approx 0.56238115491,
  \]
- with the new endcap result,
  \[
  R_{\rm cap}\approx 0.99581161464,
  \qquad
  \widehat M_{LL}\approx 0.07069609422.
  \]

To leading separated order,
\[
R_{\rm full}=R_{\rm side}R_{\rm cap}\approx 0.90230243917,
\]
\[
\widehat M_{\rm full}=
\begin{pmatrix}
0.56238115491 & 0\\
0 & 0.07069609422
\end{pmatrix}.
\]

The resulting full-wall composite breathing branch is
\[
\lambda_-\approx 6.05235326,
\qquad
\lambda_+\approx 251.08293474,
\]
with residues
\[
R_-\approx 0.00285529,
\qquad
R_+\approx 0.38643043,
\qquad
R_-+R_+=\frac{109}{280}.
\]

The corresponding physical pole-squared ratios relative to the sharp-wall TF baseline are
\[
\boxed{
\frac{\Omega_-^2}{\Omega_{-,\rm sharp}^2}\approx 1.13198989,
\qquad
\frac{\Omega_+^2}{\Omega_{+,\rm sharp}^2}\approx 1.16963464.
}
\]

The one-pole Padé reduction is still excellent:
\[
\lambda_{\rm eff}\approx 193.59552541,
\qquad
\max \text{ rel. err.}\approx 7.72\times 10^{-5}.
\]

So this is the first actual **full wall-completed monopole breathing branch** in the current Family-1 program, at least to separated leading order in the sidewall and endcap thickness parameters.

---

## Bottom line

This step narrows the remaining finish-line task substantially.

We now have:

1. **static local wall support/source law** from the Family-1 traction sector,
2. **static monopole closure** from the reduced geometry Hessian,
3. **bulk TF inertia** from the parent PDE,
4. **radial sidewall dynamic correction** from the Family-1 radial soft wall,
5. **endcap dynamic correction** from the new filled-to-endcap cap-layer reduction,
6. and a first **full wall-completed dynamic monopole branch** by combining the carried-forward sidewall result with the new cap result.

So the remaining gap is no longer “derive the dynamic wall sector.”
It is much narrower:

> derive the **fully coupled** sidewall-plus-endcap boundary-layer reduction beyond separated order, and then package the whole wall sector into one final throat-response module.
