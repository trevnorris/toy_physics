# 2PN Family-1 traction-balance / wall-stress completion: current result

## 1. What this step closes

The strict isotropic boundary-layer pullback was too small to reproduce the passed static support sector.
The next natural question was whether the solved support operator could be rewritten as an **explicit local wall-stress / traction-balance surface energy**.

The answer is yes.

The cleanest exact closure is:

- keep the monopole wall channel separate,
  
  \[
  K_{00}=\frac{4}{45},
  \]

- and solve the dipole/quadrupole support sector with one base pressure profile and one tangential wall-stress profile.

---

## 2. Minimal exact traction-balance ansatz on the support sector

On the \(\ell=1,2\) support sector, write

\[
K_{\ell m}=\langle z_{\rm base}\rangle_{\ell m}+\bigl(\ell(\ell+1)-2\bigr)\langle t\rangle_{\ell m},
\]

with low-order Family-1 profiles

\[
z_{\rm base}(\mu)=b_0+b_2\mu^2,
\qquad
 t(\mu)=t_0+t_2\mu^2+t_4\mu^4.
\]

Matching the five carried-forward support targets

\[
K_{1\perp}=\frac27,
\quad
K_{10}=\frac14,
\quad
K_{20}=\frac49,
\quad
K_{21}=\frac23,
\quad
K_{22}=\frac83,
\]

fixes the profiles **uniquely**:

\[
b_0=\frac{17}{56},
\qquad
b_2=-\frac{5}{56},
\]

\[
t_0=\frac{593}{672},
\qquad
t_2=-\frac{1553}{672},
\qquad
t_4=\frac78.
\]

So

\[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\]

\[
t(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4.
\]

In Legendre form,

\[
z_{\rm base}(\mu)=\frac{23}{84}-\frac{5}{84}P_2(\mu),
\]

\[
t(\mu)=\frac{211}{1008}-\frac{129}{112}P_2(\mu)+\frac{7}{18}P_2(\mu)^2.
\]

This is the exact same Family-1 low-order basis that kept reappearing in the earlier flare / Robin / variational steps.

---

## 3. Equivalent local wall-energy form

For

\[
\mathcal O_{\rm supp}=z_{\rm base}+\frac12\{-\Delta_S-2,\,t\},
\]

the equivalent local quadratic form on the \(\ell=1,2\) support sector is

\[
E_{\rm supp}[\Psi]
=\frac12\int d\Omega\,\Big[p(\mu)\Psi^2+t(\mu)|\nabla_S\Psi|^2\Big],
\]

with

\[
p(\mu)=z_{\rm base}(\mu)-2t(\mu)-\frac12\Delta_S t(\mu).
\]

Using

\[
\Delta_S t(\mu)= -\frac{35}{2}\mu^4+\frac{2729}{112}\mu^2-\frac{1553}{336},
\]

this gives the exact pressure profile

\[
p(\mu)=\frac{571}{672}-\frac{5141}{672}\mu^2+7\mu^4.
\]

In Legendre form,

\[
p(\mu)= -\frac{155}{168}-\frac{2005}{1008}P_2(\mu)+\frac{28}{9}P_2(\mu)^2.
\]

So the solved support sector is now an explicit local wall model:

- one base pressure profile \(z_{\rm base}(\mu)\),
- one tangential wall-stress / curvature-compliance profile \(t(\mu)\),
- and the induced pressure profile \(p(\mu)\) that follows from the anticommutator structure.

---

## 4. Exact verification

With the local wall-energy form,

\[
K_{\ell m}=\int d\Omega\Big[p(\mu)|Y_{\ell m}|^2+t(\mu)|\nabla_S Y_{\ell m}|^2\Big],
\qquad \ell=1,2,
\]

all five support-channel residuals vanish identically.

So the dipole/quadrupole support sector is no longer just a passed operator fit. It now has a direct **traction-balance / wall-stress** reading.

---

## 5. Physical reading

This suggests the strict isotropic boundary layer was missing exactly one new structural ingredient:

> an axisymmetric tangential wall-stress / curvature-compliance profile in the same low-order Family-1 basis.

The monopole wall mode remains separate,

\[
K_{00}=\frac{4}{45},
\]

which is consistent with the earlier PDE observation that a finite monopole pole/stiffness is its own channel.

So the constructive picture is now:

1. strict isotropic pullback alone = insufficient,
2. add one explicit tangential wall-stress profile \(t(\mu)\),
3. keep the monopole wall channel separate,
4. the full passed support sector closes exactly.

That is the first explicit local traction-balance completion of the strict isotropic layer that is exact on the solved 2PN support sector.
