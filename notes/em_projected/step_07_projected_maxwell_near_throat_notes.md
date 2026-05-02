
# Near-Throat Projected Maxwell Notes

## What this extension does

This note pushes the **projection-first Maxwell** derivation toward the regime
that is most relevant to the moving-throat PDE: **close to the throat mouth**,
before the usual far-field zero-mode suppressions are imposed.

The main idea is simple:

- keep the exact projected Maxwell equation,
- work on a **finite throat / half-line domain** instead of only the whole-line brane reduction,
- and take explicit **narrow-kernel limits** to see what the projected theory looks like locally near the mouth.

That exposes the terms that reduction-first Maxwell throws away too early.

---

## 1. Exact projected Maxwell law on a throat domain

Start from the weighted localized bulk equation

\[
\partial_M\!\bigl(Z(w)F^{MN}\bigr)
+\frac{1}{\xi}\partial^N\!\bigl(H(w)B\bigr)
=
\mu_0 J^N,
\qquad
B\equiv \partial\!\cdot\!A.
\]

For a brane index \(\nu\in\{0,1,2,3\}\), project over a throat domain
\(D \subset w\) with a normalized kernel \(W(w)\):

\[
\langle Q\rangle_W := \int_D W(w)\,Q(w)\,dw,
\qquad
\int_D W(w)\,dw=1.
\]

Then the exact projected inhomogeneous law is

\[
\partial_\mu \langle ZF^{\mu\nu}\rangle
+
\bigl[WZF^{w\nu}\bigr]_{\partial D}
-
\langle W' ZF^{w\nu}\rangle
+
\frac{1}{\xi}\partial^\nu \langle H B\rangle
=
\mu_0 \langle J^\nu\rangle.
\]

Equivalently, recombining the boundary and \(W'\) term,

\[
\boxed{
\partial_\mu \langle ZF^{\mu\nu}\rangle
+
\langle \partial_w(ZF^{w\nu})\rangle
+
\frac{1}{\xi}\partial^\nu \langle H B\rangle
=
\mu_0 \langle J^\nu\rangle.
}
\]

This is the best starting point near the throat because it keeps the mixed
sector visible.

---

## 2. The key near-throat term

Define

\[
M^{\mu\nu}(w):=Z(w)F^{\mu\nu}(w),
\qquad
Q^\nu(w):=Z(w)F^{w\nu}(w),
\qquad
G(w):=H(w)B(w).
\]

Then the projected law is

\[
\partial_\mu \langle M^{\mu\nu}\rangle
+
\langle \partial_w Q^\nu\rangle
+
\frac{1}{\xi}\partial^\nu \langle G\rangle
=
\mu_0 \langle J^\nu\rangle.
\]

The crucial extra piece is therefore

\[
\boxed{
\langle \partial_w(ZF^{w\nu})\rangle.
}
\]

That is the projected signature of the mixed \(A_w/F_{\mu w}/J^w\) block.
It vanishes only if the usual far-field assumptions are made:

\[
F_{\mu w}\approx 0,
\qquad
J^w\approx 0,
\qquad
\partial_w A_\mu \approx 0.
\]

So **close to the throat**, projection-first electrodynamics is not just “ordinary
Maxwell with a funny coupling.” It contains a leading transverse-flux-gradient
channel.

---

## 3. Two useful local limits

### 3.1 Symmetric narrow kernel around an interior slice

Take an even normalized kernel \(W_\sigma\) centered at an interior point
\(w_0\). If odd moments vanish and the second/fourth moments are \(m_2,m_4\),
then for a smooth quantity \(Q\),

\[
\langle Q\rangle_\sigma
=
Q_0
+\frac{m_2\sigma^2}{2}Q_2
+\frac{m_4\sigma^4}{24}Q_4
+O(\sigma^6),
\]

\[
\langle \partial_w Q\rangle_\sigma
=
Q_1
+\frac{m_2\sigma^2}{2}Q_3
+\frac{m_4\sigma^4}{24}Q_5
+O(\sigma^6),
\]

where \(Q_n\equiv \partial_w^n Q|_{w_0}\).

So the projected Maxwell law becomes

\[
\partial_\mu
\left[
M^{\mu\nu}_0
+\frac{m_2\sigma^2}{2}M^{\mu\nu}_2
+\cdots
\right]
+
\left[
Q^\nu_1
+\frac{m_2\sigma^2}{2}Q^\nu_3
+\cdots
\right]
+
\frac{1}{\xi}\partial^\nu
\left[
G_0
+\frac{m_2\sigma^2}{2}G_2
+\cdots
\right]
=
\mu_0
\left[
J^\nu_0
+\frac{m_2\sigma^2}{2}J^\nu_2
+\cdots
\right].
\]

**Result:** away from the mouth, the first projection-width correction is
even, \(O(\sigma^2)\).

---

### 3.2 One-sided mouth kernel

Now take the throat mouth to be at \(w=0\) and use the normalized one-sided
kernel on \(w\ge 0\),

\[
W_\ell(w)=\frac{1}{\ell}e^{-w/\ell}.
\]

For a smooth \(Q(w)\),

\[
\langle Q\rangle_\ell
=
Q_0+\ell Q_1+\ell^2 Q_2+\ell^3 Q_3+\cdots,
\]

\[
\langle \partial_w Q\rangle_\ell
=
Q_1+\ell Q_2+\ell^2 Q_3+\ell^3 Q_4+\cdots.
\]

If one keeps the boundary split, then

\[
[WQ]_0^\infty = -\frac{Q_0}{\ell},
\qquad
-\langle W'Q\rangle = \frac{Q_0}{\ell}+Q_1+\ell Q_2+\cdots,
\]

so the apparent \(1/\ell\) singularity cancels exactly, leaving the same
finite derivative expansion.

Thus the near-mouth projected Maxwell law is

\[
\partial_\mu
\left[
M^{\mu\nu}_0+\ell M^{\mu\nu}_1+\ell^2 M^{\mu\nu}_2+\cdots
\right]
+
\left[
Q^\nu_1+\ell Q^\nu_2+\ell^2 Q^\nu_3+\cdots
\right]
+
\frac{1}{\xi}\partial^\nu
\left[
G_0+\ell G_1+\ell^2 G_2+\cdots
\right]
=
\mu_0
\left[
J^\nu_0+\ell J^\nu_1+\ell^2 J^\nu_2+\cdots
\right].
\]

**Result:** right at the mouth, asymmetry turns the first projection-width
correction into **\(O(\ell)\)** rather than \(O(\sigma^2)\).

That is one of the cleanest mathematical signs that the mouth is genuinely
special.

---

## 4. Near-throat zero-mode effective parameters

Now impose the zero-mode ansatz

\[
A_\mu(x,w)=a_\mu(x),
\qquad
A_w=0,
\qquad
F^{w\nu}=0,
\qquad
J^\nu(x,w)=j^\nu(x)S(w).
\]

Then the projected equation becomes

\[
\langle Z\rangle\,\partial_\mu f^{\mu\nu}
+\frac{\langle H\rangle}{\xi}\partial^\nu(\partial\!\cdot\!a)
=
\mu_0 \langle S\rangle\,j^\nu,
\]

so the local projection-first effective parameters are

\[
\mu_{\rm eff}=\mu_0\frac{\langle S\rangle}{\langle Z\rangle},
\qquad
\xi_{\rm eff}=\xi\frac{\langle Z\rangle}{\langle H\rangle}.
\]

### 4.1 Symmetric interior kernel

\[
\mu_{\rm eff}^{(\mathrm{sym})}
=
\mu_0\frac{S_0}{Z_0}
\left[
1+
\frac{m_2\sigma^2}{2}
\left(
\frac{S_2}{S_0}-\frac{Z_2}{Z_0}
\right)
+O(\sigma^4)
\right],
\]

\[
\xi_{\rm eff}^{(\mathrm{sym})}
=
\xi\frac{Z_0}{H_0}
\left[
1+
\frac{m_2\sigma^2}{2}
\left(
\frac{Z_2}{Z_0}-\frac{H_2}{H_0}
\right)
+O(\sigma^4)
\right].
\]

### 4.2 Mouth kernel

\[
\mu_{\rm eff}^{(\mathrm{mouth})}
=
\mu_0\frac{S_0}{Z_0}
\left[
1+
\ell
\left(
\frac{S_1}{S_0}-\frac{Z_1}{Z_0}
\right)
+O(\ell^2)
\right],
\]

\[
\xi_{\rm eff}^{(\mathrm{mouth})}
=
\xi\frac{Z_0}{H_0}
\left[
1+
\ell
\left(
\frac{Z_1}{Z_0}-\frac{H_1}{H_0}
\right)
+O(\ell^2)
\right].
\]

So the first mismatch is controlled by **local shape slippage**:

- symmetric interior:
  \[
  \frac{S_2}{S_0}-\frac{Z_2}{Z_0},
  \qquad
  \frac{Z_2}{Z_0}-\frac{H_2}{H_0},
  \]
- mouth:
  \[
  \frac{S_1}{S_0}-\frac{Z_1}{Z_0},
  \qquad
  \frac{Z_1}{Z_0}-\frac{H_1}{H_0}.
  \]

### Exact special cases

If

\[
H=Z,
\]

then

\[
\boxed{\xi_{\rm eff}=\xi}
\]

exactly, not just perturbatively.
The updated script now keeps the half-line mouth integral
\(\int_0^\infty W_\ell(w)H(w)\,dw\) symbolic and performs the substitution
\(H\mapsto Z\) inside SymPy, so the equality is attached to the actual mouth
projection rather than to two separately typed identical expressions.

If

\[
S=C\,Z,
\]

then

\[
\boxed{\mu_{\rm eff}=C\,\mu_0}
\]

exactly, for any projection kernel.
The same symbolic-substitution check is used for \(S\mapsto CZ\) in the
mouth-kernel source integral.

This contains the earlier global matching channel as a special case:
if \(S=Z/Z_{\rm int}\), then \(\mu_{\rm eff}=\mu_0/Z_{\rm int}\).

---

## 5. Gaussian concrete checks

Take

\[
Z(w)=e^{-w^2/\lambda^2}.
\]

### Symmetric Gaussian observer kernel

With

\[
W_\sigma(w)=\frac{1}{\sqrt{2\pi}\sigma}e^{-w^2/(2\sigma^2)},
\]

the exact average is

\[
\langle Z\rangle
=
\frac{\lambda}{\sqrt{\lambda^2+2\sigma^2}}
=
1-\frac{\sigma^2}{\lambda^2}
+\frac{3\sigma^4}{2\lambda^4}
+O(\sigma^6).
\]

### One-sided exponential mouth kernel

With

\[
W_\ell(w)=\frac{1}{\ell}e^{-w/\ell},
\qquad w\ge 0,
\]

the exact average is

\[
\langle Z\rangle
=
\frac{\sqrt{\pi}\lambda}{2\ell}
e^{\lambda^2/(4\ell^2)}
\operatorname{erfc}\!\left(\frac{\lambda}{2\ell}\right),
\]

with local expansion

\[
\langle Z\rangle
=
1-\frac{2\ell^2}{\lambda^2}
+\frac{12\ell^4}{\lambda^4}
-\frac{120\ell^6}{\lambda^6}
+O(\ell^8).
\]

The accompanying SymPy audit checks this expansion two independent ways: from
the large-argument `erfc` asymptotic form and from direct integration of the
Taylor-expanded Gaussian localizer against the one-sided exponential kernel.

Here the \(O(\ell)\) term vanishes only because this Gaussian localizer has
\(Z'(0)=0\). A generic asymmetric source or gauge profile would produce an
actual linear mouth correction.

---

## 6. What seems new and useful for the PDE

The most useful takeaways are:

### A. The exact throat-domain projection keeps the mixed sector alive

The leading extra term is

\[
\langle \partial_w(ZF^{w\nu})\rangle,
\]

so projection-first electrodynamics near the throat is naturally sensitive to
\(A_w\), \(F_{\mu w}\), and \(J^w\).

### B. The mouth is more sensitive than an interior slice

- interior symmetric observation: first correction \(O(\sigma^2)\),
- mouth-anchored observation: first correction \(O(\ell)\).

So if there is a genuinely throat-local EM effect, the mouth is the place where
it should show up first.

### C. \(H=Z\) remains the clean gauge-sector choice

The near-throat expansion keeps confirming the same thing:

\[
H=Z
\quad\Rightarrow\quad
\xi_{\rm eff}=\xi
\]

exactly, even before far-field reduction.

### D. The first source mismatch is local, not only global

The earlier global matching condition \(S=Z/Z_{\rm int}\) is still valid.
But near the throat, the expansion shows something more local:

- if the source and localization profiles have the same **local shape**, the
  first correction to \(\mu_{\rm eff}\) disappears.

That gives a sharper way to think about how close the near-throat theory is to
the far-field reduction.

---

## 7. How this connects to the moving-throat program

This extension helps on the **Maxwell/mixed** side of the PDE program.

It does **not** by itself compute the full grouped-\(P_2\) conservative bundle,
the outgoing prefactor \(P_0\), or the normalization target

\[
mhat_0^2 P_0 = \frac{54 G c_s^5}{5 a^5 c^5}.
\]

What it does is clarify the local projected EM structure that should feed those
objects close to the throat:

- the mixed transverse derivative term,
- the distinction between symmetric interior and mouth-local observation,
- the clean role of \(H=Z\),
- and the first local shape-mismatch corrections.

So this looks like a plausible missing piece for understanding why the Maxwell /
mixed block can matter much more near the throat than in the far-field brane
limit.

---

## 8. What this still does **not** solve

This is still a **reduced electromagnetic derivation**, not the completed
moving-throat theorem.

In particular, it does not solve:

1. the actual wall dynamics,
2. the full grouped real \(P_2\) overlap integrals on the true branch,
3. the final outgoing normalization,
4. or the parent-action gap that remains until a throat action \(S_\Sigma\) or
   \(S_\eta\) is promoted into the declared total action.

So the safest reading is:

> this extension sharpens the EM/mixed near-throat side of the PDE program, but
> it does not replace the need for the actual moving-throat branch data.

---

## Files

- `step_07_projected_maxwell_near_throat_sympy.py`
