# Moving-Throat PDE — Stage 246: Compensated Multimode Source Compiler Beyond Positive Family-1

## Status

**Exact within the carried positive-source Family-1 mouth-bias theorem, the carried two-moment source map \((\mathfrak g[\sigma],\mathcal S[\sigma])\), and the declared mean-preserving two-mode compensated source family with the radial attenuation closure**
\[
s(r)=\frac{r_\sigma^2}{r^2+r_\sigma^2},
\]
**used here as the minimal continuation that reproduces the sampled Session-I source point.**

This stage does **not** yet assemble the full reduced stationary barrier.
It compiles the source lane that Stage 243 only declared abstractly:

1. the exact sign-change test for the first compensated two-mode source,
2. the exact mouth-bias functional \(\mathfrak g[\sigma]\),
3. the exact shell-loading functional \(\mathcal S[\sigma]\),
4. the exact mixed-to-shell loading ratio \(\mathcal R[\sigma]\),
5. and the first explicit source packet that can be fed into Stage 247 together with the Stage-244 leakage/work lane and the Stage-245 non-rigid `U/V` packet.

---

## Purpose

Stage 243 already lifted the positive-source Family-1 corridor to the first sign-changing multimode source family,
\[
\sigma(x)=1+a\cos(\pi x)+b\cos(2\pi x),
\qquad x\in[0,1],
\]
and proved the exact sign-change test.

Stage 245 then stated the next job clearly:

> promote the mouth/core source to the first compensated multimode branch, compile its mouth-bias and shell-loading functionals, and feed those functionals into the stationary barrier compiler.

That is exactly what Stage 246 does.

The important structural point is that this stage does **not** invent a new mouth observable.
It keeps the same carried Family-1 mouth-bias functional and the same carried mixed-to-shell loading law, and only extends them from the positive corridor to the first mean-preserving sign-changing two-mode family.

---

## Provenance

This stage sits directly after:

- **Stage 243**, which declared the compensated sign-changing source lane but did not yet compile its source moments,
- **Stage 244**, which compiled the leakage/work lane on the selected-support side,
- **Stage 245**, which compiled the non-rigid `U/V` packet on the orbit side,
- the compact program’s positive-source Family-1 theorem,
- and the barrier-session Section I, which explicitly used the sign-changing source branch and recorded
  \[
  \mathfrak g[\sigma](r_{\rm eval})\approx 0.82823667,
  \qquad
  \mathcal R[\sigma](r_{\rm eval})\approx 0.21677037,
  \qquad
  \sigma_{\min}<0.
  \]

So Stage 246 is the source-side continuation of the relaxed branch.

---

## 0. Why this stage is needed

Before this step, the derivation stack had an asymmetry:

- the open-system export lane had explicit observables from Stage 244,
- the non-rigid mouth/dressing lane had explicit physical compilers from Stage 245,
- but the compensated source lane was still only a declared profile with a sign-change flag.

That was no longer enough.
The compact program already says that, at first order in source deformation, the mouth/source sector enters the Family-1 throat core only through **two** scalar source moments:
\[
\mathfrak g[\sigma]
\quad\text{and}\quad
\mathcal S[\sigma].
\]

So the next honest step is not “pick another profile by taste.”
It is to compute those two moments explicitly for the first compensated multimode family and then reconstruct the mixed-to-shell loading ratio from the carried Family-1 gain law.

That is the whole content of Stage 246.

---

## 1. Carried positive-source Family-1 source law

Work on the dimensionless mouth interval
\[
x=\frac{z}{L}\in[0,1].
\]

The compact program already fixes the two carried source kernels
\[
c(x)=\cos\!\left(\frac{\pi x}{2}\right),
\qquad
K_q(x)=\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}.
\]

The corresponding source moments are
\[
\boxed{
\mathfrak g[\sigma]=\int_0^1 \sigma(x)\,c(x)\,dx,
\qquad
\mathcal S[\sigma]=\int_0^1 \sigma(x)\,K_q(x)\,dx.
}
\]

Inside the positive-source theorem one has
\[
\sigma(x)\ge 0,
\qquad
\int_0^1 \sigma(x)\,dx=1,
\]
and therefore
\[
0\le \mathfrak g[\sigma]\le 1.
\]

The carried Family-1 mixed-to-shell loading ratio is
\[
\boxed{
\mathcal R[\sigma]
=
\frac{\bigl(\mathfrak g[\sigma]-\mathfrak r_{F1}\bigr)^2}{1+\mathfrak r_{F1}^2},
}
\]
with the fixed Family-1 geometric constant
\[
\boxed{
\mathfrak r_{F1}
=
\sqrt{\frac{12}{\pi^2}\left(\frac{37}{20}\right)^2-1}.
}
\]

So the positive-source Family-1 mouth law is
\[
\boxed{
\Pi=\Sigma_0\bigl[1-\mathcal R[\sigma]\mathcal S[\sigma]\bigr].
}
\]

Stage 246 keeps this law and changes only the source family inserted into it.

---

## 2. Mean-preserving compensated two-mode source family

Take the first compensated two-mode branch
\[
\boxed{
\sigma_{a,b}(x)=1+a\cos(\pi x)+b\cos(2\pi x),
\qquad x\in[0,1].
}
\]

Because the cosine modes integrate to zero on \([0,1]\),
\[
\boxed{
\int_0^1 \sigma_{a,b}(x)\,dx = 1.
}
\]

So this family changes source **shape** while preserving total normalized source weight exactly.

Write
\[
y=\cos(\pi x),\qquad y\in[-1,1].
\]
Then
\[
\cos(2\pi x)=2y^2-1,
\]
so the compensated source becomes the exact quadratic
\[
\boxed{
\sigma_{a,b}(y)=1-b+a y+2 b y^2.
}
\]

The candidate minima are therefore explicit:

- boundary values
  \[
  \sigma_{a,b}(1)=1+a+b,
  \qquad
  \sigma_{a,b}(-1)=1-a+b,
  \]
- interior stationary point, when \(b\neq 0\),
  \[
  y_*=-\frac{a}{4b},
  \]
  with vertex value
  \[
  \sigma_{a,b}(y_*)=1-b-\frac{a^2}{8b}.
  \]

Hence the exact minimum is
\[
\boxed{
\sigma_{\min}(a,b)=
\begin{cases}
1+b-|a|, & b\le 0\ \text{or}\ |a|\ge 4b,\\[4pt]
1-b-\dfrac{a^2}{8b}, & b>0\ \text{and}\ |a|\le 4b.
\end{cases}
}
\]

So the sign-change condition is simply
\[
\boxed{
\mathrm{signchg}(a,b)\iff \sigma_{\min}(a,b)<0.
}
\]

This is the exact extension of the Stage-243 sign-change declaration.

---

## 3. Exact mouth-bias and shell-loading compilers

The crucial point is that the two carried source moments stay exactly computable on the compensated family.

### 3.1 Mouth-bias functional

Using
\[
\mathfrak g[\sigma]=\int_0^1 \sigma(x)\cos\!\left(\frac{\pi x}{2}\right)dx,
\]
one gets
\[
\boxed{
\mathfrak g(a,b)
=
\frac{2}{\pi}\left(1+\frac{a}{3}-\frac{b}{15}\right).
}
\]

So the same mouth-bias functional that was carried through the positive-source theorem extends linearly to the compensated two-mode family.

### 3.2 Shell-loading functional

Using
\[
\mathcal S[\sigma]=\int_0^1 \sigma(x)\,K_q(x)\,dx,
\qquad
K_q(x)=\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)},
\]
one gets
\[
\boxed{
\mathcal S(a,b)
=
\frac{2\tanh(\pi/2)}{\pi}
\left(1+\frac{a}{5}+\frac{b}{17}\right).
}
\]

So the compensated two-mode branch is exactly the first source family that can tune the two carried source moments away from the positive exponential Family-1 slice while still remaining fully explicit.

---

## 4. Exact two-moment source map and its inverse

Define the normalized source moments
\[
\widetilde g:=\frac{\pi}{2}\mathfrak g-1,
\qquad
\widetilde S:=\frac{\pi}{2\tanh(\pi/2)}\mathcal S-1.
\]

Then
\[
\boxed{
\begin{pmatrix}
\widetilde g\\[4pt]
\widetilde S
\end{pmatrix}
=
\begin{pmatrix}
\dfrac13 & -\dfrac1{15}\\[8pt]
\dfrac15 & \dfrac1{17}
\end{pmatrix}
\begin{pmatrix}
a\\[4pt]
b
\end{pmatrix}.
}
\]

The determinant is
\[
\boxed{
\det M_{\rm src}=\frac{14}{425}>0.
}
\]

So the first compensated two-mode branch is **exactly invertible** as a source compiler.
Its inverse map is
\[
\boxed{
a=\frac{85}{42}\widetilde S+\frac{25}{14}\widetilde g,
\qquad
b=\frac{425}{42}\widetilde S-\frac{85}{14}\widetilde g.
}
\]

This is structurally important:

> the two-mode branch is not just “some sign-changing profile.”
> It is the first exact two-parameter source family that independently spans the two carried Family-1 source moments \((\mathfrak g,\mathcal S)\).

That is why it is the correct Stage-246 source compiler.

---

## 5. Exact mixed-to-shell loading ratio on the compensated branch

The carried Family-1 loading law gives
\[
\boxed{
\mathcal R(a,b)
=
\frac{\bigl(\mathfrak g(a,b)-\mathfrak r_{F1}\bigr)^2}{1+\mathfrak r_{F1}^2}.
}
\]

So after Stage 246 the compensated source branch enters the reduced mouth law only through the explicit source packet
\[
\boxed{
\mathcal M_\sigma(a,b)
=
\bigl(
\sigma_{a,b},\,
\sigma_{\min},\,
\mathrm{signchg},\,
\mathfrak g(a,b),\,
\mathcal S(a,b),\,
\mathcal R(a,b)
\bigr).
}
\]

### 5.1 Compensation line

For any target mouth-bias level \(\mathfrak g_c\), the corresponding exact line in \((a,b)\)-space is
\[
\boxed{
b
=
5a+15-\frac{15\pi}{2}\,\mathfrak g_c.
}
\]

In particular, on the carried Family-1 lower compensated branch
\[
\mathfrak g_-(\mathfrak r)
=
\mathfrak r-\frac12\sqrt{1+\mathfrak r^2},
\]
so
\[
\boxed{
b
=
5a+15-\frac{15\pi}{2}\,\mathfrak g_-\!\bigl(\mathfrak r_{F1}\bigr).
}
\]

### 5.2 Exact quarter-ratio test

If
\[
\mathfrak g(a,b)=\mathfrak g_\pm(\mathfrak r_{F1})
=
\mathfrak r_{F1}\pm\frac12\sqrt{1+\mathfrak r_{F1}^2},
\]
then
\[
\boxed{
\mathcal R(a,b)=\frac14.
}
\]

So the quarter-ratio theorem survives exactly on the compensated two-mode branch.
The positive-source theorem had forced the lower branch as the only admissible positive candidate.
The compensated family keeps the same quarter-ratio surface but no longer enforces positivity.

---

## 6. Radially transported source branch used for the Session-I readback

To connect the exact two-mode source family to the recorded Session-I sampled point, use the minimal radial attenuation closure
\[
\boxed{
a(r)=a_0\,s(r),
\qquad
b(r)=b_0\,s(r),
\qquad
s(r)=\frac{r_\sigma^2}{r^2+r_\sigma^2}.
}
\]

Then the full transported source is
\[
\boxed{
\sigma(x;r)=1+a(r)\cos(\pi x)+b(r)\cos(2\pi x).
}
\]

The two carried source moments become
\[
\boxed{
\mathfrak g(r)
=
\frac{2}{\pi}
\left[
1+s(r)\left(\frac{a_0}{3}-\frac{b_0}{15}\right)
\right],
}
\]
\[
\boxed{
\mathcal S(r)
=
\frac{2\tanh(\pi/2)}{\pi}
\left[
1+s(r)\left(\frac{a_0}{5}+\frac{b_0}{17}\right)
\right].
}
\]

The mixed-to-shell ratio becomes
\[
\boxed{
\mathcal R(r)
=
\frac{\bigl(\mathfrak g(r)-\mathfrak r_{F1}\bigr)^2}{1+\mathfrak r_{F1}^2}.
}
\]

### 6.1 Exact sign-change threshold for the Session-I orientation

On the Session-I orientation
\[
a_0>0,\qquad b_0<0,
\]
the branch stays in the boundary-minimum case, so
\[
\boxed{
\sigma_{\min}(r)=1+b(r)-a(r)=1-(a_0-b_0)s(r).
}
\]

Therefore the exact sign-change threshold is
\[
\boxed{
\mathrm{signchg}(r)\iff
s(r)>\frac{1}{a_0-b_0}
\iff
r<r_\sigma\sqrt{a_0-b_0-1}.
}
\]

So the transported compensated source has a sharp same-charge contact window in which the source profile necessarily becomes sign-changing.

---

## 7. Exact Session-I readback

Use the recorded Session-I source-side parameters
\[
a_0=2.2,
\qquad
b_0=-0.6,
\qquad
r_\sigma=0.8,
\qquad
r_{\rm eval}=1.00217028,
\]
together with the carried Family-1 constant
\[
\mathfrak r_{F1}\approx 1.77799353547498.
\]

Then
\[
s(r_{\rm eval})
=
\frac{r_\sigma^2}{r_{\rm eval}^2+r_\sigma^2}
\approx 0.38921266210419.
\]

So the transported mode amplitudes are
\[
a(r_{\rm eval})\approx 0.856267856629217,
\qquad
b(r_{\rm eval})\approx -0.233527597262514.
\]

The exact Stage-246 source compilers then give
\[
\boxed{
\mathfrak g[\sigma](r_{\rm eval})
\approx 0.828236674079292,
}
\]
\[
\boxed{
\mathcal S[\sigma](r_{\rm eval})
\approx 0.675847711465632,
}
\]
\[
\boxed{
\mathcal R[\sigma](r_{\rm eval})
\approx 0.216770371559385,
}
\]
\[
\boxed{
\sigma_{\min}(r_{\rm eval})
\approx -0.089795453891731,
\qquad
\mathrm{signchg}=\mathrm{True}.
}
\]

The last three values reproduce the sampled Session-I mouth/source diagnostics exactly at the precision reported in the session table.

Two structural consequences follow immediately.

1. The transported Session-I source branch is genuinely **outside** the positive-source corridor because \(\sigma_{\min}<0\).
2. The same branch still lives inside the **same carried Family-1 source law**: nothing new had to be invented beyond the compensated two-mode source compiler.

A further diagnostic is
\[
\mathfrak g(0)
=
\frac{2}{\pi}\left(1+\frac{a_0}{3}-\frac{b_0}{15}\right)
\approx 1.12893906>1,
\]
so the transported branch really does leave the positive-source theorem window at sufficiently small separation.

---

## 8. What this stage achieves physically

Stage 246 closes the source-side gap left open by Stages 243–245.

1. It keeps the exact carried Family-1 source observables \(\mathfrak g[\sigma]\) and \(\mathcal S[\sigma]\) rather than replacing them with new phenomenological knobs.
2. It proves that the first sign-changing two-mode source family is exactly invertible onto those two source moments.
3. It reconstructs the mixed-to-shell loading ratio \(\mathcal R[\sigma]\) from the same carried gain law that the compact program already uses.
4. It reproduces the sampled Session-I mouth/source diagnostics exactly with one explicit transported compensated source branch.

So after Stage 246 the compensated source lane is no longer just a profile with a sign-change flag.
It is an exact source packet that can be inserted directly into the stationary barrier compiler.

---

## 9. Result

The compensated multimode source lane now has an explicit derivation-stack compiler.

The exact packet is
\[
\boxed{
\mathcal M_\sigma(r)
=
\bigl(
\sigma(x;r),\,
\sigma_{\min}(r),\,
\mathrm{signchg}(r),\,
\mathfrak g[\sigma](r),\,
\mathcal S[\sigma](r),\,
\mathcal R[\sigma](r)
\bigr),
}
\]
with
\[
\mathfrak g[\sigma](r)
=
\frac{2}{\pi}
\left[
1+\frac{a(r)}{3}-\frac{b(r)}{15}
\right],
\qquad
\mathcal S[\sigma](r)
=
\frac{2\tanh(\pi/2)}{\pi}
\left[
1+\frac{a(r)}{5}+\frac{b(r)}{17}
\right],
\]
\[
\mathcal R[\sigma](r)
=
\frac{\bigl(\mathfrak g[\sigma](r)-\mathfrak r_{F1}\bigr)^2}{1+\mathfrak r_{F1}^2}.
\]

This is exactly the source-side packet that Stage 247 needs.

---

## 10. Immediate next step

The next stage should now assemble the first honest reduced stationary barrier compiler from the three relaxed lanes:

1. Stage-244 leakage/work packet,
2. Stage-245 non-rigid `U/V` packet,
3. Stage-246 compensated source packet.

At that point the lowered barrier branch can be written as one explicit reduced stationary compiler rather than as three separate session-side ingredients.
