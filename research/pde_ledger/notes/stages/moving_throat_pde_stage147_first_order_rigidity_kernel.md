
# Moving-Throat PDE — Stage 249: First-Order Rigidity Kernel at the Canonical Family-1 Point

## Goal

Turn the general first-order correction formulas from Stage 248 into one explicit
rigidity law for the physical normalized mouth traction
\[
\widehat T_m.
\]

---

## 1. Canonical branch traction law

On the lower compensated branch one has
\[
R_q=\frac14,
\qquad
\Sigma_0=\frac{\Pi}{1-\mathcal S_q/4},
\qquad
\widehat T_m=\sqrt{\frac{9\Sigma_0}{20}}.
\]

So at fixed canonical branch structure,
\[
\delta \Sigma_0
=
\frac{\delta\Pi}{1-\mathcal S_*/4}
+
\frac{\Pi_*}{4(1-\mathcal S_*/4)^2}\,\delta\mathcal S_q.
\]

Substituting the Stage 248 retuning formulas gives
\[
\boxed{
\delta \widehat T_m
=
\epsilon\Big[
A_T\,(\bar g_\varsigma-\mathfrak g_*)
+
B_T\,(\bar S_\varsigma-\mathcal S_*)
\Big],
}
\]
with
\[
A_T
=
-\frac{9}{40\widehat T_{m,*}}
\left[
\frac{1}{\mathfrak g'_*(1-\mathcal S_*/4)}
+
\frac{\Pi_*\,\mathcal S'_*}{4\mathfrak g'_*(1-\mathcal S_*/4)^2}
\right],
\]
\[
B_T
=
\frac{9}{40\widehat T_{m,*}}
\frac{\Pi_*}{4(1-\mathcal S_*/4)^2}.
\]

Numerically,
\[
\boxed{
A_T \approx -4.27263956256927,
\qquad
B_T \approx 0.134875005736706.
}
\]

So the canonical branch is vastly more sensitive to overlap changes than to
mixed-kernel changes:
\[
\boxed{
\frac{|A_T|}{B_T}\approx 31.6785.
}
\]

---

## 2. One effective rigidity kernel

Because \(\Sigma_\epsilon-\Sigma_*\) integrates to zero, the traction shift can be written
as a single weighted overlap:
\[
\boxed{
\delta \widehat T_m
=
\epsilon
\int_0^1
\mathcal W_*(x)\,
\bigl[\varsigma(x)-\Sigma_*(x)\bigr]\,dx,
}
\]
with centered weight
\[
\boxed{
\mathcal W_*(x)
=
A_T\Big(c(x)-\mathfrak g_*\Big)
+
B_T\Big(K_q(x)-\mathcal S_*\Big).
}
\]

So after retuning the electrochemical bias, the entire first-order non-exponential
shape sensitivity collapses to **one scalar kernel** \(\mathcal W_*(x)\).

This is the right rigidity statement:

> once the canonical lower branch is fixed, a positive mouth-layer deformation can move the
> Family-1 point only through its overlap with one known rigidity kernel.

---

## 3. Meaning

The branch is not infinitely rigid. But it is much more rigid than the earlier
branch-choice ambiguity suggested:

- the deformation space collapses to two source moments,
- the canonical retuning removes one of them as an independent datum,
- and the remaining first-order traction shift is controlled by one explicit kernel.

That is a strong reduction of the mouth-side uncertainty.
