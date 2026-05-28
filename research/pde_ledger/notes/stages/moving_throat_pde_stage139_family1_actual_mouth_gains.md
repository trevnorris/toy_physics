# Moving-Throat PDE — Stage 139: Actual Family-1 Mouth Gains

## Goal

Evaluate the explicit gain formulas of Stages 137–138 on the concrete Family-1 branch.

This answers the practical question:

> once the geometric core branch is fixed, what mouth gains does the actual throat-core
> ask for at the canonical compensation point `\Pi_*`?

---

## 1. Family-1 geometric input

From Stage 121,
\[
\mathfrak r_{F1}
=
\sqrt{\frac{12}{\pi^2}\left(\frac{37}{20}\right)^2-1}
\approx 1.77799353547498.
\]

From Stage 134,
\[
\Pi_* \approx 1.50882951349316,
\qquad
\mathcal S_q(\Pi_*)\approx 0.658075937605429.
\]

So the actual Family-1 fixed-point law is
\[
\Pi_* = M_s + M_q\,\mathcal S_q(\Pi_*).
\]
Using Stage 138,
\[
M_q=-R_q M_s,
\qquad
R_q=\frac{(\mathfrak g_c-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
\]
so the exact shell gain required at `\Pi_*` is
\[
\boxed{
M_s^*(R_q)=\frac{\Pi_*}{1-R_q\mathcal S_q(\Pi_*)},
\qquad
M_q^*(R_q)=-R_q M_s^*(R_q).
}
\]

---

## 2. Natural equal-normalized mouth-source branch

On the simplest natural branch,
\[
\mathfrak g_c=1.
\]
Then
\[
\boxed{
R_q^{\rm nat}
=
\frac{(1-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}
\approx 0.145454452260421.
}
\]

So the actual Family-1 mouth gains selected at the canonical compensation point are
\[
\boxed{
M_s^{\rm nat,*}
=
\frac{\Pi_*}{1-R_q^{\rm nat}\mathcal S_q(\Pi_*)}
\approx 1.66854252965624,
}
\]
\[
\boxed{
M_q^{\rm nat,*}
=-R_q^{\rm nat}M_s^{\rm nat,*}
\approx -0.242696939724365.
}
\]

So the natural equal-normalized core source law does **not** exactly reproduce the
outlet-consistent canonical ratio, but it is already on the correct sign branch and not
far away from it.

---

## 3. Exact compensated branch

On the exact compensated family,
\[
\mathfrak g_c
=
\mathfrak r_{F1}-\frac12\sqrt{1+\mathfrak r_{F1}^2}
\approx 0.758035078944663,
\qquad
R_q=\frac14.
\]

Then the actual gains are
\[
\boxed{
M_s^{\rm comp,*}
=
\frac{\Pi_*}{1-\mathcal S_q(\Pi_*)/4}

\approx 1.80594111095636,
}
\]
\[
\boxed{
M_q^{\rm comp,*}=-\frac{M_s^{\rm comp,*}}{4}
\approx -0.451485277739090.
}
\]

This is exactly the Stage 135 one-parameter canonical branch.

---

## 4. Quantitative comparison

The natural and canonical gains differ only moderately:
\[
\frac{M_s^{\rm comp,*}}{M_s^{\rm nat,*}}-1
\approx 0.0823464663669,
\qquad
\frac{|M_q^{\rm comp,*}|}{|M_q^{\rm nat,*}|}
\approx 1.86028418097.
\]

So the shell gain changes by about `8.23%`, while the mixed gain magnitude must increase
by about a factor `1.86` to land exactly on the canonical compensated branch.

This is consistent with the earlier source-family result that the Family-1 branch was already
much closer to the lower compensated branch than to the forbidden upper one.

---

## Result

The actual Family-1 mouth gains are now explicit.

- natural equal-normalized core source:
  \[
  M_s\approx 1.66854,
  \qquad
  M_q\approx -0.24270;
  \]
- exact compensated canonical branch:
  \[
  M_s\approx 1.80594,
  \qquad
  M_q\approx -0.45149.
  \]

So the remaining ambiguity is no longer “what are the gains?” It is only whether the real
mouth core stays on the natural branch or shifts modestly toward the lower compensated one.
