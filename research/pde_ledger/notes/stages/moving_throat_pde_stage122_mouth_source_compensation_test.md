# Moving-Throat PDE — Stage 122: Concrete Mouth-Source Branch and Compensation Test

## Goal

Compute the normalized mouth-coupling ratio
\[
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}
\]
on the first explicit local mouth-source branch, and test whether that branch already lands on the canonical compensation family.

## 1. Natural equal-normalized local mouth source

Take one local mouth source profile that couples to the shell and the mixed channel with equal strength **after normalization by their quadratic channel weights**:
\[
\boxed{
g_s=g_m\sqrt{K_s},
\qquad
g_q=g_m\sqrt{K_q}.
}
\]
This is the simplest concrete local branch with **no channel favoritism** at the mouth.

Then
\[
\boxed{
\mathfrak g_{\rm nat}=1.
}
\]

This should be read as a concrete branch ansatz, not as a theorem of the full PDE: it is the cleanest “same mouth source, same normalized loading” core closure.

## 2. Exact compensation family

Stage 119 already showed that the compensated canonical outgoing branch requires
\[
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2,
\]
equivalently
\[
\boxed{
\mathfrak g_{\rm comp}^{\pm}(\mathfrak r)
=
\mathfrak r\pm\frac12\sqrt{1+\mathfrak r^2}.
}
\]

Inserting the Stage 121 geometric value
\[
\mathfrak r=\mathfrak r_{F1}=\frac{\sqrt{4107-100\pi^2}}{10\pi}
\]
gives the two exact compensated mouth-coupling values
\[
\boxed{
\mathfrak g_{\pm}^{F1}
=
\frac{2\sqrt{4107-100\pi^2}\pm 37\sqrt3}{20\pi}.
}
\]

Numerically,
\[
\boxed{
\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathfrak g_+^{F1}\approx 2.79795199200529.
}
\]

## 3. Test of the natural branch

The natural equal-normalized mouth branch gives
\[
\mathfrak g_{\rm nat}=1.
\]
That is **not** exactly on the compensation family at \(\mathfrak r=\mathfrak r_{F1}\), since
\[
1\neq \mathfrak g_-^{F1},
\qquad
1\neq \mathfrak g_+^{F1}.
\]

Equivalently, the exact compensation defect evaluated on the natural branch is
\[
\mathcal C_{\rm nat}
:=
1+\mathfrak r_{F1}^2-4(1-\mathfrak r_{F1})^2
=
\frac{-12321+80\pi\sqrt{4107-100\pi^2}}{100\pi^2}
\approx 1.74016524722739,
\]
so the natural branch misses the compensated surface.

## 4. How far away is the nearest compensated branch?

The **near** compensated branch is the minus branch, since
\[
|\mathfrak g_{\rm nat}-\mathfrak g_-^{F1}|
<
|\mathfrak g_{\rm nat}-\mathfrak g_+^{F1}|.
\]
Numerically,
\[
\boxed{
\Delta \mathfrak g_- = 1-\mathfrak g_-^{F1}\approx 0.241964921055337,
}
\]
while
\[
\boxed{
\Delta \mathfrak g_+ = \mathfrak g_+^{F1}-1\approx 1.79795199200529.
}
\]

So the natural equal-normalized mouth source is **not** automatically compensated, but it is much closer to the lower compensated branch than to the upper one.

## 5. Exact traction renormalization factors

Since Stage 119 gave
\[
\mathfrak g \propto \frac{1}{\mathcal T_m},
\]
the traction amplitude required to move from the natural branch \((\mathfrak g=1)\) to the compensated branches is
\[
\boxed{
\frac{\mathcal T_m^{(\pm)}}{\mathcal T_m^{\rm nat}}
=
\frac{1}{\mathfrak g_\pm^{F1}}.
}
\]

Numerically,
\[
\boxed{
\frac{\mathcal T_m^{(-)}}{\mathcal T_m^{\rm nat}}
\approx 1.31920016339112,
\qquad
\frac{\mathcal T_m^{(+)}}{\mathcal T_m^{\rm nat}}
\approx 0.357404273860789.
}
\]

So the nearest compensated branch is reached by a **31.9% traction enhancement** relative to the simplest equal-normalized mouth-source branch.

## Result

The first truly concrete mouth-core test is now complete:

- the actual throat geometry fixes
  \[
  \mathfrak r_{F1}\approx 1.778;
  \]
- the simplest equal-normalized local mouth source gives
  \[
  \mathfrak g_{\rm nat}=1;
  \]
- this does **not** lie exactly on the compensation surface;
- but the nearest compensated branch is only a modest traction renormalization away.

So the surviving ambiguity is no longer broad outlet algebra. It is whether the real moving-throat mouth source is exactly equal-normalized, or instead slightly biased toward the lower compensated branch.
