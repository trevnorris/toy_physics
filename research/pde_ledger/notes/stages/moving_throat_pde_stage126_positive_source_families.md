
# Moving-Throat PDE — Stage 126: Explicit Positive Source Families and the Family-1 Compensation Point

## Goal

Evaluate explicit positive localized source laws and determine whether the physical lower compensated branch can be reached without exotic sign-changing mouth profiles.

## 1. Self-matched derivative profile

The most natural positive axial source on the first D/N interval is the normalized derivative profile itself,
\[
\boxed{
\sigma_{\rm match}(z)=k\cos(kz),
\qquad
k=\frac{\pi}{2L},
}
\]
since
\[
\int_0^L k\cos(kz)\,dz = 1.
\]

Its mouth-bias factor is exact:
\[
\mathfrak g_{\rm match}
=
\int_0^L \sigma_{\rm match}(z)\cos(kz)\,dz
=
k\int_0^L \cos^2(kz)\,dz
=
\frac{\pi}{4}.
\]
So
\[
\boxed{
\mathfrak g_{\rm match}=\frac{\pi}{4}\approx 0.785398163397448.
}
\]

This already lies much closer to the physical lower compensated branch than the naive point-source value \(\mathfrak g=1\).

## 2. Family-1 comparison

On the explicit geometric branch from Stages 121–122,
\[
\mathfrak g_-^{F1}\approx 0.758035078944663.
\]
Therefore the self-matched source overshoots the exact compensated lower branch by
\[
\Delta\mathfrak g_{\rm match}
=
\frac{\pi}{4}-\mathfrak g_-^{F1}
\approx 0.0273630844527852.
\]

Since traction scales as \(\mathfrak g^{-1}\),
\[
\boxed{
\frac{\mathcal T_m^{(-)}}{\mathcal T_m^{\rm match}}
=
\frac{\pi/4}{\mathfrak g_-^{F1}}
\approx 1.036097385480999.
}
\]
So the exact Family-1 compensated branch is only a **3.61% traction enhancement** away from the self-matched derivative profile.

## 3. Exact positive family interpolating through the compensated point

Introduce the convex positive family
\[
\boxed{
\sigma_\xi(z)
=
(1-\xi)\,k\cos(kz)+\xi\,\frac1L,
\qquad 0\le \xi\le 1.
}
\]
This is normalized and nonnegative on \([0,L]\). Its mouth-bias factor is
\[
\mathfrak g_\xi
=
(1-\xi)\frac{\pi}{4}+\xi\frac{2}{\pi}.
\]

Because
\[
\frac{2}{\pi}<\mathfrak g_-^{F1}<\frac{\pi}{4},
\]
there is a unique \(\xi_*\in(0,1)\) such that
\[
\mathfrak g_{\xi_*}=\mathfrak g_-^{F1}.
\]
The exact solution is
\[
\boxed{
\xi_*
=
\frac{\frac{\pi}{4}-\mathfrak g_-^{F1}}{\frac{\pi}{4}-\frac{2}{\pi}}
=
\frac{-37\sqrt3-5\pi^2+2\sqrt{4107-100\pi^2}}{5(8-\pi^2)}
\approx 0.183918405511538.
}
\]

So only an **18.4% admixture** of the fully washed positive profile \(1/L\) into the self-matched derivative profile reaches the exact lower compensated branch.

## Result

The explicit mouth-source bias is now much narrower than before:

- the point-source branch \(\mathfrak g=1\) is too high,
- the exact self-matched D/N derivative source gives \(\mathfrak g=\pi/4\),
- and the true Family-1 compensated branch is reached by a small \(18.4\%\) broadening of that already-natural positive source law.

So the canonical branch no longer looks like a delicate coefficient fit. It sits inside a simple exact family of positive localized mouth sources.
