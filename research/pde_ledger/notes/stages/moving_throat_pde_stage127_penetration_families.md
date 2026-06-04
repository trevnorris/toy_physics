
# Moving-Throat PDE — Stage 127: Geometric Mouth-Penetration Families

## Goal

Translate the lower compensated mouth-source bias into simple geometric penetration scales.

## 1. Uniform slab source

Take the positive localized source profile
\[
\boxed{
\sigma_x^{\rm slab}(z)=\frac{1}{xL},
\qquad
0\le z\le xL,
\qquad
0<x\le 1,
}
\]
and zero elsewhere on \([0,L]\).

Then
\[
\mathfrak g_{\rm slab}(x)
=
\int_0^{xL}\frac{dz}{xL}\cos\!\left(\frac{\pi z}{2L}\right)
=
\boxed{
\frac{2}{\pi x}\sin\!\left(\frac{\pi x}{2}\right).
}
\]

Solving
\[
\mathfrak g_{\rm slab}(x)=\mathfrak g_-^{F1}
\]
gives the unique Family-1 slab depth
\[
\boxed{
x_*^{\rm slab}\approx 0.797839360904564.
}
\]

So a positive uniform mouth source extending over about \(80\%\) of the throat span reaches the lower compensated branch exactly.

## 2. Truncated exponential source

Take the normalized positive exponential family
\[
\boxed{
\sigma_x^{\exp}(z)
=
\frac{e^{-z/(xL)}}{xL\left(1-e^{-1/x}\right)},
\qquad
z\in[0,L],
\qquad
x>0.
}
\]

Its mouth-bias factor is
\[
\boxed{
\mathfrak g_{\exp}(x)
=
\frac{2\left(2+\pi x\,e^{-1/x}\right)}
{\left(4+\pi^2x^2\right)\left(1-e^{-1/x}\right)}.
}
\]

Solving
\[
\mathfrak g_{\exp}(x)=\mathfrak g_-^{F1}
\]
gives
\[
\boxed{
x_*^{\exp}\approx 0.662765402623160.
}
\]

So an exponentially localized mouth source with penetration depth about \(0.66L\) reaches the exact lower compensated branch.

## Result

These two geometric positive-source families both hit the same lower compensated branch with **moderate** mouth penetration depths:
\[
x_*^{\exp}\approx 0.66,
\qquad
x_*^{\rm slab}\approx 0.80.
\]

Combined with Stages 125–126, this means:

- the upper compensated branch is excluded by positivity,
- the lower compensated branch is the unique physical candidate,
- and simple positive localized source laws reach it without requiring sign-changing or finely oscillatory mouth forcing.
