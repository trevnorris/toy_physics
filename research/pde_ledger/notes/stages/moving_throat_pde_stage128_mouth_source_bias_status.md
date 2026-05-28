
# Moving-Throat PDE — Stage 128: Mouth-Source Bias Status

## What is now fixed

The outlet-core branch-selection problem has narrowed further.

1. On the explicit Family-1 throat geometry,
\[
\mathfrak r_{F1}\approx 1.77799353547498
\]
is already fixed by geometry.

2. For any **positive normalized localized mouth source** on the first D/N interval,
\[
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)\cos\!\left(\frac{\pi z}{2L}\right)\,dz
\]
satisfies
\[
0\le \mathfrak g\le 1.
\]

3. Therefore the upper compensated branch
\[
\mathfrak g_+^{F1}\approx 2.79795199200529
\]
is ruled out by positivity alone.

4. The lower compensated branch
\[
\mathfrak g_-^{F1}\approx 0.758035078944663
\]
is the unique physically admissible canonical branch under positive mouth sourcing.

## What the explicit source laws show

- Point source:
  \[
  \mathfrak g=1
  \]
  (the old naive equal-normalized branch).

- Self-matched D/N derivative source:
  \[
  \mathfrak g=\frac{\pi}{4}\approx 0.785398
  \]
  which is already only \(3.61\%\) away in traction from exact lower-branch compensation.

- Exact convex positive family:
  \[
  \sigma_\xi(z)=(1-\xi)\,k\cos(kz)+\xi/L
  \]
  hits the exact lower compensated branch at
  \[
  \xi_*\approx 0.183918.
  \]

- Uniform slab and truncated exponential penetration families also hit the exact lower branch at moderate penetration depths:
  \[
  x_*^{\rm slab}\approx 0.797839,
  \qquad
  x_*^{\exp}\approx 0.662765.
  \]

## Meaning

The branch-selection ambiguity is no longer “does some unknown source law maybe rescue the canonical outlet?”

It is now much sharper:

\[
\boxed{
\text{Under positive localized mouth sourcing, the lower compensated branch is uniquely selected and is easily reachable.}
}
\]

So the remaining PDE-facing question is not branch sign. It is the detailed shape of the actual mouth source profile \(\sigma(z)\), or equivalently the exact amount of positive mouth broadening away from the point-source limit.
