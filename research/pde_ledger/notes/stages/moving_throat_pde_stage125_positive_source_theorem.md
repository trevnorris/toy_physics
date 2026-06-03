
# Moving-Throat PDE — Stage 125: Positive Local Mouth-Source Theorem

## Goal

Replace the still-free normalized mouth-coupling ratio
\[
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}
\]
by an explicit **localized mouth-source law** and determine which compensation branch is physically admissible.

## Localized source law

Take one scalar mouth amplitude \(u\) and one nonnegative normalized axial source profile
\[
\sigma(z)\ge 0,
\qquad
\int_0^L \sigma(z)\,dz = 1,
\qquad
z\in[0,L].
\]
Assume the shell mouth load uses the same positive source density throughout the mouth layer, while the mixed channel couples through the first D/N half-wave derivative on the same throat interval,
\[
\chi_{1/2}(z)=\sqrt{\frac{2}{L}}\sin\!\left(\frac{\pi z}{2L}\right),
\qquad
\chi_{1/2}'(z)=\sqrt{\frac{2}{L}}\frac{\pi}{2L}\cos\!\left(\frac{\pi z}{2L}\right).
\]
This is the first one-lane positive localized-source closure for the mouth problem.
It is meant to capture positive localized source data on the first D/N throat
interval, not arbitrary sign-changing, multimode, or nonlocalized mouth data.

After normalizing to the point-source branch \(\sigma(z)=\delta(z)\), the mouth-bias factor becomes
\[
\boxed{
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)\cos\!\left(\frac{\pi z}{2L}\right)\,dz.
}
\]

So the mouth-source bias is the first cosine moment of the positive axial source profile.

## Exact positivity theorem

Because
\[
0\le \cos\!\left(\frac{\pi z}{2L}\right)\le 1
\qquad \text{for } z\in[0,L],
\]
every positive normalized source law satisfies
\[
\boxed{
0\le \mathfrak g[\sigma]\le 1.
}
\]

This immediately selects the physically admissible compensation branch.

On the explicit Family-1 branch, the carried branch values are
\[
\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathfrak g_+^{F1}\approx 2.79795199200529.
\]
Therefore
\[
\boxed{
\mathfrak g_+^{F1}>1
\quad\Rightarrow\quad
\text{the upper compensated branch is impossible for any positive localized source law.}
}
\]

By contrast,
\[
\boxed{
0<\mathfrak g_-^{F1}<1
\quad\Rightarrow\quad
\text{within this positive localized-source setup, the lower compensated branch is the unique physically admissible compensated Family-1 branch.}
}
\]

## Interpretation

This is a strong narrowing.

The outlet-core ambiguity is no longer “which branch do we choose?” Within a localized positive source model on the first D/N throat interval:

- the upper compensated branch is ruled out,
- the lower compensated branch is the only canonical candidate.

So the remaining question is not branch sign, but the **shape** of the positive source profile \(\sigma(z)\).
