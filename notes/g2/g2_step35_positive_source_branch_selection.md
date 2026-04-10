# Step 35 — Positive mouth sourcing selects the lower branch without fine tuning

## Goal

Step 34 reduced the compensated outlet to a one-parameter parent family labeled by
`\mathfrak r`, with the required mouth-coupling ratio fixed by

```math
\mathfrak g=\mathfrak r\pm\frac12\sqrt{1+\mathfrak r^2}.
```

The next real question is no longer algebraic. It is physical:

> does any natural localized mouth-source law actually land on that family, and if
> so, does it do it without fine tuning?

This step answers that question on the explicit Family-1 branch.

---

## Step 35A — The explicit Family-1 parent ratio and its two compensated branches

For the carried explicit Family-1 geometry,

```math
\boxed{
\mathfrak r_{F1}
=
\sqrt{\frac{12}{\pi^2}\left(\frac{37}{20}\right)^2-1}
\approx 1.77799353547498.
}
```

So the two compensated parent branches are

```math
\boxed{
\mathfrak g_-^{F1}
=
\mathfrak r_{F1}-\frac12\sqrt{1+\mathfrak r_{F1}^2}
\approx 0.758035078944663,
}
```

```math
\boxed{
\mathfrak g_+^{F1}
=
\mathfrak r_{F1}+\frac12\sqrt{1+\mathfrak r_{F1}^2}
\approx 2.79795199200529.
}
```

So the purely algebraic family still contains two branches. The physical mouth
source law has to decide which one is actually admissible.

---

## Step 35B — Positive mouth-source theorem

Let `\sigma(z)` be any nonnegative normalized axial source profile on the first
D/N throat interval,

```math
\sigma(z)\ge 0,
\qquad
\int_0^L \sigma(z)\,dz = 1,
\qquad
z\in[0,L].
```

The normalized mouth-bias factor is then

```math
\boxed{
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)
\cos\!\left(\frac{\pi z}{2L}\right)
\,dz.
}
```

Because

```math
0\le \cos\!\left(\frac{\pi z}{2L}\right)\le 1
\qquad \text{on } [0,L],
```

any positive normalized source law satisfies

```math
\boxed{
0\le \mathfrak g[\sigma]\le 1.
}
```

That immediately kills the upper Family-1 branch, because

```math
\mathfrak g_+^{F1}>1.
```

So the upper compensated branch is physically impossible for any positive localized
mouth source.

By contrast,

```math
0<\mathfrak g_-^{F1}<1,
```

so the lower compensated branch is the **unique physically admissible canonical
candidate** under positive mouth sourcing.

---

## Step 35C — The natural self-matched derivative profile is already very close

The most natural positive axial source on the first D/N interval is the normalized
first-derivative profile itself,

```math
\boxed{
\sigma_{\rm match}(z)=k\cos(kz),
\qquad
k=\frac{\pi}{2L}.
}
```

It is normalized, and its mouth-bias factor is exactly

```math
\boxed{
\mathfrak g_{\rm match}=\frac{\pi}{4}\approx 0.785398163397448.
}
```

Comparing this to the exact Family-1 lower branch,

```math
\mathfrak g_-^{F1}\approx 0.758035078944663,
```

gives only a small traction mismatch. Since traction scales as
`\mathfrak g^{-1}`,

```math
\boxed{
\frac{\mathcal T_m^{(-)}}{\mathcal T_m^{\rm match}}
=
\frac{\pi/4}{\mathfrak g_-^{F1}}
\approx 1.036097385480999.
}
```

So the exact lower compensated Family-1 branch is only a **3.61% traction
enhancement** away from the most natural self-matched derivative source law.

That is already a strong “not fine tuned” signal.

---

## Step 35D — An explicit convex positive family hits the lower branch with only an 18.4% admixture

Introduce the exact convex positive family

```math
\boxed{
\sigma_\xi(z)
=
(1-\xi)\,k\cos(kz)+\xi\,\frac1L,
\qquad 0\le \xi\le 1.
}
```

It stays positive and normalized, and its bias factor is

```math
\boxed{
\mathfrak g_\xi
=
(1-\xi)\frac{\pi}{4}
+
\xi\frac{2}{\pi}.
}
```

Solving `\mathfrak g_\xi=\mathfrak g_-^{F1}` gives the exact interpolation point

```math
\boxed{
\xi_*
=
\frac{\frac{\pi}{4}-\mathfrak g_-^{F1}}{\frac{\pi}{4}-\frac{2}{\pi}}
\approx 0.183918405511538.
}
```

So only an **18.4% admixture** of the fully washed positive profile `1/L` into the
natural derivative profile reaches the exact lower compensated branch.

Again, that is a small deformation of a very natural source family, not a delicate
coefficient fit.

---

## Step 35E — The regular canonical point occurs at finite bias

For the explicit exponential mouth-layer family, the bias law is

```math
\boxed{
\mathfrak g_\Pi
=
\frac{2\Pi(2\Pi e^\Pi+\pi)}{(4\Pi^2+\pi^2)(e^\Pi-1)}.
}
```

A direct rearrangement gives an exact positivity decomposition for `1-\mathfrak g_\Pi`,
so in particular

```math
0<\mathfrak g_\Pi<1
\qquad \text{for every finite } \Pi>0,
```

and

```math
\lim_{\Pi\to\infty}\mathfrak g_\Pi = 1.
```

So the naive equal-normalized branch `\mathfrak g=1` is only reached in the
**singular point-source limit**, not at any finite regular mouth bias.

By contrast, the lower compensated Family-1 point occurs already at the finite bias

```math
\boxed{
\Pi_*\approx 1.50882951349316,
\qquad
\mathfrak g_\Pi(\Pi_*)=\mathfrak g_-^{F1}.
}
```

So the physically selected lower branch is regular and nearby, not a singular limit.

---

## Main result of the step

Positive mouth sourcing makes the branch selection problem almost trivial.

For the explicit Family-1 branch,

```math
\boxed{
\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathfrak g_+^{F1}\approx 2.79795199200529.
}
```

But every positive normalized source law satisfies

```math
0\le \mathfrak g[\sigma]\le 1,
```

so the upper branch is ruled out and the lower branch is the unique admissible
canonical candidate.

Even better, the lower branch is **close to natural source laws**:

- the self-matched derivative profile sits only `3.61%` away in traction,
- the exact lower branch is reached by only an `18.4%` admixture of the uniform
  positive profile into that derivative source,
- and the regular canonical point occurs already at the finite bias
  `\Pi_*\approx 1.50883`.

So the next question is no longer whether a natural source family can reach the
balanced outlet. It can. The sharper question is how the **actual coupled mouth/core
branch** deforms around that regular lower compensated point.
