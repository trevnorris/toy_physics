# Step 31 — The quartic g-2 sliver becomes an exact finite-`f` isotropic DtN surface

## Goal

Step 30 fixed only the **tangent** isotropic DtN branch-selection law

```math
5b+\frac{a_0}{3}+9a_5=-\Lambda_1.
```

The moving-throat notes already give the corresponding **exact** isotropic DtN
deformation law, after the canonical even `l=2` fingerprint is enforced,

```math
\chi_Q
=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
```

So the next honest move is to stop working only at tangent level and insert the
exact Step-29 electron-point target

```math
\chi_e=\frac{1}{1+\Lambda_1 f}.
```

That is what this step does.

The main result is the exact finite-`f` isotropic anomaly surface

```math
\boxed{
3\big(S\beta^5+9\Sigma_5\big)(1+\Lambda_1 f)=3S-\Sigma_0.
}
```

Equivalently,

```math
\boxed{
3S\!\left[(1+\Lambda_1 f)\beta^5-1\right]
+\Sigma_0
+27(1+\Lambda_1 f)\Sigma_5
=0.
}
```

So the quartic g-2 sliver is no longer only a first-order branch-selection
constraint. It is one exact finite-`f` isotropic DtN surface.

---

## Step 31A — Exact electron-point target

From Step 29,

```math
\Delta_Q^{(e)}
=
-\frac{\Lambda_1 f}{1+\Lambda_1 f},
\qquad
\chi_e = 1+\Delta_Q^{(e)} = \frac{1}{1+\Lambda_1 f}.
```

Using the carried benchmark numbers,

```math
\Lambda_1 \approx 0.279605891931464,
\qquad
f \approx 0.001161409732093,
```

one has

```math
\boxed{
\Lambda_1 f \approx 3.24737004039746\times 10^{-4},
}
```

and therefore

```math
\boxed{
\chi_e \approx 0.999675368415848.
}
```

So the electron-point branch is only a very small finite deformation away from
the canonical outgoing value `\chi_Q=1`.

---

## Step 31B — Exact isotropic DtN anomaly surface

The Stage-90/91 deformation law is

```math
\chi_Q
=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0},
```

where

- `S` is the overall mouth normalization,
- `\beta` is the outgoing-argument deformation,
- `\Sigma_0` is the static isotropic throat-core shift,
- `\Sigma_5` is the extra odd isotropic `l=2` core outlet.

Matching `\chi_Q` to the exact electron target `\chi_e` gives

```math
\boxed{
3\big(S\beta^5+9\Sigma_5\big)(1+\Lambda_1 f)=3S-\Sigma_0.
}
```

Solving for the odd core slot gives the exact branch-selection surface

```math
\boxed{
\Sigma_5
=
\frac{3S-\Sigma_0}{27(1+\Lambda_1 f)}
-\frac{S\beta^5}{9}.
}
```

This is the finite-`f` replacement for the linear Step-30 tangent relation.

---

## Step 31C — Canonical-preservation surface as the zero-anomaly limit

If one instead asks for exact preservation of the canonical outgoing branch,
`\chi_Q=1`, then the same algebra gives

```math
\boxed{
\Sigma_5
=
\frac{S(1-\beta^5)}{9}
-\frac{\Sigma_0}{27}.
}
```

So the exact anomaly surface is a direct deformation of the exact
canonical-preservation surface, obtained simply by replacing `1` by
`1/(1+\Lambda_1 f)`.

That is a good consistency check: at `f\to 0`, the g-2 branch-selection law
collapses back to the canonical outgoing-preservation law.

---

## Step 31D — Three exact pure isotropic realizations

The exact finite-`f` surface already admits three clean pure realizations.

### Pure compensated argument deformation

Set `\Sigma_0=\Sigma_5=0`. Then

```math
\boxed{
\beta = \chi_e^{1/5} = (1+\Lambda_1 f)^{-1/5}.
}
```

Numerically,

```math
\boxed{
\beta_e \approx 0.999935065250674.
}
```

So the required argument deformation is tiny.

### Pure static additive core

Set `\beta=1` and `\Sigma_5=0`. Then

```math
\boxed{
\Sigma_0 = -3S\,\Lambda_1 f.
}
```

At unit scale `S=1`,

```math
\boxed{
\Sigma_0^{(e)} \approx -9.74211012119238\times 10^{-4}.
}
```

### Pure odd `l=2` core outlet

Set `\beta=1` and `\Sigma_0=0`. Then

```math
\boxed{
\Sigma_5 = -\frac{S\,\Lambda_1 f}{9(1+\Lambda_1 f)}.
}
```

At unit scale `S=1`,

```math
\boxed{
\Sigma_5^{(e)} \approx -3.60701760168546\times 10^{-5}.
}
```

So the exact electron-point outgoing defect can be carried by a very small
odd-core outlet as well.

---

## Step 31E — Exact surface reduces to the Step-30 tangent law

Write

```math
S=1+\varepsilon s,
\qquad
\beta=1+\varepsilon b,
\qquad
\Sigma_0=\varepsilon a_0,
\qquad
\Sigma_5=\varepsilon a_5,
```

and expand the exact surface together with

```math
\chi_e=\frac{1}{1+\varepsilon\Lambda_1}.
```

The first-order condition is exactly

```math
\boxed{
5b+\frac{a_0}{3}+9a_5+\Lambda_1=0,
}
```

i.e.

```math
\boxed{
5b+\frac{a_0}{3}+9a_5=-\Lambda_1.
}
```

So Step 30 is the strict tangent reduction of the present exact finite-`f`
branch-selection law.

---

## Main result of the step

The quartic g-2 sliver is now encoded by the exact isotropic DtN surface

```math
\boxed{
3\big(S\beta^5+9\Sigma_5\big)(1+\Lambda_1 f)=3S-\Sigma_0.
}
```

with exact pure-branch representatives

```math
\boxed{
\beta = (1+\Lambda_1 f)^{-1/5},
\qquad
\Sigma_0=-3S\Lambda_1 f,
\qquad
\Sigma_5=-\frac{S\Lambda_1 f}{9(1+\Lambda_1 f)}.
}
```

So the next PDE-facing question becomes sharper again:

> **which explicit moving-throat outlet model can realize this exact finite-`f` surface while preserving the already-fixed conservative even `l=2` branch?**
