# Step 43 — The adiabatic anomaly admits a pure odd normalized DtN representative

## Goal

Step 42 showed that on the anomaly-only adiabatic branch

- the parent compensation ratios stay fixed,
- the auxiliary D/N geometry stays fixed,
- and only the odd outlet coefficient detunes.

That still lives in the **microscopic** compensated-core / hybrid-outlet variables.
The next honest question is therefore:

> what does the same anomaly look like after we reduce it all the way back to the **normalized isotropic DtN branch-selection surface**?

This step shows an unexpectedly clean answer:

```math
\boxed{
\text{the entire adiabatic anomaly has an exact pure-odd normalized DtN representative.}
}
```

In that representative the even branch is frozen exactly and the whole observable branch motion sits in one isotropic odd slot.

---

## Step 43A — Exact isotropic DtN deformation surface

From the earlier isotropic DtN reduction,

```math
\chi_Q
=
\frac{3\bigl(S\beta^5+9\Sigma_5\bigr)}{3S-\Sigma_0},
```

with

- `S` = overall isotropic mouth normalization,
- `\beta` = outgoing-argument deformation,
- `\Sigma_0` = static isotropic additive deformation,
- `\Sigma_5` = isotropic odd `l=2` core outlet.

The adiabatic anomaly target is still

```math
\chi_Q = e^{-\ell} = \frac{1}{1+x},
\qquad x:=e^{\ell}-1=\Lambda_1 f.
```

So the problem is now very narrow: how is this exact target represented on the DtN surface once the adiabatic anomaly-only closure freezes the even branch?

---

## Step 43B — Fixed-even slice gives a pure odd exact law

Because Step 42 already proved that the anomaly-only branch does **not** move the parent compensation ratios or the even D/N geometry, the natural normalized DtN slice is

```math
\boxed{\beta = 1,\qquad \Sigma_0 = 0.}
```

Then the isotropic DtN surface collapses to

```math
\chi_Q = 1 + \frac{9\Sigma_5}{S}.
```

Matching this to `\chi_Q=e^{-\ell}` gives the exact pure-odd law

```math
\boxed{
\Sigma_5(\ell)
=
\frac{S}{9}\bigl(e^{-\ell}-1\bigr)
=
-\frac{Sx}{9(1+x)}.
}
```

So once the even branch is frozen, the exact electron anomaly is represented by **one** pure odd isotropic outlet coefficient and nothing else.

---

## Step 43C — This pure odd representative is independent of the microscopic loading share

On the compensated hybrid outlet we also have

```math
\chi_Q = \frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
```

Step 42 showed that `\sigma_W` drifts and `\gamma_W` retunes so that this still equals the electron target.

Now form the normalized DtN representative directly from `\chi_Q`:

```math
\Sigma_5 = \frac{S}{9}(\chi_Q-1).
```

After substituting the exact hybrid-outlet solution for `\gamma_W`, all explicit dependence on `\sigma_W` cancels, leaving

```math
\boxed{
\Sigma_5(\ell)=\frac{S}{9}(e^{-\ell}-1).
}
```

So the microscopic loading-share drift is real on the compensated-core side, but once the branch is repackaged as a **normalized outgoing DtN deformation**, the observable anomaly is carried entirely by the pure odd slot.

That is the cleanest reduced statement reached so far.

---

## Step 43D — Canonical normalized gauge and the tangent law

If we choose the canonical normalized gauge

```math
S=1,
```

then the exact finite-`f` odd slot is

```math
\boxed{
\Sigma_5
=
-\frac{\Lambda_1 f}{9(1+\Lambda_1 f)}.
}
```

Numerically at the electron point,

```math
\boxed{
\Sigma_5^{(e)}
\approx
-3.60701760168546\times 10^{-5}.
}
```

Now expand the full isotropic DtN surface with

```math
S=1,
\qquad
\beta = 1+\varepsilon b,
\qquad
\Sigma_0 = \varepsilon a_0,
\qquad
\Sigma_5 = \varepsilon a_5,
\qquad
\chi_Q = 1 - \varepsilon\Lambda_1 + O(\varepsilon^2).
```

The exact first-order constraint is the already-known tangent law

```math
5b + \frac{a_0}{3} + 9a_5 = -\Lambda_1.
```

On the fixed-even adiabatic slice `b=0`, `a_0=0`, this reduces to

```math
\boxed{a_5 = -\frac{\Lambda_1}{9}.}
```

So the adiabatic anomaly-only branch picks the pure odd tangent representative

```math
\boxed{b=0,\qquad a_0=0,\qquad a_5=-\Lambda_1/9.}
```

---

## Step 43E — What changed conceptually

There are now two equivalent descriptions of the same quartic anomaly layer.

### Microscopic compensated-core / hybrid-outlet language

The anomaly slightly softens the loading share and slightly retunes the odd outlet coefficient, while leaving the parent compensation ratios and even geometry fixed.

### Normalized isotropic DtN language

The exact same branch motion can be represented by

```math
\boxed{
\beta = 1,
\qquad
\Sigma_0 = 0,
\qquad
\Sigma_5 = \frac{S}{9}(e^{-\ell}-1).
}
```

So at the normalized outgoing-DtN level, the anomaly is **pure odd**.

That is a very useful simplification for the next moving-throat step, because it means the PDE-facing target is no longer “find a mix of even and odd isotropic deformations.” It is just:

> compute the physical odd quadrupole outlet coefficient on the true branch and see whether it lands on the pure-odd target above.

---

## Main result of the step

The adiabatic anomaly-only branch admits an exact pure-odd normalized DtN representative:

```math
\boxed{
\Sigma_5(\ell)=\frac{S}{9}(e^{-\ell}-1)=-\frac{S\Lambda_1 f}{9(1+\Lambda_1 f)},
\qquad
\beta=1,
\qquad
\Sigma_0=0.
}
```

In canonical normalized gauge,

```math
\boxed{
\Sigma_5
=
-\frac{\Lambda_1 f}{9(1+\Lambda_1 f)},
}
```

and at tangent level,

```math
\boxed{b=0,\qquad a_0=0,\qquad a_5=-\Lambda_1/9.}
```

So the improved g-2 picture is now extremely compressed:

- the adiabatic ground state freezes the wall and the parent compensation sheet,
- the microscopic core/hybrid outlet retunes only its odd slot,
- and the normalized outgoing DtN branch sees the whole anomaly as one exact pure-odd isotropic coefficient.
