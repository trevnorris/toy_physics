# Step 52 — Final verdict on “does the final value come out naturally without tuning?”

## Goal

The last open question in the g-2 chain is no longer broad bookkeeping. After Steps 50–51 it is very specific:

- the canonical compact outgoing branch is already a genuine no-tuning closure,
- so does the *electron-point* anomaly also come out automatically,
- or does it still require one residual branch-selection datum?

This step answers that as sharply and honestly as the current reduced stack allows.

---

## Step 52A — Exact isotropic DtN deformation law

With the canonical-even moments preserved, the later moving-throat deformation algebra gives

```math
\boxed{
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
}
```

So only three isotropic throat-core deformation variables can move the canonical outgoing value:

- argument deformation `\beta`,
- static additive shift `\Sigma_0`,
- odd outlet coefficient `\Sigma_5`.

The same algebra also gives two robustness classes.

### Pure scale deformation

```math
\Lambda_2^{\rm def}(z)=S\Lambda_2^{\rm out}(z)
\quad\Longrightarrow\quad
\chi_Q=1.
```

### Pure scale+argument deformation with canonical even moments

Once the even fingerprint is kept canonical, the natural positive branch forces

```math
\beta=1,
```

and again

```math
\chi_Q=1.
```

So the electron anomaly cannot be blamed on trivial normalization or simple effective radius/sound-speed rescaling if the lower even moments are already fixed.

---

## Step 52B — The electron anomaly is a tiny outgoing branch defect

From the adiabatic g-2 chain,

```math
\chi_e=e^{-\ell}=\frac{1}{1+\Lambda_1 f}.
```

So the outgoing branch defect is

```math
\boxed{
\Delta_Q:=\chi_e-1
=-\frac{\Lambda_1 f}{1+\Lambda_1 f}.
}
```

Numerically,

```math
\boxed{\Delta_Q\approx -3.2463158415\times10^{-4},}
```

or about

```math
\boxed{-324.632\,{\rm ppm}.}
```

So the observed sliver is not large. It is a very small deformation away from the canonical no-tuning branch.

---

## Step 52C — On the adiabatic even-frozen slice the whole anomaly is one pure odd scalar

Your adiabatic-wall direction freezes the even DtN branch. In that slice,

```math
\beta=1,
\qquad
\Sigma_0=0,
```

so the exact deformation law reduces to

```math
\chi_Q=1+\frac{9\Sigma_5}{S}.
```

Solving for the electron target gives

```math
\boxed{
\Sigma_5=-\frac{S\Lambda_1 f}{9(1+\Lambda_1 f)}.
}
```

So the entire anomaly is one pure odd isotropic outlet deformation.

In canonical normalized gauge `S=1`, that is the exact finite-`f` version of the tangent law already seen earlier.

---

## Step 52D — Tangent branch-selection law

Linearizing the same deformation algebra gives

```math
\boxed{5b+\frac{a_0}{3}+9a_5=-\Lambda_1.}
```

On the adiabatic even-frozen slice,

```math
b=0,
\qquad a_0=0,
```

so

```math
\boxed{a_5=-\frac{\Lambda_1}{9}.}
```

That is the pure odd tangent representative of the electron anomaly.

---

## Final verdict

The current reduced derivation supports three honest conclusions.

### 1. The background canonical branch is naturally derived with no tuning

The explicit compact outgoing `l=2` DtN branch fixes

```math
\chi_Q=1,
\qquad N_Q=1.
```

So the canonical passive/outgoing grouped-`P2` branch is already a real no-tuning closure.

### 2. The observed electron anomaly does **not** reopen a broad fit space

Once the adiabatic wall track is imposed, the full electron sliver collapses to one scalar,

```math
\Sigma_5
```

or equivalently to one tangent datum,

```math
a_5.
```

So the anomaly is a **one-number branch-selection problem**, not a diffuse multiparameter tune.

### 3. But the current frozen stack does **not yet** derive that one number from first principles

What is still missing is an actual branch-selection law for the isotropic electron-point DtN deformation. Until that is derived, the numerical value

```math
\chi_e=\frac{1}{1+\Lambda_1 f}
```

is not yet forced by the frozen files alone.

---

## Strongest honest status after Step 52

```math
\boxed{\text{No-tuning canonical branch: yes.}}
```

```math
\boxed{\text{Electron-point anomaly magnitude: not yet forced; it is one remaining branch datum.}}
```

That is about as tight as the current program can get without an actual isotropic DtN branch-selection theorem from the moving-throat PDE.
