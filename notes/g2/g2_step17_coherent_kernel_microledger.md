# Step 17 — Coherent-kernel microscopic demand-ratio ledger

## Goal

Step 16 reduced the quartic branch-selection question to one scalar criterion:

```math

delta\ln R_{\rm target}=0
\qquad\text{or}\qquad
\delta\ln R_{\rm target}=-\Lambda_1.
```

But that was still an abstract branch-composite statement.
The natural next move is to insert the **actual coherent local D/N kernel formulas** from the moving-throat notes and ask what the microscopic variables really control.

This step does that.

The main result is stronger than expected:

```math
\boxed{
\Lambda_1
=
-\delta\ln\Lambda
-2\,\delta\ln(1-\epsilon)
+\delta\ln Z_W
+2\,\delta\ln(1+\chi_0)
}
```

so the quartic universal transfer-shape drift is controlled entirely by the
**mixed/outgoing microscopic variables**.
The wall–`U` dressing coordinate `\epsilon_\eta` cancels identically from the inferred `\Lambda_1` law.

At the same time, the coherent support ratio `\zeta` does not enter either
`R_{\rm tr}` or `R_{\rm target}` at all. It only moves the total baseline `M_{\rm tr}`.

So the next branch-selection step is no longer “what does the support lane do to the target ratio?”
It does nothing to the target ratio directly.
The live question is whether the actual support lane compensates the exact tracking deficit by increasing the baseline at fixed target, or whether the completed PDE forces a genuine retargeting drift in the mixed/outgoing microdata.

---

## Inputs from the coherent local D/N kernel

From the coherent-kernel placement map in the moving-throat notes,

```math
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
```

```math
\epsilon
=\epsilon_W^{(\mathrm{split})}
=\epsilon_W\Bigl[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\Bigr],
```

```math
M_{\rm mix}
=
\frac{8 Z_W (1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
```

```math
M_{\rm supp}
=
\frac{8 \zeta Z_W (1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\zeta\epsilon)},
```

```math
M_{\rm tr}=M_{\rm mix}+M_{\rm supp}=M_{\rm mix}S(\zeta;\epsilon),
```
with
```math
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
```
and
```math
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
```

These are exactly the formulas that Stage 30 compresses the coherent local kernel into.

---

## Step 17A — The support ratio drops out of both `R_tr` and `R_target`

The script verifies

```math
\frac{\partial R_{\rm tr}}{\partial \zeta}=0,
\qquad
\frac{\partial R_{\rm target}}{\partial \zeta}=0.
```

So on the actual coherent local-kernel branch:

- `\zeta` does **not** move the tracking factor,
- `\zeta` does **not** move the demand ratio,
- `\zeta` only moves the total baseline through
  ```math
  M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon).
  ```

That already tells us something important:
whatever the support lane is doing, it is **not** directly retargeting the selected-branch demand ratio.

---

## Step 17B — Exact microscopic drift of the tracking factor

Use the natural logarithmic drift variables

```math
q_\chi:=\delta\ln(1+\chi_0),
\qquad
q_U:=\delta\ln(1+\delta_U).
```

Then the exact logarithmic drift of the tracking factor is

```math
\boxed{
\delta\ln R_{\rm tr}
=
-\frac{\delta_U\,q_\chi+\chi_0\,q_U}{1+\chi_0+\delta_U}.
}
```

So the constructive branch lowers `R_{\rm tr}` only through motion in the two microscopic placement variables `(1+\chi_0)` and `(1+\delta_U)`.
The support ratio `\zeta` is absent.

---

## Step 17C — Exact microscopic drift of the split blocking ratio

If

```math
q_W:=\delta\ln\epsilon_W,
```
then the split blocking ratio obeys

```math
\boxed{
q_\epsilon:=\delta\ln\epsilon
=
q_W-\frac{2}{11+9\delta_U}\,q_U.
}
```

So the split-blocking drift is fixed by the primitive mixed blocking drift plus the axial split drift.

---

## Step 17D — Exact microscopic drift of the demand ratio

Using the natural drifts

```math
q_\Lambda:=\delta\ln\Lambda,
\qquad
q_Z:=\delta\ln Z_W,
\qquad
q_\eta:=\delta\ln\epsilon_\eta,
\qquad
q_\epsilon:=\delta\ln\epsilon,
```

we get the exact microscopic demand-ratio law

```math
\boxed{
\delta\ln R_{\rm target}
=
q_\Lambda-q_Z-2q_\chi
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta
-2\frac{\epsilon_*}{1-\epsilon_*}q_\epsilon.
}
```

Equivalently,

```math
\delta\ln R_{\rm target}
=
\delta\ln(1-\epsilon_\eta)
+\delta\ln\Lambda
+2\,\delta\ln(1-\epsilon)
-\delta\ln Z_W
-2\,\delta\ln(1+\chi_0).
```

So the actual coherent-kernel demand ratio is controlled by

- the common wall–`U` dressing factor `(1-\epsilon_\eta)`,
- the mixed/outgoing scale `\Lambda`,
- the split blocking `(1-\epsilon)`,
- the wall-to-mixed overlap `Z_W`,
- and the common interference ratio `(1+\chi_0)`.

Again, `\zeta` is absent.

---

## Step 17E — Exact cancellation theorem with the Step-16 law

Step 16 proved abstractly that

```math
\delta\ln R_{\rm target}
=
\delta\ln(1-\epsilon_\eta)-\Lambda_1.
```

Now insert the exact microscopic Stage-30 formula above. The `\delta\ln(1-\epsilon_\eta)` term cancels identically, leaving

```math
\boxed{
\Lambda_1
=
-\delta\ln\Lambda
-2\,\delta\ln(1-\epsilon)
+\delta\ln Z_W
+2\,\delta\ln(1+\chi_0).
}
```

Or, in the logarithmic drift variables,

```math
\boxed{
\Lambda_1
=
-q_\Lambda+q_Z+2q_\chi+2\frac{\epsilon_*}{1-\epsilon_*}q_\epsilon.
}
```

This is the main result of the step.

The quartic universal transfer-shape drift is **not** controlled by

- the support ratio `\zeta`, or
- the wall–`U` dressing coordinate `\epsilon_\eta`.

It is carried entirely by the mixed/outgoing microscopic variables.

---

## Step 17F — What the cancellation means physically

This exact cancellation changes the interpretation of the remaining anomaly problem.

The branch-selection ambiguity is no longer “support versus dressing.”
The support lane does not move `R_{\rm target}` at all, and the wall–`U` dressing coordinate cancels from the inferred `\Lambda_1` law.

So the next honest question is now much sharper:

> does the physical coherent support lane compensate the exact tracking-branch deficit by increasing the available baseline at fixed demand ratio, or does the completed PDE force a real retargeting drift in the mixed/outgoing microdata `(\Lambda,Z_W,\chi_0,\epsilon)`?

That is the direct continuation point.

---

## Continuation point

The next clean move is to use the exact support-compensation theorem from the coherent local D/N branch and see whether the physical support lane is naturally a **fixed-target / load-compensation** mechanism.

If it is, then the real PDE-side continuation will favor the Step-16 coherent side

```math
\delta\ln R_{\rm target}=0,
```

rather than the direct retargeting law.
