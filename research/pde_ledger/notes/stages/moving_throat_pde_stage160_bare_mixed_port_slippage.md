
# Moving-Throat PDE — Stage 245: Bare Mixed-Port Slippage Theorem and the Collapse of the Last Tangential DtN Gate

## Goal

Stage 244 reduced the co-evolving Family-1 problem to one exact first-order obstruction:
\[
\Delta_Q
=
-\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
\]

So the next real question is no longer broad.

It is:

> what microscopic quantity inside the concrete compensated Robin–mixed core
> actually controls `\delta\gamma_W` under a tangential co-evolving mouth
> deformation?

This note gives that bridge exactly.

The main result is that, on the compensated canonical branch of the concrete
two-channel core model, the effective odd defect `\delta\gamma_W` is controlled
only by one **bare mixed-port slippage scalar**
\[
\delta\mathfrak B_W
:=
\delta\gamma_0-\frac13\,\delta\kappa_0.
\]

Equivalently,
\[
\boxed{
\delta\gamma_W
=
\frac{\delta\mathfrak B_W}{1+r_{c,*}}
=
\frac{\delta\mathfrak B_W}{1+\mathfrak r_{F1}^2}.
}
\]

So the last linear 2.5PN obstruction is not “any tangential mouth deformation.”
It is only the part of that deformation that breaks the pure-scale locking of the
bare mixed side-channel.

If the bare mixed branch stays a pure-scale deformation of the canonical compact
outgoing `l=2` port, then
\[
\delta\mathfrak B_W=0,
\qquad
\delta\gamma_W=0,
\qquad
\Delta_Q=0
\]
at first order.

That is the sharpest narrowing reached so far.

---

## 1. Exact core identities on the compensated branch

Stage 97 gave the concrete core-level outlet coefficients
\[
\kappa_W=\frac{\kappa_0}{1+r_c},
\qquad
\gamma_W=\frac{\gamma_0}{1+r_c},
\qquad
r_c=\frac{\lambda^2}{K_sK_q}.
\]

Stage 98 showed that the compensated canonical branch is characterized by
\[
\kappa_W=\frac13,
\qquad
\gamma_W=\frac19,
\]
which at the bare level means
\[
\kappa_0=\frac{1+r_c}{3},
\qquad
\gamma_0=\frac{1+r_c}{9}.
\]

Now linearize around that compensated branch:
\[
r_c=r_{c,*}+\delta r_c,
\qquad
\kappa_0=\frac{1+r_{c,*}}{3}+\delta\kappa_0,
\qquad
\gamma_0=\frac{1+r_{c,*}}{9}+\delta\gamma_0.
\]

Then
\[
\delta\kappa_W
=
\frac{\delta\kappa_0}{1+r_{c,*}}
-
\frac{\kappa_{0,*}}{(1+r_{c,*})^2}\,\delta r_c,
\]
\[
\delta\gamma_W
=
\frac{\delta\gamma_0}{1+r_{c,*}}
-
\frac{\gamma_{0,*}}{(1+r_{c,*})^2}\,\delta r_c.
\]

Using
\[
\kappa_{0,*}=\frac{1+r_{c,*}}{3},
\qquad
\gamma_{0,*}=\frac{1+r_{c,*}}{9},
\]
one gets the exact compensated-branch identity
\[
\boxed{
\delta\gamma_W-\frac13\,\delta\kappa_W
=
\frac{1}{1+r_{c,*}}
\left(
\delta\gamma_0-\frac13\,\delta\kappa_0
\right).
}
\]

This is the key algebraic bridge.

The dependence on `\delta r_c` cancels identically.
So the odd/even mismatch is insensitive to first-order drift of the static/mixed
hybridization itself.

---

## 2. Collapse under the Stage 244 canonical-even gate

Stage 244 already proved that exact first-order preservation of the canonical even
`l=2` fingerprint forces
\[
\delta\kappa_W=0,
\qquad
\delta\mathfrak g=0.
\]

So the identity above collapses immediately to
\[
\boxed{
\delta\gamma_W
=
\frac{1}{1+r_{c,*}}
\left(
\delta\gamma_0-\frac13\,\delta\kappa_0
\right).
}
\]

Using the Family-1 identification
\[
r_{c,*}=\mathfrak r_{F1}^2,
\qquad
\mathfrak r_{F1}\approx 1.77799353547498,
\]
this becomes
\[
\boxed{
\delta\gamma_W
=
\frac{\delta\mathfrak B_W}{1+\mathfrak r_{F1}^2},
\qquad
\delta\mathfrak B_W
:=
\delta\gamma_0-\frac13\,\delta\kappa_0.
}
\]

Numerically,
\[
1+\mathfrak r_{F1}^2\approx 4.161261012190819,
\]
so
\[
\boxed{
\delta\gamma_W
\approx
0.240311770175051\,\delta\mathfrak B_W.
}
\]

So the last effective odd defect is carried by one scalar only:
the **bare mixed-port slippage**
\[
\delta\mathfrak B_W.
\]

---

## 3. Pure-scale harmlessness theorem

The compensated concrete core was built from a bare mixed side-channel that, on the
canonical branch, is a pure-scale deformation of the compact outgoing `l=2` port.

At first order, pure-scale preservation means that the bare even and odd coefficients
continue to move in the canonical ratio
\[
\delta\gamma_0=\frac13\,\delta\kappa_0.
\]

Then
\[
\delta\mathfrak B_W=0,
\]
and therefore
\[
\boxed{
\delta\gamma_W=0.
}
\]

Combining with the Stage 244 defect law gives
\[
\boxed{
\Delta_Q=0,
\qquad
N_Q-1=0
}
\]
to first order on the natural source-map branch.

So the tangential mouth deformation is not dangerous by itself.
It is dangerous only if it drives the bare mixed side-channel away from pure-scale
locking.

This is the exact reduced form of the remaining theorem gate.

---

## 4. Tangential DtN susceptibility

To connect this to the actual co-evolving mouth data, define the tangential
mixed-port slippage susceptibility by
\[
\boxed{
\delta\mathfrak B_W
=
\Upsilon_\Pi\,\delta\Pi_{\rm tan},
}
\]
where the Stage 244 tangential mouth deformation is
\[
\delta\Pi_{\rm tan}
=
0.832409471081635\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S
\qquad
(\delta\mathfrak g=0).
\]

Then
\[
\boxed{
\delta\gamma_W
=
\frac{\Upsilon_\Pi}{1+\mathfrak r_{F1}^2}\,
\delta\Pi_{\rm tan}.
}
\]

This is the cleanest possible reduced DtN interface between the co-evolving mouth
system and the outlet-core normalization problem.

Everything microscopic that the full moving-throat PDE still has to compute is now
packaged into one scalar susceptibility `\Upsilon_\Pi`.

---

## 5. Final linear defect law in mouth variables

Using the Stage 244 result
\[
\Delta_Q
=
-\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W,
\qquad
N_Q-1
=
\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W,
\]
we obtain
\[
\boxed{
\Delta_Q
=
-\frac{9\sigma_*\,\Upsilon_\Pi}
{(1-\sigma_*)(1+\mathfrak r_{F1}^2)}
\,\delta\Pi_{\rm tan},
}
\]
\[
\boxed{
N_Q-1
=
\frac{9\sigma_*\,\Upsilon_\Pi}
{(1-\sigma_*)(1+\mathfrak r_{F1}^2)}
\,\delta\Pi_{\rm tan}.
}
\]

Numerically,
\[
\frac{9}{1+\mathfrak r_{F1}^2}
\approx
2.16280593157546,
\]
so
\[
\boxed{
\Delta_Q
\approx
-\frac{2.16280593157546\,\sigma_*\,\Upsilon_\Pi}{1-\sigma_*}
\,\delta\Pi_{\rm tan},
}
\]
\[
\boxed{
N_Q-1
\approx
\frac{2.16280593157546\,\sigma_*\,\Upsilon_\Pi}{1-\sigma_*}
\,\delta\Pi_{\rm tan}.
}
\]

Substituting the explicit tangential transport,
\[
\boxed{
\Delta_Q
\approx
-\frac{\sigma_*\,\Upsilon_\Pi}{1-\sigma_*}
\left(
1.80034014155495\,\delta\Sigma_0
-
2.51482073756543\,\delta\mathcal S
\right),
}
\]
\[
\boxed{
N_Q-1
\approx
\frac{\sigma_*\,\Upsilon_\Pi}{1-\sigma_*}
\left(
1.80034014155495\,\delta\Sigma_0
-
2.51482073756543\,\delta\mathcal S
\right).
}
\]

Using
\[
\delta\Sigma_0
\approx
6.42981496203006\,\delta\widehat T_m,
\]
this is equivalently
\[
\boxed{
\Delta_Q
\approx
-\frac{\sigma_*\,\Upsilon_\Pi}{1-\sigma_*}
\left(
11.5758539789133\,\delta\widehat T_m
-
2.51482073756543\,\delta\mathcal S
\right),
}
\]
\[
\boxed{
N_Q-1
\approx
\frac{\sigma_*\,\Upsilon_\Pi}{1-\sigma_*}
\left(
11.5758539789133\,\delta\widehat T_m
-
2.51482073756543\,\delta\mathcal S
\right).
}
\]

So the whole first-order 2.5PN normalization problem has now collapsed to a single
scalar:
\[
\boxed{\Upsilon_\Pi.}
\]

---

## 6. Best current theorem statement after Stage 245

Inside the explicit co-evolving Family-1 closure and the compensated Robin–mixed
outlet class, the sequence of reductions is now:

1. exact first-order canonical-even preservation forces
   \[
   \delta\mathfrak g=0,
   \qquad
   \delta\kappa_W=0,
   \]
2. the concrete outlet-core model then shows that the only remaining odd defect is
   the bare mixed-port slippage
   \[
   \delta\mathfrak B_W
   =
   \delta\gamma_0-\frac13\,\delta\kappa_0,
   \]
3. pure-scale preservation of the bare mixed side-channel kills that slippage
   exactly,
4. and the completed moving-throat PDE now has only one genuine new linear datum to
   supply:
   \[
   \delta\mathfrak B_W
   =
   \Upsilon_\Pi\,\delta\Pi_{\rm tan}.
   \]

So the next serious calculation is not another outlet-family search.
It is the much smaller question:

> compute the tangential mixed-port slippage susceptibility `\Upsilon_\Pi` on the
> actual moving-throat branch.

If that susceptibility vanishes, the first-order reduced 2.5PN defect vanishes
automatically.

## Immediate next step

The clean next derivation is now:

1. keep the concrete compensated core model,
2. perturb the bare mixed side-channel by the tangential mouth deformation,
3. compute the induced slippage scalar
   `\delta\mathfrak B_W = \delta\gamma_0 - \delta\kappa_0/3`,
4. and extract the coefficient `\Upsilon_\Pi`.

That is the exact next theorem gate.
