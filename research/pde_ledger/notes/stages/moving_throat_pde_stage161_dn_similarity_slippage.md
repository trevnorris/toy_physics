# Moving-Throat PDE — Stage 161: D/N Similarity Decomposition of the Tangential Mixed-Port Slippage Susceptibility

## Goal

Stage 160 reduced the whole first-order co-evolving Family-1 normalization problem to one scalar susceptibility:
\[
\delta\mathfrak B_W
=
\Upsilon_\Pi\,\delta\Pi_{\rm tan},
\qquad
\delta\mathfrak B_W:=\delta\gamma_0-\frac13\delta\kappa_0.
\]

That is already a sharp theorem gate, but it still treats `\Upsilon_\Pi` as a black box.

The next honest step is to translate that black box into the actual geometric language of the concrete D/N mixed-tube realization of Stage 116.

This note does exactly that.

The main result is that the tangential slippage susceptibility is not fundamental. It is the exact image of one still smaller microscopic defect:
\[
\boxed{
\Upsilon_\Pi
=
\frac{1+r_{c,*}}{9}
\bigl(\Xi_\gamma-2\Xi_L\bigr),
}
\]
where
\[
\Xi_\gamma:=\frac{\delta\ln\gamma_0}{\delta\Pi_{\rm tan}},
\qquad
\Xi_L:=\frac{\delta\ln(L_W/a)}{\delta\Pi_{\rm tan}}.
\]

So the last linear 2.5PN obstruction is controlled only by the difference between

- the odd bare-port normalization strain `\Xi_\gamma`, and
- twice the D/N-tube length-ratio strain `2\Xi_L`.

Equivalently, if the actual tangential deformation preserves the compensated D/N similarity class
\[
\gamma_0\propto (L_W/a)^2,
\]
then
\[
\Xi_\gamma=2\Xi_L,
\qquad
\Upsilon_\Pi=0,
\qquad
\Delta_Q=0
\]
at first order.

So the remaining theorem gap narrows once again: it is no longer a generic DtN susceptibility. It is one scalar **D/N similarity-slippage susceptibility**.

---

## 1. Exact similarity-defect parametrization of the bare mixed side-channel

On the compensated canonical branch of Stages 115–116,
\[
\kappa_0=\frac{1+r_c}{3},
\qquad
\gamma_0=\frac{1+r_c}{9},
\qquad
r_c=\frac{\lambda^2}{K_sK_q}.
\]

To measure deviations away from that exact pure-scale branch, introduce the exact bare similarity-defect variables
\[
\boxed{
\kappa_0=\frac{1+r_c}{3}\,(1+\varepsilon_\kappa),
\qquad
\gamma_0=\frac{1+r_c}{9}\,(1+\varepsilon_\gamma).
}
\]

Then the bare mixed-port slippage scalar becomes
\[
\boxed{
\mathfrak B_W
:=
\gamma_0-\frac13\kappa_0
=
\frac{1+r_c}{9}\,\bigl(\varepsilon_\gamma-\varepsilon_\kappa\bigr).
}
\]

So the bare slippage is exactly the difference between two dimensionless similarity defects.

Linearizing about the compensated branch `\varepsilon_\kappa=\varepsilon_\gamma=0` gives
\[
\boxed{
\delta\mathfrak B_W
=
\frac{1+r_{c,*}}{9}
\bigl(\delta\varepsilon_\gamma-\delta\varepsilon_\kappa\bigr).
}
\]

Because the basepoint already satisfies `\varepsilon_\gamma-\varepsilon_\kappa=0`, the direct hybridization drift `\delta r_c` drops out automatically at this stage.

That is the first exact sharpening.

---

## 2. Exact D/N-tube meaning of the even similarity defect

Stage 116 gave the concrete mixed-tube realization
\[
\boxed{
\kappa_0
=
\frac{4L_W^2}{\pi^2 a^2},
}
\]
where

- `L_W` is the auxiliary D/N half-wave tube length,
- `a` is the exterior mouth radius.

Combining this with the similarity-defect parametrization above gives the exact identity
\[
\boxed{
\varepsilon_\kappa
=
\frac{12L_W^2}{\pi^2 a^2(1+r_c)}-1.
}
\]

Therefore, on the compensated branch,
\[
\delta\varepsilon_\kappa
=
2\,\frac{\delta L_W}{L_{W,*}}
-
2\,\frac{\delta a}{a_*}
-
\frac{\delta r_c}{1+r_{c,*}}.
\]
Equivalently,
\[
\boxed{
\delta\varepsilon_\kappa
=
2\,\delta\ln\!\left(\frac{L_W}{a}\right)
-
\delta\ln(1+r_c).
}
\]

So the even bare similarity defect is nothing but the failure of the D/N tube-length ratio `L_W/a` to track the hybridization factor `1+r_c` in the canonical way.

---

## 3. Exact odd similarity defect and the final cancellation of `\delta r_c`

The odd similarity-defect variable is
\[
\boxed{
\varepsilon_\gamma
=
\frac{9\gamma_0}{1+r_c}-1.
}
\]

Linearizing on the compensated branch gives
\[
\delta\varepsilon_\gamma
=
\frac{9\,\delta\gamma_0}{1+r_{c,*}}
-
\frac{\delta r_c}{1+r_{c,*}}.
\]
Since
\[
\gamma_{0,*}=rac{1+r_{c,*}}{9},
\]
this is equivalently
\[
\boxed{
\delta\varepsilon_\gamma
=
\delta\ln\gamma_0
-
\delta\ln(1+r_c).
}
\]

Subtracting the even defect from the odd defect, the hybridization drift cancels exactly:
\[
\boxed{
\delta\varepsilon_\gamma-\delta\varepsilon_\kappa
=
\delta\ln\gamma_0
-
2\,\delta\ln\!\left(\frac{L_W}{a}\right).
}
\]

This is the cleanest microscopic identity in the whole stage.

It says that the only first-order source of bare mixed-port slippage is the failure of the odd bare-port normalization `\gamma_0` to track the square of the D/N tube ratio `L_W/a`.

No direct `\delta r_c` term survives.

---

## 4. Tangential susceptibility decomposition

Now define the tangential microscopic similarity susceptibilities
\[
\boxed{
\Xi_\gamma
:=
\frac{\delta\ln\gamma_0}{\delta\Pi_{\rm tan}},
\qquad
\Xi_L
:=
\frac{\delta\ln(L_W/a)}{\delta\Pi_{\rm tan}}.
}
\]

Then the exact linear bare slippage law becomes
\[
\boxed{
\delta\mathfrak B_W
=
\frac{1+r_{c,*}}{9}
\bigl(\Xi_\gamma-2\Xi_L\bigr)
\delta\Pi_{\rm tan}.
}
\]
So the Stage 160 susceptibility is
\[
\boxed{
\Upsilon_\Pi
=
\frac{1+r_{c,*}}{9}
\bigl(\Xi_\gamma-2\Xi_L\bigr).
}
\]

Using the Family-1 identification
\[
r_{c,*}=\mathfrak r_{F1}^2,
\qquad
\mathfrak r_{F1}\approx 1.77799353547498,
\]
one finds
\[
1+r_{c,*}=1+\mathfrak r_{F1}^2\approx 4.161261012190819,
\]
so numerically
\[
\boxed{
\Upsilon_\Pi
\approx
0.462362334687869\,\bigl(\Xi_\gamma-2\Xi_L\bigr).
}
\]

Define the single remaining microscopic similarity-slippage scalar by
\[
\boxed{
\Xi_{\rm slip}:=\Xi_\gamma-2\Xi_L.
}
\]
Then
\[
\boxed{
\Upsilon_\Pi
=
\frac{1+r_{c,*}}{9}\,\Xi_{\rm slip}.
}
\]

So the unknown is no longer a generic DtN susceptibility. It is one dimensionless similarity-strain mismatch.

---

## 5. Collapse of the final defect law

Stage 160 already gave
\[
\Delta_Q
=
-
\frac{9\sigma_*}{(1-\sigma_*)(1+r_{c,*})}
\Upsilon_\Pi\,\delta\Pi_{\rm tan},
\qquad
N_Q-1
=
\frac{9\sigma_*}{(1-\sigma_*)(1+r_{c,*})}
\Upsilon_\Pi\,\delta\Pi_{\rm tan}.
\]

Substituting the exact Stage 161 decomposition, all prefactors cancel:
\[
\boxed{
\Delta_Q
=
-
\frac{\sigma_*}{1-\sigma_*}
\Xi_{\rm slip}
\,\delta\Pi_{\rm tan},
}
\]
\[
\boxed{
N_Q-1
=
\frac{\sigma_*}{1-\sigma_*}
\Xi_{\rm slip}
\,\delta\Pi_{\rm tan}.
}
\]

This is a substantial simplification.

The first-order reduced 2.5PN defect no longer depends explicitly on the core hybridization ratio `r_c` at all. It depends only on

1. the hybrid branch weight `\sigma_*`,
2. the tangential mouth deformation `\delta\Pi_{\rm tan}`,
3. the D/N similarity-slippage scalar `\Xi_{\rm slip}`.

Using the Stage 159 tangential transport,
\[
\delta\Pi_{\rm tan}
=
0.832409471081635\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S,
\]
one gets
\[
\boxed{
\Delta_Q
=
-
\frac{\sigma_*\Xi_{\rm slip}}{1-\sigma_*}
\left(
0.832409471081635\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S
\right),
}
\]
\[
\boxed{
N_Q-1
=
\frac{\sigma_*\Xi_{\rm slip}}{1-\sigma_*}
\left(
0.832409471081635\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S
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
-
\frac{\sigma_*\Xi_{\rm slip}}{1-\sigma_*}
\left(
5.35223887169622\,\delta\widehat T_m
-
1.16275838754222\,\delta\mathcal S
\right),
}
\]
\[
\boxed{
N_Q-1
\approx
\frac{\sigma_*\Xi_{\rm slip}}{1-\sigma_*}
\left(
5.35223887169622\,\delta\widehat T_m
-
1.16275838754222\,\delta\mathcal S
\right).
}
\]

So the entire first-order normalization problem has now collapsed one step further:
\[
\boxed{\Xi_{\rm slip}.}
\]

---

## 6. Exact D/N similarity-preservation theorem

Stage 116 fixed the compensated D/N-tube realization by
\[
\kappa_0=\frac{1+r_c}{3},
\qquad
\gamma_0=\frac{1+r_c}{9},
\qquad
\kappa_0=\frac{4L_W^2}{\pi^2 a^2}.
\]
Therefore, on the compensated D/N branch,
\[
\boxed{
\gamma_0
=
\frac13\kappa_0
=
\frac{4L_W^2}{3\pi^2 a^2}.
}
\]

Taking the logarithmic first variation gives the exact D/N similarity condition
\[
\boxed{
\delta\ln\gamma_0
=
2\,\delta\ln\!\left(\frac{L_W}{a}\right).
}
\]
Equivalently,
\[
\boxed{
\Xi_\gamma=2\Xi_L,
\qquad
\Xi_{\rm slip}=0.
}
\]

So if the actual tangential mouth deformation preserves the compensation-selected D/N similarity class to first order, then
\[
\boxed{
\Upsilon_\Pi=0,
\qquad
\delta\mathfrak B_W=0,
\qquad
\delta\gamma_W=0,
\qquad
\Delta_Q=0,
\qquad
N_Q-1=0.
}
\]

This is the strongest first-order positive result reached so far.

It says that the reduced 2.5PN obstruction survives only if the tangential co-evolving mouth deformation drives the concrete mixed side-channel **out of its compensation-selected D/N similarity class**.

---

## 7. Best current theorem statement after Stage 161

Inside the explicit co-evolving Family-1 closure, the compensated Robin–mixed outlet class, and the concrete D/N-tube realization of the bare mixed side-channel:

1. exact canonical-even preservation still forces
   \[
   \delta\mathfrak g=0,
   \qquad
   \delta\kappa_W=0,
   \]
2. the only remaining first-order odd datum is the bare mixed-port slippage
   \[
   \delta\mathfrak B_W,
   \]
3. that slippage is exactly the image of one D/N similarity-strain mismatch
   \[
   \Xi_{\rm slip}=\Xi_\gamma-2\Xi_L,
   \]
4. and the reduced 2.5PN defect is therefore
   \[
   \Delta_Q
   =
   -\frac{\sigma_*}{1-\sigma_*}
   \Xi_{\rm slip}
   \delta\Pi_{\rm tan}.
   \]

So the next theorem gate is now as small as it can be without solving the true moving-throat branch:

> compute whether the actual tangential co-evolving mouth deformation preserves the compensated D/N similarity law
> \(
> \gamma_0 = 4L_W^2/(3\pi^2 a^2)
> \)
> to first order, or produces a nonzero similarity-slippage scalar `\Xi_{\rm slip}`.

That is the exact next microscopic closure question.
