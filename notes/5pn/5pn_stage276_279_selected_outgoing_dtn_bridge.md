# 5PN continuation — Stages 276–279

## Goal

Continue directly from the Stage-272–275 outgoing-normalization collapse without reopening the bulk PDE.
The next clean move was to connect the **actual selected moving-throat outgoing branch** to the isotropic DtN deformation variables

\[
(\beta,\Sigma_0,\Sigma_5),
\]

so that `chi_Q` stops being an abstract leftover scalar and becomes a direct continuum-kernel observable.

The key point that emerged immediately is that this connection is **not unique**. The exact isotropic DtN deformation algebra only fixes a two-parameter **gauge family** of triples `(beta, Sigma_0, Sigma_5)` once the observable pair `(N_Q, chi_Q)` is specified. So the right result is not “the” unique triple, but an exact family together with physically useful gauge choices.

---

## Stage 276 — exact selected-branch DtN gauge family

Start from the exact isotropic deformation laws already fixed at Stage 274:

\[
N_Q = \frac{3}{3S-\Sigma_0},
\qquad
\chi_Q = \frac{3\bigl(S\beta^5+9\Sigma_5\bigr)}{3S-\Sigma_0}.
\]

Solving exactly for the deformation data gives

\[
\boxed{
\Sigma_0 = 3S - \frac{3}{N_Q},
\qquad
\Sigma_5 = \frac{\chi_Q}{9N_Q}-\frac{S\beta^5}{9}.
}
\]

So the observable pair `(N_Q, chi_Q)` fixes only those two combinations. The actual selected outgoing branch therefore determines a **two-parameter DtN gauge family**, labeled by `(S, beta)`.

Three exact gauges are especially useful.

### A. Core gauge
Set
\[
S=1,
\qquad
\beta=1.
\]
Then
\[
\boxed{
\Sigma_0 = 3\Bigl(1-\frac{1}{N_Q}\Bigr),
\qquad
\Sigma_5 = \frac{\chi_Q-N_Q}{9N_Q}.
}
\]
This packages the whole defect into a static core shift plus an odd core outlet.

### B. Scale gauge
Set
\[
\Sigma_0=0,
\qquad
\beta=1.
\]
Then
\[
\boxed{
S = \frac{1}{N_Q},
\qquad
\Sigma_5 = \frac{\chi_Q-1}{9N_Q}.
}
\]
This packages the same selected branch into pure mouth normalization plus an odd core outlet.

### C. Argument gauge
Set
\[
\Sigma_0=0,
\qquad
\Sigma_5=0,
\]
with the natural positive branch. Then
\[
\boxed{
S = \frac{1}{N_Q},
\qquad
\beta = \chi_Q^{1/5}.
}
\]
So the same selected branch can also be viewed as a pure effective-argument deformation.

---

## Stage 277 — exact compensated Robin–mixed outlet dictionary

Now take the explicit isotropic moving-throat outlet class

\[
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+\rho_R
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}.
\]

On the exact compensated canonical-even branch

\[
\boxed{
\rho_R = 4\sigma_W,
\qquad
\kappa_W = \frac13,
}
\]

one gets the exact outgoing normalization factor

\[
\boxed{
\chi_Q^{\rm hyb} = \frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
}
\]

This explicit outlet admits two exact DtN gauge embeddings.

### Core gauge embedding
Choose
\[
S=1,
\qquad
\beta=1,
\qquad
\Sigma_0 = 3\sigma_W,
\qquad
\Sigma_5 = -\sigma_W\gamma_W,
\]
with the required even slots
\[
\Sigma_2 = -\frac{\sigma_W}{3},
\qquad
\Sigma_4 = -\frac{\sigma_W}{9}.
\]
This reproduces the compensated hybrid outlet exactly.

### Scale gauge embedding
Choose
\[
S = 1-\sigma_W,
\qquad
\beta = 1,
\qquad
\Sigma_0 = 0,
\qquad
\Sigma_5 = \sigma_W\Bigl(\frac19-\gamma_W\Bigr),
\]
with
\[
\Sigma_2 = \Sigma_4 = 0.
\]
This also reproduces the same outlet exactly.

So the compensated Robin–mixed outlet is the first explicit moving-throat outgoing branch that can be embedded into the Stage-276 DtN gauge family in more than one natural way.

A very useful special case is the canonical compensated outgoing point
\[
\gamma_W = \frac19.
\]
Then in the **scale gauge**
\[
\Sigma_5 = 0,
\qquad
\Lambda_2^{\rm hyb}(z)=(1-\sigma_W)\Lambda_2^{\rm out}(z),
\]
so the whole compensated outlet reduces to the robust pure-scale class and therefore
\[
\boxed{\chi_Q=1.}
\]

---

## Stage 278 — `chi_Q` as a direct continuum-kernel observable

The exact compensated branch makes the outgoing factor explicit:

\[
\boxed{
\chi_Q = \frac{1-9\sigma_W\gamma_W}{1-\sigma_W},
\qquad
\Delta_Q:=\chi_Q-1 = \frac{\sigma_W(1-9\gamma_W)}{1-\sigma_W}.
}
\]

So the last reduced outgoing defect is already a direct continuum observable of the selected branch.
It depends only on

- the static mixed loading `sigma_W`, and
- the odd mixed outlet coefficient `gamma_W`.

The relation can be inverted exactly:

\[
\boxed{
\gamma_W = \frac{1-(1-\sigma_W)\chi_Q}{9\sigma_W}.
}
\]

If one additionally imposes the natural source-map odd normalization condition from Stages 83–84, then the **required** conservative normalization becomes

\[
\boxed{
N_Q^{\rm req} = \frac{1}{\chi_Q} = \frac{1-\sigma_W}{1-9\sigma_W\gamma_W}.
}
\]

So the actual selected moving-throat outgoing branch no longer hides behind an abstract `chi_Q`. Once `(sigma_W, gamma_W)` are computed on the physical branch, both `chi_Q` and the required `N_Q` are known immediately.

The exact DtN-gauge observables are then:

- **core gauge**:
  \[
  \Sigma_0 = 3\sigma_W,
  \qquad
  \Sigma_5 = -\sigma_W\gamma_W,
  \]
- **scale gauge**:
  \[
  S = 1-\sigma_W,
  \qquad
  \Sigma_5 = \sigma_W\Bigl(\frac19-\gamma_W\Bigr)=\frac{(\chi_Q-1)(1-\sigma_W)}{9}.
  \]

So in the scale gauge, the entire outgoing defect is seen as one odd core slot riding on top of a pure mouth renormalization.

---

## Stage 279 — linear Family-1 canonical-even projection into `(beta, Sigma_0, Sigma_5)`

The exact linear defect formulas on the compensated hybrid outlet are

\[
\delta E_2 = \frac{\delta\mathcal C - 9\sigma_*\,\delta\kappa_W}{27(1-\sigma_*)},
\qquad
\delta E_4 = \frac{5\delta\mathcal C - 72\sigma_*\,\delta\kappa_W}{243(1-\sigma_*)},
\]

\[
\Delta_Q = \frac{\delta\mathcal C - 27\sigma_*\,\delta\gamma_W}{3(1-\sigma_*)}.
\]

Imposing exact first-order preservation of the canonical even `l=2` fingerprint,

\[
\delta E_2 = 0,
\qquad
\delta E_4 = 0,
\]

forces uniquely

\[
\boxed{
\delta\mathcal C = 0,
\qquad
\delta\kappa_W = 0.
}
\]

Using the explicit Family-1 transport law from Stage 141–142, that further implies

\[
\boxed{
\delta\mathfrak g = 0.
}
\]

So on the canonical-even compensated branch the only surviving linear mouth/core defect is the odd mixed-channel renormalization `delta gamma_W`, and the remaining outgoing defect is

\[
\boxed{
\Delta_Q = -\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
}
\]

Now compare with the Stage-92 linear DtN law

\[
\Delta_Q = 5b + \frac{a_0}{3} + 9a_5.
\]

A natural compensated **core gauge** choice is therefore

\[
\boxed{
\delta\beta = 0,
\qquad
\delta\Sigma_0 = 0,
\qquad
\delta\Sigma_5 = -\frac{\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
}
\]

Then indeed

\[
\boxed{
\Delta_Q = 9\,
\delta\Sigma_5.
}
\]

So after enforcing canonical-even preservation, the whole linearized Family-1 mouth/core defect projects into the isotropic DtN triple as a **pure odd core outlet**. That is the sharpest reduced bridge so far.

---

## Best current reading after Stage 279

The selected moving-throat outgoing branch is now connected to the isotropic DtN deformation variables in the strongest honest way available without solving the full PDE.

1. The exact isotropic deformation algebra determines only a **gauge family** of `(beta, Sigma_0, Sigma_5)` for a given `(N_Q, chi_Q)`.
2. The compensated Robin–mixed outlet gives the first explicit moving-throat branch inside that family.
3. On that explicit branch,
   
   \[
   \chi_Q = \frac{1-9\sigma_W\gamma_W}{1-\sigma_W}
   \]
   
   is already a direct continuum-kernel observable.
4. On the canonical-even Family-1 branch, every first-order mouth/core defect except the odd mixed-channel renormalization is forced away, so the DtN projection becomes
   
   \[
   \delta\beta=0,
   \qquad
   \delta\Sigma_0=0,
   \qquad
   \delta\Sigma_5=-\frac{\sigma_*}{1-\sigma_*}\delta\gamma_W.
   \]

So the remaining reduced outgoing theorem gap is now as narrow as it can be in this language:

> compute the physical odd mixed-channel renormalization `gamma_W` of the selected passive/outgoing branch, because that one continuum observable already determines the surviving isotropic DtN defect.
