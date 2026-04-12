# 5PN continuation notes — stages 231 through 233

These stages continue directly from the Stage-228/230 two-selector reduction, but now tie it back to the **actual microscopic coherent moving-throat variables** rather than leaving it as an abstract overlap-image test.

The payoff is strong:

> on the minimal coherent continuum branch, the Stage-230 overlap-image selectors are no longer abstract. They collapse to **two explicit microscopic defect scalars**.

Those two scalars are:

- the **dressing slippage**
  \[
  \Sigma_\eta = 2c_1-\kappa_U-\kappa_\eta,
  \]
- the **support-plane scalar**
  \[
  \Sigma_{\rm sup}
  =
  \frac{\kappa_U}{2}
  -c_K\kappa_\eta
  +c_\chi\Sigma_\chi
  +c_Z\zeta_Z,
  \]
  with
  \[
  c_K=\frac{8575}{4392}\approx 1.95241347905,
  \]
  \[
  c_\chi\approx 14.3144730533,
  \qquad
  c_Z\approx 1.99787596774,
  \]
  and
  \[
  \Sigma_\chi = c_1+\gamma_1-\kappa_U,
  \qquad
  \zeta_Z = 2\lambda_1-\kappa_\eta-\kappa_W.
  \]

So the next theorem gate is now much sharper than “extract the whole 11-component tangent.”
It is:

> does the actual moving-throat branch dynamically generate
> \[
> \Sigma_\eta = 0,
> \qquad
> \Sigma_{\rm sup}=0
> \]
> on the coherent branch-law manifold?

---

## Stage 231 — exact coherent selector bridge

Files:
- `5pn_stage231_coherent_selector_bridge.py`
- `5pn_stage231_coherent_selector_bridge_output.txt`

Using the Stage-230 selector formulas and substituting the actual minimal coherent-continuum observables
\[
 d\ln K = \kappa_\eta,
 \qquad
 d\ln M = 0,
 \qquad
 d\ln \Omega_U = \kappa_U/2,
\]
\[
 d\ln \chi_0 = \Sigma_\chi = c_1+\gamma_1-\kappa_U,
\]
\[
 d\ln Z_W = \zeta_Z = 2\lambda_1-\kappa_\eta-\kappa_W,
\]
\[
 d\ln \epsilon_W = \Sigma_\epsilon = 2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,
\]
one gets the exact collapse
\[
 S_{\rm shape} = -\Sigma_\eta.
\]

This is already a major simplification. The first reduced selector is not some hidden overlap residual. It is exactly the microscopic dressing slippage.

Even better, Stage-166 had already shown
\[
 \mathcal R_1+\Xi_1
 =
 -\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta.
\]
Combining the two gives the exact selector bridge
\[
 \mathcal R_1+\Xi_1
 =
 \frac{\epsilon_\eta}{1-\epsilon_\eta}
 S_{\rm shape}.
\]
So on the physical branch \(0<\epsilon_\eta<1\), the following are equivalent:
\[
 S_{\rm shape}=0,
 \qquad
 \Sigma_\eta=0,
 \qquad
 \mathcal R_1+\Xi_1=0.
\]

This means the Stage-230 shape selector is the overlap-image form of the same microscopic dressing mismatch that already controlled the selected-branch demand law.

---

## Stage 232 — exact support-plane selector on the coherent continuum branch

Files:
- `5pn_stage232_coherent_support_plane_test.py`
- `5pn_stage232_coherent_support_plane_test_output.txt`

The second Stage-230 selector is
\[
 S_{\rm support}
 =
 d\ln \Omega_U - d\ln \Omega_U^{({\rm pred})}(d\ln K,d\ln M,d\ln\chi_0,d\ln Z_W).
\]

On the same minimal coherent-continuum branch, this becomes one exact microscopic scalar:
\[
 S_{\rm support} = \Sigma_{\rm sup}
 =
 \frac{\kappa_U}{2}
 -c_K\kappa_\eta
 +c_\chi\Sigma_\chi
 +c_Z\zeta_Z.
\]

So the support selector depends only on

- wall-stiffness drift \(\kappa_\eta\),
- internal-\(U\) frequency drift \(\kappa_U/2\),
- coherent interference-ratio drift \(\Sigma_\chi\),
- wall-to-mixed overlap drift \(\zeta_Z\).

And it is **independent of** \(d\ln\epsilon_W\) at this reduced level.

That is the second real narrowing of the theorem gap: the support part of the overlap-image test is now a concrete microscopic plane condition instead of an abstract overlap residual.

---

## Stage 233 — microscopic two-selector theorem and selected-spectral corollary

Files:
- `5pn_stage233_microscopic_two_selector_theorem.py`
- `5pn_stage233_microscopic_two_selector_theorem_output.txt`

Combining Stages 231 and 232 with the exact Stage-230 residual vector gives the full reduced microscopic form
\[
 (0,0,S_{\rm support},S_{\rm shape}/8,0,S_{\rm shape},0)
 =
 (0,0,\Sigma_{\rm sup},-\Sigma_\eta/8,0,-\Sigma_\eta,0).
\]

So overlap-image membership on the coherent branch-law manifold is now equivalent to exactly
\[
 \boxed{\Sigma_\eta=0,\qquad \Sigma_{\rm sup}=0.}
\]

This is the sharpest current 5PN microscopic continuation point.

### Exact restoration map

Because the Stage-230 restoration map is
\[
 d\ln\Omega_U \to d\ln\Omega_U-S_{\rm support},
\]
\[
 d\ln\epsilon_W \to d\ln\epsilon_W-S_{\rm shape},
\]
\[
 d\ln\Omega_W \to d\ln\Omega_W-S_{\rm shape}/8,
\]
it becomes, microscopically,
\[
 d\ln\Omega_U \to d\ln\Omega_U-\Sigma_{\rm sup},
\]
\[
 d\ln\epsilon_W \to d\ln\epsilon_W+\Sigma_\eta,
\]
\[
 d\ln\Omega_W \to d\ln\Omega_W+\Sigma_\eta/8.
\]
So the reduced overlap-image restoration has a direct microscopic meaning.

### Raw selected-spectral branch corollary

Applying the same selector formulas to the raw selected-spectral branch from the earlier spectral extractor,
\[
 d\ln K = \sigma_K^{\rm raw},
 \qquad
 d\ln M = 0,
 \qquad
 d\ln\Omega_U = 0,
\]
\[
 d\ln\chi_0 = 0,
 \qquad
 d\ln Z_W = 2\sigma_\kappa^{\rm raw}-\sigma_K^{\rm raw},
 \qquad
 d\ln\epsilon_W = 2\sigma_\kappa^{\rm raw},
\]
one gets the exact selector miss
\[
 S_{\rm shape}^{({\rm raw})} = \sigma_K^{\rm raw},
\]
\[
 S_{\rm support}^{({\rm raw})}
 =
 2c_Z\sigma_\kappa^{\rm raw}-(c_K+c_Z)\sigma_K^{\rm raw}.
\]
Numerically,
\[
 S_{\rm support}^{({\rm raw})}
 \approx
 3.99575193547\,\sigma_\kappa^{\rm raw}
 -3.95028944679\,\sigma_K^{\rm raw}.
\]

So the raw selected branch misses the packet-null image in exactly two ways:

1. a **shape/dressing miss** carried by \(\sigma_K^{\rm raw}\),
2. a **support-plane miss** carried by a fixed linear combination of \(\sigma_K^{\rm raw}\) and \(\sigma_\kappa^{\rm raw}\).

That is the cleanest current diagnostic of what actual moving-throat co-evolution the 5PN branch still has to generate.

---

## Best current continuation point after stages 231–233

The continuation point is now very narrow.

It is no longer:
- extract the whole overlap tangent,
- or solve a generic packet-null image problem,
- or guess how the selected branch should be corrected.

It is now exactly:

> compute the physical coherent-branch slippages that control
> \(\Sigma_\eta\) and \(\Sigma_{\rm sup}\),
> and test whether the actual moving-throat branch drives both to zero.

That is the next honest theorem gate.
