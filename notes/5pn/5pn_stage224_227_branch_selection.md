# 5PN continuation notes — stages 224 through 227

These stages continue directly from the Stage-221/223 overlap realizability work.
The target was to answer the next honest question after the compensated coherent-profile family existed algebraically:

> do the **actual moving-throat coherent branch relations** force that compensation family,
> kill it, or leave it underdetermined?

The answer is now sharp:

- the current moving-throat coherent branch laws do **not** kill the compensated packet-null family,
- they do **not** force it either,
- they leave it **underdetermined by exactly four scalar conditions**.

So this is **not** a dead end. It is a precise realizability bottleneck.

---

## Stage 224 — exact overlap-image parameterization

Files:
- `5pn_stage224_overlap_image_parameterization.py`
- `5pn_stage224_overlap_image_parameterization_output.txt`

This stage takes the constructive support-restored packet-null family and rewrites it completely in the overlap-observable drift variables

\[
(d\ln K,\ d\ln M,\ d\ln C,\ d\ln\varpi,\ d\ln\Omega_U,\ d\ln\Omega_W,
\ d\ln\chi_0,\ d\ln\epsilon_\eta,\ d\ln Z_W,\ d\ln\epsilon_W,\ d\ln\delta_U).
\]

### Main result

Although the compensated family begins from five primitive overlap-side drifts,

\[
(\sigma_K,\ \sigma_\kappa,\ \ell_B,\ \ell_W,\ \ell_R),
\]

its observable overlap image has exact rank `4`.

A convenient exact coordinate system on that image is

\[
(d\ln\chi_0,\ d\ln Z_W,\ d\ln C,\ d\ln\varpi).
\]

So the packet-null family already splits naturally into

- a **2D orbit/shape layer**
  \[
  (d\ln\chi_0,\ d\ln Z_W),
  \]
- a **2D support layer**
  \[
  (d\ln C,\ d\ln\varpi).
  \]

The exact constructive image also gives the simple universal formulas

\[
d\ln\delta_U = -\frac{18}{11} d\ln\chi_0,
\]
\[
d\ln\Omega_W = \frac{41}{44} d\ln\chi_0 + \frac58 d\ln Z_W,
\]
\[
d\ln\epsilon_W = 2 d\ln\chi_0 + d\ln Z_W,
\]
\[
d\ln\epsilon_\eta = 0.
\]

The nontrivial remaining image coordinates

\[
d\ln K,
\qquad
d\ln M,
\qquad
d\ln\Omega_U,
\]

are exact linear functions of
\((d\ln\chi_0,d\ln Z_W,d\ln C,d\ln\varpi)\).
The script keeps them symbolically exact and prints useful numerical approximations.

---

## Stage 225 — exact similarity-orbit bridge in overlap variables

Files:
- `5pn_stage225_similarity_orbit_bridge.py`
- `5pn_stage225_similarity_orbit_bridge_output.txt`

This stage takes the direct coherent moving-throat branch laws from the compact PDE program and rewrites them in overlap variables.
The three exact coherent zero-defect residuals are

\[
\Sigma_\eta = d\ln\epsilon_\eta,
\]
\[
\Sigma_{\rm tr} = (1+\chi_0)\,d\ln\delta_U + (1+\delta_U)\,d\ln\chi_0,
\]
\[
\Sigma_{\rm nt} = (d\ln Z_W - 2 d\ln\Omega_W) + E_* d\ln\epsilon_W - F_* d\ln\delta_U.
\]

### Main result

On the constructive overlap image,

\[
\Sigma_\eta = \Sigma_{\rm tr} = \Sigma_{\rm nt} = 0
\]

identically.

So the full compensated overlap image already lies tangent to the exact monomial-preserving similarity orbit of the coherent moving-throat branch.

Just as important, the orbit verdict is **support-blind**:

\[
\partial_{d\ln C}\Sigma_\eta
=
\partial_{d\ln C}\Sigma_{\rm tr}
=
\partial_{d\ln C}\Sigma_{\rm nt}
=0,
\]
\[
\partial_{d\ln\varpi}\Sigma_\eta
=
\partial_{d\ln\varpi}\Sigma_{\rm tr}
=
\partial_{d\ln\varpi}\Sigma_{\rm nt}
=0.
\]

So the current moving-throat coherent/invariant laws only see the orbit/shape side, not the explicit support pair.

---

## Stage 226 — branch-selection theorem

Files:
- `5pn_stage226_branch_selection_theorem.py`
- `5pn_stage226_branch_selection_theorem_output.txt`

This is the actual answer to the fork question.

The full constructive overlap image is encoded by seven exact residuals

\[
(R_K,\ R_M,\ R_{\Omega_U},\ R_{\Omega_W},\ R_{\epsilon_\eta},\ R_{\epsilon_W},\ R_{\delta_U}).
\]

The actual moving-throat coherent branch laws are only the three monomial-preservation residuals

\[
(B_\eta,\ B_{\rm tr},\ B_{\rm nt})
=
(\Sigma_\eta,\Sigma_{\rm tr},\Sigma_{\rm nt}).
\]

### Exact equivalence proved by the script

The full overlap-image residual system is exactly equivalent to

1. the three actual branch-law residuals
   \[
   (B_\eta,B_{\rm tr},B_{\rm nt})=0,
   \]
2. plus one exact **selector quartet**
   \[
   (R_K,\ R_M,\ R_{\Omega_U},\ R_{\epsilon_W})=0.
   \]

Once those four selectors vanish, the remaining two image residuals vanish automatically:

\[
R_{\delta_U}=0,
\qquad
R_{\Omega_W}=0.
\]

### Dimension count

The script computes exact ranks:

- full overlap-state space: dimension `11`,
- actual coherent branch-law system: rank `3`, so branch-law manifold dimension `8`,
- constructive overlap image: rank `7`, so image dimension `4`.

Therefore the codimension gap is

\[
7-3 = 4.
\]

### Final verdict

The current moving-throat coherent branch laws

- do **not** forbid the compensation family,
- do **not** force it,
- and leave it **underdetermined by exactly four scalar conditions**.

That is the cleanest status statement we have had so far.

---

## Stage 227 — overlap realizability tester

Files:
- `5pn_stage227_overlap_realizability_tester.py`
- `5pn_stage227_overlap_realizability_tester_output.txt`

This stage packages the full overlap-image condition as a practical finite-dimensional tester.

### Input
An observed first-order overlap tangent

\[
(d\ln K,\ d\ln M,\ d\ln C,\ d\ln\varpi,\ d\ln\Omega_U,\ d\ln\Omega_W,
\ d\ln\chi_0,\ d\ln\epsilon_\eta,\ d\ln Z_W,\ d\ln\epsilon_W,\ d\ln\delta_U).
\]

### Output
The exact seven residuals

\[
(R_K,\ R_M,\ R_{\Omega_U},\ R_{\Omega_W},\ R_{\epsilon_\eta},\ R_{\epsilon_W},\ R_{\delta_U}).
\]

and the quadratic realizability distance

\[
D_{\rm real}^2 = \sum_i R_i^2.
\]

So once we have an actual moving-throat overlap extraction from a more expensive PDE or branch computation, there is now a direct yes/no test against the current constructive packet-null family.

---

## What changed in the roadmap

Before Stage 224–227, the continuation question was still phrased loosely as

> “Does the actual moving-throat branch generate the compensation ratios?”

Now it is much sharper.

The actual moving-throat coherent branch already guarantees the three monomial-preserving relations.
That is **not** the remaining bottleneck.

The remaining bottleneck is the exact selector quartet

\[
(R_K,\ R_M,\ R_{\Omega_U},\ R_{\epsilon_W}).
\]

So the next honest theorem gate is:

> extract the actual overlap tangent of the moving-throat branch and test whether those four remaining selector residuals vanish.

That is smaller and cleaner than the old “derive 5PN somehow” problem.
