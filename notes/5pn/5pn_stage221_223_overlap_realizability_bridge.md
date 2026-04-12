# 5PN continuation notes — stages 221 through 223

## Purpose

These stages take the support-restored packet-null family from Stages 218–220 and attack the next honest problem:

> how does an actual moving-throat overlap branch map into that packet-null family, and does the simplest coherent profile branch land inside it?

This is the realizability step that had remained open after the support-restoration rescue. The answer is now much sharper.

There are three distinct results:

1. an exact **realizability compiler** from overlap-state tangent variables into the Stage-220 packet-null carrier space,
2. a **no-go theorem** for the pure coherent-profile branch,
3. and a unique **compensated realization family** showing how the branch can still survive once the upper-leg and frequency variables are allowed to co-move.

---

## Stage 221 — exact overlap realizability map

Files:
- `5pn_stage221_overlap_realizability_map.py`
- `5pn_stage221_overlap_realizability_map_output.txt`

### Input

Start from the Stage-220 support-restored packet-null carrier variables

- `alpha_K`
- `alpha_GW`
- `alpha_GU`
- `alpha_R`
- `alpha_OU`
- `beta_B`
- `beta_varpi`

and restrict them to the moving-throat overlap dictionary supplied by the coherent finite-throat branch:

- `C = lambda_B * kappa`
- `G_U = lambda_U`
- `G_W = lambda_W * kappa`
- `R = lambda_R * kappa`

with primitive overlap-state drifts

- `sigma_K      = d ln K_geo`
- `sigma_kappa  = d ln kappa`
- `ell_B        = d ln lambda_B`
- `ell_W        = d ln lambda_W`
- `ell_R        = d ln lambda_R`

and compensation variables

- `ell_U        = d ln lambda_U`
- `ell_OmegaU   = d ln Omega_U`
- `ell_varpi    = d ln varpi`.

### Exact carrier map

The overlap-state tangent maps into the Stage-220 carrier space as

- `alpha_K     = sigma_K`
- `alpha_GW    = ell_W + sigma_kappa`
- `alpha_GU    = ell_U`
- `alpha_R     = ell_R + sigma_kappa`
- `alpha_OU    = ell_OmegaU`
- `beta_B      = ell_B + sigma_kappa`
- `beta_varpi  = ell_varpi`.

So the realizability question becomes purely linear:

> does there exist a choice of `(ell_U, ell_OmegaU, ell_varpi)` such that the induced carrier vector lies in the exact Stage-220 packet-null family?

### Exact realization solve

At the constructive point `(E_*,F_*) = (1/4, 5/6)`, solving the full packet-null equations gives

\[
ell_U
\approx
0.84421482069524\,\sigma_K
+0.13857149963136\,\sigma_\kappa
-0.002148228031649\,\ell_B
+1.02617196909895\,\ell_R
-0.885452241435934\,\ell_W,
\]

\[
ell_{\Omega_U}
\approx
0.57648223573430\,\sigma_K
+0.067558848518307\,\sigma_\kappa
-0.013251749605984\,\ell_B
+0.718150303781305\,\ell_R
-0.637339705657013\,\ell_W,
\]

\[
ell_{\varpi}
\approx
-5.45139383054918\,\sigma_K
+0.216979224836926\,\sigma_\kappa
+0.489740172938324\,\ell_B
-4.25085088235995\,\ell_R
+3.97808993425855\,\ell_W.
\]

So the overlap-state tangent is not generically free. But it does admit an exact three-parameter compensation solve.

### Induced orbit/companion drifts

The same solve fixes the dependent drifts required to stay on the Stage-220 packet-null manifold:

\[
d\ln\delta_U
\approx
0.50522670126549\,\sigma_K
-0.0056516769732235\,\sigma_\kappa
-0.039854080113249\,\ell_B
-0.965244046150368\,\ell_R
+0.999446449290394\,\ell_W,
\]

\[
d\ln M
\approx
0.46453483007811\,\sigma_K
-0.14202530222611\,\sigma_\kappa
-0.022207043148670\,\ell_B
-0.616043330635281\,\ell_R
+0.496225071557841\,\ell_W,
\]

\[
d\ln\Omega_W
\approx
-0.27660634196525\,\sigma_K
+0.51753444540572\,\sigma_\kappa
+0.0039179033515622\,\ell_B
+0.0731670124294903\,\ell_R
+0.440449529624671\,\ell_W.
\]

### Meaning

This is the first exact realizability compiler from the moving-throat overlap variables into the support-restored 5PN packet-null family.

So the old vague question

> “does the branch somehow realize the packet-null directions?”

has now become

> “does the actual branch generate the exact compensation ratios above?”

---

## Stage 222 — pure coherent-profile branch no-go

Files:
- `5pn_stage222_pure_profile_branch_nogo.py`
- `5pn_stage222_pure_profile_branch_nogo_output.txt`

### Setup

The cleanest candidate realization is the pure coherent-profile branch coming from the nonconstant finite-throat family already isolated in the moving-throat notes:

- only `kappa(theta)` and `K_geo(theta)` move,
- all couplings `lambda_B, lambda_W, lambda_R, lambda_U` are frozen,
- all frequencies are frozen.

In carrier language this is the tangent

\[
(\alpha_K,\alpha_{GW},\alpha_{GU},\alpha_R,\alpha_{OU},\beta_B,\beta_{\varpi})
=
(\sigma_K,\sigma_\kappa,0,\sigma_\kappa,0,\sigma_\kappa,0).
\]

### Exact result

Solving the full packet-null equations on this restricted branch gives only the trivial solution

\[
\sigma_K = \sigma_\kappa = 0.
\]

So the pure profile branch is an exact no-go for nontrivial packet-null realizability.

### Actual coherent family inserted

The nonconstant coherent family has

\[
\kappa(\theta)=\kappa_0\cos\theta+\kappa_1\sin\theta,
\qquad
\kappa_0=\frac{2\sqrt2}{\pi},
\qquad
\kappa_1=-\frac{4}{3\pi},
\]

and

\[
K_{\rm geo}(\theta)=K_0+T_{\rm grad}\sin^2\theta.
\]

Therefore

\[
\sigma_\kappa(\theta)=\partial_\theta\ln\kappa(\theta)\,d\theta
=
 d\theta\,\frac{3\sqrt2\sin\theta+2\cos\theta}{2\sin\theta-3\sqrt2\cos\theta},
\]

\[
\sigma_K(\theta)=\partial_\theta\ln K_{\rm geo}(\theta)\,d\theta
=
\frac{2T_{\rm grad}\sin(2\theta)}{2K_0-T_{\rm grad}\cos(2\theta)+T_{\rm grad}}\,d\theta.
\]

Two concrete points are especially useful:

- At `theta = 0`:
  \[
  \sigma_\kappa/d\theta = -\frac{\sqrt2}{3},
  \qquad
  \sigma_K/d\theta = 0.
  \]

- At the max-coupling point `theta = theta_max = -arctan(\sqrt2/3)`:
  \[
  \sigma_\kappa/d\theta = 0,
  \qquad
  \sigma_K/d\theta = -\frac{12\sqrt2\,T_{\rm grad}}{22K_0+4T_{\rm grad}}.
  \]

So the actual nonconstant overlap family does produce nontrivial profile drifts, but those drifts do **not** lie inside the packet-null family by themselves.

### Meaning

This is a real obstruction, but not a dead end.

It says only that profile motion alone cannot carry the reduced 5PN closure.
The upper-leg and frequency sector must co-evolve with it.

---

## Stage 223 — unique compensated profile family

Files:
- `5pn_stage223_profile_compensation_family.py`
- `5pn_stage223_profile_compensation_family_output.txt`

### Setup

Now keep the same coherent profile branch, but allow the upper-leg and frequency variables to compensate:

- `d ln lambda_U`
- `d ln Omega_U`
- `d ln varpi`

while still freezing `lambda_B`, `lambda_W`, and `lambda_R`.

This is the smallest compensated realization family worth testing after the pure-profile no-go.

### Exact compensation law

The Stage-221 realization compiler then collapses to

\[
d\ln\lambda_U
\approx
0.84421482069524\,\sigma_K
+0.13857149963136\,\sigma_\kappa,
\]

\[
d\ln\Omega_U
\approx
0.57648223573430\,\sigma_K
+0.067558848518307\,\sigma_\kappa,
\]

\[
d\ln\varpi
\approx
-5.45139383054918\,\sigma_K
+0.216979224836926\,\sigma_\kappa.
\]

So once the coherent profile chooses a tangent `(sigma_K, sigma_kappa)`, there is a **unique** compensation direction in `(lambda_U, Omega_U, varpi)` that restores packet-nullness.

### Induced companion drifts

The required co-moving branch also forces

\[
d\ln\delta_U
\approx
0.50522670126549\,\sigma_K
-0.0056516769732235\,\sigma_\kappa,
\]

\[
d\ln M
\approx
0.46453483007811\,\sigma_K
-0.14202530222611\,\sigma_\kappa,
\]

\[
d\ln\Omega_W
\approx
-0.27660634196525\,\sigma_K
+0.51753444540572\,\sigma_\kappa.
\]

So the realization family is not just “adjust three knobs.” It is a full co-evolving microscopic branch.

### Concrete evaluations

#### At `theta = 0`

Since `sigma_K = 0` and `sigma_kappa = -(sqrt(2)/3) dtheta`, the compensation slopes are

\[
\frac{d\ln\lambda_U}{d\theta}
\approx -0.0653232313790,
\qquad
\frac{d\ln\Omega_U}{d\theta}
\approx -0.0318475466110,
\qquad
\frac{d\ln\varpi}{d\theta}
\approx -0.102284987506.
\]

The induced companion drifts are

\[
\frac{d\ln\delta_U}{d\theta}
\approx 0.00266422607523,
\qquad
\frac{d\ln M}{d\theta}
\approx 0.0669513695361,
\qquad
\frac{d\ln\Omega_W}{d\theta}
\approx -0.243968077229.
\]

#### At `theta = theta_max`

Since `sigma_kappa = 0`, the whole family is driven only by `sigma_K`, so all compensation slopes become proportional to

\[
\frac{T_{\rm grad}}{2K_0 + 0.363636\ldots T_{\rm grad}}.
\]

The main ones are

\[
\frac{d\ln\lambda_U}{d\theta}
\approx
-1.30243641707\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}},
\]

\[
\frac{d\ln\Omega_U}{d\theta}
\approx
-0.889384359537\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}},
\]

\[
\frac{d\ln\varpi}{d\theta}
\approx
+8.41029282436\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}}.
\]

The induced drifts are

\[
\frac{d\ln\delta_U}{d\theta}
\approx
-0.779452857821\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}},
\]

\[
\frac{d\ln M}{d\theta}
\approx
-0.716674316608\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}},
\]

\[
\frac{d\ln\Omega_W}{d\theta}
\approx
+0.42674229845\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}}.
\]

### Meaning

This is the strongest positive result of the session.

The pure profile branch is dead, but the actual coherent overlap family is **not**.
It survives as a unique compensated realization family inside the exact Stage-220 packet-null manifold.

So the remaining question is no longer algebraic and no longer generic. It is now extremely concrete:

> does the true moving-throat PDE branch dynamically generate the compensation ratios required by the Stage-223 family?

---

## Best current summary after Stages 221–223

The realizability problem is now sharply split.

### What is ruled out

- A pure coherent-profile deformation with frozen upper-leg and frequency data cannot realize the 5PN packet-null condition.

### What survives

- The coherent finite-throat overlap family does admit a unique support/frequency compensation direction that restores exact packet-nullness.

### What the PDE must now decide

The completed moving-throat branch has to choose between two possibilities:

1. it does **not** co-evolve in the Stage-223 ratios, in which case the present reduced 5PN closure fails on that branch;
2. it **does** co-evolve in those ratios, in which case the reduced 5PN obstruction survives only at the deeper level of actual branch existence and stability.

So the program is now exactly where it should be:

- the algebraic bottleneck is gone,
- the realization bottleneck is explicit,
- and the next test is a genuine microscopic branch-selection test rather than another symbolic repackaging.
