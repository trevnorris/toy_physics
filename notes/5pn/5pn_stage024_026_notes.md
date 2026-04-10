# 5PN / Moving-Throat continuation — Stages 24–26

These stages continue the explicit moving-throat spectral-placement branch after the concrete axial / loaded-profile / selected-prefactor scripts.

## Stage 24 — source map, microscopic normalization, softening depth

Files:
- `5pn_stage24_source_map_microscopic_normalization_softening.py`
- `5pn_stage24_source_map_microscopic_normalization_softening_output.txt`

Main results:

1. The natural D/N source-map factor is no longer independent:
   
   \[
   \hat m_-^2 = \frac{s_-(x)}{\kappa_0^2}.
   \]

2. The selected-branch normalization product becomes
   
   \[
   N_-(x) = \frac{\beta_0\, s_-(x)^2}{\kappa_0^2 (A-x)}.
   \]

3. The selected branch is exactly parameterized by the softening depth
   \[
   x = A-\lambda_-.
   \]

4. The exact secular loading law is
   
   \[
   \alpha_0(x)=\frac{x(x+\Delta K_{ax})}{\kappa_0^2(x+\Delta K_{ax})+\kappa_1^2 x}.
   \]

5. The exact support-loading requirement is
   
   \[
   \frac{g_{B,req}^2}{\varpi^2}=\alpha_0(x)-\alpha_{mix}.
   \]

So the selected quadrupole bridge is no longer an eigenvector problem. It is a one-variable spectral-placement problem in `x`.

---

## Stage 25 — universal D/N branch geometry

Files:
- `5pn_stage25_dimensionless_normalization_and_support_frontier.py`
- `5pn_stage25_dimensionless_normalization_and_support_frontier_output.txt`

Main results:

Using the exact D/N constants
\[
\kappa_0^2=\frac{8}{\pi^2},\qquad
\kappa_1^2=\frac{16}{9\pi^2},\qquad
\eta=\frac{\kappa_1^2}{\kappa_0^2}=\frac29,
\]
and the dimensionless variables
\[
\xi=\frac{x}{A},\qquad \delta=\frac{\Delta K_{ax}}{A},
\]
we get two universal branch functions:

1. The normalization function
   
   \[
   F(\xi,\delta)=
   \frac{(9\delta+11\xi)^4}{81(1-\xi)(9\delta^2+18\delta\xi+11\xi^2)^2}.
   \]

2. The support-feasibility function
   
   \[
   G(\xi,\delta)=\frac{9\xi(\xi+\delta)}{9\delta+11\xi}.
   \]

Exact geometric facts verified in SymPy:

- `F` is strictly increasing on the stable branch,
- `F(0,delta)=1`,
- `F -> +infinity` as `xi -> 1^-`,
- `G` is strictly increasing,
- `G(0,delta)=0`,
- `G_max(delta)=9(1+delta)/(9 delta + 11)`.

So for fixed `delta`, the selected branch is controlled by a unique normalization locus `F = R_target` and a sharp support-feasibility frontier `M_mix <= G`.

---

## Stage 26 — continuum-kernel extraction and placement map

Files:
- `5pn_stage26_continuum_kernel_extraction_and_placement_map.py`
- `5pn_stage26_continuum_kernel_extraction_and_placement_map_output.txt`

Main results:

From the first explicit finite-throat continuum kernel, the reduced branch data are extracted exactly:

\[
A=\frac{K_U K_{\eta}^{eff}-c_{\eta U}^2}{\mu_\eta K_U},
\qquad
\Delta K_{ax}=\frac{\pi^2 T_w}{\mu_\eta L^2},
\]
\[
\alpha_{mix}=
\frac{(K_U c_{\eta W}+c_{UW} c_{\eta U})^2}
{\mu_\eta K_U (K_U K_W^{eff}-c_{UW}^2\sigma)},
\]
\[
\beta_0=
\frac{\mu_W}{\mu_\eta}
\frac{(K_U c_{\eta W}+c_{UW} c_{\eta U})^2}
{(K_U K_W^{eff}-c_{UW}^2\sigma)^2}.
\]

Then the full placement problem compresses to the exact dimensionless kernel ratios

\[
\epsilon_\eta = \frac{c_{\eta U}^2}{K_U K_\eta^{eff}},
\qquad
\epsilon_W = \frac{c_{UW}^2\sigma}{K_U K_W^{eff}},
\]
\[
\rho = \frac{c_{UW} c_{\eta U}}{K_U c_{\eta W}},
\qquad
Z_W = \frac{c_{\eta W}^2}{K_\eta^{eff} K_W^{eff}},
\qquad
\delta_0 = \frac{\pi^2 T_w}{L^2 K_\eta^{eff}},
\]
plus the radiative demand scale
\[
\Lambda = \frac{27 \pi^2 G c_s^5 K_W^{eff}}{20 a^5 c^5 \mu_W}.
\]

The exact placement formulas are

\[
\delta = \frac{\delta_0}{1-\epsilon_\eta},
\qquad
M_{mix} = \frac{8 Z_W (1+\rho)^2}{\pi^2 (1-\epsilon_\eta)(1-\epsilon_W)},
\]
\[
R_{target} = \frac{\Lambda (1-\epsilon_\eta)(1-\epsilon_W)^2}{Z_W (1+\rho)^2}.
\]

And the exact product law is

\[
R_{target} M_{mix} = \frac{8\Lambda(1-\epsilon_W)}{\pi^2}.
\]

This cleanly separates the placement problem into three lanes:

1. geometry lane: `delta = delta_0 / (1-eps_eta)`
2. mixed-stability/product lane: `R_target M_mix = 8 Lambda (1-eps_W)/pi^2`
3. redistribution lane: `(eps_eta, Z_W, rho)` move the branch along the fixed product curve.

---

## What these stages change in the 5PN program

The remaining bridge is no longer a vague “solve more PDE” task.
It is now a very narrow branch-placement problem:

1. compute the actual continuum ratios `(eps_eta, eps_W, rho, Z_W, delta_0, Lambda)`,
2. map them to `(delta, M_mix, R_target)`,
3. find the unique `xi` solving `F(xi,delta)=R_target`,
4. and check the support-feasibility condition `M_mix <= G(xi,delta)`.

If that succeeds, the selected quadrupole branch is admissible on the natural finite-throat continuum placement map.
