# 5PN continuation notes — Stages 326–328

These stages perform the extraction that Stage 323–325 was still waiting for.

Before this batch, the coherent branch was already split into

- an **orbit-lock packet** built from
  \[(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda),\]
- and a separate **support/normalization packet** that adds only
  \[\zeta.\]

What was still missing was the exact bridge from the **microscopic kernel state and its grouped weak-axisymmetric drifts** into those actual branch variables.

These stages close that gap.

---

## Stage 326 — exact microscopic coherent placement-state extractor

Files:
- `fivepn_stage326_328_common.py`
- `5pn_stage326_microscopic_coherent_placement_state.py`
- `5pn_stage326_microscopic_coherent_placement_state_output.txt`

The actual coherent placement state is now extracted directly from the microscopic kernel state
\[(\lambda_W,c_{\eta U},\gamma,K_U,K_\eta^{(\mathrm{eff})},K_W^{(\mathrm{eff})},\mu_W,T_U,\lambda_\phi,K_\phi^{(\mathrm{eff})}).\]

The exact formulas are

\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U},
\qquad
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\]
\[
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}},
\qquad
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
\qquad
\Lambda=\Lambda_0\frac{K_W^{(\mathrm{eff})}}{\mu_W},
\]
\[
\zeta=\frac{\lambda_\phi^2K_W^{(\mathrm{eff})}}{\lambda_W^2K_\phi^{(\mathrm{eff})}},
\qquad
\Lambda_0=\frac{27\pi^2Gc_s^5}{20a^5c^5}.
\]

The older reduced Stage-318 variable is recovered exactly as
\[
\widehat Z_W = Z_W\frac{\Lambda_0}{\Lambda}
=\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}.
\]
So the earlier reduced state was not wrong; it was just hiding the quotient structure of two actual branch observables.

The exact coherent observables then become
\[
R_{\rm tr}=\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\]
\[
\epsilon=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
\[
R_{\rm target}=\Lambda\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
\]
The support packet also compiles directly in microscopic variables through
\[
M_{\rm mix}=\frac{8Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
\qquad
M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon),
\]
with
\[
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\qquad
R_{\rm target}M_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]

So after Stage 326, the actual coherent branch is no longer expressed by an abstract reduced packet. It is a direct quotient state of the microscopic kernel.

---

## Stage 327 — exact microscopic coherent placement-drift extractor

Files:
- `5pn_stage327_microscopic_coherent_placement_drifts.py`
- `5pn_stage327_microscopic_coherent_placement_drifts_output.txt`

The grouped weak-axisymmetric drifts of the actual placement variables are now explicit.
Writing the microscopic drift variables as
\[(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1,\phi_1,\kappa_\phi),\]
we get

\[
\delta\ln\chi_0=\gamma_1+c_1-\kappa_U,
\qquad
\delta\ln\delta_U=\tau_1-\kappa_U,
\]
\[
\delta\ln Z_W=2\lambda_1-\kappa_\eta-\kappa_W,
\qquad
\delta\ln\epsilon_W=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,
\]
\[
\delta\ln\epsilon_\eta=2c_1-\kappa_U-\kappa_\eta,
\qquad
\delta\ln\Lambda=\kappa_W-\mu_1,
\]
\[
\delta\ln\zeta=2\phi_1-2\lambda_1+\kappa_W-\kappa_\phi.
\]

The older reduced drift is recovered exactly:
\[
\delta\ln\widehat Z_W
=
\delta\ln Z_W-
\delta\ln\Lambda
=
2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W.
\]

This stage also makes the packet split completely sharp at drift level:

- the **orbit side** uses only
  \[(\delta\ln\chi_0,\delta\ln\delta_U,\delta\ln Z_W,\delta\ln\epsilon_W,\delta\ln\epsilon_\eta,\delta\ln\Lambda),\]
- the **support side** adds only
  \[\delta\ln\zeta.\]

So support transport cannot hide inside the orbit packet even infinitesimally.

---

## Stage 328 — exact microscopic actual-branch packet compiler

Files:
- `5pn_stage328_microscopic_actual_branch_packet_compiler.py`
- `5pn_stage328_microscopic_actual_branch_packet_compiler_output.txt`

This stage closes the loop between the older microscopic monomial theorem and the newer actual coherent placement-map packet.

### 1. Finite orbit packet

The finite quotient packet
\[(q_{\rm tr},q_{\rm nt},q_\eta)\]
compiled from the actual placement state is exactly the finite log-ratio packet of the three direct microscopic monomials. There is no mismatch between the two languages.

### 2. Infinitesimal orbit defect packet

Feeding the actual placement drifts into the Stage-324 packet compiler gives
\[
\Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta,
\]
and from them
\[
\Theta_1,
\qquad
\Xi_1,
\qquad
\mathcal R_1.
\]
The exact theorem proved in the script is that these are **identical** to the older Stage-313/314 microscopic compatibility ledger.

So the actual coherent placement map introduces no hidden extra weak-axisymmetric obstruction. The packet is the same packet in two equivalent coordinate systems.

### 3. Direct observable zero-defect form

On the actual microscopic branch the first-order orbit-lock theorem is exactly
\[
\delta\ln R_{\rm tr}=0,
\qquad
\delta\ln R_{\rm target}=0,
\qquad
\delta\ln\epsilon_\eta=0.
\]

### 4. Final two-packet split

The actual coherent moving-throat branch now ends at the exact microscopic two-packet decomposition

- **orbit packet**
  \[(q_{\rm tr},q_{\rm nt},q_\eta)\]
  or equivalently
  \[(\Theta_1,\Xi_1,\mathcal R_1),\]
- **support packet**
  \[(\zeta;M_{\rm mix},S,M_{\rm tr}).\]

The orbit packet is blind to `lambda_phi`, `K_phi^(eff)`, `phi_1`, and `kappa_phi`; the support packet carries them all.

---

## Where this leaves the 5PN / moving-throat program

The extraction step is now done in the sharpest exact form we have had so far.

1. The actual coherent branch state is a **7-scalar microscopic quotient state**:
   \[(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda,\zeta).\]
2. The weak-axisymmetric orbit-lock theorem depends only on the first six.
3. The support/normalization theorem adds only `\zeta`.
4. The actual orbit packet compiled from those six scalars is **exactly the same** packet as the older microscopic monomial-rigidity ledger.

So there is no further hidden algebraic bottleneck between the microscopic kernel variables and the actual coherent branch observables.

The remaining honest theorem gate is now purely PDE-side:

> compute the actual finite state and grouped weak-axisymmetric drifts of
> \[(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda,\zeta)\]
> on the completed moving-throat operator and test them against the exact orbit packet and support packet above.
