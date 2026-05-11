# Moving-Throat PDE — Stage 113: Outlet-Model Status Update

## What is now explicit

The abstract isotropic branch-selection triple `(b,a_0,a_5)` is no longer just symbolic. Three explicit outlet classes have now been checked:

1. **Pure geometric Robin core**
   \[
   \Lambda_2^{\rm R}=\Lambda_2^{\rm out}+\rho_R,
   \]
   which shifts the canonical outgoing normalization as
   \[
   \chi_Q^{\rm R}=\frac{3}{3-\rho_R}.
   \]

2. **Standalone mixed `A_w/F_{\mu w}`-type hidden pole**
   \[
   \Lambda_2^{\rm mix}=\Lambda_2^{\rm out}-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5},
   \]
   which cannot preserve the already-fixed even `l=2` branch unless it vanishes.

3. **Hybrid Robin–mixed outlet**
   \[
   \Lambda_2^{\rm hyb}=\Lambda_2^{\rm out}+\rho_R-
   \frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5},
   \]
   which admits one nontrivial compensated canonical-even branch
   \[
   \rho_R=4\sigma_W,
   \qquad
   \kappa_W=\frac13.
   \]
   On that branch,
   \[
   \chi_Q^{\rm hyb}=\frac{1-9\sigma_W\gamma_W}{1-\sigma_W},
   \]
   and canonical outgoing normalization is preserved exactly when
   \[
   \gamma_W=\frac19.
   \]

## PDE-facing consequence

The remaining isotropic branch-selection question is now much sharper than before.
The completed moving-throat PDE does **not** need to decide among an unlimited space of deformations. At low frequency, the first explicit outlet audit says the actual branch must land in one of two categories:

- either a nearly pure scale/argument deformation class, which is harmless,
- or a compensated Robin–mixed outlet of the exact type above.

Pure Robin loading alone and a naive standalone hidden mixed pole are not enough.
