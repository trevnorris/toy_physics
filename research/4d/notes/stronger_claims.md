## 3) Step A: Conservative path now, roadmap to the stronger claim

### A1) Conservative claim (what you said you want first)

This version is basically what the handoff sheet says Paper 7 needs: one frozen model spec + self-contained derivation + controlled limits + minimal numerics. 

You’d structure the paper’s “math claim” like:

1. **Define the 4D bulk model** by action and conventions.
2. Derive **bulk EOM** (GNLS + Maxwell + force ledger).
3. Define **brane observables** by projection (ρ_brane, J_brane) and derive the projected continuity identity exactly.
4. Present the **Poisson hook** as a *regime/limit*, i.e. a Helmholtz decomposition of v_brane and an exact identity that becomes Poisson-like under stated approximations. 
5. Present the **EM localization reduction** as a *controlled zero-mode reduction* giving a 3+1D Maxwell equation with μ0eff.
6. Present **leakage** as an explicit source term in the projected continuity, with symmetry gates (when it vanishes, when it turns on).

That’s fully defensible and doesn’t require you to claim a closed 3D effective theory.

### A2) Roadmap to the stronger claim (explicit reduced 3+1D effective theory)

To make the stronger claim, you need to “freeze” and quantify the reduction — i.e. not just show identities, but produce a bona fide reduced 3+1D system with controlled correction terms:

1. **Mode expansion + small parameter**
   Adopt a truncation like (\psi \approx \psi_0(x,t)\chi_0(w) + \varepsilon,\psi_1(x,t)\chi_1(w) + \dots), and state what “small” means (e.g. small excitation of higher w-modes outside the throat). This is already aligned with your “hard-mode” direction: 4D is primary, brane is emergent. 

2. **Projected GNLS for ψ0, with explicit overlap constants**
   You already have the separable HO-w projected GNLS and C10 overlap constant in the appendix output. The stronger version adds:

   * the ψ1 backreaction terms (mode coupling),
   * the correction to the effective nonlinearity if ψ1 is not negligible.

3. **Projected momentum equation + “Reynolds stress” correction** (high referee value)
   This is the missing piece if you want Poisson emergence to look less “hand-wavy”.
   The SymPy derivation tool already anticipates this and prints a brane-projected momentum balance with extra stress sectors under a `--show_brane_momentum` style flag.
   In the paper, you’d define
   [
   \Pi_{ij}^{\rm brane}=\int W,\rho,v_i v_j,dw,\qquad
   R_{ij}=\Pi_{ij}^{\rm brane}-\rho_{\rm brane},v_i^{\rm brane}v_j^{\rm brane},
   ]
   and show where it enters the reduced momentum equation.

4. **Localized Maxwell KK tower estimate (optional but the real “strong claim” backbone)**
   The reduction you have is correct *as a zero-mode truncation*. To claim “the brane EM is 3D-like” robustly, you’d add the standard argument:

   * the zero mode is normalizable because (\int Z(w),dw<\infty),
   * higher modes are massive in 3+1D and give Yukawa-suppressed corrections on-brane.
     This can be done analytically (Sturm–Liouville) with Z(w)=Gaussian, plus a short scaling estimate.

5. **Leakage sector beyond parity-symmetric baseline**
   You’ll want one explicit knob that “turns leakage on” and show leading-order behavior. v5 now includes exactly that (shifted weight example).
