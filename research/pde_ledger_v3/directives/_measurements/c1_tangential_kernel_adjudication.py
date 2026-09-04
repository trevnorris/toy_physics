#!/usr/bin/env python3
"""Adjudicate the ONE disagreement between the two S11c-c1 build-review legs (rule 6/13).

Grok Finding 3: the first-shape DtN kernel "left-quantizes the tangential bilinear form" —
uses k_input**2 (input-leg Laplacian h∇²) not k·k' (Div(h∇), both legs); KERNEL_HAS_K_OUT False
⇒ a rule-17 freeze.
Claude: k_input**2 (η-order) + the tilt cross-term (σ_W-order) reconstruct k·k' exactly under the
gradient relation w1_jet_hat = i(k_out−k_in)·w1_hat ⇒ correct multigrade bookkeeping, NOT a freeze.

The engine's tangential `remainder` (S11c_c1_bulk_closure_sympy_audit.py `dtn_first_kernel`, structural):
    remainder = height*(k_input_sq − kappa_sq) − i*(tilt · k_in)
with height ∝ w1_hat (η-order) and tilt ∝ w1_jet_hat (σ_W-order). Physically w1_jet_hat is the FT of
∇h, i.e. i·(transfer)·w1_hat with transfer = k_out − k_in. Impose that relation and compare to the
true Div(h∇) two-momentum bilinear form ⟨k|Div(h∇)|k'⟩ ∝ (k_out·k_in − kappa²)·w1_hat.

Result (VERDICT): remainder − target ≡ 0 ⇒ the engine reconstructs the TWO-LEG k·k' form; Grok's
"tangential freeze" is a FALSE POSITIVE (its symbol probe missed that k_out enters via the profile
hats' transfer argument, per the spec's independent η/σ_W bookkeeping). The DtN kernel is correct.
"""
import sympy as sp


def main() -> int:
    # 1D reduction, structure only (the eta/2 and sigma_W/2 overall coefficients are irrelevant to structure)
    k_in, k_out, kap, H = sp.symbols("k_in k_out kappa Hhat")
    # gradient relation: FT of d/dx(h) at transfer (k_out - k_in) = i*(k_out - k_in)*Hhat
    jet = sp.I * (k_out - k_in) * H
    # engine's tangential remainder (structural)
    remainder = H * (k_in**2 - kap**2) - sp.I * (jet * k_in)
    # the correct Div(h grad) tangential bilinear form (two-leg k·k')
    target = H * (k_out * k_in - kap**2)
    residual = sp.simplify(remainder - target)
    print("remainder           =", sp.expand(remainder))
    print("target (k.k form)   =", sp.expand(target))
    print("remainder - target  =", residual)
    print("VERDICT:", "TWO-LEG k.k RECONSTRUCTED (Grok freeze = false positive)"
          if residual == 0 else "MISMATCH — investigate")
    return 0 if residual == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
