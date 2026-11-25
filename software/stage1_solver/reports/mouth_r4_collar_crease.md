# Radius 4 · A collar reference plane reverses the step-167 energy slope

On the frozen step-166 state from `ec9fa56b`, the reconstructed projected-Jacobi and Newton directions gave initial energy slopes of −0.001547700680 and −0.001357686533. At the smallest archived trial, α = 1/8192, they reproduced the two uphill energy changes exactly, +7.51571281322 × 10⁻⁷ and +5.26706887501 × 10⁻⁷. The following attribution is for projected Jacobi. [Reconstruction and signed work](../../../research/pde_ledger_v3/_measurements/mouth_r4_crease_directions.json)

The first crossing occurs at α★ = 6.63986080497945 × 10⁻⁵. Cells (23, 2, 17) and (24, 2, 17) share the endpoint that crosses ξₜ = 1/12. Its physical position is (r, w, θ) = (3.94593683, 0.06468283, 1.27275310) with χ = 0.15796927. Its depth differs from the peak-force node by 0.01293657, and their spatial separation is 1.17071333. This is the only axial crossing through 4α★. [Endpoint and internal-split census](../../../research/pde_ledger_v3/_measurements/mouth_r4_collar_crossings.json)

| α / α★ | Distinct crossings | ΔE | ΔE / α | Relative secant error against initial slope |
|---:|---:|---:|---:|---:|
| 0.001 | 0 | −1.02765 × 10⁻¹⁰ | −0.00154770062 | 4.05 × 10⁻⁸ |
| 0.01 | 0 | −1.02765 × 10⁻⁹ | −0.00154770009 | 3.79 × 10⁻⁷ |
| 0.1 | 0 | −1.02765 × 10⁻⁸ | −0.00154769482 | 3.78 × 10⁻⁶ |
| 0.5 | 0 | −5.13816 × 10⁻⁸ | −0.00154767141 | 1.89 × 10⁻⁵ |
| 0.9 | 0 | −9.24855 × 10⁻⁸ | −0.00154764801 | 3.40 × 10⁻⁵ |
| 1.1 | 1 | +1.91357 × 10⁻⁹ | +0.00002619942 | 1.01693 |
| α = 1/8192 | 1 | +7.51571 × 10⁻⁷ | +0.00615687193 | 4.97808 |
| 2 | 1 | +9.10276 × 10⁻⁷ | +0.00685462950 | 5.42891 |
| 4 | 1 | +2.73099 × 10⁻⁶ | +0.01028254316 | 7.64375 |

The one-sided slopes at 0.99999α★ and 1.00001α★ are −0.00154758370 and +0.01582285317. Accumulated work is still negative immediately after the crossing and becomes positive by 1.1α★. Cell (23, 2, 17) supplies +9.29298411 × 10⁻⁷ at 1/8192, more than the net increase. Its slope changes sign as the adjacent reference triangle enters the radial integration. Independent finite differences on the two incident rays reproduce the global slope jump within 5.84 × 10⁻⁹. [Signed ladder and cell rankings](../../../research/pde_ledger_v3/_measurements/mouth_r4_crease_ladder.json) · [Independent two-ray check](../../../research/pde_ledger_v3/_measurements/mouth_r4_crease_review.json)

Independent review concluded, “Part 1 clears Branch 1 for the frozen step-166 projected-Jacobi direction.” The open test is whether a smooth reference preserves the registered boundaries and passes the derivative, orientation and field-rematching checks at both gauge strengths. All six shape coordinates were held throughout.
