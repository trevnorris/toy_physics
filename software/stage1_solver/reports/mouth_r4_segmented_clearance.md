# Radius 4 · Held-shape restart from step 108

The audit uses the two independently equilibrated step-108 field states from `cc11a9f3` and the two crossing locations identified in `974127ee`. The entrance radius is 4, depth is 0.81914036 and exit radius is 3.30637727. All six shape coordinates are held.

| Step-108 measurement | ε = 10⁻⁸ | ε = 10⁻¹⁰ |
|---|---:|---:|
| Full signed depth-probe energy change | −6.12972676753 × 10⁻⁷ | −6.12972872025 × 10⁻⁷ |
| Integrated analytic work minus energy change, divided by full probe width | +7.61525 × 10⁻¹⁰ | +1.05690 × 10⁻⁹ |
| Largest local analytic slope discrepancy across terms, cells and intervals | 3.69395 × 10⁻⁸ | 3.71074 × 10⁻⁸ |
| Peak active crystal force | 0.09881032648 | 0.09881032561 |
| Inward current | 10.707728693785 | 10.707728693789 |
| Feed minus bulk current | +4.86722 × 10⁻¹³ | +3.89022 × 10⁻¹³ |

Each of the four affected cell rays crosses one collar reference plane. The two distinct crossing offsets split the full probe into three smooth intervals. Each interval passed its own derivative and signed work checks against the 10⁻⁷ force tolerance. An independent scan of 6,290 active cells detected four changed rays on the negative wide probe. Both signs of the 2⁻¹⁹ and 2⁻²⁰ probes retained each ray’s sequence of reference triangles. Endpoint shape-force measurements use this geometric test to choose the two signed depth probes before measuring their forces. [Segmented work](../../../research/pde_ledger_v3/_measurements/mouth_r4_segmented_depth.json) · [Reference-face scan](../../../research/pde_ledger_v3/_measurements/mouth_r4_endpoint_reference_faces.json)

The physical crystal evaluator and the whole-cell chart now call the same subdivision routine. Both returned a positive lower bound of +6.474133271946422 × 10⁻⁸ for the close-approach trial in melted cell (25, 11, 10), including the roundoff allowance. On the inverted control, both witnessed −4.14306920 × 10⁻⁶ on the edge. The 39-point gauge separately rejected its moving sample at −7.18900962 × 10⁻⁶. The shared certificate changed the energies, physical crystal label gradients, crystal loads on order, peak active force, five native field residuals and ten flow/balance readouts by exactly zero on both saved states. [Positive, negative and invariance controls](../../../research/pde_ledger_v3/_measurements/mouth_r4_shared_orientation.json)

Independent review concluded, “Cleared for held step 109 on the two Gate-1 states. The segmented depth work and shared whole-cell orientation certificate pass both gauges.” The primary Newton–CG continuation starts at the freshly matched step-108 state and checks zero-order-node work after each field re-equilibration. The secondary starts from the certified primary endpoint and must meet the same 0.000030 internal residual target on its own functional. The paired stationary forces and eligibility for shape release remain open. [Independent review](../../../research/pde_ledger_v3/_measurements/mouth_r4_segmented_gate_review.json)
