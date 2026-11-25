# Radius 4 · Held relaxation stops at step 167

The ε = 10⁻⁸ continuation from `05a1b666` accepted steps 109–166. All 58 lowered fixed-field energy and preserved positive orientation. Including the two field updates, total energy fell by 0.00154804675549. The smallest peak crystal force was 0.04374951 at step 123. It reached 0.10460517 at step 166 against the 0.000030 target. The strongest remaining force was axial, −0.10460517 at node (24, 3, 20), depth 0.07761939 and χ = 0.11448753. Depth stayed at 0.81914036 and exit radius at 3.30637727. [Accepted work and trial ledger](../../../research/pde_ledger_v3/_measurements/mouth_r4_segmented_stop_work.json)

| Held primary quantity | Step 108, fresh fields | Step 148, fresh fields | Step 166, frozen fields |
|---|---:|---:|---:|
| Total energy | 92185.39080106 | 92185.38944731 | 92185.38925302 |
| Peak active crystal force | 0.09881033 | 0.36399481 | 0.10460517 |
| Order optimality residual | 1.94581 × 10⁻⁹ | 3.99757 × 10⁻⁹ | 0.00126198 |
| First fully melted cross-section, depth | 0.33635069 | 0.33635069 | 0.33635069 |

Step 166 retained the order and carrier fields from step 148. Inward current was 10.70772687, the dense subsonic share was 66.75897%, and feed minus bulk current was +8.95284 × 10⁻¹³. At the fresh updates, all 52,040 free zero-order nodes at step 128 and all 51,993 at step 148 gained energy under +10⁻⁹. Their largest absolute gauge-to-physical work ratios were 8.81434 × 10⁻¹⁰ and 2.13714 × 10⁻¹⁰, respectively. The whole-cell chart lower bound at step 166 was +3.08172600 × 10⁻⁵. The smallest sampled physical inverse determinant was +0.00130119 globally and +0.43442795 in ordered material. [Saved-state field and orientation readback](../../../research/pde_ledger_v3/_measurements/mouth_r4_segmented_stop_local.json)

**Step 167 failed the energy-decrease test.** Nine larger trials exceeded the radial label range and one Newton trial at 1/16 inverted a corner in fully melted cell (25, 10, 7). All 18 smaller trials passed orientation and raised energy. At 1/8192 the increases were +5.26707 × 10⁻⁷ for Newton and +7.51571 × 10⁻⁷ for projected Jacobi. CG encountered nonpositive curvature on every accepted step and again on the unsuccessful attempt. The last accepted state is preserved. [Search fractions, signed work and refusal predicates](../../../research/pde_ledger_v3/_measurements/mouth_r4_segmented_stop_work.json)

The saved path contains 421 recrossings of ξₜ = 1/12 at 91 distinct quadrature-ray endpoints. Of those recrossings, 58 were within one grid index of the instantaneous peak force in each direction. No crystal label node recrossed a plane. Step 166 recorded no crossing. The rejected trials lack saved label maps, so their energy increases cannot yet be attributed to that plane. The open diagnostic is the one-sided work of the entrance force near the plane. [Same-plane crossing measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_segmented_stop_chatter.json)

Independent review concluded, “Verified a numerical step-167 line-search stop on the primary held shape. Step 166 is feasible and above the internal target. The paired stationary readout remains unavailable.” The secondary gauge was never admitted. All six shape coordinates stay held, and radius 3 stays blind. [Independent review](../../../research/pde_ledger_v3/_measurements/mouth_r4_segmented_held_stop_review.json)
