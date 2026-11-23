# Does the ordered radius-4 mouth support itself?

**The mouth drains through an ordered interior. Its core and collar can reach internal elastic balance with the conversion layer responding to both reservoirs. The exterior shape still carries a measurable load.** This test used one quantum of circulation through both baths, a stationary supersolid lattice and an entrance radius of 4. The sealed lateral face and reversible melting exit follow the existing preparation. No occupied breathing wave, interior melting requirement or additional material coefficient entered the calculation.

## What moved

The previous registered mouth balanced five shape coordinates at depth 0.81583 and waist/exit radius 3.57310. We released both components of material displacement across each collar together with the same five coordinates and meridional order. The coarse sharp-exit solve reached depth **0.760626** and waist/exit **3.471100**. Its largest shape-force residual was **0.00009654** and its largest free lattice-node residual was **0.00005639**. On the refined grid the last accepted shape had depth **0.793316** and waist/exit **3.494726**. Its corresponding residuals were **0.07938** and **0.10064**. A subsequent trial stopped in the receiving-order line search at residual 2.36 × 10⁻⁷. These lengths describe the restricted calculations. They do not give a mesh-converged equilibrium geometry. [Shape measurements](../../../research/pde_ledger_v3/_measurements/mouth_joint_collar_release.json)

A displacement passing continuously through the core–collar join still lowered the energy. It needed no bond breaking or slip. We released the shared nodes and both displacement components throughout the core and inner collar at the last refined shape, retaining its external boundaries. We then relaxed the material, terminal order, carrier density and phase to a common stationary point. The crystal receives the interpolated order at its material integration points. The front receives the elastic work deposited back onto those same nodes. Independent signed work showed that sampling elastic energy on a separate grid introduced an error of 0.28283 in a displacement derivative. The shared integration reduced that discrepancy below 10⁻¹² in the same test. [Bonded displacement](mouth_internal_collar.md) · [Independent work checks](mouth_bonded_transfer_review.md)

## What the held mouth transmits

| Measurement on the refined shape, in solver units | Result |
|---|---:|
| Maximum free material-node residual, 33 × 32 solid mesh | 0.00000507 |
| Carrier and conversion-layer stationarity residual | 0.000000182 |
| Sampled principal stretches in the core and inner collar | 0.33272–1.74982 |
| Sampled volume ratio after material relaxation | 0.14554–1.54670 |
| Inward current, one quantum / same lattice without winding | 11.64367 / 13.81320 |
| Intake retained / incoming-layer current denser than half the bath and slower than sound including swirl | 84.29% / 68.96% |
| Drive-induced axial stress integrated on the two sides of the exit disk | +1.08657 inside / +0.44241 receiving |

The 24 × 24 solid mesh gave inward current 11.64370, material-force residual 0.00001086 and carrier/front residual 0.000000189 on the same held geometry. The middle of the bonded join moved inward by 0.00704 and toward the entrance by 0.01324. All sampled material volume ratios remained positive. The fully melted receiving fluid has zero crystalline stiffness, so the model assigns no unique solid stretch to its material labels. The wound current at the two far reservoirs agreed within 4.7 × 10⁻¹³. At the disk itself the directly reconstructed dense, subsonic shares were 67.98% and 75.42% on its two sides. The two stress readings include carrier and order stresses with the resting reference subtracted. Their spread limits the precision of a single transmitted traction. [Reservoir and disk comparisons](mouth_supercore_joined_front.md)

## Why translation stays gated

Releasing the internal solid changes the force on the exterior shape. At the held fine state, the five signed generalized forces were **(+1.4413, −1.6862, −11.2875, −19.1078, −11.8797)** for log depth, scaled log exit radius, middle width, width asymmetry and the localized waist deformation. Their units follow the shape coordinates. A joint step of 0.00015 in their normalized downhill direction lowered energy by **0.00377214**. The coarser solid mesh gave −19.0010 in the same width-asymmetry test and −11.5220 in the localized waist test. Both meshes favor reshaping the waist and shifting where the throat is widest.

Circulation also produces unequal pressure around a single cross section. A smooth wall displacement that widens the polar regions and narrows the equator, preserving average radius to first order and leaving the exit fixed, gave force **+0.08406**. The identical unwound control gave **−0.00000021**. A separately implemented test with a different depth window found the same sign. A polar distortion is therefore an available way to lower the energy. A circular meridian excludes that motion. The opening stayed finite during the computations. Neither a self-supporting shape nor closure or material failure has been established. The next calculation must allow the polar shape and solid displacement to vary together, retain the mobile bonded core and work-conserving front, and check the forces on both sides of the exit. The radius-3 prediction and free translation remain reserved. [Independent mechanical review](mouth_collar_continuum_review.md) · [0PN roadmap](../../../docs/zero_pn_free_fall_roadmap.md)
