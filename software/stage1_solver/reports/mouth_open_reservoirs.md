# How much wave energy returns to a draining mouth from the brane and bulk?

**Most of the carrier wave launched by the bend escapes, and the two trapped-wave excitations remain distinguishable in the one-arm calculation.** A round throat connects a laterally extended brane slab to a four-dimensional bulk half-space. Outside the opening the slab faces and bulk plane reflect the carrier. The primary test uses a slab one mouth radius thick, the same carrier equation on both sides, stationary inward mass current 0.0373999726942 and peak bend displacement 0.00001. These are geometry and material assumptions. [Open-reservoir solver](../tools/mouth_open_reservoirs.py) · [Earlier guide comparison](mouth_carrier_bend_work.md)

| Drive angular frequency | Carrier energy lost per cycle | Brane share of escaped power | Bulk share |
|---:|---:|---:|---:|
| 0.456605832128 | 7.47883001722 × 10⁻¹¹ | 42.36% | 57.64% |
| 3.93100103488 | 5.43092854376 × 10⁻¹⁰ | 48.12% | 51.88% |

These are 128-cell-per-radius measurements in model units. Moving the absorbing region farther out changed both powers by less than 0.0005% on the 64-cell grid. Against an independent continuum solution at zero intake, the 128-cell powers were 0.31% and 2.13% low. Opening the reservoirs reduced the zero-intake continuum radiation by 4.11% and 4.27% relative to the previous guides. Slab thickness matters for where the energy goes. At the lower tone, the brane received 27.74% with a half-radius thickness and 45.77% with a two-radius thickness on the 64-cell grid. [Calculations and independent review](../../../research/pde_ledger_v3/_measurements/mouth_open_reservoirs.json)

**Return to the neck.** For the mixture produced by the bend, the brane returned **7.09%** of the traveling power incident on it at the lower tone and **5.62%** at the upper tone. The bulk returned **5.65%** and **3.42%**, respectively. These 96-cell readings use local straight-neck modes and account separately for energy transmitted by fields that decay away from the aperture. A separate launch of only the lowest radial dipole at the upper tone gave **64.39%** return from the brane. A single reflection coefficient misses this dependence on the incident wave. [Driven-wave decomposition](../tools/mouth_driven_return.py) · [Isolated-junction launches](../tools/mouth_aperture_return.py)

**Spectral consequence.** With trial bend mass 12 and static stiffness 6, the flowing three-frequency calculation gave a linewidth of 0.000934278953269 against a spacing of 0.00406283876392. The independent zero-intake continuum width was 0.000946966344636. Both retained two peaks with the earlier equal-area selective readout, one cold loss channel, a common shift and zero additional phase noise. The next material input is the inertia and stiffness of this particular bend at fixed intake. The loss calculated here concerns a longitudinal carrier wave. The throat’s electric sign and 0PN free acceleration remain open.

Working base `df13acaa`. Local review is recorded. Formal review and commit are pending.
