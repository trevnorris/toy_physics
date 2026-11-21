# Can the resting mouth determine its trapped-wave energy and wall tension?

**One trapped ground wave exerts outward pressure throughout the tested flare. Its pressure and a uniform wall tension cannot hold the archived shape in equilibrium, even allowing an arbitrary uniform backing pressure.** This tests an isolated mouth of radius 1 and depth scale 1 draining at 0.25, before applying a companion flow. Spherical sections in the three brane coordinates open out along the fourth spatial coordinate. The carrier solve prescribes this wall shape.

**Measurement.** In each of eight disjoint depth bands we tested the work of a smooth outward displacement. We measured the draining carrier pressure from wall stress and independently from an interior energy variation at fixed intake. Widening the wave's confining well locally lowers its mode frequency. At fixed action $J=E/\omega$, the wave pressure per unit action is $p_1=-\delta\omega/\delta V$. The measured values ranged from 0.0206545 to 0.0326838 on the finest grid. A direct change of the mode frequency checked the integrated wave force. $K=\delta A/\delta V$ is the change of wall area per change in enclosed volume for the same displacement. The proposed equilibrium requires $p_{\rm carrier}+Jp_1=P_{\rm back}+\sigma K$ in every band. [Calculation](../tools/mouth_occupied_radial_balance.py) · [Wave and confinement assumptions](mouth_wave_flare_work.md)

| Carrier mesh | Smallest possible worst pressure error, from wall stress | From interior work |
|---:|---:|---:|
| 26 | 0.000200946 | 0.000200839 |
| 34 | 0.000200963 | 0.000201015 |
| 48 | 0.000201068 | 0.000201067 |

Fitting three preselected bands exactly missed five unused bands. On mesh 48 the signed errors ran from −0.00027330 to +0.00035565 in model pressure units. The table gives the best possible fit using all eight. An independent calculation showed why adjusting the constants cannot remove the discrepancy. The band at depth 1.44–1.50 has the same mean curvature and unit-action wave pressure as a weighted combination of three other bands. Its measured carrier pressure is lower by **0.000402135**. Uniform backing, uniform tension and wave action each cancel in this comparison. [Guarded measurements and independent checks](../../../research/pde_ledger_v3/_measurements/mouth_occupied_radial_balance.json)

**Decision for 0PN.** The saved carrier field fixes inward flow. The homogeneous wave equation fixes frequency and mode shape and admits any amplitude. Neither field provides the stored energy or the total sideways wall force. A physical mouth could settle into a different flare or carry wall or backing stress that varies with depth. The test leaves both possibilities open. The earlier companion accelerations still describe specified preparations. [Occupied-mouth accelerations](mouth_occupied_companion.md)

**Next measurement.** For a mouth prepared independently of its gravity response, record the absolute wave energy and frequency, giving $J$. At the same intake, hold the entrance and impose a small relative sideways offset at each bend, measuring the total restoring reactions. The response at the other bend tests how far the stress travels. This requires a material model or independent reading of the wall reaction and a physical preparation of the wave. Prescribing the missing force at each depth to maintain the archived shape would fit this test by construction.

If the steady drain is what maintains the wave, the input for that preparation is the work the flow does on the wave per cycle and the energy the mode loses per cycle. A stable crossing of these two rates as stored energy varies would select an occupation. The present time-independent confining well supplies no such selection.

Working base `1dc45d09`. One independent reviewer reproduced the pressure bound. A second integrated the wave pressure directly and found the same discrepancy. “Scoped clear for the radial-work mismatch in the prescribed cubic mouth with one scalar ground mode, uniform backing and uniform area tension.” Formal review and commit pending.
