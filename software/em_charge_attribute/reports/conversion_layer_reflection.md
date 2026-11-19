# Does resonance spacing track the conversion layer’s midpoint?

**Broadening a stationary layer produced apparent motion under both tested inertia laws.** Translating an unchanged layer by \(0.05\) gave a phase displacement of \(0.05\) at each fixed probe frequency. The spacing of its first two resonances recovered that translation approximately. All lengths are in units of the initial rim-to-midpoint distance. [Earlier sharp-end frequency test](bend_front_attachment.md)

**Preparation.** Use one straight axial patch, rim at \(q=0\) and initial midpoint \(Y=1\). The ordered fraction \(\chi\) falls linearly from one to zero across width \(h\), centered at \(Y\). Keep the mouth radius, intake and ordered-side wave speed fixed. Through-flow corrections to propagation are neglected. The specified shear energy gives stiffness \(\mu_0\chi\). Compare transverse inertia \(\rho_0\chi^m\) within the layer for \(m=0\) (constant) and \(m=1\) (proportional to order). The ordered-side speed is \(c=\sqrt{\mu_0/\rho_0}=1\). Adopt the gapless branch for the numerical frequencies. The ramp, inertia laws and lack of shear traction from the bulk are trial assumptions. [Order-dependent stiffness](../../../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md#L98) · [Open kinetic normalization](../../../docs/moving_interface_conversion_work_proposal.md#L93)

| Layer change from \(Y=1,\ h=0.10\) | Transverse inertia | Actual midpoint shift | Shift from resonance spacing | Shift from phase at low / high probe |
|---|---|---:|---:|---:|
| Translate to \(Y=1.05\) | Constant | 0.05 | 0.049648 | 0.050000 / 0.050000 |
| Broaden to \(h=0.20\) | Constant | 0 | 0.075936 | 0.052661 / 0.077602 |
| Translate to \(Y=1.05\) | Proportional to order | 0.05 | 0.049951 | 0.050000 / 0.050000 |
| Broaden to \(h=0.20\) | Proportional to order | 0 | 0.003749 | 0.000362 / 0.003397 |

**Readouts.** For each inertia law, fix the two phase probes at the *original* resonances. For a return amplitude \(R(k)\) referred to the rim, the apparent shift is \(\Delta\arg R/(2k)\). Separately, relocate the two rim-node resonances and difference the spacing estimates \(\sqrt{2}\pi/\sqrt{k_1^2-k_0^2}\) against baseline. For \(a=Y-h/2,\ k=\omega/c,\ v=2kh/(m+1)\), the regular finite-energy solution of \((\chi U_{\chi})_{\chi}+(kh)^2\chi^m U=0\) gives
\[
R=e^{2ika}\frac{J_0(v)+iJ_1(v)}{J_0(v)-iJ_1(v)}.
\]
At long wavelengths the reflection phase divided by \(2k\) approaches \(a+h/(m+1)\). Constant inertia therefore puts the apparent reflection point at the outer edge to this order. Inertia proportional to order places it at the midpoint to this order. Higher-frequency corrections account for the different phase readings when the layer broadens.

**Returned power.** The isolated lossless shear channel returns all incident power. As a separate sensitivity check, add a passive velocity drag only inside the ramp, with dimensionless rate \(0.1\). After broadening, the low / high probes return \(0.923107/0.923116\) of incident power for constant inertia and \(0.960788/0.960788\) for inertia proportional to order (rounded). Rigid translation preserves the returned fraction. The damping example measures loss from shear only. Its destination and the actual conversion response require a material law.

**Implication.** In this patch, different displacements inferred at two fixed frequencies rule out a rigid translation as the sole change. Agreement alone cannot establish translation for arbitrary profiles. To infer motion at a draining mouth requires the transverse inertia through conversion and the frequency-dependent stress or energy transfer to other channels. The resonances alone cannot identify where the opposite rim force acts.

SymPy/SciPy used the regular Bessel solution. Mathematica integrated the wave equation from the free edge and checked a second edge cutoff. The original Wolfram output omitted frequency columns; a spectra-only run supplied them. All 96 paired scalar readouts agree to ten decimal places. Working result on base e0a94594, uncommitted pending independent review. [Commands and output](../../../docs/_measurements/conversion_layer_reflection.json) · [Analytic solution](../conversion_layer_reflection.py) · [Direct integration](../conversion_layer_reflection.wl)
