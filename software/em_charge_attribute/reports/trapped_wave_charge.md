# Does the wave that supports a mouth carry its electric sign?

**The present trapped-wave equations leave the electric sign undetermined.** Their real standing mode stores energy and pushes on the wall. Its wave action \(E/\omega\) is positive for either mouth orientation. In the fixed cavity the wave vanishes at the mouth and has zero normal derivative at the reflecting end. It carries zero cycle-averaged axial energy flux. Reflecting the cavity preserves these properties while the separate material flow remains inward. The direct source obtained by varying the wave action with respect to a weak static electric potential is zero. The static mouth model also explicitly holds the trapped-wave energy constant in its electric functional. [Wave Hamiltonian](../../../notes/inner_throat/inner_throat_hard_mode.md#L459) · [Standing-mode boundary conditions](../../../notes/lepton_mass_notes.md#L266) · [Static electric functional](../g0_closure_card_v0.md#L197)

**A more general standing wave has a signed quantity.** Give two tangent shear polarizations the same spatial profile, frequency and reflecting conditions. This is a conditional extension of the scalar cavity. Let \(x,y\) be their real amplitudes in a fixed brane basis and \(M\) their common integrated inertia. The symmetry that rotates those two amplitudes conserves \(S=M(x\dot y-y\dot x)\) in the absence of torque. In terms of the two circularly polarized actions \(I_\pm\),

\[
 E=\omega(I_++I_-),\qquad S=I_+-I_- .
\]

The difference \(S\) measures the sense and amount of rotation of the material displacement. At one fixed energy it ranges continuously from \(-E/\omega\) through zero to \(+E/\omega\). A geometric reflection in \(w\) leaves the two in-brane components of a reflected preparation unchanged in the common basis. It does not select the opposite sense of polarization. In a draining cavity the wave can also exchange angular momentum with the rim and bulk. A matched intake alone supplies no torque balance. [Shear action](../../stage1_solver/reports/pathA_35_G0_freeze.md#L201) · [Reflection conventions](../../../docs/native_light_em_and_vortex_throat_interpretation.md#L4093) · [Elastic-wave polarization and angular momentum](https://pmc.ncbi.nlm.nih.gov/articles/PMC6176628/)

**Radius and energy test.** Set \(L=1\) and \(\omega(R)=\sqrt{1+4/R^2}\) in common model units, and compare prepared cavities. The radial load is the wave's outward generalized force \(-\partial_R E\) at fixed \(I_\pm\). The radii are prescribed and the calculation assumes exact polarization degeneracy. [Cavity dispersion](../../stage1_solver/reports/pathA_26_derrick.md#L21) · [SymPy and Mathematica outputs](../../../docs/_measurements/trapped_wave_charge.json)

| Radius | Circular actions \((I_+,I_-)\) | Stored energy | Radial load | Signed polarization \(S\) |
|---:|---:|---:|---:|---:|
| 1 | (1, 0) | \(\sqrt5\) | \(4/\sqrt5\) | \(+1\) |
| 2 | (2, 1) | \(3\sqrt2\) | \(3\sqrt2/4\) | \(+1\) |
| 2 | (1, 2) | \(3\sqrt2\) | \(3\sqrt2/4\) | \(-1\) |
| 2 | (3, 0) | \(3\sqrt2\) | \(3\sqrt2/4\) | \(+3\) |

Adding equal amounts to both actions can therefore raise the stored energy while leaving the imbalance fixed. The three radius-2 preparations also show exactly why the radial load cannot select one signed unit. The direct electric reading from the wave was zero in every case. As a separate detection control, coupling the potential to rotations of the polarization yields \(Q_{\rm wave}=gS\), with \(g\) supplied. A physical identification would require the rim to lock \(S\) to opposite fixed values in mirrored mouths and an electric coupling to that quantity. The previously calculated signed geometric count offers another possibility. The wave could protect that count by keeping the mouth open. The next test is whether the count can unwind at finite energy while the mouth stays open and continues to drain. [Geometric count and conditions](stationary_mouth_charge.md) · [Flow subtraction](mirrored_mouth_probe.md)

Working calculation on base e0a94594, uncommitted pending independent review.
