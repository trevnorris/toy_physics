# Does shear at the mouth push on the brane?

Stored shear biases a locally flat conversion front toward the ordered side. Its equal and opposite generalized force acts on whatever holds the shear profile in place. The existing equations leave open how either force reaches the bending slab. [Order energy and equilibrium front](../../../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md#L93) · [Displacement traction and material reference](../../../docs/elastic_reference_limit_comparison.md#L61)

**One-mouth calculation.** Compare a relaxed patch and a patch with weak prescribed shear at the same mouth radius and background intake. Use a planar patch at the tip, with coordinate \(q\) pointing from ordered material into bulk along the local normal. Use the reflected equilibrium order profile as a specified preparation, \(\chi(q)=1/(1+e^{q/\delta})\), with width parameter \(\delta>0\). Let \(b=\mu|\mathrm{curl}_4u_d|^2/2\) be a localized, nonnegative shear energy density before multiplication by \(\chi\). Translating the shear profile by \(X\) and the order front by \(Y\) isolates the extra energy per unit three-area \(A\),
\[
\frac{E}{A}=\int\chi(q-Y)b(q-X)\,dq.
\]
All other fields and the local area are held fixed for this partial virtual work. A draining front may have a different profile. No trapped-mode amplitude or discharge law enters this calculation.

At coincident centers the signed forces along \(q\) are
\[
\frac{F_Y}{A}=\int\chi' b\,dq,\qquad
\frac{F_X}{A}=\int\chi b'\,dq=-\frac{F_Y}{A}.
\]
The relaxed state gives zero. Both forces are quadratic in the weak strain. For this logistic front, \(\lvert F_Y\rvert/E=\langle1-\chi\rangle_E/\delta\), where the average has weight \(\chi b\). This places the ratio between zero and \(1/\delta\). Two smooth preparations \(b=b_0\,\mathrm{sech}^{2p}(q/(2\delta))\) give \(1/(3\delta)\) for \(p=1\) and \(2/(5\delta)\) for \(p=2\). At equal stored energy the narrower profile exerts \(6/5\) as much force on the front. These strains can be prepared through the normal derivative of one tangential displacement component. SymPy integrated in the order fraction. Mathematica differentiated the energy in coordinates attached to the shear profile. [Commands and measured results](../../../docs/_measurements/conversion_shear_work.json)

**What reaches the brane.** If \(z\) is the slab bend displacement, the translation contribution for these fixed-shape profiles is \(F_z=F_X\,dX/dz+F_Y\,dY/dz\). Moving both equally gives zero work from this term. Moving the shear pattern with the slab while fixing the front gives the slab the profile reaction. Both are illustrative attachments. On this flat patch the bare curl energy's direct displacement traction for the chosen shear is tangential. The force conjugate to advancing the order front can exist even when that cut term has zero normal component.

To predict the bend, the model must specify or measure how the shear-bearing material and conversion front each move under the same slab displacement, including where the rim or outer support takes the reaction. If the strain changes shape, its additional work must also be measured. The frozen mouth already uses an external holding reaction. A charge strength across radii remains open. [Frozen support](../g0_closure_card_v0.md#L57) · [Bulk momentum balance](rounded_end_bulk_flux.md)

Working result on base e0a94594, pending independent review.
