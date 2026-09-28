I will review only the files in this packet, starting with its layout and the three named governing documents, and I will not run scripts or follow paths outside the packet.**NEEDS REVISION.** The outgoing-prescription stop is an honest limit of this first stage. The revision is the independent LU binding, which does not assemble a scalar 5×5 inverse.

## Outgoing prescription

The supplied branch and mode data do not determine a real-axis singular integral for the coupled strong-symbol inverse. Stopping with the explicit inverse and the regular spectral density is the approved fallback.

What is already fixed, and already kept, is the regular piece only:

- Shared-physics reconstruction is \(\hat f(s)=\int e^{-is\xi}f\,d\xi\) and \(f=(1/2\pi)\int e^{+is\xi}\hat f\,ds\) (`S11c_d_SHARED_PHYSICS.md` lines 234–235).
- `EdgeReduction` builds `fourier_mass` from \(\int dq\int du\,\exp(-a u^2+iqu)\) (`source-definitions.md` lines 44–50). The strong symbol is the coefficient of \(e^{ik_n z}\) (`source-definitions.md` lines 167–173).
- The worker’s density is \(e^{ik_n(z-z_p)}P^{-1}/\) `fourier_mass` (`S11c_d_reference_kernel.py` lines 251–263). That is the regular integrand in this convention. The time factor \(e^{-i\omega t}\) is saved from `wave_phase` (`source-definitions.md` lines 84–85) and is not folded into that density.
- `constantFourierMass` is the Abel/profile integral (`source-definitions.md` lines 101–102). It is stored beside the density and is not used as the \(k_n\) measure. The profile Abel regulator stays a separate saved operand.

What is not fixed is the indentation of real \(k_n\) singularities of \(\det P_{\mathrm{strong}}\). `construct_end` chooses discrete end modes by sheet membership, bulk-decay certification, signed current, and \(\mathrm{orientation}\cdot\mathrm{Im}(k)\) (`source-definitions.md` lines 422–429). Those are end-pencil mode labels. The constructor plan’s hand substitution at the saved point makes the rest-acoustic \(q^2\) negative for every real \(k_n\) (`S11c_d_FORM_constructor_plan.md` lines 50–66), so that scalar bulk branch does not even display a radiating cut to copy. The amendment already says a discrete census on a decaying-bulk input does not cover a radiating continuous spectrum (`S11c_d_SCATTERING_FORM_AMENDMENT.md` lines 151–161). The plan allows the inverse-symbol ingredient to stand when that coupled prescription is not determined (`S11c_d_FORM_constructor_plan.md` lines 94–98). The worker’s status `REGULAR_SPECTRAL_DENSITY_BUILT_OUTGOING_PRESCRIPTION_PENDING` and `outgoingGreenOperatorCompleted=false` (`S11c_d_reference_kernel.py` lines 334–355) match that boundary. No extra contour, \(i\varepsilon\) insertion, or pole census belongs in this worker.

## Blocking finding

Independent direct-LU binding does not produce a numeric 5×5 matrix. In `direct_solve` (`S11c_d_reference_kernel.py` lines 302–305) each right-hand side `mp.eye(5)[:, j]` is correct, and `mp.lu_solve` returns a column. The entry is then read as `columns[j][i]`. In mpmath a single index selects a row, so that value is a 1×1 matrix, not the scalar `columns[j][i, 0]`. `mp.matrix` therefore does not receive five scalar rows, and `native` cannot turn the result into an 80-digit scalar matrix. The required comparison with the cofactor inverse (implementation gate, `S11c_d_reference_kernel_implementation.md` lines 41–46) cannot succeed on a valid invertible symbol.

Smallest repair: in that row builder, use `columns[j][i, 0]`. The layout `entry(i, j) =` component `i` of the solve against column `j` is the right orientation once the entry is a scalar.

## What holds

- Cofactors are the transpose: minor drops source row `j` and source column `i`, with sign \((-1)^{i+j}\) (`S11c_d_reference_kernel.py` lines 115–131). For zero-based indices that sign agrees with the classical adjugate. Left and right products match \(A\,\mathrm{adj}(A)=\mathrm{adj}(A)\,A=(\det A)\,I\). The strong symbol is the saved 5×5 `REFERENCE` matrix, not the six-component weak matrix (`source-definitions.md` lines 360–385).
- Symbol units follow `measure(row_i) - known(field_j)` (`source-definitions.md` lines 182–186). The worker’s row reconstruction and inverse units `field_i - row_j` (`S11c_d_reference_kernel.py` lines 183–203) match that formula in the order `(u1, u2, u3, theta, eW)`.
- Generic carriers, exact `xreplace`, unexpanded physical residuals, and probe operands are written before the numeric guards (`S11c_d_reference_kernel.py` lines 100–112 and 322–333). The decoder matches the existing `SavedCodec` (`source-definitions.md` lines 618–624). Containment is checked before the SymPy import and matches the guard’s 2 GiB, zero swap, 32 tasks, one CPU, nice at least 15, and one native thread (`s11c_guarded_run.py` lines 30–47 and 179–186; worker lines 54–75).

This does not clear a Green operator, FORM, A11, or A12. The worker has not run.

## Optional

- `kernelEntryUnits` adds the unit of \(k_n\) to the inverse units (`S11c_d_reference_kernel.py` lines 251–259). That is the unit of \(\int(\mathrm{density})\,dk_n\). The saved density is still the unintegrated quotient by the dimensionless Fourier mass. Both factors are stored.
- The mutation flips the first numerically nonzero entry in row-major order (lines 312–314). A tiny leading entry can fail the `1e-20` movement test while a larger entry would respond.
- `reduction['reductionState']` is not shown as a pickle key in the excerpts. The attributes read there do exist on `EdgeReduction`. A wrong wrapper name fails closed at restore.
