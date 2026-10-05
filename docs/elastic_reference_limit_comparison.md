# Comparing elastic reference limits during brane formation

Status: new paper comparator for substantive Claude assessment. The user authorized comparing retained-reference and newly relaxed-reference limits, not adopting either law. No scientific worker, symbolic restoration or numerical calculation is involved. All displayed new identities are proposed paper derivations.

The question is how newly ordered material contributes to throat support. Two materials can have the same current shape and local stiffness but different stored stress because they remember different unstressed states. This comparison holds physical displacement, velocity, density, order and geometry the same, and varies only the proposed elastic reference. It does not reset material positions, choose a return loop, or assign the desired electric force sign.

## Existing physics and the extra comparator assumption

The [stage-006 material specification](../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md), P1/P5/P6, supplies an order-gradient energy and the rotational shear term

\[
f_B(\chi)+\frac{\kappa}{2}|\nabla\chi|^2
+\chi f_{\rm sh},\qquad
f_{\rm sh}=\frac{\mu}{2}|\mathcal C u_d|^2,
\quad \mu=\mu_R^{(4)}>0.
\]

Here \(u_d\) is displacement, distinct from the material velocity in that source. \(\mathcal C\) denotes its declared linear curl-like measure with the same norm; no component normalization is changed or newly certified here. The [S11c-a definitions](../research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md), §§1a/1c, already provide the effective brane material map \(x(X,t)=X+u(X,t)\) and inertia \(\rho_{\rm br}^0|\partial_tu|^2/2\). Their reduction from the four-dimensional order model at a forming throat remains unjoined. The [September interpretation](native_light_em_and_vortex_throat_interpretation.md), §§3.2/5.3/5.7, preserves that distinction and the open rotational-reference and angular-momentum obligations.

To represent the two alternatives, introduce a **proposed reference datum** \(B_*\), of the same tensor type and dimensions as \(B=\mathcal C u_d\), and compare the energy density

\[
e_{\rm sh}(\chi,B,B_*)=\frac{\chi\mu}{2}|B-B_*|^2,
\qquad r=B-B_*.
\]

This is an explicit conditional extension of the original shear term, which it recovers for \(B_*=0\). The original sources do not supply \(B_*\), its physical carrier, formation rule or work budget. It is not an already-derived microrotation, new gauge redundancy or automatically free internal variable. In particular, we do **not** minimize over \(B_*\) at every instant: that would make \(B_*=B\) and erase this stiffness. A finite nonlinear, frame-indifferent material completion is not being claimed.

Work in a fixed Cartesian linearization patch, with constant \(\mu,\kappa\) and smooth prescribed fields. The norm and measure remain the same in both cases. For formal component variation only, write \(B_A=C_{Ai\alpha}\partial_i u_{d\alpha}\), with fixed coefficients representing \(\mathcal C\); \(A\) enumerates its independent components in that norm. This leaves the original component convention symbolic rather than substituting a new three-dimensional curl for a four-dimensional one.

## Two initialization limits on the same material history

Let \(t_b\) label re-ordering of a material parcel, with its prescribed physical shear measure \(B_b\) at formation. The required continuation or trace of the displacement into that state is part of the comparator's supplied data, not a derived bulk shear mode. The initialization cases are

| Limit | Reference after formation | What is retained |
|---|---|---|
| Retained reference | \(B_*^R=B_{\rm old}\), the specified pre-disordering reference | A material memory datum; its physical carrier through the disordered interval is an open input. |
| Newly relaxed reference | \(B_*^N=B_b\) at formation, then held to that formation datum during later small perturbations | The actual current configuration as the initial reference of this selected elastic measure. This is not zero physical displacement or an energy-free reset of existing ordered material. |

For this local comparison the reference components are carried with the parcel in the stated fixed linearization frame after initialization. This is a comparator transport prescription, not a finite-rotation constitutive law. Neither case adds a relaxation rate. Values must be specified consistently before testing support or a force; they cannot be fitted to the desired sign.

Set \(r_b=B_b-B_{\rm old}\). At the same order fraction \(\chi\) and formation configuration, the proposed comparison is

\[
e_R=\frac{\chi\mu}{2}|r_b|^2,\qquad e_N=0,
\qquad \partial_Be_R=\chi\mu r_b,\quad\partial_Be_N=0.
\]

The relaxed case is unstressed **in this selected shear contribution** at that configuration. It need not have zero compression, interface stress, velocity or total energy. The retained case may have a nonzero conjugate shear load; this does not yet say whether its physical mouth force points inward or outward.

For the same later perturbation \(\delta B\), with both references fixed to their respective histories,

\[
e_R(B_b+\delta B)=\frac{\chi\mu}{2}|r_b|^2
+\chi\mu r_b\!\cdot\delta B+\frac{\chi\mu}{2}|\delta B|^2,
\qquad e_N(B_b+\delta B)=\frac{\chi\mu}{2}|\delta B|^2.
\]

Thus the candidate local Hessian in \(B\) is \(\chi\mu I\) in both cases, and in displacement gradients it is \(\chi\mu C^TC\). This is only an incremental comparison at fixed \(\chi\), reference, geometry and coefficients. It does not prove equal full coupled spectra, equilibrium, confinement or stability. With the same already-supplied inertia, it explains why being born unstressed need not remove later shear-wave response. Background reference gradients, interface motion and coupled order perturbations remain relevant to a complete operator.

## First variation and the boundary terms actually available

For the **partial potential**

\[
F_{\rm part}=\int_\Omega\left[f_B(\chi)+\frac{\kappa}{2}|\nabla\chi|^2
+\frac{\chi\mu}{2}|r|^2\right]d^4X,
\qquad P_{i\alpha}=\chi\mu\sum_A C_{Ai\alpha}r_A,
\]

independent variations at fixed spatial domain and fixed remaining physical fields give the proposed identity

\[
\begin{split}
\delta F_{\rm part}={}&\int_\Omega\left[
\left(f_B'(\chi)-\kappa\Delta\chi+\frac{\mu}{2}|r|^2\right)\delta\chi
-\partial_iP_{i\alpha}\,\delta u_{d\alpha}
-\chi\mu r_A\,\delta B_{*A}\right]d^4X\\
&+\int_{\partial\Omega}\left[
\kappa\partial_N\chi\,\delta\chi
+N_iP_{i\alpha}\,\delta u_{d\alpha}\right]dA.
\end{split}
\]

The boundary term \(N_iP_{i\alpha}\) is the load conjugate to this displacement variation on a cut. It is **not yet the complete physical Cauchy traction or normal mouth force**. The rotational stress's admissibility and angular-momentum balance are still those of the original open substrate problem. Likewise, the displayed order derivative is only the contribution from \(F_{\rm part}\), not the full \(\mu_\chi\) with density/core/mixed terms.

During a mechanical perturbation after formation, \(\delta B_*=0\). The separate variation in \(B_*\) records the work associated with changing a physical reference; it is not a licence to vary that datum freely in the old equations. A variation of the actual throat shape would also change domain, order, fields and potentially reference transport. Those variations must be derived together before this cut load can be used in an electric-mouth force.

## Conversion and reference work remain explicit

For constant \(\mu\), the partial stored-energy rate along a material path is

\[
D_te_{\rm sh}=\frac{\mu}{2}|r|^2D_t\chi
+\chi\mu r\!\cdot D_tB
-\chi\mu r\!\cdot D_tB_*.
\]

The three pieces respectively record changing order, changing physical deformation and changing reference. The derivative of the integrated energy on a material volume also includes \(\int e_{\rm sh}\nabla\cdot v\,d^4X\). \(D_tB\) is the actual rate of the measure; it is not silently replaced by \(\mathcal C D_tu_d\) in an inhomogeneous advected field. No material kinetic energy, compression, gradient boundary work, native memory or external/core work is included in this partial rate.

At fixed formation deformation and reference, increasing \(\chi\) has shear-storage derivative \(\mu|r_b|^2/2\) in the retained case and zero in the newly relaxed case. This is a conditional constitutive comparison, not an imposed frozen-wall drain solution or net front power. A real moving front changes additional quantities and must obey the original material/order balances. Return remains bulk re-ordering; the total balance rate is not augmented by the separately labelled relaxation adjunct.

Initializing a reference while the parcel is fully disordered, \(\chi=0\), costs zero **in this gated term alone**; any memory/reference degrees of freedom may have other energy. Resetting a reference in already ordered material at fixed \(B,\chi\) from nonzero residual \(r\) to zero would instead change this stored energy by \(-\chi\mu|r|^2/2\). That energy cannot disappear from a total balance. Its destination is not automatically heat, radiation or wave amplification. Neither initialization limit supplies a sustaining power source just by its name.

## What this comparison can decide next

The proposed result to assess is a narrow one: the two limits can differ in initial stress and formation work while retaining the same selected local incremental stiffness. This identifies where a support prediction becomes sensitive to reference formation. It does not choose which limit the medium realizes.

Before using either case as an actual throat model, the source/physics gaps are concrete: the carrier and transport of reference data; compatibility of reference data between parcels and with the rotational symmetry/torque balance; the native four-dimensional-to-brane displacement and inertia map; the core coupling that converts these loads into a normal mouth stress; and a complete balance for conversion, material transport, memory and boundary work. A field \(B_*\) is not necessarily the curl of a globally compatible reference displacement. If compatibility fails, defects or additional state may be needed; none is inserted here.

An oriented \(\pm w\) throat must still produce like-orientation repulsion and opposite-orientation attraction. The two scalar storage comparisons neither prove nor forbid that requirement. No orientation sign is inserted into the positive stiffness to force it. The complete reflected state must transform the reference data along with displacement and geometry, and the electric source/work convention must be derived from the same core. No extra microrotation or retired polar-P construction is adopted.

This is a paper comparison with prescribed histories, not a front solver, numerical mode, core equilibrium, power integral, electric-force calculation or leakage result. All prior scientific inputs and results remain unchanged. Any next concrete calculation requires its applicable assessment/readiness; mixed numerical recovery, selected centre-drive implementation and leakage workers remain parked.
