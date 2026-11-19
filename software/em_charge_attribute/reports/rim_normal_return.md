# Does an inward drain supply a net load to the slab along \(w\)?

**The feeding field supplies a positive normal momentum current in the saved one-mouth benchmark. Its destination remains open.** A second calculation shows that a return force on the slab can cancel the entire far bend while leaving the same slope near the mouth. The bulk can also carry the compensating momentum elsewhere. Conservation constrains the sum over all recipients. [Two-face exterior](two_face_rim_bend.md) · [Native normal exchange and support](../../../research/pde_ledger_v3/steps/O2_steady_brane_balance.md#L106)

**What reaches the mouth.** Take \(+w\) from the feed toward the receiver. On the radius-one aperture of the fixed slab at inward number rate \(0.25\), integrate the stationary carrier tensor over its three-dimensional opening. In solver units the excess integrand is \(j_w^2/n+(n^5-1)/5-(n/4)\partial_w^2\log n\); the uniform bath pressure has been subtracted. On the finer grid (bath extent 24, six lateral cells per radius), the carried part was \(+0.043839159857\), excess pressure \(-0.00597325082679\) and the one-sided density-gradient stress \(-0.0297312380831\). The total was \(+0.00813467094708\). Recovering the normal curvature from the field equation instead gave \(+0.00718186190367\). Both methods gave a positive reading on all three grids; their difference limits the numerical precision. This measures one cut of the feeding carrier. [Commands, guarded results and input hashes](../../../docs/_measurements/rim_normal_return.json)

**Where the reaction goes.** Enclose the core and the receiving region whose load on the exterior slab is sought. With each material cut and mechanical contact counted once, the force exerted on those slab degrees of freedom is \(I_{\rm feed}-\Phi_{\rm other,out}+F_{\rm other,applied}-dP_w/dt\). Use excess fluxes in the same lab frame and close every side boundary. Include carried momentum and transmitted stress in the cuts, and any wall or bound-wave momentum in the storage. A stationary or cycle-averaged bounded state removes the storage rate. A later return from outside this volume must also be included when predicting the far bend. The current receiving benchmark imposes an inlet field different from the feeding solution. The electric mouth functional already supplies a signed bare load \(s g_{\chi h}\) conjugate to its scalar \(h\); its local restoring term changes that load by \(-k_m h\). Turning this into a mechanical force on the midplane requires the displacement map and the reaction of the source or its support. Adding the bare load to the feed current without that accounting could count the same transfer twice. [Inlet mismatch](../../stage1_solver/reports/mouth_inlet_trace.md) · [Electric mouth functional](../puncture_deflection_electric_sign_result.md#a-field-action-kernel-and-ledger)

**A test of cancellation.** In the static spherical midplane model with stiffness \(\kappa_m>0\), let \(F\) pull the slab upward at radius \(a\), and let a signed reaction \(-R\) act back on the slab at \(b>a\). Both loads are explicit diagnostic inputs. With the distant slab at rest the bend is

\[
m(r)=\frac{1}{4\pi\kappa_m}
\begin{cases}
F/r-R/b,&a\le r<b,\\
(F-R)/r,&r\ge b .
\end{cases}
\]

The slope measures \(Q_m(r)=-4\pi\kappa_m r^2m'(r)\). It reads \(F\) inside the return radius even for complete cancellation, and \(F-R\) outside. For a distributed return this quantity accumulates the signed load out to each radius. The discriminating test is the **outer** reading around an unforced draining mouth, together with the independent signed receiver and support momentum balance on the same preparation. The fixed-face archive lacks that pair. The normal response also leaves the in-plane force and inertia needed by 0PN open.

SymPy matched the radial fields at the return surface. Mathematica integrated the spherical loads against the radial Green function. Both checked the normal-stress identity and complete-return cancellation. Working result on base e0a94594, pending independent review.
