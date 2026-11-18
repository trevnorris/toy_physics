# A mouth-level equation for a weak sideways push at 0PN

**What the rim supplies.** Linearize about one steady inward drain in resting brane. A sideways translation has rim displacement \(q(w,t)\cos\theta\), with \(X=q(0,t)\) the geometric center displacement to first order. Absorb the angular factor \(4\pi/3\) in the line coefficients \(\mu,T,H\). The wall equation is \(\mu\ddot q-(Tq')'+Hq=S\), where \(S\) includes the responding carrier load. Boundary work gives the force required at the mouth for the wall channel as \(-T(0)q'(0,t)\). This is a usable force if the slope is known. The inward rate fixes a different quantity. The full control-volume force also accounts for the stress and momentum carried through the receiving and edge cuts and the change of stored momentum. Internal reactions count once. [Wall action and forcing](../../../research/pde/paper/pde.tex#L773) · [Rim projection](one_mouth_translation_equation.md#work-exchanged-when-the-opening-shifts) · [Moving-cut balance](one_mouth_moving_momentum.md)

**Eliminating the fields.** Choose a cut past the mouth and call its lateral displacement \(Y\). At real Laplace rate \(s>0\), with zero initial perturbations, let the resolved region require endpoint forces \(AX+BY\) and \(CX+DY\) for imposed displacements \(X,Y\). A continuation that requires force \(K_e(s)Y\) exerts the opposite force on this region. Eliminating \(Y\) gives

\[
\widehat F_{\rm push}(s)=\left[A-B(D+K_e)^{-1}C\right]\widehat X(s).
\]

This is the free-center law when the resolved response and termination together include all channels, \(D+K_e\) is invertible and the push is the only independent incremental drive. An interface carrying several independent traces needs the corresponding matrix and conjugate forces. Eliminating them still leaves a scalar total response along the chosen push. The feeding-side stress alone does not give these blocks. [Open surface and trapped-field inputs](../../em_charge_attribute/u1_body_dynamics_inputs.yaml#L48)

**Why the missing response matters even in slow fall.** For the homogeneous passive wall equation let \(p_s(w)\) be its profile for unit rim displacement, with \(-T(L)p_s'(L)=K_e(s)p_s(L)\) at the cut. Subtract the static equation multiplied by \(p_s\) from the driven equation multiplied by \(p_0\) and integrate. If \(K_e(s)=K_e(0)+M_e s^2+o(s^2)\) and \(p_s\) tends regularly to \(p_0\), the rim's slow force per acceleration after static subtraction is

\[
M_{\rm wall+end}=\int_0^L\mu p_0^2\,dw+p_0(L)^2 M_e .
\]

This wall-channel result assumes no independent distributed drive or linear drag term. The second term is inertia transmitted across the cut, weighted by how far that end moves. No bare rim mass was introduced.

**Nonuniqueness check.** Keep the near-mouth equations fixed. As diagnostic choices for the unclosed lateral response beyond the cut, use two passive field layers with a traction-free reflecting far end. Each has length \(b\), wave speed \(c\), and characteristic impedance \(Z_i=\mu_i c\). Both have zero static translation force. Their slow inertia is \((b/c)(Z_1+Z_2)\). Scaling the layer densities and tensions together changes that inertia at fixed speed, geometry and static force. Even keeping that inertia fixed, swapping unequal layers changes the force response at order \(s^4\). These examples modify only the perturbation sector and illustrate the freedom left by the supplied equations. Neither is a proposed interior solution of the active drain. [SymPy transfer and Mathematica boundary-work checks](../../../docs/_measurements/rim_end_load.json)

**Physical outcome.** The native equations determine the form of the boundary force but leave its magnitude under a push open. For the same draining preparation, a causal relation between boundary motion and total lateral reaction, computed independently or measured, suffices without a full interior solution. To predict only a regular slow acceleration, its inertial coefficient suffices once static forces, drag and other loads have been accounted for. The existing prescribed-motion runs measure the feeding contribution. No physical push trajectory or trapped-wave following rate follows from them. [Existing force and motion records](wave_fourth_reading.md)

Working result on base cda9a8eb, pending independent review.
