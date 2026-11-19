# Can a larger mouth compensate for more trapped wave energy?

**Axial balance.** Use the existing closed cylindrical cavity with independent radius \(R\) and reflecting length \(L\). Its trapped standing mode has energy \(E_{\rm w}=J\omega\) at wave action \(J\), with \(\omega^2=\omega_g^2+v_r^2/R^2+v_z^2/L^2\). The time-averaged outward force on the reflecting end is \(F_{\rm w}=-\partial_L E_{\rm w}|_{R,J}=\beta E_{\rm w}/L\), where \(\beta=v_z^2/(L^2\omega^2)\). The other end feels the opposite force. Each prepared state keeps its own action during the virtual displacement. Different preparations can have different actions, energies and sizes. [Cavity functional](../../stage1_solver/reports/pathA_26_derrick.md) · [Endpoint work](bend_front_attachment.md)

At free axial equilibrium the wave push balances the restoring force. With positive wall tension \(\sigma\), nonnegative volume-energy density \(p\), and the optional side curvature coefficient \(\kappa_b\), that balance requires

\[
 F_{\rm w}(R)=4\pi\sigma R^2+\frac{4\pi}{3}pR^3+\frac{16\pi}{9}\kappa_b .
\]

Cap tension contributes zero axial derivative at fixed \(R\). With the same material and pressure parameters in every preparation, doubling \(R\) increases the required integrated push by \(12\pi\sigma R^2+(28\pi/3)pR^3\). The pressure term here denotes the positive coefficient of cavity volume in the reduced energy. This is a necessary condition for a free, conservative equilibrium. It neither asserts a stable draining branch nor predicts the size selected by adding energy.

**Radius-independent term.** Side curvature by itself supplies a constant axial pull in this geometry. With tension and volume cost removed, the same trapped mode exerts a positive radial drive \(Jv_r^2/(R^3\omega)\), while side curvature contributes zero radial force at fixed \(L\). With \(J,v_r>0\) those two terms alone cannot select a finite equilibrium radius. Near-constant axial force over a limited range would require the constant term to dominate there and a separately consistent radial balance. [Wall-curvature calculation](order_rim_radius.md)

**Electric interpretation.** A wave can support a time-averaged exterior disturbance if the wall transmits a net load to a brane field. An axial balance between opposite internal forces does not specify that transmission. Even an exactly constant force into the existing passive bend leaves the earlier like-mouth attraction. An independent candidate is a conserved, orientation-signed field flux fixed by one supported mouth while its radius and wave energy vary. Its coupling, unit strength and force on a second mouth would have to follow from the brane and rim equations. [Radius comparison](unequal_mouth_electric.md) · [Two-mouth sign](two_mouth_thickness.md)

SymPy and Mathematica agreed on the radius increment, the wave force and the radial control. Working calculation on base e0a94594, uncommitted pending independent review. [Checked expressions and runs](../../../docs/_measurements/wave_load_vs_radius.json)
