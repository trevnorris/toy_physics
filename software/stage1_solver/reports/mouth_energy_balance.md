# Energy required by the ideal 0PN mouth force

**Deliverable:** account for exterior energy flux, pressure work and stored flow energy, then identify the additional energy condition at the opening. This uses the same ideal constant-density, inviscid potential exteriors as the [momentum calculation](mouth_conversion_drag.md). It does not specify the physical conversion law or inertia.

**Energy crossing a translating cut.** Let $\mathbf v$ be laboratory material velocity, $\mathbf V$ the cut velocity, $e$ specific internal energy, and $p$ pressure. The outward energy rate is

$$
\Phi_E=\int\left[\rho(e+|\mathbf v|^2/2)(\mathbf v-\mathbf V)+p\mathbf v\right]\cdot\mathbf n\,dA.
$$

Pressure work uses material velocity. For steady translation, integrating the ideal radial source/sink plus a uniform stream and decaying dipole gives $\Phi_E=q(H+U^2/2)$, where $q$ is signed outward mass rate and $H=e+p_\infty/\rho$. The unknown internal-energy and pressure levels remain explicit. Equal mass transferred between the 3D brane and 4D bulk does not imply equal $H$.

**Power available and power delivered.** Set $q_f=-D$, $q_b=D$, $D>0$, and let $U,W,V$ be the signed lab speeds of distant brane, distant bulk and mouth. The steady single-holder branch gives

$$
\mathcal A=D\left[H_f-H_b+\frac{U^2-W^2}{2}\right],\qquad
F=D(U-W),\qquad P_{\rm mech}=FV,
$$
$$
\mathcal R=\mathcal A-FV
=D\left[H_f-H_b+\frac{(U-V)^2-(W-V)^2}{2}\right].
$$

Here $\mathcal A$ is net energy arrival from the exterior cuts. $FV$ is power delivered to the holder during constant-speed motion. Applying it to a free center also needs that center's dynamics. The remainder $\mathcal R$ is frame invariant. Its required boundary balance is

$$
\mathcal R+P_{\rm supply}=\dot E_{\rm unresolved}+P_{\rm other,out}.
$$

These terms cover storage and transfers omitted from the two ideal cuts. Order energy already included in $H_f-H_b$ cannot also be counted as supplied power. Each additional supplier or recipient must be identified.

**Diagnostic, not a physical failure.** With equal heads, $W=0$ and $V=U$, the remainder is $-DU^2/2$. That steady preparation requires extra power or declining stored energy; a passive stationary interface cannot sustain it from these inputs alone. A nonzero physical head difference could change this conclusion. Towing through both resting media instead gives zero mechanical power and $\mathcal R=D(H_f-H_b)$. Thus zero drag does not establish an energy-free conversion process.

**Stored energy was retained.** In either exterior of dimension $n=3,4$, write the lab dipole potential as $Ux+B(t)x/r^n$ and use a fixed-radius translating shell $a<r<R$. With $q,U,\rho,e$ fixed and $\Omega_n$ the area of the unit sphere, the excess kinetic energy above the radial and uniform backgrounds is

$$
E_{\rm dip}=\frac{\rho\Omega_n(n-1)B^2}{2n}(a^{-n}-R^{-n}).
$$

The unsteady Bernoulli pressure includes $-\rho\dot B\cos\theta/r^{n-1}$. It contributes $-\rho\Omega_n\dot B[U-(n-1)B/r^n]/n$ to the energy flux, giving $\Phi_E(R)-\Phi_E(a)=-\dot E_{\rm dip}$. Omitting that pressure term fails the shell balance. The physical mouth condition still has to select $B(t)$; this energy is not an assigned mouth mass.

**Conclusion and evidence.** The exterior calculation supplies a required energy budget, but no available physical conversion power. The missing datum is the mouth-frame total energy drop per transferred mass, together with its storage/export/supply relation. The [adopted energy ledger](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#6-energy-input-and-forcepower-pairing) explicitly leaves these open. Maintained incident streams and pressure levels also do not establish the energy budget of an isolated falling pair. SymPy and Mathematica outputs and guard validation are in the [run manifest](../../../docs/_measurements/mouth_energy_balance.json). Base `fd667ea9291b4430a96afe3fd932926232f5507d`; working, uncommitted pending review. Interior energetics and free release remain deferred.
