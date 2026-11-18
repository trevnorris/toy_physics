# Does the existing stress determine the mouth's energy transfer?

**Deliverable:** derive the energy current available from the existing equations and determine whether it closes the two-drain power requirement without a throat-interior solve.

**Result.** The specified neutral carrier action gives an exact energy-current formula. It does not select the two-sided energy transfer of a physical moving mouth. The full live brane stress, material/reference response and conversion energy partners remain open. Thus the required pair power is still a condition to be supplied, rather than a derived available source.

**Known-sector derivation.** With $n=|\psi|^2$, $c=\hbar^2/(2m)$ and the supplied neutral GNLS action, the canonical energy density and current are

$$
\mathcal E=c|\nabla\psi|^2+V_{\rm conf}n+U(n),\qquad
\mathbf J_E=-2c\,\operatorname{Re}(\partial_t\psi^*\nabla\psi),
$$
$$
\partial_t\mathcal E+\nabla\cdot\mathbf J_E=n\,\partial_tV_{\rm conf}.
$$

This follows from the bulk field equation with general $U$, without an assigned mouth mass. It concerns the known neutral carrier, not a completed live brane energy law. [Parent action and field equation](../../../research/pde/paper/pde.tex) · [Live energy obligations](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#6-energy-input-and-forcepower-pairing)

On a cut moving with velocity $\mathbf W$, the outward energy rate is $\int(\mathbf J_E-\mathcal E\mathbf W)\cdot\mathbf N\,dA$. In amplitude and phase variables the current contains both $-\hbar n\,\partial_t\theta\,\mathbf v$ and $-\hbar^2(\partial_tn)\nabla n/(4mn)$. Consequently a normal material current and a stationary density snapshot are insufficient for a changing mouth. One needs the energy current itself, or the corresponding time and spatial field traces, on both sides and the remaining cuts.

**What the held solver supplies.** For a stationary carrier $\psi=\phi e^{-i\mu t/\hbar}$, the formula reduces to $\mathbf J_E=\mu\mathbf j_N$. The existing held-feed implementation chooses $\mu(U)=\mu_0+U^2/2$ in its solver units while fixing inward number rate $J_N$. Its prescribed continuum energy increment is therefore $J_NU^2/2$. This is energy supplied by that feed preparation. It is not a measured moving-mouth conversion response or a check of the falling pair's power. [Implemented preparation](../tools/mouth_feed_weak_flow.py) · [Receiving-side mismatch](one_mouth_outer_holding_force.md)

**Stress work and its recipient.** Local traction $\mathbf t$ does power $\int\mathbf t\cdot\mathbf v_{\rm app}\,dA$, where $\mathbf v_{\rm app}$ is its application velocity. Splitting off center velocity gives

$$
P_{\rm traction}=\mathbf F\cdot\mathbf V+
\int\mathbf t\cdot(\mathbf v_{\rm app}-\mathbf V)\,dA.
$$

Equal opposite tractions with different application velocities can have zero resultant and nonzero power. The integrated mouth force therefore supplies only the center-work term. Furthermore, for a rigidly translated confinement profile, $\int n\partial_tV_{\rm conf}\,dV=-\mathbf F_{\rm profile}\cdot\mathbf V$, with $\mathbf F_{\rm profile}=\int n\nabla V_{\rm conf}\,dV$. That source is reciprocal profile work; counting it again as conversion fuel duplicates energy transfer.

**Check against the pair requirement.** The [ideal pair](two_drain_energy.md) requires a net inner-region interaction contribution $2P_{\rm mech}$. Evaluating whether it is supplied requires the receiving energy flux, local stress-work velocities, and order/reference/storage response of the same moving mouth. The equations above do not determine these missing data from the incoming stress or intake. An independently specified effective boundary energy relation could suffice; constructing the entire interior is not required by the balance. No such relation was selected here.

**Evidence and stop.** SymPy's direct field balance and Mathematica's Noether construction agree, including moving-cut transport, the stationary feed identity and reciprocal profile work. [Scripts, sources and guarded outputs](../../../docs/_measurements/mouth_stress_energy.json). Base `fd667ea9291b4430a96afe3fd932926232f5507d`; working, uncommitted pending independent review. Physical energy supply and free acceleration remain open.
