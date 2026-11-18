# Can the bulk change the inverse-square force on a held mouth?

**0PN working result.** Symmetry allows the source's brane flow to change the outward flux of sideways momentum at order $1/d^2$ even when the first-order change in total intake is zero. The draining model leaves the sign and magnitude of this term open, including whether it vanishes. The target is held at the same radius and intake in the paired and isolated preparations.

**Complete exterior balance.** Take $+x$ toward the companion, distance $d$. Compare on a common fixed boundary of radius $R$ enclosing only the target. Take large $d$ at fixed $R$. Write $\Delta$ for paired minus isolated, $\Phi$ for outward sideways momentum per time (carriage and stress), and $H_x$ for the holder's force on the mouth. On the mouth and enclosed material together,

$$
\Delta H_x=\Delta\dot P_x+\Delta\Phi_{\rm br,x}+\Delta\Phi_{\rm bulk,x}-\Delta F_{\rm other,x}.
$$

The exterior cuts together cover the entire boundary. Opposite sides of an internal transfer cancel in this sum. An independent force on a feed, wall or support, or a change in momentum storage, still enters the balance. For a stationary comparison with unchanged other applied forces, write $C_{\rm br}=-\partial_U\Phi_{\rm br,x}$ and $B_{\rm bulk}=\partial_U\Phi_{\rm bulk,x}$ at zero drive. If the local brane speed $U=K_s/d^2+o(d^{-2})$ is the sole leading incident drive and the response is differentiable, then

$$
-\Delta H_x=(C_{\rm br}-B_{\rm bulk})K_s/d^2+o(d^{-2}).
$$

A positive $B_{\rm bulk}$ reduces the force toward the source. Other independent incident fields of this order would require their own response. [Shell observable](../../../research/pde_ledger_v3/steps/S15_two_drain_region.md) · [Signed near-mouth budget](../../../docs/drain_driven_initial_acceleration.md#from-the-s15-shell-reading-to-a-0pn-release-measurement)

**Angular check.** Consider a fixed round bulk cut and a baseline rotationally symmetric within the three brane directions. The two sides in $w$ may differ. At polar angle $\theta$ to the $w$ axis put $s=\sin\theta$ and $\mu=\hat{\mathbf r}\cdot\mathbf e_x$, where $\hat{\mathbf r}$ points radially within the brane. At each $\theta$, write the outward mass flux density as $j_0+Uj_1\mu$, the axial material velocity as $v_r\mu+U[A\mu^2+B(1-\mu^2)]$, and the pressure as $p_0+UP_1\mu$, through first order. These coefficients can vary with $\theta$. For material carriage and scalar pressure the outward momentum flux density is $jv_x+p s\mu$. Averaging over in-brane directions gives $\langle\mu\rangle=0$ and $\langle\mu^2\rangle=1/3$. Hence the integrated mass-rate derivative vanishes, while the momentum-rate derivative is

$$
\frac{4\pi R^3}{3}\int_{\rm bulk}\sin^2\theta\,
\big[j_1v_r+j_0(A+2B)+P_1\sin\theta\big]\,d\theta .
$$

Sum over all fixed bulk sectors and include any further momentum currents of the chosen model in the total. If the brane/bulk partition moves on the cut, vary its area as well or measure the undivided total. The expression is a projection of boundary data. It supplies no dynamics selecting those data. The native model specifies in-plane carriage at conversion and face-normal bulk loading. The load amplitude and its response about an active drain remain open. [Native premises 3 and 4](../../../research/pde_ledger_v3/steps/O2_steady_brane_balance.md) · [Scope of the rest-bulk equations, §8.5](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md)

**Where cancellation is established.** In the ideal steady, constant-density potential receiver, carriage plus Bernoulli pressure cancels a decaying velocity dipole exactly. A uniform bulk drift $W$ contributes $\dot M_t W$. The ideal companion's direct 4D discharge gives $W=O(d^{-3})$, so this benchmark has zero inverse-square bulk coefficient. The saved inlet has no matched physical bulk field. A separate attempt at direct matching to a flat potential receiver failed its phase check. Both algebra systems checked the angular projection and the ideal cancellation. [Guarded calculations](../../../docs/_measurements/mouth_leading_bulk_export.txt) · [Receiver comparison](one_mouth_outer_holding_force.md)

**Requirement for release.** Read the total signed bulk export on the same exterior under opposite weak sideways brane flows, with intake inward and reservoir preparation fixed. The difference divided by $2U$ estimates $B_{\rm bulk}$. Repeat at half strength to check the linear regime. An independent total holding-force reading can supply $C_{\rm br}-B_{\rm bulk}$ directly under the stated balance conditions. Predicting initial acceleration also requires the same mouth's independently calibrated acceleration per applied force. The physical bulk coefficient and free acceleration remain unmeasured. Working calculation at base fd667ea9, pending review.
