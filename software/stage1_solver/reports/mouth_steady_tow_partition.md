# Steady towing of the rigid reference mouth at 0PN

The two reservoirs exert nearly equal, opposite sideways forces on the same draining mouth. Their separately measured responses predict a uniform tow through resting reservoirs in two finite boxes, including the change of sign of its small residual. This calculation assumes sealed intact slab faces, an open downstream exit and a rigid cubic flare. It measures the coherent carrier force on the entire prescribed assembly. Physical flare equilibrium, force retained by a free mouth and infinite-domain drag remain open. [Guarded measurements and source hashes](../../../research/pde_ledger_v3/_measurements/mouth_steady_tow_partition.json)

### Separate brane and bulk drives

For each box, the radius is 1 and intake is 0.25. The first two force coefficients below give the medium's force on a held mouth per unit positive stream speed. The tow columns give the external force per unit positive mouth speed through resting reservoirs. Both quantities are steady, linear responses. The larger box preserves the original inner mesh.

| Box radius / roof depth | Brane drive | Bulk drive | Raw sum | Tow predicted with frame change | Separate lab tow |
|---|---:|---:|---:|---:|---:|
| 9 / 10 | +0.2495011362 | −0.2493825283 | +1.1860791 × 10⁻⁴ | −6.3624938 × 10⁻⁶ | −6.3625266 × 10⁻⁶ |
| 18 / 20 | +0.2498891653 | −0.2498691390 | +2.0026228 × 10⁻⁵ | +9.1780324 × 10⁻⁸ | +9.1774037 × 10⁻⁸ |

The stationary solve driving both reservoirs at velocity −1 reproduces the negative of the raw sum to within $3.2\times10^{-15}$. This checks linear superposition. The independent tow uses a moving-coordinate volume source, zero lateral drive at both reservoirs, and the conservative momentum force. Its largest-box reading was archived before the separate drives were calculated. Two interior virtual displacements there predict tow coefficients $9.0922432\times10^{-8}$ and $9.1338585\times10^{-8}$, against $9.0921637\times10^{-8}$ and $9.1338201\times10^{-8}$ in the independent tow. Both tests use the same linearized carrier equations. [Earlier force audit](mouth_finite_wall_force_audit.md)

### Why the raw sum overstates the difference

For the base field $f(r,w)$, changing from the common negative stream to the lab tow adds the phase perturbation $i x f$ per unit speed. In the sideways angular harmonic its radial coefficient is $i r f$. Substitution in the steady momentum stress gives the base mass current $\mathbf j_0=\operatorname{Im}(f^*\nabla f)$ through each port, after angular integration. With outward normals and the actuator sign, the tow prediction is the two positive-stream coefficients added to $\int_{\rm ports}\mathbf j_0\cdot\mathbf n\,dS$. This integral vanishes for a stationary mass-balanced continuum solution. The nodal finite-element boost was evaluated from the base field alone. Independently evaluating the analytic product with the one-sided base gradients gives net mass traces of −0.00012497018 and −0.00001993444 in the two boxes. These are errors in the remote gradient readout of the balanced weak solution. They explain why subtracting the two large forces without transforming the readout gives a spurious tow coefficient. [Frame diagnostic](../tools/mouth_tow_frame_check.py)

### Downstream momentum

Under brane-only drive, the signed bulk-side hemisphere fluxes per unit stream speed at radii 3, 6 and 12 are respectively $−2.034\times10^{-10}$, $−1.599\times10^{-10}$ and $−1.656\times10^{-10}$. The outlet disk and surrounding shelf have opposing one-sided traces of order $1.8\times10^{-6}$. Closing each hemisphere with both gives a residual near $+2.17\times10^{-8}$. This closure residual exceeds the hemisphere signal and supplies a diagnostic scale only. The measurements resolve neither a finite exported fraction nor an exact zero. The next physical comparison uses independently supported mouths and the same three drives. The existing prescribed-family results can already test how the intake and inertia change separately with radius. [Stream and cut calculation](../tools/mouth_steady_tow_partition.py)

Independent review: “Scoped clear for these Step 3 finite-box comparisons and unresolved physical wall allocation.” The reviewer independently integrated the four base mass traces and compared both interior windows with the archived lab tow.
