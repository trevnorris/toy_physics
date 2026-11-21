# Can the native 0PN substrate locate attachment anchors?

**The existing calculation leaves both anchor spacings possible.** I tested whether weak angular variations of shear order could grow into distinct sites on the conversion layer of the same radius-one mouth at inward intake 0.25. They decay on both meshes. The model supplies the shape and depth of the conversion layer but has no identified pair of attachment ends on which to impose an offset of one eighth of a material pitch. Its single-attachment holding force is therefore undetermined. [Native preparation and force observables](mouth_offset_brane_readout.md)

The test used the archived carrier and half-head order profile, with the carrier, reservoir values and conversion drive fixed. The local order at the comparison point on the finer mesh was 0.5000000000000019. The stationary state is stored as a function of brane radius and fourth-coordinate depth. To check that the angular symmetry of the saved state had not hidden an instability, I allowed order perturbations with angular dependence on the full brane sphere. Linearization of the scalar order equation gives the quadratic form for angular degree \(\ell\),

\[
 Q_\ell[v]=Q_0[v]+\nu\ell(\ell+1)\int W\,|v|^2/r^2\,dV,\qquad
 W=\exp[-\arg(\psi)/\nu],\quad \nu=0.4 .
\]

Here \(Q_0\) includes the derivative of the imposed conversion drive. We used the same weighted mass in each mode, fixed the perturbation to zero at the distant reservoirs, and enforced zero amplitude at the brane axis for nonspherical modes. The positive base quadratic form and the positive angular cost exclude a growing small angular pattern of any degree on these meshes. [Native driven order equation](../tools/mouth_conversion_supply.py) · [Calculation](../tools/mouth_anchor_order_spectrum.py)

| Slowest decay rate, per model time | Mesh 34 | Mesh 48 |
|---|---:|---:|
| Spherically symmetric, \(\ell=0\) | 0.3080592 | 0.3080172 |
| Dipole, \(\ell=1\) | 2.1051202 | 2.1050814 |
| Quadrupole, \(\ell=2\) | 4.8901336 | 4.8900969 |

These are relaxation rates of the order field. The result concerns infinitesimal perturbations of this particular stationary state. Defects created during formation, finite-amplitude transitions, and feedback from other material fields remain open.

**Which force fits?** With zero sideways step, the two existing crossed-strand trials predict brane-subtracted holding forces of 0.0342657283508 and 0.165540524726 model force units at anchor spans of 0.25 and 1 pitch, respectively. Their anchor positions, strand stiffness and unloaded length were inputs. The native carrier and order fields contain neither site identities nor bond connectivity. The older tagged-contact and internal-orientation models also assign their positions or graph. An angularly symmetric average can hide microscopic pair correlations. An independent comparison needs a substrate calculation that resolves connected constituents, with its pitch and force scale fixed on a separate material patch. Then one can locate the relaxed anchors, clamp the measured offset and subtract a matched control. Choosing the force nearest an uncalibrated trial would leave the free-fall input open. [Two-site trial readings](mouth_two_site_clamp.md) · [Contact model](mouth_tagged_contact.md) · [Material constitutive inputs](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#32-open-material-and-constitutive-inputs)

Working base `49dca23e`. [Measurements and review scopes](../../../research/pde_ledger_v3/_measurements/mouth_anchor_order_spectrum.json)
