# Does doubling the mouth radius change its electric strength?

**Yes, for the current fixed-source mouth.** Its orientation-dependent force attracts like orientations at all three radius pairings. Holding the mouth bend fixed gives the desired repulsion, but doubling both radii makes that force almost four times stronger. The model still needs a mouth law that selects both the sign and a common strength. [Mouth functional](../g0_closure_card_v0.md#32-functional-sign-and-conjugate-natural-condition) · [Boundary work conventions](../puncture_deflection_electric_sign_result.md#b-full-ensemble-definitions)

**Calculation.** For outward force \(F\) at center separation \(d\), form \(F_{\rm odd}=(F_{++}+F_{--}-F_{+-}-F_{-+})/4\). Positive means the orientation-dependent force repels like mouths. The entries below give \(B=4\pi\kappa\lim_{d\to\infty}d^2F_{\rm odd}\) in common model units. Use the same brane stiffness \(\kappa=7/8\), annulus width \(1/20\), local restoring coefficient \(k_m=1\), and field normalization for physical radii \(a=1,2a=2\). The datum is fixed within each row. Gravity and all orientation-even forces cancel in this comparison.

| What each mouth maintains | Small–small | Small–large | Large–large |
|---|---:|---:|---:|
| Unit mean bend across its annulus | 129.165936 | 254.114839 | 499.933291 |
| Unit integrated load as the replacement boundary condition | 0.844795 | 0.879780 | 0.916213 |
| Original unit source doing work on the bend | −0.844795 | −0.879780 | −0.916213 |

For unit applied generalized load, solve \(-\kappa\nabla^2h+k_m\eta h=\eta\) with the model’s volume-normalized annulus \(\eta\), regular center and decaying exterior. The mean bend is \(S=\int\eta h\,d^3x\). The transmitted load is \(z=1-k_mS\), giving \(z(a)=0.919127347040\) and \(z(2a)=0.957190292880\). Propagation between distant mouths contributes \(z_i z_j/(4\pi\kappa d)\) to the reciprocal response. Differentiating the energy including the work that maintains each boundary condition gives the following coefficients. Here \(\phi,q,j\) are the held bend, integrated imposed load, and source magnitude, respectively. The shared \(g\) is the original source strength when retained, or zero when replaced.

| Boundary condition | \(B_{ij}\) |
|---|---|
| Held bend | \((z_i\phi_i/S_i)(z_j\phi_j/S_j)\) |
| Held integrated load | \(z_i z_j(q_iq_j-g^2)\) |
| Fixed source | \(-z_i z_j(j_i+g)(j_j+g)\) |
| Fraction \(\lambda\) of core source work retained | \(z_i z_j[(1-2\lambda)q_iq_j-\lambda g(q_i+q_j)-g^2]\) |

**What the mouth would have to supply.** Equal positive coefficients require \(\phi(2a)/\phi(a)=0.508297494451\) for the held bend, or \(q(2a)/q(a)=0.960234713909\) for the replacement load. With \(g=0\) the fractional-work branch also requires \(\lambda<1/2\) for repulsion. With common \(g>0\) and \(0\le\lambda<1\), the determinant \(B_{11}B_{22}-B_{12}^2=-(1-\lambda)^2g^2z_1^2z_2^2(q_1-q_2)^2\). Equality of all three positive coefficients would require equal \(q\), then equal \(z\), which the annular solution excludes. At \(\lambda=1\) the sign is attractive. The exterior determines these requirements. The physical rule that fixes the bend or load remains open.

**Checks.** In the sharp spherical limit without a local restoring term, equal bends give the exact far-field ratio \(1:2:4\). A direct two-sphere field solution at \(d/a=20,40,80\) approaches that limit. Equal integrated loads give \(1:1:1\) there by prescription, whereas equal load per area gives \(1:4:16\). For the finite annulus, SymPy matched the radial solution analytically and Mathematica integrated the radial equation. The earlier constant curvature contribution is a wall-end force. Its transfer into the exterior bend remains a separate physical requirement. Working calculation on base e0a94594, pending independent review. [Guarded runs and outputs](../../../docs/_measurements/unequal_mouth_electric.json) · [Wall and end-support work](common_face_work.md)
