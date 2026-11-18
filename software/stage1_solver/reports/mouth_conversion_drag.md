# Does keeping the conversion velocity forbid attraction without towing drag?

**Deliverable:** check momentum compatibility of the ideal 0PN attraction and towing responses with local-velocity carriage at conversion, including stress reactions and the required radius scaling.

**Working result.** There is no contradiction at the level of the integrated steady momentum balance. The calculation leaves the mechanical reaction at the mouth open. It establishes neither a matched physical solution nor free acceleration. The adopted conversion rule fixes transported momentum; it separately permits face-normal bulk traction with unknown amplitude. [Adopted premises 3–4](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#2-adopted-premises-and-their-accounting-representation)

Take one axis and equal steady inward and outward material mass rate $D>0$. Let $U_{\rm br},U_{\rm bulk}$ be the distant lab-frame streams and $V$ the constant mouth speed. The previously calculated ideal exterior momentum rates, including pressure, supply $L=D(U_{\rm br}-V)$ inward and $E=D(U_{\rm bulk}-V)$ outward. These inputs require the separate constant-density, inviscid potential exteriors; the 3D and 4D enclosing spheres are not a common conversion surface. [Exterior calculation and conditions](one_mouth_outer_holding_force.md#conditional-exterior-force-calculation)

Let $C$ be the carried sideways momentum crossing the actual conversion region, evaluated in the mouth frame, and $T$ the separate mechanical transfer between the two material subsystems. The same $C$ occurs on both sides under the conversion premise. Its first variation retains changes in current, velocity and measure. Let $F_f,F_b$ denote forces assigned to the mouth, and $K_f,K_b$ signed momentum transferred to separately supported surroundings. With no steady storage, the subsystem balances give

$$
0=L-C-T-F_f-K_f,\qquad 0=C+T-E-F_b-K_b,
$$
$$
F_f+F_b=D(U_{\rm br}-U_{\rm bulk})-K_f-K_b.
$$

Thus the ideal single-holder branch ($K_f=K_b=0$) has zero force while towing through both resting media, and force $DU_{\rm br}$ when held in a brane stream with the distant bulk at rest. Conversion itself adds no second impulse. The receiving subsystem must transmit $F_b=C+T-E-K_b$ mechanically. These are **required reactions**, not a calculation of the physical traction law. Momentum taken by another supported structure cannot be assigned to the mouth.

**Normal pressure can have a sideways projection.** For an upper face $w=h(\mathbf x)$ with bulk above, $\mathbf n\,dA=(-\nabla h,1)d^3x$. Pressure on the material below gives $F_x=\int p\,\partial_xh\,d^3x$. On a fixed projected spherical shell, write $h=h_0(r)+\epsilon h_1(r)\mu$ and $p=p_0(r)+\epsilon p_1(r)\mu$, with $\mu=x/r$. Keeping the shape and pressure live gives

$$
\left.\partial_\epsilon F_x\right|_0=\frac{4\pi}{3}\left\{
[p_0r^2h_1]_{\rm inner}^{\rm outer}
+\int_{\rm inner}^{\rm outer}r^2(h'_0p_1-p'_0h_1)\,dr\right\}.
$$

The endpoint term is retained. A flat unchanged face has zero sideways pressure force. For example, $h_0=k_hr^2/2$, $p_1=k_pr$, $h_1=0$ on $0<r<a$ gives $4\pi k_hk_pa^5/15$. This demonstrates an allowed normal-force moment; those profiles are not a solved fluid state. The lower face needs its own orientation, and a moving projected boundary needs its boundary variation. A force on brane material counts as mouth force only when its transmission to the opening is established.

**Across radii.** If the physical far monopole carries the full source intake, $U_s=D_s/(4\pi\rho_\infty d^2)$. With the above ideal force and an independently determined scalar inertia $M_t$, the leading slow-response acceleration is $D_tD_s/(4\pi\rho_\infty M_td^2)$. A common positive $D(a)/M(a)=\lambda$ gives $G=\lambda^2/(4\pi\rho_\infty)$ across this family. This is a requirement on measured responses, not a prescription for inertia. Other leading drives, memory or cross-acceleration require their own terms. The ideal 4D inverse-cube contribution is a geometric correction, not automatically a post-Newtonian one.

**Evidence and open inputs.** SymPy and Mathematica check the balances and pressure moment using polar and Cartesian geometry respectively; commands, literal readings and guarded outcomes are in the [run manifest](../../../docs/_measurements/mouth_conversion_drag.json). The physical pressure/shape response, its reaction on the geometric mouth, $D(a)$ and independent inertia remain open. No interior model or numerical free-release prediction was added. Working on base `fd667ea9291b4430a96afe3fd932926232f5507d`; uncommitted pending review.
