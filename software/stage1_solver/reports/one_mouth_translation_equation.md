# What makes one draining opening move at 0PN?

Let $X$ be the volume centroid of the three-dimensional opening in the brane. A sideways push $f_x$ does work $f_x\,\delta X$ and leaves other deformations free. The radius-one benchmark prescribes inward number rate $J=0.25$. Use the same draining preparation throughout. [Geometric probe](../src/stage1_solver/mouth_probe.py) · [Preparation](one_mouth_hold_release_requirements.md)

**The equation.** Follow a region enclosing the rim. With push $f_x$ and holder force $H_x$ as the only separately applied forces, lab-frame momentum conservation gives

$$
 f_x+H_x=\frac{dP_{\Omega,x}}{dt}+T_{\rm recv,x}+E_x-T_{\rm feed,x}.
$$

| Term | What determines it |
| --- | --- |
| $T_{\rm feed,x}$ | Sideways momentum entering from the brane through the opening, including carriage and stress. The saved carrier contribution at incident speed $U=+0.01$ is $+0.00249999985299$ in solver force units. |
| $T_{\rm recv,x}$ | Momentum continuing past the opening, using the same normal into the drain. Equal inward number current leaves the receiving traction undetermined. |
| $E_x$ | Net outward momentum at the remaining cuts, including edge and slab-face tractions. |
| $P_{\Omega,x}$ | All material, field and localized momentum inside the region. How this inventory changes when the rim translates is unknown. |

Each current is relative to its moving cut. Each traction enters once. Measuring the total actuator reaction for a given motion would give the whole right-hand side at the mouth. [Moving balance and receiving stresses](../../../docs/_measurements/mouth_moving_momentum.txt) · [Carrier integral](../../../docs/_measurements/held_crossflow_budget.txt)

**What prevents a prediction.** For a finite initial acceleration $A$ from rest, the saved finite-grid feeding load changes at order $At^2$. At the instant of release it balances none of a nonzero step push, whatever $A$ is. Two checked field couplings illustrate the ambiguity without assigning a bare mass. A brane deformation that follows the rim immediately requires additional force $c^2A$. An elastically coupled deformation requires $c^2A[1-\cos(\omega t)]$ for the same constant acceleration. Here $c$ is the normalized field shift per unit static translation and $\omega$ its relaxation-free oscillation frequency in the elastic example. Both give the same static deformation and slow inertia. Their onset forces differ. [Guarded onset check](../../../docs/_measurements/mouth_rim_push_compatibility.txt) · [Two force histories in SymPy and Mathematica](../../../docs/_measurements/mouth_rim_force_equation.txt)

**Decision.** Start with the existing sideways deformation mode. At a spherical opening a shift $\delta X$ gives normal rim displacement $\delta X\cos\theta$. The geometry equation has this mode, but the saved drain supplies neither its normalized material momentum nor its rim traction during motion. Derive the work exchanged through that boundary with the other deformations free. An independent total-force and displacement measurement could calibrate it. The first distinction to resolve is whether a finite part of the brane deformation follows the rim immediately or responds after it moves. [Geometry equation and missing inputs](one_mouth_hold_release_requirements.md#translation-mode-and-the-measured-response)

Once this coupling closes the total balance, use opposite weak pushes in resting brane to predict the two half-strength displacements. Independently obtain the total holder force $H(U)$ in settled flow. For a regular linear onset with acceleration-per-force slope $\chi$ at rest and the same flow maintained, removing the holder predicts $a_{\rm free}=-\chi H(U)+o(U)$. A delayed response requires the displacement history on the gravity time scale. The physical pushes remain unrun. Working note on base fd667ea9, pending review. [Held-state identity and step response](../../../docs/_measurements/mouth_center_step_response.txt)

## Work exchanged when the opening shifts

In the quadratic geometry action let the displacement be $\eta=q(w,t)\cos\theta$, with the geometric mouth at $w=0$, $w$ increasing into the drain, and $q(0,t)=X(t)$ measured from rest. The integral of $\cos^2\theta$ over the sphere is $B=4\pi/3$. Use the parent's coefficients with the background measure absorbed. Varying the action gives the signed work **on the opening by this deformation**,

$$
\delta W_{\rm geom\to rim}=B\,T_w(0)\,\partial_w q(0,t)\,\delta X .
$$

The opening does the opposite work on the deformation. The gradient and stiffness are evaluated at the mouth. The perturbation of the carrier can drive this gradient through the field equation. A transmitted force counted through this boundary term must enter the total balance only once. The paper's implemented port basis omits the translation dipole; the formula above projects its structural action, using the physical centroid normalization. [Geometry action](../../../research/pde/paper/pde.tex#L773) · [Implemented force channels](../../../research/pde/paper/pde.tex#L1020) · [Geometric normalization](../../../docs/_measurements/mouth_rim_dipole.json)

**How much inertia can the deformation supply?** Consider the passive geometry contribution with no distributed driving perturbation, zero initial perturbations and a frequency-independent remote condition that makes the endpoint work pairing cancel between two solutions. Let $p_s(w)$ be its response to unit mouth displacement at Laplace frequency $s$, and $D(s)=-B\,T_w(0)p_s'(0)$ the required force per displacement. Subtracting the static and driven equations and integrating once gives

$$
\frac{D(s)-D(0)}{s^2}=B\int \mu_\eta(w)\,p_0(w)\,p_s(w)\,dw .
$$

For a regular slow limit with a finite weighted norm, the coefficient of acceleration is $B\int\mu_\eta p_0^2\,dw$. It measures how much of the deformation moves per unit shift of the opening. The same coefficient can be read entirely at the mouth from the force response. If the dipole equation also has distributed source $S_s(w)\cos\theta$, the drive force for this channel becomes $D(s)X_s-B\int p_s S_s\,dw$ under the same remote condition. A responding carrier, open termination or dissipative brane can therefore change the physical force law. [SymPy and Mathematica derivations, commands and outputs](../../../docs/_measurements/rim_dipole_work.json)

**Check and decision.** For a uniform elastic segment of length $L$ with a free far end and zero local restoring term, both engines give a slow force per acceleration of $4\pi\mu_\eta L/3$ and a rapid limit of zero. This control shows why a sudden push can miss inertia relevant to slow motion. No segment parameters were assigned to the drain. The radius-one, $J=0.25$ archive contains the carrier field, current, spacings, speed and aperture mask. It supplies neither the moving rim traction nor the passive mode norm. The next discriminating datum is the **total actuator force and geometric displacement history for this same drain** under weak sideways motion at several slow drive rates. Resolve any static restoring force and velocity-dependent drag separately before inferring an acceleration coefficient. The physical four pushes and two-body free fall remain open. No new commit pending review. [Archive inventory and elastic control](../../../docs/_measurements/rim_dipole_work.json) · [Physical source absent from the existing wall solver](../src/stage1_solver/p2_driven_absorber.py#L1)

## Slow sideways force on the same draining mouth

**Available measurement.** The isolated radius-one mouth drains at prescribed rate $0.25$ in a resting lateral reservoir. On the grid with reservoir distance 24 mouth radii, six lateral cells per radius and 24 depth intervals, the archived solver translates the intake and restoring profile together. It subtracts the static load at the same position. It records the incoming carrier momentum and stress. Its output called “compensating” is their negative sum. The separate harmonic force-drive example likewise obtains its force from the carrier response. Neither supplies an independent total actuator reaction. [Actual load readout](../tools/lateral_mouth_acceleration.py#L163) · [Harmonic force construction](../tools/mouth_force_drive.py#L51)

**Prediction from two drive rates.** Use the carrier force per acceleration at $\omega=0.005$ and $0.01$ to extrapolate linearly in $\omega^2$. Reserve the separately calculated zero-frequency tangent for comparison. The boundary-law readout predicts **0.00186796136018** against **0.00186796183684**. Reconstructing flux from the interior carrier field gives **0.00186631018191** against **0.00186631065853**. These are finite-grid carrier contributions in solver mass units. The two spatial readouts differ more than either frequency extrapolation error. At $0.005$ the acceleration-dependent carrier force per unit displacement amplitude is only $-4.67039907442\times10^{-8}$. This sets a scale for resolving that particular contribution against the static load. [Source readings](../../../docs/_measurements/lateral_mouth_acceleration.json) · [SymPy and Mathematica checks](../../../docs/_measurements/mouth_slow_force_readout.json)

**What a total-force readout must supply.** Constrain just the geometric center. Record the signed force exerted by that constraint and the actual center displacement after transients settle. Obtain $K$ from opposite small static shifts of the same mouth. For $X=d\cos\omega t$, write the force's first harmonic as $d[C_\omega\cos\omega t+S_\omega\sin\omega t]$. The apparent inertia is $(K-C_\omega)/\omega^2$ and the velocity coefficient is $-S_\omega/\omega$. At a single rate the restoring and acceleration forces share a phase. A slow local inertia requires a resolved plateau as frequency decreases, with the same coefficients at half displacement amplitude. Finite-band readings can include delayed boundary forces. [Phase and work checks](../../../docs/_measurements/mouth_slow_force_readout.json)

**Stop.** The missing input is the sideways reaction of the moving opening, including momentum retained there and transferred onward. The full translation equation could return this as a constraint force, or an independent paired force and position measurement could supply it directly. Further carrier frequency sweeps cannot determine it. Physical total force, inertia and free fall remain open. Working on base fd667ea9, uncommitted pending review.
