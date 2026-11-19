# Does one electric flux spacing survive if the field can redistribute?

**Yes, under a specific mouth law.** In a spherical collar with a sufficiently strong preference for one circular polarization, one signed flux sector lies strictly below every other sector even when the electric field and both material amplitudes can vary freely. At the earlier parameter values, three-dimensional relaxation favors a unit in the trials, but a neutral configuration invalidates the earlier lower bound for that sector. A separate symmetry test asks whether the existing mouth can supply the required handed preference. [Guarded outputs and independent review](../../../research/pde_ledger_v3/_measurements/mouth_free_flux.json) · [Earlier fixed-flux-profile comparison](mouth_coherence_energy.md)

### Which fields were freed

Keep a mouth of fixed radius \(R\), matched inward intake, and a zero-potential comparison reservoir. The angular connection \(A\) represents electric displacement in the previously proposed dual description. On every sphere outside the excised mouth its curvature \(B=*dA\) has integrated flux \(2\pi n\), with electric charge \(q=nq_0\). Vary its distribution with full energy \(\gamma\int_{r\ge R}|B|^2d^3x\), including the exterior beyond the collar. With physical displacement \(\mathbf D=q_0 B/(2\pi)\) the coefficient is \(\gamma=q_0^2/(8\pi^2\epsilon)\). The pure monopole costs \(\pi\gamma n^2/R\). The squared radius times the **fluctuating** radial field has zero angular mean on each sphere. At a fixed multipole trace on the outer collar surface, the cost-minimizing vacuum continuation has both radial and tangential field. For example a dipole trace \(r_b^2\delta B_r=e\cos\theta\) costs \(2\pi\gamma e^2/(3r_b)\) outside radius \(r_b\).

Within \(R<r<\Lambda R\), the two relative-coherence fields still have material energy density \(K\sum_h|D_h\Psi_h|^2-K\sum_h g_h|\Psi_h|^2/r^2+u(\sum_h|\Psi_h|^2)^2/2\). Here \(D_h\) contains the varying exterior connection and the fixed native polarization connection of degree \(C_h=-2hs\) for helicity \(h\) and orientation \(s\). The native curvature remains uniform on each round section. Material edge traces are free. No extra penalty for angular fluctuations is needed for the existence result below. Radius, integral flux spacing, curvature binding and collar geometry remain inputs. Energy ordering across charge sectors assumes access to compensating charge. An isolated mouth needs a physical current to change sectors.

### What happened at the original parameters

At \(\Lambda=5/4,\ g_{\rm strong}=3/2,\ \gamma u/K^2=1/5\), a neutral, equal-gain dipole trial has total energy \(-0.1869334140 B_3\), including its positive field cost. The trial connection is continuous at the join to the vacuum solution and has finite energy. Smoothing its radial derivative there gives arbitrarily close energies. The old neutral lower bound with uniform flux was \(-B_3/9\). For \(|n|\ge2\), the monopole cost and the material quartic alone still put the energy at or above \(-B_3/9\). A relaxed unit trial undercuts it. Here and below \(B_3=2\pi K^2g_{\rm strong}^2(1-1/\Lambda)/(uR)\) is an energy unit.

| Attained energy divided by \(B_3\) | Equal helicity gains | Weaker gain \(=4g_{\rm strong}/5\) |
| --- | ---: | ---: |
| Neutral | −0.188121765417 | −0.117007995196 |
| Favored unit | −0.207628194173 | −0.207628194173 |
| Other unit | −0.207628194173 by symmetry | +0.012820801130 |

The final variational search used 2,000 real coefficients, both complex material fields, angular modes of every azimuth through the retained orders, radial freedom, and Maxwell multipoles through order eight. Two initializations reproduced the tabulated energies. Increasing from 1,200 to 2,000 coefficients changed the favored-unit energy by \(8.6545\times10^{-8} B_3\). These values give upper bounds to the sector infima, so the unrestricted neutral–unit ordering at this parameter choice remains open.

### A comparison over all admissible fields

Use the same equations with a wider collar \(\Lambda=5/3\), gains \(g_{\rm strong}=9/5,\ g_{\rm weak}=1/10\), and \(1/2\le m=\gamma u/K^2\le51/100\). These are example material constants, chosen to test existence. The same ones apply to both mouth orientations and all radii. At \(m=1/2\) the Wolfram comparison gives the following bounds, rounded for display, in units of \(B_3\).

| Sector | Bound |
| --- | ---: |
| Unit aligned with the native degree of the favored helicity | Upper −0.193519733506 |
| Neutral | Lower −0.184367425815 |
| Opposite unit | Lower +0.165123456790 |
| Any magnitude two or greater | Lower −0.180666796158 |

For \(D=d-iA_{\rm eff}\) and curvature \(F_{\rm eff}=dA_{\rm eff}\), the one-form \(\bar\Psi D\Psi\) is globally defined on each closed sphere. Integrating \(d(\bar\Psi D\Psi)=\overline{D\Psi}\wedge D\Psi-i|\Psi|^2F_{\rm eff}\) gives zero. Expanding the two squares \(|D_1\Psi\pm iD_2\Psi|^2\) in an oriented orthonormal frame then bounds the integrated angular gradient energy from below by either sign of \(\int |\Psi|^2 F_{\rm eff}\), including when the curvature changes sign across the sphere. For the neutral and opposite-unit sectors choose opposite signs for the two fields. Completing a square in the Maxwell fluctuation at \(m=1/2\) leaves a nonnegative product of the two densities and the square of their mean difference. The tangential field energy discarded in this lower bound is nonnegative. A second square in the radial covariant derivative accounts for the energetic cost of tracking the curvature at every radius. Its multiplier vanishes at both collar edges. The upper bound uses one admissible charged field. The lower bounds include arbitrary mixtures, zeros, and radial and angular flux redistributions. [First-order surface identity](https://doi.org/10.1007/s11040-025-09533-z)

In detail, let \(t=r/R\), \(c=3/(1+\Lambda+\Lambda^2)\), and choose
\[
p_a(t)=\frac{a(t-1)(t-\Lambda)(t+\Lambda+1)}{(1+\Lambda+\Lambda^2)t^2}.
\]
Use \(p_a,-p_a\) for the two neutral channels with \(a=g_{\rm strong}-1\). Then \((t^2p_a)'=a(ct^2-1)\) and the residual weak-density coefficient \(2-g_{\rm strong}-g_{\rm weak}-2t^2p_a^2\) stays positive. A common radial multiplier bounds the large-flux sectors after dropping their nonnegative angular gradient. The opposite unit already separates with the angular estimate alone. SymPy independently used a different radial bound and a different unit trial at \(m=1/2\). Increasing \(m\) preserves the fluctuation bound and changes the retained mean-field cost linearly. Checking both endpoints gives a strict margin of at least \(33154674220343/6262353981253440\) in units of \(B_3\) across the stated interval. The physical margin scales as \(1/R\). With both native degrees set to zero the neutral sector wins. These are bounds on fixed-sector infima, with no assertion about trapping times or uniqueness of profiles.

### What selects the sign

On the certified branch \(n=-s h_{\rm favored}\). Mirroring a mouth reverses its charge **if the same helicity is favored relative to its local projected normal**. Reversing that preference as well preserves the charge. For a neutral patch, proper in-plane rotations permit a Hermitian gain matrix \(G=g_0 I+i g_1 J\), with \(J\) the real quarter-turn. A tangent reflection that preserves the material state takes \(J\) to \(-J\). Such a state must have \(g_1=0\), even with a nonzero inward flow that preserves that reflection. The example has fractional circular contrast \(17/19\). This value is sufficient in the tested regime and has not been measured.

A single fixed spin vector projects with opposite signs on the two hemispheres. A local spin texture \(\mathbf L=\chi\hat{\mathbf n}\) could carry a uniform handed state, if the material couples to it. Holding that texture fixed in laboratory brane coordinates as the projected normal reverses gives the same preferred charge. Reversing both retains \(\chi\) and gives opposite charges. Two energetically equivalent states \(\chi\) and \(-\chi\) leave either sign available to an isolated mouth of the same orientation. The native degenerate wave action and a circularly prepared excitation supply no equilibrium preference on their own. [Trapped-wave control](../../em_charge_attribute/reports/trapped_wave_charge.md) · [Existing mirror matching](../../em_charge_attribute/reports/mirror_internal_matching.md)

### Could an interaction between mouths align their handedness?

As a separate far-separation test, give each mouth a local handed state \(\chi_i=\pm1\) and use the conditional unit law \(n_i=-s_i\chi_i\). Introduce a passive real field \(\phi\) with positive energy \(\kappa\int(|\nabla\phi|^2+\phi^2/\ell^2)d^3x/2\) and coupling \(-\lambda\sum_i\chi_i\phi(\mathbf X_i)\). A reflection-invariant realization needs \(\phi\) to transform with \(\chi_i\). Hold the two centers and source magnitudes fixed, subtract identical isolated self-energies and include the electric interaction. Minimizing \(\phi\) gives
\[
 U(d)=\frac{\chi_1\chi_2}{d}\big(C s_1s_2-A e^{-d/\ell}\big),
\qquad C=\frac{q_0^2}{4\pi\epsilon},\qquad A=\frac{\lambda^2}{4\pi\kappa}.
\]
Two mouths pointing the same way prefer aligned handedness, hence like charges, only if \(A e^{-d/\ell}>C\) when the sign preparations can be compared against a common reservoir. Their outward force on that branch is \([C-A(1+d/\ell)e^{-d/\ell}]/d^2\), which is then **inward**. Taking \(\ell\) infinite preserves the incompatibility between strict alignment preference and net repulsion. Conversely, a pair prepared with fixed like charges repels at sufficiently large separation for finite \(\ell\). This is a test of Coulomb plus this particular alignment interaction. It includes neither a drain-gravity term nor a process that lets an isolated particle change charge.

### Can forming a pair select the sign?

Take two mouths pointing to opposite sides. Assume the angular impulse delivered per unit solid angle is radial in brane coordinates, with magnitude \(l_i\) and sign \(s_i\chi_i\). Its integral over each sphere vanishes for either handedness. A stronger formation hypothesis asks two corresponding patches to receive equal and opposite impulses, compared in a common laboratory frame. Let \(M\) map the unit radial directions of the first mouth to those of the second by an orthogonal transformation. The mean squared residual impulse is
\[
 (l_1-l_2)^2+\frac{l_1 l_2}{3}\|I-\chi_1\chi_2 M\|_F^2.
\]
For equal nonzero strengths, local cancellation on patches facing the same direction requires equal handedness, giving opposite charges. Matching antipodal patches requires opposite handedness, giving equal charges. An isolated birth from a neutral state with only two charged mouths requires opposite charges by electric conservation itself. Even if every neutral pair matches the same directions, two independently born pairs can put either like or unlike charges on their upward mouths. A simultaneous reversal in either birth leaves its local balance intact. Both engines checked the two-birth cases. Identifying the physical torques and the correspondence between patches requires a formation law. [Earlier pair-formation calculation](../../em_charge_attribute/reports/forming_pair_circulation.md)

**Physical decision.** The field-redistribution test permits a radius-independent unit. To obtain universal same-side repulsion, independently formed mouths must share a rule for the handedness. The remaining choice is what carries that information in the material or during formation. One discriminator is the difference in small-amplitude energy for right- and left-circular probes on an uncharged patch with matched intake, as a function of position and mouth orientation, including a preparation with zero imposed spin. The current contact and standing-wave equations supply no such difference. Any proposed source of handedness also owes its effect on the far-brane light modes and its stress in the 0PN momentum balance. [Joint requirements](../../../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md#energy-selection-of-one-electric-flux-spacing)
