## Referee report: one-mouth 0PN release record

### Held-equilibrium identity and its derivative

For the prepared family with the mouth held at rest under total translation support \(H_a(U)\) and incident label \(U\), replacement of that support by the work-calibrated load \(f=H_a(U)\) without impulse leaves the free initial data unchanged and must produce zero release acceleration:

\[
g_a\bigl(H_a(U),U\bigr)=0.
\]

In the symmetric zero-flow reference, \(H_a(0)=0\). Assume \(H_a\) differentiable and \(g_a\) jointly differentiable at the origin. Differentiate the identity in \(U\) and evaluate at \(U=0\):

\[
\partial_f g_a(0,0)\,H_a'(0)+\partial_U g_a(0,0)=0.
\]

Define

\[
\chi_a:=\partial_f g_a(0,0),\qquad C_a:=-H_a'(0).
\]

Then

\[
\partial_U g_a(0,0)=\chi_a C_a,
\]

and the first-order expansion of the release acceleration in the two independent drives is

\[
g_a(f,U)=\chi_a(f+C_a U)+o(|f|+|U|).
\]

The signed comparison in the candidate is exactly the odd-part form of the same relation:

\[
\lim_{u\to0}\frac{g_a(0,u)-g_a(0,-u)}{2u}
=-\Bigl(\lim_{h\to0}\frac{g_a(h,0)-g_a(-h,0)}{2h}\Bigr)
\Bigl(\lim_{s\to0}\frac{H_a(s)-H_a(-s)}{2s}\Bigr),
\]

i.e. \(\partial_U g_a(0,0)=-\chi_a H_a'(0)=\chi_a C_a\). The minus sign is required by the definition \(C_a=-H_a'(0)\) once \(H_a\) and \(f\) share the same work calibration \(H_a\delta X\) and \(f\delta X\).

### Holding-force sign, release-time acceleration, chain rule, signed readout

All four elements are justified for one and the same prepared branch under the hypotheses already fixed by the reviewed support-removal record: smooth instantaneous response, zero acceleration under the total supports, common positions and velocities, independent initial data and history held fixed when \(f\) is varied, and constraint stresses and field derivatives recomputed on release.

- **Sign of \(C_a\).** With the common work convention, \(H_a(U)\) is the support load that must be reapplied to keep \(\ddot X(0^+)=0\). The combination that appears as an effective force beside \(f\) is therefore \(C_a U=-H_a'(0)U\). The candidate correctly leaves the measured sign of \(C_a\) to the net force; positivity is not claimed a priori.
- **Release-time acceleration.** The observable is \(\ddot X(0^+)\), read as \(\lim_{t\downarrow0}2[X(t)-X(0)]/t^2\) from rest. That is the quantity controlled by the earlier record. The archived independent example with vanishing linear residual and nonzero quadratic remainder confirms that a clean linear piece can sit beside a quadratic remainder; it does not by itself supply \(\chi_a\) or \(C_a\) for the mouth.
- **Chain rule.** It is only the \(U\)-derivative of \(g_a(H_a(U),U)=0\) at the origin. It relates three response coefficients of one jointly differentiable family; it is not a substitute for a constitutive law at the opening.
- **Signed readout.** Antisymmetric limits extract the odd linear response and suppress even contamination. They remain valid only on a reproducible branch with jointly differentiable \(g_a\).

### Held state versus free state

The held state prepares the settled fields and fixes \(H_a(U)\). The free state is the same material and geometric data at the same instant, with the support replaced by the calibrated actuator \(f\) and with no impulse. The identity \(g_a(H_a(U),U)=0\) is the bridge: every held equilibrium is a free initial condition that produces zero acceleration when \(f=H_a(U)\). Differentiating that bridge yields the linear release law. Field derivatives and constraint stresses are allowed to jump with the support removal; the kinematic state \((X,\dot X)\) is not.

### Singular or history-dependent response

If the onset is singular, \(\lim_{t\downarrow0}2[X(t)-X(0)]/t^2\) fails to exist as a finite number, \(\chi_a\) is not defined as a finite response coefficient, and \(1/\chi_a\) cannot be read as an effective inertia. If preparation history fails to select a single reproducible branch, \(g_a(f,U)\) is not a well-defined function of the two labels, joint differentiability fails, and the chain-rule comparison is void. Both failures are already flagged in the candidate; they are not repaired by the formal identity.

### Meaning of the effective inertia

When \(\chi_a\) is finite and nonzero, \(1/\chi_a\) is the reciprocal of the release-time force-to-acceleration response of the translational coordinate \(X\) for this preparation family. It is a signed response coefficient at \(t=0^+\). It is not identified with the momentum content \(P_X\), not guaranteed positive without measurement, and not claimed to govern the motion at later times.

### Scope of the physical evidence

Verified from the given source material:

- The one-body balance on the moving cut is \(\dot P_X+R_x=f+L_x\), equivalently \(\dot P_X=f+L_x-R_x\), with \(L_x\) and \(R_x\) as stated and independent rim forces counted once.
- The calibrated dipolar rim load does work \(f\delta X\) with \(X\) the geometric centroid.
- The support-removal linear expansion is conditioned on smooth instantaneous response and the common-state hypotheses above; the archived example has zero linear residual and nonzero quadratic remainder.
- The periodic held neutral benchmark supplies the feeding-side current \(L_x=mJ(U)U\) and does not measure the mechanical reaction.
- Local continuations of a common feeding trace and inward current admit smooth pass-through or a localized restoring contribution; neither selects a global physical throat.
- The finite-grid isolated onset archive prescribes \(X(t)\) and reads the brane-side load; it contains no independent free-center acceleration.

Open physical inputs (correctly stated in the candidate):

- Total held-force response \(H_a(U)\) for the release protocol’s incident label \(U\).
- Instantaneous force-to-center response \(\chi_a\) of the same mouth on the same branch.
- Independent matching of the benchmark feed speed \(U_{\rm feed}\) to the release label \(U\).
- A constitutive opening law that would determine both \(H_a\) and \(\chi_a\).

Until those inputs exist, the two linear accelerations remain conditional.

### Concrete discriminant

Fix the feeding-side current at the benchmark value \(L_x=mJ(U)U\). Continue the same trace in two ways:

1. Smooth transmission of the whole carrier load through the cut. The mechanical reaction on the mouth can vanish, so \(H_a(U)=0\) near \(U=0\) and \(C_a=0\).
2. Localized restoring coupling at the rim. The same \(L_x\) coexists with \(H_a'(0)\neq0\), hence \(C_a\neq0\).

The chain-rule identity then predicts different free release accelerations \(\chi_a C_a U\) at identical \(L_x\). Measuring \(H_a(U)\) (or the flow-only release) distinguishes the continuations; the benchmark current alone does not. A second discriminant is singular onset: if the free trajectory behaves as \(X(t)-X(0)\sim t^{\alpha}\) with \(\alpha\neq2\), the finite-\(\chi_a\) reading fails while the balance \(\dot P_X=f+L_x-R_x\) can still hold.

### What must be fixed

No change to the balance, the held-equilibrium identity, the definition \(C_a=-H_a'(0)\), the chain-rule comparison, or the signed readout is required. Two scoping tightenings are required so that the text cannot be read as claiming more than the evidence supports:

1. State explicitly that \(\ddot X(0^+)=\chi_a(f+C_a U)\) is the first-order expansion of the release map under joint differentiability and smooth onset, for the prepared family at rest, and is not a closed equation of motion for later times or a derivation of \(\chi_a\) or \(C_a\) from \(L_x\).
2. State explicitly that \(1/\chi_a\) is the reciprocal release-time response coefficient of this family, and that the native prescribed-motion archive and the periodic \(L_x\) benchmark do not furnish \(\chi_a\) or \(H_a\).

No other correction is forced by the source facts.

### Verdict

The holding-force sign convention, the use of release-time acceleration, the chain-rule relation, and the signed readout are justified for the same prepared branch under the smoothness, common-state, and reproducibility hypotheses inherited from the reviewed support-removal record. The candidate’s momentum budget matches the known one-body balance. The open-input section correctly separates verified observations from the missing held-force and free-response data. No error in the record changes what can be claimed: the linear release law remains conditional on those inputs, and neither the feeding-side benchmark nor the prescribed-motion archive closes it. With the two scoping tightenings above, the one-mouth record stands as a consistent conditional identity linking held force slope and free acceleration response; it does not yet supply a measured \(\chi_a\) or \(C_a\).
