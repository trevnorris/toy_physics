## Referee report: one-mouth 0PN release record

### Derivative of the held-equilibrium identity

Hold the geometric center at rest under total signed support \(H_a(U)\) in settled sideways flow \(U\), with the mouth on the prepared inward-draining branch and with independent initial fields and history held fixed when the drives are varied. Write

\[
g_a(f,U)=\ddot X(0^+)
\]

for the release acceleration after the support is replaced, without impulse, by a known external force \(f\). Held equilibrium is the identity

\[
g_a\bigl(H_a(U),U\bigr)=0.
\]

In the symmetric zero-flow reference, \(H_a(0)=0\). If \(g_a\) is differentiable at the origin, differentiate the held identity in \(U\) and evaluate at \(U=0\):

\[
\partial_f g_a(0,0)\,H_a'(0)+\partial_U g_a(0,0)=0.
\]

Define

\[
\chi_a=\partial_f g_a(0,0),\qquad C_a=-H_a'(0).
\]

Then

\[
\partial_U g_a(0,0)=\chi_a C_a,
\]

and the first-order release law on that preparation family is

\[
\ddot X(0^+)=\chi_a(f+C_a U)+o(|f|+|U|).
\]

The odd finite-difference form of the same chain rule is exactly the displayed comparison in the candidate:

\[
\lim_{u\to0}\frac{g_a(0,u)-g_a(0,-u)}{2u}
=
-\Biggl(\lim_{h\to0}\frac{g_a(h,0)-g_a(-h,0)}{2h}\Biggr)
\Biggl(\lim_{s\to0}\frac{H_a(s)-H_a(-s)}{2s}\Biggr),
\]

i.e. left-hand side \(\partial_U g_a(0,0)\) equals \(-\chi_a H_a'(0)=\chi_a C_a\).

---

### Verdict (scoped)

**The holding-force sign convention, the use of release-time acceleration, the chain-rule relation, and the signed readout are justified for one and the same prepared inward-draining branch, under the smoothness and common-history hypotheses already stated in the reviewed support-removal record.** The linear release identity does not by itself determine the physical values of \(\chi_a\) or \(C_a\). Those remain open inputs. No error in the candidate changes what may be claimed at linear order on that scoped family. One wording fix is required in the rim-force gloss so that the balance matches the source counting rule without ambiguity.

---

### What is verified

1. **Balance form.** The written budget \(\dot P_X=f+L_x-R_x\) is the same statement as the source balance \(\dot P_X+R_x=f+L_x\). Momentum assigned to the moving cut, incoming feed \(L_x\), and exterior storage/pass-on \(R_x\) are partitioned consistently with “each transfer once,” provided the calibrated rim load \(f\) is not also inserted as a second independent contribution inside \(R_x\).

2. **Held–free relation.** Held and released states are values of one response map \(g_a\) on a common preparation: same positions, velocities, independent initial data, and history; constraint stresses and field derivatives recomputed on release. Held: \(f=H_a(U)\) and \(g_a=0\). Released: support replaced by prescribed \(f\), readout \(\ddot X(0^+)=g_a(f,U)\). The linear link between them is only the chain rule above.

3. **Sign of \(C_a\).** Defining \(C_a=-H_a'(0)\) makes the combination \(f+C_a U\) the first-order vanishing drive on the held locus. That is a convention fixed by the identity \(g_a(H_a(U),U)=0\); the physical sign of \(C_a\) is then whatever the net force measurement returns. It is not fixed by the periodic neutral benchmark.

4. **Release-time acceleration.** Using \(\ddot X(0^+)\) matches the conditioning of the earlier support-removal record (smooth instantaneous response, zero acceleration under the total supports). The archived independent example with vanishing linear residual and nonzero quadratic remainder is compatible with a clean linear statement plus higher-order remainder; it does not contradict the first-order claim.

5. **Effective inertia.** If \(\chi_a\) is finite and nonzero, \(1/\chi_a\) is the effective inertial coefficient of the *release response* of that draining state. It is the reciprocal linear map from the drive combination to initial center acceleration. It is not identified, without further constitutive input, with the literal mass content of \(P_X\).

6. **Open inputs (correctly scoped).**
   - Periodic held neutral benchmark: supplies \(L_x=mJ(U)U\); does not measure the mechanical reaction.
   - Local continuations of a common feeding trace and inward current: smooth transmission and localized restoring coupling are both admissible; neither selects a global physical throat.
   - Finite-grid isolated onset archive: prescribes \(X(t)\) and reads the brane-side load; supplies no independent free-center acceleration.
   Until total held-force response and instantaneous force-to-center response (or an equivalent opening constitutive law) exist for the same mouth, \(\chi_a\) and \(C_a\) stay conditional. The candidate states this correctly.

---

### What must be fixed

**Rim-force wording only.** The phrase that an independent rim force “enters with the opposite sign in \(R_x\)” can be read as putting the calibrated load into both \(f\) and \(R_x\), which would violate the source rule that independent rim forces are counted once. Replace that gloss by: the calibrated signed dipolar load appears once, as \(f\), doing work \(f\,\delta X\); \(R_x\) excludes that same calibrated contribution and retains pass-on and exterior storage. The displayed balance equation itself need not change.

No change is required to the held identity, the definitions of \(\chi_a\) and \(C_a\), the first-order release law, or the signed comparison.

---

### Singular or history-dependent response

The chain-rule identity and the signed comparison assume:

- differentiability of \(g_a\) at \((0,0)\);
- that the held curve \((H_a(U),U)\) and the force-step and flow-only releases used in the odd limits lie on one preparation family.

**If the response is singular** (\(\chi_a=0\), \(\chi_a\) infinite, or a non-differentiable cusp/deadband at the origin), \(\chi_a\) and \(C_a\) as derivatives need not exist, \(1/\chi_a\) ceases to be a finite effective inertia, and the linear release equation is not justified.

**If the response is history-dependent**, two preparations that share the same instantaneous mouth geometry and the same \(U\) (or the same \(f\)) can return different \(g_a\). Then \(g_a(H_a(U),U)=0\) still constrains the held family, but flow-only releases prepared differently need not satisfy \(\partial_U g_a(0,0)=\chi_a C_a\). Switching on a remote flow at the release instant is already excluded by the candidate as a different history; that exclusion is necessary and correct.

---

### Concrete discriminant / counterexample

**Discriminant (operational).** On one fixed inward-draining preparation branch:

1. Measure signed holding force \(H_a(s)\) for small settled \(s\), form \(C_a^{\mathrm{hold}}=-\lim(H_a(s)-H_a(-s))/(2s)\).
2. Measure force-step accelerations in resting brane, form \(\chi_a\).
3. Measure flow-only release accelerations \(g_a(0,\pm u)\), form the left-hand odd limit.
4. Test equality of that left-hand limit with \(\chi_a C_a^{\mathrm{hold}}\).
5. Halve the drives and repeat to confirm the linear window.
6. Repeat (3) after an alternate history that reaches the same geometric state and the same \(U\) (e.g. transient load then settle, versus settle under support alone).

Failure of (4) at weak drive, or disagreement between histories in (6), falsifies the use of a single differentiable \(g_a\) on a common family. Agreement of (4) with a nonzero quadratic remainder in displacement (as in the archived example) leaves the linear claim intact and confines error to higher order.

**Illustrative counterexample.** Suppose a rate-independent deadband: the mouth does not accelerate until \(|f+C U|>f_*\), and \(H_a(U)=-C U\) still holds inside the band. Then \(g_a(H_a(U),U)=0\), the odd holding limit still defines a finite \(C\), but \(\chi_a=\partial_f g_a(0,0)=0\) and the flow-only odd limit vanishes, while finite steps \(|h|>f_*\) give a nonzero response. The formal product identity at the origin holds in the trivial sense \(0=0\cdot C\), yet there is no useful finite effective inertia and no linear release law for practical drives. A hysteretic holding loop with \(H_a^{\uparrow}(0)\neq H_a^{\downarrow}(0)\) likewise makes the signed limit for \(H_a\) history-slice dependent and breaks equality with a single flow-only release series.

---

### Scope of physical evidence vs claim

| Item | Status |
|---|---|
| Control-surface balance and single counting of \(f\) | Verified against source (after wording fix) |
| Chain rule \(\partial_U g=\chi_a C_a\) on a smooth common preparation | Verified by differentiation |
| Signed odd-limit comparison | Verified as the same identity |
| Linear residual zero / quadratic remainder in prior archive | Compatible; does not prove values of \(\chi_a,C_a\) |
| \(L_x=mJ(U)U\) from periodic neutral benchmark | Verified as feed input only; no reaction |
| Selection of physical throat / transmission law | Open |
| Independent free-center \(\chi_a\) from prescribed-motion archive | Not available |
| Numerical values and physical sign of \(C_a\) | Open; must come from net-force measurement on the same branch |

**Claim that remains allowed:** for a mouth held on a stated inward-draining branch, with smooth instantaneous response and common preparation history, the initial release acceleration is \(\ddot X(0^+)=\chi_a(f+C_a U)\) to first order, with \(C_a=-H_a'(0)\), and the signed comparison is the experimental check of that relation.

**Claim that is not allowed:** that present archives already fix \(\chi_a\), the sign or magnitude of \(C_a\), a unique throat transmission law, or a general finite-time translation equation beyond the release instant.

---

### Required correction (minimal)

In the balance paragraph, replace the “opposite sign in \(R_x\)” gloss with an explicit once-only statement for the calibrated load \(f\), aligned with the source. Leave the held-equilibrium derivative, the definitions \(\chi_a=\partial_f g_a(0,0)\) and \(C_a=-H_a'(0)\), the release law, the signed comparison, and the open-input inventory unchanged.

---

### Separated lists

**Verified observations**

- Balance \(\dot P_X+R_x=f+L_x\) with work \(f\,\delta X\) at the geometric centroid.
- Held identity \(g_a(H_a(U),U)=0\) and its derivative \(\partial_U g_a(0,0)=\chi_a C_a\) under differentiability and common history.
- Signed odd-limit comparison as the finite-difference form of that derivative.
- Conditioning on smooth response and recomputed constraints at release, consistent with the prior support-removal record.
- Periodic benchmark value of \(L_x\); non-uniqueness of local transmission; prescribed-motion archive without free-center acceleration.

**Open physical inputs**

- Measured \(H_a(U)\) (hence sign and magnitude of \(C_a\)) for the physical mouth.
- Measured \(\chi_a\) from independent force-to-center response on the same preparation.
- Constitutive selection among local continuations of feed and inward current.
- Behavior when the response is singular or history-dependent (outside the stated hypotheses).
- Any extension from release-time linear response to a closed evolution equation for \(X(t)\) at later times.

**End of report.**
