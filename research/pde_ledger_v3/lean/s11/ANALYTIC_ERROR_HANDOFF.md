# Handoff to the S11c calculation session: tail and Abel error bounds

This addresses **question 1**, the requested tail/regulator assessment and
bounded proof proposal. T1–T4 now has passing Lean proofs and local checks;
Claude and Grok independently cleared the approved fixed packet.
**The bounded T1–T4 contract is complete. This is a conditional mathematical result,
not yet a numerical error certificate for the actual scattering solve.**

The useful result is that explicit tail/moment bounds and a stable outgoing
operator realization can be converted into controlled solution errors. The
remaining task is to establish those premises for the actual reduced operator.
No S11c source, running calculation or pinned export was changed for this work.

For complex amplitudes g on the real line, a measurable multiplier b satisfying
|b(x)| <= B, and any measurable kept domain s, Lean proves the actual integral
bound

\[
 \left|\int b(x)g(x)\,dx
       -\int_s b(x)e^{-a|x|}g(x)\,dx\right|
 \le B\left[\int_{s^c}|g(x)|\,dx
             +a\int |x|\,|g(x)|\,dx\right].
\]

The premises are a >= 0, integrable g, and a finite first absolute moment.
The general theorem permits any measure on R; Lebesgue measure is the physical
specialization. It includes the cutoff s={|x|<=R}, zero regulator and zero
amplitude. The proof establishes the integrability needed for the actual
integrals. It neither projects onto the real part nor drops constant/step terms.

For a constant-plus-step profile contribution, take b=f_-+Delta f H and g to
be the appropriate folded trial/test amplitude in dimensionless position.
This application still needs the actual Fourier-pairing identity and moment
estimates for that amplitude. Bounds for a few Gaussian test fields do not
supply a uniform operator estimate over the scattering domain.

Executed native statements plus inspection of the existing `hat`/`prescribe`
Jacobian and measure mapping identify the Fourier convention and Abel
half-line transforms. With s=L_W(k-k'), q=k-k', the convolution measure
L_W/(2*pi) gives physical regulator width h=a/L_W, hence a/10 for the approved
input. The constant kernel and step kernel are

\[
 P_h(q)=\frac{h}{\pi(h^2+q^2)},\qquad
 \frac12\bigl(P_h(q)-iQ_h(q)\bigr),\quad
 Q_h(q)=\frac{q}{\pi(h^2+q^2)}.
\]

The step's odd/PV part and its factor of one half must remain. This is exact
symbolic checking of selected original source statements, not a Lean CAS bridge.
The native Gaussian Fourier mass 2*pi was checked separately; the raw native
Abel `delta_mass` Meijer-G expression was not evaluated or certified. The Abel
parameter here is a position-space damping parameter, not an outgoing-frequency
prescription.

Lean also constructs the perturbed inverse on explicit Banach spaces. Given
a supplied boundedly invertible A:X->Y, a bounded perturbation E, and

\[
 \|A^{-1}\|\le\kappa,\qquad \|E\|\le\varepsilon,
 \qquad \kappa\varepsilon<1,
\]

it proves that A+E has an actual inverse and that

\[
 \|(A+E)^{-1}\|\le\frac{\kappa}{1-\kappa\varepsilon},\qquad
 \|(A+E)^{-1}-A^{-1}\|
 \le\frac{\kappa^2\varepsilon}{1-\kappa\varepsilon}.
\]

For u=A^{-1}f and u_tilde=(A+E)^{-1}(f+delta f), it further proves

\[
 \|\widetilde u-u\|_X
 \le\frac{\kappa}{1-\kappa\varepsilon}
       \bigl(\|\delta f\|_Y+\varepsilon\|u\|_X\bigr).
\]

A fixed bounded observation/channel map multiplies this bound by its operator
norm. Errors in the observation map or flux normalization need their own bounds;
those errors are not included in the fixed-map theorem.

To apply these results, the calculation side needs to supply:

1. **The spaces and domain:** a fixed outgoing/graph-space realization A:X->Y
   for the intended spectral region, including its boundary conditions.
2. **Uniform tail and moment estimates:** identify the actual folded amplitudes
   and bound their tails and first moments uniformly on the relevant trial/test
   unit balls, with legitimate distributional pairing and integration order.
3. **One full operator error bound:** control coordinate/momentum truncation,
   quadrature, Abel and boundary/channel approximation in the same X->Y norm.
   Local differential, ordinary integrable and distributional terms may require
   different estimates before their errors can be combined.
4. **A stability margin:** establish kappa and check kappa*epsilon<1, with source
   and observation errors included. A uniform conclusion requires uniform
   estimates throughout the named spectral region. Threshold/pole/coalescence
   neighborhoods cannot silently inherit a margin proved elsewhere.

Thus sufficiently small full-operator norm error, together with bounded inverse
stability on the stated domain, is a sufficient route to convergence of the
solve. Weak/distributional regulator convergence or agreement between cutoffs
alone does not establish these hypotheses. The earlier assessment discusses
possible ways to supply them; its stronger Fourier-space estimates and profile
tail formulas were not added to this Lean increment.

Evidence: four modules plus the audit root built freshly; 41 standard-axiom
audits passed; twelve deliberately false statements each reduced to False;
sixteen positive executions (fifteen distinct statements) passed; native sign/scale/omission controls passed.
These controls test the selected mathematical claims, not every possible source
mutation. All source, dependency and input/output object hashes matched at
review authorization. Earlier reviewed proof/evidence and D4B objects remained
unchanged. No certified numerical epsilon or kappa is claimed.

The source entry points are [Tail.lean](S11AnalyticError/Tail.lean),
[Abel.lean](S11AnalyticError/Abel.lean) and
[Stability.lean](S11AnalyticError/Stability.lean). See the
[T1–T4 contract](ANALYTIC_ERROR_COVERAGE.md),
[fidelity record](ANALYTIC_ERROR_FIDELITY.md) and
[initial assessment](S11C_ANALYTIC_ASSESSMENT.md).
The approved 24-file review packet has aggregate SHA256
`79a8697aebca1ba1c1ac29e694083daed6c1df512259d07fc59444b27e2b8a48`.
This handoff is a new summary outside that frozen packet and does not change it.
Questions 2 (variable coefficients/interfaces) and 3 (pole/projector hypotheses)
are not solved by T1–T4.
