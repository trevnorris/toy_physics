# S11c analytic gaps: assessment and bounded proof proposal

2026-09-16. Requested by the user after consultation with the calculation
session. Governed by [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).
**Status: mathematical assessment and proposed scope, not a completed Lean
contract or a certified error budget for the running calculation.**

Prioritize tail/regulator bounds (question 1). Queue variable-coefficient bulk
identities (question 2) and nonlinear-pencil pole hypotheses (question 3).
The concrete projector counterexample below should reach the calculation
session before its pole work. Existing D4B reviews remain independent; this
document changes neither their frozen packet nor S11c sources or exports.

## 1. What the actual reduced operator permits us to claim

The governing [S11c-d specification](../../directives/S11c_d_SHARED_PHYSICS.md),
§§1c, 2 and 3b, supplies short-range positive-order profile jets, a two-ended
interface, and a frequency-dependent differential/nonlocal operator. Read-only
inspection of [the engine](../../scripts/S11c_d_mixing_scattering_sympy_audit.py)
covered `ReducedPencil`, `BoundedActionQuadrature`, `EdgeReduction.subtraction`
and `EdgeReduction.prescribe`. Output, input and middle momenta are distinct;
normal derivative terms and both full asymptotic channel pencils matter.

The current class assumes integrable consumed profile derivatives, not a
quantitative decay envelope. This gives qualitative Fourier/tail convergence
for those derivatives, **not a uniform rate over the profile class**, and does
not by itself prove convergence of the complete scattering operator. For an
unequal-asymptote zero jet, retain

\[
 f=f_-+\Delta f\,H+f_{\rm loc}.
\]

Transforming the remainder requires the additional stated premise
\(f_{\rm loc}\in L^1\). The constant and step cannot be discarded as tails.
Their delta/PV contributions encode the asymptotic interface.

The [domain report](../../_measurements/S11c_d_quadrature_domain_report.md)
already distinguishes finite-domain quadrature from complete-operator tails
and the Abel limit. Its denominator sign records are not uniform lower bounds
on an infinite spectral domain. Agreement of two quadratures or two cutoffs
is useful evidence, but neither is a certified bound on the omitted integral.

Separate the errors before combining them: finite-box quadrature, coordinate
tails, each momentum tail, the Abel limit, boundary/channel truncation, and
inversion conditioning. Parent-theory retained-order remainders are additional.

### Useful bounds, with explicit assumptions

For an integrable complex amplitude \(g\), real transfer \(s\), and \(R\ge0\),

\[
 \left|\int_{|x|>R}e^{-isx}g(x)\,dx\right|
 \le T_g(R):=\int_{|x|>R}|g(x)|\,dx.
\]

This is uniform in real transfer. If
\(M_p=\int(1+|x|)^p|g(x)|dx<\infty\), \(p>0\), then
\(T_g(R)\le M_p/(1+R)^p\). An envelope
\(|g(x)|\le C e^{-b|x|}\), \(b>0\), gives
\(T_g(R)\le2Ce^{-bR}/b\). Complex transfer requires a compatible exponential
moment and a bounded imaginary part; the real-transfer estimate cannot simply
be reused at complex poles.

There is an immediate concrete application. The
[approved input](../../_measurements/S11c_d_variable_profile_development_input.json)
has \(w=(1+\tanh\xi)/2\), \(m=\operatorname{sech}^2\xi/3\), \(L_W=10\).
Elementary integration gives, for \(R\ge0\),

\[
 T_{w-H}(R)=\log(1+e^{-2R})\le e^{-2R},\qquad
 T_m(R)=\tfrac23(1-\tanh R)\le\tfrac43e^{-2R},\qquad
 T_{w'}(R)=1-\tanh R\le2e^{-2R}.
\]

These bound the corresponding profile Fourier integrals. Higher consumed jets,
kernel factors, coordinate Jacobians and retained prefactors must also be
accounted for before claiming an action/operator bound. Gaussian test fields
permit explicit weighted tails, but scattering states are not Gaussian tests.
Also, \(w-H\) has a jump at the subtraction origin: spatial exponential decay
does not imply exponential momentum decay of every separated operand.

For multiple momentum variables, a union bound over the complement of the
integration box applies to an **absolutely integrable effective integrand**.
Its weighted moments must bound all remaining integrations uniformly. Preserve
the native integration order unless absolute integrability or another theorem
licenses an interchange. Handle cancellation-dependent distributional pieces
before applying absolute-value bounds.

### The Abel regulator has a specific meaning here

The engine uses \(e^{-a|\xi|}\), \(a>0\), on the constant/step transform.
This is separate from its Gaussian Fourier-normalization device and from an
outgoing-frequency boundary value. With \(s=L_W(k-k')\), physical width is
\(h=a/L_W\). Including the convolution measure, the constant kernel is

\[
 P_h(q)=\frac{h}{\pi(h^2+q^2)},
\]

and the step kernel is \((P_h-iQ_h)/2\), where
\(Q_h(q)=q/[\pi(h^2+q^2)]\). The odd PV contribution is essential.

A small first proof can avoid a large distribution library: for bounded
\(|b|\le B\), integrable \(g\) with finite first absolute moment, and
\(a,R\ge0\), set

\[
 I=\int b(x)g(x)dx,\qquad
 I_{R,a}=\int_{|x|\le R}b(x)e^{-a|x|}g(x)dx.
\]

Then the triangle inequality and \(0\le1-e^{-t}\le t\) give

\[
 |I-I_{R,a}|\le B\left[T_g(R)+a\int |x|\,|g(x)|dx\right].
\]

Take \(b=f_-+\Delta fH\) and identify \(g\) with the actual folded test/trial
amplitude in dimensionless position. This controls constant and step together.
The Fourier pairing identity and the required moment bounds remain explicit
application obligations; a few Gaussian witnesses do not discharge them for
all scattering fields. For an operator claim, the same right side must be
bounded uniformly on the chosen unit balls of trial and test spaces.

Two useful diagnostics guide the numerical application:

- Finite momentum boxes need a buffer around an Abel peak. For
  \(|k'|\le K-d\), \(d>0\), the mass of \(P_h(k-k')\) omitted outside
  \([-K,K]\) is at most \(2h/(\pi d)\). At \(k'=K\), the retained mass tends
  to one half, not one; shrinking the regulator does not repair that edge.
- Weak distributional convergence is insufficient for inverse stability.
  On the whole line, the Poisson convolution multiplier is \(e^{-h|\zeta|}\),
  so \(\|P_h*-I\|_{L^2\to L^2}=1\) for every \(h>0\). On inputs with
  \(|D|^\alpha g\in L^2\), \(0<\alpha\le1\), Plancherel and
  \(1-e^{-t}\le t^\alpha\) instead give
  \(\|P_h*g-g\|_2\le h^\alpha\||D|^\alpha g\|_2\).
  The conjugate Poisson/PV difference has the same bound because its extra
  Hilbert multiplier has modulus at most one. These are whole-line statements,
  not finite-box quadrature estimates. Their full formalization is deferred.

### What is sufficient for the scattering solve?

Choose a fixed outgoing problem \(A:X\to Y\) on stated Banach spaces, with a
bounded inverse of norm at most \(\kappa\). The spaces must incorporate the
domain of the differential operator and the outgoing boundary conditions.
Real continuum scattering does not supply a bounded inverse on unweighted
\(L^2\) automatically; a weighted outgoing estimate or a stable, consistently
truncated boundary problem is an additional obligation.

For \(\widetilde A=A+E\), \(\|E\|_{X\to Y}\le\varepsilon\),
\(\kappa\varepsilon<1\), and sources \(\widetilde f=f+\delta f\),

\[
 \|\widetilde A^{-1}-A^{-1}\|
 \le\frac{\kappa^2\varepsilon}{1-\kappa\varepsilon},\qquad
 \|\widetilde u-u\|_X
 \le\frac{\kappa}{1-\kappa\varepsilon}
       (\|\delta f\|_Y+\varepsilon\|u\|_X).
\]

This follows from the Neumann series and the resolvent identity. A bounded
channel extraction map transfers the solution bound to amplitudes; errors in
that map and in flux normalization must also be included. Uniform estimates
need a named spectral region with controlled inverse norm, modal gaps and
nonzero channel currents. Thresholds, poles and coalescences require separate
treatment. Resolving quadrature near a narrow peak does not establish these
stability margins.

One sufficient route to \(\varepsilon\): bound local coefficient errors in
the relevant Sobolev/graph norm, and bound the regular nonlocal error kernel
after separating local and distributional terms. Schur bounds
\(\sup_x\int\|E(x,y)\|dy\le A_0\),
\(\sup_y\int\|E(x,y)\|dx\le B_0\) imply the unweighted \(L^2\) operator
bound \(\sqrt{A_0B_0}\); weighted spaces require the corresponding conjugated
kernel. This is sufficient, not an assertion that the actual kernel satisfies
absolute Schur bounds. If it does not, retain the cancellations and find a
suitable operator estimate before choosing the formal application.

## 2. Proposed first Lean increment and stopping rule

This closes the gap between **explicit analytic error hypotheses** and their
consequences. It does not manufacture missing estimates for the full operator.

| Obligation | Bounded deliverable |
|---|---|
| T1: tails and Abel pairing | Prove the integrated tail inequality and the combined first-moment Abel/truncation bound above, including complex amplitudes, zero regulator and zero amplitude cases. Use an explicit tail majorant; sharper fractional rates are deferred. |
| T2: stability transfer | Reuse Mathlib's Banach-algebra/continuous-linear-map inversion machinery to prove the perturbed inverse and source-to-solution bounds for Banach spaces, then composition with a fixed bounded observation map. |
| T3: fidelity and coverage | Identify only the native constant/step half-line transforms, Fourier measure, origin and width \(a/L_W\). Partition application obligations into local differential, ordinary integrable kernel/profile, and distributional constant/step terms. Every retained contribution must have an applicable estimate or remain explicitly unresolved; no generic sample establishes coverage. |
| T4: controls and review | Reject wrong Abel scaling, omitted tail/first-moment terms, an unjustified strengthened bound, and inversion claims without the smallness condition. Include nonzero passing examples and an explicit failure at \(A=1,E=-1\). Obtain two independent fidelity reviews before closure. |

First instantiate the compact identification at the supplied LAB_HELD /
RHO4_CONSTANT profile. Before calling this a **numerical certificate for S11c**,
the calculation side must supply uniform folded-amplitude bounds (or suitable
operator estimates), coordinate/momentum tail constants, quadrature guarantees,
the outgoing domain, and a stability margin. Distinguish a proved theorem with
unverified premises from a fully instantiated error bound.

Stop after T1–T4 and the explicit application-obligation record. Do not prove a
general limiting-absorption theorem, a complete scattering solver, every native
integral, pole existence or parent-theory remainder estimates in this increment.
No new Lean implementation or production calculation was started for this
assessment. A library search found existing normed-operator and Banach-algebra
inversion machinery; no build is needed for this document-only proposal.

## 3. Queued: variable coefficients and interfaces

The completed [D3 bulk contract](D3_BULK_COVERAGE.md) uses constants. For smooth
spatial coefficients and \(G_{ij}=\partial_i u_j\), direct differentiation of
its same \(L=-[a(\operatorname{div}u)^2+b\operatorname{tr}(G^2)+c\|G\|^2]/2\)
gives, with the contract's variational sign,

\[
 \mathrm{EL}_j=(a+b)\partial_j\operatorname{div}u+c\Delta u_j
 +(\partial_j a)\operatorname{div}u
 +\sum_i(\partial_i b)\partial_j u_i
 +\sum_i(\partial_i c)\partial_i u_j.
\]

Thus \(b=-a,c=0\) is generally no longer bulk-null when \(a\) varies.
If \(\operatorname{div}J=(\operatorname{div}u)^2-\operatorname{tr}(G^2)\),
then \(a\operatorname{div}J=\operatorname{div}(aJ)-\nabla a\cdot J\).
Across interfaces, formulate the identity weakly with specified traces and
retain the flux/jump terms; multiplying distributions by unspecified traces is
not a boundary condition. The same warning applies to variable D4 \(\beta\):
\(\beta P=\operatorname{div}(\beta K)-\nabla\beta\cdot K\).
Propose a later compact product-rule/traction contract, not a general interface
classification. These displayed extensions have not yet been checked in Lean.

## 4. Queued, with an immediate correction: nonlinear-pencil projectors

The specification §3b calls
\(P=(2\pi i)^{-1}\oint L(z)^{-1}L'(z)dz\) a Riesz projector for isolated
poles, extending to multiple poles. **That formula is not generally an
idempotent for nonlinear pencils.** Direct counterexample on a circle about 0:

\[
 L(z)=z^2,\quad L^{-1}L'=2/z,\quad P=2,\quad P^2=4\ne P.
\]

It has algebraic multiplicity two, one-dimensional kernel, and a second-order
inverse pole whose residue is zero. These are different invariants.

For a finite-dimensional analytic regular pencil at an isolated semisimple
zero, let \(V,W\) be full right/left nullspace basis matrices and assume
\(D=WL'(z_*)V\) invertible. The residue is
\(R=VD^{-1}W\). This is the precise semisimple normalization in
[Schumacher, §7, Corollary 7.5](https://arxiv.org/html/2412.15985v1#S7).
For a contour enclosing only that zero, the contour object is \(RL'(z_*)\);
direct multiplication gives idempotence and rank \(\dim\ker L(z_*)\).
For a simple zero this reduces to the stated scalar normalization.

At defective nonlinear zeros, use root-chain/principal-part data or an
appropriate linearization with its own spectral projection and explicit maps
back to physical fields. For finite matrices the trace of the logarithmic
derivative contour counts determinant zeros with multiplicity, not the rank
of a universally defined projector. Even enclosing several distinct simple
zeros does not in general produce an idempotent on the original field space.

For the actual nonlocal pencil, first specify fixed domain/codomain spaces,
an analytic sheet away from branch cuts, Fredholm index zero, an invertible
point, and the hypotheses ensuring finite-dimensional singular data. An
operator trace/determinant requires its applicable additional assumptions.
[Beyn–Latushkin–Rottmann-Matthes](https://arxiv.org/abs/1210.3952)
treat holomorphic Fredholm pencils and contour-based recovery of spectral data;
their framework is a reference, not verification of these premises here.

With an analytic Fredholm homotopy and an admissible contour, the smallness
condition \(\sup_\Gamma\|L_{\rm ret}^{-1}R_{\ge2}\|<1\) protects contour
invertibility and the appropriate enclosed algebraic count. It does not alone
preserve the number of distinct poles, kernel dimensions or semisimplicity,
nor supply an individual pole-displacement rate. A linearized Riesz rank needs
that linearization to be established. The specification should distinguish
these claims before general multiple-pole outputs are interpreted physically.

All assessments above are read-only with respect to S11c. They do not assert
that a pole or scattering numerical result is wrong; they identify missing
hypotheses and one false general interpretation of a contour formula.
