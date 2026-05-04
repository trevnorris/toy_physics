# Parent Throat Action Promotion — Master Note

## What this bundle adds

This bundle continues the parent-action branch after the wall-action audit.
It does two things.

1. It proposes the **minimal nonlinear parent-complete throat action** in the
   current throat coordinates and derives its exact quadratic reduction.
2. It lifts that promoted action into the **weak-axisymmetric grouped-`P2`**
   language already used by the later moving-throat program.

## Files

- `step_14_parent_throat_action_candidate_notes.md`
- `step_15_parent_throat_action_weak_axisym_notes.md`
- `step_13_parent_throat_action_master_sympy.py`
- `step_14_parent_throat_action_candidate_sympy.py`
- `step_15_parent_throat_action_weak_axisym_sympy.py`

## Executive result

The project files already fixed the core status split:

- the present confinement-only parent action gives a **wall source/force**, not
  an autonomous throat PDE;
- the quadratic wall action is already a **consistent linear effective closure**.

The new derivation shows that a **minimal gauge-fixed nonlinear throat action**
can bridge those two statements cleanly.

The promoted action

\[
S_{\rm total}=S_\psi[\psi,A;\Sigma]+S_{\rm EM}[A]+S_\Sigma[R]
\]

with

\[
\mathcal L_\Sigma
=
\frac12\mu_\Sigma(R,w)R_t^2
-
\frac12 T_{w,\Sigma}(R,w)R_w^2
-
\frac12 T_{\Omega,\Sigma}(R,w)|\nabla_\Omega R|^2
-
U_\Sigma(R,w)
\]

has an exact quadratic limit equal to the audited linear wall action, with

\[
K_\eta
=
U_{\Sigma,RR}(R_0,w)
-
\frac{d}{dw}\Bigl(T_{w,\Sigma,R}(R_0,w)R_0'\Bigr)
+
\frac12 T_{w,\Sigma,RR}(R_0,w)(R_0')^2.
\]

The updated SymPy audit now carries the relevant linear and quadratic
integration-by-parts boundary terms explicitly before reading off the bulk
background equation and \(K_\eta\), and it discharges the corresponding
boundary densities on a concrete decaying Gaussian example by explicit
`sp.limit` evaluation. The same boundary operator is also tested on
\(\arctan w\), whose endpoint discharge is nonzero, so the zero-boundary
Gaussian checks are not just a broken boundary reader returning zero for
everything. So the quadratic recovery no longer hides the decay/compact-support
step.

So the moving-wall PDE can be promoted from effective closure to parent status
without breaking the existing reduction stack.

## Most useful new consequence

The promotion also gives a parent-level origin for the already-known
weak-axisymmetric grouped `P2` line:

\[
(20,21,22)\sim\Bigl(1,\frac12,-1\Bigr),
\qquad b=3a.
\]

The current scripts now derive this from the real-harmonic identity
\[
\int Y_{20}\,\bigl(Y^{\rm real}_{2m}\bigr)^2\,d\Omega
=
(-1)^m\!\int Y_{20}Y_{2m}Y_{2,-m}\,d\Omega,
\]
so the \(m=1\) sign is no longer inserted by hand.

But the same derivation also shows something important and honest:

> a pure wall-only anisotropy cannot close the live even gates nontrivially.

That statement should be read narrowly. In the present reduced wall-only gate
model it is the linear-algebra fact that the \(2\times 2\) coefficient matrix
for \((\delta K,\delta M)\) has nonzero determinant. The script now derives
that matrix as the zero-support / zero-mixed specialization of the full weak-
axisymmetric gate system, so the obstruction is genuinely tied back to the
later bundle equations rather than to a disconnected toy solve. It also writes
\(\delta M\) and \(\delta K\) as the actual overlap-generated wall slots
\[
\delta M=\int \delta\mu\,\beta^2\,dw,
\qquad
\delta K=\int \left[\delta T_w(\beta')^2+(\delta K_\eta+6\delta T_\Omega)\beta^2\right]dw,
\]
before specializing the gates, so the obstruction is explicitly tied to wall
kernel data instead of only to abstract slope symbols. It is still an
obstruction in this truncated system rather than a deeper standalone no-go
theorem.

The wall-only obstruction check is now mutation-guarded at the solve level:
the script perturbs the \(K_1\) gate before solving the two wall-only equations
and verifies that the resulting \((\delta K,\delta M)\) solution is no longer
the zero pair. This tests the gate solve itself rather than only substituting
the already-solved zero result back into the equations. It also perturbs the
\(dK\) coefficient in the wall-only \(K_1\) gate and verifies that the
two-equation determinant shifts, so the guard is not only an RHS
invertibility check.

So promoting `S_\Sigma` is **necessary** for parent completeness, but it is not
by itself the whole theorem. The support and Maxwell/mixed sectors still have
to participate in the actual branch realization.

## Why this matters

This promotion changes the program in a useful way.

- `K` and `M` stop being pure closure knobs and become branch integrals of a
  parent throat action.
- The old `(a,L)` sector remains the lowest axisymmetric truncation, not a
  discarded toy law.
- The grouped real `P2` wall lanes remain literal mode families of the same
  parent throat field.
- The later isotropic and weak-axisymmetric tests still apply exactly, but now
  their wall-side inputs have a parent-level source.

## Practical next move

The next clean derivation is no longer “invent more algebra.”
It is:

1. combine the promoted `S_\Sigma` wall block with the already-derived support
   and Maxwell/mixed reduced blocks,
2. export a frozen actual branch,
3. compute the resulting wall/support/transfer packet,
4. test the isotropic full-bundle surface and the weak-axisymmetric packet
   without post-target retuning.

So this bundle does not finish the theorem. It removes one very specific status
objection and converts it into a sharper branch-realization problem.
