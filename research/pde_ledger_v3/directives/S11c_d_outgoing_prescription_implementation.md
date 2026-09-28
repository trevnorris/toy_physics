# Outgoing reference prescription: bounded implementation for review

The user said **Continue** after the [next-stage plan](S11c_d_outgoing_prescription_plan.md).
The candidate implementation is
[`S11c_d_outgoing_prescription.py`](../_measurements/S11c_d_outgoing_prescription.py).
No scientific operation or native packet restoration has run during preparation.
Execution requires the independent method/implementation gate. External packet
submission is a separate explicit permission; it has not yet occurred.

## Actual new construction

The worker consumes the accepted inverse, determinant, all 25 cofactors,
source/unit/context artifacts, and the accepted own REFERENCE modal and
pairing packets. It restores no scientific producer object and imports no
producer module. The old pairing checkpoint has a stale `/tmp` root; its
preserved project-local packet has the exact accepted byte count and SHA256,
also joined to the modal checkpoint's pairing hash. The original checkpoint
and accepted files are not rewritten.

The old EndResolventAudit already computed numerical local residues in earlier
constant-end diagnostics. Its report and inventory remain prior evidence,
including their scope and fingerprint-only heavy output limitation. The new
worker does not rerun that constructor, its Cauchy circles, bank inversions,
or mode census. Here the exact singular coefficients are derived directly
from the **new accepted scalar inverse**, and are independently compared with
the saved modal derivative/basis blocks. This supplies the join required by
the new reference-kernel representation.

1. Verify every exact route/hash before restoration. Join checkpoints,
   reference case, original/current input parameters, and the saved symbolic
   modal pencil to its original pairing packet. Select zero background for
   the reference, keeping the current physical coefficients unchanged.
2. From the saved physical acoustic wave row, compute `q_phys²(k_n)`. Require
   its negative to be an even real quadratic with strictly positive constant
   and quadratic coefficients. Thus this worker is explicitly restricted to
   the already approved evanescent real-normal-momentum slice. Select the
   outward-decaying real-axis value and require an exact full-matrix identity
   between the saved modal symbol on that branch and the bound accepted
   strong symbol. Check row/field/normal/bulk units. No complex frequency is
   inserted into a real-axis Piecewise expression.
3. Reuse the complete saved candidate list. Select real normal momenta on the
   physical sheet with the saved bulk-decay and current-normalization evidence.
   Require inherited finite polynomial coverage, zero saved certificate
   residuals, resolved reality/sheet classifications, denominator exclusions,
   and both unique normal lifts for every saved root disk. These are joins to
   the saved certificates, not a repeated census or a new global theorem.
   Require saved exact pole positions, whole mode blocks, and algebraic
   multiplicity equal to nullity, bounded here to one or two. The observed
   metadata has two eligible nullity-two current-normalized blocks. The code
   discovers the actual selection rather than supplying expected locations,
   residues or signed currents.
4. Bind the accepted determinant and cofactors at the unchanged physical
   point and differentiate only these saved scalar expressions in `k_n`.
   For an inherited multiplicity `m`, verify determinant derivatives of
   orders below `m` vanish and the `m`th does not; verify cofactor derivatives
   of orders below `m-1` vanish. The simple inverse-pole residue is then the
   actual matrix `m * adj^(m-1)(k_a) / det^(m)(k_a)`. This follows by matching
   the first nonzero Taylor coefficients; nullity two does not imply a
   second-order inverse pole.
5. Reuse saved full right/left bases, total normal derivative, frequency
   derivative and signed current. Verify coordinate, symbol, projected
   derivative and flux-frequency joins. Compare the exact scalar-derived
   residue to `R (L* P_k R)^(-1) L*`, using explicit one/two-dimensional
   scalar block inversion. There is no new physical nullspace/SVD construction,
   root search or full-symbol LU solve. SVD only checks the small projected
   derivative's conditioning. Both Laurent residuals and the derivative
   projector residual are retained.
6. Derive the whole block's spatial direction from the saved signed current,
   requiring one sign across that block. A one-sided reversal of this direction
   must change the current check and the singular contribution while the source
   inverse and residue remain fixed. This tests the prescription connection;
   it is not independent physical verification of the current.

The exact symbolic joins are required to be zero. Numerical comparisons to
the saved machine-precision modal data use `1e-8` in the unchanged reference
unit frame, including a relative residue comparison. This does not claim
50-digit modal data merely because exact expressions are evaluated with
`evalf(50)`. The projected derivative must have smallest singular value above
`1e-10 * max(1, largest)`. There is no rescan or relaxed fallback.

## Concrete prescription candidate and its review question

For each real block at `k_a`, let `A_a` be the actual computed residue and
`s_a` its outward current direction along increasing normal coordinate. The
candidate boundary distribution is the full accepted inverse interpreted by
symmetric principal value, plus `i*pi*sum(s_a*A_a*delta(k-k_a))`. The sign
corresponds to the denominator `k-k_a-i*0*s_a` under the source's `exp(+ikx)`
reconstruction.

The worker writes a concrete matrix of limits of integrals of the **actual
saved spectral density**, excluding `(k_a-epsilon,k_a+epsilon)` at each saved
exact real pole. It adds the explicit terms
`exp(i*k_a*(z-zp))*i*pi*s_a*A_a/fourier_mass`. It retains the continuous/regular
integrand rather than replacing it by a finite sum over modes. The exclusion
parameter carries normal-momentum units and is distinct from the retained
profile Abel regulator. At separated source/observation positions, the local
Fourier contour orientation is checked on both spatial sides: upper closure
for positive separation and clockwise lower closure for negative separation.

**Method review must decide whether this is the justified spatial outgoing
prescription on this fixed input.** The derivation uses the radiation-selected
real-axis branch, full block currents, and accepted real-root coverage; it
does not assume a complex-frequency `omega+i0` path or prove equivalence to a
retarded resolvent of the dissipative pencil. Review the transfer of saved
root coverage through the actual symbol/coordinate join, the sign rule, and
the distribution/convergence domain (including any coincident-position
contact term). Identify a concrete missing premise or smallest repair if
these ingredients do not suffice. Do not turn a global stability theorem,
frequency-pole campaign, new mode census or A12 witness into a prerequisite
without a source-grounded reason for this specific construction.

The candidate currently reports
`FIXED_INPUT_OUTGOING_PRESCRIPTION_CANDIDATE_BUILT`, with full outgoing Green
acceptance and full FORM false. Passing algebraic tests is not automatic
method/result acceptance. If a specific domain premise remains missing, the
exact residue and candidate prescription artifacts stand as partial work.
Nothing in this worker establishes other physical inputs, all four cases,
the two-asymptote response, frequency analyticity, A11 or A12.

## Persistence and launch gate

Each journaled call saves complete arguments and a start receipt before the
operation, then its actual return before later guards. Exact derivatives,
residues, coordinate/unit/branch joins, current checks, mutations and the
candidate integral remain saved even if a later stage fails. All consumed
files are checked again after construction. No accepted file is rewritten.

The worker checks exact cgroup/native bounds before importing SymPy or NumPy.
Launch is one existing `s11c_guarded_run.py` wrapper around the existing
normalization supervisor: 900 seconds, 2 GiB, zero swap, one CPU, nice 15,
32 tasks and one native thread. Failure is preserved without automatic retry.
The host completion hook must be armed first and target the current session.
The gate receipt must pin this exact worker and input manifest and record the
user's Continue scope decision plus the adjudicated method reviews.

Both fresh read-only review reports must finish before edits/adjudication,
with no peer report sharing or scientific execution. Their literal verdicts
are preserved. Optional wording does not trigger automatic re-review.
No scientific job has launched and no gate receipt is created by this proposal.

## Bounded corrections after both independent reviews

The immutable reviewed version remains in the fixed packet. Claude's literal
verdict is **NEEDS REVISION**; Grok's is **CLEAR FOR THIS BOUNDED STAGE**.
The [disposition](../_measurements/S11c_d_outgoing_prescription_review_disposition.md)
records the actual findings and local corrections without claiming a new
independent CLEAR. Both reviewers support the spatial current-oriented method
on this fixed slice; neither requires a new mode or frequency campaign.

The corrected implementation adds these execution conditions:

- Join `fixed_inverse * fixed_det - fixed_adj` exactly to zero, and restore
  and join the accepted regular spectral-density artifact to the same inverse,
  Fourier phase and saved mass.
- Represent the selected positive radical by `r`, with `r² = a + b*k²` and
  strictly positive saved `a,b`. For the determinant and each adjugate entry,
  save the rational-coordinate expression and exact back-substitution identity.
  Reduce its denominator along this relation and require a nonzero constant
  times a nonnegative power of `r`. Any unsupported denominator stops the job.
- For all 25 actual inverse entries, save the exact leading Laurent behavior
  at both tails via `k=±1/t`, `r=sqrt(a*t²+b)/t`, with `t>0`. Require a rational
  expression in the two coordinates and an exact back-substitution identity.
  Every nonzero entry must have a positive integer leading power of `t`.
  Because the positive-radical expression is analytic near `t=0`, these
  decaying Laurent tails also have integrable derivatives, giving convergence
  of the oscillatory tails at nonzero source/observation separation.
- Save the explicit pointwise domain `Ne(z,zp)`, the full tail/domain evidence,
  and reciprocal-momentum/radical units. The polynomial contact part is zero
  only after all tails pass. If an entry is constant or growing, preserve its
  actual expansion and stop before emitting a kernel; a general contact-term
  extension is not silently constructed. No value or complete distributional
  extension on the diagonal is claimed, including possible logarithmic behavior.
- For every excluded real-normal lift, require the exact opposite-sheet join
  `q_out(K) + scale*Q_native = 0`. A tolerance-based flag alone is insufficient.
  Save each candidate's inherited operands/classification before its guards.
- Save the actual zero reference-grade bindings separately from the unchanged
  development input file. Current-sign mutation and contour-sign identities
  remain wiring checks, not independent evidence for the physical current.

These are bounded operations on accepted saved expressions/certificates within
the existing 900-second/2-GiB ceiling. No new root, full-symbol inverse, current,
complex-frequency continuation or producer constructor is introduced. Actual
result acceptance still depends on inspecting saved evidence after execution.
