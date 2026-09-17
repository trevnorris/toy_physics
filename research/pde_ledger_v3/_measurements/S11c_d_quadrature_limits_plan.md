# S11c-d quadrature and limit checks

The complete bounded action check is published at `009667f7`. Source and
assembled actions agree to 2.49e-16, but changing the two quadrature grids
changes the action by up to 0.04670 in the declared unit frame. Continue step 5
of the numerical action plan using the saved operands; do not rebuild the
reduced source or the test-field substitutions.

1. Inventory every source denominator that remains coordinate dependent after
   numerical binding, including the three bulk-momentum sites and all three
   Abel transfer pairs. Retain original and bound operands, real/imaginary
   parts and squared moduli. Record positive-regulator and real-momentum
   assumptions; sampled distances alone do not establish global coverage.
2. Transform the already computed Abel even kernel to the physical normal
   momentum coordinate using SymPy's integral change of variables. Derive its
   width and bounded antiderivative from that operand. Compare ordinary Gauss
   grids with the exact bounded mass at both centered and off-node peaks.
   Compare a rule split around the computed width with the same exact mass.
   Include all transfer-pair records; a generic momentum sample misses their
   diagonal concentration. This checks resolution at positive regulator, not
   a pointwise substitute for the distributional weak limit.
3. For every actual nested profile transform, derive its transfer coordinate
   from its phase. Bind only the approved profile instance, retain its original
   leg assignment, and evaluate zero and nonzero transfers. Compare split
   Gauss integration with independent adaptive integration at two spatial
   cutoffs. Record the absolute profile-tail integrals and integration error
   estimates. The zero-jet subtraction remains intact on both half-lines.
4. Use these results to select the next full-action refinement. Bound memory
   before increasing tensor sizes; split or change integration coordinates
   around computed concentration loci when needed. Verify the change of
   variables, Jacobian and source action joins. Reuse all saved direct and
   assembled test operands. Refine spatial/momentum domains and regulator only
   after resolving each fixed-regulator quadrature.
5. Establish the needed action-level tail and Abel weak limits before boundary
   matching. Elementary kernel/profile tests supply numerical operands and
   resolution evidence; they do not establish the full physical operator limit
   or an S-matrix. Retain exceptional domains and failures explicitly.

## Applying the Lean T1–T4 handoff

The calculation session has received
[`ANALYTIC_ERROR_HANDOFF.md`](../lean/s11/ANALYTIC_ERROR_HANDOFF.md).
Its local verification passes; independent statement-fidelity reviews are
pending in the other session. Receipt is not final acceptance of that packet,
a new Lean proof task, or a numerical scattering certificate. Preserve the
reviewed revision and its eventual disposition before claiming final reliance.
The reception checkpoint records the current handoff and fixed-packet hashes.

T1 supplies the actual complex integral estimate
`|full pairing - cutoff/Abel pairing| <= B (tailMass + a firstMoment)` for a
bounded measurable multiplier and an integrable folded amplitude with finite
first absolute moment. T2 constructs the perturbed inverse under compatible
fixed Banach-space operator bounds and `kappa * epsilon < 1`, with quantitative
solution and fixed-observation errors. The compact native check retains physical
width `a/L_W`, the constant contribution and both even and odd/PV step terms.
It is a selected source identity, not a formalization of all native integrals.

Advance step 5 through the following application obligations, reusing accepted
assembly and quadrature operands:

1. **Identify the realization and norm.** Give the complete retained reduced
   operator a fixed outgoing/graph-space map `A: X -> Y` on a named spectral
   domain, with boundary conditions and component units. Local derivative
   orders 0–3 need corresponding domain/embedding estimates; the 80 nonlocal
   operands retain their full source/field/profile/measure census. Do not infer
   an unweighted L2 outgoing inverse. A finite discretization has different
   spaces until explicit lifts/restrictions and complement control relate it
   to this realization.
2. **Identify actual folded amplitudes.** Construct the weak pairing for each
   term family from the accepted reduced operands. Keep local differential,
   ordinary integrable and constant/step distributional terms distinct. Prove
   the Fourier pairing/integrability and any integration-order change being
   used. Bound the actual amplitude tails and first moments uniformly over
   the relevant trial/test unit balls, not just the two Gaussian witnesses.
   Retain the fixed step origin, all three Abel transfer pairs, their measures,
   the odd/PV contribution and the half factor. If a proposed space lacks the
   required moments, the T1 application is unresolved; do not impose a moment
   bound by sampling or silently change the scattering domain.
3. **Assemble one error budget.** Bound spatial and momentum truncation,
   quadrature/discretization, Abel damping and boundary/channel approximations
   in the same `X -> Y` operator norm, after the justified identifications in
   step 1. Keep coefficient and derivative-domain bounds for local terms and
   full folded-kernel estimates for nonlocal terms. Record each bound, domain,
   unit and source identity before combining them into epsilon. Refinement
   differences and native action agreement remain numerical evidence, not
   certified upper bounds without a separate argument. This budget concerns
   the retained operator; omitted physical higher-order terms are separate.
4. **Establish stability and observables.** Construct or justify the outgoing
   inverse and an actual bound kappa; a finite-matrix condition number alone
   does not control the full-space inverse or unresolved complement. Check
   the strict product margin on the stated domain and retain source errors.
   Apply the fixed-observation theorem only to a justified bounded channel map.
   Approximate observation maps, modal normalization and current/flux ratios
   need their own errors and nonzero incident-current bounds. Threshold, pole
   and coalescence neighborhoods retain their own domains, not a generic margin.

The current wider-box triple run supplies finite quadrature evidence for this
application and continues unchanged. Its result will select any further
fixed-regulator refinement. The theorem gives a clear certificate target; it
does not establish numerical epsilon or kappa, permit omitting unresolved
convergence checks, or solve the variable-coefficient/interface and nonlinear
pole questions. The latter is governed separately by the approved
`S11c_d_NONLINEAR_POLE_CONTRACT.md` addendum.

Preserve the full symbolic input and the builder report contract suffix. Keep
durable runs in repository _scratch/s11c/, publish validated .out files using
DataLad/git-annex, and commit each step. Use one CAS job and a silent local
completion/error watcher, with no recurring checks or model polling.
