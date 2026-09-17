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
The bounded T1–T4 work subsequently completed at `ad365b5f`: its recorded
Claude and Grok fidelity reviews are clear and the closure record has no
required finding remaining. The original reception checkpoint preserves its
earlier review-pending state; the variable-coefficient reception checkpoint
records the accepted closure revision and unchanged T1 proof sources. This
establishes the bounded conditional result, not a numerical scattering
certificate or a reason to enlarge the Lean task.

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

## Applying the variable-coefficient/interface handoff

[`VARIABLE_COEFFICIENT_HANDOFF.md`](../lean/s11/VARIABLE_COEFFICIENT_HANDOFF.md)
supplies VC1–VC4 for the specified local D3 family, D4 odd density and a scalar
normal-slice interface identity. Local verification passes; its separate
fidelity reviews remain pending. This is no finding of a missing term in the
complete closed S11c operator and establishes no general transmission problem.
The new reception checkpoint pins this packet without changing it or launching
more Lean work. Apply its hypotheses to step 5's realization and subsequent
boundary matching as follows.

1. **Match the actual action and indices.** Lean uses derivative rows,
   `G_ij = partial_i u_j`; native S11c-b stores `grad_u[field][direction]`, so
   the map is `G_ij = grad_u[j][i]`. Preserve the transposed contraction of
   the b-gradient in the D3 identity. Match the density sign, coefficient
   factors, wave-amplitude normalization, retained grades and `1/L_W` coordinate
   scaling before applying a component formula. Native `operator_from_density`
   uses `epsilon * diff(density, grad_u[a][i])`; it must not be equated to the
   negative-action Lean momentum merely by its array shape.
2. **Keep profile dependence during variation.** Source inspection finds
   `construct_energy` forms coefficient-times-invariant densities before
   `build_operator` calls `operator_from_density`. That function takes density
   derivatives and then `total_derivative`, which differentiates the live
   W/MU background and retained profile jets. The inspected path therefore
   includes the product-rule mechanism. S11c-d's `ReducedActionAssembly`
   collects already-derived probe actions and checks their reconstruction;
   it does not manufacture an operator by putting profiles into an old
   constant-coefficient EL formula. These are source observations, not a new
   all-coefficient fidelity proof or a validation of every retained term.
3. **Retain weighted divergence corrections and currents.** A selected
   constant-coefficient representative is not automatically interchangeable
   with its profile-weighted divergence equivalents. For any such replacement,
   derive the coefficient-gradient correction and the boundary current from
   the original density; preserve both. Native uniform basis selection uses
   Euler signatures before coefficient binding, so its representative choices
   must be matched explicitly when applying the VC identities. No error in
   those choices is inferred here. A D4 odd identity requires an actual D4
   action/normalization identification; five reduced field slots do not imply
   four spatial derivative directions. The existing Q9 dependency trace found
   no consumed Q9-family input in S11c-d; this handoff does not alter that trace.
4. **Derive the complete boundary pairing.** The VC scalar theorem keeps
   `(pi_minus - pi_plus) dot h` for the normal from minus to plus. Match action
   sign and distinguish momentum flux from stiffness traction and from the
   S11b-derived energy current used in scattering. For a multidimensional
   application supply tangential integration, side regularity, common test
   traces, admissible variations and any surface action/source. Do not infer
   field continuity or the existence of arbitrary solution traces. Our smooth
   profile is not replaced by a discontinuous thin interface; artificial
   domain-decomposition cuts introduce no new physical surface action. The
   reduced local orders 0–3 and full nonlocal closed response need their actual
   Green/boundary pairing, beyond the first-gradient local VC families.

Before using these identities as a S11c operator/matching certificate, test the
actual mapped action and retained order against the profile-aware derivation,
including an omitted-gradient, derivative-index and reversed-interface-sign
control. Reuse the accepted operands; generic local identities alone do not
identify the nonlocal closure or its boundary data. A demonstrated upstream
disagreement would require the user's repair checkpoint; no new discrepancy or
production rerun is established by this reception. Preserve the running job's
sources and validate it on its original inventory when the watcher completes.

Preserve the full symbolic input and the builder report contract suffix. Keep
durable runs in repository _scratch/s11c/, publish validated .out files using
DataLad/git-annex, and commit each step. Use one CAS job and a silent local
completion/error watcher, with no recurring checks or model polling.
