# Tail, Abel and inverse-stability contract T1–T4

Authorized by “Let's commit and move on to the next item (possibly one of the
questions the other session posed?)”, following the bounded proposal in
[S11C_ANALYTIC_ASSESSMENT.md](S11C_ANALYTIC_ASSESSMENT.md). Governed by
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). Status: bounded T1–T4 complete. Local verification
passes and Claude and Grok independently returned CLEAR. This is not a
numerical scattering certificate.
D4B is complete at `2eac8635`; preserve its historical proof and evidence.

| Item | Finite deliverable |
|---|---|
| T1 | For complex amplitudes on the real line, prove the norm bound for an omitted integral and the combined truncation/Abel error `B (tail + a firstMoment)`. State measurability, integrability, bounded multiplier and finite first-moment premises explicitly. Include arbitrary measurable kept domains, with `abs x <= R` as the physical specialization, zero regulator and zero amplitudes. |
| T2 | For Banach spaces and a supplied continuous linear equivalence A, prove that `norm(A^-1) <= kappa`, `norm E <= epsilon`, and `kappa epsilon < 1` give an actual inverse for A+E. Prove its norm, inverse-difference, source-to-solution and bounded-observation estimates. |
| T3 | Identify the actual native constant/Heaviside Abel transforms, Fourier measure, fixed origin and physical width a/L_W. Record the partition into local differential, ordinary integrable and distributional terms, and list unresolved full-operator application premises. All kernel, scattering-domain and inverse-stability bounds remain explicit obligations, not inferred from witnesses. |
| T4 | Include meaningful paired mathematical mutations and positive controls for omitted tails/moments, Abel scale/sign, lost conditioning factors and removed smallness. Check standard axioms and source/object correspondence; obtain two independent non-author fidelity reviews of a new fixed packet. |

The integral theorems may quantify over any measure on the real line; Lebesgue
measure is the physical application. This also permits finite nonzero controls
without unrelated integration machinery. No integral is used as evidence for
the action unless its integrability premises hold. The Abel factor is
`exp(-a abs xi)` in dimensionless position. The supplied multiplier is bounded
and measurable; a constant-plus-step is allowed. Keeping the folded amplitude's
first moment explicit avoids a claim of uniform L2 operator-norm convergence.

The inverse theorem concerns bounded maps between stated Banach spaces; the
actual differential operator requires its graph/outgoing domain to be supplied.
The estimate does not posit an unweighted L2 scattering resolvent. Coordinate,
momentum, quadrature, channel truncation and Abel errors must be controlled in
the same operator norm before applying the inverse result.

Coverage is the exhaustive kept/complement split for every measurable domain,
all admissible complex amplitudes and nonnegative regulator values, and every
operator perturbation satisfying the stated strict margin. Boundary/failing
margin cases have counterexamples, not inverse guarantees. No generic sampling
or count of native integrals establishes the missing application estimates.

The recorded run passed four modules plus the audit root, 41 standard-axiom
audits, twelve mathematical rejections and sixteen positive controls. Each
false statement produced exactly one unsolved `False` in `contract_control`;
no compiler/instrument failures were accepted. The native check passed eight
identities and four sign/scale/omission controls. See
[ANALYTIC_ERROR_VERIFICATION.txt](ANALYTIC_ERROR_VERIFICATION.txt) and
[ANALYTIC_ERROR_FIDELITY.md](ANALYTIC_ERROR_FIDELITY.md).

Both substantive independent reviews cleared the authorized fixed packet.
Nonblocking wording/count findings were resolved in documentation; proofs and
instruments are unchanged. Current source/object correspondence passes. See
[ANALYTIC_ERROR_FIDELITY_REVIEW.md](ANALYTIC_ERROR_FIDELITY_REVIEW.md).
The sixteen positive executions contain fifteen distinct statements. The
listed physical application premises remain open.

Stop after T1–T4 and the application-obligation record. No general
limiting-absorption theorem, Schur/pseudodifferential library, complete solver,
every-output CAS bridge, variable-coefficient/interface proof, nonlinear-pole
formalization, parent-theory bound or production rerun is included. Questions
2 and 3 remain queued. Preserve all S11c files/exports and other-session work.

Use one Lean worker, `-j1 -M4096`, sequential processes with 600-second limits
and process-group cleanup. Long verification uses a silent local completion
hook. The user authorized this exact review packet and subsequently requested a
commit after review clearance, then question 2 as a separate bounded increment.
That later instruction supersedes the original no-commit wording. It does not
authorize transfer of a future question-2 packet.
