# Full D4 constant-coefficient bulk contract

Status: **complete** for D4C.1–D4C.4. Recorded local verification passes and
both independent fidelity reviews are CLEAR. Authorized 2026-09-17 as the next
S11 increment. Governed by `../FORMALIZATION_POLICY.md`; prior D4 density
classification and D4B odd-density proofs remain unchanged. Review provenance
and dispositions are in [D4_BULK_FIDELITY_REVIEW.md](D4_BULK_FIDELITY_REVIEW.md).

## Claim and domain

For every constant real coefficient vector `v=(a,b,c,beta)`, use the supplied
action `L=-Q/2`, where

`Q(G)=a (tr G)^2+b tr(G^2)+c tr(G G^T)+beta P(G)` and
`P=F01 F23-F02 F13+F03 F12`, `F=G-G^T`.

The gradient convention is `G_ij=partial_i u_j`. Fields have four real
components on real spacetime with four spatial coordinates and one time
coordinate. Coefficients are unrestricted, including zero and negative values;
they are fixed parameters, not varied fields. Backgrounds are smooth and
variations smooth with compact support. The action has no time derivatives or
inertia. The relative action is the actual integral of the change in density.

## Bounded obligations

| Item | Deliverable and coverage |
|---|---|
| D4C.1 | Identify the full already classified SO(4)-invariant density, derive its actual jet momenta, and prove the actual first variation and equivalence of stationarity with the local EL equation. Reuse the S10 analytic lemmas and D4B odd action. |
| D4C.2 | Prove `EL=(a+b) grad(div u)+c Delta u` on all smooth fields. Prove equality of bulk operators iff both `c` and `a+b` agree. This covers the full family, not selected waves alone; waves may prove necessity. |
| D4C.3 | Prove the full variationally null family is `(t,-t,0,beta)`, with arbitrary `t,beta`. The response map `(a,b,c,beta)->(c,a+b)` has image dimension two and kernel dimension two. Retain explicit even and odd divergence currents and nonzero density/momentum witnesses. The compact homogeneous map is `rho=0,mu=c,B=a+b+c`; zero wavevector is retained. |
| D4C.4 | Check a compact exact native basis/action/EL identification, plus passing and deliberately false normalization, coefficient, count, null-family, current and modal controls. A failed instrument, warning, missing import or timeout is not mathematical rejection. |

The even current is `J_i=sum_j(u_i partial_j u_j-u_j partial_j u_i)`;
`div J=(div u)^2-tr(G^2)`. The odd current is the reviewed
`K_i=(1/2)sum_j u_j M_ij`, with `div K=P`.
Consequently a null density is a linear combination of these divergences.
Zero bulk variation does not imply pointwise zero density or no boundary effects.

## Fidelity, evidence and reviews

Reuse the reviewed D4 complete density classification (four-dimensional SO
space, three-dimensional O space, one-dimensional odd space) at `d6119f75`,
and the D4B actual odd first variation/current at `2eac8635`. The new proof
specializes the D3 even-family variation argument to four spatial dimensions
without editing reviewed sources. No new rotation certificate is needed.

The compact native check selects existing Q9/coordinate helpers without
importing the production module or running its driver. It identifies the
computed full D4 basis with the trace/orientation basis by an invertible exact
change of basis and compares the actual `L=-Q/2` EL normalization. An all-zero
odd bulk expression alone is insufficient: nonzero density/momentum controls
must retain the odd sign and factor. This translation check is not a Lean
kernel proof of the CAS software.

Recorded run2 freshly compiled 26 unchanged imports, six new modules and the
audit root. All 67 standard-axiom audits pass; sixteen paired mathematical
mutants reduce to `False` and twenty positive executions pass. There are
eighteen distinct positive source statements: two normalization positives
are deliberately repeated for separate sign/factor pairs. The native check
passes twelve identities, two mutated native EL implementations, and seven
target checks (three wrong-formula rejections and four positives). See
`D4_BULK_VERIFICATION.txt` and the author validation record for exact provenance. Build reuse requires
recorded transitive source/pin/generator, command and input/output object guards.
Historical sources and evidence remain unchanged; any freshly rebuilt imports
receive a new current object record. The preservation manifest is
`_measurements/S11_lean_d4_bulk_preserved_inputs.json` relative to the ledger.

The user explicitly approved the fixed 52-file packet. Claude and Grok each
returned a substantive CLEAR review with no required correction. Their source
reviews did not independently execute the builds or hash checks. The author
revalidated the packet, sources, dependency pins, instruments, logs and objects;
no shared-source drift was found. The closure adds only status and evidence
clarifications; the approved snapshot and all canonical proofs remain unchanged.
See `D4_BULK_FIDELITY_REVIEW.md` and
`_measurements/S11_lean_d4_bulk_closure_validation.json` relative to the ledger.

## Exclusions and stopping rule

No D5, variable coefficients, interfaces, general null-Lagrangian theorem,
new spectrum/stability/scattering/pole results, S11c calculations or exports,
production reruns or systematic CAS bridge. This is not a new physical discovery.

Stop when D4C.1–D4C.4, source/axiom checks, meaningful controls and both fidelity
reviews clear, and compact documentation records their scope and provenance.
No commit is requested for this increment.
