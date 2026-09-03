# Independent physics review — S11c-c SHARED PHYSICS spec (decision leg)

You are one of two independent reviewers of a **physics specification**. This spec will be read by two blind
computer-algebra engines (a SymPy engine and a from-scratch Wolfram engine) that each build the S11c-c physics
from it. An error in this spec is the worst kind: both engines would faithfully compute the same wrong thing and
their agreement would hide it. Your job is to find every such error **before** any engine is built.

⛔ This is a DOCUMENT review, not a script review. There is no code to run. But a prose claim that you
"re-derived X and it works" is worth nothing unless you show the computation. Where a question is settled by a
short symbolic derivation, **write a small script (SymPy or by hand with explicit algebra), save it and its
literal output to a named absolute path, and cite the path** — otherwise your derivation claim is discarded.
`sympy` is available. Save any script under `/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/`
or `~/.s11_build/`.

## The artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c_SHARED_PHYSICS.md`

## What S11c-c is (in one paragraph, so you can judge completeness)
S11c is a staged family computing the non-uniform (curved-brane) transverse↔thickness coupling. S11c-a produced
the tilted-face interface shape derivatives (T-a..T-i); S11c-b produced the variable-coefficient slab operator and
its off-diagonal coupling kernel, but left the bulk pressure `δp_s` and the flat-face response kernels `Λ_I(ω)`
**symbolic** in every face/flux slot (it did NOT solve the bulk). **S11c-c closes the bulk**: solve the perturbed
curved two-face outgoing bulk acoustic problem to get the nonlocal Dirichlet-to-Neumann / impedance operator
(`δp_s = Z·v_bulk,s`), compose it with the interfacial mass balance and face closure (S11b's B0c, on curved
faces), and fold the eliminated-bulk result back into S11c-b's slab operator to produce the coupled nonlocal
self-energy operator. This is the curved-face generalization of S11b's flat-face B0b (`Z`, three regimes,
per-face inertial loading) and B0c (permeable face response).

## The sources of truth — read these and form YOUR OWN view of what S11c-c must require, BEFORE opening the spec
Read, in this order, and build your own list of the objects, setup, controls, and hazards a correct S11c-c spec
must contain. THEN read the artifact and compare. (Reading order is method, not a blindness control — you have
all of it; use it to derive independently first.)
1. `research/pde_ledger_v3/directives/S11c_decisions.md` — the ratified family decision list (N-series). N1
   (staged family / chain wiring / own comparator+exports), N5 (⛔ no global dispersion), N8 (T7 + reconciliation
   schema inherited), **N11a** (`v_bulk_normal_0` inherited as a standing rest-frame limit; every result
   conditional on `|q·v_bulk_normal_0/ω|≪1`; ⛔ the convective bulk operator is NOT an S11c task), N12
   (multigrade), N13/N7 (confinement/falsification are S11c-e), N14 (fresh names).
2. `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` §1, §1b, §2, and §9 tasks B0a/B0b/B0c — the
   FLAT-face bulk closure S11c-c generalizes: the impedance definition `Z ≡ δp_face / (bulk OUTWARD normal
   velocity)`, the radiation condition, the branch object `q_out` and its sound-cone branch points, the three
   regimes (bulk normal wavenumber² positive/negative/zero=grazing), per-face inertial loading, the permeable
   face response, the degenerate-loci LOCUS PROTOCOL (§8), and the two measured traps (a real part of `Z` in the
   propagating regime is bulk RADIATION, not interfacial transfer; branch re-selection can turn a sink into a
   source).
3. `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — the tilted-face substrate. In particular T-i
   (`S11CA_CLOSURE_SHAPE_DERIV`, §4) which explicitly EXCLUDES the bulk solve ("NOT B0c's `δp=Z·v_bulk`; no bulk
   DtN, impedance, or pressure-response solve belongs in T-i") — i.e. the seam S11c-c fills; and T-a/T-a′/T-a″/
   T-b/T-c/T-c′/T-d/T-e that S11c-c consumes.
4. `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` §0, §1c, §3b — the slab operator with symbolic
   `δp_s`/`Λ` face slots that S11c-c closes; the θ-row = `evolution_mass_balance − Σ closure_shape_deriv`.
5. `research/pde_ledger_v3/steps/S11c_a_interface_shape_derivatives.md` and
   `research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md` — what the two prior sub-steps actually
   established (per-engine-verified; S11c-b's cross-engine residual is DEFERRED — S11c-c imports it as such), and
   the pinned reconciliation schema (S11c-a record, "reconciliation schema" section).

## What to check — report a finding only where it catches a way the physics could be wrong
Derive your own answer first, then check the spec against it. In particular:

1. **Leaked answers (rule 5).** Does the spec state any value, coefficient, sign, regime outcome, parity outcome,
   branch selection, or cancellation that the engine is supposed to COMPUTE? The spec is meant to be a pure
   obligation-to-compute (there is no withheld acceptance value in S11c-c). Two SUPPLIED interpretation caveats
   are inherited from S11b (Re Z in the propagating regime = bulk radiation; the per-face inertial-loading sign
   convention). Judge whether each is a legitimate inherited SETUP caveat (as S11b's own spec carried) or whether
   it leaks the answer for the S11c-c curved object. Flag any true leak.

2. **Frozen varying quantities (rule 17 — the recurring root cause).** The bulk normal wavenumber / branch object
   `q_out`, the three-regime discriminant, the face tilt, and every background coefficient VARY. Does the spec
   keep them LIVE and differentiated, or does it (anywhere) let an engine freeze one to proceed? Is the §5d
   regime/branch-liveness control adequate and non-tautological? Is there any place a "required freeze" is
   smuggled in as a step rather than surfaced as a finding? Is the rest-frame limit (`N11a`) handled as a scope
   CONDITION (not a frozen term), and is the exclusion of the convective operator correct — or does closing the
   bulk actually REQUIRE the convective/drain term at the order requested (i.e. is dropping `v_bulk_normal_0`
   from the operator legitimate at first shape order, given the smallness domain)?

3. **Recipe vs object (rule 3).** Does the spec NAME the objects (the curved DtN operator, the permeable
   response, the self-energy operator) and fix the expansion order, or does it over-specify a derivation PATH
   (a particular solve method) that manufactures a well-definedness question? Conversely, is it UNDER-specified
   anywhere — a supplied law missing, a coupling that cannot appear because the setup omitted it (the `∇·u`
   fell-out-four-times failure mode)?

4. **Correct consumption of the substrate.** Is the "seam" identified correctly — is it true that S11c-b left
   `δp_s`/`Λ` symbolic in the θ-row `closure_shape_deriv` and elsewhere, so §3c's fold-back is the right
   assembly? Are the S11c-a T-objects consumed correctly (normal, conormal, shifted trace, kinematic balance)?
   Check the two geometric claims: (a) "the two exterior half-spaces do not share a bulk region (the slab
   interior is not bulk)"; (b) "the two faces tilt oppositely at first shape order (`∂h_s/∂x = (s/2)∂W_bg`), so
   the curvature correction carries a definite parity under `s→−s`." Are these right? Does face parity actually
   couple the two half-spaces, or only structure the (δW, ζ_c) response? Is the impedance still per-face?

5. **The DtN/self-energy assembly.** Is `δp_s = Z·v_bulk,s` the correct bulk relation to close, with
   `v_bulk,s` the bulk OUTWARD normal velocity at the curved face (T-a normal, T-a′ conormal)? Is the
   composition with `v_bulk,s = V_s + J_s/ρ_m` and the §1d closure the correct B0c generalization? Does folding
   the closed loads into `S11CB_SLAB_OPERATOR`'s marked face slots correctly yield the self-energy, and is the
   instruction to touch ONLY the term-origin-marked face/flux/closure slots (not re-derive or alter non-face
   rows) sound? Is anything DOUBLE-COUNTED (e.g. a face load counted both in T-i's closure shape-derivative that
   S11c-b already folded AND again in §3c)?

6. **N-series compliance.** N1 (own spec/two blind engines/frozen T7 comparator/`S11c_c_exports.py`; SymPy
   chains from `S11c_b_exports.py`, WL blind), N5 (⛔ no `ω(k)` — is the DtN emitted as an operator/kernel, not a
   spectrum?), N8 (T7 + the S11c-a reconciliation schema; ⛔ no pre-registered fold — and a branch/regime the
   engines key differently must be a computed residual, not a fold), N11a (above), N12 (every object multigraded
   `(ε,η,σ_W)`; is the truncation "first order in ε, first shape order in η and σ_W" the right order to capture
   the leading curvature correction?), N13/N7 (confinement + falsification correctly deferred to S11c-e), N14
   (fresh `S11CC_*` names, no imported-key reuse).

7. **Control adequacy.** Are the §5 controls genuine and non-tautological? Specifically: does the §5a two-route
   representation-invariance (direct Eulerian boundary-perturbation vs material face-flattening) give two
   GENUINELY INDEPENDENT constructions of the SAME object (so their residual is meaningful), and is the
   one-sided corruption a real independence test? Does §5b form-ablation test physics (form) not arithmetic
   (coefficient)? Is the §5c uniform-limit correctly scoped as a REGRESSION that cannot validate the
   first-shape-order curvature (which vanishes at η,σ_W→0)? Do the §5a/§5b/§5d controls actually have TEETH on
   the first-shape-order curvature and on the three-regime/branch structure — the physics the uniform limit
   cannot reach?

8. **Tractability / split (N2).** S11b took eleven directive revisions and S11c-a's comparator needed delegation
   because the surface was too large. Is the S11c-c surface (curved-bulk DtN solve + permeable response + loci +
   dissipation audit + self-energy fold, over 2 anchorings × 2 densities × 2 faces × 3 regimes × 2 parities) too
   large for one spec + one engine pair to build and independently review? If so, propose the finer split (which
   objects form a reviewable sub-unit) — but only if you can name a concrete reason the current single-spec
   scope is unreviewable, not merely large.

9. **Scope leaks.** Anything nonlinear, spectral, or leakage/confinement-related wrongly IN scope? Anything the
   curved-bulk closure genuinely NEEDS that is missing or wrongly deferred?

## Physics filter
Report a finding only if it catches a way the physics could be wrong, a leaked answer, a frozen varying quantity,
a missing/incorrectly-consumed input, a tautological or toothless control, or a genuine tractability/split
blocker. Do NOT report style, wording, or "an engine could misread this" unless a concrete misread produces wrong
physics. For each finding: quote the spec line, quote the source-of-truth line it contradicts (or the derivation
that shows it wrong), and state the concrete wrong-physics consequence. If you believe the spec is correct and
complete on a given axis, say so explicitly — a clean axis is a finding too.

## Output
A ranked list of findings, most-severe first, each with: the spec location, the contradicted source (or your
saved derivation path), the wrong-physics consequence, and a concrete proposed correction. End with an overall
verdict: is the spec safe to build two blind engines against as-is, safe after the listed folds, or does it need
a re-author / finer split?
