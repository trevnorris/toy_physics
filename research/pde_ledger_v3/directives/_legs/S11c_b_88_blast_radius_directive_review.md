# Decision-list review — S11c-b #88 blast-radius BUILD DIRECTIVE

## Artifact (the thing you are reviewing)
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_88_blast_radius_build_directive.md`

This is an **orchestrator-written build directive** for a diagnostic SymPy instrument. You
are reviewing the **directive**, before any build runs. Your job is to catch defects in the
directive that would make the resulting instrument compute the wrong object, be ambiguous
to a builder, leak a withheld answer, or answer the wrong question — so that no build round
is wasted. A productive review finds real defects; "looks fine" is weak output.

## Context you need (read these in the repo)
- The engine the directive governs: `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
  (esp. lines cited by the directive: 616, 611-613, 936-969, 972-1003, 1025-1028, 1073-1090,
  1178-1215, 1228-1322, 1459-1511, 1850-1874, 2122-2137, 713-725).
- The DOF→row map: `research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py`
  `extract_slab` (~L760-799).
- The spec: `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` §1d (the
  total-divergence quotient does NOT lift to variable coefficients) and §3a.
- The settled #86 result the directive builds on:
  `research/pde_ledger_v3/directives/_measurements/S11c_b_86_reference_result.md`
  (corrected §3a basis = 40 = 10 uniform + 15 ∂W_bg + 15 ∂μ_R,bg; per-source 15 exact,
  nullity 0; reduces to the engine's committed frozen 26).

## What the #88 question actually is
"Does correcting the §3a basis (the frozen-spurion defect, engine count 26 → correct 40)
disturb the operator rows whose cross-engine verdicts were already adjudicated — i.e. would
completing the basis change any strong DOF row (u-momentum, θ, thickness) at the retained
grades, not only the coupling kernel?" The instrument must expose, per row, whether the
physically-correct operator differs from the emitted (frozen) one at retained grade — and
print enough structure (drivers + a non-absorbability rank) that the orchestrator can decide
whether a row's frozen-based verdict is invalidated. The instrument must NOT itself decide.

## Check, and report each as a numbered finding (file:line — problem — fix)
1. **Premise fidelity (verify against source; do not take the directive's word).** For each
   engine claim the directive makes — (a) the frozen quotient never differentiates the
   spurion so the strong rows freeze the Hessian (`:611-613`, `:616`, `:936-969`,
   `:1459-1511`); (b) the committed Hessian-retaining pattern is `operator_dx`/`background_dx`
   (`:1850-1874`, `:2122-2137`); (c) the retained-grade rule is η≤1 ∧ σ_W≤1 via
   `first_shape_series` (`:713-725`); (d) the DOF→row map — read the cited code and confirm
   it says what the directive says. Flag any misquote or line drift.
2. **Is the Hessian actually retained at the truncation, not dropped?** The directive asserts
   the second background derivative enters as `sigma_W · w1_profile_d{i}d{j} / L_W` (a single
   σ_W) so it survives η,σ_W ≤ 1. **Verify this in CAS** — write a short script that takes a
   term `grad_W[a]·(DOF)`, applies the Hessian-retaining dx, applies `first_shape_series`, and
   shows whether the Hessian term survives or is truncated. A prose argument is not accepted;
   save the script + literal stdout to named absolute paths and report them. If the Hessian is
   actually σ_W² (truncated), the whole instrument is measuring nothing — this is the highest-
   value thing to check.
3. **Is the named object the right one to answer #88?** The directive computes, per strong
   row, `EL_correct − EL_frozen` at retained grade, plus a rank-gain non-absorbability metric.
   Argue whether "the correct operator differs from the emitted operator on row X (non-
   absorbably)" correctly implies "row X's frozen-based cross-engine verdict is invalidated."
   Is PY-side correct-vs-frozen sufficient, or must the instrument also compute the WL side?
   (The directive claims PY-side suffices for invalidation because PY-frozen ≠ correct means
   the emitted PY row is wrong there regardless of WL.) Contest or confirm this logic.
4. **Non-absorbability metric adequacy.** The directive uses two signals: second-jet atoms
   absent from the frozen row, and a rank of DOF-monomial coefficient-vectors before/after
   adjoining the residual. Is the rank construction well-posed (are free coefficients correctly
   treated as scalars, not monomials; is the monomial basis the right space)? Could a residual
   be non-absorbable yet the rank test miss it, or absorbable yet the atom test falsely flag it?
   Propose the strongest single non-absorbability test if these are inadequate.
5. **Driver decomposition correctness.** DRIVER_B (Hessian on the selected 8) and DRIVER_C
   (omitted invariants under frozen dx) are meant to attribute the residual. Are these two
   drivers well-defined and non-overlapping? Does DRIVER_B + DRIVER_C reconstruct the total
   residual, or are there cross terms the directive misses (e.g. Hessian acting on the omitted
   invariants)? If the decomposition is incomplete, say what is missing.
6. **Builder-misreadable ambiguity (the #86 failure mode).** #86's reference build was
   defective because the directive's uniform-family carrier phrasing was read as per-source
   duplication. Hunt this directive for any phrasing a builder could implement two ways —
   especially: the completed density construction (which candidates, which coefficients), the
   live substitution vector per source, the μ_R,bg Hessian normalization, the monomial space
   for the rank. Name each ambiguity and give the disambiguating sentence.
7. **rule-5 / leak / fabrication discipline.** Does the directive leak the answer (the per-row
   disturbance, the settled family verdicts, an expected sign/count/nonzero-ness)? Is the
   acceptance genuinely "compute-and-print" with the only assert being the Control-A harness
   invariant, or does it smuggle a target the builder would iterate toward? Is Control A
   (freeze-dx + drop-omitted ⇒ residual 0) a real, non-tautological invariant (produced by an
   independent route — the un-corrected construction vs the engine's committed density), or is
   it zero by construction and therefore vacuous? Flag either failure.
8. **Shared-blind-spot / missing-control risk.** Is there a way the instrument returns a clean
   or misleading answer while the physics is wrong — e.g. a live-substitution that silently
   zeroes a whole family, a truncation applied in the wrong order, the frozen operator built
   with a subtly different density than the engine emits? Name the control that would catch it.

## Method
Read the cited engine source and the spec FIRST; form your own view of what the correct
variable-coefficient operator and its retained-grade truncation are; only then judge the
directive against it. Where a claim is checkable in CAS (finding 2 especially, and any premise
you doubt), **write the script and report its literal stdout** — a prose re-derivation is
discarded (this ledger's standing rule). Save scripts + stdout to named absolute paths under
your home or `/tmp` and report the paths. Copy anything you run to a scratch dir; do not modify
the working tree.

## Physics filter
Report a finding only if it catches a way the instrument would answer the #88 question wrongly,
be ambiguous to a builder, leak the withheld answer, or rest on a false engine premise. Do not
report style, or "it would be wrong on a different physics problem."
