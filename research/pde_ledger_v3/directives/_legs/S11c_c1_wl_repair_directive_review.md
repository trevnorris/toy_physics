# Directive design review — S11c-c1 Wolfram engine REPAIR directive (rule-7 gate, run retroactively)

## What you are reviewing

`research/pde_ledger_v3/directives/S11c_c1_wl_repair_directive.md` — an **orchestrator-written** repair directive
for the blind Wolfram engine of S11c-c1. It should have had these two decision legs **before** the repair build
ran; it did not (the orchestrator skipped the gate). The repair has since been built and 2-leg re-reviewed
(`_measurements/S11c_c1_wl_repair_directive.md`), so this review runs retroactively: your job is to judge whether
the **directive itself** is sound — faithful to the confirmed findings and the physics spec, leaking no expected
value, mis-specifying no fix, smuggling no new physics — and to flag anything that, if the directive was wrong,
could have propagated a defect the re-review might not have caught.

## Sources of truth — form your own view FIRST, then read the directive
- The confirmed defects the repair must fix: `directives/_measurements/S11c_c1_wl_build_review.md` (3 MUST-FIX,
  code-verified: R1 the `DTN_OPERATOR` composition freezing the input leg; R2 the energy bulk operand being the
  face quantity × a free `farFieldPhase` rather than a genuine far-field Poynting flux of φ; R3 the energy face
  operand using a reconstructed impermeable traction rather than the response `t_s`) + four NITs.
- The sole physics authority: `directives/S11c_c1_SHARED_PHYSICS.md` (§3a `:247-261` the two-momentum operator;
  §3b `:321-330` the three-object dissipation audit — the bulk operand is the outgoing Poynting flux of φ at
  `|w|→∞`, the face operand is the true-area traction pairing built from the §3b `t_s`; §6 the locus protocol +
  three script clauses). Siblings `S11c_a_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md` as needed.
- `CLAUDE.md` rules 2, 3 (name the object not the recipe), 5 (no expected value), 6, 11, 16, 17.

## The questions this review must answer (report a finding for any "no")
1. **Fidelity to the findings.** Does the repair directive's R1/R2/R3 each correctly and completely target the
   confirmed defect, at the CONSTRUCTION (re-entering at the composition operator / the bulk φ solution / the
   response `t_s`), never at a result? Does R2 correctly demand the bulk operand be the outgoing Poynting flux of
   the half-space φ at `|w|→∞` (spec §3b), with `farFieldPhase` removed — not a re-labelled face quantity? Does
   R3 correctly demand the face operand BIND the response `t_s` (carrying `Λ_X`)?
2. **Leak / rule 5.** Does the directive state, imply, or pre-register ANY computed value, sign, parity, residual,
   or expected result? Scrutinize the "construction invariant" probes (`COMPOSITION_HAS_Q_INPUT=True`,
   `FACE_HAS_TS_FROM_RESPONSE=True`, "a one-sided corruption moves only its operand"): are these structural
   invariants (legitimate), or do any of them encode an expected physics value?
3. **Name the object, not the recipe (rule 3).** Where the directive specifies a fix (especially R2's far-field
   Poynting from φ), does it name the OBJECT the emit must be and let the engine compute it, or does it over-
   specify a derivation recipe / under-specify so the builder could satisfy it wrongly? Is R2 specified tightly
   enough to be correct but not so tightly it dictates a particular non-unique construction?
4. **New physics / scope.** Does the directive smuggle any new physics, or correctly stay within fixing the 3
   findings + 4 NITs? Does it correctly protect the 2-leg-sound core (the kernel, flat symbol, T-a, inverse, μ_θ)
   as byte-identical, and correctly NOT touch it?
5. **The four NIT fixes (R4)** — are they correctly specified (real `SOURCE_EQUATIONS`; computed δW/ζc parity
   blocks; computed layer-potential route; genuine `REAL_ADMISSIBLE` test), each an object-not-recipe?
6. **Script clauses / method.** Does the directive correctly bind §6's three clauses (print not assert; emit
   A/B/A−B before any guard; no VERDICT; typed boolean not native), the structural rule, and the run discipline?
7. **Anything a builder could satisfy while producing a wrong or weak fix** — an under-specified invariant, a
   fix that re-enters at a result, a control the directive lets stay `A−A`.

## Method and output
This is a DOCUMENT review — quote the directive text and the finding/spec text it must honor (file:line) for
every claim (rule 2). State severity: MUST-FIX (the directive would produce a wrong/weak fix or leak an answer)
vs NIT. If the directive is sound, say so and name the two or three fixes you checked most closely. Since the
repair is already built, if you find a directive defect, also state whether the built code (per the repair record)
would have inherited it or whether the re-review would have caught it. Report only findings that catch a way the
physics or a control could be wrong.
