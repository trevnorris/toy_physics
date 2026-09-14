# Mathematica c2 reference-pressure repair — awaiting user approval

The [source-pinned audit](S11c_wolfram_repair_audit_report.md) establishes a
duplicate displacement in `mathematica/S11c_c2_N6_mathematica_audit.wl`.
`faceLaws` constructs reference-face pressure plus its normal Taylor term.
`responseFamilies` constructs physical boundary pressure. `buildContraction`
and the `buildCase` closed-slot guard feed that physical response and its
normal continuation directly into the reference slots. The guard compares
two uses of the same image and therefore cannot establish physical trace meaning.

This plan is a reviewable proposal. The audit authorization does not authorize
editing that producer. REFERENCE normalization and further scattering pause at
this repair decision. The corrected SymPy chain remains accepted on its recorded
domains; this Wolfram finding does not itself require rerunning it.

1. Preserve the audited native source/output pins and failing witness. Derive
   the affine pressure-evaluation map from the actual native face law, keeping
   its reference-pressure and normal-jet meanings explicit. Compute its ordered
   inverse against the native physical response before closing reference slots.
   Retain physical-face pressure as its own operand; do not alter c1's pressure
   meaning or subtract the measured downstream residual.
2. Retain the complete eta/sigma rectangle, including the height map acting on
   first-shape response and both ordered three-leg products. Keep output, input
   and intermediate momenta distinct. Derive normal continuation from the
   reference-pressure field. Document all matching and inverse denominators;
   identities on their domain do not resolve thresholds or spectral defects.
3. Route the computed reference images through both `buildContraction` and
   the native closed-slot construction. Preserve each carrier/source/cross
   channel's meaning and external-work orientation. Add independent physical
   trace reconstruction controls to the existing algebraic closure guards.
   Extend the focused audit to nonconstant profiles, both coordinate routes,
   anchorings, density representatives and faces. Retain raw/discarded terms
   and test missing, doubled and reversed trace-map mutations at the input.
4. Develop on one bounded case, then regenerate the four-case native c2 N6
   output once after the repair and controls are ready. Check exact source
   census, grades, dimensions, channel reconstruction and affected fingerprints.
   Inspect the existing c2 ablation harness's source selectors for compatibility;
   refresh necessary native control evidence with explicit source pins. Preserve
   pre-repair artifacts. Cross-engine comparison/review remains a separate task.
5. Publish changed `.out` evidence through DataLad/git-annex, code and concise
   inventories through Git, and verify payload hashes after the checkpoint.
   Resume the pending SymPy REFERENCE normalization only after recording the
   repaired audit disposition. Then return to variable-profile matching and
   the remaining S11c-d program.

The earliest changed producer is the self-contained **Mathematica c2 N6** script.
Its main transcript, source-dependent c2 native controls and audit checkpoint
must be refreshed. No change to b/c1 actions or pressure solves is indicated;
their sources and transcripts should remain unchanged and hash-verified.
There is no Mathematica full-c2 export or d engine in this audited entry point.
The existing SymPy b/c1/c2 exports and d normalization artifacts do not consume
this Wolfram producer, so they need no value regeneration from this repair.
Historical cross-engine c2 conclusions cannot serve as clearance of new values.

Use `math -script <path>`, durable repository run directories and serial heavy
CAS execution. Check other sessions before each launch; the license ceiling is
two Mathematica scripts including internal kernels. Commit each substantive
step, with no push, scheduler, review legs, comparator or full Wolfram d build.
