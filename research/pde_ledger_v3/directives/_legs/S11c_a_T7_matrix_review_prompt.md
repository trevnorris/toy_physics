# Independent review — S11c-a T7 feasibility/adjudication MATRIX

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_T7_adjudication_matrix.md`
(+ its measurements twin `…/directives/_measurements/S11c_a_T7_adjudication_matrix.md`).

## What this matrix is
Step 0 of a cross-engine reconciliation. Two independently-built engines compute the S11c-a interface
shape-derivatives; a comparator is supposed to join them. Before building it, this matrix **audits
where the two engines diverge in CASE STRUCTURE and semantics, and classifies each divergence** into a
structural bucket (AGREE-STRUCTURE / AXIS-DIVERGENCE / DECOMPOSITION / BESPOKE / CONSEQUENCE). It
deliberately does NOT decide which engine is physically correct (that is a later step) and asserts no
cross-engine result. Your job: is the classification **accurate, complete, and free of a mislabel that
would let a genuine physics disagreement be waved through as "serialization" — or the reverse**?

## Sources of truth (form your OWN view from these BEFORE trusting the matrix)
- PY tag stream: `/home/trevnorris/.s11_build/S11c_a_sympy_engine.out`
  (verify `sha256sum` == `6386471555b1e99d0aeb0f716eea30f839d59be50c0cedd4677ea7b376b79129`).
- WL tag stream: `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`
  (`sha256` `82062bd36cfb07b1f18631077f0c63ac1cbce7834967686f680fa9f30019e4ec`).
- PY engine `scripts/S11c_a_interface_geometry_sympy_audit.py`; WL engine
  `mathematica/S11c_a_interface_geometry_mathematica_audit.wl`.
- Spec (closed): `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — §4 (T-a…T-i laws),
  §5 (controls), §7 (emit grammar). Use it to judge whether a divergence is a case the spec REQUIRES,
  FORBIDS, or leaves open — but do NOT try to finalise which engine is right; only judge the matrix's
  STRUCTURAL classification.

## Required method — COMPUTE, do not read
⛔ A prose re-derivation is worth nothing. **Write your OWN census script** (Python; parse both `.out`
streams; do NOT run the orchestrator's `_measurements/s11ca_t7_census/s11ca_axis_census.py` as your
check — you may inspect it afterwards to spot a parser bug, but your independent number must come from
code you wrote). Save your script AND its literal stdout to named absolute paths and report those paths;
without them your findings are discarded. ⛔ Write your census script + any scratch UNDER
`/home/trevnorris/.s11_build/` (outside the repo) — do NOT create files in the working tree, and do
NOT modify either engine or any committed file.

Your census must, per joined tag, compute each engine's **leaf-case count and case-axis set** (flatten
PY's quantity-nesting to leaves; WL flattens axes into the pipe-string key). Then:
1. **Reproduce or refute** the matrix's per-tag counts and axis-sets (§ "AXIS-SET census" table /
   `axis_census.stdout`). Report every tag where your independent census disagrees with the matrix.
2. **Check every bucket assignment.** In particular hunt for:
   - a tag classified **AGREE-STRUCTURE** (Family A) whose stored LEAF representation actually differs
     (Family F: PY graded `(background, ε·deriv, total)` tuple + `epsilon_shape` vs WL derivative
     COEFFICIENT under `EXPRESSION` with order in `MULTIGRADE_EPSILON_ETA_SIGMAW`). Sample the leaf of
     each Family-A tag on BOTH sides and say whether the case grid agreeing hides a content divergence.
   - a **DECOMPOSITION** (FIELD/ORIGIN itemisation, Families D/E) that is NOT in fact a lossless
     decomposition of the other engine's aggregate — i.e. a real extra-content finding mislabeled as
     candidate serialization.
   - an **AXIS-DIVERGENCE** (Family B density, Family C virtual-dof) that the matrix under- or
     over-states: is the density axis really absent on the side claimed? is WL's virtual-work really
     the full physical×virtual matrix and PY really the diagonal? Confirm at the producer lines
     (PY `:916`,`:919-924`,`:1466`; WL `:1039-1051`,`:846-874`).
   - a **CONSEQUENCE** (control family H) that is actually an INDEPENDENT divergence, not merely
     inherited from the primaries.
   - any **joined tag missing from the matrix**, or any bucket that should be split.
3. **Do NOT run the step-1 tests.** The matrix names a RHO4-vs-RHOBR residual and a virtual-work spec
   question as step-1 work; you are reviewing the MAP, not resolving the physics. If you think the
   matrix mis-NAMES a test or its decidability, say so — but do not compute the residual yourself.

## Physics filter
Report a finding only if it changes what the matrix CLASSIFIES or what a downstream step may CLAIM — a
mislabeled divergence, a wrong count that moves a tag between buckets, a missed tag, a leaf-semantics
divergence hidden under an agreeing grid, or a rule-2/rule-5 defect. Do not report cosmetic wording.

## Rule-2 / rule-5 audit
- rule 2: does any matrix claim (a count, an axis-set, a producer line) lack a command behind it in the
  twin? Flag it.
- rule 5: does the matrix (or twin) leak an EXPECTED cross-engine RESULT — i.e. pre-decide which engine
  is correct, or state what the density/virtual-work residual WILL be? It must not. Flag any leak.

## Output
A short verdict: is the matrix sound and complete as a step-0 classification, or the exact rows/buckets
to fix, each with your independent census number and the `.out`/source evidence. Name the absolute
paths of your census script + its stdout. Focus on catching a genuine physics disagreement mislabeled
as serialization (or an agreeing grid hiding a leaf-semantics divergence) — that is the failure that
would corrupt the whole reconciliation.
