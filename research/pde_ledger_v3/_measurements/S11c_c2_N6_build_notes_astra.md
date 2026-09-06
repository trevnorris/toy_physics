S11c-c2 N6 build transcript

Deliverable: /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py
Transcript: /tmp/S11c_c2_N6/transcript.jsonl
Backend checks: /tmp/S11c_c2_N6/backend_tests.log
Triage log: /tmp/S11c_c2_N6/derived_or_declared.log

SUPPLIED, and not verified within this build: the §5c object definition;
§3a close; §3c weak increment, slot-linearity and extract; the imported c1
face response; the S11c-a/b material builders. This computes N6 given them.

WITHHELD: every expected value or acceptance criterion for the N6 residual
and the independence controls. The computed differences are for external
adjudication. No residual value controls the diagnostic's exit status.

Construction: Eulerian pressure coefficients come from the imported slab.
Material coefficients come from a native material face-only adapter. Material
mu is varied termwise from the pulled-back quadratic scalar; its constitutive
alias binds the face circuit and its amplitude binds the imported c1 source.
There is no transformation of a carrier, increment, or weak row. Both routes
use the same imported response-kernel circuit and one Eulerian extract.

Numeric table columns identify weak block, retained grade, integral signature,
and face (0 denotes the computed sum of the two faces). Each entry is a pair
[numerator, denominator] modulo the corresponding sample prime. Sample rows
join the preceding N6_SAMPLE_PROVENANCE catalog. Metadata gives independent
unit inference, target units, epsilon support, degree bounds, and the exact
conditional false-negative bound. The bound is conditional on a good selected
prime for the tested nonzero coefficients and the stated hash-sampler model;
it concerns the declared formal weak-kernel coefficients. No integral or
Fourier phase is evaluated.

The final run is the baseline invocation made by the available legacy
reduction/derived_or_declared.py, whose copy is in the local review baseline
(the current checkout has no reduction/derived_or_declared.py). The outside-
repo triage_runner.py streams the child output and disables the tool's internal
timeout. It leaves the classification and parser logic unchanged. The legacy
parser accepts colon-tagged text; its handling of the required JSONL stream
is recorded in derived_or_declared.log.

Computation locations in the deliverable:
- constitutive(): quadratic source pullback and termwise mu variation;
- face_factory(): native face adapter, mu alias binding, normal-only injection;
- source_terms()/kernel_coefficients(): imported c1 source and ordered kernel;
- template()/build_increment(): carrier-first, single weak extraction;
- slot_guard()/closure_guard(): independent full-circuit and carrier evaluation;
- residual(): literal differences of the separately compiled operands;
- pit(): exact modular numerators/denominators, branch samples, bounds, emission.
