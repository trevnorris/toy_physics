# S11c Mathematica audit of the four SymPy upstream repairs

The user approved queuing this audit on 2026-09-14 and specified a license
limit of **two concurrent Mathematica scripts**. This is a repository work
queue, with no scheduled task or automatic restart. The new authorization
covers focused Mathematica checks of these four repairs; it supersedes the
earlier SymPy-only lane restriction for this bounded audit. It does not
authorize a new upstream repair, a comparator, review legs, or the full
S11c-d Wolfram build. The retained solver/export contract stays unchanged.

Finish validation and checkpointing of the already-completed LEFT normalization
first. Run this audit before starting REFERENCE normalization or further
variable-profile scattering work. A confirmed defect requires a concrete
repair and dependency-regeneration plan and a stop for user approval. If the
audit is unresolved, report the exact remaining test and the implication for
continuation; do not silently turn an inconclusive check into clearance.

The history check found no changes to the corresponding S11c Mathematica
sources or transcripts between the parent of a74da30a and c870c32a. The four
SymPy repair records explicitly exclude Wolfram runs. This establishes an
unclosed audit obligation, not the presence of any of these bugs in Mathematica.
Earlier cross-engine checks do not establish agreement with the repaired chain.

| Issue and SymPy checkpoint | Independent Mathematica check to construct |
| --- | --- |
| Conservative inertia sign — a74da30a | Derive the time variation of the supplied kinetic action and compare it with the native assembled mechanical rows in their declared sign convention, including the material-constraint route. Keep the stored-energy and kinetic operands separate. |
| Mechanical-load sign/normalization — c643112a | Derive the external-work contribution in the same action convention as the stored mechanical row; compare each face's native load routing and its power pairing. Preserve physical generalized force as a distinct operand. Test source mutations that can detect a routing/sign error. |
| Reference-face versus physical-face pressure trace — a05b05e3 | Trace the native c1 pressure operand into c2 closure. Construct the affine face evaluation and its inverse for reference pressure independently, with the normal jet and ordered two/three-leg products retained across the full eta/sigma rectangle. Check for missing or duplicate face displacement. |
| Thickness coordinate in the kinetic action — 3b52afcb | Compose the declared physical displacement from e_W = deltaW/W_0 before time variation. Independently compare the native action and resulting inertia across both density representatives and nonuniform backgrounds. Equal-frequency or reference-only checks do not cover the original failure. |

Start from the actual Wolfram entry points and trace the helpers/data they use:

- `mathematica/S11c_b_brane_operator_mathematica_audit.wl`
- `mathematica/S11c_c1_bulk_closure_mathematica_audit.wl`
- `mathematica/S11c_c2_N6_mathematica_audit.wl`

Inventory their precise construction and output scopes before choosing a
harness. A similarly named tag is not evidence that it covers the repaired
assembly. Pin sources, relevant authorities and existing outputs before any
calculation. Read the corresponding repair reports for failure locations and
coverage requirements, but derive checks from the shared action/trace
definitions and native Wolfram operands. Do not transplant corrected SymPy
expressions as Mathematica answers or treat an expression diff as the audit.

Develop one bounded case using the existing producers/helpers, then cover both
anchorings and density representatives where the property depends on them.
Prefer symbolic coefficient/rank identities where practical. Otherwise state
each tested component, grade, parameter domain and exceptional denominator
locus explicitly. Do not let one generic witness establish a union of loci,
full degenerate subspace, or exceptional-stratum claim. Emit both independent
operands and their residual, with dimensions and perturbative order. Preserve
raw, retained and discarded terms separately.

Use focused extractions and small independent variations before committing to
full regeneration. Record each issue as PRESENT, ABSENT_ON_AUDITED_DOMAIN, or
UNRESOLVED, with evidence and limitations. The initial queue status is
NOT_YET_AUDITED. An inability to compute a required object is UNRESOLVED; never
substitute an empty result. If a discrepancy is found, identify the earliest
affected producer and which descendants need changed values versus refreshed
source pins. Report and stop before modifying that producer.

Launch each audit with the user-specified command `math -script <path>`, where
`<path>` is the actual Wolfram script path. Record that exact invocation with
the run's durable logs and source pins.

Execution is serial by default. The hard maximum is two concurrently running
Mathematica scripts across sessions, including scripts outside this audit.
Check occupied Wolfram processes/license use before launch. Internal parallel
kernels must not bypass this limit. The separate memory discipline still
forbids overlapping heavy CAS jobs, including SymPy and Mathematica; two
license slots are not permission to run two heavy jobs together. Use this
repository's `_scratch/s11c/` for durable logs, source snapshots and caches.
Do not terminate another session's kernel to acquire a slot.

Deliver focused Wolfram audit code, bounded `.out` evidence under
`mathematica/out/`, a source-pinned per-issue checkpoint, and a concise report
under `_measurements/`. Put published `.out` files in DataLad/git-annex and
ordinary code/plans/inventories in Git; verify payload hashes after saving.
Commit each substantive step under the user's standing authorization. Preserve
failed attempts and historical outputs. No push or automation is part of this
queue. Update the next-work pointers after the audit disposition so the
REFERENCE/scattering continuation cannot bypass an unresolved repair decision.
