# Independent physics review — S11b SymPy interface-coupling-law engine

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py`
— a SymPy engine. Running it flushes a tagged stdout stream AND writes a second product,
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_exports.py` (the accumulated chain record
the next step imports). Both products are in scope.

## What to check
The engine computes the S11b **brane–bulk interface coupling law** and emits its objects, then serializes an
export record. Two independent burdens, and you must discharge BOTH:

**(A) The physics is genuinely COMPUTED, not typed.** Every emitted physical result must be reached by a
computation in the file — a solved root, a symbolic determinant, a discriminant, a passivity/energy form
built from the action and ansatz. A payload that is a hand-authored CAS object with no data dependence on
the derivation is the exact defect this review exists to catch, and it is invisible to code-reading.

**(B) The export-chain wiring is correct.** Specifically:
- imports `LEDGER` from `scripts/S11_exports.py` and carries it forward flat;
- binds `c_s0`→`LEDGER['c_s0']`, `mu_R`→`LEDGER['mu_R']`, `rho_br`→`LEDGER['rho_br']` to the IMPORTED
  objects — ⛔ never re-declared under a new identity, ⛔ no second inertia knob minted for `rho_br0`;
- `F9` collision handling is **three-valued** and its comparison is **TOTAL over the imported row shapes**
  (elementwise into containers — tuples/booleans/strings/relations), so an `F9b` equal-merge is decided by
  object comparison, ⛔ never by a subtraction residual that would raise on non-expression rows and silently
  fall through; an `F9c` write uses the `s11b_` prefix and leaves the imported row untouched;
- `BUILD_INPUT_DIGESTS` sha256-pins EXACTLY three inputs: this engine's own source, `S11_exports.py`,
  `S11b_SHARED_PHYSICS.md`;
- the `D3` in-run round-trip reconstructs the written module and compares against the live object;
- `_RELATIONALS` reviver is present and `_restore` revives relationals/inequalities;
- the export is frozen (`MappingProxyType` outer + per-record, then `del _LEDGER`);
- ⛔ NO verdict: no `PASS`/`FAIL`, no `assert residual == 0` standing in for a printed residual, no prose
  conclusion as an `emit` payload. A boolean must be a symbolic test's OBJECT, never a typed literal.

## What you are handed
- The engine script and its second product `S11b_exports.py` (both paths above).
- **The physics source of truth:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md`
  — the spec whose §§0–13 define the setup, field content, governing equations, ansätze, supplied laws,
  derivation routes and tag obligations. This is what the engine is supposed to compute.
- The chain-rule contracts (for the wiring burden, ⛔ not physics answers): `F1`/`F9` at
  `research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` (F9 at :171-186, :157, :203);
  `D1`/`D3` at `research/pde_ledger_v3/S9_REWRITE_PLAN.md:210,:217`; class vocabulary at
  `S9_REWRITE_PLAN.md:41`.
- `CLAUDE.md` governs (esp. rule 2: a script PRINTS operands + residual, never asserts a conclusion).
- ⚠ The build directive `directives/S11b_sympy_build_directive.md` exists, but your job is **independent
  verification of the physics AND the wiring**, ⛔ NOT directive box-ticking — a script can satisfy its
  directive and still type a payload or misread the spec. Derive the physics yourself.

## Required method — this is a SCRIPT
Derive independently. Write your OWN derivation script for the load-bearing objects (a coupling coefficient,
a passivity/energy form, a dispersion root) **before** opening the artifact, and save BOTH the script and its
literal stdout to named absolute paths — report those paths. ⛔ A prose "I re-derived it and got X" is a
typed sentence with no computation behind it and will be discarded.

⛔⛔ **A FORM ABLATION IS MANDATORY, ⛔ NOT OPTIONAL — it is the only thing that has ever caught the worst
defect.** Change the STRUCTURE of a load-bearing object in a copy of the engine — flip a sign AND an
off-diagonal, collapse two independent symbols into one, alter a governing equation's form — then re-run and
report the LITERAL diff of the emitted tags. A COEFFICIENT rescale tests arithmetic only; only a FORM change
tests physics. If the ablation produces byte-identical output for a tag, that tag was TYPED, not computed —
report it by name. Measured precedent: a script `emit`ted a physics conclusion as a typed sentence and its
ablation output was byte-identical; eight fidelity legs missed it.

For every emitted claim, ask **WHICH LINE COMPUTED THIS?** — give the line number or report it as uncomputed.

**Two-route / one-sided corruption:** where the engine claims two independent routes to one object and emits
their residual, a zero residual proves nothing alone. Corrupt ONE route at a time and report which objects
move; if breaking route A also moves route B's object, they were never independent and the residual is
decoration. ⛔ Do not verify a supplied object against the artifact's own identity — that is circular.

Report any `assert` that PRECEDES the value it guards — an assert-before-emit hides a form defect because a
perturbation strong enough to flip the check kills the process, so you see only PASS-or-crash.

## Physics filter
Report a finding only if it catches a way the PHYSICS could be wrong (a typed payload, a misread governing
equation, a false-equal merge that collapses two distinct quantities, a circular/tautological residual, a
wiring defect that corrupts the exported objects). ⛔ Do NOT report "the script would be wrong on a
different input."

## Ablation sandbox
Copy the engine to `/tmp` and ablate the COPY. ⛔ Never modify the working tree. Save every ablation script
AND its literal stdout to named absolute paths, and report those paths. (This is pure Python — no CAS kernel,
no licence seat, no timeout constraint.)
