# Independent physics review — S11c-b SymPy engine (a SCRIPT)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
(repo root `/var/projects/toy_physics`).

A Codex-built SymPy engine. It constructs the S11c-b variable-coefficient brane (slab) operator and the
off-diagonal transverse→{θ,e_W,u_L} coupling kernel, and emits `PY_S11CB_*` tags. It chains on
`scripts/S11c_a_exports.py` (imports the accumulated `LEDGER`) and writes `scripts/S11c_b_exports.py`.

## What you are handed
- The script above (the artifact).
- The physics authority to derive against and check faithfulness to:
  `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (sole authority),
  and the specs it inherits by reference (`S11c_a_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md`).
- ⛔ You are NOT handed the sibling Wolfram engine's output, and must not seek it — the cross-engine check is
  a separate downstream comparator; your job is whether THIS engine faithfully computes the spec.
- There is no expected numeric answer: the spec withholds every computed value (the coupling grade/sign, the
  new invariants, the admissibility residual, the basis count). ⛔ Do not treat any value in the script as a
  target to confirm.

## Required method — this is a SCRIPT; derive independently and ABLATE
Write your OWN small derivation script for the load-bearing physics BEFORE opening the artifact in depth, and
save BOTH your script and its literal stdout to named absolute paths (report them). A prose re-derivation is
discarded (a typed sentence with no computation behind it is the exact defect this rebuild removes).

⛔⛔ **A FORM ABLATION IS MANDATORY, not optional — it is the only thing that has ever caught the worst
defect.** Change the STRUCTURE of a load-bearing object in a `/tmp` COPY of the script and re-run the affected
task, reporting the LITERAL diff. A COEFFICIENT rescale tests arithmetic; only a FORM change tests physics.
Do at least these ablations (each on a /tmp copy, re-running only the affected task, literal stdout):

1. **The coupling is gradient-driven.** Set a background FIRST JET to zero at its source (e.g. zero the
   supplied `∂W_bg`/`σ_W·∂ξ w₁` map, or `∂μ_R,bg`) and re-derive the coupling kernel. If the off-diagonal
   `transverse→{θ,e_W,u_L}` block does NOT collapse when the background gradient is removed, the kernel is not
   actually built from the gradient (it may be a hand-typed or substitution artifact). ⚠ `∇W_bg→0` is the
   engine's own uniform limit — use it to test that the kernel VANISHES there (a forbidden gradient-
   INDEPENDENT surviving term is a real defect), NOT as a validation of the coefficient.
2. **The energy is CONSTRUCTED, not substituted.** Check the variable-coefficient energy basis (§3a): is it
   built from the symmetry group with the background jets as spurions, EMITTING new gradient-of-background
   invariants — or is it just `W_0→W_bg` substituted into the uniform `U` (which would omit the undifferentiated-
   `u` spurion couplings, spec §1a/§1c/`N15`)? Ablate: remove one spurion channel from the enumeration and
   confirm a corresponding new-invariant / kernel term disappears; if nothing moves, the basis was substituted.
3. **Admissibility is the background-order (ε⁰) balance, not ε→0 of the wave operator.** Confirm
   `PY_S11CB_ADMISSIBILITY_OPERATOR_OPERAND` is the background-order generalized force (body force + per-face
   traction) sourced by the profile, in the same pairing as `𝒮_hold⁰` (spec §2b/§3d) — NOT the wave operator
   at zero perturbations (which is identically 0, making the test vacuous). Ablate: if the operand is
   structurally the `ε→0` limit of the §3b operator, it will be identically zero — report that.
4. **The off-diagonal block is EXTRACTED from the operator, not asserted.** Verify the kernel is read off the
   divergence-form operator's curl↔div structure (§3c), not hand-typed. Delete the operator-construction and
   check the kernel emission breaks (a hand-typed kernel would survive).

Also report, for the general script defects:
- Any `assert` that PRECEDES the value it guards (an assert-before-emit hides a form ablation — a perturbation
  that flips the check kills the process, so the leg sees only PASS-or-crash). Report every such assert.
- Any emitted payload that is a HAND-TYPED CAS object with no data dependency on the derivation (delete the
  `Solve`/`diff`/block-extraction and the output does not move) — the "hand-typed CAS object is still
  hand-typed" defect.
- Any TAUTOLOGICAL residual: `A − B` that is zero by construction for any input (e.g. `q:=A/B` then `A − q·B`).
- Any conclusion emitted as an unconditional literal, or a check whose expected value lives inside the artifact.
- Any VERDICT/PASS/FAIL terminal token, or a native boolean serialized as a residual operand (§6 forbids both).
- For each emitted §4 object, ask: WHICH LINE COMPUTED THIS? Give the line number or report it as uncomputed.
  Cross-check with `reduction/derived_or_declared.py` on the deliverable (triage, not a verdict).

## Physics filter
Report a finding only if it catches a way the physics could be WRONG (a coupling that isn't gradient-built,
a substituted energy missing the N15 channels, a vacuous admissibility test, a hand-typed/tautological
emission, an assert hiding an ablation). Do not report "the script would be wrong on a different input." A
leg that finds nothing is weak evidence — say which load-bearing objects you ablated and what the literal
diffs were.

## Operational constraints (IN THIS PROMPT, both legs identical)
- ⛔ Copy the script to `/tmp` and ablate the COPY. ⛔ Never modify the working tree.
- The engine writes `S11c_b_exports.py` (26MB) on a full run — ⛔ do not run the full package repeatedly;
  re-run only the affected task for an ablation, and direct any run to a /tmp copy / scratch path so the
  repo's `scripts/S11c_b_exports.py` is not overwritten.
- Save every ablation script AND its literal stdout to named absolute paths, and report those paths.

## Output
A short list of findings (each: the script quote + line, the ablation you ran + its literal diff, why it
could make the physics wrong, a one-line repair), or the explicit ablated-and-clean list. Read-only on the
working tree.
