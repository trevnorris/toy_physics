# Independent physics review — S11c-b SymPy engine, REPAIR round 1 (a SCRIPT)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
(repo root `/var/projects/toy_physics`). A Codex-repaired SymPy engine. Its first build had four defects
(B1–B4) found by FORM ablation; those were just repaired against the corrected spec. **Verify the repairs
actually hold, and that nothing ablation-clean regressed.** There is no expected numeric answer (the spec
withholds every value); do not treat any value in the script as a target.

## What you are handed
- The script above.
- Physics authority: `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`
  (the corrected §3a/§3c/§3d/§5a), and the specs it inherits (`S11c_a_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md`).
- The prior finding record (what B1–B4 were):
  `/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_b_sympy_build_review.md`.
- ⛔ NOT the sibling Wolfram engine's output — the cross-engine check is a separate downstream comparator.

## Required method — SCRIPT; derive independently and ABLATE (a prose re-derivation is discarded)
Write your OWN derivation for the load-bearing physics before studying the artifact; save it + literal stdout
to named /tmp paths and report them. Copy the script to /tmp and ablate the COPY; re-run only the affected
task; report the LITERAL diff. ⛔ Never modify the working tree; ⛔ do not run the full package repeatedly (it
writes a 28MB export — direct any run to a /tmp copy / scratch path).

**Verify each repair by FORM ablation (the key checks):**
1. **B1 (BLOCKING) — admissibility is now the background-order (ε⁰) balance, NOT vacuous.** Confirm
   `PY_S11CB_ADMISSIBILITY_OPERATOR_OPERAND` is the first variation of the FULL-field (§3a) energy at 𝔅⁰ with
   the profile's own gradients (`∇(W_bg+δW)`) retained — a genuine functional of the background (has data
   dependence on `W_bg`/its jets, including the generated second spatial derivative). ⛔ It must NOT be
   identically zero and must NOT be `zero_wave(euler_derivative(wave_density,…))`. Ablate: inject/alter a
   background profile jet at the source and confirm the operand MOVES (nonzero data dependence); confirm it is
   not the ε→0 limit of the §3b wave operator.
2. **B2 — coupling kernel is extracted from the §3b operator via the §3c weak variational restriction**
   (solenoidal/irrotational trial+test), and its off-diagonal block reduces to §1d's decoupled zero in the
   uniform limit. Ablate: delete the operator construction → the kernel must break; zero the background jets →
   the off-diagonal block must collapse (no gradient-independent survivor left in the emitted block).
3. **B3 — controls are non-tautological.** The adjointness residual is the pairing-based operator-block
   residual, ⛔ not `∂²U/∂u_T∂e_W − …` (Clairaut). The independence control's corruption reaches the kernel and
   admissibility (⛔ not `A−A` identical calls). Ablate/inspect: confirm a genuine one-sided corruption moves
   the corrupted operand for kernel and admissibility (or the object is omitted only on computed structural
   absence).
4. **B4 — no dormant dimension crash** from the `uniform_coefficient` fallback.
5. **No regression:** the energy-basis spurion construction (N15), the §3b operator, and the gradient-driven
   coupling must still be ablation-clean (the repair was supposed to keep them unchanged).

Also report any general defect: an `assert` before the value it guards; a hand-typed payload with no data
dependence; a tautological residual; a VERDICT/PASS/FAIL or native-boolean residual; and for each §4 object,
WHICH LINE COMPUTED IT.

## Physics filter
Report a finding only where the physics could be WRONG (a still-vacuous admissibility, a kernel not reducing
to the uniform zero, a still-tautological control, a regression). A leg that finds nothing is weak evidence —
state which objects you ablated and the literal diffs.

## Output
Findings (script quote + line, ablation + literal diff, why wrong, one-line repair), or the explicit
ablated-and-clean list with the repair-verification diffs. Read-only on the working tree.
