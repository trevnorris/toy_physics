# Independent physics review — S11c-b Wolfram engine (a SCRIPT, Mathematica)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
(repo root `/var/projects/toy_physics`).

A Codex-built **blind** Wolfram engine: it imports nothing and re-derives every object from the specs. It
constructs the S11c-b variable-coefficient brane operator and the off-diagonal transverse→{θ,e_W,u_L}
coupling kernel, and emits `WL_S11CB_*` tags.

## What you are handed
- The script above (the artifact).
- The physics authority to derive against and check faithfulness to:
  `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (sole authority),
  and the specs it inherits (`S11c_a_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md`).
- ⛔ You are NOT handed the sibling SymPy engine's output, and must not seek it — the cross-engine check is a
  separate downstream comparator; your job is whether THIS engine faithfully computes the spec. (This is the
  blind engine; its whole value is independent construction.)
- There is no expected numeric answer: the spec withholds every computed value (the coupling grade/sign, the
  new invariants, the admissibility residual, the basis count). ⛔ Do not treat any value in the script as a
  target to confirm.

## Required method — this is a SCRIPT; derive independently and ABLATE
Write your OWN small derivation (a short `.wl` or `.py`) for the load-bearing physics BEFORE opening the
artifact in depth, and save BOTH your script and its literal stdout to named absolute paths (report them). A
prose re-derivation is discarded.

⛔⛔ **A FORM ABLATION IS MANDATORY.** Change the STRUCTURE of a load-bearing object in a `/tmp` COPY of the
`.wl` and re-run the affected task, reporting the LITERAL diff. A COEFFICIENT rescale tests arithmetic; only
a FORM change tests physics. Do at least these four (each on a /tmp copy, re-running only the affected task):

1. **The coupling is gradient-driven.** Zero a background FIRST JET at its source (the supplied `∂W_bg` /
   `σ_W·∂ξ w₁` map, or `∂μ_R,bg`) and re-derive the coupling kernel. The off-diagonal transverse→{θ,e_W,u_L}
   block must COLLAPSE when the background gradient is removed (its uniform limit is S11b's decoupled zero,
   spec §1d/§5c). A gradient-INDEPENDENT survivor at jet→0 is a defect — but ⚠ first check whether such a
   survivor is a **total divergence** (representational, routed to §5c adjudication by §1d) before calling it
   a physics defect.
2. **The energy is CONSTRUCTED, not substituted.** Is the variable-coefficient basis (§3a) built from the
   symmetry group with the background jets as spurions, EMITTING new gradient-of-background invariants — or is
   it just `W_0→W_bg` substituted into the uniform `U` (which omits the undifferentiated-`u` N15 couplings,
   spec §1a/§1c/§3a)? Ablate: remove one spurion channel and confirm a new-invariant / kernel term disappears.
3. **⭐ Admissibility is the background-order (ε⁰) balance, NOT ε→0 of the wave operator.** This is the highest-
   priority check. Confirm `WL_S11CB_ADMISSIBILITY_OPERATOR_OPERAND` is the background-order generalized force
   (body force + per-face traction) the PROFILE sources — the first variation of the FULL-field
   energy-and-geometry at 𝔅⁰ with the profile's own gradients retained (e.g. the `κ_W` bending term at
   `W=W_bg`) — in the same pairing as `𝒮_hold⁰` (spec §2b/§3d). ⛔ If it is instead the wave operator
   evaluated at zero perturbations (`ε→0` of the §3b bilinear operator), it is identically zero and the
   can-fail N12 test is vacuous. Ablate: check whether the operand is structurally the `ε→0` limit of the
   §3b operator (→ it will be identically zero — report that with the literal output).
4. **The off-diagonal block is EXTRACTED from the operator, not asserted or built by a parallel route.**
   Verify the kernel is read off the divergence-form operator's curl↔div structure (§3c). Ablate: break the
   operator construction and confirm the kernel emission breaks (a parallel-route or hand-typed kernel
   survives).

Also report the general script defects: any premature guard that hides an ablation (a check that crashes
rather than emits); any emitted payload that is a HAND-TYPED object with no data dependency on the derivation
(delete the construction and the output does not move); any TAUTOLOGICAL residual (`A−B` zero by construction
for any input — e.g. a control whose "baseline" and "corrupted" operands are identical, or a self-adjointness
residual that is `A−A` by mixed-partial commutation); any conclusion emitted as an unconditional literal; any
`VERDICT`/`PASS`/`FAIL` token or a native/evaluated boolean serialized as a residual operand (§6 forbids
both — a boolean must be the unevaluated relational retaining its operands). For each emitted §4 object, ask
WHICH LINE COMPUTED THIS.

## Physics filter
Report a finding only if it catches a way the physics could be WRONG (a coupling not gradient-built, a
substituted energy missing N15, a vacuous admissibility, a hand-typed/tautological emission, a guard hiding
an ablation). Do not report style. A leg that finds nothing is weak evidence — say which load-bearing objects
you ablated and the literal diffs.

## Operational constraints (IN THIS PROMPT, both legs identical — Mathematica)
- ⛔ Wrap EVERY kernel run in `timeout 600`. A 600s hit is a FAILED ablation — report it and move on. ⛔ Never
  raise the timeout.
- ⛔ Run at most ONE Wolfram kernel at a time (the licence has TWO seats and the sibling leg may hold one).
- ⛔ If a kernel's resident memory passes 6 GB, kill it by PID and report — an orphaned kernel OOMs the machine.
- ⛔ Copy the `.wl` to `/tmp` and ablate the COPY. ⛔ Never modify the working tree; ⛔ never write under
  `mathematica/out/`.
- Save every ablation script AND its literal stdout to named absolute paths, and report those paths.

## Output
A short list of findings (each: the script quote + line, the ablation you ran + its literal diff, why it
could make the physics wrong, a one-line repair), or the explicit ablated-and-clean list. Read-only on the
working tree.
