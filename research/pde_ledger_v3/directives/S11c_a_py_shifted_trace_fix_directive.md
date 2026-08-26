# S11c-a PY fix — the traced-bulk shifted trace (§3c) — directive (rev 2, folded from 2 legs)

## Why this exists
T7 cross-engine reconciliation found the two blind engines agree exactly on the pure face-geometry
shape-derivatives and diverge on every object that traces a bulk field. The divergence was adjudicated
against the spec by independent derivation: the SymPy (PY) engine mis-applies the §3c trace-linearisation
law for the traced bulk fields. This directive re-scopes the PY engine to the spec's §3c (now clarified —
see §3c of `S11c_a_SHARED_PHYSICS.md`).

⚠ SUPPLIED / UNFALSIFIABLE WITHIN THIS BUILD: the §3c reading and the background-state premises are verified
physics, supplied as premises (§2d/§2a/§3b/§1/§3c). ⛔ Do NOT tune the output to match any recorded value or
a sibling engine, and ⛔ do NOT infer a "correct answer" to iterate toward — compute and print; the
cross-engine diff and the ablation review happen on the orchestrator's side, where a residual is a finding,
not a build failure (rule 5).

## The two PROPERTIES that must hold (state what must be TRUE; reach everything by computation)
P1 — NO FABRICATED BACKGROUND JETS. Every background face value or normal derivative that enters the §3c
   term `δh_s·∂_w f⁰` must be obtained by DIFFERENTIATING a member of the supplied background state `𝔅⁰`
   (§2d) — never introduced as a free `PREMISE` symbol. The engine currently injects free premises for the
   background normal jets of the bulk velocity, pressure, current, and density
   (`d_w_v_bulk_0_*`, `d_w_delta_p_0_*`, `d_w_j_0_*`, `d_w_rho_4D_0_*`). Replace each with the actual `∂_w`
   of the corresponding supplied background; where that background is zero or `w`-independent the jet is
   reached as zero by computation, not asserted.
P2 — PERTURBATION TRACED AT THE BACKGROUND FACE. Every traced perturbation (velocity, pressure, current,
   density) is evaluated at the background face `h_s⁰ = sW_bg/2` (§3c), so its first-shape-order dependence
   on the background face position — which carries `η` through `W_bg` (§2a) — is reached by composition,
   not frozen at the flat reference `w=sW₀/2`.
EXEMPTION: the conormal derivative's generic traced-field probe (`trace_grad_f`, `d_w_trace_grad_f`) is an
   operator on an ARBITRARY field, not a supplied physical background; its shifted-evaluation term stays. Do
   not delete it. P1 governs only the supplied physical fields (velocity, pressure, current, density).

## Scope — every object that traces a bulk field, and everything derived from it
Direct: `RELATIVE_FLUX`, `TRACTION`, `KINEMATIC_BALANCE`, `EVOLUTION_MASS_BALANCE`, `EVOLUTION_TERM_ORIGINS`,
`VIRTUAL_WORK_SHAPE_DERIV`, `CLOSURE_SHAPE_DERIV`, `FACE_SHIFT`, and the pressure/velocity/current/density
traces they consume (the `FaceSource` traces, `build_face_shift_raw`, `uniform_face_shift_reference`).
Derived (recompute automatically when the engine re-runs, but they MUST be re-emitted): the representation,
form-control, control-independence, and uniform-limit packages that consume any of the above.
⛔ Do NOT touch the pure face-geometry objects (`FACE_NORMAL`, `FACE_MEASURE_SHAPE_DERIV`, `FACE_VELOCITY`)
or the conormal probe — they are correct. Every branch, both faces, both DOFs, both density representatives,
both routes, and the controls must be recomputed for the objects in scope.

## Script obligations — non-negotiable (rule 2)
1. The script PRINTS computed objects; it never states a conclusion. Every emit payload is a CAS object.
2. PRINT operands and the residual; do NOT `assert` a value before emitting it.
3. Interpretation belongs to the step record, not the script.
⛔ The only place physical symbols may be combined by hand is in constructing the supplied §3c law and the
`𝔅⁰`/ansatz; every traced object must be REACHED BY COMPUTATION from `f = f⁰ + ε δf` composed with `w → h_s`
and from `∂_w`(supplied background), never typed as a result. A control re-enters at the field composition.

## Acceptance — VALUE-FREE and PROVENANCE-BASED (rule 5); correctness is verified DOWNSTREAM, not here
- The engine runs to completion and re-emits every in-scope object AND its derived packages, each with
  operands + residual per case (all branches/faces/DOFs/representatives/routes/controls).
- P1 holds by PROVENANCE, not by name: no in-scope traced physical object may contain a `PREMISE`-classified
  symbol standing for a background normal jet; every background jet in it is `∂_w` of a supplied `𝔅⁰` member.
  (The orchestrator verifies this over the deliverable with the reduction provenance tooling and by form
  ablation — ⛔ do NOT self-certify by grepping for the old names; renaming a premise does not satisfy P1.)
- P2 holds by construction: the traced perturbation objects depend on the background face position
  (a `W_bg`/`w₁` factor is reached by composition, not typed as a monomial).
- ⛔ Agreement with any recorded value or sibling engine is NOT a build acceptance item. The orchestrator
  diffs cross-engine over all cases and the review legs ablate the deliverable; do not target any value.

## Target
`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py` (HEAD 2d0f0055). The clarified
§3c is already in `S11c_a_SHARED_PHYSICS.md`; derive from it.
