# Independent FROM-SPEC adjudication — S11c-a T-a′ CONORMAL_DERIV cross-engine representation

## The question
Two blind CAS engines each computed the first-shape-order (ε) shape derivative of the S11c-a conormal
operator (object T-a′, tag `S11CA_CONORMAL_DERIV`). They emit it in DIFFERENT representations. Derive the
object yourself FROM THE SPEC, then determine ONE of:
  (A) the two representations are the SAME physics — give the exact spec-grounded relation identifying them,
      and show by CAS that the two forms are equal under it; OR
  (B) one engine deviates from the spec-derived object — name which engine and the precise term/order, with
      the CAS residual that exhibits it.
Both outcomes are acceptable; a genuine divergence is the most valuable output. ⛔ Do NOT force agreement,
and ⛔ do NOT assume a divergence — TEST the identification either way.

## What you derive (from the spec ALONE, in CAS)
Object T-a′ (spec §4, "T-a′ · Conormal"): the face operator `n̂_s·∇₄`, INCLUDING evaluation on the supplied
graph. Governing spec sections you must read and use:
- §3a — the face map `R_s^α`, the Eulerian graph heights `h_s^α`, `F_s^α ≡ w − h_s^α`, the outward unit
  normal `n̂_s^α` (unit normal to `F_s^α=0`, orientation `s(n̂·ŵ)>0`), the graph measure `a_s^α=√(1+|∇_x h_s|²)`.
- §2a — the supplied non-uniform background `W_bg` and its `η`/`σ_W`/`w1` profile structure.
- §3c — the shifted-trace / evaluation-on-graph law `δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s ∂_w f⁰(x,h_s⁰)`; evaluate
  the perturbation trace at the background face `h_s⁰ = sW_bg/2`.
Full spec: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`.

Write your OWN sympy (preferred; no licence seat) or wolframscript derivation that builds `n̂_s·∇₄` from §3a,
applies it to a generic bulk field, evaluates on the graph, and takes the first ε-order shape derivative for
ONE branch/face/DOF case (e.g. LAB_HELD, minus face, δW). SAVE the script AND its literal stdout to named
absolute paths and report them. ⛔ A prose derivation is discarded — script + literal stdout or nothing.

## The two engine representations to adjudicate (shown NEUTRALLY — derive FIRST, compare LAST)
Case LAB_HELD | FACE_MINUS | DOF δW, the first ε-order shape derivative:
- Engine A expresses it in `trace_grad_f_{1..4}` and `trace_grad_f_{1..4}_dw` (components of the traced
  gradient of a generic bulk/test field and their normal jets), with an overall `W_0` power. Sample terms:
  `W_0*e_W*trace_grad_f_4_dw/2  −  W_0*e_W_d1*trace_grad_f_1/2  +  W_0*e_W*sigma_W*trace_grad_f_1_dw*w1_profile_d1/4  + …`
- Engine B expresses it in `conormalBackground` / `conormalPerturbation` (the background/perturbation parts of
  the conormal field) and their normal jets (up to THIRD normal order), with an overall `W_0^2` power and
  `eta_bg`. Sample terms:
  `W_0^2*∂_{x1}∂_w(conormalBackground)*e_W_d1*eta_bg*w1_profile/4  −  W_0^2*∂_{x1}∂_w∂_w(conormalBackground)*e_W*eta_bg*sigma_W*w1_profile*w1_profile_d1/8  + …`
Here `e_W ≡ δW/W_0` and `e_W_di ≡ ∂_{x_i} e_W`; `w1_profile` is the §2a background profile; `sigma_W`,`eta_bg`
are the §2a contrast/scale. (These symbol names are given only to let you read the two forms — do not treat
them as an answer.)

## Determine and report (with CAS evidence)
1. From your independent derivation: what IS `δ(T-a′)` to first ε-order, in the natural variables?
2. Is there a spec-grounded relation identifying Engine A's `trace_grad_f` representation with Engine B's
   `conormalBackground/Perturbation` representation — e.g. the conormal defined via a generic test field, or
   the evaluation-on-graph (§3c) bringing the higher normal orders and the extra `W_0`? If yes, GIVE it and
   show by CAS that the two forms coincide.
3. If no such relation reconciles them, state which engine's form disagrees with the spec-derived object and
   the precise term/order, with the CAS residual.

## Method / discipline
- Derive from the spec independently BEFORE studying the engine forms.
- Physics filter: report a finding only if it catches a way the physics or spec-fidelity is actually wrong.
- wolframscript (only if you use it): wrap EVERY kernel run in `timeout 600`; ONE kernel at a time (2-seat
  licence); copy any artifact to /tmp and work on the COPY; ⛔ never modify the working tree.
- Save every script AND its literal stdout to named absolute paths; report those paths.
