# Independent FROM-SPEC adjudication — S11c-a T-f PROJECTION dynamic-window representation

## The question
Two blind CAS engines each computed the first-shape-order (ε) shape derivative of S11b's projection identity
(object T-f, tags `S11CA_PROJECTION_SHAPE_DERIV` / `_STATIC_OPERAND` / `_DYNAMIC_OPERAND` / `_RESIDUAL`).
They represent the smooth slab WINDOW differently. Derive the window and the projection shape-derivative
FROM THE SPEC, then determine ONE of:
  (A) the two window representations are the SAME object — give the exact identification and show by CAS that
      the two projection shape-derivatives coincide under it; OR
  (B) one engine's window deviates from the spec — name which engine and the precise dropped/added dependence,
      with the CAS residual.
Both outcomes are acceptable; a genuine divergence is the valuable output. ⛔ Do NOT force agreement and ⛔ do
NOT assume a divergence — it is entirely possible a single-argument form is equivalent to a two-argument one
(e.g. if one argument is fixed by the other through the slab thickness). TEST the identification either way.

## What you derive (from the spec ALONE, in CAS)
Object T-f (spec §4, "T-f · Dynamic-window projection"): the shape derivative of S11b's projection identity
using the anchored window `Ω^α`, plus the static-flat projection operand and their residual. Governing spec:
- §1b — the projection object = integrate the mass-conservation law `∂_tρ_4D + ∇₄·j = 0` against the window
  `Ω` and integrate by parts in `w`; use the DYNAMIC anchored window; do NOT replace it by a static window
  before differentiation.
- §3c — "Supply one fixed smooth two-argument window function `𝒪` and define `Ω^α(x,w,t) ≡ 𝒪(G_+^α, G_-^α)`",
  with `G_s^α ≡ s F_s^α` and `F_s^α ≡ w − h_s^α` (§3a). `𝒪 ≈ 1` when both arguments are negative, → 0 outside.
- §3a — the face maps/heights `h_s^α`, thickness `W^α = h_+^α − h_−^α`.
Full spec: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`.

Write your OWN sympy (preferred; no licence seat) or wolframscript derivation that builds `Ω^α = 𝒪(G_+, G_-)`
from §3c, forms the projection integral of §1b, integrates by parts in `w`, and takes the first ε-order shape
derivative for ONE branch/DOF/density case. SAVE the script AND its literal stdout to named absolute paths.
⛔ A prose derivation is discarded — script + literal stdout or nothing.

## The two engine representations to adjudicate (shown NEUTRALLY)
- Engine A: window `O_window(window_argument_plus, window_argument_minus)` — a function of TWO arguments —
  with `Subs(Derivative(O_window(·,·), (arg, n)), …)` for its partial derivatives w.r.t. EACH argument;
  integrated in the normal variable `w` via `Integral(…, (w, …))`.
- Engine B: window `windowFunction[normalCoordinate − widthProfile[x]/2 − waveOrder*(…)]` — evaluated on a
  SINGLE argument; integrated via `Inactive[Integrate][…, {normalCoordinate, …}]`.

## Determine and report (with CAS evidence)
1. From §3c: how many arguments does the window `𝒪` take, and what are they in terms of `G_±^α = ±F_s^α`?
2. Does Engine B's single-argument `windowFunction[·]` faithfully represent the spec window `𝒪(G_+, G_-)`,
   or does it drop/alter a dependence the spec requires? Show by CAS whether the two projection
   shape-derivatives are equal under a declared window-argument identification, or leave a residual.
3. If one engine deviates from the spec window, name it and the precise term, with the CAS residual.

## Method / discipline
- Derive from the spec independently BEFORE studying the engine forms.
- ⚠ The window `𝒪` is an UNSPECIFIED smooth function — keep it symbolic (an undetermined function of its
  argument(s)); do not substitute a concrete bump. The question is about its ARGUMENT STRUCTURE and the shape
  derivative, not its shape.
- Physics filter: report a finding only if it catches a way the physics or spec-fidelity is actually wrong.
- wolframscript (only if you use it): wrap EVERY kernel run in `timeout 600`; ONE kernel at a time (2-seat
  licence); copy any artifact to /tmp and work on the COPY; ⛔ never modify the working tree.
- Save every script AND its literal stdout to named absolute paths; report those paths.
