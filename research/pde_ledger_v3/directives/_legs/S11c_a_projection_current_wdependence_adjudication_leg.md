# Independent spec adjudication — S11c-a T-f: is the projection's perturbation current w-dependent, or face-frozen?

## Your task
S11c-a builds the interface shape-derivatives of a toy superfluid brane in two INDEPENDENT computer-algebra
engines (a SymPy engine and a blind Wolfram engine), each re-deriving the same physics from one shared spec.
On the T-f object (`S11CA_PROJECTION_*`, the dynamic-window projection identity and its shape derivative),
the two engines represent the **perturbation bulk current `δj`** differently inside the projection. Your job:
determine, **from the spec alone**, which representation is faithful to the T-f projection object — or
whether both are, or neither, or whether the spec is genuinely ambiguous on this point. Derive it yourself;
do not assume either engine is correct, and do not try to reverse-engineer an "expected" answer from how this
is phrased.

## A measured fact you may take as given (so you do not re-litigate it)
A separate CAS computation (independent of this adjudication) established that these two representations are
**NOT algebraically equivalent** — they yield a genuinely different projection object; the difference is a
window-shape-weighted term carrying the current's normal-derivative profile `∂_w δj_w`, nonzero for a generic
`w`-dependent current. So this is a real fork, not two spellings of one object. Your task is **only** which
representation the spec's T-f object requires. (Do NOT treat "they differ" as evidence that either specific
engine is the wrong one — that is exactly what you must decide from the spec.)

## Read the source of truth FIRST, before the engine forms
Open and read these, form your own view of what the T-f projection object requires, and only then read the
engine treatments below. Quote the **exact** spec sentences you rely on — do not paraphrase a premise you
are going to lean on.
- `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`
  - **§1b** — the current law `j = ρ_4D v_bulk`, the conservation law `∂_tρ_4D + ∇₄·j = 0`, and the sentence:
    "The projection object is the result of integrating this conservation law against `Ω` and **integrating
    by parts in `w`**." Note `∇₄·j = ∂_i j_i + ∂_w j_w` includes the **normal**-derivative term `∂_w j_w`.
  - **§3a** — the slab geometry: `w` is the normal coordinate; the slab spans a finite thickness `W_bg`
    about the face `h_s`.
  - **§3c** — the shifted-trace law `δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s ∂_w f⁰(x,h_s⁰)` (which applies to
    **traced** bulk fields evaluated **at the face**), the two-argument dynamic window, and the sentence
    beginning "In this scope the traced bulk velocity, the perturbation pressure, and the bulk current have
    zero background …".
  - **§4 T-f** — the `S11CA_PROJECTION_*` objects: the projection identity, its dynamic-window and
    static-flat routes, and its shape derivative.

## The exact engine forms (re-read the RAW files yourself — a mis-copied form has manufactured a false finding here before)

**Engine P — SymPy** (`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`):
- lines ~175-181: `j_bulk = (delta_j_bulk_1 … delta_j_bulk_4)` are **plain symbols** (no `w`-argument);
  `grad_j_bulk[i][j] = delta_j_bulk_i_dj` are **in-plane** jets; `dw_delta_j_bulk[...]` are **normal** jets.
- the projection is built in `projection_terms` (lines ~1114-1171). Read it. It uses `j_bulk[i]` (constant),
  `current_divergence = Σ_i grad_j_bulk[i][i]` (in-plane jets **only**), and
  `WINDOW_NORMAL_CURRENT = -∫ j_bulk[3]·window_normal dw`. Check for yourself whether the **normal** jets
  `dw_delta_j_bulk` (the `w`-variation of `δj`) enter the projection at all.
  ⇒ Engine P represents each current component as its value at the background face, **independent of `w`**.

**Engine W — Wolfram** (`research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl`):
- lines ~442-448:
  ```wolfram
  currentWWave[coordinates_List, normal_] := currentWPerturbation @@ Append[coordinates, {normal, time}];
  currentXWave[index_][coordinates_List, normal_] :=
    Symbol["currentXPerturbation" <> ToString[index]] @@ Append[coordinates, {normal, time}];
  ```
  so each perturbation current component is a **bulk field `δj_a(x, w, t)` that depends on the normal
  coordinate** `normal` (= `w`).
- line ~183: "Projection integrals stay opaque to every outer algebraic transform" — the projection integral
  `∫ … dw` is kept formally un-evaluated; the `w`-dependent current sits **inside** it.
  ⇒ Engine W represents each current component as a full `w`-dependent bulk field integrated across the slab.

## The crux question to settle from the spec
The conservation law integrated in §1b is `∂_tδρ + ∇₄·δj`, whose normal term `∂_w δj_w` is integrated by
parts in `w` against the window. Does the T-f projection object therefore **retain the perturbation current's
`w`-dependence across the slab** (Engine W), so that `∂_w δj_w` contributes a real `-∫(∂_wΩ)δj_w dw` term — or
is the perturbation current a **face-traced quantity** (§3c), evaluated at `h_s⁰` and hence `w`-independent
(Engine P), so its normal variation does not enter? Is `δj` in the projection a *bulk* field being integrated
over `w`, or a *traced* field at the face?

## Output
- A verdict: **Engine P faithful** / **Engine W faithful** / **both** / **neither** / **spec ambiguous**.
- The exact spec sentences (quoted) that force your verdict.
- If the spec does not unambiguously determine it, say so plainly — a one-engine correction is a spec
  question first, and a closed spec is not necessarily a correct one. Distinguish "the spec clearly requires
  X" from "the spec is silent / underspecified and both readings are defensible."
- This is a spec-reading adjudication: read the spec first, form your own view, then read the engine forms.
  Do not assume; do not guess an expected answer.
