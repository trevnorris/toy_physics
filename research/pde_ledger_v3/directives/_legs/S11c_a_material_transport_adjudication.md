# From-spec adjudication — does the T-e MATERIAL_ADVECTED density trace carry a material-transport term?

Two independent CAS engines (SymPy + Wolfram) build S11c-a blind. A build review raised a concern about
the **MATERIAL_ADVECTED × RHOBR_CONSTANT** shifted-face density trace (T-e, `S11CA_FACE_SHIFT`). Both
engines currently AGREE on the answer below, so a cross-engine residual cannot settle whether it is right —
you must decide it **from the spec**. Do not assume the engines' agreement is correct; if the spec requires
a term both omit, that is a shared-blind-spot finding.

## The spec facts (source of truth: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`)
- **§2c**: `MATERIAL_ADVECTED: Q_bg^M(x,t) ≡ Q_bg(χ(x,t))` for every background profile, including
  `ρ_4D,bg⁰`. "distinct physical anchorings … an engine may not substitute another parameterization."
- **§3a**: `ρ_4D,bg^{0,M} ≡ ρ_4D,bg⁰(χ(x,t))` (material coords χ), vs `ρ_4D,bg^{0,L} ≡ ρ_4D,bg⁰(x)` (lab).
- **§3c** trace law: `δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s ∂_w f⁰(x,h_s⁰)`; "the supplied density background
  depends on the in-plane anchor, not on `w`"; "Evaluate the perturbation trace at the background face
  `h_s⁰ = sW_bg/2` … retain its first-shape-order dependence." "Which trace terms then survive is computed
  from these premises, not stated here."

## What both engines currently emit (verbatim, MATERIAL_ADVECTED × RHOBR_CONSTANT, FACE_MINUS, DOF=DELTA_W)
```
PY : -W_0*delta_rho_4D_face_minus_dw*eta_bg*w1_profile/2 + delta_rho_4D_face_minus
WL : -W_0*eta_bg*(d/dw rhoBulkPerturbation)*w1_profile/2 + rhoBulkPerturbation(x1,x2,x3,(-W_0/2,time))
```
i.e. the perturbation trace re-centred to the background face, **with no term carrying the background
density's in-plane variation** — identical in form to the LAB_HELD case.

A "branched" alternative (evaluating the background at χ and differentiating through the material map)
would instead add a material-transport term of the shape `~ η·rho_br·(u·∇w1)/(W_0(1+η·w1)^2)` (u = the
material displacement DOF, w1 = the thickness profile).

Separately: the engines' T-e `EXACT_TRACE_SOURCE` (the un-differentiated trace) shows the MATERIAL density
background at **lab** coordinates `ρ(x)`, not `ρ(χ)`.

## The questions (derive each from the spec; show computation)
1. **Does the §3c T-e shifted density trace for MATERIAL_ADVECTED carry a material-transport term** (the
   in-plane variation of the background density under the shape/advection), or does §3c's law — a *normal*
   shift `δh_s·∂_w f⁰` with `∂_wρ⁰=0` — correctly exclude it? Derive the first-shape-order shifted density
   trace for the MATERIAL branch from §2c/§3a/§3c yourself and state whether the transport term appears.
2. **Where does the material advection legitimately enter?** Check whether the advective displacement
   `u·∇W_bg` enters the face displacement `δh_s` (the MATERIAL-branch face shift, §3a/§3b) and/or the
   projection/evolution — and whether putting it *also* into the T-e density trace would double-count.
   Both engines route the advection through `δh_s` (MATERIAL branch) and the projection integrand; is that
   the correct home for it?
3. **Is the engines' omission correct, or a shared blind spot?** State plainly: are both engines right to
   emit the density trace with no material-transport term, or does the spec require that term (both wrong)?
4. **The EXACT_TRACE_SOURCE frame.** Should the un-differentiated MATERIAL density trace source display
   `ρ(χ)` per §3a rather than `ρ(x)`? Does that frame choice change the emitted *shape-derivative* operand
   (question 1), or is it a display-only fidelity point with no effect on the computed T-e object?

## Method
- Derive from the spec; do not decide by "the engines agree." A prose derivation is discarded — write a
  runnable SymPy (or Wolfram) snippet, run it, paste its literal stdout, and give absolute paths (save
  under /tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/).
- Distinguish the shape-derivative OPERAND (physics-bearing) from the EXACT display (fidelity).
- Physics filter: a finding is a way the emitted T-e density operand is physically wrong per the spec.

## Return
For each of the 4 questions: a spec-grounded verdict with file:line citations and your snippet's literal
stdout. End with: is the current T-e density operand correct, and if not, exactly what term is missing and
which spec clause requires it.
