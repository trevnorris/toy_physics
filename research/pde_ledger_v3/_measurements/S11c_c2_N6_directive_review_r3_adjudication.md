# N6 build-directive review — round 3 adjudication (2026-09-06)

r2-folded directive + §5c (`ead5588a`) → 2 decision legs r3. **Legs DISAGREED** (M1: the disagreement is the
measurement): **Codex-sol = CLEAR TO BUILD; Grok = NOT clear.** Reports: `scratchpad/{codex,grok}_N6_dirreview_r3.txt`.

## The disagreement + its RESOLUTION on evidence (verified by orchestrator reading — rule 13)
- **Grok:** the μ_θ binding is a genuine **spec ambiguity**: the self-energy increment is δp-slot-only; `S_P` carries no
  `θ`, so a bulk θ-shift on the carrier is a no-op/annihilating; N4 can enter ONLY via the μ_θ binding; Eulerian-μ_θ vs
  material-μ_θ are **two different operators**, either of which passes all 6 DoD clauses; DoD 2 (optional N4) vs DoD 4
  (mandatory reverse channel) are inconsistent. ⇒ must pin μ_θ in §5c before build.
- **Codex:** CLEAR — the mandatory reverse-channel DoD forces the map; implementation freedom, not ambiguity.

**Decisive code (read):**
1. `face_generalized_force_rows` (`brane:2135-2179`): μ_θ is **bound into the virtual work BEFORE differentiation**
   (`bind_mu_theta_operand` → `diff(thickness_density, delta_v_u/e_W)`). ⇒ the face δp-rows get their θ/N4 dependence
   **through μ_θ**, ⛔ not through a θ-shift applied to the carrier rows. Grok's mechanical claim CONFIRMED.
2. **The parent face N6** `task_rep_invariance` (`S11c_a_interface_geometry_sympy_audit.py:1481-1496`) applies **NO
   field-redefinition `T`**: it builds each face/geometry quantity (`RELATIVE_FLUX, TRACTION, CLOSURE_SHAPE_DERIV,
   VIRTUAL_WORK_SHAPE_DERIV, VIRTUAL_CONSTRAINT, EVOLUTION_MASS_BALANCE`) Eulerian (`RAW_PRIMARY`) vs
   `route="MATERIAL"` and **differences them directly** (`object_difference`). The material→Eulerian map lives INSIDE
   the material face builders (covector inverse-transpose), ⛔ not as a `T` on the increment.

⇒ **Grok is correct; Codex's CLEAR-TO-BUILD is overruled.** The entire `T_{M→E}`-on-the-increment framing — carried
through r1/r2/r3 — is the **wrong type**. Codex cleared a directive built on that wrong frame.

## The correct object (parent-mirroring)
```text
route 1 (Eulerian):  I_E     = extract( close(SLAB)   − SLAB   )
route 2 (material):  I_{M→E} = extract( close(SLAB_M) − SLAB_M )    ⛔ NO separate T on the increment
```
`SLAB_M` = the c2 face-slot carrier built with **native material-route face ingredients** (traction /
`closure_shape_deriv` / virtual work / `V_s`, already Eulerian-mapped via the material builders' covector
inverse-transpose) + **material μ_θ**, folded into the **same δp symbols**; `close` = the same imported same-`(α,ρ)`
c1 response, opaque. The representation difference lives in the material-vs-Eulerian face/μ_θ construction, exactly as
the parent's `task_rep_invariance` does. The N4 advection channel enters via **material μ_θ** (the μ_θ binding is a
SPEC DECISION, ⛔ not a builder choice).

## The "author keeps breeding defects" pattern → CHANGE THE AUTHOR (rule 15; user-approved 2026-09-06)
My hand-written route-2 formula was wrong **3 rounds running** (r1 extract-then-map → r2 material_pullback-on-rows →
r3 T-on-carrier). ⇒ ⛔ stop re-inventing the formula. **astra authors the route-2 construction spec**, anchored to the
parent `task_rep_invariance`; reviewed **fresh Claude + Grok** (astra = Codex-family; G1); review-until-clear; then I
fold the cleared route-2 into §5c + directive, and the combined directive gets its Codex + Grok decision pass.

## Landed at r3 (both legs; ⛔ not re-opened): F2 advection representative (RHO4_CONSTANT absent), F3 tilt factory,
audit edit LAB−MATERIAL, tractability shape, leak fences, dimensions.

## Disposition
astra authors route-2 spec → fresh Claude + Grok (review-until-clear) → fold into §5c+directive → Codex+Grok directive
pass → build. Prior baseline `ead5588a`.
