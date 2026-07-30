# Stage 044-v2 — redo 044 with a DYNAMICAL-Σ sleeve (un-freeze S_hold): PREP / DECISION ANCHOR

> **Status: PAUSED ledger thread; the decision was MADE this session (2026-07-24) but not yet built.** Paused
> behind the dimension-rewrite/integration-test detour (see
> `manifests/DIMENSION_REWRITE.md`). Resume after the integration
> work reaches a milestone. This note is the durable record of the decision + the mechanically-confirmed locality.

## The decision (user, 2026-07-24)
Stage 044 (VII-2a) committed `S_hold = −∫dt∫_{Γ_Σ}d³A λ_Σ(r_B−½)` — a holonomic pin of the `r_B=½` contour to a
prescribed body-attached mid-surface. That FREEZES the sleeve SHAPE (fixed-Σ MVP). The user ruled it out — "nothing
is static" applies to the shape too — and directed:
1. **Un-freeze the shape** at the ACTION level (promote Σ to a dynamical free boundary).
2. **Commit the bending/anchoring knobs** (Tier-2, user's explicit choice over the minimal no-new-knob route):
   `κ_bend·∫H²` (Willmore/Helfrich; dim `E/L`), `κ_anchor` (body-anchoring; `E/L³`), collar-tension (`E/L²`) — 3 new
   `[POSTULATE]` knobs (currently zeroed slots in the G0 card §2.2; well-posed/regularizes cap pinch-off).
3. **REDO stage 044** (user: "if we're changing the premise… revisit it") — NOT a downstream patch. The reopen is
   **044-local**; 044's other results survive.

## The un-freeze mechanism (author the 044-v2 directive around this)
Replace `S_hold` (frozen level-set pin) with a dynamical sleeve-surface action `S_Σ`:
- leading tension = the **DERIVED** `σ_wall = √(2a_Bκ_B)/6` (register R20 — NO new knob) + the 3 committed bending
  knobs `κ_bend∫H² + κ_anchor + collar-tension`.
- the `λ_Σ` multiplier is REMOVED, replaced by the free-boundary Euler–Lagrange law (Young–Laplace + bending
  balance); the collar-Kirchhoff conormal jump `[κ_χ∂_n r_B]=λ_Σ` becomes the tension/curvature balance
  `[κ_χ∂_n r_B] = σ_wall·H − κ_bend(…) + …`.
- The wall SHAPE is now dynamical; the full free-boundary SOLVE stays **SIM-DEFERRED** (structure-complete-at-floor,
  like the rest of the ledger). ⚠ Honest: a pure-tension throat wants to pinch off (Rayleigh–Plateau); bending
  regularizes it and what holds the throat open is the drain-flow/amplitude balance — a sim question, nothing static.
- ⛔ **Constraint:** `S_Σ` lives entirely in the `r_B`/Σ sector with **ZERO coupling to the charge field `h`** — this
  preserves stage032's `internal_inconsistency=none` (032 keeps `S_hold` off `h`).
- Template for the bending functional: the pathA_26 Willmore side energy `λ_W·(16π/9)·L` (`H=2/(3a)`,
  `software/stage1_solver/tools/pathA_26_derrick.wl:314`) — a real `∫H²dA` instance to adapt.

## Locality — CLEAN (confirmed by fresh reads AND mechanically by the stage044 manifest S_hold-map)
- **S_hold-DEPENDENT (what the un-freeze edits):** the `s_hold` claim; `action_roster` (6→5 summands, S_hold→S_Σ);
  `field_set_union` (drops the `λ_Σ` aux, RETURNS the pinned wall mode as a physical DOF); the `dimensional_homogeneity`
  term count (24→23 + the new S_Σ terms); NUANCE — the wall-Hessian translational zero mode `r₀'` becomes LIVE
  (un-freezing correctly revives it).
- **S_hold-INDEPENDENT (survives untouched):** the wall dedup `r_B≡χ_B` (`λ_χ=4a_B⇒ℓ=δ`; proven on the UNCONSTRAINED
  wall, so robust by construction), `Z_χ`, the three wave speeds (`c_±²=(3±√2)/2`, `c_γ²`, `c_s0²`), `drain_named_both`
  (S_hold does NOT toggle the drains), P-retirement.
- **Stage 032 independent:** there `S_hold` is only a string scope label (`"r_B-1/2 only"`), no load-bearing consume —
  032 needs at most a one-line "superseded-MVP" pointer; its charge result (four `A_X`, `R1_REQUIRED(bc_selection)`)
  is untouched. (032 is NOT the origin of the frozen functional — 044 is; 032 only referenced the scope abstractly.)

## Consequence + process
- The 3 new knobs bump the parameter count → deferred committed re-count to **046** (alongside `Z_χ`); add register
  rows for `κ_bend`/`κ_anchor`/collar-tension.
- ⚠ **Process updated 2026-07-29/30** ([[feedback-physics-not-ceremony]]): ⛔ ~~full new-derivation gauntlet + GLM-5.2
  tertiary + Codex→Grok→Codex bookend~~ ⇒ **ONE design-review pass on a fresh reviewer, TWO mutually
  independent fresh review legs (each fidelity + per-tooth ablation — a stage's math is physics-bearing,
  amended 2026-07-30), physics leg BLOCKING.** Author 044-v2 directive → one design-review pass
  → dual-engine build re-deriving the verdict with `S_Σ` → the existing 044 teeth (now certifying, as COMPUTED
  results, that the dedup/`Z_χ`/wave-speeds survive) + NEW teeth (`S_Σ` dim-homogeneity, the free-boundary EL law, the
  no-`h`-coupling firewall, bending-term properties) → deliverables → commit. 044's note/card: fixed-Σ MVP → dynamical
  Σ, the frozen form documented as a SUPERSEDED departure (surviving-solution-only). Then **045** builds the
  non-variational block on the corrected, fully-dynamical parent action.
