# Directive — REPAIR the S11c-a spec (targeted; ⛔ NOT a re-author)

**You (Codex) authored this spec; this is a targeted repair of it.** Two review legs found **§1 fidelity
CLEAN** (every inlined S11b equation verified term-by-term) and the structure sound. ⛔ **KEEP all of that.**
Fix exactly the five items below and leave everything else byte-stable. Reviewed by a fresh Agent + Grok.

## Deliverable
Edit `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` in place (currently HEAD `7b4a50de`).
⛔ Do not re-author, restructure, or re-derive the sound parts.

## ⛔ KEEP — both legs confirmed correct, ⛔ do not touch
§1 fidelity (all inlined S11b equations); the obligation-to-compute discipline / no rule-3 typed results
(⛔ do not introduce "n_∥=O(η)"-style asserted terms/orders); the two anchoring maps and `V_{n,s}`/`δ_vx_s`
construction; the `(η, L_W, σ_W)` slope/contrast bookkeeping; T-a…T-i coverage with `δW` and `ζ_c` on every
object; the reserved-to-S11c-b admissibility residual; scope deferrals; exit semantics; F9 reservations; the
form/independence control structure (except the T-h addition in FIX 3).

## FIX 1 (the core — this corrects a FLAWED INSTRUCTION in the authoring directive)
⚠ The authoring directive's item 6 told you to "retain nonzero background face fields `p_s⁰, J_s⁰`." **That
instruction was wrong.** It forced §3b's affinity/traction into **total** pressure
(`𝒜_s^α = μ_s^α − p_s^α/ρ_m`, `t_s^α = −(p_s^α + Λ_X𝒜)n̂`), which (a) makes the background flux
`J_s⁰ = −Λ_A(0)p_s⁰/ρ_m ≠ 0` on a flat rest wall — and the kinematic identity then makes `J_s⁰` a
**bulk-normal drain**, which S11b forbids as an active term (`S11b_SHARED_PHYSICS.md:111`; decision N11a:
the drain is a *standing scope limit*, the coupling is *slab-side geometry*); and (b) breaks the
uniform-limit regression (S11b's traction has no `p⁰`, so §5c's `(η,σ_W)→0` residual is `−p_s⁰ n̂`, not a
geometry regression); and (c) runs a static pressure through the memory kernel `Λ_A(ω)`, a different
constitutive law from the inlined S11b one.

⭐ **The fix — restore S11b's PERTURBATION-pressure law in §3b, and make the leading-order background
impermeable and at rest (a supplied premise, not a computed result):**
```
leading-order background:  V_s⁰ = 0 ,  J_s⁰ = 0 ,  𝒜_s⁰ = 0
   (μ_s⁰ = 0 since θ⁰ = e_W⁰ = 0 and U is quadratic; δp_s⁰ = 0 by definition of the perturbation)
𝒜_s^α ≡ μ_s^α − δp_s^α/ρ_m ,
t_s^α ≡ −(δp_s^α + Λ_X(ω) 𝒜_s^α) n̂_s^α  +  t_hold,s⁰ ,
```
where `δp_s^α` is the **perturbation** bulk face pressure (background value zero). ⛔ Any static background
face pressure `p_s⁰` the externally-held curved brane carries enters **only** as the mechanical background
traction `t_hold,s⁰` in the support bundle `𝒮_hold⁰`; ⛔ never inside `𝒜`, ⛔ never through `Λ_A(ω)`.
⭐ In §2d, ⛔ do not declare `p_s⁰, J_s⁰` "not set to zero"; the leading-order face flux is zero by the
closure, and a static `p_s⁰` (if present) sits in `𝒮_hold⁰`. The true-area measure `a_s` stays in the exact
laws; because `J_s⁰=0`, its shape-derivative multiplies the **perturbation** flux/traction (drain-free) — the
`δa_s·J_s⁰` term simply vanishes, so item 6's valid surface-measure point survives while the drain
reactivation is removed.

## FIX 2 (Grok F1 — a binding identity living in a word)
`Σ^α` (the LHS of the tilted mass balance, §3b) is used but never **defined** by an equation. Supply, from
the §3a face maps, the graph thickness and slab mass:
```
W^α(x,t) ≡ h_+^α(x,t) − h_−^α(x,t) ,      Σ^α ≡ ρ_4D^α W^α ,      ρ_4D^α = ρ_4D,bg⁰(1+θ) .
```
so T-g/T-h/T-i cannot disagree on whether `u·∇W_bg` and `u·∇ρ_4D,bg⁰` appear (the N4 channel).

## FIX 3 (Grok F3 — independence-test coverage)
§5a runs the two-route + one-sided-corruption package only for T-g, T-c, T-d, T-i. **Add T-h** — it is the
one object where N3 (area/tilt) and N4 (advection, `u·∇Σ_E⁰`) sit in the same equation and can cancel or
double-count. Include the advection-omission corruption (`Σ_E(x(X,t),t) → Σ_E(X,t)`) on T-h's direct route.

## FIX 4 (Grok minor — the freeze as a restriction)
The freeze is shown as two density maps but never stated as a restriction of the triple. Add: `W_bg` varies
in every branch; the two representatives are the two factorizations of `ρ_br = ρ_4D W`; an independent
`ρ_4D(x)` at frozen `W` (pure N4, no N3) is **not** a background here.

## FIX 5 (Agent F2 minor — terminology)
"externally held" (the §2d admissibility premise, applies to **both** anchoring branches at background order)
collides with "LAB_HELD" (a §2c anchoring). Rename the premise (e.g. "support-stabilised") or state
explicitly that it applies identically to both §2c anchorings.

## Report back (⛔ under 20 lines)
Which fixes you applied and where; anything in the fix list that looked wrong (you are the author — flag it,
⛔ do not silently deviate). ⛔ No physics conclusions.
