# Independent physics review — the S11c-a spec (Codex, after a targeted repair)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`

⚠ This Codex-authored spec was just **repaired** (a targeted edit, not a re-author). A prior review found
its **§1 inherited setup fidelity CLEAN** (two independent term-by-term audits against S11b), and **§1 is
byte-identical after the repair** — so ⛔ do not re-audit §1 line-by-line; a quick spot-check suffices. The
repair touched **§2b, §2d, §3a, §3b, and §5a**. ⛔ Do not assume the repair is correct: a repair can be
physically wrong or introduce a new defect. Derive your own view and let it disagree.

## What it is
S11c is the non-uniform transverse coupling at a brane–bulk interface — *is light's confinement
unconditional?* The prior step **S11b** (uniform) is CLOSED (transverse coupling identically zero on a flat
wall). **S11c-a** un-freezes the in-plane background and computes the first-order shape-derivative of every
S11b interface law. Two blind engines (SymPy + Wolfram) read this spec; an error corrupts both.

## Derive your own view FIRST of the CONSTITUTIVE + BACKGROUND setup, then read the changed sections
The repair's crux is the **constitutive form and the leading-order background**. Before reading the spec,
settle for yourself, from the sources, what is correct:
- **Which pressure enters the affinity/traction?** S11b's law is `𝒜_± = μ_s − δp_±/ρ_m`,
  `t_± = −(δp_± + Λ_X𝒜_±)n̂_±` with `δp` the **perturbation** bulk face pressure
  (`research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md:223,370`). A *total*-pressure form is a
  different constitutive claim.
- **The drain.** `v_bulk_normal_0` is NOT an active operator term (`S11b_SHARED_PHYSICS.md:111`); decision
  N11a (`research/pde_ledger_v3/directives/S11c_decisions.md:N11`) makes it a standing rest-frame scope
  limit, the coupling being slab-side geometry. At background order with faces at rest, the kinematic
  identity `n̂·v_bulk⁰ = J_s⁰/ρ_m` means a nonzero background flux `J_s⁰` **is** that drain.
- **The uniform limit** must regress to the S11b traction (which carries no background pressure `p⁰`).

Then read: `S11b_SHARED_PHYSICS.md` (as needed), `S11c_decisions.md` (N4, N6, N11b, N11a, N12, N15),
`V3_STEP_PLAN.md:1179`, and the rule-2 twin `directives/_measurements/S11c_a_SHARED_PHYSICS.md`.

## What to check in the repaired sections
1. **§3b / §2d — the constitutive + background repair (the load-bearing one).** Does the affinity/traction
   now use the **perturbation** `δp_s^α` (⛔ not total `p_s^α`)? Is the leading-order background supplied as
   impermeable and at rest (`V_s⁰=0, J_s⁰=0, 𝒜_s⁰=0`) — and is that *self-consistent* (does `J_s⁰=0`
   actually follow from `μ_s⁰=0`, `δp_s⁰=0`, and the closure)? Is any static background pressure `p_s⁰`
   confined to the mechanical support traction `t_hold,s⁰` (⛔ never inside `𝒜`, ⛔ never through `Λ_A(ω)`)?
   ⚠ Does removing the background flux **break** anything the earlier review wanted kept — e.g. the
   true-face-area measure `a_s` on the perturbation flux/traction (it should survive; only the `δa_s·J_s⁰`
   term should vanish)?
2. **§3a — `Σ^α`/`W^α`.** Are the graph thickness `W^α ≡ h_+^α − h_−^α` and slab mass `Σ^α ≡ ρ_4D^α W^α`
   now **defined by equations** (not left as a bare symbol), so T-g/T-h/T-i cannot disagree on the N4
   advective terms `u·∇W_bg`, `u·∇ρ_4D,bg⁰`?
3. **§5a — independence coverage.** Is **T-h** now in the two-route + one-sided-corruption package, with an
   advection-omission corruption? (T-h is where N3 tilt and N4 advection co-occur and can cancel/double-
   count.)
4. **§2b — the freeze.** Is it now stated as a restriction of the triple (`W_bg` varies in every branch; the
   representatives are the two factorizations of `ρ_br=ρ_4D W`; a frozen-`W` pure-N4 background is excluded)?
5. **Repair side effects.** Did any change introduce a **new** defect, an inconsistency between a changed
   section and an unchanged one, a pre-stated result (rule 3: no asserted term/order/channel-survival), or a
   dangling reference (e.g. a symbol renamed in one place but not another)? Check `𝒜_s⁰`/`p_s⁰`/`J_s⁰`
   usages are consistent across §2d, §3b, and the admissibility premise.
6. **Unchanged parts (spot-check only).** §1 fidelity (byte-identical, settled); the anchoring maps; the
   `(η,L_W,σ_W)` bookkeeping; the T-a…T-i coverage; scope deferrals to S11c-b…e; exit semantics; F9
   reservations. Confirm the repair left these intact.

## Physics filter
Report a finding only if it catches a way the physics or the set-up could be wrong. ⛔ "Nothing survives" is
weak evidence: say what you checked and derived, especially on the §3b constitutive form and the `J_s⁰=0`
self-consistency.

## Output
1. Findings, most-serious first — source quote, spec line, and the concrete failure if uncorrected. State
   first whether the **§3b/§2d constitutive+background repair is physically correct**.
2. Any **new** defect the repair introduced, or any cross-section inconsistency.
3. Is the spec now sound enough to hand to the two engines? If not, exactly what remains.
