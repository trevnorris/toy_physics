# Final check — the S11c-a spec after repair-2 (three targeted edits)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`

This Codex-authored spec's **constitutive+background core was already confirmed physically correct by two
prior review passes** (perturbation-`δp` affinity/traction; impermeable-rest background `J_s⁰=0`; drain off;
§1 inlined S11b setup verified term-by-term clean). ⛔ Do not re-audit those — a spot-check suffices. This
final check targets **three small edits** just applied (repair-2) and hunts for any **side effect** they
introduced. It is the last review before the two engines are built, so an error slipping through corrupts
both.

## The three edits to confirm (each closes a defect a prior leg found)
Compare against the prior version with:
`git -C /var/projects/toy_physics diff 9365c39c -- research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`

1. **§3a — density branch map (was: `ρ_4D^α = ρ_4D,bg⁰(1+θ)`, un-branched).** Now
   `ρ_4D^α ≡ ρ_4D,bg^{0,α}(1+θ)` with `ρ_4D,bg^{0,L} ≡ ρ_4D,bg⁰(x)`, `ρ_4D,bg^{0,M} ≡ ρ_4D,bg⁰(χ(x,t))`, and
   `μ_s^α ≡ μ_θ^α/ρ_br,bg^{0,α}`. **Confirm:** the N4 disagreement is now closed — `Σ^α` and `𝒜_s^α` are
   uniquely determined per branch, so T-g/T-h/T-i cannot disagree on `u·∇W_bg` / `u·∇ρ_4D,bg⁰` on the
   material branch. **Side-effect hunt:** is the new `ρ_4D,bg^{0,α}` / `ρ_br,bg^{0,α}` / `μ_θ^α` notation
   *consistent* with §2c's branch map and §2a's `μ_s(y)` definition — no clash, no now-dangling `ρ_4D,bg⁰`
   or `μ_s^α` used un-branched elsewhere?
2. **§3b — `t_hold` removed from the constitutive traction.** Now
   `t_s^α ≡ −(δp_s^α + Λ_X(ω)𝒜_s^α)n̂_s^α` (S11b's law, no `+ t_hold,s⁰`). `t_hold,s⁰` remains only in
   `𝒮_hold⁰` for the reserved S11c-b admissibility residual. **Confirm:** the §5c uniform-limit regression
   now recovers S11b's traction with no `p⁰` residual; no `O(εη)` prestress channel rides in T-d.
   **Side-effect hunt:** is `t_hold,s⁰` still *used* (in `𝒮_hold⁰` / the S11c-b residual) and not left
   dangling? Does any object (T-d, `δ_v𝒲_bulk`) still reference a `t_hold`-bearing `t_s`?
3. **§5a — T-h's corruption now mutates its own Eulerian law.** Now "replace `∇_x·(Σ^α v)` by
   `Σ^α ∇_x·v`, dropping `v·∇_xΣ^α`, only in the direct-route source law." **Confirm:** this forces a nonzero
   two-route residual for T-h (the N3-area / N4-advection independence test now bites). **Side-effect hunt:**
   is the corruption well-defined and one-sided (direct route only)?

## Also (quick)
- Any **new** cross-section inconsistency, dangling symbol, or pre-stated result (rule 3) introduced by the
  three edits?
- Spot-check that the parts NOT in the diff (§0, §1, §2a, §2c, §2d, §3c, §4, §5b, §5c, §6, §7, §8) are
  unchanged and still coherent with the edited lines.

## Method
Read the three edited passages against their sources (S11b for the constitutive law; `S11c_decisions.md` N4,
N11a, N15 for the branch/drain/scope obligations) and the surrounding sections they touch. Quote both sides
for any finding. ⛔ "Nothing survives" is weak evidence — say what you checked.

## Output
1. For each of the three edits: **closed / not-closed**, with the reason. Any **side effect**?
2. Any new defect or inconsistency anywhere from the three edits.
3. **Is the spec now sound to hand to the two engines?** (yes / no + exactly what remains).
