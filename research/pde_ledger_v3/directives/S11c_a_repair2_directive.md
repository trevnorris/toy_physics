# Directive — REPAIR-2 the S11c-a spec (targeted; ⛔ NOT a re-author)

**You (Codex) authored and repaired this spec.** A review leg (Grok) found 3 residual defects, all
introduced or left by the previous *repair directive's* imprecise wording — ⛔ **not** by your
transcription. This directive relays Grok's precise diagnosis and fixes. Apply exactly these three, leave
everything else byte-stable. Reviewed by a fresh Agent + Grok.

## Deliverable
Edit `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` in place (currently HEAD `9365c39c`).

## ⛔ KEEP — both review legs confirmed the constitutive core is CORRECT
⛔ Do **not** reopen: the perturbation-`δp` affinity/traction; the impermeable-rest background
`V_s⁰=J_s⁰=𝒜_s⁰=0`; the `a_s` true-area measure; §1 (byte-identical, fidelity clean); the anchoring maps;
the `(η,L_W,σ_W)` bookkeeping; scope deferrals; exit semantics; F9 reservations. Only the three fixes below.

## FIX A (Grok F1, must — the N4 branch map is missing from the density identity)
⚠ §3a defines `Σ^α ≡ ρ_4D^α W^α` with `ρ_4D^α = ρ_4D,bg⁰(1+θ)` — but the RHS density has **no branch
superscript and no `χ`**, so `ρ_4D^L = ρ_4D^M`, dumping all branch dependence onto `W^α`. On
MATERIAL_ADVECTED × RHOBR-CONSTANT the two readings disagree at linear order in `u` by
`ε·rho_br·u·∂W_bg/W_bg` — the exact N4 disagreement FIX 2 was written to close.

⭐ **Grok's fix — write the §2c branch map into the density identity** (and the same on `μ_s^α`):
```
ρ_4D^α(x,t) ≡ ρ_4D,bg^{0,α}(x,t) (1+θ) ,
   ρ_4D,bg^{0,L} ≡ ρ_4D,bg⁰(x) ,        ρ_4D,bg^{0,M} ≡ ρ_4D,bg⁰(χ(x,t)) ;
μ_s^α ≡ μ_θ^α / ρ_br,bg^{0,α} ,   with ρ_br,bg^{0,α} the §2c-branched ρ_br,bg⁰ .
```
so `Σ^α` and `𝒜_s^α` are uniquely determined per branch and T-g/T-h/T-i cannot disagree on
`u·∇W_bg`, `u·∇ρ_4D,bg⁰`.

## FIX B (Grok F2, must — `t_hold,s⁰` is parked out of `𝒜` then put back into the differentiated traction)
⚠ §3b keeps `p_s⁰` out of `𝒜`/`Λ_A(ω)` (correct), but writes
`t_s^α ≡ −(δp_s^α + Λ_X𝒜_s^α)n̂_s^α + t_hold,s⁰`, and this `t_s^α` feeds `δ_v𝒲_bulk^α`. On a flat face
`t_hold·δ_vx = −p⁰δζ` — the S11b-absent `p⁰` term through a different door: it breaks the §5c regression and
puts an `O(εη)` prestress/clamp channel into T-d that **both engines agree on** (so the cross-engine check
cannot catch it).

⭐ **Fix — take `t_hold,s⁰` OUT of the constitutive traction:**
```
t_s^α ≡ −(δp_s^α + Λ_X(ω) 𝒜_s^α) n̂_s^α          (S11b's law, no t_hold).
```
`t_hold,s⁰` remains **only** in the background support bundle `𝒮_hold⁰`, consumed by the **S11c-b**
admissibility residual (`N15` — the background energetics/prestress operator is S11c-b's). ⛔ S11c-a does
**not** compute a prestress virtual-work channel here; it is deferred to S11c-b. ⭐ Design call
(orchestrator): the prestress×curvature coupling is part of S11c-b's variable-coefficient brane operator, so
S11c-a's traction is exactly S11b's constitutive law. `δ_v𝒲_bulk^α` uses this `t_hold`-free `t_s^α`.

## FIX C (Grok F3, should — T-h's corruption is a no-op on its Eulerian law)
⚠ §5a's T-h corruption is `Σ_E(x(X,t),t) → Σ_E(X,t)` — the **T-g** substitution. T-h's law is **Eulerian**
(`∂_tΣ^α + ∇_x·(Σ^α v) = −Σ_s a_s^α J_s^α`); that substitution never writes `x(X,t)`, so it is a no-op and
the N3/N4 independence test on T-h does not bite.

⭐ **Grok's fix — give T-h a source-level mutation of ITS OWN law on the direct route only:** drop the
advective divergence term (e.g. replace `∇_x·(Σ^α v)` by `Σ^α ∇_x·v`, dropping `v·∇Σ^α`), so the corrupted
route loses the N4 advective content and the two-route residual is forced nonzero.

## Report back (⛔ under 15 lines)
Which of A/B/C you applied and where; flag anything in the fix set that looked wrong. ⛔ No physics
conclusions.
