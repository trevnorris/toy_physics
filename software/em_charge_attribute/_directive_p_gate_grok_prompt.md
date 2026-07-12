# Grok compute-verification pass — the native-`Pᵃ` constraint-class gate directive (Codex already `DIRECTIVE_CLEAN`)

You are the **math-verification bookend** on a directive Codex has iterated to `DIRECTIVE_CLEAN` over five rounds. Do NOT execute the gate or design its scripts. **⚠️ Do NOT write or run any scripts (script execution fails in this environment and will kill your run) — verify the Dirac/constraint-algebra claims by ANALYTIC hand-reasoning.** The brackets here are simple enough for that (e.g. `{P²−1, P·π} = {P_aP_a, P_bπ_b} = 2P_aP_bδ_ab = 2P²`). **Your distinct value: independently check the constraint-algebra claims the directive relies on** — the kind of gap a single reviewer misses. This is a toy 4D superfluid analog program (brane = ordered/shear phase = our 3D space; throats = `±w` punctures = charges); falsification is the goal, a computed no-go is first-class.

## Read (in the repo)
- **`software/em_charge_attribute/directive_native_p_constraint_gate.md`** — v5, the directive under verification.
- `docs/em_medium_first_generative_plan.md` §4 — the gate spec (3/3 `SOUND_BUT_REFRAME` panel).
- `docs/conceptual_foundation.md` §3 (MacCullagh shear), §4 (throat/`±w`), open items #1 (easy-axis), #12 (route-a/c).

## The gate in one line
Decide by a FULL nonlinear Dirac–Bergmann analysis whether the native little-arrows sector (hard O(4) sigma `PᵃPᵃ=1` + MacCullagh shear + `χ_B`, most-general symmetry-allowed two-derivative couplings) develops an **ADDITIONAL first-class local Gauss constraint** a conserved charge sources (emergent `U(1)` = native EM), or does not. The pre-existing radial pair is second-class — NOT the test.

## VERIFY these by analytic hand-reasoning (show the brackets/chains; NO scripts)
1. **The central radial-pair claim:** for the hard sigma model with `λ(P·P−1)`, confirm the Dirac chain `π_λ≈0 → C₁=P²−1 → C₂=P·π`, and **compute `{C₁,C₂}`** — is it `2P²≈2≠0` (invariantly second-class), independent of any shear coupling? This is the anti-tautology cornerstone; verify it holds and that shear couplings cannot promote THIS pair.
2. **The "gauging ADDS a new chain, not promotes the radial pair" claim:** on the frozen Stückelberg / gauged-hard-unit control (2-component `φ²=1`, `Dφ=(∂−A)φ`, `A₀` retained), compute the constraint chain and confirm the **MIXED** structure — a NEW first-class Gauss (`π⁰≈0 → ∂·E−ρ`) coexisting with the second-class `φ²=1` radial pair.
3. **The Maxwell controls:** confirm free (ungauge-fixed) Maxwell gives `π⁰≈0` primary → Gauss secondary → first-class + conserved current; and that a **non-conserved** external current fails by **inconsistent preservation** `∂ₜG+{G,H_T}∝−∂_μj^μ≠0` (NOT by the equal-time matrix going second-class).
4. **The MacCullagh + `∇·u≈0` Dirac system:** is `½ρ(∂_tu)²−½μ(∇×u)²` with a local multiplier `σ` enforcing `∇·u≈0` a consistent constraint system, and does curl-only stiffness + that constraint correctly leave the 2 transverse modes (no stray longitudinal)?
5. **Basis finiteness sanity:** under the frozen `SO(3)_brane×Z₂(w→−w)×brane-parity×T×u-shift`, ≤2-derivative, ≤quartic, mod-IBP/derivative-free-redef closure rules — is the cross-coupling enumeration plausibly finite and unique, or is there a hole that makes it infinite/ambiguous?
6. **Verdict-logic + source-test soundness:** any compute-level way the decision table (branch order) or the source-free-first-then-one-linear-`j`-coupling test yields a wrong/ambiguous classification?

## Output
- **VERDICT:** `DIRECTIVE_CLEAN` (math-confirmed) or `DIRECTIVE_NEEDS_FIXES`.
- For each of the 6, a one-line **result** (the actual bracket / chain you derived by hand) + PASS/FAIL against the directive's claim.
- Any catch: the specific fix + why. Concrete over polite. Do not rewrite the directive.
