# Directive pathA_33 — Gate 4: Quadrupole normalization (the `54/5`) → EARNED outgoing fingerprint vs CALIBRATED magnitude

**STATUS:** ✅ DONE — executed 2026-06-26, `QUAD_CALIBRATED` (full tri-review CLEAN after a remediation). Gate 4
of the moving-throat PDE completion ladder (`research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md`),
building on Gate 3 (`pathA_32`, grouped-P2/ℓ=2 isotropy = `ISOTROPY_CALIBRATED`, committed `6711167a`). This is the
gate where the **gravity coupling `G` first enters** and the **outgoing/radiative sector** first appears. Push
memory: `project-moving-throat-pde-push`. Passed the full gauntlet (Codex→GLM→Codex design-review → execute →
tri-review).

**RESULT (2026-06-26) — `QUAD_CALIBRATED`.** On Gate-3's verified isotropic branch, the outgoing ℓ=2 sector delivers
the **EARNED predictive surplus cleanly SEPARATED from the CALIBRATED magnitude**. The outgoing fingerprint
`1/9, 4/81, 1/27` is **DERIVED** from the DtN `Ŷ₂^out=−3/Λ₂^out` Hankel series (NOT hardcoded); the prefactor algebra
follows from `P(ω)=D₀N(ω)/D^cons(ω)²` (the `−2D₂N₀` term); `P0_target_scaling=a⁻⁵`; the derived `χ_Q=1` closes
`m̂₀²P₀=54Gc_s⁵/5a⁵c⁵ ⟺ γ_quad^eff=2G/5c⁵`. **The headline result — the earned/calibrated split:** `54/5 = 2·27/5`
— the **`27` is EARNED** (`derived_in_gate`, from the outgoing fingerprint) while the **`2/5` + `G` + the assembled
magnitude are CALIBRATED** (`external_bridge_input`). Classified by a 4-way **PROVENANCE** partition (NOT `G→λG`
invariance — `54/5` is G-invariant yet calibrated, so invariance alone would mislabel it). `G` is `GENUINE_BLOCKED`:
the PDE delivers the FORM/branch, not `G`. **Dimensional milestone:** the gate carries a genuine, `μ̂₀`-FREE,
able-to-fail dimensional check (`[P₀^phys]=(c_s/a)²·N₀/D₀` must be dimensionless from the sourced `[N₀]=L⁻¹M`,
`[D₀]=L⁻¹T⁻²M`); it caught that the handoff's `P₀=N₀/D₀` silently drops a `(c_s/a)²` factor (a natural-units trap).

**Hard-won (process) — v1 was tri-review-REJECTED then remediated.** The adversarial-with-ablation leg caught two
pass-by-construction sub-controls that fidelity + arbiter both passed: (a) the "non-vacuous" dimensional gate was
**STILL tautological** — it back-solved the FREE carrier `μ̂₀` to force homogeneity (free-carrier tautology); and
(b) the per-probe `self_ablation` field was a **constant, not a re-run**. Both remediated — `μ̂₀`-free dim gate + a
new corrupt-`[N₀]` probe (§3d′) + a real two-verdict `self_ablation` — and re-verified EARNED: corrupting `[N₀]` now
fires `FAIL_DIMENSIONAL` while corrupting `[G]` does not.

**Commit state:** pathA_33 is NOT yet committed (prior committed: Gate 3 = `6711167a`).

**One-line goal.** On Gate-3's verified **isotropic passive/outgoing quadrupole branch**, add the **outgoing
ℓ=2 sector** (the radiative moments `N_{A,n}` + the outgoing-wave boundary condition) and **separate the EARNED
predictive surplus from the CALIBRATED magnitude**: (i) DERIVE the outgoing fingerprint *shape*
`Ŷ₂^out(ω)=1+a²ω²/9c_s²+4a⁴ω⁴/81c_s⁴+i·a⁵ω⁵/27c_s⁵+O(ω⁶)` from the genuine outgoing-wave BC (the rational
coefficients are G-INDEPENDENT and the sharp able-to-fail surplus); (ii) establish the prefactor algebra
`P₀=N₀/D₀`, `P₂`, `P₄` and the target scaling `P0_target_scaling=a⁻⁵`; (iii) prove the target equivalence
`m̂₀²P₀ = 54Gc_s⁵/5a⁵c⁵ ⟺ γ_quad^eff = 2G/5c⁵` and verify its **dimensional consistency**; while (iv)
**labeling `G` and the `54/5` magnitude as the absorbed calibration knob** (`GENUINE_BLOCKED` on overall
normalization). (Scaffold §8.3; handoff §10.3, §11, §12.)

**What this gate is NOT.** It is *not* the numerical branch data `(D_n,N_n)` from the solved nonlinear PDE
(handoff §13.2 = **Gate 6, the wall** — Gate 4 carries the port/branch scalars SYMBOLICALLY, as Gate 3 carried
`B̃_n,Z̃_n`); *not* a derivation of `G` or the literal `54/5` (CALIBRATED — `G` is a labeled knob, the PDE
delivers the FORM/branch); *not* the cross-ℓ unification (Gate 5) or the PN match-back (Phase 4, the downstream
decisive falsifier — `research/4d_*pn*`, do NOT re-derive, memory `project-pn-gravity-ladder`).

**VERDICT LADDER (the script must COMPUTE, not assert, which rung fires; `QUAD_CALIBRATED` is the DEFAULT —
`PASS` would require *deriving* the magnitude, which is `GENUINE_BLOCKED`):**
- `QUAD_PASS` — the full EARNED form (below) holds **AND** the magnitude `54/5` / `G` is **derived from the
  branch**, not calibrated. **Not expected** — `G` is a calibration input (`GENUINE_BLOCKED`); recorded as a
  rung only so the ladder is honest.
- `QUAD_CALIBRATED` *(default, expected)* — the **EARNED predictive surplus is all confirmed**: the outgoing
  fingerprint rational coefficients `1/9, 4/81, 1/27` are **DERIVED** from the genuine outgoing-wave BC (§2.1);
  the prefactor algebra `P₀/P₂/P₄` holds (§2.2); the target scaling `P0_target_scaling=a⁻⁵` (§2.3); the equivalence
  `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵ ⟺ γ_quad^eff=2G/5c⁵` is proven (§2.4) and **dimensionally consistent** (§2.5); the
  outgoing branch is genuinely radiative (passivity EMERGES, not imposed — §2.7) — **while** the magnitude (`G`,
  hence `54/5`) is the labeled calibration knob. The script emits the explicit **earned-vs-calibrated
  partition** (§2.6). Honest toy-model outcome; **does not block the push**.
- `QUAD_FAIL_<reason>`:
  - `FAIL_FINGERPRINT` — the ℓ=2 outgoing-wave BC does NOT yield the rational coefficients `1/9, 4/81, 1/27`
    (the branch / outgoing-sector is wrong). **The sharp earned able-to-fail** (§2.1).
  - `FAIL_NOT_OUTGOING` — the fingerprint was produced by an **imposed phenomenological dissipation** rather
    than a genuine outgoing-wave BC (scaffold §5.3: passivity/outgoing must EMERGE, not be inserted) (§2.7/§3b).
  - `FAIL_SCALING` — the TARGET scaling `P0_target_scaling` ≠ `a⁻⁵` (§2.3; the actual-branch `P₀=N₀/D₀` scaling
    is deferred to Gate 6 unless symbolic `a`-scalings of `N₀,D₀` are supplied).
  - `FAIL_PREFACTOR_ALGEBRA` — `P₀=N₀/D₀`, `P₂`, `P₄` do NOT follow from the low-frequency expansion of the
    prefactor object `P(ω)=D₀·N(ω)/D^cons(ω)²` (NOT `N/D`) (§2.2).
  - `FAIL_EQUIVALENCE` — `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵` and `γ_quad^eff=2G/5c⁵` are NOT the same target under the stated
    `m̂₀`/`γ_quad^eff` definitions (§2.4).
  - `FAIL_DIMENSIONAL` — **`[P₀^phys] = (c_s/a)²·[N₀]/[D₀]` is NOT dimensionless** (computed from the
    independently-sourced port dimensions `[N₀]=L⁻¹M`, `[D₀]=L⁻¹T⁻²M`, `[c_s]=LT⁻¹`, `[a]=L`) (§2.5). **Genuine,
    `μ̂₀`-INDEPENDENT, NON-VACUOUS gate** (the `μ̂₀`-homogeneity check is tautological — `μ̂₀` is a free carrier — so
    it is demoted to a diagnostic); able-to-fail via §3d (drop `(c_s/a)²`) and §3d′ (corrupt `[N₀]/[D₀]`).
  - `FAIL_ENGINE_DISAGREE` — the two engines disagree on a derived expression (distinct from `ENV_BLOCKED`).
- `QUAD_ENV_BLOCKED` — missing Mathematica / license / timeout / infrastructure. Not a physics verdict — re-run.

---

## §0. Scope, framing, and what is able-to-fail (vs what is definitional)

**Toy-model contract (standing).** Calibrations allowed but every claim *earned*, not asserted; goal = a
legitimate brane/bulk/soliton structure supporting EM+GR (memories `project-calibrated-pde-goal`,
`feedback-calibrate-predict-methodology`).

**THE organizing principle of this gate — EARNED surplus vs CALIBRATED magnitude (user directive 2026-06-26).**
The `54/5` quadrupole normalization is the ledger's single sharp target, but `G` is a **calibration input**
(`GENUINE_BLOCKED` — memory `project-pathA-build` Gate 4). The honest, falsifiable value of this gate is the
**predictive surplus that is G-INDEPENDENT**, cleanly separated from the magnitude the calibration absorbs. The
script MUST emit an explicit partition and drive the verdict from it (§2.6):
- **EARNED (`derived_in_gate`, able-to-fail predictive surplus):** the outgoing fingerprint *shape* (the rational
  coefficients `1/9, 4/81, 1/27`); the `P₀/P₂/P₄` algebra; the `P0_target_scaling=a⁻⁵`; the
  `m̂₀²P₀ ⟺ γ_quad^eff` equivalence (incl. the derived `χ_Q`); the dimensional consistency; the genuine-outgoing
  (emergent-passivity) character.
- **CALIBRATED (`external_bridge_input`):** the value of `G`; the PN target `2/5` (→ the assembled `54/5`
  magnitude); the model-level identity to Einstein's `2G/5c⁵` (set by the PN bridge, not derived here).

⚠️ The classifier is **PROVENANCE**, not `G→λG` invariance: the pure number `54/5` is G-invariant yet
calibrated, so invariance alone would mislabel it. The verdict is driven by the provenance partition (§2.6),
with `G→λG` emitted only as a secondary diagnostic.

**Nothing is static (keystone).** `P₀` is the **ω→0 limit of the genuinely ω-dependent outgoing response**
`Y₂^out(ω)`, and the fingerprint is that response's low-frequency Taylor/Frobenius expansion — NOT a frozen
static gain. The imaginary `i·a⁵ω⁵/27c_s⁵` term is the **radiating** part; it must arise from a real
outgoing-wave BC, not an inserted damping (§2.7; memory `project-model-mechanics-corrections` §1; scaffold §5.3
"let passive/outgoing emerge").

**The tautology to avoid (load-bearing for this gate).** The fingerprint coefficients `1/9, 4/81, 1/27` are
**stated in handoff §12.1** — so "verifying" them by comparing to the literal is **theatre**. The genuinely
able-to-fail content is that they **EMERGE from the named DtN object** `Ŷ₂^out=−3/Λ₂^out` (§1a/§2.1), and that a
**wrong BC yields a DIFFERENT tuple**: incoming `h₂^{(2)}` flips ONLY the sign of the imaginary `v₅ω⁵` term (even
coefficients unchanged); standing `j₂` removes the imaginary term / uses a different static slot (§3a). This is
the exact analogue of Gate 1 deriving the D/N pole ladder from the BC (not asserting it). §4 enforces: no
hardcoding the target coefficients.

**Definitional vs earned.** The prefactor *formulas* `P₀=N₀/D₀` etc. (handoff §12.2) and the `N_{A,n}` port
expressions (§10.3) are **definitions** — the earned content is that they (a) follow from the low-frequency
expansion of `P(ω)=D₀·N(ω)/D^cons(ω)²`, (b) combine with the BC-derived fingerprint to give the
`P0_target_scaling=a⁻⁵` form, and (c) close the `γ_quad^eff` equivalence dimensionally.

---

## §1. The object to build

Work at the **linearized isotropic reference** level (the Gate-1→3 regime; Gate-3's lane collapse gives a single
isotropic branch `D_{20,n}=D_{21,n}=D_{22,n}=D_n`). Carry the conservative scalars (`D_n` from Gate 3) and the
outgoing port scalars **symbolically** (the numerical values need the nonlinear solve — §7).

**(1a) Outgoing ℓ=2 fingerprint, DERIVED via the explicit DtN object (handoff §12.1).** Build the ℓ=2 exterior
**Dirichlet-to-Neumann (DtN)** object from the genuine outgoing-Hankel branch — name it explicitly so the
derivation is not gameable: with dimensionless frequency `z = aω/c_s`,
`Λ₂^out(z) = z·h₂^{(1)′}(z)/h₂^{(1)}(z)` (logarithmic-derivative DtN of the outgoing spherical Hankel
`h₂^{(1)}`), and the **static-slot-normalized** response `Ŷ₂^out(z) = −3/Λ₂^out(z)` (the `−3` is the ℓ=2 static
DtN slot `−(ℓ+1)`, so `Λ₂^out(0)=−3` and `Ŷ₂^out(0)=1`). Low-frequency (`z≪1`) expand to obtain
`Ŷ₂^out = 1 + u₂ω² + u₄ω⁴ + i·v₅ω⁵ + O(ω⁶)`; the **derived** coefficients must come out
`u₂=a²/9c_s²`, `u₄=4a⁴/81c_s⁴`, `v₅=a⁵/27c_s⁵`. (These are G-independent — the EARNED surplus. They must EMERGE
from the `−3/Λ₂^out` expansion, NOT be inserted — §2.1, §4.)

**(1b) Outgoing moments + conservative moments (handoff §10.3, §10.4).** For each lane `A` and port `r`, with
`P_{A,r}=Ω_{U,A,r}²g_{W,A,r}+R_{A,r}g_{U,A,r}` and resonance denominator `Δ_{A,r}`:
`N_{A,0}^{(r)}=P_{A,r}²/Δ_{A,r}²`, `N_{A,2}^{(r)}=2P_{A,r}(P_{A,r}S_{A,r}−Δ_{A,r}g_{W,A,r})/Δ_{A,r}³`,
`N_{A,4}^{(r)}=(Δ_{A,r}²g_{W,A,r}²−2Δ_{A,r}P_{A,r}²−4Δ_{A,r}P_{A,r}S_{A,r}g_{W,A,r}+3P_{A,r}²S_{A,r}²)/Δ_{A,r}⁴`;
sum over ports → `N_{A,n}`. On the **isotropic branch** (Gate 3) `N_{20,n}=N_{21,n}=N_{22,n}=N_n` and
`D_{20,n}=D_{21,n}=D_{22,n}=D_n`. The port scalars `(Ω_{U},Ω_{W},g_U,g_W,R,Δ,S)` and `D_n,N_n` are **carried
symbolically** (Gate-6 supplies numbers).

**(1c) Prefactors (handoff §12.2).** `P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`,
`P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³`. ⚠️ These are the Taylor coefficients of the **prefactor function**
`P(ω) = D₀·N(ω)/D^cons(ω)²` (with `N(ω)=N₀+N₂ω²+N₄ω⁴`, `D^cons(ω)=D₀+D₂ω²+D₄ω⁴`) — **NOT** of `N(ω)/D(ω)`
(plain `N/D` gives `P₂=(D₀N₂−D₂N₀)/D₀²`, missing the factor of 2 on `D₂N₀`; the squared denominator is what
produces it). `P₀=N₀/D₀` is the ω→0 value of `P(ω)`. On the isotropic branch the leading odd (radiative)
coefficient is controlled only by `P₀` (handoff §11/§12.2). **⚠️ Units convention (see §2.5):** these `§12.2`
formulas are in **natural units** (`c_s=a=1`); throughout the normalization target (§1d) and `Γ₅`, `P₀` denotes
the **units-restored physical** prefactor `P₀ ≡ (c_s/a)²·(N₀/D₀)` (dimensionless), since the bare `N₀/D₀` carries
`[N₀/D₀]=T²` (`[N₀]=L⁻¹M`, `[D₀]=L⁻¹T⁻²M`). The `(c_s/a)²` is the factor natural units drops.

**(1d) The normalization target + equivalence, via the bridge variable (handoff §12.3; scaffold §8.3).** The
target is `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵`, equivalently `γ_quad^eff=2G/5c⁵`. The equivalence is **NOT direct** — it runs
through the bridge variable `Γ₅` and the odd (radiative) fingerprint coefficient:
- `Γ₅ = χ_Q · P₀ · a⁵/(27c_s⁵)` (the odd `v₅=a⁵/27c_s⁵` slot times the prefactor times `χ_Q`),
- `γ_quad^eff = m̂₀² · Γ₅`,
- where `χ_Q` is the **canonical-branch factor DERIVED from the outgoing DtN fingerprint** (§2.4) — explicitly
  `χ_Q = (DtN-derived radiative coeff `v₅`) / (a⁵/27c_s⁵)`, i.e. the ratio of the derived `v₅` (the imaginary
  `iω⁵` coefficient of `Ŷ₂^out=−3/Λ₂^out`) to the canonical slot `a⁵/27c_s⁵`. It comes out `+1` on the canonical
  outgoing (`h₂^{(1)}`) branch, `−1` on the incoming (`h₂^{(2)}`) branch — a RESULT, **not assumed** (typed
  `χ_Q=1` ⇒ tautology, §4).
- **Source map `m̂₀` (dimensional split, §2.5):** `m̂₀ = μ̂₀ · m̃₀`, where `m̃₀ = 1+O(a²/r²)` is the
  **dimensionless** profile (leading point-particle order `m̃₀→1`) and `μ̂₀` is a **dimensionful scale carrier**
  whose dimension is fixed by homogeneity (§2.5). (Writing `m̂₀=1+O(a²/r²)` alone is dimensionless and cannot
  carry the target's units — the prior-art mismatch.)
- `γ_quad^eff` is extracted from `(K̄₀,K̄₂)` (scaffold §8.3); the physical `P₀=(c_s/a)²·N₀/D₀` is dimensionless
  (§2.5; the bare `N₀/D₀` is `T²`).
The script must (i) state these definitions, (ii) prove `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵ ⟺ γ_quad^eff=2G/5c⁵` THROUGH the
`Γ₅`/`χ_Q` chain with `χ_Q` derived, (iii) check dimensional consistency (§2.5).

---

## §2. Able-to-fail sub-checks (each independently computed)

**§2.1 FINGERPRINT DERIVATION (the EARNED able-to-fail core).** Derive `Ŷ₂^out` from the **named DtN object**
`Λ₂^out(z)=z·h₂^{(1)′}(z)/h₂^{(1)}(z)`, `Ŷ₂^out(z)=−3/Λ₂^out(z)`, `z=aω/c_s` (§1a) — by symbolic series expansion
of the actual spherical-Hankel logarithmic derivative, NOT by referencing the literal coefficients and NOT by an
unconstrained "outgoing solve" (an executor must not be free to pick a different Hankel ratio and back into the
answer). Emit the `h₂^{(1)}` expansion, `Λ₂^out(z)` to `O(z⁵)`, and `Ŷ₂^out=−3/Λ₂^out`. The **derived**
coefficients must equal `u₂=a²/9c_s²`, `u₄=4a⁴/81c_s⁴`, `v₅=a⁵/27c_s⁵`; any mismatch ⇒ `FAIL_FINGERPRINT`.
(Counterfactual §3a: incoming `h₂^{(2)}` flips the sign of the imaginary `v₅ω⁵` term while keeping the even
coefficients; standing `j₂` removes the imaginary term / uses a different static slot — proving the radiative
rational is outgoing-BC-selected, not universal.)

**§2.2 PREFACTOR ALGEBRA.** Independently series-expand the prefactor function `P(ω)=D₀·N(ω)/D^cons(ω)²` (NOT
`N/D`) with `N=N₀+N₂ω²+N₄ω⁴`, `D^cons=D₀+D₂ω²+D₄ω⁴` (symbolic `D_n,N_n`) and confirm the `ω⁰/ω²/ω⁴` coefficients
equal the §12.2 `P₀,P₂,P₄` formulas exactly (incl. the `−2D₂N₀` term, which only arises from the squared
denominator). A mismatch ⇒ `FAIL_PREFACTOR_ALGEBRA`. (Self-check: expanding plain `N/D` would give `−D₂N₀`, not
`−2D₂N₀` — emit both so the squared-denominator object is provably the right one.)

**§2.3 `1/a⁵` TARGET SCALING (`P0_target_scaling`).** The RHS `54Gc_s⁵/5a⁵c⁵` carries `a⁻⁵`; confirm the
**target** value of `P₀` (the magnitude the calibration must hit) scales as `a⁻⁵`. ⚠️ This is the **target**
scaling, NOT the actual branch scaling of `P₀=N₀/D₀`: with `D₀,N₀` carried symbolically (their numerical
`a`-dependence comes from the Gate-6 solve), Gate 4 cannot derive the actual-branch `a`-power. Either (a) emit
`P0_target_scaling = −5` from the RHS form, **and** (b) IF explicit symbolic `a`-scalings of `N₀,D₀` are
supplied, check the actual `P₀` `a`-power matches; otherwise mark the actual-branch scaling **DEFERRED (Gate 6)**.
`P0_target_scaling ≠ −5` ⇒ `FAIL_SCALING`.

**§2.4 EQUIVALENCE (through the bridge variable).** Prove `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵ ⟺ γ_quad^eff=2G/5c⁵` via the
chain `Γ₅=χ_Q·P₀·a⁵/27c_s⁵`, `γ_quad^eff=m̂₀²·Γ₅` (§1d) — NOT as a direct identity. **`χ_Q` must be DERIVED from
the outgoing DtN fingerprint** (`Ŷ₂^out=−3/Λ₂^out`): emit the computed `χ_Q` and confirm it equals `1` on the
canonical outgoing branch as a RESULT, not an assumption (a typed `χ_Q=1` ⇒ tautology, §4). With `m̂₀=μ̂₀·m̃₀`
(§1d) and `γ_quad^eff` extracted from `(K̄₀,K̄₂)` (scaffold §8.3), show one form implies the other. A failure (incl.
a derived `χ_Q≠1` on the canonical branch, or the chain not closing) ⇒ `FAIL_EQUIVALENCE`.

**§2.5 DIMENSIONAL CONSISTENCY (elevated — units RESTORED, NON-VACUOUS, operate on REAL expressions).**
⚠️ **Anti-vacuous mandate (2026-06-26 audit lesson — `feedback-dimensional-consistency-check`):** the Gate-1–3
dim checks were a hand-typed `{M,L,T}` tuple ledger DISCONNECTED from the computed objects, so they could not
fail. This gate's check MUST instead **carry units on the ACTUAL assembled `sympy` expressions** — give
`a,c_s,c,ħ,m,L0,…` real dimensions and run the dimension function ON the computed `P₀=N₀/D₀`, the DtN fingerprint
`Ŷ₂^out`, the `Γ₅` chain, and both sides of `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵` — NOT on a typed tuple table. It MUST carry an
**able-to-fail perturbation** (§3d) and be **dual-engine** (both engines run it).

**Provenance of dimensions:** `[P₀],[c_s],[c],[a]` are DERIVED from the parent-action conventions (handoff
§2.1–§2.2) and the §10 `P₀=N₀/D₀` definition; `[G]` enters through the **external PN bridge** (standard GR
`[G]=L³M⁻¹T⁻²`), NOT from the parent action (G does not appear in the GNLS+Maxwell parent action). Note
`c_s⁵/c⁵=(c_s/c)⁵` is dimensionless. Split `m̂₀=μ̂₀·m̃₀` (`[m̃₀]=1`); `[μ̂₀]` is fixed by demanding homogeneity.

**Where the teeth are (GENUINE, `μ̂₀`-INDEPENDENT gate — supersedes the `μ̂₀`-homogeneity framing, which is
tautological).** ⚠️ **`[μ̂₀]` is a FREE dimensional carrier** (it is not in the parent action; it is a fitted
scale): for ANY `[P₀]`, homogeneity `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵` can be solved by setting
`[μ̂₀]²=[G][c_s⁵]/([a⁵][c⁵][P₀])`. So a "homogeneity_pass" check that back-solves `[μ̂₀]` can **NEVER fail** — it is
DEAD CODE (this is the exact tautology the adversarial audit caught; do NOT make it verdict-bearing). **The genuine
verdict-bearing dimensional gate is `μ̂₀`-INDEPENDENT:**
> **`dimensional_ok` ⟺ `[P₀^phys] = (c_s/a)²·[N₀]/[D₀]` is DIMENSIONLESS** (the zero `{L,M,T}` vector), computed
> from the INDEPENDENTLY-SOURCED port dimensions `[N₀]=L⁻¹M`, `[D₀]=L⁻¹T⁻²M` (parent action / `pathA_22a`) and
> `[c_s]=LT⁻¹`, `[a]=L`.

This is genuinely able-to-fail: `[P₀^raw]=[N₀]/[D₀]=T²`, and `[(c_s/a)²]=T⁻²`, so `[P₀^phys]=1` **only if** both
the port dimensions are right AND the `(c_s/a)²` is present. Corrupt `[N₀]` or `[D₀]` → `[P₀^phys]≠1` → FAIL (§3d′);
drop the `(c_s/a)²` → `[P₀^phys]=T²≠1` → FAIL (§3d). The `μ̂₀`-homogeneity (and the resulting
`[μ̂₀]=L⁻¹T⁻¹M⁻¹ᐟ²`, both sides `=L⁻²M⁻¹T⁻²`) is emitted ONLY as a **non-able-to-fail diagnostic**, explicitly
flagged as such. **Expected resolution:** `[P₀^phys]=1` (dimensionless), so `dimensional_ok=True` — but it is now a
real gate that the §3d/§3d′ mutations flip.

**⚠️ The hidden `(c_s/a)²` factor (corrects a prior slip — this IS the dropped factor the check exists to catch).**
`N₀` and `D₀` are NOT like-dimensioned: the prior SymPy audit `reports/pathA_22a_dimensional_skeleton.md` gives
`[N₀]=L⁻¹M`, `[D₀]=L⁻¹T⁻²M`, so the **raw** `P₀^raw=N₀/D₀` has `[P₀^raw]=T²`. The handoff's `P₀=N₀/D₀` (§12.2) is
in **natural units** (`c_s=a=1`); the units-restored **physical** prefactor carries the factor natural units
drops: `P₀ ≡ (c_s/a)²·(N₀/D₀)` (`[(c_s/a)²]=T⁻²` ⇒ `[P₀]=1`). Use this normalized `P₀` in the target AND in `Γ₅`
(the whole prefactor function takes the same `(c_s/a)²` overall normalization; the `P₀:P₂:P₄` relative structure
is unchanged). This dropped `(c_s/a)²` is exactly the natural-units trap from
`feedback-dimensional-consistency-check`.

**Implementation (the genuine gate is `μ̂₀`-free, so there is no re-solve trap).** `dimensional_ok` is computed
purely from `[P₀^phys]=(c_s/a)²·[N₀]/[D₀]` being the zero vector — `[μ̂₀]` does NOT enter the gate at all (it is a
free carrier), which structurally removes the back-solve tautology. Emit: `[N₀],[D₀],[P₀^raw]=T²,[P₀^phys]`, the
`dimensional_ok` boolean (computed), AND — labeled **"diagnostic, non-able-to-fail (μ̂₀ is a free carrier)"** — the
solved `[μ̂₀]=L⁻¹T⁻¹M⁻¹ᐟ²` and the `m̂₀²P₀` LHS/RHS table. Both engines run it.
Pairs with `feedback-dimensional-consistency-check` (D-dependent `G`; the anti-vacuous regression; the dropped
natural-units factor).

**§2.6 EARNED-vs-CALIBRATED partition (the verdict driver).** ⚠️ **`G→λG` invariance alone is NOT a valid
classifier** — the pure number `54/5` is G-independent yet is the *calibrated magnitude*, so a `g_independent`
test would mislabel it `earned`. Partition every quantity by **PROVENANCE** into four categories:
- `derived_in_gate` (EARNED — computed here from the BC/algebra/dimensions): the fingerprint coefficients
  `1/9, 4/81, 1/27` (and the `27` in `Γ₅`), the `P₀/P₂/P₄` algebra, the `P0_target_scaling=−5`, the equivalence
  chain `Γ₅`/`χ_Q`, the dimensional homogeneity, the derived `χ_Q=1`, emergent passivity;
- `external_bridge_input` (CALIBRATED): `G`, the PN normalization target `2/5` (→ the assembled `54/5`
  magnitude), the Einstein-`2G/5c⁵` identity;
- `deferred_branch_data` (Gate 6): the numerical `D_n,N_n` and port scalars, the actual-branch `a`-scaling;
- `convention` (definitional): unit choices, normalization slots.

The verdict is computed from this provenance partition + the §2.1–§2.5/§2.7 gate booleans: `QUAD_CALIBRATED`
requires ALL `derived_in_gate` items confirmed AND the magnitude correctly classed `external_bridge_input`;
`QUAD_PASS` additionally requires `G`/the magnitude to move from `external_bridge_input` to `derived_in_gate` (not
expected — `GENUINE_BLOCKED`). The classification per quantity must be COMPUTED (provenance traced: does it come
from a BC/algebra solve here, from an external constant, or from deferred numbers?), not typed; the `g_independent`
(`G→λG`) test is emitted as ONE diagnostic but is explicitly **not** the classifier (§3f probes this). A quantity
landing in the wrong category ⇒ `FAIL` of the partition.

**Self-evident decomposition (emit this).** The target magnitude factorizes as `54/5 = 2·27/5` — the **earned
`27`** (the DtN fingerprint's radiative slot `a⁵/27c_s⁵`, `derived_in_gate`) times the **calibrated `2/5`** (the
PN-bridge target `2G/5c⁵`, `external_bridge_input`). Equivalently `γ_quad^eff = 2χ_Q·G/5c⁵` with the derived
`χ_Q=1`. Emitting this decomposition makes the provenance partition self-evident: the `27` is ours, the `2/5`
(and `G`) is the calibration.

**§2.7 NOT-STATIC / passivity EMERGES.** Confirm the fingerprint's radiating (imaginary) term arises from a
genuine outgoing-wave BC, and `P₀` is the `ω→0` limit of the one outgoing response `Y₂^out(ω)` — not a static
gain with damping appended. Building the radiative part by inserting a phenomenological dissipation ⇒
`FAIL_NOT_OUTGOING` (§3b). Static↔dynamic consistency: the `ω→0` limit and the finite-ω expansion are two limits
of one object.

---

## §3. Mandatory counterfactual self-test (able-to-fail probes baked into the script)

Each probe must show the verdict **changes**; a probe that leaves it unchanged ⇒ `FAIL_NOT_ABLE_TO_FAIL`. Every
probe's gate boolean must be **computed from the mutated input** (NOT a typed literal) and carry a self-ablation
(FAIL fires WITH the mutation, vanishes WITHOUT) — the Gate-3 lesson (memory
`feedback-negative-verdict-short-circuit`).

⚠️ **`self_ablation` must be an ACTUAL RE-RUN, not a constant (adversarial-audit defect to avoid).** Each probe's
`self_ablation` field must **re-execute the gate logic with the mutation REMOVED** and record the verdict that
unmutated run actually produced — it must NOT return the global `baseline_verdict`/`fail_suppressed=(baseline==…)`
constant (that is the typed-assertion the gauntlet forbids). The required evidence is two real computed verdicts
per probe: `with_mutation` (= the FAIL) and `without_mutation` (= NOT the FAIL).

**§3a WRONG-BC.** Replace the outgoing (`h₂^{(1)}`) DtN with a **standing-wave** (regular `j₂`) and an
**incoming** (`h₂^{(2)}`) DtN; recompute the fingerprint and assert the **tuple changes in the predicted way**:
incoming flips ONLY the sign of the imaginary `v₅ω⁵` term (the even `u₂,u₄` are unchanged → `χ_Q=−1`); standing
`j₂` gives `Λ_stand(0)=+ℓ=+2` (NOT the outgoing `−(ℓ+1)=−3`), so `Ŷ_stand(0)=−3/2≠1` and the imaginary term is
absent (not the outgoing normalized tuple). This proves the radiative rational `+1/27` is outgoing-BC-selected. If
neither wrong BC perturbs the tuple, §2.1 is vacuous.

**§3b IMPOSED-DISSIPATION.** Build the radiative part by inserting a phenomenological damping term instead of the
outgoing BC; assert this is **detected** (`FAIL_NOT_OUTGOING`) — the imaginary term must come from the BC, not an
inserted sink. Conversely, removing the genuine outgoing BC must remove the imaginary `v₅` term.

**§3c WRONG-SCALING.** Perturb the **target RHS** `a`-exponent (e.g. force the RHS `∝a⁻⁴`); assert
`P0_target_scaling` recomputes to `−4≠−5` and §2.3 fires `FAIL_SCALING`. (Mutating the *actual-branch* `P₀=N₀/D₀`
scaling is not the test here — that is DEFERRED to Gate 6 unless symbolic `N₀,D₀` `a`-scalings are supplied, in
which case mutate those instead.)

**§3d DIMENSIONAL-BREAK (drop `(c_s/a)²`).** Drop the `(c_s/a)²` factor (use raw `N₀/D₀`) → `[P₀^phys]=T²` (not
dimensionless) → the `μ̂₀`-free `dimensional_ok` gate (§2.5) fires `FAIL_DIMENSIONAL`; with `(c_s/a)²` restored,
`[P₀^phys]=1` does NOT fire. ⚠️ Do NOT mutate `[G]` (a free `[μ̂₀]` absorbs it — vacuous) and do NOT route the
gate through `μ̂₀`-homogeneity at all (tautological) — the teeth are `[P₀^phys]` being dimensionless.

**§3d′ DIMENSIONAL-BREAK (corrupt port dimension).** Corrupt `[N₀]` (or `[D₀]`) away from `L⁻¹M` / `L⁻¹T⁻²M`
(e.g. `[N₀]→L⁻¹M·L`) → `[P₀^phys]≠` dimensionless → `dimensional_ok` fires `FAIL_DIMENSIONAL`; the correct
sourced port dimensions do NOT fire. (This is the breadth the §3d-only check lacked — it makes the gate sensitive
to wrong port-response dimensions, not just the `(c_s/a)²` drop.)

**§3e EQUIVALENCE-BREAK.** Perturb a coefficient in `γ_quad^eff` or `m̂₀`; assert §2.4 fires `FAIL_EQUIVALENCE`.

**§3f PARTITION-MISLABEL.** Reclassify an `external_bridge_input` (e.g. `G`, or the PN `2/5` target, or the
assembled `54/5` magnitude) as `derived_in_gate`; assert the provenance partition check (§2.6) catches it —
**by provenance, not by `G→λG` invariance** (the `54/5` is G-invariant yet calibrated, so an invariance-only
test would MISS this; that is exactly the trap §2.6 guards). Conversely, reclassify a `derived_in_gate` item
(the fingerprint `27`) as `external_bridge_input` and confirm that is caught too. Guards the central
earned-vs-calibrated claim from being mislabeled in either direction.

**§3g WRONG-PREFACTOR-OBJECT.** Replace the prefactor object `P(ω)=D₀N(ω)/D^cons(ω)²` with plain `N(ω)/D(ω)`;
recompute `P₂` and assert it comes out `(D₀N₂−D₂N₀)/D₀²` (missing the factor of 2) ≠ the §12.2 formula ⇒
`FAIL_PREFACTOR_ALGEBRA` fires; and that the correct object `D₀N/D^cons²` does NOT fire (self-ablation). Makes the
§2.2 prefactor check a real §3-framework mutation, not only a static comparison.

---

## §4. Anti-tautology firewall + reduction certificate

**Firewall.** The script MUST: (1) **derive** the fingerprint coefficients from the **named DtN object**
`Λ₂^out(z)=z·h₂^{(1)′}/h₂^{(1)}`, `Ŷ₂^out=−3/Λ₂^out` (symbolic Hankel expansion) and EMIT `Λ₂^out(z)` + the
expansion — **forbidden:** hardcoding `1/9, 4/81, 1/27` and "checking" against them, or an unconstrained
"outgoing solve" that lets the executor pick a Hankel ratio (quote the derivation from `−3/Λ₂^out`, not the
literal); (2) **derive `χ_Q` from the DtN fingerprint** (forbidden: typed `χ_Q=1`); (3) compute the
earned-vs-calibrated partition by **PROVENANCE** (`derived_in_gate` / `external_bridge_input` /
`deferred_branch_data` / `convention`), NOT by `G→λG` invariance alone and NOT typed; (4) expand the prefactor
object `P(ω)=D₀N(ω)/D^cons(ω)²` (not `N/D`); (5) carry the §3 able-to-fail probes as computed mutations. If the
fingerprint, `χ_Q`, or the partition is asserted rather than derived ⇒ `FAIL_TAUTOLOGICAL`.

**Reduction certificate (state precisely).** FROZEN-INPUT: Gate-3's isotropic linearized reference; the
symbolic conservative scalars `D_n`; the symbolic outgoing port scalars `(Ω_U,Ω_W,g_U,g_W,R,Δ,S)`/`N_n`; the
calibration knob `G`. COMPUTED: the BC-derived fingerprint `Ŷ₂^out`; the `P₀/P₂/P₄` algebra; `P0_target_scaling=a⁻⁵`;
the `m̂₀²P₀ ⟺ γ_quad^eff` equivalence; the dimensional table; the earned-vs-calibrated partition; the §3
counterfactuals. DEFERRED (§7): the numerical `(D_n,N_n)` from the solved nonlinear branch + the actual branch
selection / passivity emergence on that solution (Gate 6); `G` + the `54/5` magnitude (calibrated). **The
conservative claim:** *given* Gate-3's isotropic branch, the outgoing ℓ=2 sector delivers the G-INDEPENDENT
fingerprint shape + `P₀/P₂/P₄` form + `1/a⁵` scaling + the dimensionally-consistent `γ_quad^eff` equivalence —
the predictive surplus — with the `54/5`/`G` magnitude a labeled calibration; no claim that the magnitude is
derived (that needs Gate 6 + the calibration).

---

## §5. Dimensional consistency (standing pre-numbers step — ELEVATED for this gate)

**Restore units** (memory `feedback-dimensional-consistency-check`). This gate has a **flagged prima-facie
dimensional mismatch** in `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵`: a prior-art audit, assuming textbook 3+1D `[G]` and
`[P₀]=energy/force`, found LHS≠RHS. Resolve it (§2.5): derive `[G],[m̂₀],[P₀]` from the parent-action conventions
(handoff §2.1–§2.2), note `c_s⁵/c⁵` is dimensionless, and either (a) show the formula IS homogeneous under the
model's conventions (identifying where the prior-art guess went wrong — likely `[P₀]` or `[G]`), or (b) report
`FAIL_DIMENSIONAL` if it genuinely cannot be reconciled. Emit a full homogeneity table (every symbol's
dimension + the LHS/RHS check) via SymPy with units restored.

---

## §6. Deliverables

**Dual-engine required** (memory `feedback-dual-engine-required`). **⚠️ Gate-3 lesson — the SECOND engine must
INDEPENDENTLY assemble the HEADLINE quantities** (the BC-derived fingerprint coefficients AND the `P₀/P₂/P₄`
expansion AND the dimensional check), not merely re-check booleans or do the auxiliary probes (Gate-3 v1's `.wl`
was vacuous `x−x` on the headline). The cross-engine comparator must compare the independently-derived fingerprint
coefficients and prefactor expressions, not just pass/fail flags.
- `software/stage1_solver/tools/pathA_33_quadrupole_normalization_sympy.py`
- `software/stage1_solver/tools/pathA_33_quadrupole_normalization.wl`
- `software/stage1_solver/reports/pathA_33_quadrupole_normalization.md` (verdict + which-rung + narrative)
- `software/stage1_solver/reports/pathA_33_results.yaml` (**YAML not JSON**) emitting, at minimum:
  - the DtN object `Λ₂^out(z)` + `Ŷ₂^out=−3/Λ₂^out` expansion and the **derived** fingerprint coefficients
    `u₂,u₄,v₅` from BOTH engines independently (with the `=a²/9c_s², 4a⁴/81c_s⁴, a⁵/27c_s⁵` match booleans,
    computed — the SECOND engine derives them too, not just re-checks);
  - the `P₀,P₂,P₄` formulas + the independent expansion of `P(ω)=D₀N(ω)/D^cons(ω)²` match, AND the `N/D`
    self-check showing `−D₂N₀` vs the correct `−2D₂N₀` (§2.2);
  - `P0_target_scaling` (the RHS `a`-power `=−5`) + the actual-branch scaling status (computed or DEFERRED) (§2.3);
  - the `m̂₀=μ̂₀·m̃₀` / `γ_quad^eff` definitions + the **derived `χ_Q`** + the equivalence proof through
    `Γ₅`/`χ_Q` (§2.4);
  - the dimensional table computed **FROM the real assembled expressions** (recursive walker, NOT a typed
    `{M,L,T}` tuple ledger): `[N₀],[D₀],[P₀^raw]=T²,[P₀^phys]`, and the **verdict-bearing `μ̂₀`-FREE gate
    `dimensional_ok ⟺ [P₀^phys] is dimensionless`** (computed); SEPARATELY the `μ̂₀`-homogeneity
    (`[μ̂₀]=L⁻¹T⁻¹M⁻¹ᐟ²`, LHS/RHS table) emitted **labeled "diagnostic — non-able-to-fail (μ̂₀ free carrier)"**;
    the dim check runs on BOTH engines (§2.5);
  - the **earned-vs-calibrated PROVENANCE partition** (per-quantity class ∈ {`derived_in_gate`,
    `external_bridge_input`, `deferred_branch_data`, `convention`}, computed by provenance; the `G→λG` invariance
    emitted as a separate diagnostic, NOT the classifier) + the **`54/5 = 2·27/5` decomposition** (earned `27` ×
    calibrated `2/5`) (§2.6);
  - the §3 counterfactual results (§3a wrong-BC: incoming flips `v₅` sign / `χ_Q=−1`, standing `Λ(0)=+2`→`Ŷ(0)=−3/2`;
    §3b imposed-dissipation; §3c target-scaling; §3d drop-`(c_s/a)²`→`FAIL_DIMENSIONAL`; **§3d′ corrupt-`[N₀]/[D₀]`
    →`FAIL_DIMENSIONAL`**; §3e; §3f; §3g wrong-prefactor-object→`FAIL_PREFACTOR_ALGEBRA`), each with **two real
    computed verdicts: `with_mutation` (the FAIL) and `without_mutation` (NOT the FAIL)** — the `self_ablation` is
    an ACTUAL re-run, not the global-baseline constant;
  - `verdict`, `which_rung`, `engine_agreement` (max symbolic/numeric delta; tol `<1e-10` symbolic / `<1e-8`
    numeric — `FAIL_ENGINE_DISAGREE` if exceeded), `dim_homogeneity_table`.

A pde_ledger feed note
`research/pde_ledger/notes/stages/moving_throat_pde_pathA_33_quadrupole_normalization_result.md` (Codex-authored,
Claude-verified).

---

## §7. Deferred / out of scope (record, do not silently drop)

- **The numerical branch data** `(D_n,N_n)` (and the port scalars) from the solved nonlinear PDE, the actual
  branch selection, and the on-solution emergence of passivity/outgoing (handoff §13.2) — **Gate 6 (the wall)**.
  Gate 4 carries these symbolically.
- **`G` and the literal `54/5` magnitude** + the Einstein-`2G/5c⁵` identity — **CALIBRATED** (`GENUINE_BLOCKED`;
  the PDE delivers the FORM/branch, not `G` — memory `project-pathA-build`).
- **Cross-ℓ consistency** (ℓ=0/1 return ↔ ℓ=2 quadrupole — the reconciliation's gate) — **Gate 5**.
- **PN match-back** (`research/4d_*pn*` + `research/1pn_orbital_dynamics`) — the decisive downstream falsifier
  (Phase 4); do NOT re-derive the audited PN ladder (memory `project-pn-gravity-ladder`).
- **BC provenance** — inherits Gate 1's `BC_DEPENDENT` for the mouth D/N (a labeled calibration input).

---

## §8. Process

1. **Codex design-review** of THIS directive (`codex exec -c model_reasoning_effort=xhigh`, backgrounded, never
   `timeout`-wrapped) — leading-question / can't-fail-gate / tautology trap hunt; **special focus: the §2.5
   dimensional resolution and whether the §2.1 fingerprint derivation is genuinely BC-driven (not literal-check)**
   (memory `feedback-directive-design-review`). Fold; Codex confirm-pass to GREEN.
2. **GLM** single fresh-perspective pass (user-run, out-of-band). Fold; Codex iterate to GREEN again
   (`feedback-review-ordering-codex-then-glm`).
3. **Codex executes** dual-engine (Claude reviews, Codex codes — `feedback-claude-reviews-codex-codes`); iterate
   to exit 0; Mathematica needs `--sandbox danger-full-access`; ≤2 concurrent `math -script`; scripts wrapped
   `timeout 600`.
4. **Tri-review** (clean agents): arbiter re-run (reliability + refresh outputs); transliteration-fidelity
   (code-vs-equations, **incl. that the SECOND engine independently assembles the fingerprint/P₀**, AND that the
   **§2.5 dim check is NON-VACUOUS** — it carries units on the REAL expressions + has able-to-fail teeth in
   `[P₀]`, not a typed tuple ledger; the 2026-06-26 audit lesson, memory `feedback-dimensional-consistency-check`);
   adversarial-with-ablation (force each verdict-bearing flag/threshold incl. `FAIL_DIMENSIONAL` via §3d + verify
   the FIX too — both directions; memory `feedback-negative-verdict-short-circuit`).
5. **Run the gauntlet via per-gate AGENTS**, not inline (`feedback-offload-review-gauntlet`); orchestrator keeps
   two carve-outs: verify flagged math folds, and one final directive read pre-exec.
6. User gate before commit (`feedback-sequential-audit-chunks`; commit only when the user asks).
