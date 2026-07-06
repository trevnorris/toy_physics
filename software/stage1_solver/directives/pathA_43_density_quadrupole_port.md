# Directive pathA_43 — The density-mode ℓ=2 radiative-port numerator `N₀`: density-native, or was the vector channel load-bearing? (gravity consolidation gate A1)

**Status:** DRAFT v2 (Codex design-review xhigh + Codex confirm-GREEN + GLM-5.2 tertiary `SOUND_WITH_CONCERNS` all folded — see
changelog §10; pending Codex confirm-to-green on v2 → user gate before execution). **Scope provenance:** `_scratch/A1_setup_scope.md` +
`_scratch/A1_setup_scope_codex_response.md` (the Codex xhigh setup red-team — corrected four sub-questions, crux verdict
picture (ii), minimal able-to-fail, red flags). **Conceptual source (read first):** `docs/conceptual_foundation.md` §3/§4;
`notes/ledger_v2_rebuild_plan.md` §4 (Phase A); `STATUS.md` ▶ RESUME HERE.

**This is the ONE genuinely-new derivation of Phase A (finish gravity).** The gravity sector is otherwise earned
(`pathA_21c` drain force, `pathA_28/29` return/radiation, `pathA_30–34` reduced-closure Gates 1–5, the separate audited PN
corpus `research/4d_*pn*`). This gate re-homes the ℓ=2 quadrupole radiative-port **numerator `N₀`** from the old-ledger
vector `U,W`/`A_w` scaffold onto the density/`c_s` lane — or emits a first-class falsification that it cannot.

---

## 0. What this gate is worth — read BEFORE executing (honesty preregistration)

Unlike `pathA_40` (a near-mechanical CALIBRATED), **this gate has a genuine binary and IS a place falsification can fire.**
The diagnostic (this session, 3-agent + 2-agent fan-outs) established that the ENTIRE exterior of the gravity quadrupole is
already density-native — the DtN fingerprint `1/9,4/81,1/27`, `χ_Q=1`, the prefactor algebra, the `(c_s/a)²` normalization,
the dim slots, and the static-impedance cancellation are all on the `c_s` lane, grep-confirmed free of `A_w` (§1). **The one
object still hosted on the vector port is the ℓ=2 source numerator `N₀`** (`pathA_33` defers it to Gate 6; `pathA_34` makes
the vector binding `N₀ = (Ω_U²g_W + R_mix·g_U)²/(Ω_U²Ω_W²−R_mix²)²` load-bearing). This gate computes `N₀` NATIVELY on the
density lane.

**What this gate EARNS (each independently able-to-fail):** that a density-lane ℓ=2 source/port numerator (a) EXISTS and is
nonzero, **derived from a continuity/interface source equation** (not a relabel — its provenance ANCESTRY carries no vector
coordinate), (b) carries `[N₀]=L⁻¹M` (sourced, μ̂₀-free), (c) lands in the `a⁻⁵` `P0_target` slot **with every `a`-power
traced to continuity/interface geometry + sourced dimensions** (not to a Gate-6-deferred overlap scalar), (d) has the
outgoing radiative sign matching the earned fingerprint (`+i z⁵/27`, `χ_Q=1`), and (e) is **vector-independent** — no
`A_w/F_{μw}/J^w/U/W` in its ancestry/taint set, invariant under ablating the vector *equations from the source graph*. If all
hold → the density lane hosts the quadrupole port, the EM scaffold retires, the diagnostic sliver closes.

**What this gate does NOT earn (honest boundary — do NOT overclaim):**
- **The MAGNITUDE.** `54/5`, the `2/5` Burke–Thorne factor, and `G` stay **CALIBRATED / external** (`G` = `GENUINE_BLOCKED`;
  the PDE gives the FORM/branch, not Newton's `G`). This gate earns EXISTENCE + FORM + vector-independence, NOT the number.
- **The exact branch coupling.** The literal magnitude of the wall→bulk-density coupling requires the nonlinear throat
  interior and is **Gate-6/sim-deferred**. This is a *reduced-closure* gate (like `pathA_32–34`): it decides
  existence/nonzero/dimension/scaling/sign/vector-independence and closure-slot-compatibility — not the solved-branch number.

**⚠ The seductive failure mode here is a rigged POSITIVE** (`feedback-negative-verdict-short-circuit`, inverted). The danger
is NOT a rigged negative — it is **relabeling the old vector 2-port with density NAMES so `N₀≠0` by fiat.** Codex red-flag
(binding): "Replacing `U/W/A_w` names with density labels while preserving the old port formula by fiat." Therefore
**symmetric scrutiny: a `DENSITY_PORT_HOSTED` verdict earns the SAME extra scrutiny as a FAIL.** The gate must SHOW the
coupling is **continuity/interface-derived from physical density variables** (the `pathA_32` wall ℓ=2 mode + the `pathA_29`
bulk-density Helmholtz mode), not a rename — and prove vector-independence by an ABLATION that actually recomputes, not by
inspection.

---

## 1. Earned imports `E` (do NOT re-derive — provenance-tagged; host-agnostic unless noted)

**The already-density-native exterior (REUSE verbatim):**
- **Outgoing DtN fingerprint + `χ_Q=1`** (`pathA_33` `build_fingerprint`; ledger stages 104/105): `z = aω/c_s`,
  `h₂ = j₂+i y₂`, `Λ₂^out = z·d/dz ln h₂`, `Ŷ₂^out = −3/Λ₂^out = 1 + z²/9 + 4z⁴/81 + i·z⁵/27 + O(z⁶)`. Grep-confirmed ZERO
  `A_w/Maxwell/F_{μw}/J^w`. `χ_Q=1` from matching the ω⁵ coefficient of `Ŷ_Q^ret` (renormalized one-pole) to `Ŷ₂^out`.
- **Prefactor algebra** (ledger stage 023 §8): `P₀ = N₀/D₀`, `P₂ = (D₀N₂ − 2D₂N₀)/D₀²`,
  `P₄ = (D₀²N₄ − 2D₀(D₂N₂+D₄N₀) + 3D₂²N₀)/D₀³`. Abstract scalar algebra — host-agnostic.
- **Dimensional normalization** (`pathA_33` Gate-4, μ̂₀-free able-to-fail dim check): `P0_physical = (c_s/a)²·N₀/D₀`, with
  sourced slots `[N₀]=L⁻¹M`, `[D₀]=L⁻¹T⁻²M`, `P0_target_scaling = a⁻⁵`.
- **Static impedance cancellation** (ledger stage 011): `D₀ = K − B₀ − Z₀`; the ℓ=2 compatibility surface eliminates the
  wall stiffness `K`, so `Z₀` drops from the compatibility VARIATION. ⚠ **Precision (Codex red-flag):** `Z₀` still enters
  `P₀` directly through `D₀` — it cancels from the *eliminated compatibility surface* only, NOT from `P₀`. Do not over-read.

**The calibrated closure (magnitude is CALIBRATED — reuse the FORM, not as a derivation):**
- `m̂₀² χ_Q N_Q = 1`, `N_Q = K̄₀/K̄₀^target` (conservative quadrupole normalization defect); with `χ_Q=1`, `m̂₀=1` → `N_Q=1`
  → `m̂₀² P₀ = 54 G c_s⁵/(5 a⁵ c⁵)` ⟺ `Γ̄₅ ≡ m̂₀² P₀·a⁵/(27c_s⁵) = 2G/(5c⁵)` (ledger stages 100/106). `54/5 = 2·27/5`
  (the **27 earned** from the fingerprint; the **2/5·G calibrated**).
- **Invariant consistency target** (from `research/4d_2_5pn`, the corpus's single open item `Γ₅=2G/(5c⁵)`):
  `Γ̄₅ = 9 K̄₂^{5/2}/K̄₀^{3/2}` and `K̄₄ = 4K̄₂²/K̄₀`, with ledger moments `K̄₀=54Gc_s⁵/(5a⁵c⁵)`, `K̄₂=6Gc_s³/(5a³c⁵)`,
  `K̄₄=8Gc_s/(15ac⁵)` (the `K̄₄=4K̄₂²/K̄₀` relation is a bonus closure cross-check — verified to hold by hand).

**The density-lane machinery to EXTEND (the port must be built FROM these physical variables):**
- **`pathA_29` bulk-density outgoing channel:** a linearized scalar continuity/Helmholtz mode `Φ_ℓ(w,r)` on the density
  sound speed `c_s`: `∂_w²Φ + (ω/c_s)²Φ = 0` on the finite slab `0<w<d`; 3D radial `g''+(2/r)g'+κ²g`, `κ²=(ω/c_s)²−m²`;
  `m=0` zero-mode Green function `g ~ 1/r`, whose FLOW `−g' ∝ 1/r²` gives the falloff exponent `p=2` (the `1/r²` is on the
  flow, not on `g` — Codex NIT). The ℓ=0/1 source-moment template: `M0 = ∫S_leak d³x`, `D1_i = ∫x_i S_leak + ∫j_i`.
  ⚠ **`pathA_29` has NO derived ℓ=2 source** (`T2_applied: false`; `Q2`/`raw_kernel2` are inert reuse metadata) — the ℓ=2
  source coupling is exactly what THIS gate computes. `Q₂` in `pathA_28/29` is a FREE anchor symbol, NOT derived density
  physics (Codex red-flag — do not treat it as done). ⚠ **Executor note (GLM NIT-5):** `pathA_28`'s `source_provenance()`
  labels `Q₂`'s `coupling_flag="derived"` — this refers ONLY to the DtN radiation KERNEL `i·a⁵ω⁵/(27c_s⁵)` (reused from
  `research/4d_2_5pn`), NOT to `Q₂`'s value or its wall→bulk coupling. Do not misread it as "the ℓ=2 coupling is already
  derived."
- **`pathA_32` wall ℓ=2 mode:** the grouped-P₂ reduced mode — mass `M₂`, angular stiffness
  `K₂ = ∫[T_w β₂'² + (K_η+6T_Ω)β₂²]`, SO(3)-covariant, on the linearized isotropic reference `R₀(w)`. ⚠ its radial/support
  scalars are CALIBRATED and its **outgoing `N` coefficients are DEFERRED** — so `pathA_32` supplies the wall mode STRUCTURE
  but did NOT compute the outgoing coupling. That coupling is new here.

**Branch / provenance tags (do NOT swap or relabel):**
- The old-ledger `N₀` is VECTOR-hosted: `N_{A0} = P_A²/Δ_A²`, `P_A = Ω_U²g_W + R·g_U`, `Δ_A = Ω_U²Ω_W² − R²`, where
  `U` = brane-like gauge coordinate, `W` = the mixed `A_w/F_{μw}/J^w` coordinate, `R` = U–W mixing, `g_U,g_W` = wall→port
  couplings (ledger stage 023 §4 `Δ_A`, §5 `P_A`/`N_{A0}`; §2 is the Lagrangian — GLM NIT-4). **This gate may NOT rename these to density labels** — it must construct the density
  port from `{q₂ (wall), Φ₂ (bulk density), c_s, a, R₀, β₂, …}` with a coupling DERIVED from projected continuity /
  interface matching, then PROVE the result carries none of `{A_w, F_{μw}, J^w, U, W}`.
- `54/5`, `G`, `2/5` are CALIBRATED anchors — out of the derived set (Codex red-flag: A1 does not derive them).

---

## 2. The crux

The ℓ=2 outgoing radiative-port numerator `N₀` is the moving throat's coupling to the outgoing quadrupole `c_s` wave. The old
ledger hosts it on a vector 2-port (`Ω_U²g_W + R·g_U`). **Decisive question:** what is the density-lane analog once `A_w` is
removed —

- **(i) single-channel density source** — the port collapses to the ℓ=2 projection of the wall motion onto the single bulk
  Helmholtz mode `Φ₂` (the 2-port mixing `R` was an `A_w` artifact with no density analog); OR
- **(ii) a genuine density two-port** `(q₂ wall mode, Φ₂ bulk-density mode)` with a continuity-DERIVED wall–bulk coupling and
  mixing (the SAME `P₀²/Δ₀²` shape but on physical density coordinates, not `A_w`); OR
- **(c) irreducibly vector-hosted** — the nonzero numerator term genuinely needs the `A_w` coordinate `W` (or the mixing `R`
  has no continuity origin); removing the vector coordinates kills or corrupts `N₀`.

**Codex crux verdict (folded):** primary target = **picture (ii)** (the `pathA_32` wall mode + the `pathA_29` bulk-density
mode, coupled), with (i) acceptable only as a Schur-compressed limit of that interface. Picture (c) is a LIVE first-class
negative. The distinction is decidable by reduced algebra; the exact branch magnitude is Gate-6/sim-deferred.

---

## 3. The gate

### 3A. Build the density-lane ℓ=2 port from physical density variables (the derivation — NOT a relabel)

Construct the density-lane ℓ=2 source/port numerator `N₀^{den}` from the physical reduced-closure objects: the `pathA_32`
wall ℓ=2 coordinate `q₂` (with `M₂,K₂`), the `pathA_29` bulk-density Helmholtz mode `Φ₂(w,r)` at `c_s`, and a wall→bulk
coupling **derived from projected continuity / interface matching** across the brane at `w=0` (the same continuity that
produced `M0,D1` at ℓ=0/1, now carried to ℓ=2). Emit `N₀^{den}` symbolically together with its explicit host-variable set.

**Requirement (acceptance, not route):** the derivation must (1) originate the wall→bulk coupling from a continuity/interface
condition or reduced-closure projection — NOT from renaming `g_W/g_U/R`; (2) yield `N₀^{den}` whose host-variable set is a
subset of `{q₂, Φ₂, c_s, a, R₀, β₂, T_w, K_η, T_Ω, m, ρ, K, …}` and contains NONE of `{A_w, F_{μw}, J^w, U, W, Ω_U, Ω_W,
R_mix, g_U, g_W}`. **The route (Schur complement / Green-function DtN matching / Lagrangian projection) is Codex's to design**
— the directive fixes only the origin discipline and the outputs. Report whether the density port is picture (i) or (ii).

**Ancestry / taint machinery (BINDING — Codex design-review BLOCKERs 1,2; the anti-relabel teeth):** the derivation must be
built as an explicit **source graph** — every symbol/equation tagged with its SOURCE (`continuity_interface`,
`pathA_29_bulk`, `pathA_32_wall`, `calibrated_anchor`, `vector_port`, …) — so `N₀^{den}` carries a computed **ancestry/taint
set** (the union of source tags of every factor that actually enters it), reported ALONGSIDE the free-symbol host set. A
symbol-name check is NOT sufficient: a renamed copy of the old vector formula has density-looking free symbols but a
`vector_port` ancestor. Both the origin check (§3B.0) and the vector-independence ablation (§3B.5) operate on this ancestry
graph, not on symbol names.

**⚠ Tag verification — the `continuity_interface` tag must be EARNED, not self-asserted (BINDING — GLM NIT-1; closes the last
anti-rig hole).** Self-reported tags are defeatable by a devious rig: write a NEW equation algebraically equivalent to the
vector port `P_A = Ω_U²g_W + R·g_U`, tag it `continuity_interface`, wire it as `N₀^{den}`'s ancestor — then §3B.0 passes,
§3B.5 removes nothing, and the §5 relabel control (which feeds the *old* structure) does not fire. To prevent this, an equation
may carry the `continuity_interface` tag ONLY if it is derived as a **structured ℓ=2 projection of the SAME
projected-continuity operator** that produced `pathA_29`'s ℓ=0/1 moments (`M0 = ∫S_leak d³x`, `D1_i = ∫x_i S_leak + ∫j_i`) —
i.e. the ℓ=2 source is an `∫Y₂*(θ,φ) S_leak d³x`-type moment of the IDENTICAL continuity/interface condition, NOT an equation
the executor merely labels "continuity." The report must EXHIBIT this ℓ=0→1→2 continuity-operator lineage; a
`continuity_interface` tag without the exhibited lineage is invalid → `FAIL_NOT_DENSITY_DERIVED`. The §5 **fake-continuity
control** tests exactly this (a vector-port-equivalent equation mis-tagged `continuity_interface` must be caught by the lineage
check + the dual-engine method-provenance trace, not pass).

### 3B. The able-to-fail checks (each a COMPUTED pass/fail on `N₀^{den}`)

0. **Origin (`origin_derivation_ok`) — Codex BLOCKER 1:** `N₀^{den}` traces, through the source graph, to a
   continuity/interface source equation (the wall→bulk coupling has a `continuity_interface` ancestor), AND its ancestry/taint
   set contains NO `vector_port` tag. A `N₀^{den}` that (a) has no continuity/interface ancestor (a bare relabel) or (b)
   carries a `vector_port` ancestor → `FAIL_NOT_DENSITY_DERIVED`. This check runs FIRST and is the primary anti-relabel gate.
1. **Nonzero:** `N₀^{den} ≠ 0`, established FROM the derived coupling (not asserted). A computed `N₀^{den} ≡ 0` (selection
   rule / vanishing overlap) is `FAIL_PORT_VANISHES` — a first-class finding (the vector channel was load-bearing).
2. **Dimension:** `[N₀^{den}] = L⁻¹M`, from SOURCED variable dimensions, μ̂₀-free and able-to-fail (the `pathA_33/34` lesson:
   a corrupted SOURCED dimension MUST fail; free carriers must NOT rescue homogeneity).
3. **Scaling + provenance (`a_scaling_provenance_ok`) — Codex BLOCKER 3:** `N₀^{den}/D₀` lands in the `a⁻⁵` slot
   (`P0_physical = (c_s/a)²·N₀^{den}/D₀` carries `P0_target_scaling = a⁻⁵`), AND every `a`-power in `N₀^{den}` traces to
   continuity/interface geometry + sourced dimensions. **If any `a`-power rides on a Gate-6-DEFERRED overlap/coupling scalar,
   that scalar must be PROVEN dimensionless AND `a`-independent** (`pathA_32/33` explicitly defer branch `a`-scaling + port
   scalars to Gate 6 — `pathA_32_results` `outgoing N_A,n deferred`, `pathA_33_results` `port_scalars: gate6_branch_solve`).
   **The `a`-independence proof must be STRUCTURAL, not asserted (GLM NIT-2):** dimensionless ≠ `a`-independent (e.g. `(a/d)²`
   is dimensionless yet `a`-dependent). A deferred scalar counts as `a`-independent ONLY if shown to be a ratio of quantities
   with identical `a`-dimension, or a pure function of dimensionless ratios (`d/a`, profile norms) held FIXED by the reference
   configuration `R₀(w)` — otherwise it is uncertifiable → case-5 `PORT_INCONCLUSIVE`.
   A computed `a`-power `≠ −5` in the slot → `FAIL_PORT_MALFORMED(scaling)`. An `a`-power that CANNOT be certified because it
   rides on an un-provable deferred scalar → `PORT_INCONCLUSIVE_SIM_DEFERRED` (NOT `DENSITY_PORT_HOSTED`).
4. **Radiative sign:** the outgoing sign of the port matches the earned fingerprint (`+i z⁵/27`, `χ_Q=1`) — outgoing/damping,
   not anti-damping/incoming.
5. **Vector-independence (`vector_independence_ok`) — Codex BLOCKER 2:** the SOURCE-GRAPH ablation — remove the vector
   *equations/coordinates* (`{A_w, U, W, Ω_U, Ω_W, R_mix, g_U, g_W}` and every equation tagged `vector_port`) from the input
   graph BEFORE the derivation runs, then RE-DERIVE `N₀^{den}` through the same pipeline. It must be UNCHANGED. (This is
   stronger than post-hoc symbol substitution: it re-runs the derivation on the vector-free graph.) If the re-derivation
   changes or kills `N₀^{den}`, the port was vector-hosted → subsumed under `FAIL_NOT_DENSITY_DERIVED` (§4) with the
   ablation delta as the certificate.

### 3C. Closure-consistency overlay (consistency, NOT magnitude derivation)

Confirm `N₀^{den}` slots into the earned host-agnostic algebra: `P₀ = N₀^{den}/D₀`, `P0_physical = (c_s/a)²·N₀^{den}/D₀`,
and the closure `m̂₀²χ_Q N_Q = 1` with `χ_Q=1` reproducing the earned/calibrated split — with `G`, `2/5`, and the literal
`54/5` held CALIBRATED/external (state each as an input, never as a derived output). Cross-check the invariant
`K̄₄ = 4K̄₂²/K̄₀` on the density-mode moments (a closure consistency witness feeding A3's 2.5PN match-back;
`Γ̄₅ = 2G/(5c⁵)` is the calibrated target, NOT derived here).

---

## 4. Verdict grammar — a STRICT ordered decision table (mutually exclusive, complete)

Evaluate in THIS order; the first matching case is the verdict:

1. `FAIL_NOT_DENSITY_DERIVED` — the origin check (§3B.0) OR the vector-independence source-graph ablation (§3B.5) fails:
   `N₀^{den}` has no continuity/interface ancestor (a bare relabel), OR its ancestry/taint set carries a `vector_port` tag, OR
   re-deriving on the vector-free source graph changes/kills it. **The port is not genuinely density-derived — the vector
   channel was load-bearing → the EM question reopens.** Carries the ancestry/taint set + the ablation delta. [first-class
   falsification; the primary anti-relabel / anti-rig verdict — subsumes the old `FAIL_VECTOR_LOAD_BEARING`]
2. `FAIL_PORT_VANISHES` — the (genuinely density-derived) `N₀^{den}` computes to `≡0` (selection rule / vanishing overlap).
   The density lane cannot source the ℓ=2 port → the EM question reopens. Carries the vanishing overlap/selection-rule
   certificate. [first-class falsification]
3. `FAIL_PORT_MALFORMED` — `N₀^{den}` is nonzero and density-derived but fails ≥1 of {`[N₀]=L⁻¹M`, a COMPUTED `a`-power `≠−5`
   in the slot, outgoing sign}. A **characterized departure** (the density lane hosts *a* port, but not the earned one);
   carries which check failed. (An uncertifiable — not wrong — `a`-power routes to case 5, not here.)
4. `DENSITY_PORT_HOSTED` — ALL checks pass: §3B.0 origin (continuity ancestor, no `vector_port` taint), nonzero,
   `[N₀]=L⁻¹M`, `a⁻⁵` slot with certified `a`-provenance, outgoing sign, vector-independent under source-graph ablation, and
   slots into the closure. **The ℓ=2 quadrupole port is density-native; the EM scaffold retires.** Magnitude CALIBRATED
   (`G,2/5,54/5` external); exact branch coupling Gate-6/sim-deferred. **Expected landing IF the physics is clean — but earns
   the SAME extra scrutiny as a FAIL** (§0 anti-rig): the report must SHOW the coupling is continuity-derived (ancestry/taint
   set clean) and the source-graph ablation genuinely re-derives.
5. `PORT_INCONCLUSIVE_SIM_DEFERRED` — the terminal fallback: existence-or-obstruction cannot be certified at the
   reduced-closure level (e.g. an `a`-power or the coupling magnitude genuinely needs the Gate-6-deferred nonlinear throat
   interior). **Non-dodge certificate (BINDING — Codex BLOCKER 4):** an inconclusive verdict MUST carry a per-predicate table
   — for EACH of {origin, nonzero, dimension, scaling, sign, vector-independence}: its evaluated result, WHY it is undecidable
   at reduced closure, and the EXACT deferred object. **Inconclusive is FORBIDDEN if the reduced computation has already
   shown** zero coupling (→ case 2), vector ancestry / ablation-dependence (→ case 1), a malformed dimension / wrong computed
   `a`-power / wrong sign (→ case 3), or a missing continuity origin (→ case 1). It may only fire where a predicate is
   genuinely *undecided*, never where a FAIL is already reachable.

**Rider attachment (every non-FAIL verdict carries, each CITING its producing computation — no decorative constants):**
- `port_picture: {i single-channel | ii two-port(q₂,Φ₂)}` with the derived coupling expression;
- `closure_slot: P0_physical=(c_s/a)²N₀^{den}/D₀; a⁻⁵ ✓; K̄₄=4K̄₂²/K̄₀ ✓` (← §3C);
- `magnitude: CALIBRATED{54/5, 2/5, G=GENUINE_BLOCKED}; exact_branch_coupling: SIM_DEFERRED(Gate-6)` (scope-boundary caveat,
  fixed statement, labelled as such);
- `feeds: v2-ledger Part-II gravity; 2.5PN open item Γ̄₅=2G/(5c⁵) (research/4d_2_5pn), A3 match-back`.

---

## 5. Controls (able-to-fail, rig-locked; `feedback-negative-verdict-short-circuit`)

**Discipline (BINDING — anti-rig):** every control fixture MUST mutate the input *equations / physical source objects* and
then **recompute** the verdict through the SAME predicates the production run uses. NO fixture may set its own answer — no
`n0_vanishes=True`, `vector_hosted=True`, typed verdict boolean, or expected-label short-circuit (the Gate-3/4/pathA_38
rigged-positive pattern). A `DENSITY_PORT_HOSTED` that cannot be flipped by any mutation is a rig.

**The verdict-flipping counterfactuals (each MUST flip, recomputed from mutated source facts):**
- **Vector-hosted control (the decisive discriminator):** feed the OLD vector port `P_A = Ω_U²g_W + R·g_U` (genuine `A_w/U/W`
  structure, tagged `vector_port` in the source graph) → the §3B.0 origin check MUST detect the `vector_port` ancestor AND the
  §3B.5 source-graph ablation MUST change/kill it → recompute `FAIL_NOT_DENSITY_DERIVED`. This proves the gate distinguishes
  density-derived from vector-hosted, and is the anti-relabel teeth.
- **Relabel-rig control (anti-fiat — targets the §0 red-flag / BLOCKER 2 directly):** feed the SAME old vector structure with
  every symbol RENAMED to a density-looking label (`Ω_U→ω_wall`, `g_W→g_ρ`, `R→r_mix`, …) but with NO continuity/interface
  ancestor wired in → the ancestry/taint machinery MUST still trace the `vector_port` origin tag (the rename does not change
  the source-equation tag) and recompute a SPECIFIC `FAIL_NOT_DENSITY_DERIVED` (NOT `DENSITY_PORT_HOSTED`). This is the
  rigged-positive short-circuit; a name-only check that passes this fixture is the rig.
- **Fake-continuity control (GLM NIT-1 — the deterministic-tag teeth):** feed an equation algebraically EQUIVALENT to the
  vector port `P_A` but MIS-TAGGED `continuity_interface` (the devious re-tag: origin check would pass, ablation would remove
  nothing) → the §3A lineage check MUST reject the tag (no exhibited ℓ=0→1→2 `∫Y₂*S_leak`-type lineage from the `pathA_29`
  continuity operator) and recompute `FAIL_NOT_DENSITY_DERIVED`. This proves the `continuity_interface` tag is EARNED (lineage
  + dual-engine method-provenance), not self-asserted — closing the last relabel hole.
- **Ablation-isolation control (GLM NIT-3 — the ablation has independent teeth):** feed a derivation that uses `vector_port`
  equations as INTERMEDIATE steps whose contribution SURVIVES into `N₀^{den}` (a Schur step that does NOT fully eliminate the
  vector dependence), constructed so the FINAL free-symbol host set looks vector-free → the §3B.5 source-graph ablation MUST
  still detect the surviving intermediate dependence (removing the vector equations changes `N₀^{den}`) and recompute
  `FAIL_NOT_DENSITY_DERIVED`, EVEN THOUGH a naive final-expression origin check might have passed. (Contrast: a Schur step that
  GENUINELY eliminates the vector dependence — ablation leaves `N₀^{den}` unchanged — is the legitimate picture-(i) limit and
  must NOT fire.)
- **Zero-coupling control:** set the DERIVED wall→bulk continuity coupling to 0 (mutate the interface condition, not a flag) →
  MUST recompute `FAIL_PORT_VANISHES`.
- **Dimensional control:** corrupt a SOURCED dimension (e.g. `[Φ₂]`, `[β₂]`, or the coupling) by one power of `L` → MUST
  recompute `FAIL_PORT_MALFORMED(dimensional)`; a free-carrier corruption must NOT (the `pathA_34` sourced-vs-free lesson).
- **Sign control:** swap the outgoing→incoming DtN boundary (`+i → −i`, incoming) → MUST recompute `FAIL_PORT_MALFORMED(sign)`
  / mismatched radiative sign vs the earned `+i z⁵/27`.
- **Scaling control:** perturb the `a`-power of the coupling to a WRONG COMPUTED value (break `a⁻⁵`) → MUST recompute
  `FAIL_PORT_MALFORMED(scaling)`.
- **Deferred-scalar control (BLOCKER 3 teeth):** route an `a`-power through a synthetic Gate-6-deferred overlap scalar whose
  `a`-independence CANNOT be proven → MUST recompute `PORT_INCONCLUSIVE_SIM_DEFERRED` (not `DENSITY_PORT_HOSTED`), with the
  per-predicate non-dodge table naming that scalar as the deferred object. Conversely a deferred scalar PROVEN dimensionless +
  `a`-independent must NOT force inconclusive (proves the certificate discriminates).

**Baseline PRODUCTION check (NOT a control):** the real derived density port → compute `N₀^{den}`, its host set, all §3B
checks, the verdict, and the §3C closure slot. If the computed result differs from the §9 prior, **preserve the result and
trigger extra scrutiny** — a computed FAIL is NOT a failure of this check; it is a first-class finding.

---

## 6. Dual-engine (BINDING; `feedback-dual-engine-required`, Stage-4 lesson)

Every computed claim (the derived coupling, `N₀^{den}`, the dimension, the `a⁻⁵` scaling, the sign, the vector-independence
ablation, each control) is **dual-engine**: SymPy + a genuinely INDEPENDENT Mathematica route, engine-agreement asserted per
quantity. **The second engine must DERIVE, not echo** (the Stage-4 vacuous-reassertion pattern; the `MATHEMATICA_MIRROR_POLICY`
this ledger will inherit). The independence must live in the METHOD — e.g. SymPy assembles the reduced-closure Lagrangian and
takes the wall→bulk coupling by Schur elimination / projection; Mathematica derives the SAME coupling by an independent route
(Green-function/DtN interface matching, or a different elimination order), and independently recomputes the ablation.
Mathematica MUST NOT consume SymPy's `N₀^{den}` expression or its pass/fail booleans. **Method-provenance evidence (Codex
NIT):** each engine emits its own derivation trace (the intermediate coupling/expression it built, before agreement is
checked) so the review can VERIFY the Mathematica route derived `N₀^{den}` independently rather than transliterating SymPy's
result — engine-agreement on a value the second engine merely re-typed is the Stage-4 vacuous-reassertion rig.

The orchestrator independently re-runs both engines as arbiter + a transliteration-fidelity audit + a clean
adversarial-with-ablation agent (per `docs/development_pipeline.md`).

---

## 7. Deliverables

- `directives/pathA_43_density_quadrupole_port.md` (this file, post-gauntlet).
- `tools/pathA_43_density_quadrupole_port_{sympy.py,.wl}` — dual-engine (the density-port derivation + the §3B checks + the
  §3C closure slot + all §5 controls; both exit 0; engine-agreement asserted; print-only/assert-zero idioms so the reshape
  into the v2 ledger is minimal).
- `reports/pathA_43_density_quadrupole_port.md` (verdict on line 1) + `reports/pathA_43_..._results.yaml` — the earned imports
  `E`, the derived `N₀^{den}` + host set + port picture, the §3B check results, the vector-independence ablation delta, the
  §3C closure slot + `K̄₄=4K̄₂²/K̄₀` cross-check, the controls, the CALIBRATED/SIM-DEFERRED scope tags, provenance.

---

## 8. Review plan (the gauntlet; `docs/development_pipeline.md`)

1. **This directive** → Codex design-review (xhigh) → fold → GLM-5.2 tertiary (via opencode; resume `-c "continue"` if it
   times out) → Codex confirm-to-green → **user gate** before execution.
2. **Execution:** Codex codes + runs + iterates to exit 0 + dual-engine; Claude reviews only.
3. **After the run:** orchestrator arbiter re-run (both engines) + transliteration-fidelity audit + clean
   adversarial-with-ablation agent (the §5 counterfactuals — vector-hosted, zero-coupling, dimensional, sign, scaling,
   relabel-rig, each recomputed from mutated source facts with no direct-boolean short-circuit — are the adversarial focus;
   the anti-relabel and the rigged-`DENSITY_PORT_HOSTED` are the primary targets) → verdict → **user gate**.
4. **Never alter the calibrated process unilaterally** (`feedback-never-alter-calibrated-process`): if a step fails, HALT +
   bring to user. Math-level disputes (a constant, a tolerance, a value-class) → Claude+Codex resolve + agree; escalate to
   user only if a fix changes the conceptual nature.

---

## 9. Scope boundary + honest priors

- **This directive = the density-lane ℓ=2 port numerator `N₀^{den}` ONLY** — existence, form, dimension, scaling, sign,
  vector-independence, closure-slot compatibility. It does NOT derive the magnitude (`54/5,2/5,G` CALIBRATED), does NOT solve
  the moving-throat PDE from first principles (the particle limit stays a controlled reduction; the exact branch coupling is
  Gate-6/sim-deferred), and does NOT re-derive the PN ladder (`research/4d_*pn*` is cited, not redone).
- **Feeds:** the v2 ledger Part II (gravity); the `research/4d_2_5pn` corpus's single open item `Γ₅ = 2G/(5c⁵)` (this gate
  supplies the density-native port the 2.5PN "conditional theorem" was conditioned on, at the reduced-closure level); A3's
  invariant consistency check (`K̄₄=4K̄₂²/K̄₀`, `Γ̄₅=9K̄₂^{5/2}/K̄₀^{3/2}`).
- **Honest prior (recorded before executing) — PLAUSIBLE BUT WEAK (Codex NIT):** `DENSITY_PORT_HOSTED` is *plausible* — the
  physical frame demands gravity radiation ride the medium's own density ripple at `c_s`, and the bulk Helmholtz mode `Φ₂` is
  the natural outgoing channel (its DtN fingerprint is already earned), so the wall's quadrupolar motion *should* couple to it
  by continuity. **But the prior is WEAK, not confident:** the ONLY concrete `N₀` that exists today is still vector-hosted
  (`pathA_34` `N0_from_port`, `free_N0_symbol_used_in_verdict:false`), `Q₂` is merely a live anchor symbol (`pathA_28`), and
  `pathA_32` DEFERRED the outgoing coupling — so no one has actually computed this. A vanishing overlap (`FAIL_PORT_VANISHES`),
  a relabel/vector-ancestry (`FAIL_NOT_DENSITY_DERIVED`), or an uncertifiable `a`-provenance (`PORT_INCONCLUSIVE`) are all live
  outcomes that would reopen the EM question or bound the claim. A clean `DENSITY_PORT_HOSTED` earns extra scrutiny, not
  celebration (`feedback-falsification-is-the-goal`: a clean "it all works" is suspicious). **Calibrate-predict honesty**
  (`feedback-calibrate-predict-methodology`): the magnitude is ABSORBED (calibrated on `G`, zero surplus from this gate);
  report the earned-vs-calibrated split plainly; the surplus lives in the held-out predictions + the 2.5PN/4PN
  consistency, not here.

---

## 10. Changelog

- **v1 → v2 (GLM-5.2 tertiary folded; verdict was SOUND_WITH_CONCERNS, no BLOCKERs).** GLM independently hand-verified both
  invariant relations (`K̄₄=4K̄₂²/K̄₀`, `Γ̄₅=9K̄₂^{5/2}/K̄₀^{3/2}=2G/(5c⁵)`) and verified every cited claim against the repo
  (17-row table, all ✓). 5 NITs folded: (NIT-1, the material one) the `continuity_interface` tag was self-reported and
  defeatable by a devious re-tag (a vector-port-equivalent equation mis-tagged continuity) → added a BINDING **tag-verification
  lineage requirement** (a `continuity_interface` tag is valid ONLY with an exhibited ℓ=0→1→2 `∫Y₂*S_leak`-type lineage from the
  `pathA_29` continuity operator) + a **fake-continuity control** (§3A, §5); (NIT-2) `a`-independence proofs for deferred
  scalars must be STRUCTURAL, not asserted (dimensionless ≠ `a`-independent) (§3B.3); (NIT-3) added an **ablation-isolation
  control** proving §3B.5 has teeth independent of the origin check (vector eqns as surviving intermediates → FAIL) (§5); (NIT-4)
  `stage 023` citation fixed §2,5→§4,5 (§1); (NIT-5) executor note that `pathA_28`'s `Q₂` `coupling_flag="derived"` refers to the
  DtN kernel, not `Q₂`'s coupling (§1). Pending Codex confirm-to-green on v2.
- **v0 → v1 (Codex design-review xhigh folded; verdict was NEEDS_REVISION).** BLOCKERS: (1) origin discipline was not a formal
  verdict predicate → added `origin_derivation_ok` (§3B.0) + an ordered `FAIL_NOT_DENSITY_DERIVED` verdict (§4 case 1) BEFORE
  `DENSITY_PORT_HOSTED`, requiring a derivation trace from continuity/interface source equations. (2) name-ablation could be
  gamed by a renamed copy of the old vector formula → added ANCESTRY/TAINT machinery (a tagged source graph; `N₀^{den}` carries
  a taint set reported alongside the host set), upgraded §3B.5 vector-independence to a SOURCE-GRAPH ablation (remove vector
  *equations* before re-deriving, not post-hoc symbol substitution), and made the relabel-rig control route to a SPECIFIC
  `FAIL_NOT_DENSITY_DERIVED` (§5). (3) the `a⁻⁵` scaling could overclaim (branch `a`-scaling is Gate-6-deferred in
  `pathA_32/33`) → added `a_scaling_provenance_ok` (§3B.3): every `a`-power must trace to continuity/interface geometry +
  sourced dims; an `a`-power riding an un-provable deferred scalar → `PORT_INCONCLUSIVE`, not `HOSTED` (+ a deferred-scalar
  control, §5). (4) `PORT_INCONCLUSIVE` lacked a non-dodge certificate → now MUST carry a per-predicate table (result / why
  undecidable / deferred object) and is FORBIDDEN where any FAIL is already reachable (§4 case 5). NITs folded: `pathA_29`
  zero-mode wording (`g~1/r`, flow `1/r²`, `p=2`) (§1); honest prior softened to PLAUSIBLE-BUT-WEAK (§9); dual-engine
  method-provenance trace evidence required (§6). No issue found with the `54/5`/`2/5`/`G` calibrated boundary (Codex confirmed
  correct). Pending Codex confirm-to-green → GLM-5.2 tertiary.
- **v0 (Claude draft from the Claude↔Codex setup scope).** Built from `_scratch/A1_setup_scope.md` + the Codex xhigh setup
  red-team (`_scratch/A1_setup_scope_codex_response.md`): corrected four sub-questions (density source object; port anatomy
  after removing `A_w`; able-to-fail negative; closure consistency not magnitude), crux verdict picture (ii), the minimal
  able-to-fail computation (the reduced ℓ=2 density-lane numerator), and the red flags (no `Q₂`-as-derived, no relabel-by-fiat,
  no `T2`-as-derived, no `54/5`/`G` claim, `Z₀`-still-enters-`P₀`, no free-carrier dim pass, no ungrounded selectors).
  Pending Codex design-review.
