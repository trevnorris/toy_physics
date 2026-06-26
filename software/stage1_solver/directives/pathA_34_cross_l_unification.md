# Directive pathA_34 — Gate 5: scalar/dipole side-conditions + CROSS-ℓ unification (the ℓ=0/1 return ↔ ℓ=2 quadrupole consistency)

**Ladder position:** Gate 5 of the moving-throat-PDE completion ladder
(`research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md`). Gates 1–4 DONE & CALIBRATED
(`pathA_30/31/32/33`). This is the **first decisive cross-ℓ gate** — the survival test for the model's one concrete
GR-departure prediction (the bounded ℓ=0/1 `c_s`-radiation residual `∝ε₀`).

**Template:** mirror `directives/pathA_33_quadrupole_normalization.md` (the hardened pattern: COMPUTED-not-asserted verdict
ladder, dual-engine HEADLINE, computed counterfactual probes with **real two-verdict `self_ablation`**, the `μ̂₀`-FREE /
free-carrier-independent dimensional check). The dim check here must be **genuine from the start** (NOT the Gate-1–3 vacuous
typed-tuple debt — see `pathA_dimcheck_retrofit_gates1to3.md` for what NOT to repeat).

**Teeth design provenance:** the verdict-bearing teeth below were settled in a Claude+Codex (xhigh) design discussion
(`software/stage1_solver/_scratch/pathA_34_gate5_teeth_design_discussion.md` +
`software/stage1_solver/_scratch/pathA_34_gate5_teeth_codex.log`), DECISION STATUS = CONVERGED; then design-reviewed by Codex
(xhigh, `software/stage1_solver/_scratch/pathA_34_design_review_codex.log` = SOUND-WITH-CONCERNS) and the 8 issues folded.
The central refinement: a naïve "one `S_return` parameter set exists" test is **TAUTOLOGICAL** — the port-kernel symbols
`Ω_{U,A,r}, Ω_{W,A,r}, R_{A,r}, g_{U/W,A,r}` (handoff §10.2–10.3) are **already-reduced ℓ=2 lane data**, not ℓ-independent knobs;
letting them float per-ℓ trivially passes. The sharedness lives **one level up**, in the projection-blind kernel **generator** /
`S_return` functional / projection rule fixed **before** ℓ-projection. Hence the decisive object is **"fix a single
projection-blind generator + admissible branch on ℓ=2 (keeping `QUAD_CALIBRATED`); then ℓ=0/1 are TARGET-BLIND PREDICTIONS,"**
guarded by a **provenance/rank audit** (no new `ε₀`, `ε₁`, or ℓ-gain knob introduced after the ℓ=2 branch is fixed).

---

## VERDICT LADDER (the script must COMPUTE, not assert, which rung fires; `CROSS_L_RESIDUAL_PREDICTION` is the DEFAULT pass)

**PASS — `CROSS_L_RESIDUAL_PREDICTION`:** a single projection-blind generator + admissible branch (i) keeps ℓ=2
`QUAD_CALIBRATED`, and (ii) makes the **target-blind** ℓ=0/1 return give **finite, nonzero, bounded** effective residuals
matching pathA_29 in **form / sign / order** (magnitude deferred):
- `A0 ≃ i·a·(ω/c_s)·M0 · ε0_eff/(1+ε0_eff)`  (leading, ω→0; raw outgoing order `ω¹`)
- `A1 ≃ i·a³·(ω/c_s)³·D1 · ε1_eff/(2(1+ε1_eff))`  (leading; raw outgoing order `ω³`)
  with `ε0_eff, ε1_eff` **computed from the shared generator** (not free flags), finite and `>0` (bounded DC transfer
  `T_ℓ(0)=1/(1+ε_eff)∈(0,1)`).

**FAIL rungs** (the **physics/logic** rungs are each able-to-fail; §3 probes prove each fires. The **infrastructure** rungs
`FAIL_ENGINE_DISAGREE` / `ENV_BLOCKED` fire from the cross-engine comparator / environment, not a §3 mutation — see §4):
- `FAIL_DECOUPLED` — the construction admits **block-separable** ℓ=0/1 vs ℓ=2 knobs (ℓ=2 still hits `P0` while `ε0` is dialed
  arbitrarily) ⇒ not one return law. (§2.1 rank/provenance audit.)
- `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` — the generator is shared, but the ℓ=2 constraints leave a **null direction** that
  changes the scalar/dipole DC return `T0(0)/T1(0)` **without** changing `P0=N0/D0`, and no PDE branch selector removes it ⇒
  ℓ=0/1 is not actually predicted. (§2.1.)
- `FAIL_EPSILON_MISMATCH` — the target-blind ℓ=0/1 prediction has **wrong sign / wrong order / not boundable** as a finite
  positive `ε_eff` (e.g. an anti-localizing / no-zero-mode / wrong-DC-transfer return) ⇒ disagrees with pathA_29's residual
  class. (§2.3.)
- `FAIL_OVERCANCEL` — the ℓ=2-compatible return forces `T0(0)=T1(0)=1` / `ε_eff=0` (clean cancellation), **erasing** the
  advertised leading residual and the GR-departure prediction. **Suspicious, not a success** (`feedback-falsification-is-the-goal`). (§2.3/§3d.)
- `FAIL_NO_CONSISTENT_RETURN` — a single shared generator is proven, but **no admissible branch** satisfies both ℓ=2 `P0` and
  the ℓ=0/1 residual class. (A genuine no-go; welcome.) (§2.4.)
- `FAIL_QUAD_REGRESSION` — the cross-ℓ rebuild **breaks Gate 4** (incoming/standing rather than outgoing DtN, wrong prefactor
  algebra `P(ω)=D₀N/D^cons²`, or loss of the `1/9, 4/81, 1/27` fingerprint / `χ_Q=1`). (§2.5.)
- `FAIL_DIMENSIONAL` — a verdict-bearing residual/transfer expression is not dimensionally homogeneous under sourced
  dimensions (free-carrier-INDEPENDENT check; §2.6/§3f).
- `FAIL_TAUTOLOGICAL` — the cross-ℓ relation, the shared-generator claim, or any fingerprint is **asserted** rather than
  derived; or a probe leaves the verdict unchanged. (§4.)
- `FAIL_ENGINE_DISAGREE` — the two engines disagree on a derived expression (distinct from `ENV_BLOCKED`).

---

## §0. Scope, framing, and what is EARNED vs CALIBRATED/DEFERRED

This gate is the **track-3 falsification** that `pathA_28` deferred: `pathA_28_cancellation_condition.yaml` states its own
`R0=−M0`/`R1=−D1` is "algebraic bookkeeping… an identity (x−x)… `cancellation_possible` is a parameter/literal flag, not a
computed quantity… real falsification lives in track 3: whether an admissible return can **actually** deliver R0=−M0, R1=−D1."
**Gate 5 makes `ε0_eff, ε1_eff` (hence how close to `R0=−M0` an admissible return gets) a COMPUTED quantity** from the shared
generator — that is the upgrade.

- **EARNED (`derived_in_gate`, target-blind, able-to-fail):**
  - the **radiative-order fingerprints** `ω^(2ℓ+1)` for ℓ=0 (`ω¹`) and ℓ=1 (`ω³`), **derived** from the respective outgoing
    DtN `Λ_ℓ^out(z)=z·h_ℓ^{(1)′}(z)/h_ℓ^{(1)}(z)` (Hankel series), exactly as Gate 4 derived ℓ=2 from `Ŷ₂^out=−3/Λ₂^out`;
  - the **FORM, SIGN, and order** of the residuals `A0, A1` above;
  - the **non-over-cancellation**: that the shared generator forces `ε_eff≠0` (a bounded DC transfer), so the GR-departure
    survives — this is the predictive content. **⚠️ EARNED only if `ε_eff≠0` FOLLOWS from the derived generator/rank result
    (§2.1).** If non-zero/positivity is merely imposed as an admissibility *domain* (assumed, not derived), the nonzero-ness
    and magnitude are reclassified **DEFERRED**, and if the generator leaves it free the verdict is
    `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` — do NOT claim it as earned.
- **CALIBRATED / DEFERRED (not this gate's job):**
  - the **magnitudes** of `ε0_eff, ε1_eff` (in pathA_29 these are free slab-partition params; magnitude equality to pathA_29 is
    **NOT required** here — only form/sign/order/boundedness);
  - the overall gravity normalization (`G`, `54/5`) — inherited `GENUINE_BLOCKED` from Gate 4;
  - the full **nonlinear** closure (Gate 6) — here everything is linear-response about the stationary isotropic branch.

**Calibration boundary (unchanged):** "thorough, calibrations allowed, nothing asserted" — the gate delivers the *form/branch*
of the cross-ℓ relation; magnitudes stay labeled.

**⚠️ Scope precision (GLM Concern-3 — what "shared generator" is at the LINEAR level).** The brane↔bulk return law `S_return` is
part of the **stationary background** (the nonlinear throat solution), NOT the linear perturbation; the linearized PDE rides on
top of it. So Gate 5 does **not** solve/derive `S_return` — it **carries the return structure SYMBOLICALLY** (exactly as Gate 4
carried `D_n, N_n` symbolically). The claim under test is a **linear-algebraic consistency check**: *given a symbolic return
structure made consistent with ℓ=2, does the SAME structure projected onto ℓ=0/1 produce residuals of the pathA_29 class, with
`ε_eff` fixed (not free)?* It is NOT a derivation of the return law itself — that is the Gate-6 nonlinear closure. Any PASS claim
is therefore a *consistency* claim about a symbolic return, not proof that such a return is dynamically realized.

---

## §1. The object to build (requirement-level; Codex designs the route — `feedback-claude-reviews-codex-codes`)

1. **Identify the projection-blind generator `G`.** The generator is the **full linearized parent system** — the BdG matter
   sector (handoff §9.2) + the localized Maxwell sector (§9.3) + the geometry projection that carries the harmonic decomposition
   (§9.4: "projecting onto harmonics gives separate modal equations for `l=0`, grouped real `l=2`, and higher multipoles") —
   **plus its shared, ℓ-INDEPENDENT medium parameters** (`μ_η, T_w, K_η, T_Ω, c_s, …`). **⚠️ "Shared generator" means THIS parent
   PDE + shared medium parameters — NOT a shared port-kernel structure.** The already-reduced grouped-`P2` port symbols
   (`Ω_{U,A,r}, Ω_{W,A,r}, R_{A,r}, g_{U/W,A,r}`; `Δ_{A,r}=Ω_U²Ω_W²−R_{A,r}²`, `P_{A,r}=Ω_U²g_W+R g_U` — §10.2–10.3) are the
   **ℓ=2-specific** reduction of `G`; ℓ=0 and ℓ=1 are DIFFERENT reductions of the SAME `G` (see step 1′). The ℓ-dependence must
   enter ONLY via (i) the angular projection — the `T_Ω·ℓ(ℓ+1)` angular-stiffness term (`6T_Ω` at ℓ=2, `2T_Ω` at ℓ=1,
   **VANISHES at ℓ=0**), (ii) the reduction route (step 1′), (iii) the DtN `Λ_ℓ^out`, (iv) the source moments — **NOT** via
   independent per-ℓ generator parameters. **Emit the generator's parameter set and how each ℓ-sector's data is produced from it.**

1′. **The ℓ=0 sector is STRUCTURALLY different — handle it via the collective-variable reduction, not a port kernel (GLM
   Concern-1).** At ℓ=0 the angular-stiffness term `T_Ω·ℓ(ℓ+1)` vanishes (operator `μ_η ∂_t² − ∂_w(T_w ∂_w) + K_η`), and the
   axisymmetric ℓ=0 sector reduces to the **collective variables `(δa, δL)`** via the two-mode truncation of handoff §8.2 — which
   is exactly the reduction Gate 2 (`pathA_31`) already built — **NOT** the BdG/gauge/mixed-channel port kernel of §10. So the
   ℓ=0 return must be produced from the §8.2 / Gate-2 collective reduction (a genuine projection of the same `G`), while ℓ=1/ℓ=2
   ride the harmonic port-kernel-style reduction. **The §2.1 rank audit must explicitly check whether the ℓ=0 return IS
   determined by that collective reduction, or whether it would require a separate port-kernel-like reduction that does NOT exist
   in the handoff — the latter ⇒ `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` for the ℓ=0 channel.** (Do not silently force ℓ=0 into a
   port-kernel mold; that would turn Gate 5 into an ℓ=1-only test by construction.)
2. **Fix `G` + admissible branch on ℓ=2 and confirm Gate 4 still holds.** Re-derive (or import + re-verify) the ℓ=2 outgoing
   prefactor and confirm `QUAD_CALIBRATED` (fingerprint `1/9, 4/81, 1/27`, `P(ω)=D₀N/D^cons²`, `χ_Q=1`, `P0` scaling `a⁻⁵`).
   Any break ⇒ `FAIL_QUAD_REGRESSION`.
3. **Project `G` onto ℓ=0 and ℓ=1.** Derive the outgoing DtN `Λ_ℓ^out(z)=z·h_ℓ^{(1)′}(z)/h_ℓ^{(1)}(z)` from `h_0^{(1)},
   h_1^{(1)}` and the **static-slot-normalized response** `Y_ℓ^out = −(ℓ+1)/Λ_ℓ^out` (so `Y0=−1/Λ0`, `Y1=−2/Λ1`,
   `Y2=−3/Λ2` — the ℓ=2 case **must** reproduce Gate-4's `−3/Λ₂`, and the `−(1+1)=−2` is the origin of the ℓ=1 factor-of-½
   in `A1`; the generic `−(ℓ+1)` is verdict-bearing, not just `−3`). Derive the **return moments** `R0, R1` the generator
   produces; confirm the raw outgoing leading orders `ω¹` (ℓ=0) and `ω³` (ℓ=1). **Fix the time convention EXPLICITLY** (e.g.
   `e^{−iωt}`) and verify `h_ℓ^{(1)}` is the OUTGOING solution under it — the incoming choice flips the sign of the imaginary
   radiative term. In Gate 5 this sign trap routes to `FAIL_QUAD_REGRESSION` (if it breaks the ℓ=2 outgoing BC, §2.5/§3e) or
   `FAIL_EPSILON_MISMATCH` (if it flips the ℓ=0/1 residual sign, §2.3/§3c); tri-review verifies. (It is the analog of Gate-4's
   own `FAIL_NOT_OUTGOING`, which is NOT itself a Gate-5 rung.)
4. **Impose the side-conditions** `R0=−M0` (kills raw `O(ω¹)`) and `R1=−D1` (kills raw `O(ω³)`), with
   `M0(ω)=∫_brane S_leak d³x`, `D1_i(ω)=∫x_i S_leak + ∫j_i` (incl. carried odd wake) — `pathA_28_cancellation_condition.yaml`.
5. **Compute the residuals `A0, A1`** that survive when the **shared-generator** return is imposed (the DC transfer is
   `T_ℓ(0)=1/(1+ε_eff)`; **whether `ε_eff≠0` is forced — i.e. whether the return CANNOT be perfect — is the §2.1/§2.3 question,
   derived from the generator/rank result, NOT assumed here**: if the generator forces `ε_eff=0` ⇒ `FAIL_OVERCANCEL`; if it
   leaves `ε_eff` free ⇒ `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`), and compare to pathA_29's class
   (`reports/pathA_29_results.yaml`: `A0_residual`, `A1_residual`, `T_ell`, `T_at_DC`) in **form/sign/order/boundedness** — NOT magnitude.

---

## §2. Able-to-fail sub-checks (each independently computed)

**§2.1 CROSS-ℓ RANK / PROVENANCE AUDIT (the decisive able-to-fail core).** Prove the construction is NOT block-separable:
(a) **a DERIVED, SECTOR-SPECIFIC generator map.** The map `G → {P0, T0(0), T1(0)}` must be **derived from named PDE / projection
/ admissibility equations**, NOT asserted — and each sector via its OWN reduction (consistent with §1 step 1′, NOT all through
the ℓ=2 port kernel): **`P0` (ℓ=2)** via the §9.4 projection + the §10.2–10.3 grouped-`P2` port-kernel reduction; **`T1(0)`
(ℓ=1)** via the §9.4 harmonic (port-kernel-style) reduction at ℓ=1; **`T0(0)` (ℓ=0)** via the **§8.2 / Gate-2 collective
`(δa,δL)` reduction** (the ℓ=0 sector has no `T_Ω` angular stiffness and no §10 port kernel — see §1 step 1′). **If any sector's
map cannot be exhibited non-asserted** — in particular if the §8.2 collective reduction cannot determine `T0(0)`, or no ℓ=0/1
reduction of `G` exists beyond the postulated flat slab — **the verdict is `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (or
`FAIL_TAUTOLOGICAL` if a map is asserted), NOT a deferred or PASS result.** This is the honest landing if the generator cannot be
reconstructed for a sector.
(b) **provenance** — trace every parameter feeding the ℓ=0/1 return back to the generator `G` fixed on ℓ=2; assert **no new
`ε0`, `ε1`, or ℓ-gain parameter** is introduced after ℓ=2 is fixed (any introduced-after parameter ⇒ `FAIL_DECOUPLED`).
(c) **rank** — after applying the ℓ=2 `P0` constraint to `G`, **DIRECTLY COMPUTE the NATIVE residual freedom (the null space)**
in the ℓ=0/1 DC return `T0(0), T1(0)` from `G` itself — do **not** rely on the §3b injection probe to reveal it (the probe tests
the *detection mechanism*; the native null space must be computed independently). If a **native null direction** moves `T_ℓ(0)`
at constant `P0=N0/D0`, it is acceptable **only if a branch selector DERIVED from named PDE/admissibility equations kills it** —
with a self-ablation: *selector derived ⇒ null killed; selector removed ⇒ `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`*. An
asserted/hand-imposed selector does NOT count (⇒ `FAIL_TAUTOLOGICAL`). The PASS condition is that the computed native null space
is empty (or fully killed by a derived selector), so `T0(0), T1(0)` (hence `ε0_eff, ε1_eff`) are **fixed functions of the
ℓ=2-constrained generator** (a genuine prediction), with the null-space computation + the selector's provenance emitted.

> **⚠️ IMPLEMENTATION FIREWALL for §2.1 (added after the v1 REJECTION — the rank audit was found rigged-to-`UNDERDETERMINED`).**
> The rank audit must operate on the **actual constraints collected from the named equations** (the §8.2/Gate-2 ℓ=0 collective
> reduction, the §9.4 ℓ=1 harmonic reduction, the §10 ℓ=2 port kernel — with ℓ=2 `N0` the DERIVED `n0_port=P²/Δ²`, NOT a free
> symbol), and it must be **able to return BOTH "determined" and "underdetermined"** depending on whether those constraints pin
> the return. A coordinate/constraint setup in which `UNDERDETERMINED` is structurally guaranteed (e.g. carrying `Z0_ret,Z1_ret`
> as free coordinates against a constraint row that is zero in their columns by construction) is **`FAIL_TAUTOLOGICAL`** — the
> verdict must be EARNED, and the freedom in `Z` must EMERGE as "no collected named constraint touches it + pathA_29 establishes
> `Z` is a premise," not be assumed. The `derived_selector` must be a **real symbolic equation that is substituted and the null
> space re-computed** (so the verdict genuinely flips), NOT a boolean flag. **§0's "carry the return structure symbolically"
> permission applies to the return FORMS, NOT to the verdict-determination mechanism** — every verdict-bearing control (rank audit
> + each §3 probe) must compute from genuinely-mutated quantities and be able-to-fail, dual-engine.

**§2.2 ℓ=0/1 FINGERPRINT DERIVATION (earned).** Derive the response `Y_ℓ^out = −(ℓ+1)/Λ_ℓ^out` (with
`Λ_ℓ^out=z·h_ℓ^{(1)′}/h_ℓ^{(1)}`) and its Hankel low-`z` series; assert (a) the normalization factor is `−(ℓ+1)` — `Y0=−1/Λ0`,
`Y1=−2/Λ1`, `Y2=−3/Λ2` (the ℓ=2 value **must** match Gate 4); (b) the raw outgoing leading orders are `ω^(2ℓ+1)` (ℓ=0→`ω¹`,
ℓ=1→`ω³`); (c) the imaginary/radiative part **emerges from the outgoing BC under the fixed time convention** (not an imposed
dissipation; the incoming convention flips its sign — cf. Gate-4 `FAIL_NOT_OUTGOING`). A wrong `−(ℓ+1)` factor, order, or
coefficient ⇒ `FAIL_EPSILON_MISMATCH` (order/coefficient branch) or `FAIL_QUAD_REGRESSION` (if the ℓ=2 fingerprint regresses).

**§2.3 RESIDUAL MATCH TO pathA_29 (form/sign/order/boundedness — NOT magnitude).** Assert the computed `A0, A1` reduce, at
leading order, to `i·a·(ω/c_s)·M0·ε0_eff/(1+ε0_eff)` and `i·a³·(ω/c_s)³·D1·ε1_eff/(2(1+ε1_eff))` with `ε_eff` finite, `>0`,
`T_ℓ(0)∈(0,1)` (bounded). Wrong sign/order/unbounded ⇒ `FAIL_EPSILON_MISMATCH`; `ε_eff→0` (clean cancellation) ⇒
`FAIL_OVERCANCEL`; no admissible branch reaches the class ⇒ `FAIL_NO_CONSISTENT_RETURN`.

**§2.4 EXISTENCE OF AN ADMISSIBLE BRANCH.** Confirm at least one admissible branch (normalizable transverse zero mode /
localizing family, per pathA_29's `p=2` survival) simultaneously satisfies ℓ=2 `P0` and the ℓ=0/1 residual class; if provably
none ⇒ `FAIL_NO_CONSISTENT_RETURN`.

**§2.5 GATE-4 NON-REGRESSION.** Re-run the Gate-4 headline on the shared generator; any deviation in fingerprint / prefactor
object / `χ_Q` / outgoing BC ⇒ `FAIL_QUAD_REGRESSION`.

**§2.6 DIMENSIONAL (μ̂₀-FREE, able-to-fail, dual-engine — genuine from the start).** Walk `dim_of`/`dimOf` over the **assembled
residual/transfer expressions** (`A0, A1, T_ℓ, ε_eff` and the generator data). **Declare a source-of-truth dimension table** for
`{M0, D1, R0, R1, A0, A1, T_ℓ, ε_ℓ, N0, D0}` plus the base symbols (`a, c_s, ħ, m, K, ρ0, …`), each sourced INDEPENDENTLY from
the action / EOS `P=Kρ⁵` / the moment definitions (`M0=∫S_leak d³x`, `D1=∫x_i S_leak+∫j_i`) — **NOT inferred from the residual
expressions** (back-solving any of these from `A0/A1/T_ℓ` is forbidden; it re-creates the free-carrier tautology). The
verdict-bearing booleans must be **free-carrier-independent**: `T_ℓ` and `ε_ℓ` must come out **dimensionless**, and `A0, A1`
homogeneous with their independently-sourced amplitude dimensions. Able-to-fail: corrupt a **sourced** input dimension →
`FAIL_DIMENSIONAL` (§3f); corrupting an unconstrained carrier must leave it unchanged.

**§2.7 PROVENANCE PARTITION (earned vs calibrated/deferred).** Classify each emitted quantity into
{`derived_in_gate` (the `ω^(2ℓ+1)` fingerprints, the `−(ℓ+1)` normalization, the residual form/sign/order, **and the forced
`ε_eff≠0` ONLY IF it follows from the §2.1 generator/rank result**), `external_bridge_input` (the `R0=−M0`/`R1=−D1`
side-conditions from pathA_28; the `M0/D1` moments), `deferred_branch_data`, `convention`}. **For the `ε_eff` magnitude, keep
three cases DISTINCT (GLM NIT-5), do not conflate:** (i) **computed from the generator but NOT cross-checked to pathA_29's
magnitude** ⇒ `derived_in_gate` *flagged `magnitude-not-cross-checked`* (the nonzero-ness/value IS derived; only the literal
pathA_29 comparison is deferred); (ii) **not computed because the generator is underdetermined** ⇒ `deferred_branch_data` (and
the verdict is `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`); (iii) the literal pathA_29 `ε0/ε1` partition magnitudes ⇒
`deferred_branch_data`. COMPUTED by provenance trace, not asserted.

---

## §3. Mandatory counterfactual probes (baked into the script; each a real two-verdict `self_ablation`)

⚠️ **`self_ablation` must be an ACTUAL RE-RUN** (the Gate-4 adversarial-audit defect to avoid): each probe re-executes the gate
logic with the mutation REMOVED and records the verdict (`with_mutation` = the FAIL, `without_mutation` = NOT the FAIL). A probe
that leaves the verdict unchanged ⇒ `FAIL_TAUTOLOGICAL`.

- **§3a DECOUPLE-KNOBS.** Inject an independent harmonic-projector scale on ℓ=0/1 (a parameter not present in `G`); assert
  §2.1 fires `FAIL_DECOUPLED` (ℓ=2 still hits `P0` while `ε0` floats) — and that WITHOUT the injected knob it does not.
- **§3b NULL-DIRECTION (TWO ablation states — GLM Concern-4).** **Inject** a *shared-generator* perturbation that preserves
  `N0/D0` (constant `P0`) while moving `T0(0)/T1(0)` — a direction in the generator's `P0`-null space, **distinct from a
  decoupled per-ℓ knob (§3a)**. Report BOTH ablation states: **(a)** remove the injection (tests that §2.1(c) DETECTS an
  artificial null direction → fires `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` with injection, reverts without); **(b)** the
  **baseline run with NO injection** — this is where §2.1(c)'s **directly-computed native null space** is read: if the native
  null space is non-empty and no derived selector kills it, the verdict is `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`
  **regardless of the probe**. The probe must NOT mask a native null direction — state (b) is the load-bearing one; the
  injection (a) only certifies the detector works even when the native case is clean.
- **§3c WRONG-SIGN / ANTI-LOCALIZING RETURN.** Flip the DC transfer sign / use a non-normalizable (delocalizing) return;
  assert §2.3 fires `FAIL_EPSILON_MISMATCH`.
- **§3d PERFECT-RETURN (over-cancel).** Force `ε_eff→0` (`T_ℓ(0)=1`); assert §2.3 fires `FAIL_OVERCANCEL` and that the
  physical bounded-residual branch does NOT.
- **§3e BREAK-GATE-4.** Swap the ℓ=2 outgoing BC for incoming/standing, or the prefactor object `D₀N/D^cons²`→`N/D`; assert
  §2.5 fires `FAIL_QUAD_REGRESSION`.
- **§3f CORRUPT-DIMENSION.** Corrupt a **sourced** input dimension (e.g. `[M0]→[M0]·L`); assert §2.6 fires `FAIL_DIMENSIONAL`;
  corrupting an unconstrained carrier leaves it unchanged (free-carrier-independence proof).
- **§3g ASSERT-NOT-DERIVE.** Replace a derived fingerprint/`ε_eff` with a hardcoded literal; assert §4 fires
  `FAIL_TAUTOLOGICAL`.
- **§3h NO-CONSISTENT-RETURN.** Mutate the admissible-branch family so that hitting ℓ=2 `P0` is provably incompatible with a
  bounded positive `T0(0)/T1(0)` (e.g. restrict to a branch whose `P0`-matching point forces `ε_eff<0`/unbounded); assert §2.4
  fires `FAIL_NO_CONSISTENT_RETURN`, and that the unrestricted family does not. (This is the probe BLOCKER-1 requires so that
  rung is demonstrably reachable.)

---

## §4. Anti-tautology firewall + the shared-generator certificate

The script must: (1) DERIVE the ℓ=0/1 fingerprints from the named DtN objects `Λ_ℓ^out`, not assert them; (2) emit the
**shared-generator certificate** — the generator parameter set + the provenance trace showing ℓ=0/1 carry NO parameter not fixed
by the ℓ=2 constraint + the rank check (no null direction); (3) carry the §3 probes as computed mutations with real
`self_ablation`; (4) keep ℓ=2 `QUAD_CALIBRATED` as a live re-check, not an import-on-faith. If the shared-generator claim, the
cross-ℓ relation, or any fingerprint is asserted rather than derived ⇒ `FAIL_TAUTOLOGICAL`. **The DEFAULT verdict is the PASS;
the script must be able to reach every PHYSICS/LOGIC FAIL** — the §3 probes (§3a–§3h) demonstrate reachability for
`FAIL_DECOUPLED`, `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`, `FAIL_EPSILON_MISMATCH`, `FAIL_OVERCANCEL`, `FAIL_QUAD_REGRESSION`,
`FAIL_DIMENSIONAL`, `FAIL_TAUTOLOGICAL`, and `FAIL_NO_CONSISTENT_RETURN`. The **infrastructure rungs** `FAIL_ENGINE_DISAGREE`
and `ENV_BLOCKED` fire from the cross-engine comparator / environment, not a §3 mutation, and are explicitly **excluded** from
the "probe-reachable" requirement.

---

## §5. Dimensional consistency (standing pre-numbers step — genuine, dual-engine)

Units-restored homogeneity on the REAL assembled expressions via the `dim_of`/`dimOf` walker (the Gate-4 pattern;
`pathA_33_*` + `src/stage1_solver/dimensional_check.py`), free-carrier-independent, with an able-to-fail perturbation (§3f),
run by BOTH engines. Emit a full homogeneity table. (Do NOT reinvent the Gate-1–3 vacuous typed-tuple ledger.)

**⚠️ DIMENSIONAL FIREWALL (standing, user-mandated — no exceptions).** **No LLM — not GLM, not Codex review prose, not this
directive's hand-worked examples — establishes a dimensional fact.** Every `[·]` exponent in the source-of-truth table (§2.6) is
declared from the action/EOS/moment definitions, and every dimensional VERDICT (`dimensional_ok`, `FAIL_DIMENSIONAL`) is
**COMPUTED in the Codex-written test code** by the dual-engine `dim_of`/`dimOf` walker on the real assembled expressions, with the
§3f able-to-fail probe firing on a corrupted *sourced* input. The tri-review's dim-check item (§8) verifies the code computes it
(not that any prose asserted it). Any dimensional claim that exists only in review prose is treated as **unverified** until the
dual-engine walker reproduces it. (LLM Hankel/sign/order checks in reviews are useful leads but are re-derived in-code too.)

---

## §6. Deliverables

- `tools/pathA_34_cross_l_unification_{sympy.py,.wl}` — dual engine; both independently assemble the HEADLINE (the shared
  generator → ℓ=0/1 outgoing fingerprints + residuals `A0,A1` + `ε_eff` + the ℓ=2 non-regression check), not just scalar
  pass/fail booleans; cross-engine comparator on the residual/transfer expressions (symbolic + numeric tolerance →
  `FAIL_ENGINE_DISAGREE` if exceeded).
- `reports/pathA_34_results.yaml` — the computed verdict (one ladder rung), the shared-generator certificate (parameter set +
  provenance trace + rank check), `Λ_0^out/Λ_1^out` + fingerprints, `A0/A1/T_ℓ/ε0_eff/ε1_eff`, the pathA_29 form/sign/order
  comparison, the provenance partition (§2.7), the dim-homogeneity table, and the §3 probes each with two real verdicts.
- `reports/pathA_34_cross_l_unification.md` — prose writeup (verdict + what is EARNED vs CALIBRATED/DEFERRED + the honest gap).

---

## §7. Deferred / out of scope (record, do not silently drop)

- **ε0/ε1 magnitude** — uncalibrated; form/sign/order only. **Literal pathA_29 magnitude reconciliation** — deferred (pathA_29
  is a postulated flat slab; the ledger's ℓ=2 data are the mixed-sector port kernel — the reconciliation note's honest gap).
- **Full nonlinear closure (Gate 6)** — everything here is linear-response about the stationary isotropic branch.
- **If the verdict is `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (GLM Concern-2 — the most likely outcome):** that is a **structured,
  informative result**, not merely an open item — it localizes *where* the linear level cannot fix the ℓ=0/1 return partition
  (no BCs without the nonlinear solution). The cross-ℓ *consistency* test then **defers to Gate 6** (where the nonlinear solve
  determines the return law and unifies the two reductions — the reconciliation note §3 calls this unification "the nonlinear
  closure"). Gate 5's contribution in that case = turning "open item #9" into a *characterized* underdetermination (which channel,
  which null direction, what a Gate-6 selector must supply), not a guess.
- **PN match-back** — the decisive downstream falsifier (`research/4d_*pn*`); do NOT re-derive the PN ladder.
- **BC provenance** — inherits Gate 1's `BC_DEPENDENT` (labeled calibration input).

---

## §8. Process (the standing gauntlet — run the review rounds via per-gate AGENTS, `feedback-offload-review-gauntlet`)

1. **Codex design-review (xhigh)** of THIS directive → fold corrections.
2. **GLM** single fresh-perspective pass (user-run, out-of-band) → fold → **Codex confirm** to GREEN.
3. **Dual-engine execute** (Codex `--sandbox danger-full-access` for the `.wl`; `codex exec … -c model_reasoning_effort=xhigh`,
   backgrounded, never `timeout`-wrapped; scripts `timeout 600`).
4. **Tri-review** (clean agents): orchestrator arbiter re-run (both engines) + transliteration-fidelity + **adversarial-with-
   ablation**. **Explicit tri-review items:** (a) the second engine independently assembles the HEADLINE (no `x−x`); (b) every
   verdict-bearing control computes-from-inputs + has an able-to-fail probe with a REAL two-verdict `self_ablation`; (c) the
   dim check is free-carrier-INDEPENDENT (corrupting a sourced input fires `FAIL_DIMENSIONAL`); (d) the **shared-generator
   certificate** is genuine — the ℓ=0/1 prediction carries no parameter not fixed by ℓ=2 (the adversarial leg must itself try
   the §3a decouple-knob and confirm the verdict flips).
5. **User gate**, then commit (only when asked; stage explicit paths). On completion: flip task #99, update `STATUS.md` ⭐⭐
   + `decisions/13` §0 + the completion ladder Gate-5 entry.

---

## Cross-refs

- Teeth discussion: `software/stage1_solver/_scratch/pathA_34_gate5_teeth_design_discussion.md` +
  `software/stage1_solver/_scratch/pathA_34_gate5_teeth_codex.log` (CONVERGED); design-review:
  `software/stage1_solver/_scratch/pathA_34_design_review_codex.log` (SOUND-WITH-CONCERNS, 8 issues folded).
- Framing: `research/pde_ledger/notes/stages/moving_throat_pde_open_item_reconciliation.md` (§2 mapping, §3 honest gap, §4.3
  cross-ℓ gate).
- Targets: `software/stage1_solver/reports/pathA_28_cancellation_condition.yaml` (`R0=−M0`, `R1=−D1`, `A0/A1` orders);
  `software/stage1_solver/reports/pathA_29_results.yaml` (`A0/A1_residual`, `T_ell`, `T_at_DC`, `ε0/ε1`).
- Machinery / template: `software/stage1_solver/directives/pathA_33_quadrupole_normalization.md` +
  `tools/pathA_33_quadrupole_normalization_{sympy.py,.wl}`; handoff §9.4, §10.2–10.3, §12.
- Memories: `project-moving-throat-pde-push`, `feedback-falsification-is-the-goal`, `feedback-negative-verdict-short-circuit`,
  `feedback-dimensional-consistency-check`, `feedback-decisive-test-not-tautological`, `feedback-claude-reviews-codex-codes`.
