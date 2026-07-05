# Directive pathA_42 — the charge-coupled scalar: does it BREAK the model or is it naturally hidden? (map the departure; magnitude SIM-GATED)

**Status:** DRAFT v3.1 (Codex design-review + confirm-GREEN + GLM-5.2 tertiary `SOUND_WITH_NITS` + Codex final ADJUDICATION all
folded — see §9. State machine = GREEN (closed). Radiation-normalization physics ADJUDICATED against pathA_38: `K_h` is NOT pinned by
the static-Coulomb calibration (that goes through `Q_E`+`N₀`, not `K_h`) ⇒ current-ledger branch = exponent-3, radiation `SIM_GATED`
on the free `M_h` even at `c_E=c_γ`; the GLM-F6 `FALSIFIABLE_TENSION` is conditional on a future `K_h`-pinning. Final Codex re-check =
directive ALIGNED/GREEN (the only stale-phrase grep hit is the archived GLM review transcript in `_scratch`, not the directive). ▶
READY FOR USER GATE (execution).)
**Prior:** pathA_39 Stage 4 = `FIELD_SCALAR_VECTOR_DEPARTURE` (EM field content = 2 transverse photons + `u_L` density + `h`-branon;
4 physical DOF; deleting the `h` block → exact Maxwell). pathA_40 cone-lock = `CONE_LOCK_CALIBRATED` (`c_E=c_γ` calibrated, not
derived; `Δr=2`). pathA_41 NG5 = `SECOND_MEDIUM_DRIFT{ρ_B0,χ_c,C_hu}` (there: `c_E`→`SIM_DEFERRED(Route-A)`, `M_h`→`CALIBRATED_GEOMETRY_INPUT`).
**Scope chosen by user (2026-07-05):** the "focused gate — earn the map" option, after "attempt to pin `M_h`/`c_E`" hit a HARD WALL.

**This gate answers task #110:** the magnetism sector landed as a *characterized departure* — a propagating charge-coupled scalar
inseparable from the EM field. Does it show up in a way that would CLEARLY BREAK the model, or is it naturally hidden? This gate
EARNS the falsifiable **map** of where the break-risk lives; it does NOT (cannot) pin the decisive magnitude.

---

## 0. What this gate is + honesty preregistration (READ FIRST)

**The decisive magnitude is SIM-GATED and this gate does NOT try to pin it.** Reachability read = `HARD_WALL`: `M_h` (branon
inertia; pathA_41 `CALIBRATED_GEOMETRY_INPUT`) and `c_E` (branon speed, `K_h=M_h c_E²`; pathA_41 `SIM_DEFERRED(Route-A)`) cannot be
projected from earned objects — `c_E` is hand-inserted into the pathA_38 dynamic radial operator (`pathA_38_results.yaml:306`) and
`M_h` is symbolic/imported (`pathA_39_scalar_admixture_results.yaml:64`). No earned quadratic parent action carries the `h`-sector
TIME-KINETIC weight (only the in-plane NORMALIZATION weight `N₀=8/(3ℓ)` exists). Pinning needs the sim-deferred nonlinear throat
solve (pathA_40 `ROUTE_A_MISSING_NONLINEAR_THROAT`). **Therefore any numeric magnitude is honestly `SIM_GATED`. Do NOT launder one.**

**⭐ THE #1 RIG RISK (Codex — guarded MECHANICALLY, not by prose):** "assuming the missing kinetic metric by notation" — reading the
zero-mode NORMALIZATION `N₀`/`K_parallel` as the dynamical INERTIA `M_h`, or the Green-inserted `c_E` as a DERIVATION of `c_E`, then
"deriving" `c_E=c_γ` or a radiation number by convention. **HARD OUTPUT INVARIANT (§4.A):** no numeric `M_h`, `c_E`, `K_h`,
`P_h/P_EM`, or EP magnitude may be emitted UNLESS an earned `h_time_kinetic_parent_action` provenance object exists in the ledger
(it does NOT today). Negative fixtures inject `M_h:=N₀`, `M_h:=K_parallel`, `c_E:=c_γ from Green form` → the pipeline MUST emit
`LAUNDERING_REFUSED`, not a number and not a bare prose `SIM_GATED`.

**What the gate DOES earn (all able-to-fail):** (a) EP-safety of the `h`-branon ON THE DECOUPLED FLOOR (charge-only ⇒ Coulomb-
degenerate ⇒ not Gate-L's `FAIL_BENDING_MASSLESS_FIFTH_FORCE`) — subject to the universality CONJUNCTION (safety additionally needs
`q_h/Q_E` universal) and the honest caveat that `C_hu` mixing REINTRODUCES a mass coupling; (b) the radiation SCALING from an
explicit dipole/flux derivation — exponent-3 on the CURRENT ledger (pathA_38 does not pin `K_h`; radiation `SIM_GATED` on the free
`M_h` even at `c_E=c_γ`), exponent-1 CONDITIONAL on a future `K_h`-pinning (SCALING earned, NUMBER `SIM_GATED`); (c) the
preferred-frame cost of the hiding limit; (d) the `q_h/Q_E` universality constraint; (e) the `u_L`+mixing
mass-coupled EP structure (magnitude `SIM_GATED` via `a_L` AND `C_hu`); (f) the Gate-L bending-mode connection. **A clean "it's all
fine" is suspicious** ([[feedback-falsification-is-the-goal]]): each channel must be able to land a break.

---

## 1. Substrate (imported; do NOT re-derive — cite the earned reports)

Linearized EM sector (from `pathA_39_scalar_admixture_results.yaml:82-91`), 2-field scalar block `Φ=(u_L, h)` + transverse `u_T`:
- **transverse photon:** `½ρ_br(∂_t u_T)² − ½μ_R(curl u_T)²`, 2 pols at `c_γ²=μ_R/ρ_br` (`pathA_39_magnetic_force.md:92`).
- **`h`-branon:** `½M_h(∂_t h)² − ½K_h k² h²`, `K_h=M_h c_E²`, NO mass term ⇒ gapless, speed `c_E`. Charge coupling `q_h·h`,
  `q_h=2Q_E tanh(b/ℓ)/b` (COMPUTED, pathA_38). **`h` bare source has mass residue 0** (`source_mass_vector=[q_M,0]`,
  `pathA_39_scalar_admixture_results.yaml:90-91`).
- **`u_L` density:** `½ρ_br(∂_t u_L)² − ½B_eff k² u_L²`, `B_eff=ρ_B0²/χ_c`, `c_s²=B_eff/ρ_br`. Couples to charge `q_L=Nu·a_L·sCharge`
  AND to mass `q_M` ⇒ the EP channel; `a_L` SIM-DEFERRED.
- **mixing:** `−C_hu k² u_L h`; the static 2×2 exchange matrix carries a charge-mass amplitude `A_qm ⊃ C_hu·q_h·q_M`
  (`pathA_39_scalar_admixture_results.yaml:107`); mixed eigenpoles carry mass residues (`:121-133`). Off-cone residual
  `det M|cone = −C_hu²k⁴` (pathA_40). `C_hu` SIM-DEFERRED/free (stability `C_hu²<B_eff K_h`).
- **`h` charge residue floor:** `4Q_E²tanh²(b/ℓ)/(M_h b²)`, `a_L`- AND `C_hu`-INDEPENDENT (`pathA_39_scalar_admixture_screen.md:24`).
- **preferred-frame:** pathA_39 records `PREFERRED_FRAME_UNLESS_cE_EQUALS_cGAMMA`, residual `(c_E²ρ_br−μ_R)/ρ_br`
  (`pathA_39_magnetic_force.md:63-65`); pathA_40 leaves `c_E=c_γ` calibrated not derived.
- **Gate-L failure mode:** pathA_35 required the embedding mode `u_w` GAPPED (`Ω_w`) to avoid `FAIL_BENDING_MASSLESS_FIFTH_FORCE`
  (`pathA_35_gateL_light.md:132`). The pathA_38 Coulomb `h` is the **same embedding-direction family, a DISTINCT ledger object**
  (gapless); no earned action bridges gapped-`u_w` → gapless-`h` (that bridge = the deferred throat program).

---

## 2. The decisive computations (five channels; each able-to-fail)

### 2.1 `h`-branon EP-safety — via the FULL static 2×2 exchange matrix (Codex BLOCKER 1)
Compute the static scalar exchange from the FULL `2×2` matrix over `(u_L,h)`: `A_qq` (charge-charge), `A_qm` (charge-mass), `A_mm`
(mass-mass), NOT just the decoupled `h` pole.
- **Decoupled floor (`C_hu=0`):** `A_qm`'s `C_hu` term vanishes ⇒ `h` couples to charge only ⇒ its static 1/r² is degenerate with
  Coulomb ⇒ absorbable into the calibrated `Q_E` ⇒ NOT a mass-coupled (composition-dependent) EP fifth force. Earn `h_EP =
  EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR`.
  **⚠ CONJUNCTION (GLM F1):** the "absorbable into `Q_E`" step is valid ONLY if `q_h/Q_E = 2tanh(b/ℓ)/b` is UNIVERSAL across species.
  If `b/ℓ` is species-indexed (§2.3 leaves this open), then `q_h ∝̸ Q_E`, and for neutral bulk matter (`ΣQ_E=0` but `Σq_h≠0`) the `h`
  scalar gives an unscreened composition-dependent long-range force = an EP fifth force ON the decoupled floor. So this token covers
  ONLY the mass-coupled channel; FULL decoupled-floor EP safety is the CONJUNCTION `EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR ∧
  universality=EARNED_SAFE`. The build MUST NOT report unqualified "decoupled-floor EP-safe" while `universality≠EARNED_SAFE`.
- **With mixing (`C_hu≠0`):** `A_qm ⊃ C_hu·q_h·q_M` ⇒ the `h`-like eigenmode INHERITS a mass coupling ⇒ a composition-dependent
  force reappears. Earn this REINTRODUCTION; its magnitude is `SIM_GATED via C_hu` (and `q_M`).
- **Able-to-fail controls (§4):** (i) inject nonzero `h` BARE mass residue → `FIFTH_FORCE_TRIGGERED` even at `C_hu=0` (proves the
  decoupled safety is earned from residue-0, not asserted); (ii) fixture `q_L=0, C_hu≠0, q_M≠0` → MUST flag `MIXED_SCALAR_EP_RISK`
  (proves the mixing caveat is computed, not decorative).

### 2.2 Radiation scaling — EXPLICIT dipole/flux derivation, under the RIGHT normalization (Codex BLOCKER 4 + GLM F2/F6)
Do NOT assert the exponent. DERIVE it in BOTH engines from a stated setup: (i) a nonrelativistic oscillating point-charge source;
(ii) the far-zone retarded solution of the `h` field (speed `c_E`) vs the transverse EM field (speed `c_γ`); (iii) the scalar
energy-flux / stress tensor vs the EM Poynting comparator.
- **⭐ THE NORMALIZATION IS DECISIVE (GLM F2/F6, ADJUDICATED against pathA_38 — Codex):** the exponent DEPENDS on what is held fixed.
  (a) **Bare-fixed** (`q_h`, `M_h` independent): `P_h/P_EM ∝ (c_γ/c_E)³`. (b) **Static-calibrated** (IF the static coupling `q_h²/K_h`
  were pinned to the measured Coulomb, fixing `K_h=M_h c_E²`): `P_h/P_EM ∝ (c_γ/c_E)¹` (order-one at `c_E=c_γ`). **ADJUDICATION (Codex,
  pathA_38 evidence):** the pathA_38 static Coulomb is `U_int ∝ q_eff²/N₀ · 1/(4πR)` — calibrated through `Q_E` and the zero-mode
  NORMALIZATION `N₀=8/(3ℓ)`, NOT through the dynamic stiffness `K_h`. So `K_h`/`M_h` are NOT pinned by the static calibration ⇒
  **the CURRENT-LEDGER physical branch is exponent-3, and the radiation magnitude stays GATED on the free `M_h`** (via the residue
  floor `q_h²/M_h = 4Q_E²tanh²(b/ℓ)/(M_h b²)`). The exponent-1 branch is CONDITIONAL — reachable only if a FUTURE ledger fact pins
  `K_h` independently of the `Q_E` source calibration. The build must DERIVE this pinned-`K_h` test (default NO on the current ledger)
  and emit BOTH exponents + which branch is active + the pathA_38 evidence. Do NOT call exponent-1 "the physical case" for the
  current ledger.
- **Able-to-fail:** ablations that corrupt the Green/flux speed MUST change the derived exponent; and a fixture applying the WRONG
  normalization (absorption vs bare-fixed) MUST produce the OTHER exponent (proves the two are distinguished in code, not conflated —
  GLM F2). If the derivation cannot be mechanized in both engines, set `radiation=AUDIT_ONLY_NOT_EARNED` (§3).
- **The two ends (current ledger, exponent-3, `M_h` free):** even at `c_E=c_γ` the magnitude is GATED on the free `M_h` ⇒
  `radiation=SIM_GATED` (NOT order-one-earned ⇒ NOT `FALSIFIABLE_TENSION` on the current ledger — Codex overturns the GLM-F6
  unconditional claim); `c_E≫c_γ` ⇒ suppressed but preferred-frame (§2.5). Emit NO numeric `P_h/P_EM`.
  **`radiation=FALSIFIABLE_TENSION` is reachable ONLY CONDITIONALLY** — if a future ledger pins `K_h` (exponent-1 branch) AND the
  model commits to `c_E=c_γ` (§3 + §4.K, both gated on the pinned-`K_h` fact).
- **Vacuum-Cherenkov deferral (GLM F3):** if the eventual pinned `c_E<c_γ`, a threshold Cherenkov channel opens (charges with
  `c_E<v<c_γ` radiate `h`), distinct from dipole radiation — carry subtag `CHERENKOV_DEFERRED` on `radiation=SIM_GATED` as a
  throat-solve re-open obligation (NOT a sixth channel).

### 2.3 Universality constraint (Codex Finding 1 — no free EARNED_SAFE)
`q_h/Q_E = 2tanh(b/ℓ)/b`. From the ledger, is `b/ℓ` a GLOBAL constant or a PER-SPECIES datum? pathA_41 records `b`,`ℓ` as
CALIBRATED_GEOMETRY inputs, NOT a universality theorem. **Production default:** `universality = SIM_GATED_REQUIRED_UNIVERSALITY` (a
condition the throat program MUST satisfy) UNLESS an earned global source schema is cited; if the ledger provides a species-indexed
source, `universality=FALSIFIABLE_TENSION` (§3 rule 2). It may NOT become `EARNED_SAFE` merely by stating a required condition.
**Able-to-fail:** a species-varying `b/ℓ` fixture → non-universality flag fires.

### 2.4 `u_L` + mixing mass-coupled EP channel — structure earned, magnitude SIM-GATED via `a_L` AND `C_hu` (Codex BLOCKER 1)
The `u_L` mode couples to charge (`q_L`) AND mass (`q_M`) ⇒ the genuine EP channel; and `C_hu` mixing transfers mass residue into the
charge-sourced eigenmodes (§2.1). Earn the STRUCTURE (couples to mass; sign; magnitude `∝ a_L²` from `q_L` PLUS the `C_hu·q_h·q_M`
mixing term). **Report the magnitude `SIM_GATED via a_L AND C_hu`; do NOT assert "gravity-weak"** (no earned argument it is
gravitational-strength — a HYPOTHESIS only).

### 2.5 Preferred-frame channel (Codex BLOCKER 3 — new)
The radiation hiding limit (`c_E→∞`) is NOT free: `c_E≠c_γ` ⇒ `PREFERRED_FRAME_UNLESS_cE_EQUALS_cGAMMA` (pathA_39), residual
`(c_E²ρ_br−μ_R)/ρ_br`. Emit `preferred_frame` ∈ {`EARNED_SAFE` iff an earned `c_E=c_γ` exists (it does NOT — pathA_40 calibrated),
`PREFERRED_FRAME_TENSION` iff suppression requires `c_E/c_γ≫1`, `SIM_GATED` otherwise}. **Binding coupling to §2.2/§3:**
`NATURALLY_HIDDEN` is IMPOSSIBLE if the only radiation-suppression path requires `c_E/c_γ≫1` (that path pays preferred-frame
tension). At the calibrated `c_E=c_γ`, radiation is `SIM_GATED` on the current ledger (gated on the free `M_h` — §2.2 adjudication);
it becomes order-one/`FALSIFIABLE_TENSION` ONLY under a future pinned-`K_h` fact (exponent-1 branch).

---

## 3. Verdict grammar (first-match; a CLOSED per-channel state machine — Codex Finding 2 + confirm-pass)
**Canonical per-channel state tokens (each channel takes EXACTLY one; no other token is legal — a build emitting an out-of-set
token is itself a `NO_GO_CONSISTENCY`):**
- `h_EP` ∈ { `EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR` (charge-only at `C_hu=0`, residue-0 holds — covers ONLY the mass-coupled
  Gate-L `FAIL_BENDING_MASSLESS_FIFTH_FORCE` channel; FULL decoupled-floor EP safety ADDITIONALLY requires `universality=EARNED_SAFE`,
  because non-universal `q_h/Q_E` gives an unscreened composition-dependent force even at residue-0 — see §2.1 + GLM F1) |
  `FIFTH_FORCE_TRIGGERED` (nonzero bare mass residue ⇒ composition-dependent force even at `C_hu=0`; a BREAK) | `NO_GO` (residue-0
  earned fact contradicted, or static-Coulomb mismatch) }.
- `radiation` ∈ { `SIM_GATED` (exponent DERIVED, magnitude gated on `c_E` AND on the free `M_h` — the CURRENT-LEDGER state, since
  pathA_38 does not pin `K_h`; gated even at `c_E=c_γ`) | `FALSIFIABLE_TENSION` (order-one with no preferred-frame-free hiding —
  reachable ONLY CONDITIONALLY: a future ledger pins `K_h` (exponent-1 branch) AND the model commits to `c_E=c_γ`; NOT the current
  ledger, Codex overturned the unconditional GLM-F6 claim) | `AUDIT_ONLY_NOT_EARNED` (the §2.2 dipole/flux derivation could not be mechanized in
  both engines ⇒ the scaling is an audit assumption, may NOT support an earned claim; rides on the mapped verdict as a flag, is NOT
  a break) | `EARNED_SAFE` (magnitude PINNED and shown negligible vs bounds — requires the deferred throat solve; NOT reachable on
  the current walled ledger) }, with the boolean SUBTAG `CHERENKOV_DEFERRED` (carried on `SIM_GATED`: if the pinned `c_E<c_γ` a
  threshold Cherenkov channel opens — a throat-solve re-open obligation, GLM F3).
- `universality` ∈ { `SIM_GATED_REQUIRED_UNIVERSALITY` (production default — `b/ℓ` global is a REQUIRED condition, not yet earned) |
  `FALSIFIABLE_TENSION` (a species-indexed source is cited/found ⇒ forced non-universality) | `EARNED_SAFE` (an earned global source
  schema proves `b/ℓ` universal — not available today) }.
- `u_L_EP` ∈ { `SIM_GATED` (structure earned; magnitude gated via `a_L` AND `C_hu`) | `NO_GO` (stability/consistency break) |
  `EARNED_SAFE` (magnitude PINNED and shown negligible vs EP/fifth-force bounds — requires the deferred throat solve; NOT reachable
  on the current walled ledger) }, with
  the boolean SUBTAG `MIXED_SCALAR_EP_RISK` (fires when `C_hu·q_h·q_M` reintroduces a mass coupling; it is a FLAG carried ON
  `u_L_EP=SIM_GATED`, NOT a separate state — it escalates to a rule-2 `FALSIFIABLE_TENSION(channel=u_L_EP)` ONLY if an earned bound
  shows the mixing-EP appreciable WITHOUT gating, which needs the walled number ⇒ not reachable now).
- `preferred_frame` ∈ { `EARNED_SAFE` (an earned `c_E=c_γ` exists — NOT available, pathA_40 calibrated) | `PREFERRED_FRAME_TENSION`
  (the only radiation-suppression path requires `c_E/c_γ≫1` ⇒ a COST that blocks `NATURALLY_HIDDEN`; rides on the mapped verdict,
  NOT a rule-2 break by itself) | `SIM_GATED` (`c_E` free, no earned tie) }.

**Pipeline-level GUARD PREDICATES (booleans OUTSIDE the per-channel machine; each fires `NO_GO_CONSISTENCY` via rule 1):**
`LAUNDERING_REFUSED` (guard A tripped on the production path — a laundered magnitude was attempted) and `STABILITY_VIOLATED`
(`C_hu²≥B_eff K_h`). These are provenance/stability invariants, not channel states; rule 1 references them explicitly.

**Overall verdict (first match; references ONLY the five channels' tokens above + the two guard predicates):**
1. Any channel = `NO_GO`, OR guard predicate `LAUNDERING_REFUSED` (production path), OR guard predicate `STABILITY_VIOLATED`
   (`C_hu²≥B_eff K_h`) → **`NO_GO_CONSISTENCY`** (+ the conflicting relation / unsat core).
2. Else any channel in a `FALSIFIABLE_TENSION` state (`h_EP=FIFTH_FORCE_TRIGGERED`, `radiation=FALSIFIABLE_TENSION`,
   `universality=FALSIFIABLE_TENSION`, or an escalated `u_L_EP`) → **`FALSIFIABLE_TENSION(channel=…)`** (name every triggering channel).
3. Else if `h_EP=EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR` AND ≥1 remaining channel is gated/cost/audit
   (`radiation∈{SIM_GATED,AUDIT_ONLY_NOT_EARNED}`, `universality=SIM_GATED_REQUIRED_UNIVERSALITY`, `u_L_EP=SIM_GATED`,
   `preferred_frame∈{PREFERRED_FRAME_TENSION,SIM_GATED}`) → **`SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED`** — the expected honest
   landing (carries every gated/cost/audit/subtag flag: `MIXED_SCALAR_EP_RISK`, `PREFERRED_FRAME_TENSION`, `AUDIT_ONLY_NOT_EARNED`,
   `SIM_GATED_REQUIRED_UNIVERSALITY` as applicable).
4. Else `h_EP=EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR` AND `radiation=EARNED_SAFE` AND `universality=EARNED_SAFE` AND `u_L_EP=EARNED_SAFE`
   AND `preferred_frame=EARNED_SAFE` (no `SIM_GATED`/cost/audit anywhere) → **`NATURALLY_HIDDEN`** (extra scrutiny — essentially
   unreachable while the magnitude is walled; a landing here is suspicious and must be justified).
Every verdict carries: the Gate-L bending-mode connection (EARNED answer), the HARD_WALL statement (`M_h`/`c_E` not pinned + the
missing object + where it must come from), and the full five-channel table with subtags.

---

## 4. Controls (able-to-fail; each = SOURCE MUTATION → RECOMPUTED INVARIANT → EXPECTED TRANSITION)
- **A. Laundering guard (the #1 rig — HARD OUTPUT INVARIANT, headline):** the pipeline may emit a numeric magnitude ONLY IF an
  earned `h_time_kinetic_parent_action` provenance object exists. **General rule (GLM F5a):** ANY emitted numeric expression
  containing `M_h`, `c_E`, or `K_h` as a free symbol — INCLUDING the charge-residue floor `4Q_E²tanh²(b/ℓ)/(M_h b²)` and `P_h/P_EM`
  and any EP-magnitude — triggers the guard. Negative fixtures inject `M_h:=N₀`, `M_h:=K_parallel`, `c_E:=c_γ from Green form`, AND
  `K_h:=N₀·c_γ²` (reading the calibrated cone-lock + normalization as a derivation — GLM F5b) → each MUST recompute
  `LAUNDERING_REFUSED` (NOT a number, NOT a bare `SIM_GATED`). Baseline production (no such provenance) → the magnitude stays symbolic
  + `SIM_GATED`.
- **B. `h`-EP decoupled able-to-fail (headline):** nonzero `h` bare mass residue → `h_EP` recomputes `FIFTH_FORCE_TRIGGERED` even at
  `C_hu=0` → overall `FALSIFIABLE_TENSION`/`NO_GO`.
- **C. Mixing-EP able-to-fail:** fixture `q_L=0, C_hu≠0, q_M≠0` → `MIXED_SCALAR_EP_RISK` fires (mass residue enters the `h`-like
  eigenmode via `A_qm`).
- **D. Radiation-scaling able-to-fail:** corrupt the Green/flux speed → the DERIVED exponent changes (proves §2.2 is derived).
- **E. Preferred-frame able-to-fail:** set the ledger's radiation-suppression path to require `c_E/c_γ≫1` → `preferred_frame =
  PREFERRED_FRAME_TENSION` and `NATURALLY_HIDDEN` becomes unreachable.
- **F. Universality able-to-fail:** species-varying `b/ℓ` fixture → non-universality flag.
- **G. Stability:** `C_hu² ≥ B_eff K_h` → `NO_GO_CONSISTENCY`.
- **H. Static-Coulomb-match:** perturb the `h` static coupling off the calibrated `Q_E` → static mismatch flag.
- **I. `u_L` mass-coupling sign able-to-fail (GLM F7):** fixture `q_M → −q_M` → the computed `A_mm`/`A_qm` mass-coupling sign FLIPS
  in the `u_L` EP structure (proves the sign is derived, not stamped by convention).
- **J. Radiation-normalization able-to-fail (GLM F2):** a fixture applying the WRONG normalization (bare-fixed vs static-calibrated)
  → the derived exponent switches (3↔1), proving the two normalizations are distinguished in code, not conflated.
- **K. Conditional-`FALSIFIABLE_TENSION` able-to-fail (GLM F6, adjudicated):** a fixture supplying an explicit pinned-`K_h` ledger
  fact (exponent-1 branch) AND committing `c_E=c_γ` → `radiation=FALSIFIABLE_TENSION` (order-one, no hiding) — proves that break-state
  is reachable UNDER THE CONDITION. The SAME fixture WITHOUT the pinned-`K_h` fact (current ledger) → `radiation=SIM_GATED` even at
  `c_E=c_γ` (proves the tension is genuinely conditional on `K_h`-pinning, not asserted).
- **Baseline production run:** the real earned ledger → the five-channel table + overall verdict (a non-`SIM_GATED` computed
  magnitude would be preserved + trigger scrutiny — but must NOT appear given the wall + guard A).

---

## 5. Dual-engine (BINDING — split genuine DERIVATION from EXTRACTION+AUDIT; Codex Finding 3)
- **Genuinely dual-engine DERIVABLE (SymPy + independent Mathematica):** (a) the static 2×2 exchange matrix `A_qq,A_qm,A_mm` +
  the decoupled-vs-mixed EP structure (§2.1); (b) the radiation power-ratio EXPONENT from the explicit dipole/flux setup (§2.2); (c)
  the stability bound `C_hu²<B_eff K_h` + the `NO_GO` algebra; (d) the algebraic PROJECTION checks.
- **EXTRACTION + dual-engine algebraic PROJECTION (parse the frozen ledger fact, then both engines verify the projection — NOT an
  independent physical derivation of the fact):** the mass-residue-0 fact `source_mass_vector=[q_M,0]` (Codex Finding 3). (Provenance
  note, GLM F8: this fact is itself DERIVED in pathA_38 from the computed `grav_even_overlap=0`, `pathA_38_results.yaml:148` — the
  "extraction+projection" label is stricter than re-derivation and does NOT weaken the underlying provenance.)
- **EXTRACTION + PROVENANCE AUDIT (classify, cited `:line`):** the Gate-L failure-mode connection; the HARD_WALL provenance (`c_E`
  inserted, `M_h` calibrated-geometry); the `b/ℓ` global-vs-species status; the `u_L` magnitude SIM-GATED tag.
- **Bookkeeping (single-engine, cited):** the five-channel table, verdict formatting.

---

## 6. Scope boundary + ledger bookkeeping (do NOT modify earned specs)
- **This gate pins NOTHING new.** It does NOT execute Route-A; it CONFIRMS the wall. Therefore it does NOT trigger the NG5 forward
  re-open and does NOT rerun pathA_40 (no `SIM_DEFERRED`→`DERIVED` flip). State this explicitly.
- `SCALAR_DEPARTURE_MAPPED` is a CHARACTERIZATION of the pathA_39 departure; non-retroactive (does not alter pathA_38/39/40/41
  earned verdicts). The `h` decoupled EP-safety is a NEW earned sub-result about the departure, not a change to any sector spec.
- Everything sim-gated carries a re-open note (row → the throat solve → "re-adjudicate on solve"), mirroring pathA_41.
- **Astrophysical/cosmological re-open row (GLM F4):** on throat-solve pinning of `M_h`/`c_E`, the scalar coupling `g̃²=q_h²/M_h =
  4Q_E²tanh²(b/ℓ)/(M_h b²)` and speed `c_E` MUST be checked against stellar-cooling (RGB/HB/BNS) and cosmological (BBN/CMB) bounds on
  a massless charge-coupled scalar — an omitted future obligation the map records now (not a current-ledger channel).

---

## 7. Deliverables + review plan
- `tools/pathA_42_charge_coupled_scalar_{sympy.py,.wl}` (dual-engine: §2 core + §4 controls; exit 0; engine-agreement asserted).
- `reports/pathA_42_charge_coupled_scalar.md` (verdict line 1 = the overall verdict) + `_results.yaml` (five-channel table, the EP
  adjudication incl. the mixing caveat, the radiation exponent derivation + SIM-GATED magnitude, the preferred-frame cost, the
  universality constraint, the `u_L`/mixing structure, the Gate-L connection, the HARD_WALL provenance, controls).
- **Gauntlet:** this directive → Codex confirm-to-green → GLM-5.2 tertiary → Codex confirm-to-green → **user gate** → execute
  (Codex codes, Claude reviews) → arbiter re-run + transliteration-fidelity + adversarial-with-ablation (fresh agents) → **user
  gate** → commit. Never alter the calibrated process unilaterally.

---

## 8. Honest expectation
The likely landing is `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED`: the `h` static/EP sector is safe against the Gate-L mass-coupled
fifth force on the decoupled floor (a real positive — CONDITIONAL on the separately-gated `q_h/Q_E` universality requirement; NOT an
unqualified all-clear), but the departure is genuine and its size is gated — radiation (exponent-3 on the current
ledger, `SIM_GATED` on the free `M_h` even at `c_E=c_γ`; hidden only at the preferred-frame `c_E≫c_γ`), the `u_L`+`C_hu` mixing EP
channel, and the `b/ℓ` universality requirement all ride on the deferred throat solve. Report the earned map + the SIM-GATED
magnitude plainly; do not oversell "hidden" and do not launder a number. A `FALSIFIABLE_TENSION` (forced non-universality; or — only
CONDITIONALLY, if a future ledger pins `K_h` — the `c_E=c_γ` radiation) or a `NO_GO` (stability) must stay reachable.

---

## 9. Changelog
- **v1 → v2 (Codex design-review xhigh folded; verdict was `NEEDS_REVISION`).** BLOCKERs: (1) `h`-EP now uses the FULL static 2×2
  exchange matrix — `C_hu·q_h·q_M` mixing REINTRODUCES a mass coupling to the `h`-like eigenmode, so decoupled-floor safety is
  separated from the `SIM_GATED via a_L AND C_hu` mixing caveat (+ fixture C). (2) Laundering guard is now a HARD OUTPUT INVARIANT
  (§0/§4.A) with negative fixtures (`M_h:=N₀/K_parallel`, `c_E:=c_γ from Green`) → `LAUNDERING_REFUSED`. (3) NEW preferred-frame
  channel (§2.5) — the `c_E→∞` hiding limit costs Lorentz invariance; `NATURALLY_HIDDEN` impossible when hiding needs `c_E≫c_γ`. (4)
  Radiation scaling now REQUIRES an explicit dipole/far-zone/flux derivation in both engines (§2.2) with speed-corruption ablations,
  else downgraded to audit. Findings: (1) universality defaults to `SIM_GATED_REQUIRED_UNIVERSALITY`/`FALSIFIABLE_TENSION_OPEN`, not
  free `EARNED_SAFE`; (2) verdict grammar rebuilt on explicit five-channel tuples; (3) mass-residue-0 relabeled "extraction + dual-
  engine algebraic projection." NITs: `M_h` provenance aligned to pathA_41 (`CALIBRATED_GEOMETRY_INPUT`); "`h` = same
  embedding-direction family, distinct ledger object" (bridge unearned); canonical spelling `FALSIFIABLE_TENSION(channel=…)`.
- **v2 → v2.1 (Codex confirm-to-green fold; §3 was `NOT_GREEN` — state machine incomplete).** §3 rebuilt as a CLOSED per-channel
  state machine: every token defined + mapped (`FIFTH_FORCE_TRIGGERED`, `AUDIT_ONLY_NOT_EARNED`, `SIM_GATED_REQUIRED_UNIVERSALITY`,
  `PREFERRED_FRAME_TENSION`, `MIXED_SCALAR_EP_RISK` as a subtag-on-`u_L_EP=SIM_GATED`); §2.3 `FALSIFIABLE_TENSION_OPEN` → canonical
  `universality=FALSIFIABLE_TENSION`; §2.2 audit-fallback → `radiation=AUDIT_ONLY_NOT_EARNED`; out-of-set token ⇒ `NO_GO_CONSISTENCY`.
  Second confirm-pass fix: added the defined `EARNED_SAFE` states for `radiation` and `u_L_EP` (each reachable only with the pinned
  magnitude ⇒ post-throat-solve) so §3 rule 4 (`NATURALLY_HIDDEN`) references only legal tokens — a defined-but-unreachable-while-
  walled landing.
- **v2.2 → v3 (GLM-5.2 tertiary `SOUND_WITH_NITS` folded; 8 findings).** F1 (HIGH): `h`-EP decoupled safety is the CONJUNCTION of
  mass-residue-0 AND `q_h/Q_E` universality — non-universal `q_h/Q_E` gives an unscreened composition-dependent fifth force even at
  residue-0; token renamed `EARNED_SAFE_ON_DECOUPLED_FLOOR`→`EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR`, §2.1 conjunction added, §8
  softened. F2/F6 (radiation normalization — the central physics): the exponent is normalization-dependent (bare-fixed 3 vs static-
  Coulomb-calibrated 1); §2.2 now requires deriving BOTH + identifying the PHYSICAL one (is `K_h` pinned by calibration?) + the
  calibrated-`c_E=c_γ`→`FALSIFIABLE_TENSION` reachability + controls J/K; **routed to the final Codex confirm to adjudicate the
  `K_h`-pinned question against pathA_38.** F3: `CHERENKOV_DEFERRED` subtag (`c_E<c_γ` threshold). F4: §6 stellar-cooling/cosmological
  re-open row. F5: guard A strengthened (residue-floor + general `M_h`/`c_E`/`K_h`-free-symbol rule + `K_h:=N₀c_γ²` negative fixture).
  F7: control I (`q_M→−q_M` sign-flip). F8: §5 provenance note (mass-residue-0 derived from `grav_even_overlap=0`).
- **v3 → v3.1 (Codex final ADJUDICATION of the GLM F2/F6 normalization physics; state machine confirmed GREEN).** Codex checked
  pathA_38: the static Coulomb is `U_int ∝ q_eff²/N₀·1/(4πR)` — calibrated via `Q_E` + the zero-mode NORMALIZATION `N₀`, NOT via the
  dynamic stiffness `K_h`. So `K_h` is NOT pinned ⇒ **current-ledger branch = exponent-3, radiation `SIM_GATED` on the free `M_h`
  even at `c_E=c_γ`** — Codex OVERTURNS the unconditional GLM-F6 `FALSIFIABLE_TENSION` claim (GLM conflated `N₀` with `K_h`). §2.2/§3/
  §4.K revised: exponent-1/`FALSIFIABLE_TENSION` is CONDITIONAL on a future pinned-`K_h` ledger fact (default NO); control K now tests
  both branches. §0/§8 corrected. The honest landing returns to a genuinely-gated departure, not an order-one Lorentz-vs-radiation
  conflict.
- **v1** — authored from `_scratch/pathA_42_scalar_scope_codex.md` + the `HARD_WALL` reachability read.
