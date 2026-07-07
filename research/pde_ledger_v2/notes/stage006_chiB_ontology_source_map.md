# I-3 (ledger_stage006) source map — two-phase material-state ontology (order field χ_B)

> Running-start prep for the NEXT stage (I-3), captured 2026-07-07 from a source-mapping fan-out so a fresh session can
> author the directive without re-discovery. Verify against the cited sources before finalizing the directive.
> Companion: `part1_medium_atomic_split.md` (row 006). Build-order id **006**, Part I.
> **⚠ FRESH-AUTHORED stage — NOT a reshape.** There is NO existing script pair; the audit is authored from PROSE sources.
> The directive scope shifts from "faithful extraction of a harness" to "author a real dual-engine audit of a postulated
> action" — so the directive must SPECIFY the action/postulates to encode (§1 below), not point at a report's checks.
> **Headline verdict:** `ACTION_SPECIFIED_CLASSIFIED` (structure). **Scope (user-approved, `part1_medium_atomic_split.md:16,27`):**
> dimensional homogeneity + recovery reduction (assert-zero) + θ-as-Maxwell-φ no-go (able-to-fail). Does NOT earn light.

## File inventory (READ-ONLY prose sources — content inlined into the note by the orchestrator)
- `notes/brane_bulk_handoff.md` — the χ_B closure: density split, conservation, order-state balance PDE, projected law, energy ledger, Task-A postulate list.
- `docs/conceptual_foundation.md` §1–§4 — the physical picture (one medium, two phases; why a double-well is needed since GNLS `U(ρ)∝ρ⁵` is single-well).
- `notes/rung_W_reframe.md` — the GLM pre-check that PASSED rung W (a double-well χ_B wall is stable AND light-permitting); concrete well forms; the input-count/drift.
- `software/stage1_solver/reports/pathA_36_c5_phase_potential.md` — the COMPUTED θ-as-φ NO-GO (Dirac–Bergmann bracket, sign argument, `BY_TUNING` control) that I-3 carries as the able-to-fail FATAL_FLAW.
- Recovery target (frozen firewall, OLD ledger): `research/pde_ledger/paper/stages/stage_243.tex` / `stage_244.tex` (`stage_244.tex:36-39`, the canonical `S_leak` "Gaussian one-mode leakage law") — assert-zero against THIS, not just the handoff restatement.

## Leg 1 — the postulated χ_B action + order-state balance PDE (SPECIFY, all labeled POSTULATED)

**Free-energy functional** (`brane_bulk_handoff.md:420-451`):
```
F = ∫ d⁴X [ ½ n|u|²  +  U(n)  +  f_B(n,χ_B)  +  (κ_B/2)|∇₄χ_B|²  +  χ_B f_shear  +  f_throat  +  f_mix ]
```
- `U(n)` = I-1's KEPT stiff GNLS/EOS bulk energy `(K/4)ρ⁵` (single-well).
- `f_B(n,χ_B)` = **POSTULATED double-well** (the crux new object). Candidate forms: φ⁴ `½a(χ_B²−1)²` (`rung_W_reframe.md:203`) or Landau `½α(T−T_c)χ_B²+¼βχ_B⁴` (α<0). **See trip-up T2** on the [0,1]-vs-±1 minima convention.
- `(κ_B/2)|∇₄χ_B|²` = **POSTULATED** interface/surface-tension gradient cost.
- `χ_B f_shear` = **POSTULATED** multiplicative gate: shear/light exists only where χ_B≠0; explicit `½χ_B μ_R(∇×u)²` (`rung_W_reframe.md:60`). Ties to I-2's `c_γ²=μ_R/ρ_br`.
- `f_throat`, `f_mix` = mouth/wall + mixed-gauge couplings (ontology-only for I-3; largely deferred).

**Density split** (`brane_bulk_handoff.md:172-186`): `n` = total conserved 4D constituent density; `n_B=χ_B n` = brane-ordered density.
**Total conservation (exact):** `∂_t n + ∇₄·(n u)=0`.
**Order-state balance PDE (the Γ_B drain/return term)** (`:203-208`):
```
∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B ,     Γ_B = Γ_return − Γ_drain ,   J_χ=0 in the simplest closure
```
Sign: `Γ_B<0` brane→bulk (de-structuring/drain); `Γ_B>0` bulk→brane (re-ordering). **Throat = phase-conversion** (site of `Γ_drain>0`, χ_B:1→0) — an **admittance/outlet, NOT a suction force** (`:88-90`); drain is μ-gradient-driven `J_repair~−M_n∇μ` (`:96-100`), not a new pairwise attraction. Optional dissipative dynamics `D_tχ_B=−M_χμ_χ+Γ_throat+Γ_return` with `μ_χ=δF/δχ_B`; energy-ledger requirement `P_order=∫μ_χ n Γ_B` — do not silently dissipate (`:479-499`).

## Leg 2 — the recovery reduction (EARNED; assert-zero)

Limit `χ_B=1, Γ_B=0, J_χ=0` (`brane_bulk_handoff.md:305-311`, Test 1 `:820-836`). The general projected law (`:264-292,676-707`)
`∂_tρ_B+∇₃·J_B=S_flux+S_convert` (with `ρ_B=∫W χ_B n dw`, `S_convert=∫W n Γ_B dw`, `S_flux=−[W(χ_B n u^w+J_χ^w)]_∂+∫W'(w)(…)dw`)
must reduce to the **old projected-leakage law** `∂_tρ+∇₃·J=S_leak`, `S_leak=−[W n u^w]_∂+∫W'(w)n u^w dw`.
**Assert-zero form:** `(S_flux+S_convert)│_{χ_B=1,Γ_B=0,J_χ=0} − S_leak ≡ 0` symbolically, against the canonical frozen `S_leak` (stage_244.tex). EARNED relative to the imposed χ_B split + declared `W(w)`.

## Leg 3 — the θ-as-Maxwell-φ no-go (able-to-fail FATAL_FLAW; carried dead-end)

Claim being killed: a complex `χ_B=|χ_B|e^{iθ}` whose phase θ serves as the MacCullagh/Maxwell scalar φ resolving C5. Verdict = **`FATAL_FLAW`**, do NOT carry into a rung-φ.
Sign argument (`pathA_36_c5_phase_potential.md:12-24,44-70`; `rung_W_reframe.md:153-184`): a **stable order-parameter phase has POSITIVE gradient energy** → `K_θ=−κ_phase<0` in the Lagrangian convention (couples to its conjugate density via Josephson `∂_tθ=−μ`, NOT to `∇·u`). **Maxwell's electric square needs the OPPOSITE sign** `K_θ=J²ρ_B0²/ρ_br>0` (`B_eff=0`, `m_θ²=0`). The only sign-flip (negative gradient stiffness) = a finite-k **Lifshitz instability = the pathA_25 wall** (`k_Rstar²=40Kmρ0⁴/ħ²`), itself already falsified (`RC…FAIL_LIGHT_STARVED`). So θ-as-φ re-enters a killed wall.
**What the audit computes + able-to-fail:** Dirac–Bergmann of the longitudinal `(u_L,θ)` system: `Φ₁=π_θ−Jkρ_B0 u_L`, bracket `{Φ₁,Φ₂}=k²(J²ρ_B0²+κ_phase ρ_br)/ρ_br`. Provenance-fixed → bracket≠0 → two **second-class** constraints → **1 stray longitudinal DOF** → `FAIL_CAUCHY_STRAY_LONGITUDINAL` (or `FAIL_C5_LONGITUDINAL_ZERO_MODE` for `B_eff=0`). **Would-PASS fixture (the teeth):** tuned `K_θ=J²ρ_B0²/ρ_br, B_eff=0, m_θ²=0` → bracket=0, first-class, 0 DOF → `C5_RESOLVED_MAXWELL_BY_TUNING`, flagged **`BY_TUNING`, NOT `WITH_PROVENANCE`** (the frozen χ_B defs do not force `K_θ>0`, its value, or removing the `ρ_B0²/χ_c` term). Detuning any of the three fires FATAL (controls all fire, `:120-138`).
**Simplest fresh-authored version (if full Dirac–Bergmann is heavier than needed):** compute the SIGN of a stable double-well phase's gradient stiffness (positive gradient energy → `K_θ<0` in L) vs the Maxwell-required electric sign (`K_θ>0`), assert OPPOSITE → FATAL; the PASS-fixture is precisely a negative-stiffness (Lifshitz/pathA_25) phase. Transverse control stays positive (2 photons `c_γ²=μ_R/ρ_br`; `ε≠ρ_br` → `FAIL_TRANSVERSE_DISTURBED`) — ties to I-2.

## Dimensional targets + consumed I-1/I-2 dictionary
χ_B is dimensionless `[1]`; `F`=energy over `d⁴X` → every functional term is a 4D free-energy density `M L⁻² T⁻²`. Consumed (CITED, not re-derived): `{L,T,M}`; `n,ρ0=L⁻⁴`; `m_GNLS=M`; `K=M L¹⁸T⁻²`; `c_s0=√(5Kρ0⁴/m)=LT⁻¹` (DERIVED); `c_γ=LT⁻¹`, `c_γ²=μ_R/ρ_br`. Forces: `[a]=[α]=[β]=M L⁻²T⁻²`, `[κ_B]=M T⁻²`; balance PDE → **`[Γ_B]=T⁻¹`**, `[J_χ]=L⁻³T⁻¹`; relaxation → `[μ_χ]=M L⁻²T⁻²`, `[M_χ]=L²T M⁻¹`. Firewall ablations (pattern `pathA_36…:141-152`): drop `ρ_B0`/omit the ∇ in κ_B|∇χ_B|²/χ_c multiply-vs-divide → break homogeneity.

## Verdict tokens + honest scope
Headline `ACTION_SPECIFIED_CLASSIFIED` (structure). Recovery sub-token EARNED (propose `RECOVERY_REDUCTION_VERIFIED`). No-go carried `FAIL_CAUCHY_STRAY_LONGITUDINAL` (+ `FAIL_C5_LONGITUDINAL_ZERO_MODE`), positive control `C5_RESOLVED_MAXWELL_BY_TUNING` flagged `BY_TUNING`. Drift `DRIFT(6)`/`SECOND_MEDIUM_DRIFT`. Dual-engine `ENGINE_AGREE`. **Landing:** POSTULATED microstructure (labeled) + EARNED recovery + carried no-go; classifies the action's structure and encodes the dead-end; does NOT earn light/Maxwell (that is the point of the no-go).

## New knobs for the parameter register (χ_B row absent today; add on build)
| Param | dim | Class | Note |
|---|---|---|---|
| `χ_B` | `1` (∈[0,1]) | ACTION (postulated OP) | independent scalar OP, NOT χ_B=|P_∥|² (load-bearing, T4) |
| `a`/`α`,`β` (well) | `M L⁻²T⁻²` | ACTION | POSTULATED double-well (GNLS is single-well) |
| `κ_B` (gradient) | `M T⁻²` | ACTION | surface-tension-like |
| `Γ_B`=`Γ_return−Γ_drain` | `T⁻¹` | ACTION/GAP | conversion rate; global-return-constrained `R_0=−M_0,R_1=−D_1` |
| `M_χ` (mobility), `J_χ`(=0), `α_aniso` | `L²T M⁻¹` / `L⁻³T⁻¹` / `M L⁻²T⁻²` | ACTION | relaxation / transport / P-orientation anisotropy |
| θ-no-go branch (DEAD): `θ,J,ρ_B0,K_θ,χ_c,B` | pathA_36 | carried no-go | do NOT admit as live knobs. Note `ρ_B0,χ_c` already appear as pathA_41 Part-VI drift — cross-reference |
Minimum live drift for real χ_B = **1 field + 4 constants + 1 gate = 6 inputs** (`rung_W_reframe.md:140`) → `DRIFT(6)`.

## Modeling choices the directive MUST pin (trip-ups)
1. **n (number, L⁻⁴) vs mass density** — `½n|u|²` needs mass; pin the density variable + its dim or the free-energy-density check is ill-posed (`brane_bulk_handoff.md:188`).
2. **Double-well form + range** — χ_B∈[0,1] (minima 0,1) vs `½a(χ_B²−1)²` (minima ±1). Pick one (e.g. `f_B=½a χ_B²(1−χ_B)²` for [0,1]) so the recovery limit χ_B=1 sits at a minimum.
3. **Same-density degeneracy** — both minima of `f_B(n,χ_B)` at the SAME n (else it's a liquid-vapor density interface = no light, already rejected). Label POSTULATED.
4. **χ_B independent of P** — the T1-escape rests on this; postulate explicitly as the single load-bearing assumption.
5. **Conservative vs dissipative order dynamics** — choose; decide if the energy ledger `P_order` is in-scope (likely a labeled adjunct, since audited scope = dim+recovery+no-go).
6. **`J_χ=0` default** for a clean recovery assert-zero; flag J_χ≠0 as deferred.
7. **Recovery target = the canonical frozen `S_leak`** (stage_243/244.tex), not the handoff restatement.
8. **Fix the sign convention explicitly** (energy-vs-Lagrangian) so "stable phase → K_θ<0 in L" vs "Maxwell electric K_θ>0" is unambiguous.
9. **Global return `R_0=−M_0,R_1=−D_1`** = global-moment closures, labeled postulates — do NOT assert-zero them locally in I-3.
10. **Throat=phase-conversion is ontology, not a solve** — encode the claim; defer the moving-throat branch solve. Gate-L translation-Goldstone hazard (`δw=u_w`) is out-of-scope, flag deferred.

## Process (fresh-authored variant)
Author directive (SPECIFY the action + 3 legs + labeled postulates + consumed dictionary) → Codex xhigh design-review → fold to `DIRECTIVE_CLEAN` (no GLM on Parts I–VI) → pre-exec user gate → Codex builds the two scripts (fresh audit, independent `.wl`) → dual-engine both exit 0 → arbiter re-run → full tri-review (arbiter+fidelity+adversarial-with-ablation) → remediate → bump 5→6 → **update + Codex-verify the parameter register** → note/card/`\input{stages/stage_006}`/registration → PDF → commit + docs/memory sync. Target stem: `ledger_stage006_two_phase_chiB_ontology` (confirm slug on authoring).
