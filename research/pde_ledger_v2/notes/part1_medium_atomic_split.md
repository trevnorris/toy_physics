# Part I — The medium: FINALIZED atomic-stage split (user-approved 2026-07-07)

> Per-Part user gate (blueprint §7) satisfied 2026-07-07. This is the frozen split for Part I.
> Build order continues after the pilot (001 gravity, 002 gravity, 003 light) → Part I = build-order
> ids **004–007**; final Part-ordered linear renumber deferred to the pre-swap pass (manifest).
> Governing: `notes/ledger_v2_blueprint.md` §2 (granularity), §3 (Part I row), §5 (reshape spec),
> §6 (verification protocol). Source maps: two Explore fan-outs (continuum foundation + brane structure),
> 2026-07-07.

## The four stages (dependency order: substrate → speeds → phase structure → frozen brane action)

| Build id | Stage | Source gate(s) | Headline verdict token | Scope | Audit / reshape note |
|---|---|---|---|---|---|
| **004** | **I-1** GNLS parent action + dimensional foundation | `pathA_19` | `RETAIN_L_T_M` | EARNED: [L][T][M] dictionary (17 algebraic checks), derived `c_s0`, `ξ_h=√2·ħ/(m·c_s0)`, pin null relation `a=ħ/(m·c_s0)`; POSTULATED: parent action, `m_GNLS=M`; **labeled non-derivation:** `m_defect` NOT emergent (`INFLOW_MASS_SOURCE_MISSING`) | SymPy logic already exists+passes in the shared harness `src/stage1_solver/dimensional_check.py` (`run_patha19_foundation`, `--patha19-foundation`); reshape **extracts** it to a standalone `_sympy_audit.py` + re-authors an independent `.wl` |
| **005** | **I-2** EOS sound speed + light/sound ratio as a free knob | `pathA_20` (`pathA_20b` reference-only) | `C_GAMMA_RATIO_UNDERDETERMINED` | EARNED: `c_s²=5Kρ⁴/m` (rel. to EOS), `c=c_γ` wave-sector ceiling, three velocity scales; POSTULATED: EOS `P=Kρ⁵`; **CALIBRATED: `λγ=c_γ/c_s` unpinned (free), `λγ³` tail carried** | SymPy in harness (`run_patha20_velocity_constants`, `--patha20-velocity`); extract standalone `_sympy_audit.py` + independent `.wl`; optional able-to-fail teeth = pathA_20b forced-equality negative control |
| **006** | **I-3** Two-phase material-state ontology (order field `χ_B`) | Item A: `notes/brane_bulk_handoff.md` + `docs/conceptual_foundation.md` | `ACTION_SPECIFIED_CLASSIFIED` (structure) | POSTULATED (labeled): one conserved medium, two phases, `χ_B∈[0,1]`, order-state balance (`Γ_B` drain/return), throat=phase-conversion, χ_B gates shear/light; EARNED: recovery reduction (χ_B=1,Γ_B=0 → old projected-leakage law); NO-GO: `θ`-as-Maxwell-`φ` = FATAL_FLAW (carried dead-end) | **fresh-authored** dual-engine audit (prose-sourced; no existing pair). Scope (user-approved): dimensional homogeneity + recovery reduction (assert-zero) + `θ`-as-`φ` no-go (able-to-fail) |
| **007** | **I-4** Shear-surface G0 freeze — frozen medium action + drift ledger + DOF | `pathA_35` G0 | `T0_SHEAR_FROZEN(d9520d3819c3)` + `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)` | KEPT: GNLS + T0 polar-OP; POSTULATED/CALIBRATED: the "11" (4 constants {ρ_br, μ_R, λ_Pu, Ω_w} + 1 function g_ℓ + 6 structural postulates); EARNED: flat-brane linear DOF=8 (dual-engine, able-to-fail), dimensional firewall | reshape existing dual-engine pair (`pathA_35_G0_sympy.py` + `pathA_35_G0.wl`, ENGINE_AGREE). Carry the 2026-07-04 erratum: the "11" STANDS; irreducible cross-sector drift `{ρ_B0,χ_c,C_hu}` is a Part-VI (pathA_41) item |

## Settled decisions (user, 2026-07-07)

- **pathA_23_stage0 (EM-native brane-elastic fork) DROPPED to history** — NOT a Part I stage. Reconciles a
  blueprint §3-vs-§4 contradiction: §3's Part I table listed it as a source, but §4's DROPPED list says
  "pathA_23 (EM-native) … stays as history, not a stage." Its MacCullagh *menu* was superseded by the
  pathA_35 G0 freeze (which made the choice). The two-phase *ontology* (I-3) comes from the χ_B material-state
  closure (`brane_bulk_handoff.md` + `conceptual_foundation.md`), a physically distinct artifact from
  pathA_23_s0's brane-elastic displacement action. → Part I = **4 stages**, no menu-vs-freeze overlap.
- **I-3 audit scope** = dimensional homogeneity + recovery reduction + `θ`-as-`φ` no-go (fresh-authored, real
  executable stage — meets the uniform-audit standard + goal-1 context capture).

## Notes / cross-Part dependencies

- **pathA_19/pathA_20 are NOT truly `.wl`-only.** The SymPy verification already exists and passes inside the
  shared `dimensional_check.py` harness; the reshape *extracts* per-stage standalone audit scripts (which we
  author fresh regardless). Not "missing verification."
- **I-4's `T0_SHEAR_FROZEN(d9520d3819c3)` is load-bearing provenance already cited by the built `ledger_stage003`
  (Part III)** and downstream Part IV — I-4 is its formal home. No conflict; Part I formalizes what the pilot cites.
- **I-3 vs Part 0 (conceptual prose):** no redundancy. Part 0 = plain-language vision; I-3 = the formalized,
  dimensionally-audited postulated action (free-energy functional + order-state balance PDE + recovery reduction).
- **I-2 running-start:** the detailed pathA_20/20b atomic source map (claim-set, scope, carried residuals, reshape notes) is in
  `research/pde_ledger_v2/notes/stage005_pathA20_source_map.md` — read it to author the I-2 directive without re-running discovery.
- **I-3 running-start:** the χ_B ontology source map (the postulated action + the 3 audit legs {dim, recovery-reduction, θ-as-φ no-go} +
  labeled postulates + the 10 modeling choices to pin) is in `research/pde_ledger_v2/notes/stage006_chiB_ontology_source_map.md`. ⚠ I-3 is
  **FRESH-AUTHORED** (no existing script pair) — the directive SPECIFIES the action, it does not extract a harness. Read the source map first.

## Progress

- **I-1 `ledger_stage004` DONE** (committed `d9544a62`) — dual-engine (SymPy 49 / Mathematica 50 PASS) + full tri-review CLEAN
  (`FIDELITY_CLEAN` + `ADVERSARIAL_CLEAN`). `RETAIN_L_T_M` earned; carried gaps first-class.
- **I-2 `ledger_stage005` DONE** — dual-engine (SymPy 80 / Mathematica 81 PASS) + full tri-review CLEAN (`FIDELITY_CLEAN` +
  `ADVERSARIAL_CLEAN`, 14-mutation ablation matrix — every tooth genuine, verdict computed + reversible) + post-remediation
  `REVERIFY_CLEAN`. `C_GAMMA_RATIO_UNDERDETERMINED` — sound speed `c_s²=5Kρ⁴/m` derived from EOS; `c=c_γ` ceiling; `λγ` a free
  calibration input (reversible to `EQUALS` only with an inserted source equation). Consumes I-1 dictionary + `EOS_FROM_GNLS_FACTOR`.
  Registration at count 5, PDF rebuilt.
- **▶ NEXT = I-3 `ledger_stage006`** (two-phase χ_B material-state ontology; fresh-authored dim + recovery + θ-as-φ no-go audit).
  Then I-4 `ledger_stage007` (pathA_35 G0 freeze). Then Part II (gravity, the swing).

## Per-stage process (unchanged, calibrated)

Codex xhigh design-review of the reshape directive → fold to CLEAN (no GLM tertiary on Parts I–VI) →
pre-exec user gate → Codex builds only the two scripts (`--sandbox danger-full-access`, background, `< /dev/null`,
xhigh) → dual-engine both exit 0 → orchestrator arbiter re-run → full tri-review on fresh agents (arbiter +
fidelity + adversarial-scoped-to-reshape-integrity, with ablation) → remediate → bump counts → **update the parameter
register (`../parameter_register.md`) → Codex-verify the register update → fold — MANDATORY every stage (blueprint §5
step 7)** → commit + docs/memory sync.
Orchestrator authors notes/cards/LaTeX/registration; Codex codes.
