# STATUS — where the Path-A program is (single front door)

**This file is the canonical "you are here."** It is a thin pointer, not a copy — the detail lives in the linked docs.
Updated at every milestone (same moment `software/stage1_solver/decisions/13` §0 is updated). Last update: **2026-07-03.**

> **New to the model / need the physical picture? Read `docs/conceptual_foundation.md` FIRST.** It is the plain-language,
> native-terms statement of what the medium, the brane, the four sectors (gravity=drain, magnetism=swirl, electric
> charge=puncture direction ±w, light=brane shear), and the particle/defect (an extended throat — there are NO point particles)
> physically ARE. This `STATUS.md` is the *program state*; that doc is the *conceptual vision*.

---

## ▶ RESUME HERE (2026-07-03, post-`/compact`) — magnetism is next

**⭐ 3 OF 4 FORCE-SECTORS EARNED, all from ONE brane+bulk.** Gravity (`pathA_29` localization + the PN ladder), light (`pathA_36`),
and electric charge (`pathA_38` = `THROAT_ELECTRIC_LOCALIZED_COULOMB`, just banked — commit `797f4d88`) — each a distinct thing the
same medium's brane+bulk does, all sharing **one localization principle** (a normalizable transverse zero mode → `1/r²`), split by
`w→−w` **parity** (even drain = gravity, odd embedding-Goldstone = charge → that parity is *why* charge ≠ mass). **The clean
end-to-end derivation chain = `docs/conceptual_foundation.md` → ⭐ THE FOUR-SECTOR CHAIN** (read it — it maps brane+bulk → all four
sectors with each link's status).

**▶ NEXT = MAGNETISM (close Maxwell).** Draft carry-over spec: `software/stage1_solver/directives/pathA_39_magnetism_close_maxwell.md`.
The gate: a MOVING throat → brane swirl/vorticity → **B**, the Lorentz force, the curl equations, with the right `E↔B` normalization —
largely DOWNSTREAM of `pathA_38` (charge/E) + `pathA_36` (light/the EM wave). **First task = resolve its §2 setup sub-questions**
(what the swirl field is + how a moving throat sources it; does E+B close full Maxwell; the Lorentz force; `λγ` consistency), THEN the
gauntlet (Codex design-review → GLM → dual-engine → tri-review), reusing the pathA_38/pathA_29 machinery.
**⚠️ REMEMBER the pathA_38 lesson:** a clean PASS can be **pass-by-construction** — pathA_38's first run hardcoded `1/(4πR)` and was
structurally unable to FAIL; only the **adversarial-with-ablation** leg (force bad inputs → the classifier must EMIT a `FAIL`, not
raise) caught it. Make every gate empirically able-to-fail ([[feedback-negative-verdict-short-circuit]]).

**▶ THEN = THE CONSISTENCY KNIT → "fully formed."** After all four sectors: (a) the sharp `λγ = c_γ/c_s = 1` cone-lock test (still
OPEN — do all four live in ONE parameter set?); (b) the NG cross-consistency gauntlet; (c) assemble the end-to-end chain into the
central `pde_ledger` (the calibrated PDE delivering GR+EM). This closes *"four sectors that each work"* → *"one self-consistent
brane+bulk model."* (Completes the SPEC; the full nonlinear sim stays deferred — a no-go is still possible.)

---

## ⭐⭐ GRAVITY-SECTOR REFERENCE (2026-06-26) — still valid

**DYNAMICAL-GRAVITY SECTOR = BUILT & GR-MATCHED; speed-of-gravity / aberration worry RESOLVED.** The conservative PN two-body ladder
**1PN→4PN + 2.5PN radiation is already derived, audited, and GR-matched** (calibrated / controlled-reduction) in `research/4d_*pn*`.
**DO NOT re-derive it.** (memory `project-pn-gravity-ladder`.)

**END GOAL = a fully CALIBRATED PDE delivering GR + EM.** Calibration is fine; first-principles is NOT required; **existing-in-any-shape
is the win.** The central spec is `research/pde_ledger/` (the 253-stage audited ledger of the target moving-throat PDE); every
calibration is a constraint that feeds it. (memory `project-calibrated-pde-goal`.)

**pathA_28 (monopole/dipole radiation) — DONE = `MONOPOLE_DIPOLE_RETURN_CONDITIONAL`.** A VERIFIED CONSTRAINT-SPEC (dual-engine;
arbiter PASS + fidelity CLEAN; adversarial CONCERNS = it is a constraint-spec, not a falsifiable test). **Handoff:** to avoid
GR-forbidden monopole/dipole gravitational radiation, the brane↔bulk return must deliver `R0 = −M0` (net mass-rate `M0 = ∫S_leak`)
and `R1 = −D1` (net dipole/momentum-rate `D1 = ∫x_i S_leak + ∫j_i`, including the carried odd wake). Artifacts under
`software/stage1_solver/` (tools / reports `pathA_28_monopole*`).

**⭐ TRACK 3 GATE-1 (brane↔bulk return, `pathA_29`) — DONE + VERIFIED (2026-06-25) = `RETURN_RESIDUAL_PREDICTION`.** The keystone
pathA_28 handed off — earned on the **4th execution (v3b)** after a full tri-review (orchestrator arbiter re-run reproduced +
adversarial CLEAN + fidelity FAITHFUL). **Given the drain premise (`Z<0` = gravity IS the inflow): 1/r² Newtonian gravity SURVIVES
the finite slab** — both admissible DC-sink return completions (de-structuring/absorbing + Bloch stack) genuinely solve a normalizable
`m=0` transverse zero mode → `p=2` via a real 3D-radial `dsolve` (counterfactual-guarded: a wrong `1/r⁵` → nonzero residual,
rejected). **And the drain comes bundled with an UNAVOIDABLE bounded monopole/dipole `c_s`-radiation residual `∝ ε0 = 1−𝒯₀(0)`**,
tied to the gravity strength — **the falsifiable departure from GR** (GR forbids monopole/dipole GW via Birkhoff/mass-conservation;
the drain breaks brane mass conservation). NOGO is genuinely reachable via a derived delocalizing warp (`p=3`). **Sharpens but does
NOT close `pde_ledger` open-item #9** (records the residual-radiation prediction); the **gravity-range (1/r²) item passes** for the
localizing flat-slab family. Artifacts: `software/stage1_solver/{directives/pathA_29_brane_bulk_return.md (v3), tools/pathA_29_*,
reports/pathA_29_*}`. **⭐ THE ACTIVE PUSH (user, 2026-06-25: "push until a wall"): complete the full nonlinear moving-throat PDE / brane↔bulk return
closure.** Central spec = `research/pde_ledger/` (253 stages, ALL algebra done; single OPEN item = **"actual branch realization"** — the
solved nonlinear PDE must RETURN the grouped-`P2` / quadrupole-normalization data `m̂₀²P₀ = 54Gc_s⁵/(5a⁵c⁵)`). **This is a ~6-gate
ladder, NOT one step** — master checklist = `research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md`.

- **Step 1 (reconcile the two open-item framings) ✅ DONE** = `moving_throat_pde_open_item_reconciliation.md`: pathA_29's "#9 /
  `R0=−M0`" (ℓ=0/1 forbidden channels) and the ledger's "branch realization" (ℓ=2 quadrupole) are **complementary multipole sectors of
  one return law**, not the same and not independent; a new cross-ℓ consistency gate falls out (Gate 5).
- **Gate 1 (`pathA_30`, frozen-wall D/N unit test) ✅ DONE = `DN_UNITTEST_BC_DEPENDENT`** (tri-review CLEAN: arbiter re-run both engines
  + fidelity + adversarial). The geometry lift's frozen reduction genuinely gives the const-coeff Helmholtz resonator (DtN
  `−(ω/c_s)tan`, half-shifted ladder, Robin counterfactual, dual-engine agreed; Bogoliubov `k⁴` deferred under `kξ≪1`). **D/N is a
  banked calibration input** (`bc_derivation_emitted=False`); earning it → `PASS` is an optional later upgrade. 2 non-blocking NITs
  logged in the ladder. Artifacts: `software/stage1_solver/{directives,tools,reports}/pathA_30_*`.
- **Gate 2 (`pathA_31`, scalar breathing `η₀₀`) ✅ DONE = `BREATHING_CALIBRATED`** (tri-review CLEAN after a remediation). The
  distributed ℓ=0 reduction genuinely reproduces the legacy `(a,L)` collective closure: derived `L₀`-harmonic profiles, `M_AB/K_AB`
  from operator projection (not `∂²E_geom`), two genuinely-independent HF routes match, structure gate computed + able-to-fail,
  dual-engine agreed. The 2-mode truncation is clean (`V₂`-overlap `o_1=0.993` at the physical `β_L0=1.85`) across the order-unity
  wall-stiffness window (`K_η/T_w≲2.6`), failing only for sharp walls. CALIBRATED ⇐ `μ_η,T_w,K_η` are calibration inputs. **Caveat:**
  the overlap doesn't guard profile-correctness (HF mismatch does); clean truncation shown *for assumed order-unity wall stiffness*.
  v1 was tri-review-REJECTED (HF `x−x`, hardcoded counterfactual, gamed threshold) → remediated. Artifacts: `pathA_31_*`.
- **Gate 3 (`pathA_32`, grouped-`P2` / ℓ=2 sector) ✅ DONE = `ISOTROPY_CALIBRATED`** (full tri-review CLEAN after a remediation). The
  distributed lift's ℓ=2 grouped-`P2` sector satisfies the isotropy / lane-degeneracy theorem at the linearized isotropic reference:
  the three grouped lanes `{20,21,22}` collapse to common conservative coefficients (raw-`D` lane defects = 0, the **PRIMARY** gate),
  the reduction is SO(3)-covariant (angular Gram = `I₅`; computed `−Δ_S²` eigenvalue `λ_m=6` per harmonic; the angular stiffness
  `K₂=∫[T_w β₂'² + (K_η+6T_Ω)β₂²]` — the term that dropped at ℓ=0 — is now alive), and the gate is genuinely able-to-fail (8
  counterfactual probes, each flips under ablation). **The quadrupole sector first appears here.** CALIBRATED ⇐ the wall constants
  `μ_η,T_w,K_η,T_Ω,β₂` + the symbolic radial scalars are calibration inputs. **Two-tier gate:** raw-`D` lane collapse is the
  verdict-bearing PRIMARY test; the published `a₂=b₂=a₄=b₄=0` (normalized `u`-defects) is a necessary-but-not-sufficient cross-check
  (a per-lane order-independent prefactor cancels in the normalized response). v1 was tri-review-REJECTED — the
  adversarial-with-ablation leg caught two pass-by-construction defects fidelity missed (dual-engine gaming = vacuous `x−x` on three
  byte-identical Mathematica copies; 5/8 probes typed their FAIL booleans + a tautological able-to-fail aggregate) → remediated
  (honest per-harmonic `.wl` + per-lane `D` cross-engine; each probe computed-from-mutated-input with a self-ablation) → re-verified
  (arbiter byte-for-byte; fidelity + adversarial re-ablation EARNED). Artifacts: `pathA_32_*`.
- **Gate 4 (`pathA_33`, quadrupole `54/5` normalization) ✅ DONE = `QUAD_CALIBRATED`** (full tri-review CLEAN after a remediation). On
  Gate-3's isotropic outgoing branch, the ℓ=2 sector delivers the **EARNED predictive surplus cleanly SEPARATED from the CALIBRATED
  magnitude**: the outgoing fingerprint `1/9, 4/81, 1/27` is **DERIVED** from the DtN `Ŷ₂^out=−3/Λ₂^out` Hankel series (not hardcoded),
  the prefactor algebra follows from `P(ω)=D₀N(ω)/D^cons(ω)²` (the `−2D₂N₀` term), `P0_target_scaling=a⁻⁵`, and the derived `χ_Q=1`
  closes `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵ ⟺ γ_quad^eff=2G/5c⁵`. **The headline result — the earned/calibrated split:** `54/5 = 2·27/5` — the
  `27` is **EARNED** (`derived_in_gate`, from the fingerprint), the `2/5` + `G` + assembled magnitude are **CALIBRATED**
  (`external_bridge_input`), classified by a 4-way **PROVENANCE** partition (NOT `G→λG` invariance, since `54/5` is G-invariant yet
  calibrated). `G` is `GENUINE_BLOCKED` — the PDE delivers the FORM/branch, **not** Newton's `G`. **Dimensional milestone:** a genuine,
  `μ̂₀`-FREE, able-to-fail dim check (`[P₀^phys]=(c_s/a)²·N₀/D₀` must be dimensionless from sourced `[N₀]=L⁻¹M`, `[D₀]=L⁻¹T⁻²M`) caught
  that the handoff's `P₀=N₀/D₀` silently drops a `(c_s/a)²` factor (natural-units trap). v1 was tri-review-REJECTED — the
  adversarial-with-ablation leg caught two pass-by-construction sub-controls fidelity + arbiter both passed: (a) the dim gate was STILL
  tautological (it **back-solved the FREE carrier `μ̂₀`** to force homogeneity); (b) the per-probe `self_ablation` was a **constant, not a
  re-run** → remediated (`μ̂₀`-free dim gate + corrupt-`[N₀]` probe §3d′ + real two-verdict `self_ablation`) → re-verified EARNED
  (corrupting `[N₀]` now fires `FAIL_DIMENSIONAL`; `[G]` doesn't). Artifacts: `pathA_33_*`.
- **Gate 5 (`pathA_34`, scalar/dipole + cross-ℓ unification) ✅ DONE 2026-06-26 = `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (EARNED).**
  The ℓ=0/1 brane↔bulk return (pathA_28/29 — `R0=−M0`, `R1=−D1`) and the ℓ=2 quadrupole (`P₀`, derived via the port kernel) are
  sectors of one return law, but at LINEAR order the cross-ℓ link runs through the return admittances `Z0_ret,Z1_ret`, which pathA_29
  establishes are a **premise** (`Z_is_premise:true`). A rank audit over 11 genuine generator dofs computes `return_moving_nullity=2`
  (ℓ=2 determined, ℓ=0/1 NOT) → underdetermined; a **real derived selector equation** collapses it to nullity 0 →
  `CROSS_L_RESIDUAL_PREDICTION` (so the gate is genuinely able-to-fail). **Honest landing: the linear theory cannot pin the ℓ=0/1
  return — Gate 6 must supply a selector fixing `Z0_ret,Z1_ret` (the first concrete, proven Gate-6 input).** ⚠️ **v1 tri-review-REJECTED**
  — the ADVERSARIAL leg caught the verdict mechanism pass-by-construction (rank audit rigged-to-UNDERDETERMINED via a zero-padded
  constraint; flag-driven probes; `.wl` headline-only). Remediated (genuine 11-dof rank audit from the §8.2/§9.4/§10 reductions + real
  selector-equation substitution + computed-from-mutation probes + dual-engine verdict machinery) and re-verified EARNED (arbiter
  reproduces under the new toolchain; fidelity CLEAN; adversarial re-ablation EARNED). Earned content (`−(ℓ+1)/Λ_ℓ` fingerprints
  `ω¹/ω³/ω⁵`, residual forms vs pathA_29, free-carrier-independent dim check) all genuine. Artifacts: `pathA_34_*`.
- **⭐ NEXT (chosen 2026-06-26) = the BRANE sector** — the gravity arm is at its sim-ready boundary (Gates 1–5 done; Gate 6 = the
  WALL, its 4D field solve sim-deferred), so the program moves to the brane: it's the **crux**, it **gates EM/light** (light rides
  the brane), and it's **algebra-tractable / NOT sim-gated** (the same gauntlet produced pathA_25's real no-gos). **Progress (2026-06-26):** `pathA_35` directive BUILT + gauntleted; **G0 freeze DONE+VERIFIED** = `T0_SHEAR_FROZEN(d9520d3819c3)` +
  `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)` (tri-reviewed; a DOF-count hole remediated & re-verified on a fresh clean agent). **⭐ Gate L
  EXECUTED + REMEDIATED = `FAIL_COUPLE_STRESS_NOGO`** (the predicted §2.6 four-way no-go): config A (massless-`Pⁱ`) fails a-ii (3 hidden
  `P` spin waves) + a-iii (C5 longitudinal zero mode) + b (gyrostat-unbounded, transverse minor `k²(Γ_P μ_R k² − λ_Pu²)<0`); the
  slaved-rigid escape clears a-ii/b but **still dies on C5** (`φ` absent; `φ=u_w` collides with the `u_w` gap; an independent `φ` = fresh
  G0). Traction PASSES (`ARROWS_SUPPLY_TRACTION`). **⭐ SUPERSEDED (2026-07-03) by the material-state pivot + pathA_36 — the Gate-L C5 leg became the crux.**
  The **material-state closure** (brane = ordered PHASE of the medium, order field `χ_B`; from `notes/brane_bulk_handoff.md`) passed
  **rung W** (GLM reframe `notes/rung_W_reframe.md`: a double-well `χ_B` wall is stable AND light-permitting), collapsing the crux to
  **C5**. **`pathA_36`** (directive gauntleted → SOUND; dual-engine SymPy+Mathematica; full Dirac–Bergmann) tested whether the `χ_B`
  PHASE supplies the C5 `φ`. **RESULT = `FAIL_CAUCHY_STRAY_LONGITUDINAL` / `BY_TUNING`, `ENGINE_AGREE`; tri-reviewed** (arbiter
  reproduced; Codex + adversarial `EARNED`; fidelity `FIDELITY_ISSUES` = the decisive sign was asserted-not-derived in-script, physics
  triply-confirmed). **⭐ ACHIEVEMENT: the medium DEMONSTRABLY carries LIGHT — 2 transverse photons at `c_γ²=μ_R/ρ_br` (first time).**
  C5: a *stable* order-parameter phase has the WRONG-SIGN gradient stiffness (that positivity = brane stability); Maxwell needs the
  opposite ("electric") sign; the only escape (negative stiffness) = a Lifshitz instability = **the SAME wall as pathA_25**. So θ-as-φ
  is a **provable vacuum no-go**. **⭐⭐ NEXT ACTION (post-compact) = the C5-WITH-A-THROAT gate.** THE REFRAME (do NOT bolt on an
  abstract gauge field): the longitudinal "third wave" is NOT a spurious photon — it is the `c_s` density/gravity mode + (static, around
  a charge) the **ELECTRIC FIELD = brane deformation around a puncture** (`conceptual_foundation.md` §4). pathA_36 analyzed VACUUM (no
  throats) — the one mode whose character comes entirely from charge, with no charge present → **CHARGE (the throat) is the missing
  ingredient, not a gauge field.** Build a C5-with-a-throat directive (pipeline: directive → Codex design-review → GLM → dual-engine →
  tri-review): does the electric field + Gauss's law EMERGE from the throat–brane coupling; is the longitudinal sector the `c_s` gravity
  channel (two speeds, not 3-modes-of-one-field); any longitudinal RADIATION coupling to charge (falsifiable, cf pathA_29). **If it
  works → the forces are CONNECTED, not just coexisting = a massive win.** This SUPERSEDES the pathA_36 vacuum-gauge remediation (vacuum
  result = documented waypoint). pathA_35 Gate-L + pathA_36 deliverables UNCOMMITTED. Read first: `docs/conceptual_foundation.md` §3 v6
  + §2 v5 + §4; `notes/brane_bulk_handoff.md`; `notes/rung_W_reframe.md`; `reports/pathA_36_c5_phase_potential.md`.
- **⭐⭐ SUPERSEDED 2026-07-03b — the C5-with-a-throat gate (pathA_37) was BUILT as an electric-FLOW gate, then the FLOW framing was
  RETIRED (charge is NOT a brane flow).** pathA_37 v1–v3 (deformation) → GLM energy-scaling no-go → v4 flow → v5 counterflow (all
  Codex-SOUND), then a 2nd GLM pass + a user-driven conceptual pass established: **charge = the interaction of two throats' 4D BODIES
  beyond the mouth** (deformation→wrong energy `1/r⁶`; one-fluid flow→charge=mass; counterflow→needs a 2nd component the `T=0`
  condensate lacks in flat 3D — but the counterflow no-go is a **flat-3D artifact**: the throat is 4D). **Gravity correction:** strictly
  ONE-WAY (medium de-structures in-throat, no return through it; return = global `S_leakage`). **Four-way taxonomy:** gravity=one-way
  drain, light=brane shear, magnetism=brane swirl, charge=4D throat-body interaction. **⭐ NEXT ACTION = the charge crux is RE-HOMED to
  the THROAT/BULK sector** (Gate-T, with mass/geon/`r_e`/annihilation): does a throat's 4D body mediate a **brane-localized `1/r²`**
  interaction (naive 4D-bulk → wrong `1/R³`; needs the SAME localization that gave pathA_29 gravity its `1/r²`) with `±w` = sign.
  **pathA_37 flow-gate framing = RETIRED — do not execute as-is.** Canonical picture: `docs/conceptual_foundation.md` §3/§4 ⭐ v7.
  (Deferred future note, do NOT compute: same localization may turn `1/r²`→`1/r` at galaxy scale — a rotation-curve/MOND-like hint.)
- **⭐⭐ pathA_38 EXECUTED (2026-07-03) = `THROAT_ELECTRIC_LOCALIZED_COULOMB` — EARNED PASS (tri-reviewed).** The charge crux is
  RESOLVED at the reduced/spec level: the throat's 4D body mediates a **brane-localized `1/r²` Coulomb** via the **gapless
  transverse-embedding / orientation-lock Goldstone `h`** (wall displacement into `w`, `±w` arrows locked to the normal; NOT the free
  internal director = gapped → Yukawa, the pathA_36/25 wall). **DERIVED** = {`p=2` localization (normalizable `sech²(w/ℓ)` zero mode,
  `N₀=8/(3ℓ)`; delocalized→`p=3` counterfactual computed and wired to the classifier), like-repel/unlike-attract SIGN STRUCTURE
  (`G₀>0`, general theorem), charge≠mass by `w→−w` parity}; **CALIBRATED** = {charge magnitude `Q_E`; nonzero-monopole modeling};
  **SIM-DEFERRED** = {`FAIL_OPERATOR_PARITY_MIXING`, `FAIL_SOURCE_NOT_COMPACT`, nonlinear throat interior}. **Gauntlet:** directive
  `SOUND_AS_IS` (Codex r1→r4 + GLM tertiary, all concerns folded; commit `2c5fe4f7`) → dual-engine (SymPy+Mathematica, `ENGINE_AGREE`).
  ⚠️ **The FIRST execution was a pass-by-construction RIG** (hardcoded `1/(4πR)`, gate structurally unable to emit a FAIL) — CAUGHT by
  the fresh **adversarial-with-ablation** leg, REMEDIATED to genuinely derive localization from the transverse-mode projection (à la
  pathA_29's `g''+(2/R)g'+…=0`) + empirically able-to-fail (9/9 forced-bad inputs → correct `FAIL_*`), then re-tri-reviewed CLEAN
  (`ADVERSARIAL_SOUND` + `FIDELITY_CLEAN`, two minor stylizations closed). Deliverables: `tools/pathA_38_throat_body_electric_{sympy.py,.wl}`
  + `reports/pathA_38_*`. **⭐ RESULT: electric charge's `1/r²` FORCE LAW + sign + charge≠mass now share gravity's localization
  PRINCIPLE — 3 of 4 sectors EARNED (gravity pathA_29, light pathA_36, charge pathA_38). NEXT = magnetism (close Maxwell, largely
  downstream of pathA_38) + the `λγ=c_γ/c_s=1` cross-sector consistency knit.** Canonical: `docs/conceptual_foundation.md` §3/§4 ⭐ v7.
- **Gate-6 Tier-A is DEFERRED as an optional cap** — it would likely just re-confirm that the `Z0_ret,Z1_ret` selector / `L/a`
  self-selection need the sim-deferred full solve (pathA_29 `Z_is_premise`; self-selection "requires dynamics"); its 4D field solve
  is sim-deferred regardless. The full road to a simulation-ready PDE (gravity arm done; brane → EM/light → throat → integration
  remaining) = **`docs/development_plan.md`**.
- **✅ Gate-1–3 dimensional checks — RETROFITTED (2026-06-26, parallel session).** The previously **VACUOUS** typed-tuple ledgers
  in `pathA_30/31/32` are now real-expression + able-to-fail + dual-engine (per `directives/pathA_dimcheck_retrofit_gates1to3.md`,
  the real `pathA_18` harness) and **everything checks out** — no OOM surprises. (Orchestrator confirm-read optional.)
- **⭐⭐ SCOPE EXPANDED (2026-06-26):** the full nonlinear **simulation is OUT OF REACH** (solver tractability, not hardware), so the
  goal is now to bring **EVERY sector (gravity / EM-light / brane / throat) to "simulation-ready"** via the algebra/calibrate gauntlet
  and **defer the sim** as future work. The win = a complete, internally-consistent, calibrated, *simulation-ready* PDE — it
  **completes the spec, does NOT prove the theory** (sim-dependent Qs stay posed; a no-go is still possible). **Master work breakdown
  = `docs/development_plan.md`** (memory `project-simulation-deferred-complete-pde-strategy`). The Gate-5/6 gravity ladder below is now
  one *arm* of that broader plan; Gate 6's full 4D field solve is sim-deferred (its algebra-tractable reduced closure is still in-scope).

Process lesson banked (memory `project-model-mechanics-corrections` + `feedback-offload-review-gauntlet`): every compute gate runs the
full Codex→GLM→Codex directive gauntlet + dual-engine + tri-review, with the **review rounds offloaded to a gauntlet-runner agent** so
the orchestrator context survives. pathA_30 took Codex(8)→reconfirm→GLM(Bogoliubov `k⁴`+5)→Codex post-GLM(2)→close-out before a clean
execution.

**Model-mechanics reminders** (memory `project-model-mechanics-corrections`): nothing is static; **three speeds** — `c_s` = speed
gravitational changes propagate (∝ρ²); `v_r` = field *strength*, not a speed; `c_γ` = light. Gravity = the **flow between drains**;
changes propagate at `c_s`, and uniform motion tracks the **current** position → **no aberration**. Throat-soliton has **no sloshing**
(`J_w=0` expected; AC→DC retired); gravity is a separate background de-structuring drain.

**COMMIT STATE (2026-06-26):** pathA_28 (`8cf6f1f1`), pathA_29 (`145c8426`), pathA_30 Gate-1 (`f460fc63`), pathA_31 Gate-2
(`765db5f0`), the prior docs-sync (`810e01f7`), **pathA_32 Gate-3 (`6711167a`)**, and **pathA_33 Gate-4 (`75ad0d0d`)** are all
committed. This **documentation sync** — relabeling the superseded verdict-count + EM / brane-existence material below as *history*,
pinning the Gate-4 hash, and disambiguating the two clashing "Gate 3 / Gate 4" namings — is **uncommitted, pending the user's go.**
**Commit only when the user asks; stage explicit paths.**

**Process discipline (unchanged):** Codex codes / Claude reviews; **dual-engine** (Mathematica: Codex needs `--sandbox
danger-full-access` — workspace-write CAN'T run it; OR the orchestrator runs `math` directly as arbiter); **review ordering** = iterate
Codex to GREEN → one GLM pass → fold → Codex to green; crux execution prompts get the full gauntlet; reports-only; `codex exec … -c
model_reasoning_effort=xhigh` backgrounded, never wrap the session in `timeout`.

---

## Earlier track — the gravity-sector verdict count (`pathA_22b`; the context the push grew out of)

> **⚠️ Naming caution — "Gate 3 / Gate 4" are OVERLOADED across two programs.** In this section, `pathA_22b` **Gate 3 / Gate 4** are
> the *verdict-count* gates (compute `χ_Q`; classify the `g_mhat²/g_G` ratio). The ⭐⭐ LATEST section above uses **Gate 3 / Gate 4**
> for the *moving-throat ladder* (ℓ=2 isotropy `pathA_32`; the `54/5` quadrupole `pathA_33`). Different programs — the dates
> disambiguate (`pathA_22b` Gate 4 = 2026-06-22; moving-throat Gate 4 = 2026-06-26).

Before the moving-throat PDE push, the gravity sector was scored by a single **verdict-count** equation
`P0·χ_Q·g_mhat²·λγ⁵/g_G = 54/5` under calibrate-predict discipline (every value DERIVED or a declared calibration gap — no silent
knobs). Terminal state of that count: `P0` extracted; **`χ_Q ≈ 0.712` COMPUTED** (`pathA_22b` Gate 3, outgoing-DtN, even-`nw`
converged); `λγ` a genuine model gap (`BETA_GENUINE_GAP`); and — **`pathA_22b` Gate 4 (2026-06-22) → `GENUINE_BLOCKED`** — the ratio
`g_mhat²/g_G` is **not derivable** from the action (no target-blind source-map kernel; `α_J` doesn't cancel; all 22 `m̂` sites in
`pde.tex` are target-facing, dual-reviewed). So the model does **not** derive its own gravity coupling: `g_G` is calibrated on
Newton's `G`, and the **`54/5` quadrupole is an ABSORBED calibration anchor, not a prediction**; the count closes only via the
**EM-sector anchor** (pins `λγ`) — load-bearing. The falsifiable payoff is the **held-out surplus** (g−2, 5PN, ringdown,
multi-defect), riding the shared *derived* `χ_Q` + `P0/D0` bundle + `c_s`. **This verdict-count framing is what motivated the EM
re-founding (history below) and, later, the more fundamental moving-throat PDE push (⭐⭐ LATEST, above).**

> **⚠️ Open reconciliation (do NOT silently merge the two `χ_Q`s):** the moving-throat Gate 4 (`pathA_33`) independently **derived
> `χ_Q = 1`** in the clean outgoing-DtN Hankel context (incoming = −1), whereas `pathA_22b` Gate 3 **computed `χ_Q ≈ 0.712`**
> numerically in the older minimal-combination context. These are different definitions/computations bearing the same name, not a
> contradiction to paper over — reconciling them is a tracked documentation item (alongside the Gate-1–3 dimensional-check retrofit).

## Earlier track — EM re-founding → "why the brane exists" (history; routes CLOSED / PAUSED, 2026-06-23→24)

> **This was the active frontier in 2026-06-23→24; it has since handed off to the moving-throat PDE push (⭐⭐ LATEST, above).** It is
> recorded here as **history** — its falsification results are first-class and load-bearing for the current picture. `decisions/13` §0
> top block carries the same "this is the history that led here" framing. Conceptual home: `docs/conceptual_foundation.md` (v3) +
> `docs/medium_requirements_and_prior_art.md`. Live EM-track ledger: `software/stage1_solver/reports/pathA_25_STATUS.md`.

**Why it started.** Pinning `λγ` exposed that the EM sector had **drifted** from the single-medium concept (canonical EM = a
fundamental gauge field on a flat metric, *decoupled*; `reports/pathA_cgamma_of_rho_derivation.md`, Type-4). The re-founding
hypothesis: **our 3D space = an elastic brane; LIGHT = the brane's in-plane SHEAR (MacCullagh rotational elasticity), on the brane
not the bulk** → bulk stays a pure fluid → magnetism preserved.

**`pathA_23` (EM medium-native, Stages 0–3b) — DONE, CONDITIONAL.** Stage-2 **CRUX = `FAIL_UNSPECIFIED_SUBSTRUCTURE`**: the medium's
record does not determine the brane shear modulus `μ_br` (verified trilemma — `μ_br=0`→no light; Cauchy→stray longitudinal "second
photon"; MacCullagh curl-only→needs postulated gyrostats), so **brane-shear EM is not derivable medium-natively and `λγ` stays a free
input.** **User chose Path 1: postulate MacCullagh, CONDITIONAL.** Leak findings: the brane↔bulk interface traction is *generically
transverse* near a draining defect (Stage 1), bounded only under an unmotivated `ε_leak≪1` price (Stage 3), but flat light *free-slips*
and the leak is **curvature-localized / far-field-vanishing** (Stage 3b) — throat/defect leak relocated, not retired. Full picture =
`decisions/15` (§11 MacCullagh, §13 λγ, §14 costs, §15 crux, §16 leak, §17–18 brane/defect).

**`pathA_24` (little-arrows domain-wall brane) — ❌ FALSIFIED.** `T1_FAIL_NO_STABLE_WALL` (`2fa91886`, tri-reviewed GENUINE): a static
polar-vector wall on a **connected `S³` vacuum manifold (π₀=0)** spreads to infinite width and unwinds with zero barrier — no stable
wall. The three-way no-win (emergent-`w` / stable-wall / light-capable-core) demonstrated; independently corroborated by a 5-agent
prior-art survey (Davies–George–Volkas; ³He).

**`pathA_25` (GNLS polar-smectic candidate) — ❌ DENSITY ROUTE CLOSED.** G0 froze the structure (`SECOND_MEDIUM_DRIFT_AT_FREEZE`, 5
calibration inputs). Then **B4 = `FAIL_NOT_CODIM1`** (the smectic onset is a rank-2 triad, not a codim-1 lamella — the GNLS cubic
`U'''=15Kρ0²>0` lowers the triad) **and R/C = `RC_DENSITY_SMECTIC_LIGHT_NOGO`** (the only driver that opens a genuine codim-1 density
window, `Cpin`, pins `P` out-of-plane → `P_∥=0` → starves in-plane light). Both tri-reviewed, dual-engine, earned. **Conclusion: a
density modulation can form a brane only by killing light → the brane must be a *light-confining shear surface*, not a density
layering.** (Ledger: `pathA_25_STATUS.md`.)

**Pivot → throat-soliton / 4D-light + drain sector (PAUSED 2026-06-24).** The brane = a light-confining shear surface; light is a
**4D shear field** (3+1D face on the brane via a gapped `w`-mode → 2 polarizations; 4D at throats; trapped = mass = geon). **pathA_26
Derrick + open-stability = `THROAT_DRAIN_DESTABILIZED`, interpreted NOT-a-kill** (a conservative throat-soliton exists generically;
instability only at unphysical large drain = black-hole regime). The remaining EM-track step — the **drain-sector derivation** — is
**PAUSED** (memories `project-light-is-4d-throat-hypothesis`, `project-brane-existence-defect-structure`); the program's live effort
moved to the gravity-sector moving-throat PDE (⭐⭐ LATEST). Deferred/parked: `pathA_23` Stages 4–6; why fluid flows *into* the mouth
but leaks *back* into the brane (`decisions/15` §9).

**Methodology (locked, applies to all tracks):** Codex codes + runs + dual-engine; **AI prose never establishes a math fact** —
orchestrator arbiter re-run + transliteration-fidelity audit + adversarial review on clean agents; user gate per gate.

## Map — what you want → which doc holds it

| You want… | Look here |
|---|---|
| **The conceptual vision — what the medium / brane / 4 sectors / defect physically ARE (read first)** | `docs/conceptual_foundation.md` |
| **⭐⭐ The master development plan — ALL sectors → "simulation-ready" (the full scope; read for "what's left")** | `docs/development_plan.md` |
| **⭐ How we work — the dev pipeline / review gauntlet (read before running ANY gate)** | `docs/development_pipeline.md` |
| **⭐ The gravity arm — moving-throat PDE ~6-gate master checklist** | `research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md` |
| **⭐ Moving-throat ladder gate artifacts (Gates 1–4)** | `software/stage1_solver/{directives,tools,reports}/pathA_30..33_*` |
| **⭐ Gate-5 cross-ℓ framing (ℓ=0/1 return ↔ ℓ=2 quadrupole)** | `research/pde_ledger/notes/stages/moving_throat_pde_open_item_reconciliation.md` |
| Full current state + resume-after-`/compact` pointer | `software/stage1_solver/decisions/13_emergent_constants_derivation.md` §0 |
| Every value classified DERIVED / INPUT / gap / benchmark | `software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md` (Gate-4 closeout = §6) |
| The defect-regime + held-out-surplus roadmap | `docs/defect_interaction_map.md` (CURRENT STATUS banner at top) |
| EM-track history — requirements list + prior-art survey + GNLS polar-smectic candidate (routes closed) | `docs/medium_requirements_and_prior_art.md`; ledger `software/stage1_solver/reports/pathA_25_STATUS.md` |
| EM-track history — physical picture + MacCullagh template + leak findings | `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md`; directive `directives/pathA_23_em_medium_native.md` |
| Per-gate derivation detail (`pathA_22b` *verdict-count* Gates 0–4 — NOT the moving-throat ladder) | `software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md` |
| The directive that ran `pathA_22b` Gate 4 (the `g_mhat²/g_G` ratio = `GENUINE_BLOCKED`) | `software/stage1_solver/directives/pathA_22b_gate4_ratio_or_blocked.md` |
| The calibrate-predict methodology | `software/stage1_solver/decisions/09_calibrate_predict_methodology.md`; `docs/pathA_preregistration.md` |

## The `pathA_22b` verdict equation (reference — the *earlier* verdict-count framing, not the moving-throat ladder)

```
P0 · χ_Q · g_mhat² · λγ⁵ / g_G  =  54/5
 ✓     ✓     gap1     gap2  cal-on-G        (✓ = derived; gap1 g_mhat absorbs 54/5; gap2 λγ ← EM anchor)
G = (a·c_s²/m_GNLS)·g_G ,  m̂0 = (c_s/(a²·√m_GNLS))·g_mhat ,  c = λγ·c_s
```

> Here `χ_Q ≈ 0.712` (the `pathA_22b` Gate-3 number). The moving-throat ladder's Gate 4 (`pathA_33`) later derived **`χ_Q = 1`** in
> the outgoing-DtN Hankel context — a different computation of the same-named factor (see the open-reconciliation note above), and the
> moving-throat ladder reaches the same `54/5` via the earned/calibrated split `54/5 = 2·27/5`, not via this verdict-count product.
