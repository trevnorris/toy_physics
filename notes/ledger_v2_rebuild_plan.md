# Ledger v2 Rebuild — multi-phase plan & post-`/compact` handoff (2026-07-05)

> **⭐ READ THIS FIRST after `/compact`.** It is the complete plan + the evidence that led to it. Front door = `STATUS.md` ▶ RESUME
> HERE (thin pointer to here). Memory = [[project-pde-ledger-fold-in-plan]] (same content, condensed) + [[project-brane-existence-defect-structure]].

---

## 0. TL;DR — where we are and what's next

- **All 4 force-sectors EARNED + the consistency knit DONE** (cone-lock `pathA_40`, NG5 `pathA_41`, charge-coupled scalar `pathA_42`).
  The four sectors live in `software/stage1_solver/` as `pathA_*` gates; the gravity ladder is `pathA_21c/28/29/30–34` + the separate PN
  corpus `research/4d_*pn*`.
- **DECISION (user, 2026-07-05): REBUILD the central ledger FROM SCRATCH.** Do **not** fold the earned work into the old 253-stage
  `research/pde_ledger/`. Build a **fresh, sector-organized ledger at `research/pde_ledger_v2/`**, reuse the old ledger's proven
  *machinery*, carry over only useful *content*; the old ledger becomes a reference **quarry** and is later overwritten into git history.
- **Sequence (3 phases):** **A** finish gravity in `stage1_solver` → **B** build `research/pde_ledger_v2/` → **C** redteam it.
- **✅ Phase A, step A1 DONE (2026-07-06) = `pathA_43` `DENSITY_PORT_HOSTED` (EARNED).** The ℓ=2 quadrupole radiative-port numerator
  `N0` is density-native — a genuine `(q2` wall,`Φ2` bulk-density`)` two-port, provably vector-free, `[N0]=L⁻¹M`, `a⁻⁵`, outgoing sign
  `+i z⁵/27`, `χ_Q=1`; **the EM `A_w`/`U,W` scaffold retires + the diagnostic sliver closes.** Structure reduced-closure-earned;
  magnitude + `G/2/5/54/5` CALIBRATED/SIM-DEFERRED. Full gauntlet + tri-review CAUGHT A RIG → remediated → re-tri-review + hardening all
  CLEAN. Artifacts `{directives,tools,reports}/pathA_43_density_quadrupole_port*`.
- **✅ Phase A, step A2 DONE (2026-07-06) = `pathA_21c` dual-engine SymPy companion (`f236ec08`).** Added
  `tools/pathA_21c_force_from_noether_stress_tensor_sympy.py` — an independent SymPy engine peering the existing `.wl` (the report had
  long *claimed* a Python engine that never existed). Reproduces all 9 dim + 9 alg checks (drain→1/r² `FORCE_ATTRACTIVE_DERIVED`),
  derives Ω_2=4π, Ω_3=2π², `⟨n_i n_j⟩=δ_ij/d` from explicit S²/S³ surface integrals + the attractive sign from primitive orientation
  conventions, no runtime consumption, mandatory able-to-fail mutation probe (7 corruptions incl. a derived-second-moment corruption).
  **Lighter gauntlet** (no new verdict): Codex design-review → GLM-5.2 tertiary → Codex re-green (all `DIRECTIVE_CLEAN`) → dual-engine
  both exit 0 → arbiter re-run + fidelity audit `FIDELITY_CLEAN`.
- **✅ Phase A, step A3 DONE (2026-07-06) → PHASE A COMPLETE = 2.5PN match-back.** Consolidation note
  `notes/pathA_43_2_5pn_matchback.md` (pathA_43's density-mode ℓ=2 port reproduces the 2.5PN Burke–Thorne `Γ̄₅=2G/(5c⁵)` +
  `K̄₄=4K̄₂²/K̄₀` that `research/4d_2_5pn` left as its single open item) + a dedicated dual-engine verification artifact
  `tools/pathA_2_5pn_matchback_{sympy.py,.wl}` (independently checks INV1–INV5 — BOTH Γ̄₅ forms → BT incl. the corpus's own
  `9K̄₂^{5/2}/K̄₀^{3/2}`, the cross-form agreement, the moment invariant, + literal anchors defeating a coherent-rescale rig — with an
  11-mutation able-to-fail probe + asserted caught-by matrix). Consistency over CALIBRATED moments, NOT a first-principles
  `Γ̄₅`/`G` derivation (`G=GENUINE_BLOCKED`; full 1PN→4PN = SIM-DEFERRED). Note gauntlet (Codex→GLM→Codex `NOTE_CLEAN`); script spec
  gauntlet (Codex design-review CAUGHT a coherent-rescale blind spot → INV5 anchors; GLM `SOUND`; Codex `DIRECTIVE_CLEAN`) →
  dual-engine exit 0 → **full tri-review: arbiter re-run + `FIDELITY_CLEAN` + `ADVERSARIAL_CLEAN`**. **▶ NEXT = Phase B (build
  `research/pde_ledger_v2/`).**
- **Estimated v2 size ≈ 40–80 stages** (vs 253) — see §6.

---

## 1. Why REBUILD, not fold-in (the evidence — all from this session's fan-outs)

**The old `research/pde_ledger/` was an aimless AI search, not a focused corpus.** Three-agent anatomy fan-out found:
- **253 = ~8 CORE stages + ~225 SEARCH tail + ~52 non-gravity.** The banked gravity results live in ~8 stages
  (`001` field-eq skeleton, `002` (a,L)+isotropy, `023/024` 54/5 + O(3) isotropy, `104/105` fingerprint `1/9,4/81,1/27`→`χ_Q=1`,
  `100/103/106` closure `m̂₀²χ_Q N_Q=1`). The other **~225 stages are a non-terminating branch-selection/optimization SEARCH** toward the
  nonlinear Gate-6 branch realization (support/mouth funnels; Part VI = simplex cardinality 3→4→5 optimizers) + robustness + scope-drift.
- **~52 stages are NOT gravity:** `004–020` = an imported projected-Maxwell EM block (17); `219–253` = a same-charge/barrier/scalar-photon
  EM-charge probe (35). Reusable gravity = `001–003 + 021–218`.
- **The 253-ledger was already SUPERSEDED by `pathA_30–34`** (Gates 1–5) — the same job in ~5 tri-reviewed gates. Ledger stages = derivation
  NOTES; the pathA gates = the EARNED versions.
- **The 1PN→4PN ladder is a SEPARATE audited corpus** (`research/4d_1pn_full`, `4d_2pn`, `4d_3pn`, `4d_2_5pn`, `4d_4pn`,
  `1pn_orbital_dynamics`). The 253 only MATCH-BACK to it ("do NOT re-derive"). Reaching 4PN did NOT require the 253 stages.

**The EM "contradiction" DISSOLVES (2-agent diagnostic).** The old ledger's `S_EM[A]`/`A_w` Maxwell field is **REPLACEABLE bookkeeping,
NOT load-bearing** for the gravity quadrupole:
- The `54/5`/`ω⁵` quadrupole is a **scalar/density (`c_s`) spherical-Hankel DtN** (`Ŷ₂^out=1+a²ω²/9c_s²+4a⁴ω⁴/81c_s⁴+i·a⁵ω⁵/27c_s⁵`).
  `A_w`/`Maxwell` appear ZERO times in the stages that COMPUTE it (104/105/194/195); the `A_w` label in stage 023 is a superseded placeholder.
- `54/5 = 2/5 × 27`: the `2/5` = GR Burke–Thorne `2G/(5c⁵)`; the `27` = the `c_s` density-Hankel tail. So `c⁵` (light) is a **GR-matching / λγ
  units bridge** (`P₀ ∝ c_s⁵/c⁵ = 1/λγ⁵`), NOT EM propagation. Delete-`S_EM[A]` counterfactual: the EM moment is an additive term that cancels
  from the ℓ=2 surface; the quadrupole physics is untouched.
- ⚠ ONE honest sliver: the outgoing-port numerator `N_0` is, as literally written, hosted on the `A_w`/`U,W` port coordinates. Both agents
  judge this an ATTACHMENT convenience (the note's own §8 re-hosts on the scalar/density lane), not a necessity. **Closing it = derive the ℓ=2
  port NATIVELY on the density mode = exactly Phase A step A1.**

**Bottom line:** gravity radiation rides the medium's OWN density ripple at `c_s` (as the physical frame demands); the EM scaffold is an
older, heavier bookkeeping device we simply don't carry. Our 4 earned sectors + knit are a focused conclusion; tacking them onto the 253-stage
sprawl would bury the signal. So we build fresh.

---

## 2. "From scratch" = fresh CONTENT + ORGANIZATION, **REUSED MACHINERY**

Keep the old ledger's proven infrastructure (this is what stops v2 ballooning back to 253):
- the **stage triplet**: note `notes/stages/moving_throat_pde_stage{NNN}_{slug}_sympy_audit.md` + script pair
  `scripts/{stem}.py` + `mathematica/{stem→_mathematica_audit}.wl` (print-only, assert-zero, exit-code) + paper card
  `paper/stages/stage_{NNN}.tex`;
- the **`MATHEMATICA_MIRROR_POLICY`** (independent `.wl`, not a payload mirror);
- the **LaTeX build** (`pde_ledger.tex` → `document_setup`/`macros` → `main_parts` → `parts/*` → `appendices/stage_appendix_partNN` → `stages/*`);
- the **provenance/coverage index + runners** (`STAGE_PROVENANCE_INDEX.md`, `STAGE_VERIFICATION_COVERAGE.md`, `run_all_audits.sh`, `run_one_audit.sh`).

**⚠ The reshape our `stage1_solver` scripts need** (they violate the ledger contract two ways): (1) they use `argparse --compare` + write
JSON/YAML payloads; (2) their `.wl` is a PAYLOAD MIRROR of the SymPy (exactly what `MATHEMATICA_MIRROR_POLICY` rejects). Folding each gate in
requires: strip the compare-harness/file-writing, make each engine standalone print-only/assert-zero/exit-nonzero, and **re-author each `.wl`
as a genuinely independent route.**

---

## 3. The v2 ledger structure (sector-organized; mirrors `model_field_equations.tex`)

| Part | Sector | Carries over from stage1_solver / old ledger |
|---|---|---|
| **I** | The medium | GNLS action + EOS (`P=Kρ⁵`, `c_s²=5Kρ⁴/m`) + brane/bulk two-phase structure + order field (`pathA_35` G0 freeze) |
| **II** | **Gravity** | drain→1/r² force (`pathA_21c`), return/1r²-survives-slab (`pathA_29`), monopole/dipole radiation (`pathA_28`), reduced-closure Gates 1–5 (`pathA_30–34`), the NEW density-mode ℓ=2 port (A1), 54/5 calibration (`GENUINE_BLOCKED`), PN match-CITE (`research/4d_*pn*`), Gate-6 branch realization = CITED sim-deferred |
| **III** | Light | brane shear, `c_γ²=μ_R/ρ_br`, 2 transverse photons (`pathA_35` gateL + `pathA_36`; FAIL top-lines carry the earned photon as content) |
| **IV** | Charge | throat-body Coulomb 1/r², ±w sign, gapless embedding Goldstone `h` (`pathA_38`) |
| **V** | Magnetism | the moving 4D throat-body interaction, 4 stages (`pathA_39` s0/1, s2, s3, s4); EM field = transverse-vector + charge-coupled scalar, NOT exact Maxwell |
| **VI** | The knit | cone-lock `λγ`/`c_E=c_γ` (`pathA_40`), one-medium drift NG5 (`pathA_41`), charge-coupled scalar map (`pathA_42`) |
| **VII** | Integration | unified frozen parent action + BCs + the calibration map (every postulate/constant labeled) + whole-system dim-check + **the simulation hand-off spec** (every sim-dependent quantity + the exact equation + confirming measurement = the literal deliverable) + the χ_Q reconciliation |

**Carry over (quarry):** the ~8 core gravity stages reframed onto the density mode; Gate-6 reconnaissance as CITED sim-deferred; PN
match-back (cite `research/4d_*pn*`); all the machinery (§2). **Drop (quarry for stray lemmas only):** the EM import `004–020`, the ~225-stage
search sprawl, the charge-probe drift `219–253`.

---

## 4. Phase A — Finish gravity in `software/stage1_solver/` (✅ COMPLETE 2026-07-06)

**✅ Phase A is COMPLETE — A1 (`pathA_43` density port), A2 (`pathA_21c` SymPy companion), A3 (2.5PN match-back) all DONE.**
▶ NEXT = Phase B (§5). The sub-step detail below is the (now-executed) plan, kept for provenance.

Gravity is ~90% earned (see §3 Part II). Remaining:

**A1 — the density-mode ℓ=2 radiative-port gate (the ONE genuinely-new derivation).** Derive the ℓ=2 quadrupole SOURCE/port on the
density/`c_s` mode (extending `pathA_29`'s ℓ=0/1 density-return machinery up to ℓ=2), reproducing the `54/5`, the fingerprint
`1/9,4/81,1/27`, and `χ_Q=1` **without** the `A_w` host. Retires the EM scaffold + closes the diagnostic sliver.
- **Genuinely able-to-fail:** if the density lane CANNOT source the ℓ=2 port with the right sign/normalization, the vector channel WAS
  load-bearing after all — a real finding that reopens the EM question. The gate must be able to emit that.
- **Setup sub-questions to resolve first (Claude↔Codex scope, before the directive):**
  1. What is the density/BdG-lane source coupling of the moving throat to the **ℓ=2** outgoing wave (extend `pathA_29`'s ℓ=0/1 result)?
  2. What is the density-lane analog of old-ledger stage 023's `U,W` 2-port that produces the numerator `N_0`?
  3. The able-to-fail condition: what output shows the density lane cannot host the ℓ=2 port cleanly (wrong sign/normalization/gauge-obstruction)?
  4. Confirm the reframe reproduces the SAME `54/5` / fingerprint / `χ_Q=1`.
- Then the full gauntlet: directive → Codex design-review xhigh → GLM-5.2 tertiary → dual-engine → tri-review (adversarial-with-ablation) → user gate.

**A2 — `pathA_21c` SymPy companion.** `pathA_21c` (`FORCE_ATTRACTIVE_DERIVED`, the drain→1/r² force) ships `.wl`-only
(`pathA_21c_force_from_noether_stress_tensor_crosscheck.wl`) — no `_sympy.py`. Add the SymPy side for dual-engine compliance, OR fold the
force law into a gravity consolidation.

**A3 — 2.5PN match-back ✅ DONE (2026-07-06) → PHASE A COMPLETE.** Delivered as `notes/pathA_43_2_5pn_matchback.md` + the dedicated
dual-engine verification artifact `tools/pathA_2_5pn_matchback_{sympy.py,.wl}` (full summary + tri-review in the §0 A3 bullet). It grew
one step beyond a pure note: rather than lean on pathA_43's closure overlay + reviewer arithmetic, it got its own self-contained
able-to-fail artifact checking INV1–INV5 (both Γ̄₅ forms → BT incl. the corpus's `9K̄₂^{5/2}/K̄₀^{3/2}`, cross-form agreement, the moment
invariant, coherent-rescale-defeating anchors). The plan below is the (now-executed) reasoning. Verify the density-mode quadrupole
(`54/5`) agrees with the PN ladder's 2.5PN coefficient (the "cheapest decisive falsifier"). **⚠️ A3-PREP ALREADY DONE (during A1):** an agent
characterized the PN corpus — 1PN→4PN is calibrated/GR-matched controlled-reduction with `G` a genuine gap; the 2.5PN sector
(`research/4d_2_5pn`) is a CONDITIONAL theorem whose SINGLE open item is exactly `Γ̄₅=2G/(5c⁵)` (equivalently `Γ̄₅=9K̄₂^{5/2}/K̄₀^{3/2}`)
— which `pathA_43` already supplies at reduced-closure (its closure overlay computed `Γ̄₅−2G/(5c⁵)=0` and `K̄₄−4K̄₂²/K̄₀=0`, dual-engine,
in `reports/pathA_43_*`). So do NOT re-run the corpus-characterization agent. **▶ A3 FIRST ACTION = write the CONSOLIDATION NOTE** (not a
new derivation): record that `pathA_43`'s density-mode quadrupole reproduces the 2.5PN Burke-Thorne target `Γ̄₅=2G/(5c⁵)` + the
`K̄₄=4K̄₂²/K̄₀` consistency, citing `research/4d_2_5pn`'s open item; likely a short markdown note under `notes/` or a stage in v2 Part II
(decide location). Optionally re-verify the two invariants are present in `reports/pathA_43_*` before citing. **Honest scope:** magnitude
CALIBRATED (`G` genuine gap); a full first-principles 1PN→4PN re-derivation from the throat stays SIM-DEFERRED (Gate 6); A3 is the
reachable consistency check, not a from-scratch 4PN derivation.

---

## 5. Phase B — Build `research/pde_ledger_v2/` · Phase C — Redteam

**Phase B:**
- **B1** — Claude↔Codex design pass → the **v2 blueprint** (Part/stage map + explicit carry-over list + build order + stage-decomposition
  granularity) for USER APPROVAL before any building. (Foundational directive → gets a Codex design-review.)
- **B2** — stand up `research/pde_ledger_v2/` by copying the old ledger's build/mirror/provenance scaffolding as the starting infrastructure.
- **B3** — assemble the sectors as stages (each earned gate → stage(s), reshaped scripts per §2, dual-engine, provenance-registered).

**Phase C:** run the redteam on the built v2 ledger (the hardening/solidification pass — after it's built, per the user's "redteam last").

After v2 is complete + hardened: **override the old `research/pde_ledger/`** (it lives on in git history).

---

## 6. Honest scope, size estimate, guardrails

- **Size ≈ 40–80 stages vs 253.** Converged content ≈ ~30 distinct results (gravity ~9 incl. A1; light ~3; charge ~1–3; magnetism ~4;
  knit ~3; medium ~4; integration ~5); as stages that's ~40–80 depending on decomposition granularity to hit the "no algebra left unturned"
  bar. Much smaller because the ~225-stage non-converged search + ~52 EM/charge drift do NOT carry over.
- **Guardrail (unchanged):** this completes the **SPEC, not the proof.** Gate-6 nonlinear branch realization stays SIM-DEFERRED; a no-go is
  still possible (esp. A1's able-to-fail). Every sim-dependent quantity stays precisely posed in the Part-VII hand-off spec.
- **`research/pde_ledger/` is tracked by the MAIN repo** (the inner `.git` is a dead leftover); commit v2 via the main repo. Nothing this
  program did ever overwrote the old ledger (git-verified) — the fold is/was append-only; the rebuild is a NEW folder.

---

## 7. References
- Memory: [[project-pde-ledger-fold-in-plan]] (condensed plan + all diagnostic detail), [[project-brane-existence-defect-structure]]
  (the 4-sector ladder), [[project-calibrated-pde-goal]] (the end goal), [[project-pn-gravity-ladder]] (the separate PN corpus),
  [[project-simulation-deferred-complete-pde-strategy]] (the sim-deferred guardrail), [[feedback-dual-engine-required]],
  [[feedback-claude-reviews-codex-codes]], [[feedback-directive-design-review]].
- Old-ledger step pattern (reverse-engineered this session): the stage triplet + `MATHEMATICA_MIRROR_POLICY` + the LaTeX build chain
  (`paper/pde_ledger.tex`, `main_parts.tex`, `appendices/stage_appendix_partNN.tex`) + `notes/STAGE_PROVENANCE_INDEX.md`.
- Gravity artifacts: `software/stage1_solver/{tools,reports}/pathA_{21c,28,29,30,31,32,33,34}_*`. PN corpus: `research/4d_*pn*`,
  `research/1pn_orbital_dynamics`. Old ledger (quarry): `research/pde_ledger/` (core stages `001/002/023/024/100/103/104/105/106`).
- Superseded prior handoff: `notes/consistency_knit_handoff.md` (the knit is done; that note is history).
- The model's consolidated field equations (skeleton for the v2 structure): `model_field_equations.tex` (untracked, user's Overleaf copy).
