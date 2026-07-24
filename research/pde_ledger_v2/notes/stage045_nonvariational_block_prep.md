# Stage 045 (VII-2b) — the NON-VARIATIONAL drain/return block + BCs + force partition: PREP (DRAFT, not the directive)

> **⚠ STATUS: prep/anchor, NOT a ratified directive.** Authored 2026-07-24 right after stage 044 (VII-2a) COMPLETE
> (`087565b0` deliverables + `2af3169b` tracker sync). The **directive** for 045 is authored AFTER a fresh read of the 045
> sources (below) + AFTER the drain-placement USER MINI-GATE. This note is the resume anchor + the mini-gate framing. HEAD at
> authoring ≈ `2af3169b`. ⭐ **READ-FIRST for 045**, then `RESUME_ROADMAP.md` §7 + the ratified `part7_integration_atomic_split.md`.

## Where 045 sits
Part VII is IN PROGRESS: **043 (VII-1, count-as-range) + 044 (VII-2a, conservative parent action) ✅ DONE** (2 of 7). 045 is
**VII-2b = the non-variational complement of 044's conservative action** — the SECOND half of the spine. 044 assembled the
conservative `S_assembled` and PROVED it source-free (contains no drain); 045 builds the explicit **non-variational** drain /
return / BC / force-partition machinery that sits *outside* the action. Class: SYNTHESIS, dual-engine (on the global closure
identities + dimensional homogeneity) + prose (on the BC taxonomy). Foundational → FULL new-derivation gauntlet + **GLM-5.2
tertiary**. ⚠ The G0 card is DRAFT-v0 → 045 assembles at the **structure / completeness-floor level only** (no committed closed
action; the shared throat solve stays SIM-DEFERRED).

## ⭐ THE DRAIN-PLACEMENT — DECIDED by the user (2026-07-24): the drain is the DYNAMICAL wall, NOT a frozen wall
Stage 044 named TWO candidate **non-variational** drain interfaces, proved `S_assembled` source-free, and left the selection to
045. **The user has now RULED OUT the frozen-wall drain** (2026-07-24, hard-flagged: *"the drain cannot be a frozen wall — several
instances where we tried to freeze the wall and it screwed up all the calculations"*; this is the standing "NOTHING is static"
principle — see memory `project_model_mechanics_corrections` §1/§1a). So the committed non-variational drain is:

- **✅ (a) Part-I `Γ_B` order-conversion drain — THE COMMITTED CHOICE** (DYNAMICAL wall): `∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n
  Γ_B`, `Γ_B = Γ_return − Γ_drain`, dim `T⁻¹`. An INTERNAL order-conversion at the throat (the DYNAMICAL wall converting
  ordered↔disordered material; pin P12: phase-conversion, NOT suction). The wall stays dynamical — no freeze.
- **❌ (b) G0-card `S_drain` ρ-mass-sink + remote return — RULED OUT** (frozen-wall regime; card §6 freezes `r_B` via `S_hold`
  and sinks `ρ`). This is the exact frozen-wall shortcut the user rules out. `S_drain = −Σ_i Γ_0 D_i` + `S_leakage = (Σ_iΓ_0)R_0`.
- **(c) both as regimes** — MOOT for the committed drain (its frozen-Σ half is (b), which is out). If a static reduction is ever
  wanted it must be DERIVED as the `ω→0` limit of the dynamical (a), never posed as an independent frozen solve.

**⭐ Consequence for the 045 build (do NOT just transcribe the G0 card §6 — it is frozen-wall-premised):** the drain SOURCE is
the dynamical `Γ_B` order-conversion (Part-I native), NOT the G0 `S_drain=−ΣΓ_0 D_i` ρ-sink. The G0 card §6 provides the
STRUCTURE to reuse — the finite-rank `I_ret`/`P_ret` return controllers, the `f`/`w` balances, the `𝔅`/mouth/collar/IR BCs, the
`F_var+F_flux+F_𝔅+F_rad` partition, the global-closure form — but the return/BC machinery must be **re-hosted on the dynamical
wall** (the return balances the `Γ_B` order-conversion, not a ρ-mass-sink). ⚠ Also verify at the read: 044's `S_hold`
(mid-surface GEOMETRIC pin, NOT a full freeze — 044 kept the wall inertial/dynamical) does not reintroduce a partial freeze that
trips the same failure mode; if it does, flag it. **At the 045 gate this is now a CONFIRM, not an open mini-gate** — briefly
confirm (a) with the user + the re-hosting scope, then build. (`τ_d` T-odd / "drain = arrow of time" still holds — the DYNAMICAL
order-conversion is itself irreversible via the drain, so no frozen sink is needed for the arrow.)

## What 045 BUILDS (the non-variational block — RE-HOSTED on the DYNAMICAL `Γ_B` drain, option (a))
- **The drain/return sources** — the DYNAMICAL `Γ_B` order-conversion (`Γ_B = Γ_return − Γ_drain`) as the non-variational RHS on
  the wall-order balance `∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B`. ⚠ NOT the G0 frozen-wall ρ-sink `S_drain = −Σ_i Γ_0 D_i` +
  `S_leakage` (that is RULED OUT; the G0 §6 sink/leakage/return STRUCTURE is a reference to adapt, not the drain to build).
- **The finite-rank return controllers `I_ret` / `P_ret`** — projectors enforcing EXACT global mass / momentum / energy closure
  despite the local sink.
- **The momentum `f_drain` / `f_leakage` and the energy `w_drain` / `w_leakage`** balances (non-Hamiltonian controller cells).
- **The `𝔅` sleeve-contact + mouth-Robin/Neumann + collar-Kirchhoff + two-sided-IR-return BCs** (the boundary taxonomy —
  READ the two_throat sim-spec's 16-branch `𝔅` BC).
- **The phase / moduli zero-mode quotients** (the return controller lives in the zero-mode complement).
- **The force partition `F_total = F_var + F_flux + F_𝔅 + F_rad`** — the total inter-defect force split into the variational
  piece (from 044's `S_assembled`, e.g. `−δS_hold/δX_i`), the flux piece, the `𝔅`-boundary piece, and the radiative piece.
- **The global closure identities `∫(S_drain + S_leakage) = ∫f = ∫w = 0`** — the exact-global-conservation certificates (the
  DUAL-ENGINE able-to-fail core: local drain, zero global change).

## OUT of 045 (do NOT do here)
- The shared THROAT SOLVE stays **SIM-DEFERRED** (the central reduction debt; 045 NAMES it, does not perform it).
- The calibration map + whole-system dimensional firewall = **046**; the permanent registers = **047**; the sim-handoff spec =
  **048**; the `χ_Q` characterization = **049**.
- No two-body / far-field magnitude solve (sim-deferred).

## Dual-engine vs prose (045's verification contract)
- **Dual-engine-checkable (SymPy + independent Mathematica):** the global closure identities (`∫(S_drain+S_leakage)=0`,
  `∫f=0`, `∫w=0` — exact global conservation given the local sink + remote return); the finite-rank controller structure
  (`I_ret`/`P_ret` projector idempotence / rank / zero-mode annihilation); the force-partition decomposition (the 4 pieces sum
  to `F_total`, each with the right symmetry/parity); the dimensional homogeneity of every non-variational term (units-restored
  `[L,T,M]`, able-to-fail); the zero-mode quotient (the controller lives in the moduli complement).
- **Prose (PRINTED, not asserted):** the BC taxonomy (`𝔅`/mouth/collar/IR) rationale; the non-variational bookkeeping argument
  (why these are weak balances on non-Hamiltonian cells, NOT extra variations); the drain-placement adjudication narrative (the
  gate outcome); the G0-card-is-DRAFT-v0 + throat-solve-deferred caveats.

## ⭐ Gravitomagnetism fold continues in 045 (per the earmark `dc617088`)
Flag the gravitomagnetic (velocity-dependent flow) terms in the RETURN-flow structure: frame-dragging = the boost of the
gravitoelectric drain/return flow — it lives HERE, in the moving/rotating-mass return circulation. The register/model_map
one-liner is still **047's** job; 045 flags it in prose. Keep the honest asymmetry (gravity matches GR incl. GM through 4PN;
EM departs Maxwell).

## Fresh-read sources for the 045 gate (read before authoring the directive)
1. **Part I's `Γ_B` order-balance (stages 006/007)** — ⭐ THE COMMITTED DRAIN SOURCE (option (a), dynamical wall): the
   `∂_t(χ_B n)+∇·(χ_B n u+J_χ)=n Γ_B` order-conversion + `Γ_return`/`Γ_drain`. Read FIRST — 045's drain is re-hosted here.
2. `software/em_charge_attribute/g0_closure_card_v0.md` **§6** (the `S_drain`/`S_leakage`/`f`/`w` + controllers `I_ret`/`P_ret`
   + the BC structure + the force partition) — a STRUCTURE reference ONLY (its drain is the RULED-OUT frozen-wall ρ-sink; reuse
   the controller/BC/force-partition FORM, re-hosted on the dynamical `Γ_B`). DRAFT-v0.
3. `docs/two_throat_simulation_handoff_spec.md` — READ-FIRST for the **`𝔅` BCs** (the 16-branch `𝔅` boundary-condition
   taxonomy + the R1→R4 ladder; partially-closed, NO closed action). The mouth-Robin/collar-Kirchhoff/IR-return BCs live here.
4. The **044 stage note** `notes/stages/ledger_stage044_parent_action_reconciliation.md` (what interfaces 044 NAMED + handed to
   045; the source-free `S_assembled`; the `F_var` piece from `−δS_hold/δX_i`; ⚠ verify `S_hold` is only a mid-surface pin,
   not a partial freeze — see the frozen-wall rule).
5. `research/pde_ledger_v2/notes/parameter_register.md` — the drain/return knobs already registered (`Γ_B` line 159; any
   `Γ_0`/`R_0`/controller rows); 045's new register rows start at **R98** (044 ended at R97).

## Process (full new-derivation gauntlet + GLM tertiary — 045 is foundational)
Author directive → **Codex→Grok→Codex bookend** → **drain-placement CONFIRM** (option (a) dynamical `Γ_B` is DECIDED — the
frozen-wall (b) is ruled out; briefly confirm the re-hosting scope + the `S_hold`-not-a-partial-freeze check with the user, NOT
a full open mini-gate) → build (dual-engine, both exit 0,
per-tooth ablation, independent `.wl`) → arbiter re-run → fresh-agent tri-review (fidelity + adversarial per-tooth +
transliteration screen) → **GLM-5.2 tertiary** → deliverables (note + `\StatusOpen` card into the Part-VII appendix + register
rows after R97) → deliverables fidelity-verify firewall → PDF rebuild → commit + tracker/memory sync.

**Carry the 044 process wins (all validated to work):**
- ⛔ **NO-CAN'T-FAIL-CONJUNCTS ban** — every conjunct mutation-corruptible; audit for `stored==stored` / `const==const` /
  `len==len` / `startswith`-literal / assigned-not-derived honesty tokens. (044's transliteration+adversarial legs caught 3 such
  sub-conjuncts; harden the BUILD prompt against them up front.)
- ⭐ **The guard triad** `FIRED_AT_OWN_ASSERT` + `MUTATION_DID_NOT_FIRE` + `UNKNOWN_MUTATION` (from the 044 build) is a REUSABLE
  ablation-integrity pattern — require it in the 045 build prompt (it made the ablation adversarially-genuine, not vacuous).
- ⚠ **ARCHIVE all prior `OUT_*` logs before each fresh Codex confirm** + an ISOLATED prompt forbidding `OUT_*`/archive reads
  (prevents the stale-log contamination that hit 043).
- **setsid-detached launches + the Monitor tool** (bounded poll for the done-marker; `run_in_background` waiters get reaped).
- **Thin-orchestrator delegation**: directive-draft-review-distill / bookend / build-monitor / tri-review / deliverables /
  tracker-sync all to sub-agents returning compact verdicts; distill every heavy `OUT_*`/build log via a fresh agent (keep
  transcripts out of context).
- **Deliverables fidelity-verify firewall** (the last honesty gate — caught nothing to fix in 044 but is load-bearing).
- **Mathematica ≤2 seats** (run `.wl` sequentially; the ablation runner invokes `.wl` per mutation — mind seat contention with
  an arbiter `.wl` run); **Grok may be usage-limited** → substitute a fresh independent compute agent.
- **GLM-5.2 tertiary is HELD-for-user-confirm** (standing confirm-first) — surface it at the gate: run on 045 (and optionally
  the still-owed 044 GLM) only on the user's OK.
- `/var/projects/toy_physics` NOT `toy_projects` (the attractor bit ~5× in the 044 run — stay vigilant; a `cd` traversing the
  wrong name FAILS, caught by `|| cd <correct>` fallbacks).

## Launchers + markers (unchanged; under `research/pde_ledger_v2/_scratch/`)
`run_codex.sh <PROMPT> <LOG>` (RO xhigh, `___CODEX_DONE___`) · `run_codex_build.sh` (danger-full-access xhigh,
`___CODEX_BUILD_DONE___`) · `run_grok.sh <ABS_PROMPT> <LOG>` (grok-4.5 high, `___GROK_DONE___`). Always `codex exec` at
`-c model_reasoning_effort=xhigh`.
