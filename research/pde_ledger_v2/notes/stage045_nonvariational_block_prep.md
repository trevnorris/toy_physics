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

## ⭐ THE DRAIN-PLACEMENT MINI-GATE (the USER's decision — deferred from 044; resolve BEFORE authoring the 045 directive)
Stage 044 named TWO candidate **non-variational** drain interfaces, proved `S_assembled` source-free, and left
`drain_equivalence=UNRESOLVED`, `drain_selection=DEFERRED_045_AND_USER_GATE`. **045 must surface the drain-placement to the user
as a short mini-gate before building** — it is a genuine physics/ontology call (like the Part VI cone-lock). The options:

- **(a) Part-I `Γ_B` order-conversion drain** (dynamical-wall regime): `∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B`,
  `Γ_B = Γ_return − Γ_drain`, dim `T⁻¹`. An INTERNAL order-conversion at the throat (ordered↔disordered material). **Total `n`
  (mass) is conserved exactly** — nothing leaves the medium; the throat CONVERTS phase (pin P12: phase-conversion, NOT suction).
- **(b) G0-card `S_drain` ρ-mass-sink + remote return** (frozen-wall regime — the G0 DRAFT-v0 choice; card §6): local sink
  `S_drain = −Σ_i Γ_0 D_i` (bump at `x=X_i`, `w=s_i a/2` body-side, width `a/4`) + a remote even-parity return
  `S_leakage = (Σ_i Γ_0) R_0` (at `w=±2L`) on TOTAL-`ρ` continuity `∂_tρ + ∇·(ρv) = S_drain + S_leakage`; charge-EVEN;
  dim `L⁻⁴T⁻¹`, `Γ_0*=10⁻³`. Mass is LOCALLY removed at the throat, GLOBALLY returned at the IR boundary (exact global closure).
- **(c) both as regimes** — (b) = the frozen-Σ static reduction of (a). ⚠ 044/Codex established these are DIFFERENT balance-law
  objects (internal-conversion vs sink+return); the reduction mapping is **NOT currently supplied**. Choosing (c) makes 045
  DERIVE that projection/reduction (real work; possibly a clean reduction, possibly a NO-GO if they don't reconcile).

**⭐ My recommendation to bring to the gate = (b), with the honest caveats**, because: (i) 045 builds the non-variational
complement of the **G0-card** conservative action `S_cons^G0`, and the G0 card §6 IS the (b) machinery already specified
(`S_drain`/`S_leakage`/`f`/`w`/`I_ret`/`P_ret` + BCs + force partition) — the concretely-buildable path; (ii) (b) matches the
user's "drain = the arrow of time" insight (`τ_d` T-odd, information lost at the throat — an irreversible sink, not a reversible
internal conversion); (iii) it keeps 045 consistent with what 044 assembled. **(a)** would re-derive the block from Part I's
`Γ_B` (less specified, and re-opens the wall dynamics 044 froze); **(c)** is the most ambitious (derive the (a)→(b) reduction —
could be its own able-to-fail result). ⚠ **Do NOT pre-decide — present all three + the recommendation and let the user choose
at the mini-gate.** The choice determines what the committed non-variational block is built on.

## What 045 BUILDS (the non-variational block — on the drain the user selects)
- **`S_drain` / `S_leakage`** — the mass-sink + global-return sources on the `ρ`-continuity (from the selected interface).
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
1. `software/em_charge_attribute/g0_closure_card_v0.md` **§6** (the `S_drain`/`S_leakage`/`f`/`w` + controllers `I_ret`/`P_ret`
   + the BC structure + the force partition) — the PRIMARY source for option (b) + the machinery. Note it is DRAFT-v0.
2. `docs/two_throat_simulation_handoff_spec.md` — READ-FIRST for the **`𝔅` BCs** (the 16-branch `𝔅` boundary-condition
   taxonomy + the R1→R4 ladder; partially-closed, NO closed action). The mouth-Robin/collar-Kirchhoff/IR-return BCs live here.
3. The **044 stage note** `notes/stages/ledger_stage044_parent_action_reconciliation.md` (what interfaces 044 NAMED + handed to
   045; the source-free `S_assembled`; the `F_var` piece from `−δS_hold/δX_i`).
4. Part I's `Γ_B` order-balance (stages 006/007) — the option-(a) drain (if the user leans (a) or (c)).
5. `research/pde_ledger_v2/notes/parameter_register.md` — the drain/return knobs already registered (`Γ_B` line 159; any
   `Γ_0`/`R_0`/controller rows); 045's new register rows start at **R98** (044 ended at R97).

## Process (full new-derivation gauntlet + GLM tertiary — 045 is foundational)
Author directive → **Codex→Grok→Codex bookend** → **drain-placement USER MINI-GATE** → build (dual-engine, both exit 0,
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
