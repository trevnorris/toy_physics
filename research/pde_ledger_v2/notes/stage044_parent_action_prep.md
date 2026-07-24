# Stage 044 (VII-2a) — the conservative parent action + `r_B`↔`χ_B` reconciliation: PREP (DRAFT, not the directive)

> **⚠ STATUS: prep/anchor, NOT a ratified directive.** Authored 2026-07-23 right after stage 043 (VII-1) COMPLETE, from the ratified
> Part VII split (`part7_integration_atomic_split.md`) + the earlier fresh-read of the unified-equation-set inputs (the G0 card, Part I's
> χ_B action stages 004–007, the Decision-16 P-retirement). The **directive** for 044 is authored AFTER a fresh read of the 044 sources
> (below) — this note is the resume anchor + the drain-placement gate framing. HEAD at authoring ≈ `f9ef70f6` (043 sync).

## Where 044 sits
Part VII is IN PROGRESS: **043 (VII-1, count-as-range) ✅ DONE** (`fda05f46`: continuous `[40,49]` + separate discrete `11`). The spine
was split (user Decision B) into **044 (VII-2a) = the conservative parent action + the `r_B`↔`χ_B` reconciliation** and **045 (VII-2b) =
the non-variational drain/return + BCs + force partition**. 044 is the HEAVIER, more foundational half; it is the genuinely-new synthesis
that assembles the medium/gravity/light/charge/magnetism sectors into ONE conservative action over ONE field set. ⚠ The G0 card is
**DRAFT-v0 / postulated, not a ratified stage** → 044 assembles at the **structure / completeness-floor level only** (it does NOT close a
committed parent action; it presents the assembled structure + reconciliation + the honest open pieces).

## ⭐ THE CRUX = drain-placement (a genuine modeling call → a short USER MINI-GATE before building)
`r_B` (G0 card) and `χ_B` (Part I) are the **SAME [0,1] wall order parameter** — identical double-well `aχ²(1−χ)²`, identical logistic
kink, identical `√(κ/2a)` width. So a naive sum double-counts the wall (that is the de-duplication). **But the substantive crux is NOT
naming — it is WHERE gravity's drain sits:**
- **Part I:** the drain/return is on the **ORDER FIELD** — `Γ_B` (a `χ_B: 1→0` order-conversion at the throat, a charge-source term).
- **G0 card:** FREEZES the wall (`Γ_B` dynamics OFF) and puts the drain on the **BULK CARRIER `ρ`** as a charge-even **mass-sink**
  `S_drain`.
These are two DIFFERENT mathematical objects for "gravity's drain." **044 must adjudicate which** (or how they reconcile) — this
determines the parent action's structure and what gets pushed to 045's non-variational block. ⚠ **Surface this to the user as a short
mini-gate before building 044** (it is a real physics decision, like the Part VI cone-lock adjudication). Frame the options: (a) drain =
`χ_B` order-conversion `Γ_B` (Part-I native); (b) drain = `ρ` mass-sink (G0 native, wall frozen as a static reduction); (c) they are the
same object in two regimes (frozen-Σ = the static MVP reduction of the dynamical `χ_B`) — my current lean, to be validated at the read.

## The parent action + the field set (what 044 assembles)
- **Conservative parent action** `S_cons = S_bulk + S_χ + S_scalar + S_mouth + S_hold + S_geon`:
  - `S_bulk` — medium + gravity-flow (Madelung `ℏρ∂_tθ` + quantum pressure + `U=(K/4)ρ⁵`; `θ` = flow phase, `v=(ℏ/m)∇θ`) [I,II].
  - `S_χ (+S_hold)` — the ONE wall (the de-duplicated `χ_B≡r_B`; `S_hold` = the holonomic freeze that reduces it to the static MVP) [I].
  - light — `L_Mac + L_uw` (transverse `u_T`, `u_w`), `χ_B`-gated, `c_γ²=μ_R/ρ_br` [III].
  - charge — `S_H` / reduced `S_Lh(u_L,h) + S_mix(C_hu∇u_L·∇h) + S_mouth(J_m Q_χ H)` [IV].
  - magnetism = the boost of charge — **no new field** [V].
- **Field set (confirm at the read):** `{ρ, θ, r_B≡χ_B, H/h, u_L, u_T, u_w}`. ⚠ ADD `u_w` (the gapped transverse-w brane DOF, stage007
  `L_uw`, `Ω_w` gap — a genuine operative field). `θ_B` (brane-order phase) is **eliminated** (`A_eff=ρ_br+C_J²/κ_phase` absorbs it) —
  correctly not independent.
- **Notation map (owed):** G0 `ρ` ≡ Part-I `n` (both `L⁻⁴`); G0 `r_B` ≡ Part-I `χ_B`; `θ_B` eliminated. **Three** historical wall
  descriptions collapse to one: fixed `g_ℓ(w)` profile → `χ_B` field → `r_B` frozen (stage007 superseded `g_ℓ` with `χ_B`, register R21).
- **Decision-16 P-retirement bookkeeping:** the retirement SIMPLIFIES the reconciliation (the wall is now unambiguously the SCALAR `χ_B`,
  route (a); no competing polar-vector wall). Debris to carry: (1) the stage007 freeze hash byte-INCLUDES the retired `P` terms — the
  operative action EXCLUDES them (a symbolic summand-set partition over an immutable hash, NOT byte surgery); (2) the route-(c)
  `χ_B=|P_∥|²` gate is retired-but-not-foreclosed (a named Part-VII-adjacent gate needing a new T0 freeze if ever reopened).

## OUT of 044 (→ 045, VII-2b): the NON-VARIATIONAL block
The drain/return machinery is 045, NOT 044: `S_drain`/`S_leakage` mass-sink + global-return on the `ρ`-continuity, the finite-rank
`I_ret`/`P_ret` return controllers (exact global mass/momentum/energy closure), the momentum `f_drain`/`f_leakage`, the energy
`w_drain`/`w_leakage`, the `𝔅` sleeve-contact + mouth-Robin/Neumann + collar-Kirchhoff + two-sided-IR-return BCs + the phase/moduli
zero-mode quotients, and the `F_total = F_var + F_flux + F_𝔅 + F_rad` force partition. **044 NAMES the interface (which drain object it
hands to 045) but does not build the non-variational machinery.**

## Dual-engine vs prose (044's verification contract)
- **Dual-engine-checkable (SymPy + independent Mathematica):** the coupled EL field equations; the dispersion / wave-speed eigenvalues
  (`c_±²=(3±√2)/2`, `c_γ²=μ_R/ρ_br`); the wall/`H` factorizations (`L_χ=A_χ†A_χ`, `O_⊥=A†A`, `O_⊥f₀=0` gapless zero mode); the stability
  Hessian; the dimensional homogeneity of every `S_cons` term (units-restored, able-to-fail).
- **Prose (PRINTED, not asserted):** the variational/non-variational partition rationale; the reconciliation + P-retirement narrative;
  the completeness-floor / sector-assembly argument; the G0-card-is-DRAFT-v0 caveat.
- The assembled two-body / far-field solve is **SIM-DEFERRED** (the shared throat solve = the central reduction debt) — 044 NAMES it,
  does NOT perform it.

## ⭐ Fold the gravitomagnetism pointer (per the earmark, committed `dc617088`)
In 044's parent-action gravity sector, **flag the gravitomagnetic (velocity-dependent flow) terms** — frame-dragging = the boost of the
gravitoelectric drain/return flow, the gravity twin of magnetism-as-boost. Reproduced by the GR-matched PN ladder; record the honest
asymmetry (gravity matches GR cleanly incl. gravitomagnetism, EM departs from Maxwell). (The register/model_map one-liner is 047's job.)

## Fresh-read sources for the 044 gate (read before authoring the directive)
1. `software/em_charge_attribute/g0_closure_card_v0.md` (the `r_B` wall action + `S_drain` mass-sink + `S_hold` freeze; note it is DRAFT v0).
2. Part I's `χ_B` action: stages 004–007 (`notes/stages/ledger_stage004..007_*.md` + `paper/stages/stage_004..007.tex`) — the `χ_B`
   two-phase action, `Γ_B` order-conversion drain, `χ_B f_shear` light gating.
3. `software/stage1_solver/decisions/16_retire_brane_polar_field.md` (the P-retirement + the `INSTABILITY_CONFIRMED_STRUCTURAL` kill + the
   operative `{S_GNLS, gL_Mac, gL_uw}`+χ_B action).
4. The assembled sector actions already in-ledger: Part II (gravity/PN), Part III (stage003 light `L_Mac`/`L_uw`), Part IV (030–033
   charge `S_H`/`S_mix`/`S_mouth`), Part V (034–039 magnetism = boost).
5. `research/pde_ledger_v2/notes/parameter_register.md` (the field/knob provenance — `u_w`, `θ_B`-eliminated, `A_eff`, `c_γ²=μ_R/ρ_br`).

## Process (full new-derivation gauntlet + GLM tertiary — 044 is a foundational stage)
Author directive → **Codex→Grok→Codex bookend** → **drain-placement USER MINI-GATE** → build (dual-engine, both exit 0, per-tooth
ablation, independent `.wl`) → arbiter re-run → fresh-agent tri-review (fidelity + adversarial ablation + transliteration screen) →
**GLM-5.2 tertiary** (foundational) → deliverables (note + `\StatusOpen` card into the Part-VII appendix + register rows after R92) →
deliverables fidelity-verify firewall → commit + tracker/memory sync. Carry the process wins: ⛔ NO-CAN'T-FAIL-CONJUNCTS build-prompt ban
(harden vs `startswith`-literal / stored-self-compare — the 2 nits 043 shipped with); the emphatic one-shot/busy-poll rule on every
bookend/verify gauntlet; **archive prior `OUT_*` logs before each fresh Codex confirm** (043 hit the stale-log contamination — Codex
echoed old verdicts; the fix = archive + isolated prompt forbidding `OUT_*`); the deliverables fidelity-verify firewall (it caught 043's
"no can't-fail conjuncts" over-claim); Grok may be usage-limited → substitute a fresh compute agent; `/var/projects/toy_physics` NOT
`toy_projects` (the attractor bit 4× in 043's agent prompts — stay vigilant).
