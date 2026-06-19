# Directive pathA_11 — Chunk B2b: DERIVE the Maxwell transfer on the Path-A self-consistent background

**Owner:** Codex (codes + DESIGNS the route: the Spike-1→Spike-2 adaptation to the Path-A `A0` background, and an
independent Python cross-check engine). Claude reviews afterward (standing transliteration-fidelity audit per math
module + adversarial). **Status:** post-freeze (GATE-A `frozen: YES`, commit `1703f4c`, `candidate_freeze_hash
ed3585…`); B1 extraction DONE (`398ba27`); **B2a DONE + committed (`fa5712e`)** — the Path-A closed-background
exporter and the derived BdG bundle `{ϖ_j, c_j, B0, B2, B4}` already exist; architecture = decision-12 Reading A
(derive, don't inherit).

## Scope of THIS chunk (B2b only)
DERIVE the **Maxwell transfer** `{Z0, Z2, Z4, N0, N2, N4}` on the **Path-A self-consistent closed background** —
specifically on the exported EM/gauge background `A0 = {A_00, A_r0, A_w0}` from the B2a closed-background bundle (the
frozen `homogeneous_isotropic_hooke_v1` family, at a neutral reference `τ`) — via the **basis-invariant Spike-1→
Spike-2 Green/self-energy transfer** (NO posited U/W ports). Output: the Maxwell bundle that B1's extraction module
consumes through its `direct_coefficients` lane (`patha_extraction.lane_extract`).

**OUT of scope (do NOT do here):**
- NO `R_norm(τ)=0` root-find, NO `R_norm`/`R_pole`/`P2`/`P4`/`D0`/`P0` assembly — that is **B2c**. In particular
  **do NOT run the `mt15_05` Spike-3 R_norm "preview"** (it computes `D0 = K−B0−Z0` and `R_norm` internally; that
  step is B2c). B2b stops at `{Z_n, N_n}` and their validation only. (Same discipline B1/B2a followed: prove the
  piece, don't run the physical verdict.)
- NO re-derivation of the matter-sector BdG bundle — that is B2a (done). **Reuse** the B2a bundle / wall mode `χ`.
- Stays target-blind: emit the Maxwell bundle + its validation; compute nothing that reveals a held-out residual.

## Why this chunk exists (the gap decision-12 surfaced)
The earlier M1c effective-closure run did **not** derive the Maxwell transfer on a self-consistent background: the
mixed-Maxwell sector was **posited** as U/W mixed ports (`Ω_U=3.25, Ω_W=4.35, g_U=0.18, g_W=0.13, R=0.08`), and the
Spike transfer ran on a **smoke** geometry (`mt15_03`/`mt15_04` hardcode `R0(w)=a+(R_exit−a)(3x²−2x³)`,
`Z(w)=Z_floor+Z_amp·Exp[…]`, with M1a open-throat constants). So the M1c `Z0=2.50e-6, N0=2.67e-6` and the −10.8
miss rested partly on posited ports + a smoke background. Reading A requires `{Z_n, N_n}` **DERIVED** on the genuine
Path-A `A0` background via the basis-invariant Green/self-energy transfer — the U/W ports removed entirely from the
live path. The recurring project sin is a faithful-but-wrong operator (MMS/arbiter cannot catch it) and a
tautological / can't-fail gate — every validation gate below must be able to FAIL, and the operator + transfer must
be fidelity-audited term-by-term against the canonical source.

## Canonical source of truth (transliteration target — operator FORMS must match these)
There is **no standalone Maxwell notes file**; the derivation lives in the parent compact + the reduced audit notes
+ decision-05. The new Path-A scripts must keep the operator/transfer FORMS identical to these; only the background
fields + geometry change.
- **Parent action / PDE:** `notes/moving_throat_pde_program_compact.md` — Maxwell Lagrangian L592–620, exact
  localized Maxwell PDE L674–689, mixed gauge invariants L769–786, `Z_n` self-energy structure L1500–1509, `D_n`
  moments L1517–1528, `Γ_port = 4a⁵/(27 c_s⁵)` L2559–2568.
- **Reduced Maxwell/mixed kernel:** `research/pde_audit/notes/stage_v2_09_maxwell_mixed_kernel_derivation.md` —
  conservative self-energy + `Z₀=Q/∆` L117–194; outgoing `N₀=P²/∆²` L299–376; `Γ_port` reduced form L384–420.
  (READ-ONLY reference.)
- **Grouped Z_n/N_n extraction fixture:** `research/pde_audit/notes/stage_v2_21_branch_extraction_fixture_derivation.md`
  L166–260 (`Z_n` L203–211, `N_n` L227–241, `D_n` L250–259). (READ-ONLY reference.)
- **Basis-invariance design (the decisive point):** `software/stage1_solver/decisions/05_mixed_maxwell_spike_design.md`
  (D4 crux). The live transfer is **basis-invariant**: `Z_n = Re Σ_cons(0)` (`Z = gᵀ·inv(K)·g`),
  `N_n = −Im Σ_ret/(Γ_port ω⁵)` from the outgoing self-energy. The `v2_09`/`v2_21` U/W-port formulas
  (`Q/∆`, `P=Ω_U²g_W+R g_U`, `P²/∆²`) are the canonical-PHYSICS cross-reference and the regression check ONLY — they
  must NOT appear in the live extraction path (decision-12 Maxwell sub-decision).
- **Existing Spike Term Maps (term-by-term provenance):** `mt15_03_spike1_vsh_maxwell_operator_report.md` L91–103;
  `mt15_04_spike2_transfer_n0_report.md` L88–97.

## Build (requirements — Codex designs the implementation)
1. **Reuse the B2a Path-A closed background (do NOT re-solve a smoke background).** Read the exported EM/gauge
   background `A0 = {A_00, A_r0, A_w0}` (on `r_centers × w_centers`) plus the grid/geometry/constants from the B2a
   closed-background bundle (the exporter in `patha_b2a_bdg.py`: `solve_closed_background` / `export_background`),
   at the **neutral reference `τ=1.0`** (a convention, NOT a fit — must NOT be chosen by peeking at any R-target).
   - Confirm the background is the converged Path-A **closed** solve (residual ≪ the smoke `243.39`; report it) and
     uses the **frozen geometry** `a=1.0, L=37/20=1.85`. A run that silently fell back to the `mt15_03`/`mt15_04`
     smoke `R0(w)/Z(w)` or M1a geometry must FAIL the self-consistency gate.
   - If the existing B2a bundle does not carry enough of `A0` (or the right faces/measure) for the VSH transfer,
     EXTEND the B2a exporter to emit what's needed — but it must remain the SAME closed solve under the SAME frozen
     family (no new background physics), and the extension must not change any B2a-exported BdG quantity.
2. **Adapt the Spike-1→Spike-2 chain to the Path-A `A0` background (Mathematica).** Create Path-A variants of
   `mt15_03_spike1_vsh_maxwell_operator.wls` (VSH Maxwell operator, 5 lanes `phi,a_r,a_E,a_B,a_w` at `ℓ=2`,
   `L=l(l+1)=6`) and `mt15_04_spike2_transfer_n0.wls` (Green/self-energy transfer → `{Z_n,N_n}`) — name them
   `mt15_03_patha_b2b_*.wls` / `mt15_04_patha_b2b_*.wls`, dispatched via a guard mirroring the B2a `--patha-b2a`
   idiom (a `If[MemberQ[$ScriptCommandLine, "--patha-b2b"], Get[…]; Exit[0]]` guard in the originals; originals'
   legacy smoke behavior untouched). The Path-A variants must:
   - READ the Path-A `A0` background bundle instead of the hardcoded smoke `R0(w)/Z(w)` profile, and use the frozen
     geometry. Resample/interpolate `A0` onto the VSH transfer mesh as needed.
   - Keep the **operator + transfer FORMS identical** to the canonical source above: the VSH localized Maxwell
     operator terms, the gauge deflation (W-orthogonal), the **basis-invariant** self-energy `Z_n = Re Σ_cons(0)`
     (via `Z = gᵀ·inv(K)·g`), the outgoing DtN condition `K_ret = K_cons + i Γ_port ω⁵ B_tan` with
     `Γ_port = 4a⁵/(27 c_s⁵)`, and the outgoing forward source `N_n = −Im Σ_ret/(Γ_port ω⁵)`. Only the background
     `A0` + geometry change, NOT the operator/transfer terms. **NO posited U/W ports in the live path** (keep the
     `v209Regression` U/W toy strictly as a basis-invariance regression check — gate 3 — not as the source of
     `Z_n/N_n`).
3. **Wall/BdG coupling reuse (the forward source `N_n`).** The forward source couples the gauge to the wall mode and
   the BdG response. Reuse the **B1-derived** wall mode `χ(w)` (from
   `patha_extraction.solve_l2_wall_eigenproblem` on this SAME background; `χᵀWχ=1`, `∫χ>0`) and the **B2a** BdG
   bundle on the SAME background — NOT a re-posited wall shape or re-derived spectrum. Form the gauge↔wall overlaps
   via the SAME algebra convention as `research/pde_audit/scripts/stage_v2_22a_profile_to_coefficient_adapter.py`
   (READ-ONLY reference). Convert to the low-frequency moments `{Z0,Z2,Z4}` and `{N0,N2,N4}` exactly as
   `patha_extraction`'s `direct_coefficients` lane consumes them.

   **Correction note (2026-06-19):** The implemented and fidelity-audited forward-source construction is the
   canonical decision-05 D3 Fréchet source over the B2a BdG-mode response. The earlier `χ`/`v2_22a` wording here was
   imprecise: the B2a bundle already carries the χ coupling through its overlaps `c_j`; B2b must not introduce a
   separate live B1-χ overlap path.
4. **Independent Python cross-check (dual-engine — REQUIRED).** Independently assemble the SAME VSH Maxwell operator
   and compute `{Z_n,N_n}` via the SAME basis-invariant Green/self-energy transfer on the SAME exported `A0`
   background in Python, with an independent linear-algebra stack (scipy/torch). **Note:** the existing Python
   Maxwell operator `operators.py:396` (`localized_maxwell_operator`) is **stationary axisymmetric `l=0` only** (3
   lanes `A_0,A_r,A_w`) and is **NOT a drop-in** for the `ℓ=2` 5-lane VSH open-boundary Green transfer — Codex must
   build a genuine `ℓ=2` transfer engine (mirroring how B2a built a new complex-BdG Python engine rather than reusing
   the stationary tangent). The two engines must NOT call each other — genuinely independent assembly. Compute
   `{Z_n,N_n}` in both and compare (gate 2).
   - **If, after honest effort, a fully-independent Python transfer is genuinely infeasible (not merely laborious),
     HALT and surface the specific obstruction to the user.** Do NOT silently downgrade to single-engine or to a
     trivial partial check ([[feedback-dual-engine-required]], [[feedback-never-alter-calibrated-process]]).

## Independent validation (the heart of B2b — every gate must be able to FAIL)
1. **Self-consistency (catches: running on the smoke/wrong background).** Assert `A0` used is the converged Path-A
   closed solve (residual ≪ smoke `243.39`; report it), geometry `a=1, L=1.85`. A silent fallback to the
   `mt15_03/04` smoke `R0(w)/Z(w)` or M1a geometry must FAIL here.
2. **Dual-engine agreement (catches: an engine-specific assembly/solve bug).** MMA vs the independent Python engine
   agree on `{Z0,Z2,Z4,N0,N2,N4}` to a stated tolerance — **abs AND rel** (the `Z_n,N_n` are ~`1e-6`; an abs-OR-rel
   gate would pass a wrong answer for free — this was the B2a M1 fix). Report max abs + rel diff per coefficient.
   Genuinely independent assembly → disagreement flags a real bug, not a copy compared to itself.
3. **Basis-invariance (catches: a basis-dependent / posited-port leak in the transfer).** The transfer is
   basis-invariant by construction (`Z = gᵀ·inv(K)·g`); confirm `{Z_n,N_n}` are invariant under a change of lane
   basis / gauge deflation (the `v209`-style regression, exercised on the Path-A transfer, or an equivalent
   basis-rotation check). State the wrong answer: a transfer that secretly depends on a posited U/W eigenbasis (the
   exact failure decision-05 D4 was designed to prevent).
4. **Outgoing physicality (catches: a wrong-sign / non-radiating / bound spurious transfer).** `N0 > 0` (mandatory
   outgoing flux — `patha_extraction` requires it), `Σ_cons(0)` real, `−Im Σ_ret > 0` (positive radiated power),
   all `{Z_n,N_n}` finite. State: an anti-radiating or bound spurious solution would be caught here.
5. **Convergence + truncation (catches: an unconverged / under-resolved transfer).** Sweep (a) the transfer mesh ≥2
   levels, (b) the **ω-window** used for the low-frequency expansion that yields `Z2/Z4/N2/N4`, and (c) the
   open-boundary / DtN radial truncation (domain length / absorber). Show `{Z_n,N_n}` settle (values + trend), pick
   a resolution/window by a stated tolerance, and **carry the residual error per coefficient into the bundle
   `error_budget`** (B2c §J consumes it). A resolution-, window-, or truncation-dependent `{Z_n,N_n}` must be caught
   here, not silently shipped. (This is the direct analogue of the B2a modal-truncation lesson.)
6. **τ-sensitivity (catches: a background NOT re-solved per τ, i.e. a frozen/stale bundle).** Derive `{Z_n,N_n}` at
   `τ=1` and one other (e.g. `τ=2`); show they genuinely MOVE with τ through the re-solved background (decision-12
   "confirm the τ-scaling early"). Equal bundles across τ would mean the background isn't actually re-solved. This is
   machinery confirmation for B2c's sweep, NOT a calibration — no `R_norm` computed.
7. **B1 consumer compatibility (catches: a bundle the extraction can't ingest / a cross-engine mismatch through the
   consumer).** Confirm the derived `{Z_n,N_n}` (combined with the B2a `{B0,B2,B4}` and B1 `{K,M}`) form a
   well-shaped `direct_coefficients` lane: all 11 required keys present, finite, `N0>0`, loadable by
   `patha_extraction.lane_extract`'s coefficient-loading path — and that the **MMA-derived** `{Z_n,N_n}` match the
   **Python-derived** `{Z_n,N_n}` through that loading path (cross-engine, not a self-comparison — the B2a R2
   lesson). **Do NOT invoke the `D0`/`P0`/`R_norm` assembly** (that is B2c); validate up to the coefficient dict only.
8. **No can't-fail gates.** For EVERY gate above, state in the report the specific wrong answer it would catch. No
   verbatim-copy "independent" recompute, no number compared to itself, no tolerance loose enough to pass a wrong
   operator/transfer.

## Acceptance criteria (Codex iterates until ALL pass, exit 0)
1. Reads the B2a Path-A **closed**-background bundle's `A0` at `τ=1` (not smoke), frozen geometry, confirms a small
   converged residual (reported). (Exporter extended only if needed, same closed solve, B2a BdG quantities unchanged.)
2. The adapted Spike-1→Spike-2 Path-A variants read that `A0`, use frozen geometry, and compute `{Z0,Z2,Z4,N0,N2,N4}`
   via the **basis-invariant** Green/self-energy transfer — NO posited U/W ports in the live path.
3. Dual-engine (MMA vs independent Python) agreement on `{Z_n,N_n}` to the stated tolerance — **abs AND rel**
   (report diffs per coefficient). (Or a documented HALT-to-user if a genuine independent Python engine is infeasible.)
4. Basis-invariance demonstrated (gate 3); outgoing physicality demonstrated (`N0>0`, positive radiated power; gate 4).
5. Convergence + ω-window + truncation demonstrated (≥2 levels each; values + trend reported); residual error per
   coefficient recorded in the bundle `error_budget`.
6. τ-sensitivity demonstrated (`τ=1` vs `τ=2`; coefficients move).
7. The bundle is emitted in the B1-ingestible `direct_coefficients` shape; the consumer-compatibility check (gate 7)
   passes WITHOUT computing `R_norm`/`D0`. Full `pytest` for the `patha_*` suite stays green (no regression).
8. NO `R_norm`/`R_pole`/`P2`/`P4`/`D0`/`P0`/root-find computed anywhere; the `mt15_05` R_norm preview is NOT run.
9. Constraints honored: every script wrapped `timeout 600` (exit 124 = reformulate, never raise the cap); **≤2
   concurrent `math -script` seats**; CPU sparse-direct (GPU off); reasoning-only steps consume 0 seats. Firewall
   untouched — do NOT write under or import from `research/pde_audit/simulation/`, do NOT touch
   `physical_export_permitted`; the `research/pde_audit/{scripts,notes}/` references are READ-ONLY. No `git add` / no
   commit (orchestrator commits after review). YAML/markdown for human-read notes; JSON only for the machine bundle.

## Report back
The derived `{Z0,Z2,Z4,N0,N2,N4}` at `τ=1`; the confirmed background residual norm (vs the smoke `243.39`); the
dual-engine max abs+rel diffs per coefficient; the convergence + ω-window + truncation table + chosen resolution +
the recorded per-coefficient error bars; the τ-sensitivity result (`τ=1` vs `τ=2`); the basis-invariance result; the
B1 consumer-compatibility result; for EACH validation gate the specific wrong answer it would catch; the files
created/modified; and an honest note on the dual-engine Python feasibility (was a genuine independent `ℓ=2` transfer
engine built, or did anything stay shared single-engine — name it). Any place the Spike operator/transfer had to be
adapted (and why) rather than reused verbatim.
