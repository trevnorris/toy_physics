# Decision 12 — Path-A chunk B2 architecture: DERIVE the BdG/Maxwell back-reaction bundle on the frozen background

**Date:** 2026-06-18
**Status:** DECIDED (user-authorized 2026-06-18). Reading **A** (fully derived) chosen. This is the architecture +
detailed chunk plan for B2 (the decisive `R_norm(τ)=0` calibration). It is the **resume-here record after /compact.**
**Mechanism:** Claude+Codex architecture consult (`_scratch/pathA_b2_architecture_consult*`) → user methodology
call. Builds on: GATE-A freeze `decisions/11` (`frozen: YES`, commit `1703f4c`, hash `ed358569…b1691c9`); B1
extraction module `decisions/11` §5 → committed `398ba27`.

## The problem B2 surfaced
The B1 extraction module computes `R_norm` GIVEN the wall channel `{K,M,χ}` PLUS the BdG bundle (`c_j,ϖ_j,B_n`) and
the Maxwell bundle (`Z_n,N_n`). **Codebase map finding:** no BdG eigensolver and no Maxwell mixed-port/transfer
solver exists in the Python solver tree. In the earlier M1c effective-closure run those were SUPPLIED as M1b
packets — and the M1b BdG ran on a *smoke/high-residual* background ("machinery, not final physics") while the
mixed-Maxwell ports were **posited** (`Ω_U=3.25, Ω_W=4.35, g_U=0.18, g_W=0.13, R=0.08`). So the −10.8
effective-closure miss was partly resting on posited modes.

## The decision (Reading A — derive, don't inherit)
**Derive the BdG spectrum and the Maxwell transfer numerically on the Path-A self-consistent background, under the
frozen conventions.** Rejected: Reading B (inherit/freeze numeric `ϖ/Ω/λ` scales, recompute only profiles/overlaps)
— hash-compatible and smaller, but materially inherits/posits more of `R_norm`, against the "derive don't fit"
methodology ([[feedback-calibrate-predict-methodology]]). User directive (2026-06-18): "let's not shy away from the
hard work… the answers will be found in the weeds of hard work"; "nothing we've done is sacred" (re revising the
freeze).

### Freeze compatibility (important)
The GATE-A freeze did **NOT** pin numeric `ϖ/Ω/λ` values or mode profiles — the canonical `freeze_sheet.json`
serializes only forms/conventions (family, ties, geometry, objective, channel, extraction formulas). **So Reading A
is compatible with the existing freeze hash `ed3585…` — NO re-freeze of numbers is required.** What is frozen and
must NOT change: the extraction formulas, family ties, geometry, calibration objective, channel selection,
tolerances. The §5b WORDING in `decisions/11` (which ambiguously listed `λ/Ω/ϖ` as "frozen conventions") gets a
dated CLARIFICATION note (markdown only; not in the hashed JSON → hash unaffected): `ϖ_j,Ω` are DERIVED outputs of
the BdG/Maxwell eigen/transfer solve on the Path-A background; only the coupling NORMALIZATION conventions (`λ`
definitions) and the extraction formulas are frozen; the Maxwell sector uses the Spike-2 basis-invariant Green
transfer (direct `{Z_n,N_n}`), NOT the posited U/W ports.

### Maxwell sub-decision (user-approved)
Use the **Spike-2 basis-invariant Maxwell Green/self-energy transfer** to produce direct `{Z_n,N_n}` — NOT the
posited U/W mixed-port form (`Ω_U,Ω_W,g_U,g_W`). Removes those inherited posited scales. The B1 extraction already
supports a `direct_coefficients` path (`patha_extraction.py` `lane_extract`) and a `derived_maxwell_transfer` path
(`stage_v2_22b`), so direct `{Z_n,N_n}` feed straight in; `D0=K−B0−Z0` is unchanged regardless of how `Z0` was
obtained.

## Reusable derivation engines (Mathematica chain — re-point at the Path-A background)
The reusable derivation path is the Mathematica M1b/Spike chain, NOT a Python rebuild (the Python JVP in
`patha_closed_newton.py:312` / `newton.py:68` is the *stationary* Newton residual — not a drop-in dynamic/symplectic
BdG eigensolver; the Python Maxwell operator `operators.py:396` is stationary axisymmetric, not the 5-lane
VSH/open-boundary Green transfer):
- **BdG:** `software/stage1_solver/mathematica/mt15_02_bdg_wall_derivation.wls` (assembles finite BdG matrix, solves
  eigensystem, normalizes modes, exports `ϖ_j,φ_j,c_j`) + report `mt15_02_bdg_wall_derivation_report.md`.
- **Maxwell transfer:** `mt15_03_spike1_vsh_maxwell_operator.wls` (Spike-1: VSH Maxwell operator) →
  `mt15_04_spike2_transfer_n0.wls` (Spike-2: basis-invariant Green/self-energy transfer → `{Z_n,N_n}`) →
  `mt15_05_spike3_clean_rnorm.wls` (Spike-3: clean R_norm assembly) + their reports.
- **Background source:** the Path-A frozen closed solve `patha_closed_newton.py` (state `{ψ,A,R0,μ}` on the frozen
  `homogeneous_isotropic_hooke_v1` family) + `coupled_branch.py` background export; B1 `patha_extraction.py` for
  `χ,K,M`.
- **Constraints:** Mathematica ≤2 concurrent `math -script` seats ([[feedback-mathematica-single-seat]]); dual-engine
  required wherever MMA can independently verify ([[feedback-dual-engine-required]]); `timeout 600` on every script;
  GPU off → CPU sparse-direct ([[project-gpu-disabled-machine]]).

## Baseline numbers being REPLACED (provenance — NOT Path-A frozen inputs)
From M1c effective-closure (`frozen/m1c/834835…/m1c_physical_derived_report.md`): `K=4.060384`, `B0=0.00465942`,
`Z0=2.5002e-6`, `N0=2.6745e-6`. M1b BdG `ϖ=6.4327,7.6951,10.7526`. Posited ports `Ω_U=3.25,Ω_W=4.35,R=0.08,
g_U=0.18,g_W=0.13`. **Honest prior** (`K=4.06 ≫ B0+Z0≈0.0047`, ~`1.6e7` cancellation, decision-08/decision-11 §6)
came from this effective-closure branch, NOT the Path-A harmonic solve. Under Reading A, `B0/Z0/N0` + the spectra
move with the actual Path-A background → can change both `R_norm` and the naturalness verdict. A miss is handled
κ_PV-style (decision-11 §4b), never a rescue DOF.

## The B2 chunk plan (sequential, explicit user gate per chunk — [[feedback-sequential-audit-chunks]])
**B2a — Derive the BdG spectrum.** Adapt `mt15_02` to solve the BdG eigenproblem on the Path-A self-consistent
`{ψ0,A0,R0}` (frozen harmonic family) → `ϖ_j`, mode profiles `φ_j`, normalized; overlaps with the B1 wall mode `χ`
→ `c_j`, then `B_n=Σ c_j²/ϖ_j^{2(n+1)}`. Dual-engine (MMA derive + Python cross-check) + transliteration-fidelity
(clean agent) + adversarial. Output: derived `{c_j,ϖ_j,φ_j,B0,B2,B4}` on the Path-A background.
**B2b — Derive the Maxwell transfer.** Re-run the Spike-1→Spike-2 chain (`mt15_03`→`mt15_04`, and `mt15_05` for the
clean assembly) on the Path-A `A0` background → direct `{Z0,Z2,Z4,N0,N2,N4}` via the basis-invariant Green/self-
energy transfer (no U/W ports). Same audit stack. Output: derived `{Z_n,N_n}`.
**B2c — Integrate + calibrate.** Feed derived `{K,M,B_n,Z_n,N_n}` (B1 `χ,K,M` + B2a + B2b) → B1 `patha_extraction`
→ the unique `R_norm(τ)=0` deterministic root-find on the stable-side `D0>0`. NOTE: each τ → re-solve the closed
background → re-derive the BdG+Maxwell bundle → extract `R_norm(τ)` (the bundle is genuinely τ-dependent through the
background; confirm the τ-scaling early — `K=τκ̂` is exact, but `{ψ0,A0,R0}` and hence `{B,Z,N}` shift with τ via
the wall stiffness). §J: closed grid-convergence with frozen forms, calibration-covariance into held-out
`R_pole/P2/P4`, margin-to-Schur (`D0`) error bars. Report `τ*`, naturalness (`|ln τ*|`, `K/(B0+Z0)` cancellation
ratio + digit count), held-out `R_pole/P2/P4`. **This LEAVES target-blind.**

## B2a STATUS — ✅ DONE + VALIDATED (2026-06-18)
Directives `pathA_09` (build) + `pathA_10` (remediation). Built: a NEW Path-A closed-background exporter
(`patha_b2a_bdg.py`), the adapted MMA engine (`mt15_02_patha_b2a_bdg_wall_derivation.wls`, dispatched via a 6-line
`--patha-b2a` guard in the original `mt15_02`), tests, and report (`reports/patha_b2a_bdg_derivation_report.md`).
**Derived τ=1 bundle (converged, K=30 of 100 positive modes):** `ϖ=[4.5287, 5.9597, 10.1091]`,
`B0=3.9010e-5, B2=1.7734e-6, B4=8.4073e-8`; closed background residual `2.749e-10` (vs smoke `243.39`).
**Review:** two independent transliteration-fidelity audits → operator FAITHFUL (the new engine generalizes the old
real-ψ shortcuts to correct complex-BdG forms: `(dh/dρ)ψ²` pairing, `Re(ψ̄u+ψv)` density response, conjugated lower
block). Adversarial review → no fatal flaw; caught a +3.0% B0 modal-truncation miss (3-mode sum) → remediation swept
to K=30 (modal error 1.4e-5) and exports the **converged** moments. Tautological consumer gate replaced by a genuine
cross-engine check; dual-engine gate tightened to `abs AND rel`. 23 tests pass. Target-blind preserved (no
`R_norm`/Maxwell/root-find).
**Error budget carried for B2c §J** (in the bundle `error_budget` block): modal-truncation `1.4e-5`,
spatial-confirmation (12×12) `9.3e-3`, spatial-ladder (10×10 vs 8×8) `2.0e-2` (max rel over B0/B2/B4).
**τ-finding for B2c:** doubling τ moves the *matter* background <1% (τ enters only the wall sector in the frozen
Hooke family) → R_norm(τ) leverage is dominated by `K=τκ̂` (exact) + the wall/Maxwell sectors, NOT by `B_n`.

## B2b STATUS — ✅ DONE + VALIDATED (2026-06-18)
Directives `pathA_11` (build) + `pathA_12` (remediation #1) + `pathA_13` (remediation #2). Built: a NEW Path-A VSH
Maxwell transfer engine (`src/stage1_solver/patha_b2b_maxwell.py` — primary 2nd-order FD **and** a genuinely-
independent 4th-order staggered FD engine), adapted MMA Spike variants (`mt15_03_patha_b2b` / `mt15_04_patha_b2b`,
dispatched via `--patha-b2b` guards in the originals), tests (32 pass), and report
(`reports/patha_b2b_maxwell_transfer_derivation_report.md`). Derived on the **same B2a closed background's `A0`** at
neutral `τ=1` via the **basis-invariant Green/self-energy transfer** (`Z_n=Re Σ_cons(0)=gᵀ·inv(K)·g`,
`N_n=−Im Σ_ret/(Γ_port ω⁵)`, `Γ_port=4a⁵/27c_s⁵`) — **NO posited U/W ports** in the live path.
**Converged τ=1 bundle** (grid 47×17, window 0.028, radial_scale 5.0):
`Z0=2.395e-8, Z2=-1.528e-8, Z4=5.017e-9, N0=2.158e-8, N2=-5.934e-9, N4=3.698e-9`.
**Per-coefficient error bars:** `Z0/Z2/N0/N2 ≈ 0.4–0.6%`; `Z4 ≈ 2.2%` (true ~4–5% incl. geometric tail), `N4 ≈ 1.7%`;
cross-engine (2nd- vs genuinely-different 4th-order) ≈ 4.4% ≈ the convergence error (honestly recorded, not
overclaimed). All in the bundle `error_budget` for B2c §J.
**Structural win confirmed on DERIVED numbers:** `D0 = K − B0 − Z0 = 7.7209 − 3.9e-5 − 2.4e-8 ≈ 7.7209`; `Z0` is
~3e-9 of `K` → **NO fragile cancellation** (vs the old effective-closure ~1.6e7 knife-edge). The soft spots
(`Z4,N4`) feed only the **held-out** `P4/R_pole`; the `R_norm` calibration rides on the tight `N0` (0.4%) + `K`.
**Forward source** = canonical decision-05 **D3 Fréchet** over the B2a BdG-mode response (NOT a B1-χ overlap; the
B2a overlaps `c_j` already carry the χ coupling). The B2a **rest** background has `A_r0 ≡ A_w0 ≡ 0` (expected — no
rest current), so this calibration transfer runs through scalar `A_00` + kinematic `j_E` + `q·δρ`; the
vector-current coupling is held-out non-rest/excited-defect coverage.
**Review (3 rounds, distrust-all-clean):** 2 transliteration-fidelity audits → operator FAITHFUL term-by-term and
UNTOUCHED through both remediations; `Γ_port` correctly resolves the decision-05 factor-4 flag and cancels between the
DtN build and the `N` normalization. Adversarial rounds 1–2 caught **real gate-theater** (31.5% unconverged + 0.75
can't-fail gate; a transliteration masquerading as the "dual engine"; a dead-conditional / algebraic-identity
basis-invariance gate) → remediations fixed each. Round 3 verified all fixes genuine: per-coefficient convergence
gate (max-over-coefficients FORBIDDEN; grow-one/shrink-another negative control fails), the 4th-order independent
engine actually wired in (different stencils, same full sweeps), the basis gate genuinely basis-dependent (92% branch
diff; fixed-port moves 9.88% while invariant `Z` is machine-zero). Target-blind CLEAN (no `R_norm`/…/`mt15_05`).
N-fit reconditioned (value-preserving column equilibration; cond 5.5e6→33; `N0/N2/N4` each above the fit floor).
**Lesson reinforced:** the standing fidelity audit verifies the operator FORM; the **adversarial pass is what catches
can't-fail/tautological gates and unconverged-shipped-as-converged** — needed 3 rounds here.

## NEXT STEPS (in order)
1. ~~§5b clarification note to `decisions/11`~~ — DONE.
2. ~~Scaffold + run B2a (`pathA_09`/`pathA_10`); review (fidelity + adversarial); commit~~ — DONE + committed.
3. ~~B2b — Maxwell transfer (`pathA_11`/`pathA_12`/`pathA_13`); basis-invariant Green/self-energy on the Path-A `A0`
   → derived `{Z0,Z2,Z4,N0,N2,N4}`; 3-round review; commit~~ — **DONE + committed.**
   **NOTE (user gate 2026-06-18): user paused AFTER the B2b commit to review the numbers themselves before B2c is
   launched — do NOT start B2c until the user gives the go.**
4. **B2c — integrate + calibrate (task #70).** Feed the full derived bundle `{K,M,B_n,Z_n,N_n}` (B1 `χ,K,M` + B2a
   `B_n` + B2b `Z_n,N_n`) → `patha_extraction` → the unique `R_norm(τ)=0` deterministic root-find on the stable-side
   `D0>0`. §J: propagate the B2a + B2b per-coefficient error bars (esp. `Z4/N4 ~2–5%`) into held-out
   `R_pole/P2/P4`; report `τ*`, naturalness (`|ln τ*|`, `K/(B0+Z0)` cancellation ratio), held-out surplus. **LEAVES
   target-blind** (the `R_norm` anchor is the calibration target; do NOT peek at the held-out `R_pole/P2/P4` targets).
Discipline unchanged: Codex codes / Claude reviews; dual-engine + fidelity + adversarial per chunk; commit per
validated chunk; build target-blind in the held-out sense (don't peek at `R_pole/P2/P4` targets; `R_norm` anchor is
the calibration target, allowed).
