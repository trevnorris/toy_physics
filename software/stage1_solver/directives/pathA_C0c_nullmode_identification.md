# Directive pathA_C0c — Identify the field-block near-null mode (symmetry/gauge vs genuine stiffness)

**Status:** READY (Claude-authored 2026-06-20; design-review = SOUND-WITH-FIXES → fixes applied [physics check PASSED:
`g_phase=iψ` confirmed correct, global phase UNFIXED in the coded residual; layout `psi_real,psi_imag,a0,ar,aw,r0,mu`];
confirm-pass = SOUND-WITH-FIXES → the one stale-phrase fix applied; user GATED "C0c: identify the mode" → executing.)
Follows `pathA_C0b` (verdict `DIAGNOSTIC_INCOMPLETE`, fidelity-verified): the τ≈0.029 wall is a PERSISTENT near-null
SUBSPACE in the PHYSICS FIELD BLOCK (true dense-SVD: ~5 modes, `v_min` field=1.0, `u_min` pde_rows=1.0 at both τ; NOT the
mass/μ border, NOT a fold). Step of option C, task #78.

**Date:** 2026-06-20
**Owner:** Codex (codes + iterates until exit 0; extends the C0b module/diagnostics). Claude reviews after.
**Trigger:** C0b localized the near-null space to the field block but `gauge_projection_status=not_available` (the
symmetry/gauge generators were never implemented). This directive IMPLEMENTS them and IDENTIFIES the mode.

## Why this is the decisive, cheap next step
A stationary gauged-GPE+Maxwell BVP that does not pin the global U(1) phase has an EXACT phase zero mode → its Jacobian is
singular BY CONSTRUCTION at every depth (exactly the depth-independent near-singularity C0b measured). If the field-block
near-null mode IS that gauge/phase mode (or another symmetry mode), the fix is a CHEAP gauge-fixing constraint / null-space
deflation — NOT a production solver, NOT pseudo-arclength. If instead the modes are NOT explained by symmetry generators,
it is genuine field-block stiffness (the harder question). The test is cheap: evaluate analytic generators at the converged
state and (a) check they are annihilated by the Jacobian, (b) check they span the measured near-null subspace. No new solve.

## Stance (carry the C0b Single-Arbiter discipline)
DIAGNOSIS ONLY. Do NOT alter the faithful PDE operators (`coupled_branch.py`/`operators.py`), frozen physics, or
`physical_export_permitted`. Do NOT implement the gauge-fix/deflation or a re-crawl here (that is the NEXT step, gated on
this result). The generators are evaluated on the EXISTING converged state + the EXISTING assembled Jacobian — no operator
edit. Reuse the C0b machinery (the true dense-SVD diagnostics + the saved assembled Jacobian matrices under
`runs/pathA_C0b_wall_diagnosis/matrices/`, from which the near-null singular vectors are RECOMPUTED here — do not assume the
vectors themselves were persisted). CPU; `timeout 600` per script (split if needed;
timeout → NOT_MEASURED); standalone `python3`; YAML/markdown human output, JSON only for machine artifacts; no commentary
`python3 -c`.

## Work items

### C0c-1 — construct the candidate symmetry/gauge generators at the converged state
The field layout (confirmed in design-review) is the 5 per-cell components `psi_real, psi_imag, a0, ar, aw` then `r0` (nw),
then `mu` (1). Build, at the converged shallow τ (verify it from `coupled_branch.py`, don't hardcode), the infinitesimal-
generator state-vectors:
  - **Global U(1) phase (the PRIME SUSPECT — confirmed an UNFIXED EXACT symmetry of the coded residual, even with the
    Maxwell `ξ⁻¹∇·A` gauge-control term):** `g_phase = (−psi_imag, psi_real, 0, 0, 0, 0, 0)` (acts only on the ψ lanes;
    zero on a0/ar/aw/r0/μ).
  - **r-translation, w-translation:** `g = ∂_r(state)`, `∂_w(state)` (finite-difference on the grid). NOTE: V_conf + the
    wall geometry likely BREAK these — they are PROBES, not assumed symmetries; report regardless.
  - **Dilation/scaling:** an `(r∂_r)`-type generator. Also a probe (likely broken by the geometry/EOS scale).
  - **Maxwell / residual-gauge probe:** because a field-block near-null mode could live in the `a0/ar/aw` lanes rather than
    the ψ lanes, ALSO probe a residual-gauge / Maxwell-sector direction (e.g. a `∇λ`-type pure-gauge `A` shift consistent
    with the coded gauge-control term), and ALWAYS report the per-lane energy split of each measured null mode (see C0c-3).
State explicitly which generators are EXACT symmetries of the actual BVP vs broken-symmetry PROBES.

### C0c-2 — NULL-MODE test (is the generator annihilated by the Jacobian?) — CONVERGED states only for the exact gate
For each generator `g`, compute the relative annihilation `‖J·g‖ / (σ_max·‖g‖)` using the same assembled bordered Jacobian
C0b used (and a JVP of the true residual as a cross-check). **IMPORTANT (equivariance):** an exact symmetry gives `J·g = 0`
ONLY at a converged ROOT; at a NON-converged state `J·g = G·F(x)` (the generator acting on the nonzero residual), which is
generally nonzero. So apply the exact null gate (`≤ 1e-8`) ONLY at the CONVERGED shallow τ (residual ~3.4e-10). At the
NON-converged deepest τ (residual ~5.3e-5), either mark the annihilation `NOT_MEASURED`, OR verify the equivariance
identity `J·g_phase ≈ phase_action(F(x))` instead (which, if it holds, still confirms g is the symmetry direction). Report
the value + `MEASURED/NOT_MEASURED` + which test was used, per generator per state.

### C0c-3 — OVERLAP test (do the generators SPAN the measured near-null subspace?)
RECOMPUTE the near-null right singular vectors by SVD of the SAVED assembled Jacobian matrices
(`runs/pathA_C0b_wall_diagnosis/matrices/*.npz` — C0b persisted the matrices, not necessarily the vectors). Take the
smallest ~5 (the field-localized cluster). For each, report (i) the normalized overlap `|⟨v_mode, ĝ⟩|` with each unit
generator `ĝ`; (ii) the RESIDUAL fraction of `v_mode` NOT captured by the span of the generators (`1 − ‖P_span v_mode‖²`);
and (iii) the PER-LANE energy split over `psi_real, psi_imag, a0, ar, aw` (so a ψ-lane mode = phase candidate vs an
`a0/ar/aw`-lane mode = Maxwell/residual-gauge candidate is visible). Guard against a hardcoded result: include a unit test
that plants a known generator as a singular vector and confirms the overlap reports ~1 (mirroring the C0b planted-null-mode
test).

### C0c-4 — σ validation (retire the stencil-truncation caveat)
At one τ, recompute σ_min from a DENSE / full-JVP Jacobian (not the colored stencil-radius-3 assembly) and confirm it
matches the C0b value to tolerance (or report the discrepancy + which is trusted). This closes the fidelity-flagged caveat
that C0b's SVD was of the assembled stencil approximation, not the exact Fréchet derivative.

### C0c-5 — THE VERDICT (falsifiable, WHOLE-cluster — one mode does not explain the rest)
Classify EVERY mode in the measured near-null cluster (all ~5, not just the smallest): each is `GAUGE_PHASE` /
`TRANSLATION` / `DILATION` / `MAXWELL_RESIDUAL_GAUGE` / `UNEXPLAINED_STIFFNESS`, by the EXACT gate — a mode is a given
symmetry mode iff its generator passes C0c-2 (`‖J·g‖/(σ_max‖g‖) ≤ 1e-8` at the converged τ, or the equivariance identity)
AND the C0c-3 overlap `≥ 0.9`. The overall verdict must account for the WHOLE cluster (identifying ONE phase mode does NOT
explain away the remaining 4):
  - **SYMMETRY_MODE_IDENTIFIED** — the ENTIRE near-null cluster is symmetry/gauge-explained (every mode classified to a
    generator, residual fractions all small): recommend the SPECIFIC fix (pin the phase / deflate the identified null
    space) for the NEXT step; state exactly which mode(s) and the constraint.
  - **MIXED** — some modes symmetry-explained, a residual subspace NOT (e.g. phase explains mode 0 but modes 1–4 remain):
    enumerate BOTH the explained modes (+ their fix) and the unexplained residual subspace (with its per-lane split). This
    is the likely outcome if only the phase is unfixed but other small modes are physical/stiffness.
  - **GENUINE_STIFFNESS** — NO mode is symmetry-explained (residual fractions dominant, generators not annihilated): real
    field-block ill-conditioning → the harder solver question stands. Do NOT call a Maxwell-lane-dominant unexplained mode
    "stiffness" without reporting its per-lane split (it may be a residual-gauge mode).
  - **DIAGNOSTIC_INCOMPLETE** — required evidence NOT_MEASURED (a generator couldn't be built, J·g couldn't be formed, the
    annihilation gate only available at a non-converged state, etc.).

## Acceptance criteria (PASS/FAIL; exit-0 NECESSARY not SUFFICIENT)
1. Generators built (phase + translation/dilation probes + a Maxwell/residual-gauge probe), each labeled EXACT-symmetry vs
   PROBE, mapped to the verified `psi_real,psi_imag,a0,ar,aw,r0,mu` layout (phase = `(−psi_imag,psi_real,0,0,0,0,0)`,
   confirmed against `coupled_branch.py`, not guessed).
2. C0c-2 `‖J·g‖/(σ_max‖g‖)` reported per generator (MEASURED/NOT_MEASURED), computed against the real Jacobian + a JVP
   cross-check; the exact `≤1e-8` null gate applied ONLY at the CONVERGED τ, and the NON-converged τ either NOT_MEASURED or
   checked via the equivariance identity `J·g_phase ≈ phase_action(F(x))`.
3. C0c-3 overlaps of the near-null singular vectors (RECOMPUTED by SVD of the saved `runs/.../matrices/*.npz`, not assumed
   saved) with each generator + the unexplained residual fraction + the per-lane (`psi_re/psi_im/a0/ar/aw`) energy split of
   each mode; genuineness protected by a planted-generator unit test (overlap→~1 when the planted vector IS the generator).
4. C0c-4 dense/full-JVP σ_min validation vs the C0b stencil value (match-to-tolerance or a flagged discrepancy).
5. Exactly one C0c-5 verdict (or DIAGNOSTIC_INCOMPLETE) with the falsifiable numeric `verdict_support` (the ‖J·g‖ values,
   overlaps, residual fractions, thresholds) + the recommended next-step fix if a symmetry mode is identified.
6. Faithful operators + frozen physics + export guard untouched (diff); no gauge-fix/deflation/re-crawl IMPLEMENTED here;
   chunk-1a/1b/1c gates still pass; report + machine JSON; tests added (no hardcoded overlaps/annihilations).

**Fail conditions:** hardcoded overlaps/annihilation values; claiming a symmetry mode without BOTH C0c-2 and C0c-3 passing
its gate; altering operators/frozen physics/export guard; implementing the fix or a re-crawl (out of scope); masking a
NOT_MEASURED as a result; raising the timeout cap.

## Out of scope
Implementing the gauge-fix / null-space deflation; the re-crawl with the fix; the full-budget+backtracking crawl
(SPIKE-vs-persistent — only meaningful AFTER the mode is addressed); the production-solver decision; `pathA_22`.

## Review (orchestrator, after Codex)
Fidelity agent: are the generators correct for the actual layout (esp. the phase generator's ψ-component mapping)? are
‖J·g‖ and the overlaps GENUINELY computed (planted-generator test) — not hardcoded? is the dense-σ validation real?
Adversarial agent: is the symmetry-identification SOUND (both gates met) or over-claimed; if GENUINE_STIFFNESS, is the
residual genuinely unexplained? Diff-check faithful operators untouched. Then gate the NEXT step (add the gauge-fix/
deflation + re-crawl IF a symmetry mode; else the solver question).
