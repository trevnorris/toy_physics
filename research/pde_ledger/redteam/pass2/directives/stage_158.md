---
unit_id: 158
batch: IV.6
created_at: 2026-06-08T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-08T11:50:02-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 158

Apply the finding below. After applying, append an `## Applied: F1` block under it with:
`files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a
question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Do NOT change the
set of asserted identities or their expected-zero targets, and do NOT change any printed numeric value
— only the Mathematica **derivation route** so the second engine is genuinely independent of the SymPy
one. You DESIGN the independent route (the directive states the requirement and the acceptance
criteria, not a line-by-line recipe).

After editing, RUN the script
(`math -script /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl`)
and iterate under `timeout 600` until it exits 0 with every `expectZero` printing `PASS:` and the final
`Stage 158 Mathematica audit passed.` line. A timeout (exit 124) is a FAILURE — reformulate, never raise
the cap.

Do NOT touch paper.tex, notes/, or the SymPy script.

## F1 — mathematica_transliteration (orchestrator-elevated; user-authorized re-author 2026-06-08)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl`

**Issue:** The `.wl` is a near-line-by-line port of the SymPy `.py`. For the two load-bearing checks it
uses the SAME mechanism on the SAME shifted closed forms as SymPy:
- "linear delta R law" — `.wl` `rLin = Normal[Series[rFun /. g -> gStar + dg, {dg,0,1}]]`
  (wl:33-34) mirrors `.py` `R_lin = sp.series(R.subs(g, g_star+dg), dg, 0, 2).removeO()` (py:40-41);
- "linear Delta_Q law" — `.wl` `chiLin = Normal[Series[chi, {eps,0,1}]]` (wl:74-75) mirrors `.py`
  `chi_lin = sp.series(chi, eps, 0, 2).removeO()` (py:90).
The gain/slope checks (`delta Mq`, `delta Pi`) likewise re-type the same `expand`-then-drop-cross-terms
construction, and the numerical-coefficient block re-types the same closed-form expressions. Both engines
are running the identical Taylor-expansion choreography, so the second engine echoes the first rather
than re-deriving — defeating the dual-engine independence guarantee (same defect class as 161; the
"no independent route exists for a pure-Taylor stage" framing is rejected — an independent route IS
feasible, see below).

**Required change (route only — keep every asserted identity, every expected-zero target, and every
printed numeric value byte-for-byte in VALUE):**

Make the Mathematica derivations of the linear laws algorithmically independent of "Series-expand the
same shifted closed form." Feasible independent routes (you choose; combine as needed):

1. **Symbolic linear laws via analytic differentiation at the base point** rather than series-expanding
   the shifted expression. E.g. obtain the `delta R` slope as `D[rFun, g] /. g -> gStar` (this yields
   `-1/Sqrt[1+r^2]` directly from the branch value `gStar-r = -Sqrt[1+r^2]/2`), and obtain the
   `delta Mq`, `delta Pi`, and `Delta_Q` (`chi`) linear coefficients as the relevant first partial
   derivatives evaluated at the base point (`D[Mq,{...}]`, `D[chi,eps]/.eps->0`), instead of
   `Series[...]` / `expand`-and-drop-cross-terms. Assert the SAME `expectZero` identities with the SAME
   targets.

2. **AND/OR an independent numerical finite-difference cross-check at the Stage-156 canonical point.**
   The stage already carries `rF1`, `Sigma0_can`, `S_can`, `T_can`. Verify each emitted numerical
   coefficient (`dR/dg`, `dMq/dSigma0`, `dMq/dg`, `dPi/dSigma0`, `dPi/dS`, `dPi/dg`, `dSigma0/dThat`,
   `dPi/dThat`) against a central finite-difference of the underlying closed form at that point (small
   step `h`, e.g. 1e-6 with adequate `WorkingPrecision`), asserting agreement to a stated tolerance.
   This is a genuinely independent (numerical, non-series) route AND it hardens the numeric deliverables
   that are currently only re-typed from the symbolic formula.

**Acceptance criteria (the verifier checks all of these):**
- The script exits 0; every `expectZero` prints `PASS:`; final `Stage 158 Mathematica audit passed.`
- The set of asserted identities and their expected-zero targets is UNCHANGED.
- Every printed numeric coefficient is unchanged in VALUE from the committed transcript (engine
  agreement with SymPy preserved): `dR/dg = -0.49021604438762603754…`, `dMq/dSigma0 = -1/4`,
  `dMq/dg = 2.28001126927792351405…`, `dPi/dSigma0 = 0.83240947108163457213…`,
  `dPi/dS = -1.16275838754221894078…`, `dPi/dg = 1.52843317823248362127…`,
  `dSigma0/dThat = 6.42981496203005499347…`, `dPi/dThat = 5.35223887169621835652…`.
- The `.wl` no longer derives the load-bearing `delta R` and `delta Q (chi)` linear laws by
  `Series`-expanding the same shifted closed form the SymPy script uses (the second engine is
  genuinely independent).

**Verification command:** the verifier runs `redteam exec-mathematica 158` and confirms (a) exit 0 with
all `PASS:` lines, (b) the printed coefficients above unchanged, (c) the derivation route is no longer a
Series-mirror of the SymPy `.py`.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl`
- summary: Replaced the Mathematica linear-law and numeric-coefficient derivations with analytic first-derivative evaluations at the base point while preserving the existing checks and printed coefficient values.
- deviation: none
