---
unit_id: 189
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 189

This is a CLEAN re-audited unit (0 findings, no source-code change). There is no Codex
directive to verify findings against; per the prompt, the verifier confirms (A) the
committed outputs carry the canonical banner and read PASS / `= 0` / `True`, and
(B) the audit's clean disposition holds on a fresh read (independence + 0 value-reconciliation
misalignments). Both confirmed.

## Per-finding outcomes

No findings in the original report (`findings_count: 0`). Nothing to classify as
resolved/partial/regressed. Disposition-confirmation below.

## Disposition confirmation

### (A) Committed outputs — banner + check status

**SymPy output** (`scripts/output/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.txt`):
- Canonical banner present: line 3 `STAGE 189 — TRANSFER-SHAPE / OUTGOING-PREFACTOR COMPILER`; ledger banner line 120 `STAGE 189 LEDGER`. No stale 172/206/240/239 label.
- Every residual check reads `= 0`: matrix residual (out 35-40, all-zero column vector), `rank = 2` (out 41), the three defect identities (out 42-44), direct-slope bridge (out 73), `Y(omega)`/`Pref(omega)`/`P0` compiler (out 84,94,95), weak-axisym slopes (out 100-102), outgoing expansion (out 107), normalization product + `Γ₅` (out 112-113), constant-prefactor branch (out 114-117). All `= 0`.
- The back-definition appears only as a printed form (out 55, `R_target selected-branch definition := ...`), NOT as a `= 0` assertion — demotion preserved in the committed output.

**Mathematica output** (`mathematica/output/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_mathematica_audit.txt`):
- Canonical banner present: line 3 `STAGE 189 — TRANSFER-SHAPE / OUTGOING-PREFACTOR COMPILER`; ledger banner line 75 `STAGE 189 MATHEMATICA LEDGER`.
- Every check emits an explicit `PASS:` line — 13 PASS lines (out 11,13,15,17,19,28,30,36,38,40,46,48,50,56,62,64,66,68,70,72), one per residual, each on a `= 0` / `{0,...}` residual. No FAIL, no non-zero residual.
- The back-definition appears only as a printed form (out 25), NOT fed to `expectZero` — demotion preserved.

### (B) Clean disposition holds on fresh read

- **Demoted back-definition NOT asserted.** Confirmed at source: `.py:101` `Rtarget_definition = sp.simplify(Lambda0*(1-epseta)/T2_direct)` is only `sp.pprint`ed (py:104-105), never passed to `expect_zero`. `.wl:105` `rTargetDefinition` is only `Print`ed (wl:107), never passed to `expectZero`. The iter-2 demotion held in both engines.
- **`Ξ₁` independent.** Built from raw input amplitudes on both sides: `.py:119-125` `eps1 = diff(epssplit,epsW)*epsW1 + diff(epssplit,deltaU)*deltaU1`, then `Xi1_closed = zetaZ − omegaW + 2*chi1/(1+chi0) + 2*eps1/(1−epssplit)`; `.wl:116-118` structurally identical (`eps1Bridge`, `xi1Closed`). It is supplied to the residual `direct_slope − ε λ_A Xi1_closed` (py:136 / wl:127-130), not back-derived from the slope. Non-tautological: a wrong amplitude combination would leave a nonzero residual.
- **Genuinely independent `.wl`.** Confirmed at source: the `.wl` recovers `C_obs→trf` as a finite-log Jacobian (`D[transferLogFunctions[[i]], obsVars[[j]]] /. zeroObs`, wl:65-73) and the defect slopes via exponential paths (`rTrPath = r0Obs*Exp[s*theta1]`, `D[Log[tSqPath], s] /. s→0`, wl:85-92), versus the `.py`'s hand-posited linear matrix + linearized-symbolic substitution. Different derivation routes. (The one section with shared variable choreography — Sec II `Xi1_closed` perturbation build — is independently cross-validated by Sec I's distinct route, so no transliteration finding; audit's BORDERLINE-but-INDEPENDENT call confirmed.)
- **Value reconciliation.** Audit checked 21 deliverable values against card+notes, 0 misaligned. Spot-confirmed the load-bearing constants emit correctly in both outputs: `Λ₀ = 27π²G c_s⁵/(20a⁵c⁵)` (sympy out 61-64, math out 26), `T_A² = Z_W(1+ρ)²/(Ω_W²(1−ε_W)²)` (sympy out 49-54, math out 24). MATCH.

## Exec log assessment

No fresh exec was run for this verification (clean unit, no source change). Disposition is
confirmed from committed source + committed outputs, per prompt scope.

**Output freshness:** confirmed. Output mtimes (Jun 9 14:04) are newer than script mtimes
(Jun 1 11:48) for both engines. Outputs are post-fix-era and regenerated.

## Material-change assessment

`material_change`: false. No source-code change to this unit; no derived result altered.
No downstream unit is affected.

## Side observations (non-blocking)

The card's `\stagefield{Verification}` line says "Mathematica audit: none yet" while a
passing `.wl` exists — a paper-prose lag, paper-side, out of red-team scope. Already noted
in the audit; non-blocking.

## Verdict justification

The committed outputs carry the canonical `STAGE 189` banner and every check reads PASS /
`= 0` (13 explicit PASS lines on the Mathematica side, all-`= 0` residuals on the SymPy
side; no FAIL anywhere). On a fresh source read the audit's clean disposition holds: the
demoted back-definition is printed-only (not asserted) in both engines, `Ξ₁`/`Xi1_closed`
is built from raw input amplitudes (not back-derived from the slope), and the `.wl` uses a
genuinely independent route (finite-log Jacobian + exponential-path slopes vs. the `.py`'s
linearized substitution). Value reconciliation is 0-misaligned. Outputs are fresh. Verdict: verified.
