---
unit_id: 114
batch: IV.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-06T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 114

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage114_concrete_core_schur_mathematica_audit.wl:33-43`. The built-in matrix-inverse Schur route was removed. HEAD had (per diff.patch):
```
m = {{kS, lam}, {lam, -kQ*dSym}};
c = {gS, gQ};
deltaD = FullSimplify[Apart[c.Inverse[m].c, dSym], ...];
```
It is replaced by an explicit scalar elimination of the 2×2 core system `M.(s,q)=(g_s,g_q)`:
```
coreEquations = {kS*sCore + lam*qCore == gS,
                 lam*sCore - kQ*dSym*qCore == gQ};
sFromFirst  = First[Solve[coreEquations[[1]], sCore]];      (* s = (gS - lam q)/kS  *)
qFromSecond = First[Solve[coreEquations[[2]] /. sFromFirst, qCore]];  (* eliminate s, solve for q *)
sEliminated = sCore /. sFromFirst /. qFromSecond;
qEliminated = qCore /. qFromSecond;
deltaD = FullSimplify[Apart[gS*sEliminated + gQ*qEliminated, dSym], ...];  (* mouth feedback gS*s + gQ*q, u=1 *)
```
No other lines changed: `rhoC`/`rC`/`sigmaTilde`/`sigmaC`/`kappaC`/`gammaC` (wl L46-51), `targetD` (L53), both `expectZero` asserts (L54, L59), and the four printed identifications (L63-66) are untouched.

**Assessment:**
Correct and the key check passes.

1. **Genuine independence.** The new route reaches `delta_Lambda(D)` by EXPLICIT GAUSSIAN ELIMINATION / back-substitution via `Solve` on the two scalar core equations — solve eq.1 for `sCore`, substitute into eq.2 and `Solve` for `qCore`, back-substitute, then form the mouth-feedback response `gS*s + gQ*q`. This is structurally distinct from the `.py`'s `delta_D = apart((c.T * M.inv() * c)[0], D)` (`py:30`), which forms the response via the built-in 2×2 matrix inverse `M.inv()`. The two engines now derive the Schur complement by genuinely different mechanisms (scalar elimination vs. matrix-inverse bilinear form), so a transcription error in the shared hand-written `targetD`/`targetZ` would no longer be masked by an identical computation path on the left-hand side.
   - `Inverse[m]` is GONE: `grep -nE 'Inverse|LinearSolve|PseudoInverse|\.inv'` over the `.wl` returns no matches (exit 1). It was NOT replaced by a thin inverse-wrapper — `LinearSolve`/`PseudoInverse` are absent; `Solve[]` on explicit scalar equations is the route. The matrix object `m`/`c` no longer exists in the file.

2. **No regression.** Both asserts still pass (Mathematica log L11-14: `Schur form identity = 0` / `PASS`, `low-frequency normalized outlet identity = 0` / `PASS`). The `rhoC`/`rC`/`sigmaTilde`/`sigmaC`/`kappaC`/`gammaC` identifications are emitted unchanged (the def block at wl L46-51 is byte-for-byte the pre-edit block; only the `deltaD`-construction lines above it changed). The asserts remain non-tautological: `deltaD` is now built from the eliminated `(s,q)`, while `targetD`/`targetZ` are independently hand-written closed forms, so a wrong target leaves a nonzero rational residual.

3. **Byte-identical output.** `git diff HEAD` on the committed transcript `mathematica/output/...stage114...txt` is EMPTY, and the file is absent from `git status --short` — byte-identical to HEAD, a method-only change. Deliverable lines unchanged: `delta_Lambda(D) = (dSym*gS^2*kQ - gQ^2*kS + 2*gQ*gS*lam)/(dSym*kQ*kS + lam^2)`, `rho_c = gS^2/kS`, `sigma_c = (gQ*kS - gS*lam)^2/(kS*(kQ*kS + lam^2))`, `kappa_c = (kappa0*kQ*kS)/(kQ*kS + lam^2)`, `gamma_c = (gamma0*kQ*kS)/(kQ*kS + lam^2)`.

No collateral edits beyond the authorized re-author: only the `.wl` is in the dirty set among stage-114 files; the `.py` reference engine is untouched; no numbering/prose touched.

## Exec log assessment

**SymPy:** exit=0. Reference engine unchanged. Notable lines: `delta_Lambda(D) = gₛ²/Kₛ − (Kₛ·g_q − gₛ·lam)²/[Kₛ(D·K_q·Kₛ + lam²)]`; `Schur form identity = 0`; `low-frequency normalized outlet identity = 0`.

**Mathematica:** exit=0 (re-authored engine). Notable lines:
- `delta_Lambda(D) = (dSym*gS^2*kQ - gQ^2*kS + 2*gQ*gS*lam)/(dSym*kQ*kS + lam^2)` — same value as before the edit and equals the SymPy form over a common denominator.
- `Schur form identity = 0` / `PASS: Schur form identity`
- `low-frequency normalized outlet identity = 0` / `PASS: low-frequency normalized outlet identity`
- `Stage 114 Mathematica audit passed.` / `# exit_code: 0`

**Output freshness:** The committed `.txt` is byte-identical to HEAD (intended — method-only change). The orchestrator re-ran both engines (logs dated 2026-06-06T01:55–01:56) at exit 0, confirming the re-authored route reproduces the same transcript. Freshness is established by the orchestrator re-run rather than mtime divergence, which is correct for a byte-identical method-only edit.

## Material-change assessment

`material_change`: false.

The edit changes only the derivation METHOD inside the Mathematica engine. Every emitted value (`delta_Lambda(D)`, `rho_c`, `sigma_c`, `kappa_c`, `gamma_c`, and both identity residuals) is byte-identical to HEAD. No derived result that downstream units (e.g. stage 137, which mirrors this stage) could depend on has changed. No `upstream_stale` propagation is warranted on numerical grounds.

## Side observations (non-blocking)

- The `Solve` calls operate on linear scalar equations with nonzero leading coefficients (`kS>0` assumed; the `q` coefficient `-kQ*dSym` is symbolic but `Solve` returns the general rational solution), so `First[Solve[...]]` is safe — a unique solution branch, no spurious `ConditionalExpression` reached the asserts (output shows clean `= 0`). Non-blocking.
- This unit is now eligible for the Independent-Mirror Set (same class as the IV.1-100 authorized re-author): distinct primitive (`Solve`-elimination vs. matrix inverse), distinct intermediate structure (`sCore`/`qCore` scalar unknowns vs. `M.inv()` bilinear form), byte-identical deliverables. Orchestrator's call.

## Verdict justification

The single finding F1 (mathematica_transliteration, medium) is fully resolved. The `.wl` no longer uses `Inverse[m]` or any thin inverse-wrapper (`LinearSolve`/`PseudoInverse` absent; grep exit 1); it reaches the Schur complement by an explicit `Solve`-based scalar elimination of the 2×2 core system, a genuinely distinct mechanism from the `.py`'s built-in matrix inverse. Both identity asserts still pass at exit 0 with the `rho_c`/`sigma_c`/`kappa_c`/`gamma_c`/`r_c` identifications emitted unchanged, and the committed Mathematica output is byte-identical to HEAD (method-only). All acceptance criteria are met.
