---
unit_id: 137
batch: IV.4
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T17:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 137

This is a REMEDIATION verify. The authoritative findings are the Codex recovery review
(`redteam/codex_reviews/stage_137.md`, verdict FINDINGS: R1/R2/R3), which the rewritten
directive (`redteam/directives/stage_137.md`) re-encodes as F1=R3 (matrix-Schur
reconstruction), F2=R1 (de-tautologize static limit), F3=R2 (nonzero-S_q outlet). F4 (banner
relabel) was already RESOLVED pre-pass. This file supersedes the prior 2026-05-27 verdict, which
verified the orchestrator-direct (Codex-bypassed) pass that the recovery review then flagged as
non-adversarial. I checked each fix for genuine non-tautology against the load-bearing criteria,
the edited scripts, the exec logs, the regenerated saved outputs, and the owner stage 114.

## Per-finding outcomes

### F1 (=R3) — insufficient_verification: matrix-Schur reconstruction of rho_c, sigma_c

**Classification:** resolved

**What changed:**
SymPy `scripts/...stage137...sympy_audit.py:12-28` inserts the matrix-Schur block AFTER the
hand-assigned `rho_c, sigma_c` (lines 9-10) and BEFORE `Ms` (line 30):
`M_core = sp.Matrix([[Ks, lam], [lam, -Kq*D_sch]])`, `v_coup = sp.Matrix([gs, gq])`,
`delta_Lambda_schur = sp.apart((v_coup.T * M_core.inv() * v_coup)[0], D_sch)`,
`rho_c_schur = limit(D_sch -> oo)`, `sigma_c_schur = rho_c_schur - delta_Lambda_schur|_{D=1}`,
then asserts `rho_c - rho_c_schur == 0` and `sigma_c - sigma_c_schur == 0`.
Mathematica `mathematica/...stage137...audit.wl:40-49` mirrors with
`mCore = {{kS, lam}, {lam, -kQ*dSch}}`, `vCoup = {gS, gQ}`,
`deltaLambdaSchur = Apart[vCoup . Inverse[mCore] . vCoup, dSch]`,
`rhoCSchur = Limit[..., dSch -> Infinity]`, `sigmaCSchur = rhoCSchur - (deltaLambdaSchur /. dSch->1)`,
and two `expectZero` asserts. Declarations added to `Clear[...]` (lines 28-30) and `$Assumptions`
(lines 31-35: `dSch, kappa0, gamma0` with positivity).

**Assessment:**
Genuinely non-tautological. The matrix entries are the PHYSICAL primitives `K_s, K_q, lam` and
the source vector is `(g_s, g_q)` — none of which are functions of the hand-assigned
`rho_c, sigma_c`. The reconstruction DERIVES the residues by inverting the stiffness matrix and
then asserts they equal the hand-typed values, so a wrong factor/sign in the hand-typed
`sigma_c` would make `sp.simplify(sigma_c - sigma_c_schur) != 0` and fail. I confirmed this is
the identical primitive the already-verified owner stage 114 uses
(`scripts/moving_throat_pde_stage114_concrete_core_schur_sympy_audit.py:27-30`:
`M = sp.Matrix([[K_s, lam],[lam, -K_q*D]])`, `c = sp.Matrix([g_s, g_q])`,
`sp.apart((c.T*M.inv()*c)[0], D)`) — confirming the project-canonical independent route, NOT
back-built from rho_c/sigma_c. Two engines, structurally distinct primitives
(`M_core.inv()` + `sp.limit` vs `Inverse[mCore]` + `Limit`). Exec logs confirm the new SymPy
line "rho_c, sigma_c reproduced from explicit two-channel core Schur complement (M_core)." and
the two new Mathematica PASS lines for D->Infinity and static D=1.

### F2 (=R1) — tautological_check: de-tautologize the static-limit / susceptibility check

**Classification:** resolved

**What changed:**
SymPy `...sympy_audit.py:49-72` replaces the old `static_limit - (rho_c - sigma_c)` X-X block.
New block defines `r_c = lam**2/(Ks*Kq)`, `kappa_c = kappa0/(1+r_c)`, `gamma_c = gamma0/(1+r_c)`,
`D_W_bare = 1 - kappa0*z_var**2 - I*gamma0*z_var**5`, the matrix source
`delta_Lambda_matrix = delta_Lambda_schur.subs(D_sch, D_W_bare)`, the reduced envelope
`delta_Lambda_reduced = rho_c - sigma_c/(1 - kappa_c*z**2 - I*gamma_c*z**5)`, asserts they are
equal, then a static specialization `static_limit = sp.limit(delta_Lambda_matrix, z_var, 0)`
asserted equal to `rho_c_schur - sigma_c_schur` (matrix-derived). Mathematica `...audit.wl:65-77`
mirrors with `rC, kappaC, gammaC, dWbare`, `deltaLambdaMatrix = deltaLambdaSchur /. dSch -> dWbare`,
`deltaLambdaReduced`, and `staticLimit = Normal[Series[deltaLambdaMatrix, {zVar, 0, 0}]]` (Series,
the prescribed engine divergence from SymPy's `sp.limit`).

**Assessment:**
The old `assert sp.simplify(static_limit - (rho_c - sigma_c)) == 0` and
`expectZero["Schur static limit equals rho_c - sigma_c", staticLimit - (rhoC - sigmaC)]` are GONE
(grep-confirmed absent; "Schur static limit equals rho_c - sigma_c" appears in neither
transcript). The comparison is now matrix-route (`delta_Lambda_schur` on `D_W_bare(z)`) vs
reduced-envelope (built from hand-assigned `rho_c, sigma_c` and the Stage-97/114 coefficient
maps). The residual is zero ONLY if the hand-typed residues AND the coefficient maps
`kappa_c, gamma_c, r_c` are all correct, so it is falsifiable. The static assert ties to the F1
matrix-derived `rho_c_schur - sigma_c_schur`, not to the hand-assigned pair. Confirmed coefficient
maps match the directive's prescription (`kappa_c = kappa0/(1+r_c)`, `gamma_c = gamma0/(1+r_c)`,
`r_c = lam^2/(K_s K_q)`) and owner stage 114 lines 35-39. Exec logs show the two new SymPy lines
and the two new Mathematica PASS lines.

### F3 (=R2) — tautological_check: nonzero-S_q outlet consistency

**Classification:** resolved

**What changed:**
SymPy `...sympy_audit.py:74-91` replaces the old `outlet_residual.subs(Sq_var, 0).subs(Pi_var, Ms)`
block. New block keeps `Sq_var` symbolic/nonzero: `Pi_map = Ms + Mq*Sq_var`,
`mixed_contribution = sp.simplify(Pi_map - Ms)`, `Mq_from_schur = -L*sigma_c_schur/Theta`,
asserts `mixed_contribution - Mq_from_schur*Sq_var == 0`, then a non-vacuity guard
`assert sp.simplify(mixed_contribution - (-Mq_from_schur)*Sq_var) != 0`. Mathematica
`...audit.wl:79-91` mirrors with `piMap`, `mixedContribution`, `mQFromSchur = -lM*sigmaCSchur/thetaSigma`,
an `expectZero`, and an `If[TrueQ[vacuityResidual === 0], fail[...], pass["...non-vacuous..."]]`.

**Assessment:**
The `S_q = 0` substitution is GONE in both engines (grep confirms no `subs(Sq_var, 0)` in SymPy
and no `sqVar -> 0` in Mathematica; the old "Outlet consistency ... at S_q = 0" PASS line is
absent from both transcripts). The mixed `M_q*S_q` term is now isolated with `S_q` nonzero and
compared against `-L*sigma_c_schur*S_q/Theta` where `sigma_c_schur` is rebuilt from the F1
matrix Schur (traces to `Inverse(M_core)`), NOT from `Mq` itself — so a flipped sign or wrong
factor in `M_q` makes the residual nonzero. The sign `-L*sigma_c/Theta` matches the notes:40-48
relation `Pi = (L/Theta)[rho_c U_s'(0) - sigma_c U_q'(0)]` (negative on the q-channel) as cited
and confirmed in the directive. The explicit non-vacuity guard proves `+M_q` and `-M_q` give
different residuals (SymPy `!= 0` assert; Mathematica `If` branch that PASSes only when the
flipped residual is nonzero) — so the check cannot be vacuously satisfied. Exec logs show both
new SymPy lines and both new Mathematica PASS lines, including
"outlet consistency non-vacuous (sign of M_q is exercised)".

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- "rho_c, sigma_c reproduced from explicit two-channel core Schur complement (M_core)." (F1)
- "Reduced core susceptibility matches the matrix-Schur source (full z dependence)." +
  "Static core residue matches rho_c_schur - sigma_c_schur from M_core." (F2)
- "Outlet mixed channel M_q*S_q matches the matrix-Schur reconstruction (S_q != 0)." +
  "Outlet consistency (paper Checks item 1) verified with nonzero S_q." (F3)
- No "Schur-complement static limit matches rho_c - sigma_c." and no "Outlet consistency reduces
  to Pi = M_s at S_q = 0." — the old tautological lines are gone.

**Mathematica:** exit=0, 9 PASS / 0 FAIL. Notable lines:
- Banner reads "STAGE 137 — EXPLICIT CORE-TO-MOUTH GAIN MAP" (no STAGE 120; .wl line 26 confirmed).
- 9 PASS: rho_c@D->Infinity, sigma_c@static-D=1 (F1); M_s/M_q paper card (kept); reduced core
  susceptibility, static core residue (F2); outlet mixed channel + outlet non-vacuous (F3);
  sigma_c r_c-form equivalence (kept).
- One benign warning `Limit::alimv` (Mathematica ignored the `dSch>0` assumption during
  `dSch -> Infinity`); the immediately following residual = 0 and PASS confirm the limit
  evaluated correctly. Non-blocking.
- No old "Schur static limit equals rho_c - sigma_c" / "Outlet consistency Pi = M_s at S_q = 0"
  PASS lines.

**Output freshness:** confirmed regenerated post-fix. Saved outputs live at
`scripts/output/...sympy_audit.txt` (mtime 2026-05-29 16:49:02) and
`mathematica/output/...mathematica_audit.txt` (mtime 2026-05-29 16:49:02), both newer than the
scripts (`.py` 16:41:48, `.wl` 16:42:16). Saved-output contents match the exec logs line-for-line
(modulo the log header), including all new PASS lines and the absence of the old tautological ones.

## Material-change assessment

`material_change`: false.

The verification surface was strengthened but no derived value/constant changed. `rho_c =
g_s^2/K_s`, `sigma_c = (K_s g_q - lam g_s)^2/[K_s(K_s K_q + lam^2)]`, `M_s = L*rho_c/Theta`,
`M_q = -L*sigma_c/Theta` are exactly the same expressions as before (transcript lines 7-10 / WL
16-19) — they are now DERIVED from the physical core matrix and cross-asserted rather than only
hand-typed. The fixes add independent anchors (matrix Schur, full-susceptibility comparison,
nonzero-S_q outlet, non-vacuity guard); none alter the gain-map outputs downstream units consume.
No downstream re-audit is warranted on account of stage 137.

## Side observations (non-blocking)

- The `Limit::alimv` warning in the Mathematica transcript is the standard Mathematica notice
  when assumptions reference the limit variable; the residual evaluates to 0 and PASSes, so it is
  cosmetic. Not a finding.
- SymPy declares `rc = sp.symbols('r_c', ...)` at line 7 (unused after F2 introduced the assigned
  `r_c = lam**2/(Ks*Kq)` at line 55, which shadows it). Harmless dead symbol; the directive
  explicitly said not to chase unused-variable cleanup. Not blocking.
- This verdict file overwrote a stale 2026-05-27 verdict (findings_total: 4, F3 marked
  blocked_legitimate) that verified the pre-remediation orchestrator-direct edits. That older
  state no longer matches the scripts and is fully superseded here.

## Verdict justification

All three findings (F1=R3 matrix-Schur reconstruction, F2=R1 de-tautologized static limit,
F3=R2 nonzero-S_q outlet) are resolved with genuinely non-tautological checks: the residues are
now derived by inverting the physical core matrix (entries K_s,K_q,lam; source g_s,g_q — never
functions of rho_c/sigma_c, matching owner stage 114), the static-limit check is matrix-route vs
reduced-envelope, and the outlet check keeps S_q nonzero with a sign-discriminating non-vacuity
guard. The old X-X static assert and the S_q=0 substitution are both gone from both engines.
Both exec logs exit 0 (Mathematica 9 PASS / 0 FAIL), saved outputs are freshly regenerated, the
banner reads STAGE 137 with no residual STAGE 120, and no derived value changed
(material_change=false). Verdict: verified.
