---
unit_id: 182
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-30T02:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 182

## Per-finding outcomes

### F1 — insufficient_verification (support-blindness checks were vacuous)

**Classification:** resolved

**What changed:**
- SymPy `scripts/...stage182...sympy_audit.py:138-146`: the six `sp.diff(Xi1_slip/R1_slip/Theta1_fact, lamphi1/kphi)` derivative `expect_zero` calls (former L139-144) are replaced by a loop over the three **microscopic-log direct forms** `Xi1_direct`, `R1_direct`, `Theta1_direct` computing `support_syms & form.free_symbols` and raising `AssertionError` if the intersection is non-empty.
- Mathematica `mathematica/...stage182...audit.wl:132-146`: the six `D[..., lamphi1/kphi]` `expectZero` calls (former L107-112) are replaced by a `Module`/`Do` over `{xi1Direct, r1Direct, theta1Direct}` asserting `FreeQ[form, lamphi1] && FreeQ[form, kphi]`, with `fail` on violation.
- Banner fold-in: sympy L34 and math L26 now read `STAGE 182` (were `STAGE 165`); confirmed in both exec logs (line 8).

**Assessment:**
Correct and non-vacuous. The original check was vacuous because it differentiated the abstract-Sigma forms, which never contain `lamphi1`/`kphi` regardless of the physics. The replacement asserts structural absence on the **direct** forms — the forms that are actually built from microscopic logs (`zetaZ`, `chi1`, `eps1`, ...). This genuinely fails if a support drift were ever wired into the defect construction: were `lamphi1` or `kphi` introduced into any of `zetaZ/omegaW/chi1/eta1/varepsW/deltaU1/eps1`, it would propagate into `Xi1_direct.free_symbols` and the intersection would be non-empty (SymPy) / `FreeQ` would be `False` (Mathematica), tripping the assertion. SymPy log L94-96 shows `leakage = set()` for all three; Mathematica log L79-81 shows `PASS ... support-blind` for all three. This is the honest in-scope check the directive prescribed; it does not attempt to reconstruct the upstream Stage-249 zeta-cancellation. The directive's recommended SymPy/Mathematica replacement code was applied verbatim with no deviation.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/...stage182...audit.wl:60-90` (four-slippage law) and `:114-130` (tracking/nontracking split) are re-derived rather than re-typed:
- `xi1Direct` (L61) is still built from microscopic-log expressions (`zetaZ - omegaW + 2 chi1/(1+chi0) + 2 eps1/(1-eps)`).
- `sigmaEqns` (L62-68) are the five Σ definitions in terms of logs (notes §3 physical input). `sigmaSolve = First[Solve[sigmaEqns, {mu1, gam1, kEta, kW, tau1}]]` (L69) gauge-fixes five of the eight logs — the consult-sanctioned linear Solve route, NOT a `Block`.
- `xi1DirectInSigmas` (L71-77) substitutes that solution into `xi1Direct` and `Collect`s on the Σ symbols. L78-81 asserts via `FreeQ[..., Alternatives @@ logSyms]` that **all eight logs are eliminated** — a real fail mode (leftover free logs).
- L82-87 extract each coefficient via `Coefficient[Expand[xi1DirectInSigmas], sigma...]` and `expectZero` against the paper's claimed coefficient (e.g. `2*chi0/(1+chi0)`, `2*epsW*(11+9*deltaU)/(11*(1-eps)*(1+deltaU))`). The former hand-typed `xi1Slip` literal (diff old L32-37) is **deleted**; `xi1Slip` is now assigned FROM the Mathematica-computed `xi1DirectInSigmas` (L88).
- The split (L115-122) is derived by substituting `sigmaChi -> (sigmaTr - (1+chi0)*sigmaDel)/(1+deltaU)` into the derived `xi1Slip` and `Collect`ing; its coefficients (L125-128) are likewise checked against paper targets. The former hand-typed `xi1Split` literal (diff old L73-78) is deleted.

**Assessment:**
Genuinely independent, not a rename. The load-bearing four coefficients (and the split coefficients) are now COMPUTED by Mathematica from the microscopic-log direct form via `Solve`/`Collect`/`Coefficient`, then compared against the paper's asserted coefficients. An algebra error in the SymPy literals could no longer be reproduced by the `.wl`, because the `.wl` does not consume those literals — it derives the value from `xi1Direct` and the Σ definitions. The comparison targets (`2*chi0/(1+chi0)` etc.) are the paper's claimed coefficients (the thing under test), which is the correct structure for a verification (compare a derived quantity against an asserted target). The `Coefficient[..., sigmaEta]` check (L84, target 0) independently confirms Xi_1 carries no Sigma_eta. Exec log L28-42 shows the gauge-elimination PASS, all six coefficient/constant PASS, and `Xi_1 direct - derived slippage form` PASS; split section L63-72 shows the four split-coefficient checks plus `Xi_1 split - slippage form` PASS. R_1 and Theta_1 retain their substitution-comparison structure, which is within scope (directive required independence for "at least the four-slippage-law and split coefficients"). No deviation from the consult resolution.

### F3 — tautological_check (physical-branch-drift self-subtractions)

**Classification:** resolved

**What changed:**
- SymPy L74-78: the six `def - (same def) == 0` self-subtractions are replaced by five notes §3 equality checks: `Sigma_chi = chi_1/chi_0`, `Sigma_eta = eta_1/eps_eta`, `Sigma_del = delta_U,1/delta_U`, `Sigma_eps = varepsilon_W/eps_W`, `Sigma_Z = zeta_Z - omega_W`, each subtracting an independently-defined microscopic quantity from `Sigma_*.subs(slip_subs)`.
- Mathematica L54-58: the analogous five `(sigma* /. slipSubs) - <independent def>` checks.

**Assessment:**
Non-tautological. E.g. `SigmaChi.subs(slip_subs)` expands to `gam1 + c1 - kU` while `chi1/chi0` reduces to the same via the independent definition `chi1 = chi0*(gam1+c1-kU)` (L48); the residual would be nonzero if `chi1`'s definition were wrong. These relate independently-defined quantities (the §3 boxed equalities), so they exercise deliverable (5) substantively. The directive's preferred replacement was applied verbatim. Both logs show all five PASS (SymPy L14-18; Mathematica L14-23).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Xi_1 direct - slippage form = 0` (L23) — substantive four-slippage law holds.
- `Xi_1 direct support-symbol leakage = set()` / `R_1 ...` / `Theta_1 ...` (L94-96) — F1 structural absence reported as empty set, not derivative-zero lines.
- `Sigma_chi = chi_1/chi_0 = 0` ... `Sigma_Z = zeta_Z - omega_W = 0` (L14-18) — F3 §3 equalities hold.

**Mathematica:** exit=0. Notable lines:
- `PASS: Xi_1 microscopic gauges eliminated` (L28) and `PASS: Xi_1 coeff sigma_Z/chi/eta/eps/del` + `constant term` (L29-40) — F2 independent coefficient derivation fires and passes.
- `PASS: Xi_1 direct - derived slippage form` (L42); split coefficients `PASS` (L63-72).
- `PASS: Xi_1 direct support-blind` / `R_1` / `Theta_1` (L79-81) — F1 `FreeQ` checks pass.
- `Stage 182 Mathematica audit passed.` (L93).

**Output freshness:** confirmed. sympy output mtime 2026-05-30 01:40:24 and mathematica output mtime 01:40:40 are both newer than their script mtimes (py 01:22:19, wl 01:26:17). Outputs re-generated post-fix.

## Material-change assessment

`material_change`: false.

All edits are verification-rigor improvements (replacing vacuous/tautological checks; making the second engine independent). No derived result changed: `Xi_1`, `R_1`, `Theta_1`, `Sigma_tr`, and the four-slippage law printed in the logs are identical to the pre-fix forms (the substantive identities passed before and after). No downstream unit depends on a numerical or symbolic value that moved. No `upstream_stale` propagation warranted on physics grounds.

## Side observations (non-blocking)

- The Mathematica `Print["free symbols of Xi_1: ", ...]` (wl L133) still emits garbled `List @@ xi1Slip` formatting (log L78). This was already noted by the original auditor as a non-load-bearing `Print`, not an assertion. It is not part of any finding and does not affect verification.
- F2's independence was strengthened only for the four-slippage law and the split (per directive scope); the R_1 and Theta_1 sections retain a substitution-and-compare structure that still echoes the SymPy algebra. This is within the directive's "at least" scope and not a defect, but a future pass could extend the independent-derivation route to R_1/Theta_1 if maximal independence is desired.

## Verdict justification

All three findings are resolved with the directive's prescribed edits applied verbatim and no collateral changes (the diff touches exactly the banner, the F3 drift checks, the F2 four-slippage and split blocks, and the F1 support-blindness block). F1's replacement is a genuinely non-vacuous structural free-symbol-absence assertion on the microscopic-log direct forms that would trip if a support drift were wired in. F2's Mathematica route now derives the load-bearing coefficients from the direct form via a consult-sanctioned gauge-Solve/Collect/Coefficient path (not a re-typed literal, not a Block), so it is independent of the SymPy literals and would not reproduce a SymPy algebra error. F3's replacement relates independently-defined quantities. Both engines exit 0 with every check passing, outputs are fresh, and no derived result changed downstream. Verdict: verified.
