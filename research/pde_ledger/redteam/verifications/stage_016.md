---
unit_id: 016
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 016

## Per-finding outcomes

### F1 — missing_verification_script

**Classification:** resolved

**What changed:**
Codex created the new file
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl`
(251 lines). The script declares dependent symbols `R[t,w,u,v]`, `mu[R,w]`, `Tw[R,w]`, `TO[R,w]`, `USig[R,w]` and derives every claim M1–M11 with Mathematica primitives (`D`, `Simplify`/`FullSimplify`, `Limit`, `Integrate`, `SphericalHarmonicY`, `Together`), per the directive's "no SymPy mirroring" requirement. Each `M*` block compares an independently computed residual to zero and `Exit[1]`s on failure, and (per the directive's mutation policy) each block also exercises at least one sign- or value-mutation that is expected to be nonzero.

**Assessment:**
The edit is correct and addresses the finding exactly as required. Substance checks per claim:

- **M1 (exact EL):** `elViaD` is built from the action `lagExact` via direct `D` chains in `t, w, u, v` (lines 51–55); `handExact` is the hand-written form (lines 56–60). The residual `elViaD - handExact` simplifies to 0 (output line 4); the mutation `-2*D[USig,R]` produces residual `2*Derivative[1,0][USig][...]`, which is nonzero (output line 6). Non-tautological — both sides are derived independently.
- **M2 (linear density):** `linDensity = D[lagExp, eps] /. eps -> 0` is taken from the full `lagExp` (built symbolically without presupposing the linear form). Mutation flips the sign on `UR0*eta` and gives residual `-2*eta*UR0` (nonzero). Non-tautological.
- **M3 (bulk after IBP):** `eBg = dTw0p - TwR0*R0p^2/2 - UR0` matches the directive M3 form; residual 0; mutation residual `-2*dTw0p*eta` (nonzero). Carried boundary term printed (`[-Tw(R0,w) R0' eta]`).
- **M4 (raw quadratic density):** `quadRaw = D[lagExp,{eps,2}]/2 /. eps->0`; residual 0; cross-term coefficient extracted by `D[D[quadRaw,eta],etaW]` independently and matched to `-TwR0*R0p`, residual 0; mutation residual `-2*R0p*TwR0` (nonzero). Non-tautological — the cross coefficient is extracted from quadRaw and compared, not assumed.
- **M5 (K_eta after IBP):** First proves the IBP substitution rule symbolically on independent test functions `aCheck[x], etaCheck[x]` (`[M5 product-rule]` residual 0), then applies the rule to `quadRaw` and matches the canonical form with `K_eta = URR0 - dTwR0p + TwRR0*R0p^2/2`. Residual 0. The mutation flips the middle sign (`+ dTwR0p` instead of `-`) and gives residual `dTwR0p*eta^2` (nonzero). Non-tautological — quadRaw came from `D[lagExp,{eps,2}]/2`, not from the canonical form. The printed K_eta line `R0p^2*TwRR0/2 + URR0 - d_TwR_R0p` matches the SymPy output transcript line 36 exactly (which expands to the directive's `K_eta = U_{Sigma,RR} - d/dw[T_{w,R} R0'] + 1/2 T_{w,RR} (R0')^2` at SymPy txt line 38).
- **M6 (concrete boundary discharges + IBP integrals):** Uses `Limit[..., w->+/-Infinity]` and `Integrate[..., {w,-Infinity,Infinity}]` on Gaussian/Lorentzian profiles directly. Sum-of-squares residual `linGaussian^2 + linLorentz^2 + quadGaussian^2 + quadLorentz^2` is 0 (all four limits vanish individually for the real-valued residuals to sum to zero), and the dedicated `linIBP`/`quadIBP` integrals are 0 individually. The `[M6 mutation]` `linGaussian - 1 = -1` confirms `linGaussian = 0` individually (not a vacuous sum). Non-tautological.
- **M7 (finite-endpoint discharge):** `Limit[-Bend*etaL,w->Infinity] - Limit[..,-Infinity] + 2 = 0` confirms the boundary discharge is exactly `-2`. Mutation drops the `+2` and yields residual `-2` (nonzero). Non-tautological.
- **M8 (atan probe + Lorentzian denominator):** `atanDischarge - Pi = 0` confirms π. The `Together` form `-1/(E^w^2*(1+w^2))` is printed and the `(1+w^2)*lorentzTogether + Exp[-w^2] = 0` check passes. Mutation residual `Pi` (nonzero) and denominator mutation residual `w^2/(E^w^2*(1+w^2))` (nonzero, since w^2 isn't identically 0). Non-tautological.
- **M9 (Y20 eigenvalue):** Build `y20` via `SphericalHarmonicY[2,0,th,ph]` then `ExpToTrig[FunctionExpand[...]]`, compute spherical Laplacian by direct `D` chain, residual `lapY20 + 6*y20 = 0`. Mutation residual `-1/8*(Sqrt[5/Pi]*(1 + 3*Cos[2*th]))` is nonzero (this is `-y20` after the +5 mutation). Non-tautological.
- **M10 (angular norm + stiffness):** Direct `Integrate` over `{ph, 0, 2 Pi}, {th, 0, Pi}` gives `angularNorm = 1` and `angularStiff = 6` (both printed). Sum-of-squares `(angularNorm-1)^2 + (angularStiff-6)^2 = 0` plus the mutation `angularStiff - 5 = 1` confirms each value individually. Non-tautological.
- **M11 (modal EL):** Build closed modal density with `(KMode + 6*TOMode) q^2`, compute EL via D chains, compare to hand form. Residual 0; mutation flips `6 -> 5` and gives residual `-q*TOmode` (nonzero). Non-tautological.

No collateral edits in the directive's `## Applied: F1` block beyond creating the named target file. No regressions visible. Output transcript shows STATUS: PASS at the final line.

## Exec log assessment

**SymPy:** exit=0. Notable lines from `redteam/exec_logs/stage_016_sympy.log`:
- Line 41: `K_eta(w)  = R0p**2*TwRR0/2 + URR0 - d_TwR_R0p` — matches Mathematica `[M5]` printed form exactly.
- Line 49: `mu_eta q_lm,tt - ∂_w(T_w q_lm,w) + [K_eta + l(l+1) T_Omega] q_lm = S_lm.` — for l=2 gives `+6 T_Omega`, matching Mathematica `[M10]` angular stiffness = 6.
- Line 60: `STATUS: PASS`.

**Mathematica:** exit=n/a. No `stage_016_mathematica.log` exists under `redteam/exec_logs/` (only `stage_016_sympy.log` and a 0-byte `stage_016_diff.patch`). However, the saved Mathematica output `mathematica/output/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.txt` exists, is fresh (mtime 1779397336 = post-creation of the `.wl` at 1779391209), and ends with `STATUS: PASS` at line 67 after all M1–M11 plus mutation guards printing PASS. Treating the saved output as the run transcript per the prompt's "saved outputs are already fresh — read them" instruction. Notable lines:
- Line 5 `[M1] PASS` … line 64 `[M11] PASS` (all eleven core claims pass).
- Line 28 `K_eta(w)  = R0p^2*TwRR0/2 + URR0 - d_TwR_R0p` — agrees with SymPy line 41 verbatim and with directive M5.
- Line 59 `Y20 angular stiffness = 6` — confirms the directive's `+6 T_Omega` factor.
- Line 67 `STATUS: PASS`.

**Output freshness:** Confirmed. The Mathematica `.wl` script mtime is 1779391209 (2026-05-21 13:20); the Mathematica output `.txt` mtime is 1779397336 (2026-05-21 15:02), well after the script was created. The SymPy `.txt` output mtime 1779397235 (2026-05-21 14:58) is also newer than the SymPy `.py` script. Both saved outputs are post-fix.

The 0-byte `stage_016_diff.patch` is consistent with Codex creating a single new file (an untracked `.wl`) and not modifying any pre-existing tracked file; `git diff` without `--cached` against tracked content reports no hunks.

## Material-change assessment

`material_change`: false.

The edit only adds a new Mathematica audit script that cross-confirms results the SymPy audit already produced (exact EL, `K_eta` form, Y20 angular factors, modal EL). No derived expression has changed; no downstream unit's symbolic inputs were touched. The orchestrator will still mark units > 016 as `upstream_stale: true` per its default rule, but no narrow re-audit is warranted — the new engine merely corroborates the existing SymPy outputs.

## Side observations (non-blocking)

- The Mathematica exec log `stage_016_mathematica.log` is absent from `redteam/exec_logs/` even though the saved output transcript is fresh and shows PASS. The orchestrator may want to ensure the exec-mathematica step writes a per-stage log alongside refreshing the `.txt` so future verifiers don't have to fall back to the saved output. Not a verification blocker.
- The `stage_016_diff.patch` is 0 bytes. For findings where the only edit is creating a new untracked file, the patch capture step might want to include `git diff --no-index /dev/null <new_file>` or just append the new file's path; otherwise verifiers can't see the change without reading the file directly. Again not a blocker — I read the `.wl` directly.
- `[M6]` and `[M10]` use sum-of-squares aggregation `(a-1)^2 + (b-6)^2`/etc to compress multiple zero-checks into a single `expectZero` call. This is fine for real-valued residuals (zero iff each summand is zero), and each block also has an individual-value mutation check that pins the components separately, so the aggregation is not a hidden tautology.

## Verdict justification

Codex created the required Mathematica audit script exactly at the directive's path, derived M1–M11 using Mathematica idioms (no line-for-line SymPy transliteration), and included substantive mutation guards on every claim. Each residual is computed from independently constructed expressions (e.g., `elViaD` from `D` chains vs. `handExact` from the gradient-form hand-written EL; `quadRaw` from `D[lagExp,{eps,2}]/2` vs. the canonical fluctuation density; `lapY20` from spherical-Laplacian operators vs. `-6*y20`), so the passing assertions are non-tautological. The Mathematica output's symbolic K_eta and `+6 T_Omega` factors match the SymPy output transcript at the directive's named lines (38, 44, 54). Verdict: verified.
