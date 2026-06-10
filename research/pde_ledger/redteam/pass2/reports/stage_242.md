---
unit_id: 242
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 242 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_242.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (row line 96; narrative 1283-1381; Thm `thm:app-part07-selected-twin-placement` line 1376-1379)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: *"Final front-end packet: actual selected placement coordinate, primitive ranking region, support/orbit packet split, weak-axisymmetric orbit packet, and separate outgoing finish-line datum N_Q-1."* The card + notes + appendix theorem (`thm:app-part07-selected-twin-placement`, line 1378) state: (1) the realized coordinate is `varrho_phys = 2(1-eps)/3`, `sigma_phys = 2eps/(1-eps) = 4/(3 varrho_phys)-2`, fed by `eps = epsW(1 - (2/11) deltaU/(1+deltaU))`; (2) thresholds rewrite to `eps_WLambda = 1/(2+beta^2)`, `eps_ULambda = beta/(1+beta+beta^2)`; (3) the selected branch lies **strictly inside the lowest symmetric twin window** `C_mix < Pi_tr < 2 C_mix` because the demand ratio `Pi_tr/C_mix = 4/3` (= `rho_alpha`, appendix line 1378); (4) reduced-state bridge `Zhat_W = Z_W Lambda_0/Lambda` leaves `R_target` invariant; (5) support-blindness of `eps, R_tr, R_target, eps_eta, (q_tr,q_nt,q_eta), (Theta_1,Xi_1,R_1)` w.r.t. `zeta`; (6) the finite orbit packet, the three infinitesimal compilers `dln eps/R_tr/R_target`, and the `Theta_1/Xi_1/R_1` direct-observable identities; (7) the two-packet split `R_target M_mix = C_mix` and `partial_zeta M_tr = M_mix(1-eps)/(1-zeta eps)^2`; (8) the separate finish line `N_Q = 1`. Items (1),(8) carry `\StatusOpen{}` for the nonlinear-branch realization; the algebraic identities are `\StatusExactClosure{}`.

## What the script claims to verify

Both engines assert every algebraic identity above. The SymPy `.py` uses `assert_zero(name, lhs-rhs)` (cancel/together then simplify), a Python-level strict-ratio comparison for the window (lines 94-97), an abstract-`Function(zeta)` substitution device for q-packet support-blindness (lines 135-163), and an `exp(t·d)`/`diff(log,t)|_{t=0}` parametrization for the infinitesimal compilers (lines 171-203). The Mathematica `.wl` uses `expectZero`/`expectTrue` (FullSimplify under real-domain `$Assumptions`), a `Resolve[ForAll[...],Reals]` quantifier-elimination **certificate** for the strict window (lines 154-165), direct `D[closedPacket, zeta]` on real closed forms for support-blindness (lines 185-212), and a `logDrift` total-log-differential `Σ x_i ∂_{x_i} log(expr)·dx_i` (lines 89-92) for the compilers. Neither asserts `N_Q=1` (status-open, out of algebraic scope; acceptable).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) varrho_phys, sigma_phys, eps | py 64-71 / wl 118-125 | match |
| (2) threshold rewrites | py 81-82 / wl 143-150 | match |
| (3) **strict** twin window `C_mix<Pi_tr<2C_mix` | py 94-97 (`1<ratio<2`) / wl 154-165 (`Resolve[ForAll]` strict `<`) | match |
| (4) reduced-state bridge | py 112 / wl 184 | match |
| (5) support-blindness eps/R_tr/R_target | py 113-115 (diff) / wl 185-187 (D) | match |
| (5) support-blindness q-packet | py 152-164 (abstract-fn) / wl 209-212 (direct D on closed forms) | match |
| (6) dln eps/R_tr/R_target | py 186,198,211 / wl 225,238,251 | match |
| (6) Theta_1/Xi_1/R_1 + compiler matrix | py 229-272 / wl 268-299 | match |
| (7) product law + support sensitivity | py 283-287 / wl 313-318 | match |
| (8) N_Q=1 finish line | none (status-open) | acceptable |

`paper_alignment: aligned`. The pass-1 partial (deliverable 3 reduced to a sign-free equality) is RESOLVED: both engines now test the window membership strictly.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 64,68 | varrho_phys/sigma_phys == 0 | claim 1 | yes |
| A2 | sympy | 81-82 | threshold rewrites == 0 | claim 2 | yes |
| A3 | sympy | 84-92 | gap-equalities (above/below bound) | claim 3 (positioning) | yes (gap = 1/(3r), 2/(3r) > 0) |
| A4 | sympy | 94-97 | `ratio==4/3` AND `1<ratio<2` (strict) | claim 3 (strict window) | yes |
| A5 | sympy | 112-115 | bridge + diff(.,zeta)==0 | claims 4,5 | yes |
| A6 | sympy | 152-164 | support-blind q-packet (abstract-fn) | claim 5 | yes |
| A7 | sympy | 186,198,211 | dln compilers == 0 | claim 6 | yes |
| A8 | sympy | 229-272 | Theta_1/Xi_1/R_1 + matrix/inv/det | claim 6 | yes |
| A9 | sympy | 283-287 | product law + support sensitivity | claim 7 | yes |
| M1 | math | 118-129 | M1/M2 placement/sigma == 0 | claim 1 | yes |
| M2 | math | 143-150 | threshold rewrites == 0 | claim 2 | yes |
| M3 | math | 152-165 | demand ratio + **`Resolve[ForAll] strict 1<4/3<2`** | claim 3 (strict window) | yes |
| M4 | math | 184-212 | bridge + `D[.,zeta]==0` on closed forms | claims 4,5 | yes |
| M5 | math | 216-251 | `logDrift` total-log-diff compilers == 0 | claim 6 | yes |
| M6 | math | 258-299 | Theta_1/Xi_1/R_1 + LinearSolve/Det | claim 6 | yes |
| M7 | math | 303-318 | product law + support sensitivity | claim 7 | yes |

No tautological rows. The pass-1 `Theta_1 := dln_Rtr` tautology (F3) is resolved: `Theta_1`/`thetaFromPacket` is now built from `-C_tr,* Sigma_tr` (py 214-218; wl 253-258) independently of `dln_Rtr`, then asserted equal to it (py 229; wl 268).

## Findings

None. (All four pass-1 findings — F1 transliteration, F2 insufficient strict-window, F3 Theta_1 tautology, F4 stale STAGE 225 banner — are confirmed RESOLVED below.)

## Independent-derivation check (Mathematica)

The `.wl` is now a genuine independent re-derivation, not a transliteration. The two pass-1-flagged transliterated devices are gone:

1. **Support-blind q-packet.** `.py` still uses the abstract-function substitution device (`Rtr_sb = sp.Function("Rtr_sb")(zeta)` … `sp.diff(q_tr_from_observables, zeta).subs(support_blind_subs)`, lines 135-163). `.wl` instead builds the real closed-form packet `closedObservablePacket` from `trailObservable/targetObservable/epsEta` and applies `D[closedObservablePacket, zeta]` directly (lines 202-212) — the genuine content (the closed forms contain no `zeta`). Different operation.

2. **Infinitesimal compilers.** `.py` uses a one-parameter `t`-flow `chi0*exp(t*dchi0)` then `diff(log(...),t).subs(t,0)` (lines 171-203). `.wl` uses `logDrift = MapThread[#1 D[Log[expr],#1] #2 &, {variables, driftSymbols}]` — a direct total logarithmic differential `Σ x_i (∂ log/∂ x_i) dx_i` (lines 89-92), invoked at lines 216,227,240. Different route to the same drift; no `Exp[t·d]` device.

3. **Strict window.** `.py` does a scalar rational comparison `nsimplify(Pi_tr/C_mix)` then Python `1<ratio<2` (lines 94-97). `.wl` does `Resolve[ForAll[{lambdaWin,epsilonWin}, Implies[...,1 < (4/3)·C/C < 2]], Reals]` (lines 154-165) — quantifier elimination over a 2-parameter family. Different primitive.

LOAD-BEARING strict-inclusion verdict: in BOTH engines the `4/3` prefactor is POSITED (carried forward from Stage 240, py:60 `Pi_tr=(4/3)C_mix`; wl numerator `(4/3)·C_mix`), and BOTH confirm `4/3 ∈ (1,2)` STRICTLY (py: `ratio>1 and ratio<2`; wl: `1 < ... < 2` inside `Resolve`). Neither is a tautology in the harmful sense: each would FAIL if the prefactor were set to 1 or 2 (py: the `assert ratio>1 and ratio<2` would raise; wl: `Resolve` would return `False`, tripping `expectTrue`→`fail`→`Exit[1]`). The `4/3` is not derived in THIS stage by either engine — that is correct, it is a Stage-240 forward; this stage's deliverable is the window membership, which both test strictly via genuinely different primitives (scalar compare vs quantifier-elimination certificate). Not a port.

## Engine cross-check

Both engines pass all checks and agree on the probe sample exactly: `epsilon = 17/55`, `rhoSelected/varrho_phys = 76/165`, `sigma = 17/19`, `R_tr/trailObservable = 19/25`, `R_target/targetObservable = 2749376/54071875`, `M_mix/mixedPacket = 89375/(1292 pi^2)`, `M_tr/throatPacket = 89375/(646 pi^2)`. No `engine_disagreement`. The agreement is now informative (the `.wl` reaches it by independent routes per the section above).

## Verdict justification

`verdict: clean`. I attacked the load-bearing strict inclusion specifically: the `4/3` is posited in both engines (legitimate Stage-240 forward), and the window membership is tested STRICTLY in both (`<`, not `≤`) via genuinely different primitives — a quantifier-elimination certificate in Mathematica vs a scalar rational comparison in SymPy; both would fail at the boundary. I confirmed the `.wl` de-transliterated the two pass-1-flagged devices (abstract-zeta → direct `D` on closed forms; `Exp[t·d]` → `logDrift` total differential), so the second engine is now independent. I confirmed the `Theta_1` tautology is gone (packet form `-C_tr,* Sigma_tr` is independent of `dln_Rtr`). I verified every `D[.,var]`/`logDrift` differentiates w.r.t. variables that actually appear, and that the `D[.,zeta]==0` checks are legitimate zero-derivative support-blindness claims (not absent-variable self-test traps — the content IS zeta-absence). The banner is canonical (`.wl:59` "STAGE 242"; output banner line 3 "STAGE 242"; no residual "225" in either script). Checkpoint bar MET: both engines present and independent, assertions substantive/non-tautological, paper alignment exact, the `4/3`/`2/3(1-eps)` load-bearing quantities re-derived. I read the paper card, notes, and appendix; the script's claims match.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `varrho_phys = 2(1-eps)/3` | py:61-66, wl:105-121, out epsilon-derived | notes:142-143, appendix:1293, card:13, appendix row:96 | MATCH |
| `sigma_phys = 2eps/(1-eps) = 4/(3rho)-2` | py:62,68, wl:113,127 | notes:154-158, appendix:1295 | MATCH |
| `eps = epsW(1-(2/11)deltaU/(1+deltaU))` | py:57, wl:96-99 | notes:118-120, appendix:1286-1288 | MATCH |
| `eps_WLambda = 1/(2+beta^2)` | py:81, wl:144-146 | notes:188 | MATCH |
| `eps_ULambda = beta/(1+beta+beta^2)` | py:82, wl:148-150 | notes:190 | MATCH |
| demand ratio `Pi_tr/C_mix = 4/3` (= rho_alpha) | py:94, wl:153, out PASS | notes:240,255, appendix:1378 | MATCH |
| strict window `C_mix < Pi_tr < 2 C_mix` | py:96, wl:154-165, out PASS | notes:246-248 | MATCH |
| `R_tr = (1+chi0/(1+deltaU))/(1+chi0)` | py:104, wl:170-173 | notes:276, appendix:1301-1303 | MATCH |
| `R_target = Lambda(1-eps_eta)(1-eps)^2/(Z_W(1+chi0)^2)` | py:105-107, wl:174-177 | notes:286, appendix:1306-1309 | MATCH |
| reduced bridge `Zhat_W = Z_W Lambda_0/Lambda` | py:102, wl:169 | notes:296-297 | MATCH |
| `R_target M_mix = C_mix = 8Lambda(1-eps)/pi^2` | py:283, wl:313 | notes:562-565 (+128 C_mix) | MATCH |
| `partial_zeta M_tr = M_mix(1-eps)/(1-zeta eps)^2` | py:286, wl:314-317 | notes:761-762 | MATCH |
| `dln eps`, `dln R_tr`, `dln R_target` compilers | py:184-210, wl:221-251 | notes:458-477 | MATCH |
| `Theta_1/Xi_1/R_1` direct-observable + matrix | py:218-272, wl:258-299 | notes:451-498, §5.1/6 | MATCH |
| finite q-packet `(q_tr,q_nt,q_eta)` | py:120-129 | notes:341-366, appendix:1313-1330 | MATCH |
| `S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps)`, `M_mix`, `M_tr` | py:277-281, wl:303-311 | notes:540-550, appendix:1357 | MATCH |
| `N_Q = 1` (finish line) | not script-asserted (status-open) | notes:602-606, appendix:1344, card:15 | MATCH (paper-side only) |

INTERNAL (probe-only / scaffolding, no finding): rational sample point `epsilon=17/55`, `rhoSelected=76/165`, `sigma=17/19`, `trailObservable=19/25`, `targetObservable=2749376/54071875`, `mixedPacket=89375/(1292 pi^2)`, `throatPacket=89375/(646 pi^2)` (notes:45-47 + wl:46-47 declare the sample probe-only, not on the proof path); `Lambda_0` symbol (carried, dimensionful bridge constant, drops out of the quotient); pass/fail flags; `t` flow variable (py); the gap quantities `1/(3rho)`, `2/(3rho)`.

reconciliation: complete; 17 deliverable values checked, 0 misaligned.

## Self-test notes

Variable independence: every `logDrift`/`D[.,var]` in the `.wl` differentiates w.r.t. variables that genuinely appear (`epsilonSelected` ↔ `{epsW,deltaU}`; `trailObservable` ↔ `{chi0,deltaU}`; `targetObservable` ↔ `{Lambda,ZW,epsEta,chi0,epsW,deltaU}`); the `D[.,zeta]==0` and `logDrift`-free `closedObservablePacket` zeta-checks are legitimate zero-derivative support-blindness claims (the closed forms contain no `zeta`), NOT the absent-variable self-test trap, since zeta-absence IS the deliverable. Strictness: confirmed the window `Resolve` uses strict `<` and would return `False` at prefactor 1 or 2 (the C_mix factor cancels, ratio = 4/3 ∈ (1,2)); the SymPy `assert ratio>1 and ratio<2` is the matching strict scalar test. No unbounded integrals (no parity concern). Trivial-case: the orbit-compiler `Det = c_eta ≠ 0` (epsEta ∈ (0,1)) so the matrix is invertible and the round-trips are non-vacuous.
