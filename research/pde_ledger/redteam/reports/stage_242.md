---
unit_id: 242
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: insufficient
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
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 96, 1283-1381)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: *"Final front-end packet: actual selected placement coordinate, primitive ranking region, support/orbit packet split, weak-axisymmetric orbit packet, and separate outgoing finish-line datum N_Q-1."* The card + notes + appendix theorem (Thm. `app-part07-selected-twin-placement`) state the stage proves: (1) the realized coordinate on the Stage-241 twin curve is `varrho_phys = 2(1-eps)/3` with `sigma_phys = 2 eps/(1-eps) = 4/(3 varrho_phys) - 2`, fed by `eps = epsW(1 - (2/11) deltaU/(1+deltaU))`; (2) the two Stage-241 thresholds rewrite to `eps_WLambda = 1/(2+beta^2)`, `eps_ULambda = beta/(1+beta+beta^2)`; (3) the selected branch sits **strictly inside the lowest symmetric twin window** `C_mix < Pi_tr < 2 C_mix` (because the demand ratio `Pi_tr/C_mix = 4/3`, equivalently `rho_alpha = 4/3`); (4) the reduced-state bridge `Zhat_W = Z_W Lambda_0/Lambda` leaving `R_target` invariant; (5) **support-blindness**: `eps, R_tr, R_target, eps_eta, (q_tr,q_nt,q_eta), (Theta_1,Xi_1,R_1)` are all independent of the support ratio `zeta`; (6) the finite orbit packet `(q_tr,q_nt,q_eta)` and infinitesimal compilers `dln eps, dln R_tr, dln R_target` plus the direct-observable `Xi_1` identity; (7) the two-packet split with `R_target M_mix = 8 Lambda(1-eps)/pi^2 = C_mix` and `partial_zeta M_tr = M_mix (1-eps)/(1-zeta eps)^2`; (8) the separate finish line `N_Q = 1`. Items (1),(2),(8) carry `\StatusOpen{}` flags for the nonlinear-branch realization and `N_Q-1`; the algebraic identities are `\StatusExactClosure{}`.

## What the script claims to verify

The SymPy script asserts each algebraic identity above as `assert_zero(name, lhs - rhs)`: the placement coordinate/sigma identities, the two threshold rewrites, a pair of "above/below bound" decompositions, the `R_target` bridge, three direct `diff(...,zeta)==0` support-blindness checks, a support-blind-propagation construction for the q-packet via abstract zeta-functions, the three infinitesimal compilers, the `Theta_1/Xi_1/R_1` direct-observable identities, a 3x3 orbit-compiler matrix with its determinant/inverse/zero-vector round-trip, and the two-packet product law + support sensitivity. A rational sample point is printed (probe-only, not asserted). The Mathematica script mirrors these one-for-one. Neither script asserts `N_Q = 1` (item 8) — that is a status-open finish-line datum and is reasonably out of algebraic scope.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) varrho_phys = 2(1-eps)/3, sigma_phys | `varrho_phys` / `sigma_phys` asserts (py 64-71) | match |
| (2) threshold rewrites eps_WLambda, eps_ULambda | py 81-82 | match |
| (3) **strict** twin-window inclusion `C_mix < Pi_tr < 2 C_mix` | py 85-92 "above/below bound" | **partial** — only an algebraic gap-equality, no strictness/sign asserted |
| (4) reduced-state bridge R_target invariant | py 107 | match |
| (5) support-blindness of eps, R_tr, R_target | py 108-110 (direct diff) | match |
| (5) support-blindness of q-packet | py 147-158 (abstract-fn construction) | partial — conditional, not direct |
| (5) support-blindness of (Theta_1,Xi_1,R_1) | (not separately tested; inherited from dln's having no zeta) | partial |
| (6) dln eps, dln R_tr, dln R_target | py 181,193,206 | match |
| (6) Theta_1, Xi_1, R_1 direct-observable | py 219-221 + matrix 232-262 | match except Theta_1 (tautological) |
| (7) R_target M_mix = C_mix; partial_zeta M_tr | py 273-277 | match |
| (8) N_Q = 1 finish line | none | missing (status-open; acceptable) |

`paper_alignment: partial` — the algebraic identities are faithfully ported, but the strict-window-inclusion deliverable (3), which the appendix theorem states as the support-side meaning of the selected branch, is reduced to a sign-free equality.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 64 | `varrho_phys - 2/3(1-eps) == 0` | claim 1 | yes |
| A2 | sympy | 68 | `sigma_phys - 2eps/(1-eps) == 0` | claim 1 | yes |
| A3 | sympy | 81-82 | threshold rewrites == 0 | claim 2 | yes |
| A4 | sympy | 85 | `sigma_sel - (1/r-2) - 1/(3r) == 0` | claim 3 | **no (identity, no sign)** |
| A5 | sympy | 89 | `(2/r-2) - sigma_sel - 2/(3r) == 0` | claim 3 | **no (identity, no sign)** |
| A6 | sympy | 107 | R_target bridge == 0 | claim 4 | yes |
| A7 | sympy | 108-110 | `diff(eps/R_tr/R_target, zeta) == 0` | claim 5 | yes |
| A8 | sympy | 147-158 | support-blind q-packet propagation | claim 5 | partial (conditional) |
| A9 | sympy | 181,193,206 | infinitesimal compilers == 0 | claim 6 | yes |
| A10 | sympy | 219 | `Theta_1 - dln_Rtr == 0` | claim 6 | **no (tautology: Theta_1 := dln_Rtr)** |
| A11 | sympy | 220-221 | Xi_1, R_1 identities == 0 | claim 6 | yes |
| A12 | sympy | 232-262 | orbit-compiler matrix/inv/det | claim 6 | yes |
| A13 | sympy | 273-277 | product law + support sensitivity | claim 7 | yes |
| M1-M13 | mathematica | 91-298 | one-for-one mirror of A1-A13 | same | same flaws |

## Findings

### F1 — mathematica_transliteration

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl:83-298`

**What's wrong:**
The `.wl` is a line-by-line transliteration of the `.py`, not an independent re-derivation. It reuses the SymPy script's exact decomposition, variable names, intermediate quantities, choreography, and even idiosyncratic test constructions. Three corresponding sections:

1. The "above/below bound" decomposition is ported verbatim, term-for-term:
   - py:87 `sigma_sel - (1 / varrho - 2) - 1 / (3 * varrho)`
   - wl:113 `sigmaSel - (1/varrho - 2) - 1/(3 varrho)`
   An independent engine would either test the window membership directly (e.g., `Reduce`/`Resolve` that `1 < Pi_tr/C_mix < 2`) or pick a different algebraic factoring — not reproduce the same artificial split.

2. The abstract-function support-blind trick is reproduced one-for-one:
   - py:130-149 `Rtr_sb = sp.Function("Rtr_sb")(zeta)` … `sp.diff(q_tr_from_observables, zeta).subs(support_blind_subs)`
   - wl:146-165 `RtrSb = RtrSbFn[zeta]` … `D[qtrFromObservables, zeta] /. supportBlindRules`
   Same expression `-Cstar Log[RtrSb/Rtr]`, same substitution-of-derivatives device. This is a literal port of a SymPy idiom into Mathematica syntax.

3. The infinitesimal compiler uses the identical `t`-exponential parametrization and the identical orbit-compiler matrix and round-trip vector:
   - py:166-262 `chi0_t = chi0 * sp.exp(t * dchi0)` … `recovered_direct_packet - sp.Matrix([Theta_var, R_var, -(1 - epsEta)*(Xi_var + R_var)/epsEta])`
   - wl:181-276 `chi0T = chi0 Exp[t dchi0]` … `recoveredDirectPacket - {ThetaVar, RVar, -((1 - epsEta)(XiVar + RVar))/epsEta}`

This violates the second-engine policy: both engines echo the same algebra rather than deriving the result by independent routes, so the Mathematica run provides no genuine cross-check.

**Why this matters:**
A transliterated second engine cannot catch a decomposition error in the first — any mistake in the SymPy choreography is faithfully copied into Mathematica, and both "agree" on the same wrong path. The checkpoint's dual-engine guarantee is hollow.

**Required change:**
Re-author the `.wl` to derive the Stage-242 claims by independent Mathematica-native routes (different decomposition, native primitives). See directive F1 claim manifest.

**Verification:**
New `.wl` keeps the `_mathematica_audit.wl` suffix, exits 0, and visibly uses a different derivation route (e.g., `Reduce`/`Resolve` for the window inequalities, closed-form `D[...,zeta]==0` on the actual closed forms rather than the abstract-function substitution device, direct substitution into closed forms rather than the `Exp[t d]` first-order trick where avoidable).

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py:84-92`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl:110-118`

**What's wrong:**
The paper's deliverable (3) is a **strict inclusion**: the appendix theorem (`stage_appendix_part07.tex:1378`) states the branch "lies on the symmetric lowest-twin support slice" with `rho_alpha = 4/3`, and the notes (`...stage242...md:243-249`) box `C_mix < Pi_tr < 2 C_mix`. The two script checks named "selected branch lies above mixed-only bound" / "below non-twin bound" are pure algebraic **equalities** that carry no sign content:
- py:85-88 asserts `sigma_sel - (1/varrho - 2) - 1/(3 varrho) == 0`, i.e. `sigma_sel = (mixed-only bound) + 1/(3 varrho)` — an identity true by construction, not a proof that the gap `1/(3 varrho)` is positive.
- py:89-92 asserts `(2/varrho - 2) - sigma_sel - 2/(3 varrho) == 0` — likewise.

The assertion names promise an inclusion (`lies above`, `lies below`); the assertions deliver only that the two gaps equal `1/(3 varrho)` and `2/(3 varrho)`. The strictness (positivity of those gaps, i.e. the actual window-membership claim, equivalently `1 < Pi_tr/C_mix = 4/3 < 2`) is never asserted. Under `varrho > 0` these gaps are provably positive, but the script never makes that step, so it does not exercise the deliverable.

**Why this matters:**
On a checkpoint, the window-inclusion is a load-bearing physics deliverable (it is the support-side meaning of the selected branch per Thm. `app-part07-selected-twin-placement`). As written, both engines would still pass even if the branch were *outside* the window, because the equality identity is independent of the sign of the gaps.

**Required change:**
Add an explicit strict-positivity / window-membership check that exercises the inclusion, not just the equality. See directive F2.

**Verification:**
A new assertion confirms the demand ratio lies strictly in `(1,2)` and/or the two gaps are strictly positive over the declared domain (`varrho > 0`). The check must fail if the ratio were set to e.g. 1 or 2.

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py:208,219`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl:235,242`

**What's wrong:**
`Theta_1` is *defined* as `dln_Rtr` (py:208 `Theta_1 = dln_Rtr`; wl:235 `Theta1 = dlnRtr`) and then "verified" by asserting `Theta_1 - dln_Rtr == 0` (py:219; wl:242). This is `x - x == 0` — it cannot fail regardless of the physics. The paper's claim here is the substantive `Theta_1 = d ln R_tr` direct-observable identity (notes §5.2, `app-part07-stage242-orbit-lock-inf`), but since `Theta_1` is assigned to be `dln_Rtr` with no intervening derivation, the assertion exercises nothing. (By contrast `Xi_1` and `R_1` involve a genuine cancellation, so they are not tautological.)

**Why this matters:**
The check appears in the inventory as if it confirms the `Theta_1 = dln R_tr` deliverable, but it gives no evidence. It is dead weight that inflates the pass count.

**Required change:**
Make the `Theta_1` check substantive: derive `Theta_1` from its defining packet form (`Theta_1 = -C_tr,* Sigma_tr` with `Sigma_tr = (1+deltaU)dlnchi0 + (1+chi0)dlndeltaU` and `C_tr,*` per notes §5.1) and assert that this packet form equals `dln_Rtr`. That is the actual identity the paper boxes (`Theta_1 = d ln R_tr`). Then `Theta_1` is no longer defined to be `dln_Rtr` by fiat.

**Verification:**
The new `Theta_1` assertion's LHS is built from `Sigma_tr` and `C_tr,*` (independent of `dln_Rtr`), and the residual against `dln_Rtr` simplifies to 0. The check would fail if `C_tr,*` were mis-stated.

### F4 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl:59`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.txt:11`

**What's wrong:**
The Mathematica banner reads `banner["STAGE 225 — ACTUAL TWIN-SUPPORT PLACEMENT AND COHERENT ORBIT-LOCK COMPILER"]` (wl:59), and this propagates to the saved output (txt:11 `STAGE 225 — ...`). Every canonical anchor (filename, paper card `stage_242.tex`, MANIFEST, the in-script "All Stage 242 ... passed." footer at wl:321) says 242. The "225" is a stale label, a known project-wide incomplete-renumber artifact from the EM-extension realignment. Canonical numbering is ground truth; this is the specific stale label, flagged — NOT a request for a batch renumber sweep.

The `.txt` mtimes are newer than the scripts, so the output content itself is fresh; this finding is the label only.

**Why this matters:**
The banner mislabels the stage in the saved transcript, which can mislead a downstream reader cross-referencing audit output to the wrong stage number.

**Required change:**
Change `wl:59` banner string from `STAGE 225` to `STAGE 242`. (This will be naturally subsumed if F1 re-authors the `.wl`, but state it explicitly so it is not lost.) Re-run regenerates the `.txt` with the corrected banner.

**Verification:**
`.wl:59` reads `STAGE 242 — ...`; regenerated `.txt:11` shows `STAGE 242`.

## Independent-derivation check (Mathematica)

The `.wl` is **not** an independent derivation. It is a transliteration (see F1) — the section structure (I-VII), variable names, intermediate quantities, the abstract-zeta-function support-blind device, the `Exp[t d]` infinitesimal parametrization, the orbit-compiler matrix, and the rational sample point all mirror the `.py` one-for-one. The only Mathematica-native idioms are cosmetic (banner/formatting helpers, `FullSimplify` vs `simplify`, `D[...,zeta]` vs `sp.diff`). No genuine second route is present.

## Engine cross-check

Both engines pass all checks and agree on the sample point exactly: `epsilon = 17/55`, `varrho_phys = 76/165`, `sigma_phys = 17/19`, `R_tr = 19/25`, `R_target = 2749376/54071875`, `M_mix = 89375/(1292 pi^2)`, `M_tr = 89375/(646 pi^2)`. No `engine_disagreement`. (The agreement is unsurprising and uninformative given the transliteration, F1.)

## Verdict justification

`verdict: findings`. The SymPy algebra is faithful to the paper on the placement coordinate, threshold rewrites, reduced-state bridge, direct support-blindness, infinitesimal compilers, the Xi_1/R_1 identities, the orbit-compiler matrix, and the product law — I attacked the dln-eps coefficient `2/((1+δ)(11+9δ))`, the `R_target M_mix = C_mix` cancellation, and the orbit-compiler row-3/det identities, and they all hold. The findings are: (F1) the second engine is a transliteration, defeating the dual-engine guarantee on a checkpoint; (F2) the strict twin-window inclusion — a load-bearing checkpoint deliverable per the appendix theorem — is reduced to a sign-free equality and never actually exercised; (F3) the `Theta_1 = dln R_tr` check is tautological by construction; (F4) a stale `STAGE 225` banner. No `paper_misalignment` (the script's identities match the paper; the gaps are coverage gaps, not contradictions), so `needs_user_resolution: false`. No `stop_cold`: none of the fixes change a derived constant a downstream unit depends on — they add/strengthen checks and re-author the second engine.

## Self-test notes

Variable independence: confirmed the `diff(...,zeta)` checks (py:108-110) act on closed forms that genuinely contain no `zeta`, so they return identically 0 — these correctly verify support-blindness (not a false-zero trap, since the claim *is* zero-derivative). For the F2 fix I verified `Pi_tr/C_mix = 4/3` is a pure number (C_mix cancels), so `1 < 4/3 < 2` is a concrete strict inequality, not a derivative; the trivial-case check is direct (4/3 lands strictly between 1 and 2). For F3 I confirmed `Theta_1` is presently assigned `= dln_Rtr` (so the residual is structurally 0) and that the proposed packet-form LHS `-C_tr,* Sigma_tr` is independent of `dln_Rtr`, making the new check non-vacuous. Parity/symmetry: no unbounded integrals in this stage. Paper round-trip: the F2/F3 fixes introduce no new constants beyond `C_tr,*`, `Sigma_tr` which the notes §5.1 already state, so no new paper_misalignment.
