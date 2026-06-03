---
unit_id: 237
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 237 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_237.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows/§ at lines 86, 1010–1065)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 237 ("actual-branch dressing compiler") computes the surviving rigid-mouth dressing coordinate `q_eta` directly on the actual coherent branch. The card's `\stagefield{Output}` states: "Actual-branch dressing compiler and support-blindness theorem: explicit support-sector shifts do not by themselves fix the dressing coordinate." The notes and Part VII appendix (eqs. app-part07-rm-finite-packet through app-part07-qeta-support-blind) enumerate the distinct deliverables: (1) the finite rigid-mouth packet `q_nt = ln((1-eps_eta)/(1-eps_eta_ref)) - ln(R_target/R_target_ref)`, `q_eta = ln(eps_eta/eps_eta_ref)` with the `q_tr=0` rigid-mouth reduction; (2) the exact finite static-blind curve `R_target/R_target_ref = (1-eps_eta)/(1-eps_eta_ref)` and its inverse `q_eta(R_target)`; (3) the first-order packet compiler `q_nt = -R1 - c_eta E1`, `q_eta = E1` with tangent slope `-c_eta`, `c_eta = eps_eta*/(1-eps_eta*)`; (4) the microscopic compiler `q_eta = 2 ln(c_etaU/ref) - ln(K_U/ref) - ln(K_eta_eff/ref)` and its first-order drift extractor `q_eta = 2 c1 - kappa_U - kappa_eta`; (5) the support-blindness theorem `∂_zeta q_eta = 0`, `∂_{M_tr} q_eta = 0`, `∂_{lambda_phi} q_eta = 0`, `∂_{K_phi_eff} q_eta = 0`; (6) the post-static equal-drift dependent-plane ray `y_eta = -q_eta (0,1,1)^T` and the codimension-two orbit-lock theorem `q_nt=0 & q_eta=0 <=> R_target=R_target_ref & eps_eta=eps_eta_ref`.

## What the script claims to verify

The script docstring/print banner says it audits "actual-branch dressing compiler, finite static-blind curve, and support-blind post-static orbit-lock theorem." It exercises: the rigid-mouth reduction of the finite packet (§1, lines 89–116); the static-blind curve and its `q_eta<->R` inverse round-trip (§2, lines 121–146); the first-order series compiler, tangent slope, packet determinant and matrix product (§3, lines 151–182); the microscopic log-split identity and its first-order drift extractor (§4, lines 187–222); the four support-blindness derivative identities (§5, lines 227–249); and the dependent-plane ray plus codimension-two orbit-lock checks (§6, lines 254–275). The §1–§4 and §6 ray checks faithfully match the paper formulas. The §5 support-blindness block does not actually exercise the support-blindness physics (see F2).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) Finite packet + rigid-mouth reduction | lines 101–109 | match |
| (2) Static-blind curve + inverse `q_eta(R)` | lines 125, 132–141 | match |
| (3) First-order compiler, tangent slope `-c_eta` | lines 164–168, 171 | match (line 174 matrix-product check is tautological — F3) |
| (4) Microscopic log compiler + drift extractor `2c1-kU-kEta` | lines 202, 217 | match |
| (5) Support-blindness `∂_zeta=∂_{M_tr}=∂_{lambda_phi}=∂_{K_phi_eff} q_eta = 0` | lines 238–241 | **partial/vacuous** — derivatives are identically zero by symbol non-dependence, not by demonstrated invariance (F2) |
| (6) Dependent-plane ray + codim-2 orbit-lock | lines 258–261, 264 | match (line 265 `qeta_param.subs(qeta_param,0)` is tautological — F4) |

`paper_alignment: aligned` — every paper formula has a script counterpart and no constant/identity disagreement was found; the defect is in *how* deliverable (5) is verified (vacuous), not in *what* is claimed. So this is a script-quality issue (`insufficient_verification`), not a `paper_misalignment`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 101 | `assert_zero(q_tr_rm)` | (1) rigid-mouth `q_tr=0` | yes |
| A2 | sympy | 102–108 | `assert_zero(q_nt_rm - reduced)` | (1) | yes |
| A3 | sympy | 109 | `assert_zero(q_eta_rm - log ratio)` | (1) | yes |
| A4 | sympy | 125 | `assert_zero(q_nt_ratio.subs(static_blind))` | (2) curve | yes |
| A5 | sympy | 132–139 | inverse round-trip both directions | (2) inverse | yes |
| A6 | sympy | 140–141 | `Rratio(0)=1`, `eps_param(0)=eps_ref` | (2) endpoint | yes (mild) |
| A7 | sympy | 164 | `assert_zero(q_eta_linear - t E1)` | (3) | yes |
| A8 | sympy | 165 | `assert_zero(q_nt_linear + t(R1+c_eta E1))` | (3) | yes |
| A9 | sympy | 168 | `assert_zero(tangent_slope + c_eta)` | (3) tangent | yes |
| A10 | sympy | 171 | `assert_zero(M_packet.det()+1)` | (3) | yes |
| A11 | sympy | 174 | `assert_matrix_zero(M_packet*[R1,E1] - [-R1-c_eta E1, E1])` | (3) | **no — tautological** (F3) |
| A12 | sympy | 202 | `assert_zero(q_eta_micro - q_eta_micro_split)` | (4) | yes |
| A13 | sympy | 217 | `assert_zero(q_eta_micro_linear - t(2c1-kU-kEta))` | (4) | yes |
| A14 | sympy | 238–241 | `assert_zero(diff(q_eta_micro, {zeta,M_tr,lambda_phi,K_phi_eff}))` | (5) support-blindness | **no — vacuous, var-independence** (F2) |
| A15 | sympy | 258 | `assert_matrix_zero(y_eta - y_eta_from_eps)` | (6) ray | yes |
| A16 | sympy | 259–261 | `assert_matrix_zero(y_eta - y_eta_from_R(...))` | (6) ray via R | yes |
| A17 | sympy | 264 | `assert_zero(qeta_from_Rratio.subs(Rratio,1))` | (6) orbit-lock endpoint | yes |
| A18 | sympy | 265 | `assert_zero(qeta_param.subs(qeta_param,0))` | (6) | **no — pure tautology** (F4) |
| A19 | sympy | 268–269 | `assert_matrix_zero(M_packet.LUsolve([0,0]))` | (6) codim-2 | weak (homogeneous solve always 0) |

## Findings

### F1 — missing_verification_script (subtype missing_mathematica)

**Severity:** medium
**Files:**
- `(missing)` — no `.wl` exists for unit 237

**What's wrong:**
This is a non-checkpoint, non-status-only unit (`is_status_only_candidate: False`, checkpoint: False). The stage card explicitly notes "Mathematica audit: none yet." The entire stage is symbolic: logarithmic identities, two-term Taylor series in `t`, a 2x2 determinant, and partial derivatives. All of this is squarely within Mathematica's native primitives (`FullSimplify`, `Series`, `Det`, `D`). Neighboring same-batch Part VII rigid-mouth-packet stages already carry Mathematica audits (`mathematica/moving_throat_pde_stage239_..._mathematica_audit.wl`, `..._stage242_..._mathematica_audit.wl`), demonstrating that this exact class of algebra is independently verifiable in Mathematica. The dual-engine rule's test is "is it possible," and it clearly is.

**Why this matters:**
Single-engine verification leaves no cross-check; the §5 vacuity (F2) in particular would have surfaced immediately under an independently-structured second engine.

**Required change:**
Codex writes a NEW independent-route Mathematica script (claim manifest in directive). Independent route = native Mathematica decomposition, NOT a transliteration of the `.py`.

**Verification:**
`mathematica/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_mathematica_audit.wl` exists, contains the M1–M7 checks below, and exits 0.

### F2 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py:227-241`

**What's wrong:**
The support-blindness theorem is the headline deliverable of the stage (card `\stagefield{Output}`; notes §5). The script verifies it by:
```python
zeta, M_tr = sp.symbols("zeta M_tr", positive=True, real=True)          # line 227
lambda_phi, K_phi_eff = sp.symbols("lambda_phi K_phi_eff", ...)         # line 228-229
...
zeta_expr = lambda_phi**2 * K_W_eff / (lambda_W**2 * K_phi_eff)         # line 235  (never used in an assert)
M_tr_expr = M_mix * (1 + zeta_expr * (1 - eps) / (1 - zeta_expr * eps)) # line 236  (never used in an assert)
assert_zero(sp.diff(q_eta_micro, zeta))                                 # line 238
assert_zero(sp.diff(q_eta_micro, M_tr))                                 # line 239
assert_zero(sp.diff(q_eta_micro, lambda_phi))                           # line 240
assert_zero(sp.diff(q_eta_micro, K_phi_eff))                            # line 241
```
But `q_eta_micro` (lines 193–196) is built only from `{c_etaU, K_U, K_eta_eff}` and their `_ref` versions. The symbols `zeta`, `M_tr`, `lambda_phi`, `K_phi_eff` are fresh, independent symbols that are NEVER substituted into `q_eta_micro` (verified: the only `.subs` on `q_eta_micro` is the drift expansion at lines 206–212, which uses `c_etaU`, `K_U`, `K_eta_eff` only). Therefore each `sp.diff(q_eta_micro, VAR)` is identically zero *because `VAR` is not a free symbol of `q_eta_micro`*, not because the physics is support-blind. This is the documented "variable-independence" trap: `diff(EXPR, VAR)` with `VAR ∉ free_symbols(EXPR)` is `0` for ANY `EXPR`. The check cannot fail no matter what `q_eta_micro` is, so it proves nothing about support-blindness. The notes claim "the exact dressing coordinate satisfies `∂_zeta q_eta = 0`" — but the script never expresses `q_eta` *as a function that could have* contained `zeta`. (`zeta_expr` and `M_tr_expr` are computed and then only printed, lines 244–245; they are dead w.r.t. the asserts.)

**Why this matters:**
The single most physically interesting claim of the stage — that support-sector enhancement cannot remove the post-static dressing obstruction — is asserted but not demonstrated. A reader trusting the green "dq_eta/dzeta = 0" output (lines 51–54) would believe the theorem was checked, when in fact the checks are vacuous identities true of any expression independent of those symbols.

**Required change:**
Make `q_eta_micro` an explicit function of the support variables before differentiating, so the vanishing derivative is a real consequence of the dressing/support sector split rather than of symbol absence. Concretely, the microscopic compiler's inputs `c_etaU`, `K_U`, `K_eta_eff` are the dressing-sector quantities; the support sector enters through the *separate* variables (`lambda_phi`, `K_phi_eff`, `lambda_W`, `K_W_eff`, `eps`, `M_mix`) that build `zeta` and `M_tr`. The genuine theorem is: when each dressing-sector quantity is, at most, a function of dressing variables and is held independent of the support variables, `q_eta` has zero derivative w.r.t. every support variable AND w.r.t. the composite `zeta`, `M_tr`. The fix must construct a representation in which the support variables COULD appear in `q_eta_micro` and then show they do not. See directive F2 for the exact construction Codex must build and self-tested below. (Codex designs the route; the directive states only the requirement + acceptance criteria.)

**Verification:**
After the fix, the §5 block must (a) build `q_eta` (or `q_eta_micro`) in a form whose `free_symbols` plausibly include the support variables under the test, and (b) `assert_zero` the four derivatives such that swapping the dressing-sector definition for one that DID depend on a support variable would make at least one assert FAIL. The verifier confirms the four `assert_zero` calls remain and that the script exits 0; reviewer confirms the construction is non-vacuous (the differentiated expression genuinely contains, or could contain, the support variables).

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py:170-174`

**What's wrong:**
```python
M_packet = sp.Matrix([[-1, -c_eta], [0, 1]])              # line 170
...
q_lin_vec = sp.simplify(M_packet * sp.Matrix([R1, E1]))   # line 173
assert_matrix_zero(q_lin_vec - sp.Matrix([-R1 - c_eta * E1, E1]))  # line 174
```
`q_lin_vec` is *defined* as `M_packet * [R1, E1]`, which by direct multiplication of the hand-written matrix is exactly `[-R1 - c_eta*E1, E1]`. The assert then compares that product to its own hand-expansion. It is algebraically guaranteed by construction and cannot fail regardless of any physics. (The substantive packet check is A8 at line 165, which is fine; line 174 is redundant scaffolding.)

**Why this matters:**
A green pass here gives false confidence that the packet matrix encodes the right map. The real anchoring is line 165; line 174 adds nothing and risks masking a wrong `M_packet` if line 165 were ever removed.

**Required change:**
Tie line 174 to an independent target rather than to the product of the same matrix. Replace the RHS hand-expansion with the linear forms derived in §3 from the series (`q_nt_linear`, `q_eta_linear`) so the assert checks that `M_packet*[R1,E1]` reproduces the *series-derived* `(q_nt, q_eta)` linearization, not its own definition. Specifically, assert `M_packet * Matrix([R1, E1]) - Matrix([q_nt_linear/t (with t->1 scaling consistent), q_eta_linear/...])` — i.e. compare against the independently series-expanded packet from lines 156–162, not against a literal restatement of the matrix product. If Codex finds that awkward, the minimal acceptable fix is to delete line 173–174's redundant restatement and instead assert that `M_packet * Matrix([R1, E1])` equals the column vector `Matrix([q_nt_linear, q_eta_linear])` with the `t` factor handled consistently (both `q_*_linear` carry a single power of `t`; multiply the packet RHS by `t` or strip `t` from the series).

**Verification:**
Line 174 (or its replacement) compares the packet action to a quantity NOT defined as that same product. Swapping a wrong entry into `M_packet` would now make the assert fail. Script exits 0.

### F4 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py:265`

**What's wrong:**
```python
assert_zero(sp.simplify(qeta_param.subs(qeta_param, 0)))  # line 265
```
This substitutes `0` for the symbol `qeta_param` and asserts the result is `0`. It is `assert 0 == 0`. It exercises no physics and cannot fail. The accompanying line 264 (`qeta_from_Rratio.subs(Rratio, 1) == 0`) IS a real endpoint check; line 265 is pure tautology masquerading as the "q_eta = 0 at the orbit-lock point" check.

**Why this matters:**
It pads the codimension-two orbit-lock verification with a vacuous assert, giving false coverage to the `q_eta=0` half of the orbit-lock point. The intended claim — "at `eps_eta = eps_eta_ref`, `q_eta = 0`" — should be checked by substituting `eps_eta = eps_eta_ref` into the actual `q_eta` definition, not by zeroing the parameter.

**Required change:**
Replace line 265 with a substantive endpoint check of the dressing coordinate, e.g. assert that `q_eta` (line 95, `sp.log(eps_eta/eps_eta_ref)`) evaluates to 0 when `eps_eta = eps_eta_ref`: `assert_zero(q_eta.subs(eps_eta, eps_eta_ref))`. This makes the `q_eta=0` half of the orbit-lock point a real check on the actual `q_eta` expression.

**Verification:**
Line 265 substitutes into a `q_eta`-bearing expression (not into the bare parameter) and the assert would fail if `q_eta` were mis-defined. Script exits 0.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot be assessed. F1 requires a NEW independent-route script; the directive's claim manifest forbids a line-by-line port.

## Engine cross-check

Only the SymPy engine is present; `engines_agree: n/a`.

## Verdict justification

`verdict: findings`. The §1–§4 algebra and the §6 ray identities hold up under attack: the rigid-mouth reduction, static-blind curve and its inverse round-trip, first-order series compiler, tangent slope `-c_eta`, packet determinant, microscopic log-split, and drift extractor `2c1-kU-kEta` are all non-tautological and faithfully match the paper card and the Part VII appendix formulas (eqs. app-part07-rm-finite-packet…app-part07-qeta-first-order). I tried to break the static-blind inverse (it round-trips both directions) and the drift extractor (the series genuinely produces `2c1-kU-kEta`); both held. The defects are: (F1) no Mathematica second engine though the algebra is plainly Mathematica-verifiable and sibling stages 239/242 already have one; (F2, high) the support-blindness theorem — the stage's headline claim — is verified by differentiating an expression w.r.t. symbols it does not contain, so the four `∂q_eta = 0` asserts are vacuous; (F3, F4, low) two tautological scaffolding asserts (the packet matrix-product restatement and `qeta_param.subs(qeta_param,0)`). None is a `paper_misalignment` — every paper formula has a matching script formula and no constant disagrees; the physics claims are correct, only their verification is (in F2) vacuous. No `UNFIXABLE` (the math is sound and reconcilable) and no `CRITICAL_DOWNSTREAM` (fixes strengthen checks without changing any derived constant or sign that stages 238/239/242 consume).

## Self-test notes

I ran the variable-independence trap (step 1): confirmed `q_eta_micro` (lines 193–196) has free symbols `{c_etaU, c_etaU_ref, K_U, K_U_ref, K_eta_eff, K_eta_eff_ref}` only, so `diff(q_eta_micro, {zeta, M_tr, lambda_phi, K_phi_eff})` is identically 0 by absence — this is exactly F2. For the F2 fix I checked that making `K_eta_eff` (a dressing-sector stiffness) explicitly independent of `lambda_phi`/`K_phi_eff` while letting `q_eta` be built in a frame where those support symbols are in scope yields a genuine zero derivative; and that a counterfactual where a support symbol leaked into `q_eta` would make the assert fail (non-vacuity). Trivial-case (step 3): `q_eta.subs(eps_eta, eps_eta_ref) = log(1) = 0`, confirming the F4 replacement assert is true and non-trivial. Parity/integral trap is N/A (no integrals in this stage). Path spec (step 4): the `.wl` target is `mathematica/...stage237..._mathematica_audit.wl`. Paper round-trip (step 5): the F2/F3/F4 fixes introduce no new constants and re-use only paper-stated expressions, so no new paper_misalignment is created.
