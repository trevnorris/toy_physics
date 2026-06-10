---
unit_id: 237
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
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
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at line 86 and §"Physical logarithmic chart" lines 1010-1065, theorem 1169-1174)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_mathematica_audit.txt`

## What the paper claims

Stage 237 computes the surviving rigid-mouth dressing coordinate `q_eta` on the actual coherent branch and proves it is support-blind. `\stagefield{Output}` reads verbatim: *"Actual-branch dressing compiler and support-blindness theorem: explicit support-sector shifts do not by themselves fix the dressing coordinate."* The `\stagefield{Derivation ledger}` adds: *"The derivation derives `q_η=2c_1-κ_U-κ_η`, shows the finite static-blind curve is support-blind, and identifies the post-static orbit-lock correction."* The notes + Part VII appendix (lines 1010-1174) enumerate the concrete deliverables: (1) rigid-mouth finite packet `q_nt = ln((1-ε_η)/(1-ε_ηref)) - ln(R_target/R_target,ref)`, `q_eta = ln(ε_η/ε_ηref)`; (2) the finite static-blind curve `R_target/R_target,ref = (1-ε_η)/(1-ε_ηref)`, its `q_eta` parameterization, and the inverse `q_eta(R_target)`; (3) first-order packet + tangent slope `-c_eta` with `c_eta = ε_η,*/(1-ε_η,*)`; (4) microscopic compiler `q_eta = 2 ln(c_ηU/...) - ln(K_U/...) - ln(K_η^eff/...)`; (5) first-order drift extractor `q_eta = 2c_1 - κ_U - κ_η`; (6) support-blindness `∂_ζ q_eta = 0`, `∂_M_tr q_eta = 0`, and the sharper `∂_λφ q_eta = 0`, `∂_Kφ q_eta = 0`; (7) the post-static dependent-plane ray `-q_eta(0,1,1)` and the codimension-two orbit-lock point. Card status is `\StatusExactClosure`, checkpoint False.

## What the script claims to verify

The SymPy script (and its independent Mathematica twin) asserts, section by section: M1 the rigid-mouth reduction of the finite packet (`q_tr`->0, `q_nt` reduces correctly, `q_eta` unchanged); M2 that `q_nt=0` is equivalent to the static-blind curve, plus both round-trip inverses `q_eta<->Rratio` and the `q_eta=0` endpoint; M3 the first-order linearization, the tangent slope `= -c_eta`, the packet matrix `[[-1,-c_eta],[0,1]]` with det `-1`, and (SymPy) that the matrix maps `(R1,E1)` to the linear packet; M4 the microscopic log expansion identity; M5 the first-order drift coefficient `2c1-κ_U-κ_η`; M6 the four support-blindness derivatives vanish; M7 the dependent-plane ray equalities, the finite endpoints, and the codimension-two `(R1,E1)=(0,0)` solution. All checks are `assert_zero` / `expectZero` with `Exit[1]` on nonzero. Both saved outputs report all checks PASS.

## Paper ↔ script cross-check

| paper deliverable | script check | status |
|---|---|---|
| (1) rigid-mouth finite packet | M1 (py 89-116 / wl 108-126) | match |
| (2) static-blind curve + inverse `q_eta(R_target)` | M2 (py 121-146 / wl 130-159) | match |
| (3) first-order packet + tangent `-c_eta` | M3 (py 151-185 / wl 163-182) | match |
| (4) microscopic compiler | M4 (py 196-205 / wl 186-198) | match |
| (5) drift extractor `2c1-κ_U-κ_η` | M5 (py 207-220 / wl 202-216) | match |
| (6) support-blindness ∂_ζ,∂_M_tr,∂_λφ,∂_Kφ | M6 (py 230-294 / wl 218-268) | match |
| (7) dependent-plane ray + codim-two orbit-lock | M7 (py 299-320 / wl 270-300) | match |
| card Verification "Mathematica audit: none yet" | a full passing `.wl` exists | mismatch |

`paper_alignment: partial` — all seven math deliverables reconcile; the single discrepancy is the card's stale Verification field (paper-side text lag, not a math defect).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 101-109 | `assert_zero(q_tr_rm)`, q_nt/q_eta reductions | claim 1 | yes |
| A2 | sympy | 125-141 | static-curve + round-trips + endpoint | claim 2 | yes |
| A3 | sympy | 164-177 | linear terms, tangent slope, det, matrix map | claim 3 | yes |
| A4 | sympy | 205 | `q_eta_micro - q_eta_micro_split` | claim 4 | yes |
| A5 | sympy | 220 | drift coeff `2c1-κU-κη` | claim 5 | yes |
| A6 | sympy | 252-286 | presence guard + 4 support derivatives | claim 6 | yes |
| A7 | sympy | 303-314 | ray equalities, endpoints, codim-two | claim 7 | yes |
| B1-B7 | mathematica | 119-300 | engine-native mirrors of A1-A7 | claims 1-7 | yes |

All rows non-tautological: each `assert_zero` differences two *independently constructed* expressions (e.g. M2 builds `Rratio_static_blind` directly while solving `q_nt_ratio=0`; M4 builds `q_eta_micro` from the `c²/(K K)` ratio and `q_eta_micro_split` as a separate sum-of-logs; they only agree if log-algebra holds). M6 is the hotspot and is handled below.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (stale Verification status)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_237.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_mathematica_audit.wl` (full file, 304 lines, all checks PASS)

**What's wrong:**
The card states (line 11): *"Mathematica audit: none yet."* But a complete, independent, all-passing Mathematica audit `.wl` exists (mtime 2026-06-02 22:24) with a fresh saved output (2026-06-03 08:16) reporting every M1-M7 check PASS. The card understates the verification coverage.

**Why this matters:**
The card misreports the dual-engine status of an `\StatusExactClosure` stage; a reader auditing coverage would believe the Mathematica engine is missing. This is a documentation lag, not a math error, but for a non-checkpoint exact-closure stage the dual-engine status should be accurate.

**Required change (paper-side; user-gated):**
Update the `\stagefield{Verification}` line so "Mathematica audit: none yet" reflects the present `mathematica/...stage237..._mathematica_audit.wl`. Codex must NOT auto-edit the paper card — routes to user resolution (see directive `## Resolve before fix_loop`).

**Verification:**
Card line 11 should name the existing `.wl` file rather than "none yet."

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, NOT a transliteration. Evidence:
1. **Static-blind curve (M2):** SymPy substitutes `Rratio_static_blind` directly into `q_nt_ratio` (py 124-125). Mathematica instead poses `staticSolve = Solve[staticExp == 1, Rratio, Reals]`, checks uniqueness (`Length != 1 -> fail`), and compares the *solved* root against the closed form (wl 132-138). Different solving route.
2. **Codim-two orbit-lock (M7):** SymPy uses `M_packet.LUsolve(Matrix([0,0]))` (py 313). Mathematica uses `Solve[{-R1-cEta E1==0, E1==0},{R1,E1},Reals]` with a uniqueness guard (wl 293-296). Different mechanism.
3. **First order:** SymPy `sp.series(...).removeO()` + `M_packet*Matrix([R1,E1])` matrix-map check (py 156-177); Mathematica `Normal[Series[...]]` + `Coefficient[...,t,1]` (wl 169-211). Engine-native, not a 1:1 line port.
The shared structure (same log formulas, same `support_args` set) is unavoidable — both engines must encode the identical physical definitions; that is not transliteration. No `mathematica_transliteration` finding.

## Engine cross-check

Both saved outputs report all checks = 0 / PASS for M1-M7. No residual disagreement. `engines_agree: true`. The SymPy output additionally prints intermediate forms (`c_eta`, `M_packet`, `q_eta_micro`) consistent with the Mathematica forms (e.g. both give `c_eta = ε_ηref/(1-ε_ηref)` / `-ε_ηref/(ε_ηref-1)`, identical).

## Verdict justification

Verdict `findings` with one low-severity paper_misalignment (stale card Verification field). The math is sound and exactly aligned across both engines and the appendix. I attacked the known hotspot (M6 support-blindness self-test trap) hardest and it holds: the diffs are taken w.r.t. variables that genuinely appear (enforced by an explicit presence guard in both engines), the differentiated coupling responses are live abstract functions of the support args (so the diffs produce real `Derivative` terms, not trivial zeros), and the pass-1 negative-control + leak-detector guards are present and load-bearing in BOTH engines (see Self-test notes). The only defect is documentation: the card says "Mathematica audit: none yet" while a full passing `.wl` exists. That routes to user resolution; no Codex script edit is warranted.

## Self-test notes

- **Variable-independence trap (every D[]/diff):** Checked all five differentiations. (i) py167/wl175-178 `d ln(Rratio_of_qeta)/d qeta_param`: `Rratio_of_qeta` contains `exp(qeta_param)` -> variable present -> nonzero result `-c_eta` (not vacuous). (ii)-(v) py274-281 / wl256-264 support derivatives w.r.t. `zeta, M_tr, lambda_phi, K_phi_eff`: every one of these appears inside the abstract support functions `c_etaU_support(zeta,M_tr,lambda_phi,K_phi_eff)` etc. (directly, and via `zeta_expr`/`M_tr_expr` after the composite subs), so each diff is genuinely live. CONCLUSION: every D[]/diff differentiates w.r.t. a variable that actually appears — no vacuous derivative.
- **Negative-control / leak-detector guards present & load-bearing:** YES, in both engines and stronger than a comment. SymPy: explicit presence guard `set(support_args).issubset(q_eta_support_frame.free_symbols)` raises if any support var is absent (py 252-253); `impose_support_independence` zeros only the *support-function* `Derivative` atoms, so the asserts can pass only because the diffs actually produced live Derivative terms that the physical-blindness rule kills. Mathematica: TWO guards — `Not[FreeQ[qEtaSupport,#]]` presence guard (wl 231-234) AND the explicit leak-detector `Not[FreeQ[#, Derivative]]` over `rawSupportDerivatives` (wl 243-246) which `fail`s if any raw derivative is trivially Derivative-free (i.e. the channel is dead). The `supportBlindRules` then zero `Derivative[__][cEtaUSupport][__]` etc. This is exactly the pass-1 fix (live-channel negative control + `Not[FreeQ[#,Derivative]]`), present and load-bearing, not removed/weakened.
- **Symmetry/parity & trivial-case:** no unbounded integrals in this stage; checked endpoint substitutions (qeta_param->0, Rratio->1, eps_eta->eps_eta_ref) all give literal 0 as expected.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 7 deliverable values checked, 0 misaligned (1 doc status-field lag flagged separately as F1).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `q_nt = ln((1-ε_η)/(1-ε_ηref)) - ln(R_t/R_t,ref)`, `q_eta = ln(ε_η/ε_ηref)` | py 90-95 / out 5-8; wl 109-114 / out 9-12 | notes §1 boxed (l.77-87); appendix eq:app-part07-rm-finite-packet (l.1011-1018) | MATCH |
| static-blind curve `R_t/R_t,ref = (1-ε_η)/(1-ε_ηref)` + inverse `q_eta(Rratio)` | py 124-130 / out 11-13; wl 131-145 | notes §2 boxed (l.105-138); appendix eq:app-part07-static-blind-curve / -qeta-from-rtarget (l.1020-1043) | MATCH |
| `c_eta = ε_ηref/(1-ε_ηref)`, tangent slope `-c_eta`, packet matrix `[[-1,-c_eta],[0,1]]` | py 154-185 / out 16-29; wl 163-182 | notes §3 boxed (l.176-181, 201-207) | MATCH |
| micro compiler `q_eta = 2ln(c_ηU/..)-ln(K_U/..)-ln(K_η^eff/..)` | py 196-205 / out 32-37; wl 186-198 | notes §4.1 boxed (l.237-246); appendix eq:app-part07-qeta-micro-finite (l.1045-1053) | MATCH |
| drift extractor `q_eta = 2c1-κ_U-κ_η` | py 220 / out 38; wl 213-216 | notes §4.2 boxed (l.256-262); appendix eq:app-part07-qeta-first-order (l.1056-1057); card Derivation ledger (l.13) | MATCH |
| support-blindness `∂_ζ=∂_M_tr=∂_λφ=∂_Kφ q_eta = 0` | py 283-294 / out 43-46; wl 256-264 | notes §5 boxed (l.295-308); appendix eq:app-part07-qeta-support-blind (l.1061-1064) | MATCH |
| dependent-plane ray `-q_eta(0,1,1)` + codim-two orbit-lock | py 299-314 / out 49-56; wl 272-300 | notes §6 boxed (l.330-339, 360-366); appendix eq:app-part07-pure-dressing-image + thm (l.1144-1172) | MATCH |

INTERNAL scaffolding (no finding): `B_star`/`C_star` (carried packet coeffs, set to 0 on rigid-mouth slice), `zeta_expr`/`M_tr_expr` definitions (intermediate), all PASS flags, the `_numeric_zero` sample dictionaries, det(`M_packet`)=-1 (invertibility scaffolding), abstract support functions `c_etaU_support`/`K_U_support`/`K_eta_eff_support`.
</content>
</invoke>
