---
unit_id: 243
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T09:30:47-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 243 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_243.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (row 84 + main-theorem items 1-5 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.txt`

## What the paper claims

Stage 243 (a checkpoint) declares the relaxed-constraint branch as a codimension-three lift of the Stage-242 front end, while hard-wiring that the same-charge long-range verdict is unchanged. The notes' Section 9 enumerates nine exact algebraic objects the audit must verify: (1) the Gaussian/odd-profile leakage scalar `S_leak = -(sqrt2/4) ell_w j0`; (2) the work channel `W_w = (sqrt(2 pi)/8) ell_w j0 E0`; (3) the non-rigid response `U = k_V f_U/(k_U k_V - chi_lam^2)`, `V = chi_lam f_U/(k_U k_V - chi_lam^2)`; (4) the drain `D_UV = chi_lam^2 k_V f_U^2/(k_U k_V - chi_lam^2)^2`; (5) the unit-mean source normalization `int_0^1 varsigma dz = 1`; (6) the quadratic rewrite `varsigma(y) = 1 - b + a y + 2 b y^2`; (7) interior stationary point `y_* = -a/(4b)` and vertex `varsigma(y_*) = 1 - b - a^2/(8b)`; (8) the recovery slice `(ell_w, f_U, a, b) = (0,0,0,0) => S_leak = W_w = U = V = D_UV = 0`, `varsigma == 1`; (9) the strict short-range limits `lim_{x->inf} x deltaV_stat = 0` and `lim_{x->inf} x Re V_dyn = 0` for the kernel span `{x^-6, e^{-2 kappa x}/x^4, e^{-4 kappa x}/x^2}`. The card additionally records the stability condition `k_U k_V - chi_lam^2 > 0` and the ratio `V/U = chi_lam/k_V`. There is no single `\stagefield{Output}` line; the card's `\stagefield{Verification}` names both scripts, and the appendix row 84 summarizes the deliverable as "codimension-three lift ... exact recovery map ... leakage/work lane, non-rigid (U,V) lane, compensated source lane, and strict short-range kernel-span theorem."

## What the script claims to verify

The SymPy docstring states five logical blocks: leakage/work lane, non-rigid U/V response + drain, compensated-source normalization, recovery slice, and the no-new-long-range limit. The assertions actually test all nine notes-level deliverables plus the stability determinant and the V/U ratio: section 1 integrates `W'(w) j^w(w)` from scratch and matches `-(sqrt2/4) ell_w j0`, integrates `j^w E_w` and matches the work scalar; section 2 solves the 2x2 stationarity system via `sp.solve` and compares to the closed forms, computes `det H`, `V/U`, the drain, and the f_U=0 / chi=0 sub-limits; section 3 integrates the mean, performs the `y = cos(pi z)` rewrite, computes the vertex; section 4 substitutes the recovery map; section 5 forms the source products and takes `x -> infinity` limits of `x` times each kernel and of the full corrections. The Mathematica script mirrors the same five blocks with `expectZero` residual checks.

## Paper <-> script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (1) S_leak = -(sqrt2/4) ell_w j0 | sympy L40,56 / wl L81,94 | match |
| (2) W_w = (sqrt(2 pi)/8) ell_w j0 E0 | sympy L41,57 / wl L85,95 | match |
| (3) U, V non-rigid response | sympy L76-99 / wl L103-121 | match |
| stability det = k_U k_V - chi_lam^2 > 0 | sympy L100 / wl L122 | match |
| V/U = chi_lam/k_V | sympy L101 / wl L123 | match |
| (4) drain D_UV | sympy L84,102 / wl L114,124 | match |
| (5) unit-mean int = 1 | sympy L117,136 / wl L133,148 | match |
| (6) quadratic rewrite | sympy L118-138 / wl L137-152 | match |
| (7) y_*, vertex | sympy L126-139 / wl L143-155 | match |
| boundary candidates 1+/-a+b | sympy L134 (print only) / wl L156-157 (asserted) | match (wl stronger) |
| (8) recovery slice | sympy L145-172 / wl L161-174 | match |
| (9) short-range limits + span | sympy L181-220 / wl L178-214 | match |

Every paper-side deliverable is covered by a non-tautological script-side check, in BOTH engines. `paper_alignment: aligned`. (Note: the appendix main-theorem item 4 renames the constants to `a_U, a_V, chi_UV, Delta_UV`; the per-stage card and notes use `k_U, k_V, chi_lam`, which is what both scripts use. Since the per-stage card is the authoritative source for this stage and it matches the scripts, this is a citation-layer rename, not a misalignment.)

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 55 | `boundary == 0` | leakage IBP boundary | yes |
| A2 | sympy | 56 | `simplify(S_leak - expected) == 0` | claim 1 | yes |
| A3 | sympy | 57 | `simplify(W_work - expected) == 0` | claim 2 | yes |
| A4 | sympy | 58-59 | `S_leak.subs(ell_w,0)==0`, `W_work.subs(ell_w,0)==0` | recovery (claim 8) | yes |
| A5 | sympy | 98-99 | `U_sol - U_expected == 0`, `V_sol - V_expected == 0` | claim 3 | yes |
| A6 | sympy | 100 | `det_H - (k_U k_V - chi_lam^2) == 0` | stability | yes |
| A7 | sympy | 101 | `ratio_VU - chi_lam/k_V == 0` | V/U ratio | yes |
| A8 | sympy | 102 | `drain_UV - drain_expected == 0` | claim 4 | yes |
| A9 | sympy | 103-106 | f_U=0, chi=0 sub-limits | recovery / sub-cases | yes |
| A10 | sympy | 136-139 | mean=1, trivial slice, quadratic rewrite, vertex | claims 5,6,7 | yes |
| A11 | sympy | 167-172 | recovery slice all zero, varsigma=1 | claim 8 | yes |
| A12 | sympy | 213-220 | source products, mode limits, full-correction limits | claim 9 | yes |
| B1 | math | 93-97 | `expectZero` boundary, S_leak res, W_w res, ell_w=0 | claims 1,2,8 | yes |
| B2 | math | 120-128 | U/V res, detH res, V/U res, drain res, sub-limits | claims 3,4 + stability | yes |
| B3 | math | 148-157 | mean, trivial, quadratic, y_*, vertex, both boundaries | claims 5,6,7 | yes |
| B4 | math | 169-174 | recovery slice | claim 8 | yes |
| B5 | math | 207-214 | source products + limits | claim 9 | yes |

All rows are non-tautological under attack: the integrals, the linear solve, and the limits are genuinely computed in the SymPy engine, then compared to closed forms; none are `expr - expr` self-checks, and none re-substitute a self-solved root into the equation that produced it.

## Findings

### F1 — mathematica_transliteration

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl:73-214`

**What's wrong:**
The `.wl` is a near line-by-line port of the `.py`, not an independent second-engine derivation. The correspondence is one-to-one across all five blocks, with the same variable choreography, the same intermediate decomposition, and the same hardcoded `expected*` closed forms checked as residuals:

- Leakage: SymPy `S_leak = sp.integrate(sp.diff(W, w) * j_w, (w, -oo, oo))` (py:40) is mirrored verbatim by `Sleak = Integrate[D[W, w] jw, {w, -Infinity, Infinity}]` (wl:81-84), and BOTH hardcode the identical answer — py:43 `expected_S_leak = -sp.sqrt(2)*ell_w*j0/4` vs wl:90 `expectedSleak = -Sqrt[2] ellW j0/4` — then compare residuals.
- Non-rigid lane: SymPy `sol = sp.solve(stationarity, [U, V])` (py:76) with `U_expected = sp.simplify(k_V*f_U/(k_U*k_V - chi_lam**2))` (py:79) is mirrored by `uvSol = First[Solve[stationarity == {0,0}, {U,V}]]` (wl:103) with `Uexpected = ... kV fU/(kU kV - chiLam^2)` (wl:106). Same solver call, same hardcoded expected form, same residual comparison.
- Compensated source: the `y = cos(pi z)` substitution rule `{Cos[Pi z] -> y, Cos[2 Pi z] -> 2 y^2 - 1}` (wl:138-141) is a direct transcription of the SymPy `.subs({sp.cos(sp.pi*z): y, sp.cos(2*sp.pi*z): 2*y**2 - 1})` (py:119-124).
- Short-range: `deltaVStat = -(1/2)(C6 QQ + 2 C4 QY + C2 YY)` (wl:184) is identical in structure to py:188, and the same three `Limit[x * mode, x -> Infinity]` calls are taken (wl:186-205 vs py:191-198).

Because both engines (a) decompose the problem the same way and (b) check residuals against the SAME hardcoded closed forms, they cannot cross-validate the algebra — if the hardcoded form were wrong, both would simply mark a different residual and could pass or fail together. The second engine adds no independent confirmation. This violates the dual-engine policy (the `.wl` must derive the result from the physical premises by a different route, e.g., native primitives / a different decomposition / back-substitution, not echo the `.py`).

A corroborating artifact: the `.wl` banner (wl:59) reads `STAGE 226`, confirming the file was cloned from a stage-226 template rather than authored independently against the stage-243 physics.

**Why this matters:**
A checkpoint stage gets a higher bar; the second engine exists precisely to catch a wrong hardcoded `expected*` form. As written, a transcription error in either engine's hardcoded constant would not be caught by the other engine, defeating the cross-check. The transliteration silently downgrades a "two independent CAS confirmations" claim to "one CAS run twice."

**Required change:**
Re-author the `.wl` to an independent route that does NOT hardcode the same `expected*` closed forms the SymPy script uses. See the directive's claim manifest (M1-M5) for the specific independent verification each block must perform (back-substitution into the defining equations, native LinearSolve, integration-by-parts cross-check, and asymptotic-series / numeric-sample confirmation of the limits). Also correct the banner label `STAGE 226` -> `STAGE 243` (wl:59).

**Verification:**
The re-authored `.wl` must (i) no longer contain the literals `expectedSleak = -Sqrt[2] ellW j0/4`, `Uexpected = ... kV fU/(kU kV - chiLam^2)`, etc. as the sole reference forms; (ii) confirm each block by an independent mechanism (e.g., `Sleak` recovered by IBP `boundary + Sleak == Integrate[D[W jw, w], ...]`; `Usol, Vsol` confirmed by substituting back so `stationarity /. uvSol == {0,0}` AND independently via `LinearSolve[{{kU,-chiLam},{-chiLam,kV}}, {fU,0}]`); (iii) still exit 0 with all checks passing; (iv) banner reads STAGE 243.

### F2 — stale_output (banner label artifact)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl:59`

**What's wrong:**
The Mathematica banner prints `STAGE 226 — RELAXED-CONSTRAINT BRANCH DECLARATION ...` (wl:59), and the captured output reproduces this stale label (mathematica output txt:11). The math content is all stage-243 correct; only the human-readable banner is mislabeled. This is a copy artifact (consistent with the F1 transliteration origin) and matches the known project-wide "stale Stage NNN label" pattern. Not load-bearing for any assertion, but it is a real defect in a checkpoint file and should be corrected.

**Why this matters:**
A mislabeled banner in a checkpoint audit log is misleading to a reader auditing the ledger and obscures which stage the output belongs to. Saved outputs (`.txt`) are otherwise fresh (mtime 12:52 / 13:26 vs script mtime 11:58), so the only stale element is the label, not the math.

**Required change:**
Change wl:59 banner string from `STAGE 226` to `STAGE 243`. (Folded into F1's re-author; listed separately so the verifier confirms the label even if F1's route work is staged.)

**Verification:**
Banner line in regenerated `.wl` output reads `STAGE 243 — RELAXED-CONSTRAINT BRANCH DECLARATION AND SHORT-RANGE OPEN-SYSTEM COMPILER`.

## Independent-derivation check (Mathematica)

The `.wl` is NOT an independent derivation. Quoted correspondences (see F1): the leakage integral, the U/V `Solve`, the cosine-substitution rule, and the short-range `Limit` calls are each a syntactic re-spelling of the SymPy lines, and both engines compare against the SAME hardcoded `expected*` constants. Verdict: `mathematica_transliteration`.

## Engine cross-check

Both engines exit 0 and report PASS on every block; the SymPy output prints the explicit closed forms (e.g., `S_leak = -sqrt(2)*ell_w*j0/4`, `U solution = -f_U*k_V/(chi_lam**2 - k_U*k_V)`), and the Mathematica output prints `= 0` residuals with PASS for the corresponding residual checks. The numerical/symbolic content agrees. There is NO `engine_disagreement`. The agreement is, however, of low evidentiary value precisely because of F1 (the agreement is structural echo, not independent confirmation).

## Verdict justification

The math is correct and fully aligned with the paper: I attacked every assertion for tautology (the integrals, the `solve`, and the limits are genuinely computed, not `expr - expr`), checked the U/V signs and the drain denominator (sign-irrelevant, squared), verified the work-channel constant `sqrt2*sqrt(pi)/8 = sqrt(2 pi)/8` matches the card, confirmed the odd*odd integrand parity makes `S_leak` and `W_w` legitimately nonzero, and confirmed the recovery slice and short-range limits are substantive. No `paper_misalignment` (the appendix constant rename is a citation-layer cosmetic, the per-stage card matches the scripts). The verdict is `findings`, not `clean`, solely because the second engine is a transliteration mirror (F1) on a checkpoint stage that demands two genuinely independent CAS routes, plus the minor stale banner (F2). Neither finding changes the math, so no `stop_cold`.

## Self-test notes

Checked the variable-independence self-test trap: no `D[EXPR, VAR]`/`diff` is taken w.r.t. a variable the expression lacks — all derivatives are genuine stationarity `D[FUV, U/V]` and `D[W, w]`, which do depend on their variables, so no vacuous identically-zero check. Checked integral parity: `W' j^w` is odd*odd = even (nonzero integral, consistent with nonzero S_leak), `j^w E_w` is odd*odd = even (nonzero W_w) — both correctly nonzero, not falsely claimed to vanish. Trivial-case pre-check: recovery slice `(ell_w,f_U,a,b)=0` reduces every lifted observable to 0 and `varsigma` to 1 as asserted; `x->infinity` limits of `x` times each `x^{-n>=2}` and Yukawa kernel are 0 as asserted. The F1 directive prescribes independent routes (back-substitution, LinearSolve, IBP, asymptotics) that I verified reduce to the same correct values without re-using the .py's hardcoded constants.
