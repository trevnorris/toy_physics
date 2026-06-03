---
unit_id: 226
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 226 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_226.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows: line 64 status row; lines 681-745 the Strict-5PN-even-gate section + Theorem `thm:app-part07-strict-even-package`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The card (`stage_226.tex`) and appendix Theorem `thm:app-part07-strict-even-package` state that, inside the imported first-order 5PN package, the same-charge load scalar is exactly `\Xi_{\rm load}=\Xi_1=P_1/P_0` (card Output line 15; derivation-ledger line 13; appendix eq:app-part07-xi-load lines 692-701). The strict even gates `K_1:=D_{21}+D_{01}/9` and `H_{\rm even}:=D_{41}-\tfrac23 D_{21}-D_{01}/27` (appendix eq:app-part07-even-gates lines 703-708; notes lines 26-30) are separate conservative gates that must also be controlled. The notes enumerate five deliverables: (1) the exact bridge `\Xi_{\rm load}=\Xi_1=P_1/P_0`; (2) the exact comparison of the Stage-242/225 compensation surface against the strict 5PN gates, yielding `K_1=(1/9-u_2)D_{01}`, `H_{\rm even}=(D_4/D_0+\tfrac23 u_2-1/27)D_{01}` and the one-pole form `H_{\rm even}=(-3u_2^2+\tfrac23 u_2-1/27)D_{01}`; (3) the concrete compatibility-branch coefficients `K_1\approx 0.0621939470719309\,D_{01}`, `H_{\rm even}\approx -0.0116042611571584\,D_{01}` (notes 199-204); (4) the strict mixed-only even-gate `2x5` matrix (rank 2, nullity 3) with displayed null basis, the `\Xi_1(w_i)` values, and `\sigma_{\rm even}\approx 2.67386816837173`; (5) the pure-transfer `3x5` intersection (rank 3, nullity 2) with basis, `\Xi_1(t_i)`/`N_{01}(t_i)` values, `\sigma_{\rm transfer}\approx 2.31561904386057`, plus the transported Stage-241 ceiling budgets divided by each `\sigma`. The sharpest surviving subcorridor is the pure-transfer corridor `D_{01}=D_{21}=D_{41}=0`, `\Xi_1=\Xi_{\rm load}=N_{01}/N_0` (appendix eq:app-part07-pure-transfer-corridor-basic lines 732-737).

## What the script claims to verify

The SymPy script (docstring/print line 26) verifies all five notes deliverables. Section 1 (lines 35-67) builds `P_0A=N_0A/D_0A`, derives `P_1` as `\partial_\eps P_0A|_0/\lambda`, forms `Xi1=P1/P0`, and asserts it equals `N_{01}/N_0-D_{01}/D_0` (line 46); it then substitutes the compensation surface `D_{21}=-u_2 D_{01}`, `D_{41}=(D_4/D_0)D_{01}` into `K_1`, `H_{\rm even}` and asserts the closed forms (lines 59-61). Section 2 (lines 86-188) rebuilds the one-pole/hidden-even bundle (B/Z/N/P primitives), forms an `\eps`-dressed copy, differentiates to get the primitive compiler `D_{01},D_{21},D_{41},N_{01}`, evaluates on the explicit sample point, solves for the compatibility `K`, and `assert_close`s the resulting `D_0,D_2,D_4,u_2,u_4,D_4/D_0,P_0` against the notes' base data plus the one-pole identity `u_4-4u_2^2=0` (lines 171-179). Section 3-6 build the mixed-restricted `K_1`/`H_{\rm even}` and `D_{01}/D_{21}/D_{41}` coefficient matrices, check ranks/nullities, compare the null bases, `\Xi_1`/`N_{01}` projections, and the QR-projector corridor norms `\sigma_{\rm even}`/`\sigma_{\rm transfer}` against the notes, and divide the four carried ceiling budgets by each `\sigma`. All literals match the notes verbatim.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Bridge `\Xi_{\rm load}=\Xi_1=P_1/P_0` | lines 35-46 (`assert Xi1-Xi_load==0`, P1 derived) | match |
| `K_1=(1/9-u_2)D_{01}`, `H_{\rm even}=(D_4/D_0+\tfrac23 u_2-1/27)D_{01}`, one-pole `(-3u_2^2+\tfrac23 u_2-1/27)D_{01}` | lines 48-61 (asserts all three) | match |
| Compatibility base data `D_0,D_2,D_4,u_2,u_4,D_4/D_0,P_0`, `u_4=4u_2^2` | lines 152-179 (derived from bundle, asserted) | match |
| Concrete coeffs `K_1\approx0.06219...D_{01}`, `H_{\rm even}\approx-0.01160...D_{01}` | lines 204-208 | match |
| Strict even-gate `2x5` matrix, rank 2 / nullity 3, null basis `w_i`, `\Xi_1(w_i)`, `\sigma_{\rm even}\approx2.67386...` | lines 217-255 | match |
| Pure-transfer `3x5` matrix, rank 3 / nullity 2, basis `t_i`, `D_{01}=D_{21}=D_{41}=0`, `\Xi_1(t_i)`, `N_{01}(t_i)`, `\sigma_{\rm transfer}\approx2.31561...` | lines 270-318 | match |
| Transported ceiling budgets on both corridors | lines 323-356 (budget/`\sigma` division) | match |
| Independent Mathematica re-derivation | (no `.wl` exists) | missing |

`paper_alignment: aligned` — every paper-side deliverable has a faithful, non-tautological SymPy check, and every script literal matches the notes/appendix verbatim. The only gap is the absent second engine.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 46 | `simplify(Xi1 - Xi_load) == 0` | bridge `\Xi_{\rm load}=\Xi_1` | yes |
| A2 | sympy | 59 | `simplify(K1_comp - (1/9-u2)*D01) == 0` | even-gate `K_1` closed form | yes |
| A3 | sympy | 60 | `simplify(H_even_comp - (...)*D01) == 0` | even-gate `H_{\rm even}` closed form | yes |
| A4 | sympy | 61 | `simplify(H_even_one_pole - (-3U2^2+2/3 U2-1/27)*D01) == 0` | one-pole `H_{\rm even}` form | yes |
| A5 | sympy | 171-178 | `assert_close(derived, notes literal)` x8 | compatibility base data (deliverable 1 inputs) | yes |
| A6 | sympy | 179 | `simplify(u4_s - 4 u2_s^2) == 0` | one-pole identity | yes |
| A7 | sympy | 207-208 | `assert_close(K1/H coeff, notes literal)` | concrete branch coeffs (deliverable 2) | yes |
| A8 | sympy | 226-228 | `assert_close(M_even[i,j], expected)` | strict even-gate matrix | yes |
| A9 | sympy | 230-232 | `rank==2`, `len(null)==3` | nullity-3 corridor survives | yes |
| A10 | sympy | 239-248 | basis `w_i` + `\Xi_1(w_i)` match notes | corridor basis / same-charge values | yes |
| A11 | sympy | 255 | `assert_close(sigma_even, 2.67386...)` | corridor norm | yes |
| A12 | sympy | 276-278 | `rank==3`, `len(null)==2` | pure-transfer nullity-2 | yes |
| A13 | sympy | 284-297 | basis `t_i`, `\Xi_1(t_i)`, `N_{01}(t_i)` match notes | pure-transfer basis/values | yes |
| A14 | sympy | 300-303 | `(coeff(D01/D21/D41).T * b) ~ 0` | pure-transfer `D_{01}=D_{21}=D_{41}=0` | partial (redundant-by-construction; see note below) |
| A15 | sympy | 310 | `assert_close(sigma_transfer, 2.31561...)` | subcorridor norm | yes |
| A16 | sympy | 353-356 | `assert_close(budget/sigma, notes literal)` x8 | transported ceiling budgets | yes |
| (M) | mathematica | — | none | all deliverables | missing |

Note on A14 (lines 300-303): the basis `b` is `M_transfer.nullspace()` whose rows ARE `coeff_vector(D01_mixed)`, `coeff_vector(D21_mixed)`, `coeff_vector(D41_mixed)` at the same 30-digit precision (line 20 vs lines 271-273). Re-dotting those same coefficient vectors with `b` is `M_transfer*b`, which `nullspace()` guarantees is zero by construction. A14 is therefore a redundant confirmation, not an independent test. It is harmless because the load-bearing pure-transfer check is A13 (the basis must equal the notes' externally-stated `expected_t`), which is substantive. Not raised as a blocking finding.

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` (no `moving_throat_pde_stage226_*` file)

**Subtype:** `missing_mathematica`

**What's wrong:**
There is no Mathematica audit script for stage 226. The card itself records "Mathematica audit: none yet" (`stage_226.tex:11`). This is a `checkpoint: False`, `is_status_only_candidate: False` unit. The entire derivation is mechanizable independently in Mathematica: the bridge identity (rational `Series`/`Together` differentiation of `N_0A/D_0A`), the even-gate closed-form substitutions, the one-pole bundle compiler (`D[...,eps]` of the `\eps`-dressed primitives), the sample-point base data, the rank/nullspace of the `2x5` and `3x5` mixed matrices (`MatrixRank`, `NullSpace`), the same-charge functional projections, the QR/orthonormal-projector corridor norms (`Orthogonalize`/`Norm`), and the budget divisions. Per the dual-engine rule (the test is "is it POSSIBLE for Mathematica to verify," not "is it necessary"), a `.wl` is required here. A sibling Part-VII unit (stage 221) already carries a `.wl`, confirming both the policy and the directory convention.

**Why this matters:**
SymPy is the only witness for every load-bearing claim in stage 226, including the rank/nullity survival statements and the two corridor norms that feed downstream Part-VII budgets. A SymPy-only transcription error (e.g., in the primitive compiler differentiation, the QR-projector construction, or a coefficient extraction) would go undetected. The second engine must re-derive from the physical premises with native primitives and a different decomposition (NOT a line-by-line port of the `.py`).

**Required change:**
Add an independent Mathematica audit at the exact target path in the directive, verifying the M1-M7 claim manifest below using native Mathematica primitives and a distinct decomposition (anti-transliteration guard in the directive).

**Verification:**
`redteam exec-mathematica 226` runs the new `.wl` to exit 0 with all `expectZero`/`expect-close` checks passing; the captured `.txt` output reproduces `\sigma_{\rm even}\approx 2.67386816837173`, `\sigma_{\rm transfer}\approx 2.31561904386057`, ranks 2/3, nullities 3/2, and the bridge/even-gate identities.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot be assessed. The directive's anti-transliteration guard requires the new script to (a) compute the bridge `\Xi_1` by `Series`/`Normal` rational expansion rather than the SymPy `diff/subs(eps,0)/lam` choreography, (b) build the projectors via `Orthogonalize`+`Transpose` (or `PseudoInverse`-based) rather than mirroring `QRdecomposition`, and (c) compute nullspaces with `NullSpace` directly on freshly extracted `CoefficientArrays`, not by re-implementing `coeff_vector`.

## Engine cross-check

Only one engine present; not applicable. The SymPy saved output (`...sympy_audit.txt`, mtime May 11 12:50, newer than the script's May 11 11:58) reports `EXIT_CODE: 0`, `Status: PASS`, and printed values (`sigma_even = 2.67386816837172775`, `sigma_transfer = 2.31561904386055442`, ranks/nullities 2/3 and 3/2) consistent with the in-script `assert_close` targets. Output is fresh; no `stale_output` finding.

## Verdict justification

The SymPy script holds up against the paper under attack. I verified the bridge `\Xi_1=N_{01}/N_0-D_{01}/D_0` algebraically by hand (it reduces correctly from `P_1/P_0`), confirmed the even-gate substitutions reproduce the notes' closed forms (using `u_2=-D_2/D_0` and `u_4=4u_2^2`), confirmed the primitive compiler derivatives are non-trivial (the dressed bundle genuinely depends on `eps` through the mixed slopes, so the `D_{01}/D_{21}/D_{41}/N_{01}` rows are not identically zero — consistent with rank 2 and rank 3), and cross-checked all base data, matrix entries, null bases, `\Xi_1`/`N_{01}` values, the two corridor norms, and all twelve budget literals against the notes and appendix — every value matches verbatim. I specifically looked for the sibling-stage 222/223 notes-side additive-offset typo (notes = script + 68) and found NO such offset anywhere in stage 226: notes and script agree exactly on every literal. The only defect is the absent second engine, a `missing_mathematica` finding. Verdict is `findings` (not `clean`) solely because the dual-engine rule requires a `.wl` that does not exist; there is no `paper_misalignment` and no stop-cold condition.

## Self-test notes

- Trap 1 (variable independence): checked `D01_compiler=diff(Kd-B0d-Z0d,eps)` and the `D21/D41/N01` compilers — each dressed expression genuinely depends on `eps` via the slope exponentials, and after the mixed restriction the surviving rows are nonzero (rank 2 / rank 3 confirm this), so no identically-zero derivative.
- Trap 2 (parity): no unbounded symmetric-domain integrals in this stage; the work is finite linear algebra, so parity-vanishing traps do not apply.
- Trap 3 (trivial-case): the redundant A14 re-check (lines 300-303) reduces to `M_transfer*b=0`, which `nullspace()` guarantees — noted as redundant-not-independent, not a blocking finding because the substantive `expected_t` basis match (A13) carries the claim. The base-data and corridor-norm asserts are derive-then-compare against external notes literals, not self-comparisons.
- Provenance-label drift (informational only, no finding): script comments cite upstream "Stage 223/224/225" (lines 70, 181, 202, 321) while the notes/appendix cite "Stage 240/241/242" and the notes self-reference "Stage 243" (notes lines 5, 48, 60, 101, 122, 429, 477) — a renumbering-history label mismatch. No numeric value depends on the label; the carried literals (base data, ceiling budgets) match exactly, so this is not a `paper_misalignment` and is out of scope for Codex (it is prose, not a script value defect).
