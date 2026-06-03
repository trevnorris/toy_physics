---
unit_id: 231
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 231 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_231.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (row at line 74; narrative lines 848-865; claimstatus line 881)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 231 ("Continuum-placement pullback") rewrites the Stage-247/248 selected-branch dynamic-class classifier `R_ND(xi,delta)` in the physical continuum-placement variables `(eps_eta, eps_W, rho, Z_W, delta_0, Lambda)`. The card's `\stagefield{Output}` states verbatim: "Continuum-placement version of the static-first verdict: the branch decision can be made directly from the physical placement packet, with $R_{\rm target}$ controlling the selected-class motion monotonically." The notes enumerate six deliverables: (1) the geometry/classifier functions `F`, `G`, `R_ND` with their exact monotonicities and endpoints; (2) the physical classifier `R_phys := R_ND(xi_req, delta)` with the implicit-derivative monotonicity `partial_{R_target} R_phys = (partial_xi R_ND)/(partial_xi F) < 0`; (3) the generic threshold polynomial `P_c(xi,delta)` with unique root and onset law `delta_c = 8/(9c)`; (4) the pulled-back sign-flip and denominator threshold curves `R_flip(delta)`, `R_den(delta)` on sample slices with collapse to onset; (5) the continuum placement map and product law `R_target * M_mix = 8 Lambda (1-eps_W)/pi^2`, plus equivalent kernel inequalities; and (6) the static-first theorem surviving pullback (`inf B_dyn > B_stat` for both robust and nonempty budgets). The appendix row (line 74) and narrative (848-865) summarize the same with the monotonicity theorem and static-first survival as the headline.

## What the script claims to verify

The SymPy script verifies, in seven blocks, that: the symbolic derivatives `dF/dxi`, `dG/dxi`, `dR_ND/dxi` match closed-form `expected_*` polynomials, with endpoint values `F(0)=1, G(0)=0, R_ND(0)=8/(9 delta)` and `F -> oo` as `xi -> 1^-` (block 1); the pullback compiler `dR_phys/dR_target = dR_ND/dF` is assembled symbolically and its sign is checked `<0` on a 6-point `(xi,delta)` grid along with the three single-derivative signs (block 2); the threshold polynomial `P_c`, its derivative, onset factor `9 delta^2 (9 c delta - 8)`, and onset law `delta_c = 8/(9c)` (block 3); the carried numeric thresholds `R_* = s_den/(-s_num) ~ 1.2292554`, `delta_dyn_* = 8/(9 R_*)`, `delta_den = 8/9` (block 4); numeric `R_flip`, `R_den` on the three sample slices via bisection plus collapse-to-onset above the onset deltas, with `R_flip <= R_den` (block 5); the continuum placement map and product law `R_target * M_mix = 8 Lambda (1-eps_W)/pi^2` plus the two bound translators (block 6); and the static-first inequalities `B_dyn > B_stat` as comparisons of four hardcoded literals (block 7).

## Paper <-> script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) geometry/classifier + monotonicities + endpoints | lines 67-116 (`expected_dF/dG/dR`, onsets, soft limit) | match (but notes formula for `dF/dxi` is mistyped — see F1) |
| (2) `R_phys` monotonicity `dR_phys/dR_target<0` | lines 121-155 (symbolic assembly + grid sign check) | match |
| (3) `P_c`, unique root, `delta_c=8/(9c)` | lines 160-175 | match |
| (4) `R_flip`, `R_den` curves + onset collapse | lines 198-236 (bisection vs sample rows, collapse checks) | match |
| (5) placement map + product law + kernel inequalities | lines 241-265 (product law assertion; bound translators printed only) | partial — product law asserted; kernel inequalities only printed, not asserted (acceptable, they are algebraic rearrangements) |
| (6) static-first survives pullback | lines 270-282 | partial — only compares 4 pre-baked literals; pullback itself not exercised (see F4) |

`paper_alignment: partial` — the script faithfully exercises deliverables 1-5; deliverable 6 is verified only as a literal-vs-literal inequality, and the notes carry a mistyped `dF/dxi` polynomial (F1).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 94 | `simplify(dF_dxi - expected_dF_dxi)==0` | claim 1 | yes (SymPy-derived `dF_dxi` vs closed form) |
| A2 | sympy | 95 | `simplify(dG_dxi - expected_dG_dxi)==0` | claim 1 | yes |
| A3 | sympy | 96 | `simplify(dR_dxi - expected_dR_dxi)==0` | claim 1 | yes |
| A4 | sympy | 102-105 | `onset_F==1`, `onset_G==0`, `onset_R-8/(9 delta)==0`, `soft_limit_F==oo` | claim 1 | yes |
| A5 | sympy | 146-153 | grid sign checks `dF>0, dG>0, dR<0, dRphys<0` | claim 2 | yes (sign on representative grid) |
| A6 | sympy | 163 | `simplify(dP_dxi - expected)==0` | claim 3 | yes |
| A7 | sympy | 166,169 | onset factor + `delta_c` zero check | claim 3 | yes |
| A8 | sympy | 186-188 | `assert_close(R_*, ...)`, `delta_dyn_*`, `delta_den` | claim 4 | partial (R_* derived from carried `s_minus`; the two `s_minus` literals are unanchored, see F2) |
| A9 | sympy | 219-224 | `assert_close` of `xi/R_flip/xi/R_den` to sample rows + `R_flip<=R_den` | claim 4 | yes (bisection independent of the sample literals it compares to) |
| A10 | sympy | 233-234 | collapse-to-onset `R_cap(...) == 1.0` | claim 4 | yes |
| A11 | sympy | 251 | `simplify(product_law - 8 Lam(1-eps_W)/pi^2)==0` | claim 5 | yes |
| A12 | sympy | 275-278 | `B_dyn_both > B_stat_both`, `B_dyn_nonempty > B_stat_nonempty` | claim 6 | no (compares 4 hardcoded literals, see F4) |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md:98`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py:81-85`

**What's wrong:**
The notes give the `partial_xi F` numerator polynomial (line 98) as
`81 delta^3 + 240 delta^2 xi + 72 delta^2 + 297 delta xi^2 + 189 xi^3`.
The script's `expected_dF_dxi` (lines 82-84) uses
`81 delta^3 + 189 delta^2 xi + 72 delta^2 + 297 delta xi^2 + 121 xi^3`.
Two coefficients disagree: the `delta^2 xi` term (notes `240` vs script `189`) and the `xi^3` term (notes `189` vs script `121`). The script's assertion `assert sp.simplify(dF_dxi - expected_dF_dxi) == 0` (line 94) passed (output line 14 shows SymPy's own factorization with `189*delta**2*xi` and `121*xi**3`), so the script's polynomial is the algebraically correct derivative of `F` and the notes line is a transcription error. This is `notes_contradicts_script` — the notes (paper side) disagree with the verified script.

**Why this matters:**
The notes are the authoritative derivation source for the card; a wrong closed-form derivative there is a latent error that could mislead a downstream reader or a future re-derivation, and it makes the notes <-> script alignment claim (notes line 627, "the note and script are therefore aligned") false on this formula.

**Required change:**
Routed to USER (paper-side prose). Codex must not edit notes. See directive `## Resolve before fix_loop`. Likely correct value: the script polynomial `189 delta^2 xi` and `121 xi^3`, since SymPy independently produced it.

**Verification:**
After the user authorizes the notes edit, notes line 98 should read `81\delta^3+189\delta^2\xi+72\delta^2+297\delta\xi^2+121\xi^3`, matching script lines 82-84 and output line 14.

### F2 — missing_verification_script (subtype missing_mathematica)

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/` (no `.wl` for unit 231)

**What's wrong:**
There is no Mathematica script for this unit. Every claim in this stage is well within Mathematica's native reach: symbolic differentiation/factoring of rational functions (`D`, `Factor`, `Together`) for `F`, `G`, `R_ND`, `P_c` and their derivatives; `Limit[F, xi -> 1, Direction -> "FromBelow"]` for the soft limit; `Reduce`/`FindRoot` for the threshold roots `xi_c` and the pulled-back curves `R_flip`, `R_den`; and direct evaluation of the product law and static-first inequalities. The dual-engine policy requires a `.wl` wherever Mathematica CAN independently verify (the test is "is it possible," not "is it necessary"), and it plainly can here.

**Why this matters:**
A single-engine stage has no cross-check; an error in the SymPy algebra (or in the `expected_*` closed forms it is compared against) would go undetected. Mathematica is the orthogonal second engine the framework requires.

**Required change:**
Codex writes a NEW independent-route Mathematica script `mathematica/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_mathematica_audit.wl` that re-derives the claims from the physical premises via native Mathematica primitives and a different decomposition than the `.py` (not a transliteration). Claim manifest in directive F2.

**Verification:**
`redteam exec-mathematica 231` produces output with all checks present and exit 0; results agree with the SymPy values (R_* ~ 1.2292554, delta_dyn_* ~ 0.7231116, sample `R_flip`/`R_den` rows, product law identity).

### F3 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:180-181`

**What's wrong:**
`s_minus_den = 0.411024574532864` and `s_minus_num = -0.334368725711457` (lines 180-181) are hardcoded literals with a comment "Carry the Stage 230 sign-flip and denominator thresholds." Stage 231's own card and notes never state these two numbers; they appear only as the upstream source of `R_* ~ 1.2292554` (notes line 152 states only `R_*`, not the `s_minus` pair). So the script reconstructs `R_*` from two numbers that have no anchor in this stage's paper material. This is acceptable carry-forward (it is genuinely an import from stage 230, and deriving `R_*` from them is non-tautological), but the literals are unanchored within unit 231.

**Why this matters:**
If the stage-230 values were ever revised, this stage would silently keep stale numbers. The risk is low because `R_*` is independently pinned by `assert_close(..., 1.229255438463336)` at line 186, so a wrong `s_minus` pair would be caught by that comparison.

**Required change:**
Add a one-line comment at line 180 citing the specific stage-230 source script/result that defines `s_minus_den`/`s_minus_num`, so the carry-forward is traceable. Do not change the values.

**Verification:**
A provenance comment naming the stage-230 source appears above line 180; script still exits 0.

### F4 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:270-278`

**What's wrong:**
Block 7 ("Static-first theorem survives the pullback") asserts `B_dyn_both_inf > B_stat_both` and `B_dyn_nonempty_inf > B_stat_nonempty` (lines 275-278), but all four quantities are hardcoded floats (lines 270-273). The inequality compares two pre-baked literals against two others; it cannot fail unless someone retypes a literal, and it exercises nothing about the pullback that the stage is supposed to prove. The paper's claim (notes Section 7, appendix Theorem) is that the pullback restricts `R_phys` to a subset of `[0,oo)` so the strict inequalities survive — the script does not connect the four numbers to the pullback at all; it merely re-states the four upstream values and compares them.

**Why this matters:**
This is the headline "static-first survives pullback" deliverable, and the script-side check is purely a numeric restatement of upstream constants. It gives false confidence that the pullback was exercised. It is low severity because the four values are themselves carried forward from a verified upstream stage (247) and the surviving-inequality argument is a set-inclusion logical step, not new algebra — but the check as written is not load-bearing.

**Required change:**
Add a provenance comment at lines 270-273 naming the upstream stage (247) script/result that defines these four budgets, and (optionally, if cheap) tie the survival statement to the pullback by asserting that the sampled `R_phys` values over the `probe_grid` lie within `[0, R_ND(0,delta)]` (the physical subset), so the inequality is connected to the pullback rather than free-floating. Do not change the four budget values.

**Verification:**
Provenance comments appear at lines 270-273; if the subset check is added, a new assertion bounding sampled `R_phys` appears and the script exits 0.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration is not assessable. F2 requires Codex to author an independent-route `.wl`; the directive's claim manifest mandates a different decomposition (e.g., `Reduce`/`Resolve` for the threshold roots and `Limit`/`Series` for endpoints) rather than echoing the `.py` choreography.

## Engine cross-check

Only one engine present. Not applicable until F2 lands.

## Verdict justification

Verdict is `findings`. The SymPy script is substantive and largely faithful: deliverables 1-5 are exercised with non-tautological assertions (the `dF/dxi`, `dG/dxi`, `dR_ND/dxi`, `dP_c/dxi` closed forms are checked against SymPy's own differentiation; the threshold roots are found by independent bisection and compared to sample rows; the product law is an honest symbolic identity). Attacks I tried that failed: I checked whether the grid sign checks (A5) could be vacuous — they evaluate genuine non-constant rational functions at six distinct points, so they exercise the signs. I checked whether the `R_flip <= R_den` check (A9) was tautological — it is not, since `R_flip` and `R_den` are computed by separate bisections at caps `R_*` and `1`. I checked the collapse-to-onset checks (A10) — the chosen deltas (0.80, 1.00) genuinely exceed the onset deltas (0.7231, 0.889), so onset-satisfaction is real. The verdict is not `clean` because: (F1) the notes carry a mistyped `dF/dxi` polynomial that disagrees with the verified script (user-routed); (F2) there is no Mathematica second engine although the stage is fully Mathematica-verifiable; (F3) two unanchored carry-forward literals; and (F4) the headline static-first block is a literal-vs-literal comparison that does not exercise the pullback. None of these are `UNFIXABLE` or `CRITICAL_DOWNSTREAM`: F1 is a prose typo, and the script's value is already independently pinned; F2/F3/F4 are additive/annotative. I confirm I read the card, notes, and appendix row before auditing the script.

## Self-test notes

Variable-independence trap: F2's claim manifest derivatives (`D[F,xi]`, `D[R_ND,xi]`, `D[P_c,xi]`) all depend on `xi`, so none collapse to identically zero. Symmetry/parity: no unbounded-domain integrals in this stage, so no parity trap applies. Trivial-case pre-check: the `dF/dxi` discrepancy (F1) was resolved by confirming SymPy's own factorization (output line 14) matches the script's `expected_dF_dxi`, identifying the notes as the erroneous side; I did not propose any new `assert_zero`/`assert_nonzero` whose residual I could not reduce by hand. Paper round-trip: F3/F4 prescribe only provenance comments and an optional subset bound (no constant changes), so they introduce no new paper_misalignment.
