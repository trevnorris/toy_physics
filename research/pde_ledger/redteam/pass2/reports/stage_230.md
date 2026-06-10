---
unit_id: 230
batch: VII.1
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
  notes_stage_files: [moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 230 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_230.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows 72, 832-846)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.txt` (present, fresh)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_mathematica_audit.txt` (present, fresh)

## What the paper claims

Stage 230 translates the Stage-229 selected-branch classifier `R_ND(xi,delta) = 72 delta^2 (1-xi) / ((9 delta + 11 xi)(9 delta^2 + 18 delta xi + 11 xi^2))` into a dynamic-window decision rule via an affine rigid-split compiler built on the four carried Stage-228 per-unit-Xi slopes. The card `\stagefield{Output}` reads verbatim: "Static-first theorem: the static $\Xi_1$ budget is the first kill condition; the dynamic corridor can be checked only after the static transported ceiling is satisfied." The notes enumerate the full deliverable set: (1) the share weights `w_num = R_ND/(1+R_ND)`, `w_den = 1/(1+R_ND)`; (2) the affine compiler `S_±(R_ND) = (R_ND s_±^num + s_±^den)/(1+R_ND)`; (3) the sign-flip threshold `R_* = s_-^den/(-s_-^num) ≈ 1.229255438463336`; (4) the onset threshold `delta_*^dyn = 8/(9 R_*) ≈ 0.723111617875019`; (5) the dynamic ceilings in `|eps Xi_1|`, their endpoints (`B_both(0)≈1.671064893775584`, `inf B_both≈0.967282389363822`, `inf B_nonempty≈0.990581810705233`) and sample values; (6) the universal static budgets (`B_stat^both≈0.367930328492646`, `B_stat^nonempty≈0.737619063660757`); (7) the global strict inequalities `B_dyn > B_stat`. The card's `\stagefield{Verification}` line states "Mathematica audit: none yet."

## What the script claims to verify

Both scripts verify the same seven-block program. M1/§1: the classifier formula, its onset `R_ND(0,delta)=8/(9 delta)`, soft limit `xi->1^-` = 0, and strict monotonicity in xi. M2/§2: the affine share compiler `S_±` and that both slopes decrease in `R`. M3/§3: the sign-flip threshold `R_*` as the unique nonnegative root of `S_- = 0`, plus the positive/negative sign split below/above it. M4/§4: the onset threshold `delta_*^dyn`, recovered by solving `R_ND(0,delta) = R_*`. M5/§5: the dynamic-ceiling formulas, endpoints, and four sample classifier rows. M6/§6: the static budgets and the two strict static-first inequalities `inf B_dyn > B_stat`. The assertions are numeric (`assert_close`/`expectClose`) and symbolic (`simplify(...)==0`/`expectZero`, plus `Resolve[ForAll[...]]` universal proofs on the Mathematica side).

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| share weights w_num/w_den | py L85-86, wl L79-80 (+ affine form asserts) | match |
| affine compiler S_±(R_ND) | py L87-93, wl L81-91 | match |
| sign-flip R_* ≈ 1.229255438463336 | py L111-112,118-121; wl L102-122 | match |
| onset delta_* ≈ 0.723111617875019 | py L114-115,119-121; wl L124-132 | match |
| dynamic ceilings + endpoints + samples | py L144-204; wl L134-213 | match |
| static budgets B_stat | py L209-213; wl L218-232 | match |
| static-first inequalities | py L212-213; wl L231-232 | match |
| Output: static-first theorem | py L226-228, wl L231-237 (inf B_dyn > B_stat) | match |
| card Verification: "Mathematica audit: none yet" | a passing `.wl` now EXISTS | mismatch (F1) |

`paper_alignment: partial` — all math deliverables align; the lone defect is the card's stale "none yet" Mathematica-verification claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 61 | `simplify(dR_dxi - expected) == 0` | claim 1 (classifier deriv form) | yes |
| A2 | sympy | 64-66 | onset = 8/(9δ); soft limit = 0 | claim 1 | yes |
| A3 | sympy | 82-83 | rigid slope sign guard | claims 2/3 (inputs) | yes |
| A4 | sympy | 92-93 | `S_± - affine form == 0` | claim 2 | yes |
| A5 | sympy | 97-98 | `dS_± < 0` | claim 2/5 monotonicity | yes |
| A6 | sympy | 112,115 | R_*, δ_* close to literal | claims 3,4 | partial (vs literal) |
| A7 | sympy | 118 | `S_-(R_*) == 0` | claim 3 | yes |
| A8 | sympy | 119-121 | `solve(Eq(onset,R_*),δ)` = δ_* (F2 de-taut) | claim 4 | yes |
| A9 | sympy | 124-135 | rep sign checks + δ=1 slice probe | claim 3 | yes |
| A10 | sympy | 151-152,163-172,191-198 | ell_±, ceiling endpoints, samples | claims 5 | yes |
| A11 | sympy | 212-213 | `B_*_inf > B_stat_*` | claim 7 / Output | yes |
| M1 | math | 64-71 | `Resolve[ForAll[...]]` monotonicity + strip equivalence | claim 1 | yes (stronger than py) |
| M2 | math | 90-99 | affine form + `Resolve[ForAll[D S_± < 0]]` | claim 2 | yes |
| M3 | math | 111-122 | unique root, R_*, `Resolve[ForAll]` sign split | claim 3 | yes |
| M4 | math | 130-132 | unique onset root via Solve = δ_* | claim 4 | yes |
| M5 | math | 176-213 | ell_±, endpoints, 4 samples | claim 5 | yes |
| M6 | math | 216-232 | `Resolve[ForAll[B' < 0]]` + static-first ineqs | claims 6,7 | yes |

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_230.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_mathematica_audit.wl` (entire file, present + passing)

**What's wrong:**
The card's verification line states verbatim:
> `\stagefield{Verification}{SymPy audit: \StageFile{...sympy_audit.py}.  Mathematica audit: none yet.}` (stage_230.tex:11)

But a Mathematica audit `.wl` now exists, is committed (added in commit 1dfc3fe, batch VII.1), and passes all checks (its saved output ends "All Stage 230 Mathematica audits passed."). The card therefore understates the verification coverage: it claims the second engine is absent when it is present and green.

**Why this matters:**
A reader citing this stage's verification status would conclude it is single-engine when it is in fact dual-engine. This is a paper-side prose claim about script coverage; under the v2 rules the resolution direction (update the card vs. remove the `.wl`) is the user's call, not Codex's. Given the dual-engine policy, the intended direction is almost certainly to update the card, but the auditor does not decide this.

**Required change:**
See `## Resolve before fix_loop` in the directive. Codex applies nothing until the user chooses a direction.

**Verification:**
After user resolution, if direction (a), the card line at stage_230.tex:11 should read "Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage230_..._mathematica_audit.wl}." and `git grep "none yet"` for stage 230 should be empty.

## Independent-derivation check (Mathematica)

The `.wl` is a **genuinely independent route**, not a transliteration. Three corresponding sections justify this:

1. **Monotonicity of the classifier.** The `.py` never proves `dR_ND/dxi < 0` over the domain — it only checks the derivative matches a pre-written closed form (`assert sp.simplify(dR_dxi - expected_dR_dxi) == 0`, py L61) and separately checks the *compiler* slopes `dS_± < 0` (py L97-98). The `.wl` instead runs `Reduce[0 <= xi < 1 && delta > 0 && D[rND,xi] < 0, {xi,delta}, Reals]` and then `Resolve[ForAll[{xi,delta}, (0<=xi<1 && delta>0) => rNDMonotoneReduce], Reals]` plus a `\[Equivalent]` strip-recovery proof (wl L59,64-71). This is a universal-quantifier proof over the full 2D stable strip — a different operation than the `.py`, and strictly stronger (the `.py` does not prove classifier monotonicity at all, only quotes its derivative).

2. **Compiler-slope monotonicity.** `.py`: `dS_plus = sp.diff(S_plus,R); assert dS_plus < 0` (a single-point/symbolic-sign check on a rational with constant numerator, py L95-98). `.wl`: `Reduce[D[sPlus,R] < 0 && R >= 0, R, Reals]` then `Resolve[ForAll[R, sPlusReduce \[Equivalent] R >= 0], Reals]` (wl L83,92-99) — proves the negativity region is exactly `R >= 0` rather than merely confirming a sign.

3. **Sign split / onset root.** `.py`: representative numeric substitutions at R=0,1/2,2 (py L124-129) + `solve(Eq(onset,R_*),δ)` (py L119). `.wl`: `Reduce`/`Solve` for the unique root, plus `Resolve[ForAll[R, (0<=R<rStar) => sMinus>0]]` and `Resolve[ForAll[R, R>rStar => sMinus<0]]` (wl L115-122) — universal proofs over the half-line rather than three sample points.

The shared elements — the classifier formula, the four Stage-228 slope rationals, the RQ figures, and the static budgets — are legitimate carry-forward INPUT constants, not echoed algebra; both engines are entitled to import the same upstream data. The verification choreography differs in kind (`Resolve[ForAll]` universal proofs vs. `.py` point-checks). **INDEPENDENCE CALL: independent.**

## Engine cross-check

Both engines agree to full precision. Spot comparison:
- `R_*`: py `1.22925543846333598788938878199` (output L17) vs wl `1.2292554384633359878893887819926314141068918213381142792773` (output L105) — identical to 29 digits.
- `delta_*`: py `0.723111617875019116300604645303` (output L18) vs wl `0.723111617875019116300604645303424784893...` (output L106) — identical.
- `inf B_both`: py `0.9672823893638217` (output L30) vs wl `0.9672823893638217522502832045...` (output L107) — agree.
- `inf B_nonempty`: py `0.9905818107052337` (output L31) vs wl `0.990581810705233684416290498...` (output L108) — agree.
- Sample S_+/S_-/B_both/B_nonempty rows agree across both (py output L21-24 vs wl output L59-94).

No `engine_disagreement`.

## Verdict justification

The math holds up. I attacked: (a) the F2 de-tautologization — `solve(Eq(onset,R_*),δ)` (py L119) genuinely re-derives δ_* from the onset relation `8/(9δ)=R_*` and R_*; it is not value-vs-itself because R_* is built from the carried slope rationals and δ_*'s correctness is independently re-asserted (`assert simplify(delta_solutions[0]-delta_dyn_star)==0`, py L121); the de-taut fix HOLDS. The `.wl` mirrors this with `Solve[rNDAtOnset==rStar,delta]` (wl L126). (b) Independence — the `.wl` proves classifier and slope monotonicity and the sign split via `Resolve[ForAll]` universal quantification, operations the `.py` does not perform; independent, not a port. (c) Threshold emission/reflection — R_*≈1.229255438463336 and δ_*≈0.723111617875019 are emitted by BOTH scripts (sympy output L17-18; mathematica output L105-106) and appear in the notes (lines 46,50,229,313). The sole defect is the card's stale "Mathematica audit: none yet" line, a `paper_misalignment` routed to the user gate. Verdict: `findings` (one paper_misalignment).

## Self-test notes

I checked: variable independence — `sp.diff(R_ND,xi)` and `sp.diff(S_±,R)` are over symbols the expressions actually depend on (xi, R), no identically-zero-derivative trap. Trivial-case: at R=0, `S_-(0)=s_-^den>0` (nonempty ceiling infinite, matches); at R=R_*, `S_-=0` exactly (output confirms `0`). Parity: no unbounded-domain integrals here. The single fix (F1) is paper-side prose only and prescribes no script edit, so no new paper_misalignment can be introduced by a Codex change.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation against the notes `.md` (the natural carrier; the `.tex` card is terse and the appendix is qualitative — per the augmentation guards a value living correctly in the `.md` is a MATCH). Output line refs are to the committed sympy/mathematica `.txt`.

| value | source (py/wl + output) | .tex/.md location | status |
|---|---|---|---|
| s_+^den = -0.301516097158113 | py L77 / wl L75 | notes:87 | MATCH |
| s_-^den = +0.411024574532864 | py L78 / wl L77 | notes:89 | MATCH |
| s_+^num = -0.508643465308977 | py L79 / wl L74 | notes:114 | MATCH |
| s_-^num = -0.334368725711457 | py L80 / wl L76 | notes:116 | MATCH |
| R_* = 1.229255438463336 | sympy out L17 / math out L105 | notes:46,229,539 | MATCH |
| delta_*^dyn = 0.723111617875019 | sympy out L18 / math out L106 | notes:50,313,325,549 | MATCH |
| ell_- = 0.323428979934714 | sympy out L27 / math out L47 | notes:351 | MATCH |
| ell_+ = 0.503852964869151 | sympy out L28 / math out L49 | notes:355 | MATCH |
| RQ_- = 30.199907560250075 | py L144 / wl L135 | notes:340 | MATCH |
| RQ_+ = 36.171186483269487 | py L145 / wl L136 | notes:342 | MATCH |
| RQ_req = 21.854566296358396 | py L146 / wl L137 | notes:346 | MATCH |
| B_both(0) = 1.671064893775584 | sympy out L29 / math out L51 | notes:400,485 | MATCH |
| inf B_both = 0.967282389363822 | sympy out L30 / math out L107 | notes:406,448 | MATCH |
| inf B_nonempty = 0.990581810705233 | sympy out L31 / math out L108 | notes:412,459 | MATCH |
| B_stat^both = 0.367930328492646 | sympy out L34 / math out L101 | notes:428,451 | MATCH |
| B_stat^nonempty = 0.737619063660757 | sympy out L35 / math out L103 | notes:433,462 | MATCH |
| sample R=1: S_+=-0.405079781233545, S_-=+0.0383279244107035, B_both=1.243836370541187 | py L176 / wl L185; outputs | notes:493,495,499 | MATCH |
| sample R=10: S_+=-0.489813704567990, S_-=-0.266605698416519, B_both=1.028662448947899, B_nonempty=1.213136035184892 | py L178 / wl L187; outputs | notes:514,516,520,522 | MATCH |

INTERNAL (genuine scaffolding / intermediate sample values, not stated deliverables, no finding):
- `S_+(R_*) ≈ -0.415730215182002` and `B_both(R_*) ≈ 1.211971000588856` — the R=R_* sample row (py L177 / wl L186). Notes §7.3 reports only `S_-(R_*)=0` (notes:507) as the deliverable for that point; the other two are intermediate sample-point figures. Terse-carrier guard → INTERNAL, not MISSING.
- Pass/fail flags, `assert_close`/`expectClose` tolerances (1e-12..5e-13), residual diffs in the `.txt`, `dR_ND/dxi` and `dS_±/dR` symbolic forms, `w_num`/`w_den` symbolic forms — scaffolding.

reconciliation: complete; 18 deliverable values checked, 0 misaligned. (The only paper_misalignment, F1, is a prose verification-coverage claim, not a value mismatch.)
