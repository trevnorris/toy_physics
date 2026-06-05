---
unit_id: 044
batch: III.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage044_continuum_selected_rank2_closure.md]
  paper_appendix: present
---

# Audit unit 044 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_044.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage044_continuum_selected_rank2_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 66; intro lines 9–36)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}{The continuum-selected quadratic \eqref{eq:app-stage044-quadratic} and normalization gate \eqref{eq:app-stage044-normalization-gate}.}` The stage inserts the Stage-043 physical support data (`R_phi`, `M_supp`) and source data (`R_U`, `M_mix`, `t² = λ0 = 2/9`, `δ`) into the Stage-041 support-loading theorem and shows the selected wall branch is pinned by one exact quadratic in the softening depth: `ξ² + B_cont ξ + C_cont = 0`, with `B_cont = δ − M_mix(1+t²R_U²) − M_supp(1+t²R_φ²)` (eq:Bcont) and `C_cont = −δ(M_mix+M_supp) + t²M_mix M_supp(R_U−R_φ)²` (eq:Ccont). The physical root is the `+` branch `ξ_phys = [−B_cont + √(B_cont²−4C_cont)]/2` (notes §3), chosen because it reduces to `ξ=0` at zero load. The normalization gate is `R_target = F_cont(ξ_phys)` (eq:normalization-gate), with `F_cont` given explicitly in notes §4. The notes add three more deliverables: the minimal-kernel source-tied surface (`R_φ=1` lands on the Stage-041/042 source-tied closure), the interference-matched tracking surface (`R_φ=R_U` collapses the branch to `M_mix+M_supp = G_q(ξ,δ)`), and the exact mismatch penalty `λ0 M_mix M_supp(R_U−R_φ)²` entering `C_cont`. Appendix row (line 66): "Physical softening depth fixed by a quadratic and selected normalization gate," `\StatusExactClosure{}`.

## What the script claims to verify

The SymPy/Mathematica scripts build `n_req^(cont)` from the Stage-041 theorem with `q = √λ0·R_U`, `r = √λ0·R_φ`, `λ0 = 2/9`, then (1) assert the numerator of `n_req − M_supp` equals `9·(ξ² + B_cont ξ + C_cont)` with the paper's `B_cont`/`C_cont` (the factor 9 is the cleared `2/9` denominator); (2) assert the `+`-branch `ξ_phys` vanishes at `M_mix=M_supp=0`; (3) build `F_cont` and check an independent literal `R_φ=2` slice; (4) check `n_req|_{R_φ=1}` and `F_cont|_{R_φ=1}` equal hand-written source-tied forms; (5) check `n_req|_{R_φ=R_U} = G_q − M_mix` (the tracking collapse) and `F_cont|_{R_φ=R_U}` equals the collapsed tracking form; (6) extract the `ξ⁰` coefficient of the branch numerator and confirm it equals `9·C_cont`, isolating the mismatch penalty. The Mathematica script additionally solves `branchEq==0` with `Solve`, selects the zero-load-vanishing root, and confirms it equals the explicit quadratic-formula `ξ_phys` — a genuinely independent derivation route the SymPy script does not have.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Quadratic `ξ²+B_cont ξ+C_cont=0`, eq Bcont/Ccont | `expect_zero("quadratic branch equation", branch_eq − 9·branch_expected)` (py:72 / wl:56) | match |
| Physical root = `+` branch (zero-load → 0) | `expect_zero("zero-load root", ξ_phys.subs{M_mix:0,M_supp:0})` (py:78 / wl:61); wl Solve cross-route (wl:65–77) | match |
| `R_target = F_cont(ξ_phys)`, `F_cont` form (notes §4) | `F_cont` built per notes §4; `third-slice F at R_φ=2` (py:109 / wl:103) constrains it | match (gate equation itself is structural; `F_cont` form verified, not evaluated at a numeric `ξ_phys`) |
| Source-tied surface `R_φ=1` (notes §5.1) | `source-tied n` / `source-tied F` (py:131–132 / wl:123–124) | match |
| Tracking surface `R_φ=R_U` → `M_mix+M_supp=G_q` (notes §5.2) | `tracking collapse of n_req` / `tracking F collapse` (py:138,146 / wl:137–138) | match |
| Mismatch penalty `λ0 M_mix M_supp(R_U−R_φ)²` in `C_cont` (notes §6) | `mismatch penalty in C coefficient` (py:154 / wl:148) | match |

`paper_alignment: aligned`. Every paper/notes deliverable maps to a non-tautological script check. `t² = λ0 = 2/9` (notes:55) matches the script literal `lambda0 = 2/9` (py:48 / wl:40); the `.tex` carries `t²` symbolically (lines 23,29), which the notes pin numerically — so the constant reconciles via the notes (MATCH).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 72 | `expect_zero(branch_eq − 9·branch_expected)` | quadratic (Output) | yes |
| A2 | sympy | 78 | `expect_zero(ξ_phys` at zero load`)` | `+`-branch selection (notes §3) | yes |
| A3 | sympy | 109 | `expect_zero(F_lit − F_lit_expected)` at `R_φ=2` | `F_cont` form (gate) | yes (independent literal slice) |
| A4 | sympy | 131 | `expect_zero(n_source − n_source_expected)` | source-tied surface (§5.1) | partial (R_φ→1 substitution vs hand-copy) |
| A5 | sympy | 132 | `expect_zero(F_source − F_source_expected)` | source-tied surface (§5.1) | partial (same) |
| A6 | sympy | 138 | `expect_zero(n_track − (G_q − M_mix))` | tracking collapse (§5.2) | yes |
| A7 | sympy | 146 | `expect_zero(F_track − F_track_expected)` | tracking surface (§5.2) | yes |
| A8 | sympy | 154 | `expect_zero(C_from_branch_eq − C_expected)` | mismatch penalty (§6) | yes (redundant w/ A1 but isolates `C`) |
| B1 | mathematica | 56 | `expectZero[branchEq − 9·branchExpected]` | quadratic (Output) | yes |
| B2 | mathematica | 61 | `expectZero[xiPhys` at zero load`]` | `+`-branch (§3) | yes |
| B3 | mathematica | 77 | `expectZero[xiPhysSolve − xiPhys]` | `+`-branch via independent `Solve` | yes (strong, independent route) |
| B4 | mathematica | 103 | `expectZero[fLit − fLitExpected]` at `rPhi=2` | `F_cont` form | yes |
| B5 | mathematica | 123–124 | `expectZero[nSource…]`, `expectZero[fSource…]` | source-tied (§5.1) | partial |
| B6 | mathematica | 137–138 | `expectZero[nTrack…]`, `expectZero[fTrack…]` | tracking (§5.2) | yes |
| B7 | mathematica | 148 | `expectZero[cFromBranchEq − cExpected]` | mismatch penalty (§6) | yes |

A4/A5/B5 are "partial" because the `_expected` forms are hand-transcriptions of the same master expression evaluated at `R_φ=1`; the check confirms the source-tied surface is exactly the `R_φ→1` limit (the actual §5.1 deliverable) but cannot independently fail if the master `n_req`/`F_cont` were themselves wrong (those are independently anchored by A1/A3/B1/B4). This is a consistency check, not load-bearing-tautological — no finding.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.txt` (mtime 2026-05-22 12:45) vs script (mtime 2026-06-03 15:59)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.txt` (mtime 2026-05-22 12:45) vs script (mtime 2026-06-03 15:59)
- Internal stale labels in the SymPy source: docstring `scripts/...sympy_audit.py:3` ("Stage 27 SymPy audit"); subbanners `:53,:80,:111,:134,:148` ("27.1"–"27.5"); ledger banner `:156` ("STAGE 44 THEOREM LEDGER" — header was renumbered but the §x subbanners and docstring were not).

**What's wrong:**
The saved `.txt` outputs predate the scripts and carry pre-renumber banners: the SymPy transcript opens with `STAGE 27 — CONTINUUM-SELECTED RANK-2 CLOSURE` and `STAGE 27 THEOREM LEDGER`; the Mathematica transcript opens with `STAGE 027 …` (yet closes with the post-renumber line `Stage 044 Mathematica audit passed.`, confirming a partial-renumber epoch). The git history shows the scripts were last touched only by the doc-only numbering-reconciliation commit `e2a4780`; the math content is unchanged from the run that produced these `.txt` files. So the staleness is label-only: every emitted symbolic form and every `= 0 / PASS` line in the saved outputs matches what the current scripts compute. The SymPy source itself still carries the old `27.x` subbanner labels and a "Stage 27" docstring even though its top banner was updated to `STAGE 44`.

**Why this matters:**
Nothing mathematical breaks — the residuals, `F_cont`, `n_req`, and `ξ_phys` forms in the saved outputs are correct and current. The only defect is cosmetic banner/label drift in the script/output band.

**Required change:**
None applied here. This is the known, user-deferred SCRIPT/OUTPUT-band numbering cleanup (memory: `1fa4f7a`, plan `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`, PENDING, content-keyed, handled in a dedicated separate pass — never offset-swept). The orchestrator's independent exec re-run will regenerate fresh `.txt` outputs; no Codex directive is emitted for the label drift to avoid colliding with the gated dedicated pass.

**Verification:**
After the orchestrator re-runs `exec-sympy 044` / `exec-mathematica 044`, the regenerated `scripts/output/…txt` and `mathematica/output/…txt` should show `STAGE 44`/`STAGE 044` headers, all the same `= 0`/`PASS` lines, and a newer mtime. No symbolic value should change.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration. It mirrors the SymPy choreography for the algebraic identities (unavoidable — the same `n_req`, `F_cont`, `B_cont`, `C_cont` definitions), but it adds a genuinely independent derivation route the SymPy script lacks: lines 65–77 call `Solve[branchEq == 0, xi]`, `SelectFirst` the root that vanishes at zero load, and `expectZero[xiPhysSolve − xiPhys]`. That confirms the explicit quadratic-formula `ξ_phys` agrees with Mathematica's own algebraic solver — a different engine path to the same root, not an echo of the SymPy `(-B+√Δ)/2` literal. Example corresponding sections: SymPy `xi_phys = (-B_cont + sqrt(Delta_disc))/2` (py:75) is matched in the `.wl` both by the same closed form (wl:59) AND independently by `Solve` (wl:65–77). This satisfies the second-engine policy.

## Engine cross-check

Both engines emit the same `= 0` / `PASS` for every shared check. Saved outputs agree:
- SymPy: `quadratic branch equation = 0`, `zero-load root = 0`, `third-slice F at Rphi=2 = 0`, `source-tied n/F = 0`, `tracking … = 0`, `mismatch penalty in C coefficient = 0`.
- Mathematica: `PASS: quadratic branch equation`, `PASS: zero-load root`, `PASS: xiPhys solve match`, `PASS: third-slice F at rPhi=2`, `PASS: source-tied n/F`, `PASS: tracking …`, `PASS: mismatch penalty in C coefficient`, then `Stage 044 Mathematica audit passed.`
The printed symbolic forms agree across engines (e.g. SymPy `n_req^(cont)` numerator `2·M_mix·R_U²·ξ + 9·M_mix·δ + 9·M_mix·ξ − 9·δ·ξ − 9·ξ²` matches the Mathematica `(xi(delta+xi) − mMix(delta+xi+2rU²xi/9))` after clearing the `/9`). `engines_agree: true`.

## Verdict justification

The scripts faithfully and non-tautologically verify every deliverable the paper card and notes state: the exact quadratic with the paper's `B_cont`/`C_cont`, the `+`-branch root selection (the `-` branch would give `−δ ≠ 0` at zero load, so the check genuinely discriminates), the explicit `F_cont` form (constrained off the two collapsed surfaces by an independent `R_φ=2` literal slice), the source-tied (`R_φ=1`) and tracking (`R_φ=R_U`) special surfaces, and the mismatch penalty. Attacks tried and failed: (a) the factor-of-9 in `branch_eq − 9·branch_expected` is the cleared `λ0=2/9` denominator, not a fudge — `B_cont`/`C_cont` are independently defined and the equality is non-trivial; (b) the zero-load root check is not tautological — it selects the physical branch; (c) the assumption set (`δ>0`, `M_mix,M_supp≥0`, `R_U,R_φ>0`, all real) matches the physical setup and is not strong enough to make arbitrary candidates pass; (d) `R_φ=2` literal slice rules out the worry that `F_cont` was only checked on the two collapsed surfaces. I read the `.tex` card, the notes, and the appendix row; the script's verified claim matches the paper's `Output` exactly. The sole finding is label-only `stale_output` (the known, user-deferred script/output-band numbering drift), which is informational and carries no math defect — hence `verdict: findings` with `material_change` effectively false and no directive.

## Self-test notes

Checked the three standard traps. (1) Variable independence: no `diff`/`D` derivative checks in this unit — the assertions are algebraic identity residuals (`expect_zero`/`expectZero` of `lhs − rhs`), so the zero-derivative trap does not arise. (2) Symmetry/parity: no unbounded integrals; n/a. (3) Trivial-case pre-check: substituted `M_mix=M_supp=0` into `ξ_phys` and confirmed the `+` branch gives exactly `0` (B_cont=δ, C_cont=0, √(δ²)=δ, so `(−δ+δ)/2=0`) while the `−` branch would give `−δ≠0` — the zero-load assertion is a real discriminator, not a pass-by-construction. Also confirmed the `branch_eq = 9·branch_expected` factor-9 is the cleared `2/9` denominator and `t²=λ0=2/9` matches notes:55 — no value_mismatch.

## Value Reconciliation (pass-2 augmentation)

All deliverables for this stage are **symbolic** (closed-form quadratic and normalization-function structures); the only numeric constant is `λ0 = t² = 2/9`. The `.txt` outputs are label-stale (F1) but their math content is current, so the reconciliation is based on the script source plus the saved symbolic forms, which agree.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `λ0 = t² = 2/9` | py:48 `lambda0 = Rational(2,9)` / wl:40 `lambda0 = 2/9` | notes:55 `t^2 = lambda_0 = 2/9`; .tex:23,29 carry `t²` symbolically | MATCH |
| `B_cont = δ − M_mix(1+λ0 R_U²) − M_supp(1+λ0 R_φ²)` | py:69 / wl:50 | .tex:21–24 (eq Bcont); notes:105–108 | MATCH |
| `C_cont = −δ(M_mix+M_supp) + λ0 M_mix M_supp(R_U−R_φ)²` | py:70 / wl:51 | .tex:25–30 (eq Ccont); notes:110–112 | MATCH |
| Quadratic `ξ²+B_cont ξ+C_cont=0` (numerator = `9·(…)`) | py:72 out:21 `= 0` / wl:56 out:11–12 PASS | .tex:14–18 (eq quadratic, boxed); notes:99–101 | MATCH |
| `ξ_phys = (−B_cont + √(B_cont²−4C_cont))/2` (+ branch) | py:75,77 out:22–52 / wl:59,76 out:13,16 | notes:116–117 (§3) | MATCH |
| `F_cont(ξ)` closed form (normalization gate `R_target=F_cont(ξ_phys)`) | py:87–91 out:57–88 / wl:85–90 out:23 | .tex:31–36 (eq normalization-gate, boxed `F_cont`); notes:147–157 (§4) | MATCH |
| Source-tied surface (`R_φ=1`) `n_source`,`F_source` | py:114–119 out:95–124 / wl:107–108 out:30–31 | notes:172–184 (§5.1) | MATCH |
| Tracking collapse (`R_φ=R_U`) `M_mix+M_supp=G_q`, `F_track` | py:136–146 out:131–132 / wl:128–138 out:40–43 | notes:186–214 (§5.2) | MATCH |
| Mismatch penalty `λ0 M_mix M_supp(R_U−R_φ)²` in `C_cont` | py:152–154 out:137 / wl:143–148 out:48 | notes:228,256 (§6/§7) | MATCH |

INTERNAL scaffolding (no finding, not expected in prose): `n_req^(cont)` intermediate rational form; `branch_eq` (cleared numerator); `Delta_disc` (`B²−4C`); `xiPhysSolve` (Mathematica `Solve` cross-route); `D_cont`/`dCont`; the `R_φ=2` literal slice `F_lit`/`fLit` and their `_expected` twins; `G_q`; all `*_expected` hand-transcribed comparison forms; `pass/fail` flags and `=0`/`PASS` residual lines.

reconciliation: complete; 9 deliverable values checked, 0 misaligned.
