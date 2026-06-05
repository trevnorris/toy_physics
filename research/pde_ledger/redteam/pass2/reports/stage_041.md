---
unit_id: 041
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
  notes_stage_files: [moving_throat_pde_stage041_rank2_support_completion.md]
  paper_appendix: present
---

# Audit unit 041 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_041.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage041_rank2_support_completion.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 60 + `\input` at line 200)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage041_rank2_support_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage041_rank2_support_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage041_rank2_support_mathematica_audit.txt`

## What the paper claims

Stage 041 (Part III, anchor MTDC-T6, claim status: exact closure) adds a second loading direction (the "support" loading `n` along `(1,r)`) on top of the mixed baseline loading `m` along `(1,q)`, on a diagonal baseline `diag(1,1+delta)`. The bottom-line `\stagefield{Output}` is verbatim: *"The support-loading theorem \eqref{eq:app-stage041-nreq}, its monotonicity \eqref{eq:app-stage041-nreq-derivative}, and the tracking/source-tied alternatives."* The distinct deliverables are: (1) the selected-branch determinant `D_sel` (eq. app-stage041-Dsel) is **linear in n**; (2) the exact support-loading theorem `n_req(xi,delta;m,q,r)` (boxed, eq. nreq); (3) the exact monotonicity law `∂n_req/∂m = -[δ+(1+qr)ξ]²/[δ+(1+r²)ξ-m(q-r)²]² < 0` (boxed, eq. nreq-derivative); (4) the tracking closure `r=q ⇒ n_req = G_q(ξ,δ) - m` (eq. tracking); and (5) the source-tied closure `n_req^(src)` with `q=tR_U`, `r=t`, `t²=λ₀=2/9` (eq. source-tied). The notes add, in §4, two feasibility thresholds (regularity pole `m < [δ+(1+λ₀)ξ]/[λ₀(R_U-1)²]` and positivity ceiling `m ≤ ξ(δ+ξ)/[δ+(1+λ₀R_U²)ξ]`) and a source-tied derivative.

## What the script claims to verify

Both engines build the reduced dimensionless loaded matrix `M = diag(1,1+δ) - m(1,q)ᵀ(1,q) - n(1,r)ᵀ(1,r)`, form `det(M-(1-ξ)I)`, and assert it equals the paper's `D_sel` decomposition (proving linearity in `n`). They then `solve` `D_sel=0` for `n` and assert the closed form equals `n_expected` (the boxed `n_req`); `diff` `n_expected` w.r.t. `m` and assert it equals the boxed monotonicity expression; substitute `r→q` and assert collapse to `G_q - m`; and substitute `q→t·R_U, r→t, t²→λ₀=2/9` into the *general* `n_expected` and assert the result equals an independently-written `n_src_expected` (the boxed source-tied formula). The scripts additionally print the two feasibility thresholds and assert the source-tied derivative `-(δ+(1+λ₀R_U)ξ)²/(…)²`. All eight `expect_zero`/`expectZero` checks must hit literal `0` or the run aborts.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1: `D_sel` linear in `n` (eq Dsel) | `expect_zero("determinant decomposition", Det - Det_expected)` (py:72, wl:56) | match |
| D2: boxed `n_req` (eq nreq) | `expect_zero("n_req - expected", n_req - n_expected)` (py:82, wl:65) | match |
| D3: boxed `∂n_req/∂m < 0` (eq nreq-derivative) | `expect_zero("dn/dm - expected", dn_dm - monotone_expected)` (py:94, wl:76) | match |
| D4: tracking `n_req=G_q-m` (eq tracking) | `expect_zero("tracking collapse", n_track - (G_q - m))` (py:102, wl:84) | match |
| D5: source-tied `n_req^(src)` (eq source-tied) | `expect_zero("source-tied formula", n_src - n_src_expected)` (py:121, wl:113) | match |
| Notes §4 regularity + positivity thresholds | printed `reg_threshold`, `num_zero_threshold` (py:124–129, wl:104–105) | match (notes-carried; printed not asserted, but they are display-only feasibility windows, not closure identities) |
| Notes §4 source-tied derivative | `expect_zero("source-tied dn/dm", …)` (py:140, wl:117) | match (extra over .tex, present in notes) |

`paper_alignment: aligned` — every boxed/Output deliverable has a faithful, non-tautological script check, and the only "extra" checks (source-tied derivative, thresholds) are anchored in the notes.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 72 | `simplify(Det - Det_expected)==0` | D1 (linearity in n) | yes |
| A2 | sympy | 82 | `simplify(n_req - n_expected)==0` | D2 (n_req) | yes |
| A3 | sympy | 94 | `simplify(dn_dm - monotone_expected)==0` | D3 (monotonicity) | yes |
| A4 | sympy | 102 | `simplify(n_track - (G_q - m))==0` | D4 (tracking) | yes |
| A5 | sympy | 121 | `simplify(n_src - n_src_expected)==0` | D5 (source-tied) | yes |
| A6 | sympy | 140 | `simplify(dn_dm_src - dn_dm_src_expected)==0` | notes §4 src derivative | yes |
| B1 | math | 56 | `FullSimplify[detExpr-detExpected]===0` | D1 | yes |
| B2 | math | 65 | `FullSimplify[nReq-nExpected]===0` | D2 | yes |
| B3 | math | 76 | `FullSimplify[dndm-monotoneExpected]===0` | D3 | yes |
| B4 | math | 84 | `FullSimplify[nTrack-(gQ-m)]===0` | D4 | yes |
| B5 | math | 113 | `FullSimplify[nSrc-nSrcExpected]===0` | D5 | yes |
| B6 | math | 117 | `FullSimplify[dndmSrc-dndmSrcExpected]===0` | notes §4 src derivative | yes |

Every assertion traces to a paper or notes deliverable. None is tautological: each LHS is produced by an independent CAS operation (`det`, `solve`, `diff`, substitution) while each RHS is the paper's independently-transcribed closed form, so the subtraction can genuinely fail if the algebra or the paper transcription disagree.

## Findings

### F1 — stale_output

**Severity:** low (informational; non-blocking — belongs to the user-gated deferred SCRIPT/OUTPUT-band numbering pass)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage041_rank2_support_sympy_audit.txt:3,7,...,77` (mtime 2026-05-22 12:34:49)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage041_rank2_support_mathematica_audit.txt:3` (mtime 2026-05-22 12:34:58)
- script mtimes for comparison: `…sympy_audit.py` 2026-06-03 15:59:11; `…mathematica_audit.wl` 2026-06-03 15:59:11

**What's wrong:**
Both saved outputs predate their scripts by ~12 days, and the SymPy output's content disagrees with the current script's banners. The committed SymPy `.txt` prints `STAGE 24 — RANK-2 SUPPORT COMPLETION` (txt:3) and `STAGE 24 THEOREM LEDGER` (txt:77), whereas the current script emits `STAGE 41` (py:51) and `STAGE 41 THEOREM LEDGER` (py:142). The Mathematica `.txt` prints `STAGE 024` (txt:3) while the current `.wl` banner is `STAGE 041` (wl:33) — and the same `.txt` later prints `Stage 041 Mathematica audit passed.` (txt:42), an internal inconsistency confirming it is from an older revision. Separately, the **script source itself** carries stale self-labels that are the known numbering-drift class: SymPy docstring `Stage 24 SymPy audit` (py:3), subbanners `24.1`–`24.4` (py:53,84,96,104) while the banner already reads `STAGE 41`, an inline comment `section 24.1` (py:107), and cross-refs to `Stage-23 one-direction geometry` (py:13,151). These are label-only; none touches the math.

**Why this matters:**
The captured transcripts no longer reflect the script's current banner text, and the script has internal banner/subbanner number disagreement (`STAGE 41` vs `24.1`). This is purely cosmetic — every residual line in both transcripts is `= 0` / `PASS`, the symbolic results match the current scripts, and no deliverable value is affected. The risk is only reviewer confusion, not a math defect.

**Required change:**
None applied by this audit. Per the user decision of 2026-06-04, all SCRIPT/OUTPUT-band stage-number relabeling (self-banner, docstring, subbanner, ledger, and stale committed `.txt` banners) is deferred to the dedicated content-keyed pass (`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`, PENDING) — it must NOT be offset-swept or hand-patched here. The orchestrator's independent exec re-run will regenerate fresh transcripts (which will then read `STAGE 41`/`STAGE 041`), resolving the mtime/banner staleness for the math content. No directive is written so as not to pre-empt the user-gated numbering pass. (Cross-ref `Stage-23` → the current one-direction geometry is Stage 040; that too is for the deferred pass, content-keyed.)

**Verification:**
After a fresh exec run, `scripts/output/…sympy_audit.txt` line 3 reads `STAGE 41 …` (matching py:51) and all `= 0`/`PASS` lines persist; `mathematica/output/…mathematica_audit.txt` line 3 reads `STAGE 041 …` (matching wl:33). The residual content is unchanged.

## Independent-derivation check (Mathematica)

The `.wl` follows the same physical premise → algebra choreography as the `.py` (same matrix `mMat`, same `detExpected`/`nExpected` targets, same `r→q` and `q→t rU, r→t, t²→λ₀` substitution sequence). This is *not* a transliteration in the disqualifying sense: both engines start from the **physical loaded matrix** (not from each other's intermediate algebra) and each independently runs its own CAS kernel for `Det`/`Solve`/`D`/`FullSimplify`. The "expected" RHS forms each engine compares against are the paper's boxed closed forms, transcribed independently into each language, not echoed intermediates. For an algebra-identity verification stage, this dual-CAS confirmation IS the legitimate second-engine value. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines reach identical results. Side by side, displayed forms differ only by overall-sign-cancelling factoring and by Mathematica vs SymPy pretty-print conventions:

- `det` decomposition: py txt `-δm-δn+δξ+mnq²-2mnqr+mnr²-mq²ξ-mξ-nr²ξ-nξ+ξ²`; wl txt `-(delta*m)-delta*n+m*n*q^2-2*m*n*q*r+m*n*r^2+delta*xi-m*xi-n*xi-m*q^2*xi-n*r^2*xi+xi^2`. Identical. `determinant decomposition = 0` in both.
- `n_req`, `dn/dm`, `n_req(r=q)`, `n_req^(src)`, thresholds, `dn/dm^(src)`: every corresponding `expect_zero`/`expectZero` residual is `0` / `PASS` in both transcripts. The source-tied numerator/denominator appear with a common `-`/`-` factor between the two pretty-prints that cancels in the ratio (manually verified: numerator and denominator are each negated, leaving the fraction identical to the paper's eq. source-tied with λ₀=2/9). `engines_agree: true`.

## Verdict justification

`verdict: findings` with a single low/informational `stale_output` (F1). Adversarial attacks that I tried and that **failed to break the stage**: (a) I checked whether the source-tied derivative numerator should be `(1+λ₀R_U²)` instead of `(1+λ₀R_U)` — it is correctly `(1+λ₀R_U)`, since the general `(1+qr)` term gives `qr = (tR_U)(t) = t²R_U = λ₀R_U`, so the script is right; (b) I checked whether the `.subs(t**2, lam0)` could leave odd powers of `t` — all `t` enter in even powers (q², and the factored `t²(R_U-1)²` denominator), and both transcripts show no residual `t`, so the substitution fully clears; (c) I checked the `solve(D_sel=0, n)[0]` for branch ambiguity — `D_sel` is genuinely linear in `n`, so the solution is unique and `[0]` is safe; (d) I checked the matrix construction against `diag(1,1+δ) - m(1,q)ᵀ(1,q) - n(1,r)ᵀ(1,r)` and it matches exactly; (e) I verified the source-tied closed form equals the paper's eq. source-tied under λ₀=2/9 despite the sign-cancelling factoring in the printout. I read the paper card, the notes, and the appendix row (part03 line 60); every boxed/Output deliverable maps to a non-tautological, well-anchored, engine-agreeing assertion. The math is clean. The only defect is the cosmetic stale banners/labels, which fall squarely in the user-gated deferred SCRIPT/OUTPUT-band numbering class and must not be hand-patched here. `material_change: false`.

## Self-test notes

I checked the trap list: (1) variable-independence — every `diff` is `d/dm` of `n_expected`/`n_src_expected`, both of which genuinely depend on `m` (m appears in numerator and denominator), so the derivative is non-trivial and `dn/dm ≠ 0` literal as printed; (2) symmetry/parity — n/a, no unbounded-domain integrals in this stage; (3) trivial-case — substituting the source-tied invariants reduces to the printed `n_req^(src)`, confirmed `=0` against the independent RHS; (4) path specs — n/a (no missing-script directive); (5) paper round-trip — no script fix prescribed, so no risk of introducing a new paper_misalignment; the only finding is deferred-by-policy. Conclusion: the algebra holds against the paper on all five deliverables; the sole finding is cosmetic stale-output/numbering, deferred to the dedicated user-gated pass.

## Value Reconciliation (pass-2 augmentation)

Every emitted deliverable value is symbolic (this stage pins no free-standing numeric constant other than the inherited `λ₀ = 2/9 = t²`, which is an upstream-defined input, not a Stage-041 deliverable). The reconciliation is therefore symbolic closed-form ↔ paper/notes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `D_sel` decomposition `ξ(δ+ξ)-m[δ+(1+q²)ξ]-n[δ+(1+r²)ξ]+mn(q-r)²` | py:66–72 / wl:51–56; sympy txt:9–15, math txt:9 | `stage_041.tex:22-29` (eq Dsel); notes §1 lines 64–68 | MATCH |
| `n_req = [ξ(δ+ξ)-m(δ+(1+q²)ξ)]/[δ+(1+r²)ξ-m(q-r)²]` | py:78–82 / wl:59–65; sympy txt:16–22, math txt:12 | `stage_041.tex:31-37` (boxed eq nreq); notes §2 lines 82–84 | MATCH |
| `∂n_req/∂m = -[δ+(1+qr)ξ]²/[δ+(1+r²)ξ-m(q-r)²]²` | py:90–94 / wl:70–76; sympy txt:27–34, math txt:19 | `stage_041.tex:39-46` (boxed eq nreq-derivative); notes §2 lines 96–99 | MATCH |
| `n_req(r=q) = G_q - m`, `G_q = ξ(δ+ξ)/[δ+(1+q²)ξ]` | py:98–102 / wl:80–84; sympy txt:39–45, math txt:26 | `stage_041.tex:50-55` (eq tracking); notes §3 lines 119–120 | MATCH |
| `n_req^(src) = [ξ(δ+ξ)-m(δ+(1+λ₀R_U²)ξ)]/[δ+(1+λ₀)ξ-mλ₀(R_U-1)²]` | py:117–121 / wl:98–113; sympy txt:50–56, math txt:33 | `stage_041.tex:57-62` (eq source-tied); notes §4 lines 163–165 | MATCH |
| regularity pole threshold `[δ+(1+λ₀)ξ]/[λ₀(R_U-1)²]` | py:124,126–127 / wl:104,114; sympy txt:57–61, math txt:36 | notes §4 lines 178–179 | MATCH (notes-carried) |
| positivity ceiling `ξ(δ+ξ)/[δ+(1+λ₀R_U²)ξ]` | py:125,128–129 / wl:105,115; sympy txt:62–66, math txt:37 | notes §4 lines 186–187 | MATCH (notes-carried) |
| `∂n_req^(src)/∂m = -[δ+(1+λ₀R_U)ξ]²/[δ+(1+λ₀)ξ-mλ₀(R_U-1)²]²` | py:132,136–140 / wl:106–117; sympy txt:67–73, math txt:38 | notes §2 monotonicity generalizes; notes §4 (regular branch) lines 100–103 | MATCH (notes-carried; specializes the general dn/dm<0 to the source-tied case) |

INTERNAL (scaffolding, no finding): `M`/`mMat` (reduced loaded matrix), `lam = 1-xi` (reduced eigenvalue), `Det_expected`/`detExpected`, `n_expected`/`nExpected`, `n_src_expected`/`nSrcExpected`, the eight `expect_zero`/`expectZero` residual `0`s and `PASS`/ ledger prose, and the symbol `t` (auxiliary source-ratio). `λ₀ = 2/9` is an inherited input (Stage-039 source ratio), correctly used; it appears in notes §4 line 152 (`t² = λ₀ = 2/9`) — MATCH as an input.

reconciliation: complete; 8 deliverable values checked, 0 misaligned.
