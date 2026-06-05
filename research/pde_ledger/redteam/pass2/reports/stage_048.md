---
unit_id: 048
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
  notes_stage_files: [moving_throat_pde_stage048_support_compensation_theorem.md]
  paper_appendix: present
---

# Audit unit 048 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_048.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage048_support_compensation_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 74 + preview eqs)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage048_support_compensation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.txt`

## What the paper claims

Stage 048 proves there is **no reduced-level support no-go** on the coherent tracking branch. `\stagefield{Output}` (stage_048.tex:31): *"The exact threshold \eqref{eq:app-stage048-zeta-req} and success condition \eqref{eq:app-stage048-success-iff}."* The two boxed results are the exact inverse of the support factor, `\zeta_{\rm req}=\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)}` (eq. zeta-req, line 22), and the success condition `\zeta_{\rm phys}\geq\zeta_{\rm req}` (eq. success-iff, line 28). The card defines `S_{\rm req}=M_{\rm req}/M_{\rm mix}` and rests on the claim that `S` is strictly increasing up to its softening pole so `\zeta_{\rm req}` is unique and below the pole. The notes flesh this out into six exact supporting identities: the tracking law `G_tr` and its critical load `M_crit=9(1+\delta)/(9\delta+9+2R^2)`, the monotone derivative `dG_tr/d\xi`, the gap `M_crit-G_tr`, the support factor `S(\zeta;\epsilon)=1+\zeta(1-\epsilon)/(1-\zeta\epsilon)` with derivative `dS/d\zeta=(1-\epsilon)/(1-\zeta\epsilon)^2`, the pole margin `1/\epsilon-\zeta_{\rm req}`, the branch margin `\zeta_{\rm crit}-\zeta_{\rm req}`, and the monotone softening response `d\xi_{\rm phys}/d\zeta>0`. Appendix row 74 summarizes: *"Unique \(\zeta_{\rm req}\) threshold and no reduced-level support no-go."* (The appendix's `\zeta_{\rm req}=1/3` numerics at lines 32/341/359 are the downstream Family-1 precursor evaluation, not this stage's symbolic deliverable.)

## What the script claims to verify

Both engines verify the same chain of exact symbolic identities supporting the theorem. They confirm: (1) the closed form of `G_tr`, its softening-limit critical load `M_crit`, the factored monotone derivative `dG_tr/d\xi`, the `F_tr(0)=1` and `F_tr→+∞` endpoints with the exact `(1-\xi)F_tr` softening coefficient, and the factored gap `M_crit-G_tr`; (2) `S(0)=1`, the derivative `dS/d\zeta`, divergence of `S` at the pole, and — the load-bearing check — that substituting `\zeta_{\rm req}=(S_{\rm req}-1)/(1+\epsilon(S_{\rm req}-2))` back into `S` returns `S_{\rm req}` (genuine inverse-map verification, similarly for `\zeta_{\rm crit}`); (3) the factored pole margin `1/\epsilon-\zeta_{\rm req}` and branch margin `\zeta_{\rm crit}-\zeta_{\rm req}` against their boxed closed forms; and (4) the implicit derivative `d\xi_{\rm phys}/d\zeta=M_{\rm mix}(dS/d\zeta)/(dG_tr/d\xi)` against its closed form. The Mathematica script independently *solves* `S==S_{\rm req}` for `\zeta` (rather than hardcoding the formula) and checks the solution count is 1.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `\zeta_{\rm req}=(S_{\rm req}-1)/(1+\epsilon(S_{\rm req}-2))` (Output, boxed eq. zeta-req) | sympy L92/L95 `S(zeta_req)-Sreq==0`; math L86–93 Solve+verify | match |
| Success condition `\zeta_{\rm phys}\ge\zeta_{\rm req}` (Output, boxed eq. success-iff) | sympy L107–114 / math L96–104 pole & branch margins `>0` (form) + monotone `S` | match (logical condition follows from verified strict-monotone/margin identities) |
| `M_crit`, `dG_tr/d\xi`, gap (notes §1) | sympy L52–74 / math L49–72 | match |
| `S`, `dS/d\zeta`, pole divergence (notes §2) | sympy L81–90 / math L76–84 | match |
| pole margin, branch margin (notes §2,§4) | sympy L107–114 / math L100–104 | match |
| `d\xi_{\rm phys}/d\zeta` (notes §5) | sympy L118–125 / math L106–112 | match |

`paper_alignment: aligned`. Every paper-side symbolic deliverable has a faithful, non-tautological script-side counterpart in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 52 | `dG_tr/dxi - <factored form> == 0` | notes §1 monotone derivative | yes |
| A2 | sympy | 57 | `F_tr(0)-1 == 0` | notes §3 endpoint | yes |
| A3 | sympy | 61 | `(1-xi)F_tr coeff - <form> == 0` | notes §3 softening endpoint | yes |
| A4 | sympy | 69 | `M_crit-G_tr - <factored form> == 0` | notes §1 gap | yes |
| A5 | sympy | 81 | `S(0)-1 == 0` | notes §2 endpoint | yes |
| A6 | sympy | 82 | `dS/dzeta - (1-eps)/(1-zeta eps)^2 == 0` | notes §2 derivative | yes |
| A7 | sympy | 89 | `limit_phys != oo → raise` (under eps=1/(1+nu)) | notes §2 pole divergence | yes |
| A8 | sympy | 95 | `S(zeta_req)-Sreq == 0` | **Output: zeta_req inverse** | yes |
| A9 | sympy | 96 | `S(zeta_crit)-Scrit == 0` | notes §4 zeta_crit | yes |
| A10 | sympy | 107 | `pole margin - <form> == 0` | notes §2 margin | yes |
| A11 | sympy | 111 | `branch margin - <form> == 0` | notes §4 margin / success cond | yes |
| A12 | sympy | 120 | `dxi/dzeta - <form> == 0` | notes §5 monotone response | yes |
| B1 | mathematica | 49 | `dG_tr/dxi - <form> == 0` | notes §1 | yes |
| B2 | mathematica | 53 | `F_tr(0)-1 == 0` | notes §3 | yes |
| B3 | mathematica | 64 | `(1-xi)F_tr coeff - <form> == 0` | notes §3 | yes |
| B4 | mathematica | 68 | `M_crit-G_tr - <form> == 0` | notes §1 | yes |
| B5 | mathematica | 76 | `S(0)-1 == 0` | notes §2 | yes |
| B6 | mathematica | 77 | `dS/dzeta - <form> == 0` | notes §2 | yes |
| B7 | mathematica | 84 | `pole coeff - (1-eps)/eps^2 == 0` | notes §2 pole divergence | yes |
| B8 | mathematica | 87 | `Solve count != 1 → fail` | **Output: uniqueness of zeta_req** | yes |
| B9 | mathematica | 93 | `S(zeta_req)-Sreq == 0` | **Output: zeta_req inverse** | yes |
| B10 | mathematica | 94 | `S(zeta_crit)-Scrit == 0` | notes §4 | yes |
| B11 | mathematica | 100 | `pole margin - <form> == 0` | notes §2 | yes |
| B12 | mathematica | 101 | `branch margin - <form> == 0` | notes §4 / success cond | yes |
| B13 | mathematica | 108 | `dxi/dzeta - <form> == 0` | notes §5 | yes |

All assertions are anchored. No tautologies: each RHS in an `expect_zero`/`expectZero` is an independently *typed-in* closed form, and `S` is constructed independently from `zeta_req`/`zeta_crit`, so the inverse-map checks (A8/A9/B9/B10) and Solve uniqueness (B8) genuinely exercise the boxed `\zeta_{\rm req}`.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage048_support_compensation_sympy_audit.txt:3` (banner "STAGE 31")
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.txt:4` (banner "STAGE 031")

**What's wrong:**
Both saved outputs are older than their scripts:
- sympy `.py` mtime `2026-06-03 15:59:11`; sympy `.txt` mtime `2026-05-26 02:04:51`.
- mathematica `.wl` mtime `2026-06-03 15:59:11`; mathematica `.txt` mtime `2026-05-26 02:04:56`.

The content also shows the drift: the current `.py` banner (line 31) prints `STAGE 48 — SUPPORT COMPENSATION THEOREM AUDIT`, but the committed sympy `.txt` (line 3) still reads `STAGE 31 — SUPPORT COMPENSATION THEOREM AUDIT`. The current `.wl` banner (line 26) prints `STAGE 048 — SUPPORT COMPENSATION THEOREM`, but the committed mathematica `.txt` (line 4) still reads `STAGE 031 — SUPPORT COMPENSATION THEOREM`. So the captured transcripts predate the most recent banner/label edit of the scripts. The mathematical content of the outputs (all residuals = 0, all PASS) is consistent with the current scripts — only the banner labels are stale.

**Why this matters:**
The transcripts no longer byte-match what the current scripts would emit; a reviewer comparing the committed `.txt` banner to the `.py`/`.wl` banner sees a `STAGE 31`/`STAGE 031` vs `STAGE 48`/`STAGE 048` mismatch. The actual checks are unaffected (the verifier's independent re-run will refresh these), so this is informational, not blocking.

**Required change:**
Re-run both scripts and re-commit the refreshed transcripts so the saved banners read `STAGE 48`/`STAGE 048`. No script logic change is required for the outputs; see the separate label note below for the in-source docstring label.

**Verification:**
After a fresh run, sympy `.txt` line 3 reads `STAGE 48 — SUPPORT COMPENSATION THEOREM AUDIT` and mathematica `.txt` line 4 reads `STAGE 048 — SUPPORT COMPENSATION THEOREM`; both `.txt` mtimes are newer than their scripts.

#### Adjacent (non-finding) stale in-source labels — flagged for the dedicated numbering pass, NOT fixed here

Per the project's deferred SCRIPT/OUTPUT-band numbering plan (content-keyed, never offset-swept), these in-source self-labels are recorded but intentionally NOT corrected in this audit:
- sympy docstring `/var/.../scripts/...stage048...py:3`: `Stage 31 SymPy audit.` — should be `048`.
- sympy banner `:31`: `STAGE 48` — correct stage, 2-digit (zero-pad to `048` for consistency, deferred).
- mathematica `.txt` final line `:42`: `Stage 048 Mathematica audit passed.` — already correct; the `.txt` banner line 4 (`STAGE 031`) is the stale one (covered by F1's re-run).

These are label-only, do not affect any assertion, and belong to the gated numbering cleanup, not this script-correctness audit.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration. Two independent devices distinguish the engines:
1. **Inverse map.** SymPy *hardcodes* `zeta_req = (Sreq-1)/(1+eps(Sreq-2))` (L92) and substitutes back. Mathematica instead *solves* `Solve[sEnhance == sReq, zeta, Reals]` (L86), asserts exactly one solution exists (`If[Length[...] != 1, fail...]`, L87 — a check SymPy lacks entirely), then verifies. This is an independent derivation of the same boxed `\zeta_{\rm req}`, plus an extra uniqueness assertion (B8) that strengthens the "unique threshold" claim of appendix row 74.
2. **Pole divergence.** SymPy substitutes a concrete physical parametrization `eps=1/(1+nu)` and checks `limit S = +oo` (L84–90). Mathematica instead extracts the *pole residue* `Limit[(1/eps-zeta) S, zeta→1/eps] = (1-eps)/eps^2` (L79–84) — a different, stronger characterization of the same divergence. SymPy has no pole-coefficient check; Mathematica has no `eps=1/(1+nu)` substitution.

The shared identity checks use the same typed-in RHS forms (expected for verifying the same paper formula), but the two distinct derivation routes above rule out a line-by-line port.

## Engine cross-check

Both transcripts report every check as zero/PASS and agree on every printed closed form (sign and factor):
- `G_tr`, `M_crit = 9(1+delta)/(9+9delta+2R^2)` — agree (sympy L9–10, math L5–6).
- `dG_tr/dxi`, gap, `S`, `dS/dzeta`, margins, `dxi/dzeta` — all residual `= 0` / `PASS` in both (sympy L11–45, math L7–40).
- SymPy `limit zeta→(1/eps)^- of S under 0<eps<1 = oo` (L27) ↔ Mathematica `pole coefficient for S = (1-eps)/eps^2` (L23): consistent characterizations of the +∞ pole. Engines agree.

`engines_agree: true`.

## Verdict justification

`findings`, driven by a single low-severity `stale_output` (F1): both committed transcripts predate the latest banner edit (`STAGE 31`/`STAGE 031` in the `.txt` vs `STAGE 48`/`STAGE 048` in the scripts; output mtimes are 8 days older than the scripts). The mathematics is sound and fully aligned with the paper. I attacked the inverse-map checks (A8/A9/B9) as potential tautologies — they are not, because `S` is built independently of `zeta_req`/`zeta_crit`, so a wrong inverse formula would fail. I attacked the margin/derivative checks as possibly verifying a "different identity than the Output" — they verify the exact boxed `\zeta_{\rm req}` and the strict-monotone/margin facts that make the `\zeta_{\rm phys}\ge\zeta_{\rm req}` success condition meaningful, matching the card and notes. I checked symbol domains (`positive=True` on `xi,delta,R,eps,zeta` in sympy; explicit `0<xi<1, 0<eps<1, sCrit>sReq>1` in math): the SymPy margin checks are pure algebraic-form identities, so the absent ordering assumptions don't invalidate them; the `>0` positivity lives in prose and follows from the verified factored forms under the stated domain. I confirmed the `.wl` is an independent derivation (Solve + uniqueness; pole residue), not a transliteration. No `paper_misalignment`, no tautology, no engine disagreement.

## Self-test notes

Traps checked: (1) **Variable independence** — every `diff(EXPR, VAR)` differentiates w.r.t. a symbol the expression genuinely depends on (`G_tr` and `F_tr` depend on `xi`; `S` depends on `zeta`); no identically-zero derivative masquerading as a pass. (2) **Tautology** — the inverse-map and Solve checks build `S` independently of the `zeta_req`/`zeta_crit` forms, so they can fail; the `expect_zero` RHS forms are independently typed, not echoes of the LHS construction. (3) **Domain/branch** — `Limit[..., Direction→"FromBelow"]` and the `eps=1/(1+nu)` parametrization correctly pin the physical `0<eps<1` side, resolving the `-oo*sign(eps-1)` ambiguity SymPy prints unguarded (L26) into `+oo` (L27). Conclusion: math is correct and paper-aligned; only the committed transcripts are stale (banner labels), hence the single low-severity finding routed to a verifier re-run, with the in-source docstring/zero-pad labels deferred to the dedicated numbering pass.

## Value Reconciliation (pass-2 augmentation)

All of stage 048's deliverables are **symbolic** closed forms (no pinned numeric constants). The appendix `\zeta_{\rm req}=1/3` numerics (part03 lines 32/341/359) are the downstream Family-1 *precursor* evaluation, not stage 048's output, so they are excluded as out-of-scope for this stage's reconciliation.

| value (symbolic deliverable) | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `zeta_req = (S_req-1)/(1+eps(S_req-2))` | py L92 / wl L86–88; sympy out L34, math out L26 | stage_048.tex:22 (boxed); notes:99,182 | MATCH |
| `zeta_crit = (S_crit-1)/(1+eps(S_crit-2))` | py L93 / wl L89; sympy out L35, math out L27 | notes:190 | MATCH |
| success condition `zeta_phys >= zeta_req` | (logical; covered by margin checks A10/A11/B11/B12) | stage_048.tex:28 (boxed eq. success-iff) | MATCH |
| `G_tr = 9 xi(xi+delta)/(9 delta+(9+2R^2)xi)` | py L39; sympy out L9, math out L5 | notes:36–37 | MATCH |
| `M_crit = 9(1+delta)/(9+9delta+2R^2)` | py L47; sympy out L10, math out L6 | notes:55 | MATCH |
| `dG_tr/dxi = 9(2R^2 xi^2+9delta^2+18 delta xi+9 xi^2)/(2R^2 xi+9delta+9xi)^2` | py L52–56; sympy out L11, math out L7 | notes:45–47 | MATCH |
| `M_crit-G_tr` factored gap | py L69–74; sympy out L17, math out L15 | notes:59–61 | MATCH |
| `S(zeta;eps) = 1+zeta(1-eps)/(1-zeta eps)` | py L77; sympy out L23, math out L18 | notes:75 | MATCH |
| `dS/dzeta = (1-eps)/(1-zeta eps)^2` | py L82; sympy out L25, math out L21 | notes:85 | MATCH |
| pole margin `1/eps-zeta_req = (1-eps)/(eps(1+eps(S_req-2)))` | py L107–110; sympy out L36, math out L32 | notes:103–104 | MATCH |
| branch margin `zeta_crit-zeta_req = (S_crit-S_req)(1-eps)/((1+eps(S_crit-2))(1+eps(S_req-2)))` | py L111–114; sympy out L37, math out L33 | notes:194–196 | MATCH |
| `dxi_phys/dzeta` (implicit, closed form) | py L118–125; sympy out L44, math out L38 | notes:210–211 (form) | MATCH |
| pole coefficient of S `(1-eps)/eps^2` (math only) | wl L84; math out L23 | notes:93 (supports `S→+oo`) | MATCH (supporting) |

INTERNAL / scaffolding (no prose finding expected): `F_tr` softening coefficient `(9delta+9+2R^2)^2(9delta+9+2R)^2/(81(9delta^2+18delta+9+2R^2)^2)` (an intermediate establishing the `F_tr→+oo` endpoint of notes §3, not a boxed deliverable); pass/fail flags; residual `= 0` lines; Solve solution count.

`reconciliation: complete; 13 values checked, 0 misaligned.`
