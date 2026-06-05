---
unit_id: 054
batch: III.2
auditor_model: Opus 4.8 (1M context)
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage054_robin_softening_support_lane.md]
  paper_appendix: present
---

# Audit unit 054 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_054.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage054_robin_softening_support_lane.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 86; `\input` at line 226)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage054_robin_softening_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage054_robin_softening_mathematica_audit.txt`

## What the paper claims

Stage 054 replaces the Dirichlet mouth of the lowest support lane by a finite Robin compliance and quantifies the resulting support-stiffness softening. The card states three boxed deliverables. (1) The lowest Robin root satisfies `\boxed{y\tan y=\eta, \; 0<y<\pi/2}` (eq:app-stage054-robin-root). (2) With `x=\pi^2 T_X/(L^2 K_W^{\rm eff})`, `0<x<4`, the softening factor is `\boxed{A_K(\eta)=K_W^{\rm eff}/K_{\phi,0}^{\rm eff}=1/(1-x/4+x y^2/\pi^2)}` (eq:app-stage054-AK). (3) The endpoint window/ceiling is `\boxed{1\le A_K\le 4/(4-x)}` (eq:app-stage054-AK-window). `\stagefield{Output}` reads: "Robin root \eqref{eq:app-stage054-robin-root}, softening factor \eqref{eq:app-stage054-AK}, and ceiling \eqref{eq:app-stage054-AK-window}." The notes add the supporting pieces: the stiffnesses `K_W^{eff}=K_X+\pi^2 T_X/(4L^2)`, `K_{\phi,0}^{eff}=K_X+T_X y^2/L^2`; the endpoint values `A_K=1` at `y\to\pi/2` and `A_K,max=4/(4-x)` at `y\to0^+`; strict monotonic decrease of `A_K` in `y`; the rescue criterion `\zeta_{req}\le4/(4-x)` ⇔ `x\ge4-4/\zeta_{req}`; and the exact threshold `y_{req}^2=(\pi^2/x)(1/\zeta_{req}-1+x/4)`.

## What the script claims to verify

Both engines derive the Robin characteristic equation from the BVP ansatz `\psi=A\cos ks+B\sin ks` with the Neumann bottom solved for B, then assert the mouth condition reduces to `k\tan(kL)-h=0` and its dimensionless form `y\tan y-\eta=0`. They then build `K_W^{eff}` and `K_{\phi,0}^{eff}` from `K_X,T_X`, derive the x-parametrized `A_K` independently (substituting `K_X=K_W(1-x/4)`, `T_X=xL^2K_W/\pi^2`) and assert it equals `1/(1-x/4+xy^2/\pi^2)`; assert the DN endpoint `A_K(\pi/2)=1` and soft-mouth endpoint `A_K(0^+)=4/(4-x)`; verify the closed form of `dA_K/dy` (monotonicity certificate); and assert the saturation `x`-floor `4-4/\zeta_{req}`. The Mathematica script adds one extra assertion (`A_K,max(x_floor)-\zeta_{req}==0`) confirming the floor inverts the ceiling. These assertions faithfully cover all three boxed deliverables plus the notes' supporting results.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Robin root `y\tan y=\eta` (boxed) | sympy 39/42, math 46/49 derive char eq and dimensionless form, assert == 0 | match |
| Softening factor `A_K=1/(1-x/4+xy^2/\pi^2)` (boxed) | sympy 60-62 (`A_K x-form`), math 69/77 derive from K_X,T_X and assert == target | match |
| Ceiling `1\le A_K\le 4/(4-x)` (boxed) | endpoints: sympy 69 (`DN limit`→1) + 70 (`soft-mouth`→4/(4-x)); monotonicity: sympy 78 (`dA_K/dy`); math mirrors 78/79/87 | match |
| `K_W^{eff}=K_X+\pi^2T_X/(4L^2)`, `K_{\phi,0}^{eff}=K_X+T_Xy^2/L^2` (notes) | sympy 48-53 print + `K_W identity` 59; math 57-67 | match |
| rescue floor `x\ge4-4/\zeta_{req}` (notes) | sympy 92-94 (`x floor`), math 93-98 | match |
| `A_K,max(x_floor)=\zeta_{req}` (notes, implicit) | math 99 only (not in sympy) | match (extra coverage in one engine; not a defect) |

`paper_alignment: aligned` — every boxed deliverable and every notes-level supporting result has a non-tautological script-side check, and all asserted values agree with the docs.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `expect_zero(charEq/A + h - k tan(kL))` | Robin root (dimensional) | yes |
| A2 | sympy | 42 | `expect_zero(dimensionless form)` | Robin root `y\tan y=\eta` | yes |
| A3 | sympy | 59 | `expect_zero(K_W identity)` | x-parametrization consistency | partial (algebraic identity of the chosen substitution; supports A4) |
| A4 | sympy | 62 | `expect_zero(A_K x-form)` | softening factor (boxed) | yes |
| A5 | sympy | 69 | `expect_zero(DN limit → 1)` | ceiling lower endpoint | yes |
| A6 | sympy | 70 | `expect_zero(soft-mouth → 4/(4-x))` | ceiling upper endpoint | yes |
| A7 | sympy | 78 | `expect_zero(dA_K/dy closed form)` | monotonicity (closes window) | yes |
| A8 | sympy | 94 | `expect_zero(x floor − (4−4/ζ))` | rescue floor | yes |
| B1 | math | 46 | `expectZero(charEq/a + h − k Tan[k ell])` | Robin root (dimensional) | yes |
| B2 | math | 49 | `expectZero(dimensionless form)` | Robin root | yes |
| B3 | math | 67 | `expectZero(K_W identity)` | x-parametrization consistency | partial (mirror of A3) |
| B4 | math | 77 | `expectZero(A_K x-form)` | softening factor (boxed) | yes |
| B5 | math | 78 | `expectZero(DN limit − 1)` | ceiling lower endpoint | yes |
| B6 | math | 79 | `expectZero(soft-mouth − 4/(4−x))` | ceiling upper endpoint | yes |
| B7 | math | 87 | `expectZero(dA_K/dy closed form)` | monotonicity | yes |
| B8 | math | 98 | `expectZero(x floor − (4−4/ζ))` | rescue floor | yes |
| B9 | math | 99 | `expectZero(A_K,max(x_floor) − ζ)` | floor↔ceiling inverse | yes |

A3/B3 are genuine algebraic identities of the chosen `(K_X,T_X)` substitution rather than independent physics — but they are not load-bearing on their own; the load-bearing softening claim is A4/B4, which derives `A_K` from the stiffnesses and matches the boxed target non-trivially (output line 18: `4*pi**2/(4*x*y**2 - pi**2*(x-4))` reduces to the target). No `tautological_check` finding warranted.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage054_robin_softening_sympy_audit.txt` (mtime 2026-05-26 03:00:52) vs script `...sympy_audit.py` (mtime 2026-06-03 15:59:11)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage054_robin_softening_mathematica_audit.txt` (mtime 2026-05-26 03:00:58) vs script `...mathematica_audit.wl` (mtime 2026-06-03 15:59:11)

**What's wrong:**
Both saved transcripts predate their scripts by ~8 days, and the captured content disagrees with the current scripts. The SymPy transcript line 3 banner reads `STAGE 37 — ROBIN-COMPLIANCE SOFTENING`, but the current script (`...sympy_audit.py:25`) emits `banner("STAGE 54 — ROBIN-COMPLIANCE SOFTENING")`. The Mathematica transcript line 4 reads `STAGE 037 — ROBIN-COMPLIANCE SOFTENING`, but the current `.wl:32` emits `banner["STAGE 054 — ROBIN-COMPLIANCE SOFTENING"]`. So the committed outputs were produced by an earlier revision of the banners; they are stale.

**Why this matters:**
The committed transcript is the authoritative record for this paper's reproducibility claim; a stale transcript whose top line names a different stage undermines that record and could mask a real regression if one were introduced between revisions. The mathematics in the stale outputs still corroborates all PASS lines, so this is not a correctness regression — but the transcript must be refreshed from the current scripts.

**Required change:**
Re-run both scripts and overwrite the two `.txt` transcripts with current output (the orchestrator does this on the independent re-run). No assertion logic changes.

**Verification:**
After refresh, SymPy `.txt` line 3 should read `STAGE 54 — ...` (or `STAGE 054 — ...` if F2 is also applied) and Mathematica `.txt` line 4 should read `STAGE 054 — ...`; all PASS/`= 0` lines remain.

### F2 — stale_output (numbering self-labels)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py:3` — docstring `"""Stage 37 SymPy audit: Robin-compliance softening of the lowest support lane."""`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py:25` — `banner("STAGE 54 — ROBIN-COMPLIANCE SOFTENING")` (two-digit, not 3-digit `054`)

**What's wrong:**
The SymPy module docstring carries a stale stage self-label `Stage 37` (this is stage 054; the `Stage 37`/`STAGE 037` band is the known +17 pre-renumber drift). The SymPy banner also uses two-digit `STAGE 54` where the Mathematica banner uses the canonical three-digit `STAGE 054`. These are the low-severity numbering self-label class; non-blocking. The Mathematica `.wl` banner (line 32) is already correct (`STAGE 054`).

**Why this matters:**
Cosmetic only — no math impact. Flagged per the standing in-loop Reading-2 policy (a verdict:findings stage gets its unambiguous self-labels fixed). It also makes the refreshed transcript banner self-consistent (`STAGE 054`).

**Required change:**
In `...sympy_audit.py`: change docstring line 3 `Stage 37` → `Stage 054`, and banner line 25 `"STAGE 54 — ..."` → `"STAGE 054 — ..."`. Re-run to refresh the transcript.

**Verification:**
`grep -n "Stage 37\|STAGE 54[^0-9]" ...sympy_audit.py` returns nothing; refreshed `.txt` line 3 reads `STAGE 054 — ROBIN-COMPLIANCE SOFTENING`.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. It uses Mathematica-native tooling rather than echoing SymPy's algebra: `bExpr` is obtained via `First@Solve[(D[psi,s]/.s->ell)==0, b]` (.wl:40) where SymPy uses `sp.solve(...)[0]` then a `.subs` (`.py:33`); the x-form is verified with `FullSimplify[..., Assumptions->$Assumptions]` over a declared `0<x<4` domain (.wl:69-77) and the soft-mouth endpoint via `Limit[..., Direction->"FromAbove"]` (.wl:72) versus SymPy's `sp.limit(..., dir="+")` (.py:66). The Mathematica script also adds an assertion the SymPy script lacks (B9, `A_K,max(x_floor)-zeta_req==0`, .wl:99), which an echo would not introduce. The `expectZero` harness independently strips `ConditionalExpression` wrappers (.wl:26) — a Mathematica-specific concern absent in SymPy. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on every shared deliverable. Robin char eq: both reduce `charEq/A` to `k\tan(kL)-h` (sympy out line 6; math out line 6, modulo the Mathematica `ConditionalExpression` domain wrapper) and assert the dimensionless residual == 0. Softening factor: sympy `A_K x-form = 0` (out line 19), math `A_K x-form = 0` (out line 22). Endpoints: both give `A_K(\pi/2)=1` and `A_K(0^+)=4/(4-x)` (`-4/(x-4)` form; sympy out 20-21, math out 19-20). Monotonicity `dA_K/dy closed form = 0` in both (sympy out 24, math out 27). x-floor `4-4/\zeta_{req}` in both (sympy out 32, math out 32). No `engine_disagreement`.

## Verdict justification

verdict: findings (both low-severity `stale_output`). The mathematics holds up under attack: every boxed deliverable (Robin root, softening factor, ceiling) and every notes-level supporting result (stiffnesses, endpoints, monotonicity, rescue floor) is exercised by a non-tautological, well-anchored assertion in BOTH engines, and all asserted values match the `.tex` and `.md` verbatim. Attacks tried and failed: (a) tautology probe on `A_K x-form` — defeated, since `A_K` is rebuilt from `K_X,T_X` via the x-substitution and the output shows a non-trivial intermediate form reducing to the boxed target; (b) symbol-domain probe — `y>0`, `0<x<4` assumptions match the paper's stated ranges, and the `y→π/2`/`y→0+` endpoint substitutions are within domain; (c) hidden-`simplify`-assumption probe — the Mathematica `Limit::alimv` warning is benign (the soft-mouth `Assumptions->0<x<4` does not involve the limit variable `y`); (d) engine-disagreement probe — none. The only defects are the two stale-output items: the committed transcripts predate the scripts and carry the old `STAGE 37`/`STAGE 037` banner, and the SymPy source still carries a stale `Stage 37` docstring self-label plus a two-digit banner. Neither affects correctness. No `paper_misalignment`, no stop-cold.

## Self-test notes

Checked (1) variable independence: `dA_K/dy` depends on `y` (and `x`), `y` is a real diff variable — the derivative is genuinely nonzero (prefactor `-2xy/π²`), so the closed-form match is non-trivial. (2) Trivial-case pre-check: at `y=π/2`, `1-x/4+x(π/2)²/π² = 1-x/4+x/4 = 1` ⇒ `A_K=1` ✓; at `y→0`, denominator → `1-x/4` ⇒ `A_K=4/(4-x)` ✓ — endpoint asserts are correctly true. (3) Path specs: F1/F2 are output-refresh + cosmetic self-label edits on existing files; no new-script path needed. (4) Paper round-trip: the only prescribed edits are transcript refresh and a stage-number self-label correction — neither touches any load-bearing constant, so no new `paper_misalignment` is introduced.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit against the `.tex` card and `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Robin char eq `k\tan(kL)=h` | py:38-39 / wl:45-46 (out py L6, wl L6) | .tex:13-20 (eq:app-stage054-robin-root, dimensionless); .md:68-84 | MATCH |
| Dimensionless root `y\tan y=\eta`, `0<y<\pi/2` | py:42 / wl:49 (out py L8, wl L9-10) | .tex:16-19; .md:80-84 | MATCH |
| `K_W^{eff}=K_X+\pi^2T_X/(4L^2)` | py:50 / wl:58,61 (out py L13, wl L11) | .md:97 | MATCH (notes carrier; .tex terse, OK per guard) |
| `K_{\phi,0}^{eff}=K_X+T_Xy^2/L^2` | py:52 / wl:57,62 (out py L14, wl L12) | .tex:31 (in `A_K` def); .md:101 | MATCH |
| `x=\pi^2T_X/(L^2K_W^{eff})`, `0<x<4` | py:55-58 / wl:65-66 (out py L16) | .tex:22-27 (eq:app-stage054-x); .md:109,117 | MATCH |
| Softening factor `A_K=1/(1-x/4+xy^2/\pi^2)` | py:61-62 / wl:74,77 (out py L18-19, wl L18,22) | .tex:29-34 (eq:app-stage054-AK, boxed); .md:33-34,122 | MATCH |
| `A_K(y=\pi/2)=1` | py:65,67 / wl:71,75 (out py L20, wl L19) | .tex:38 (window lower); .md:134 | MATCH |
| `A_K(y\to0^+)=4/(4-x)` | py:66,68 / wl:72,76 (out py L21, wl L20) | .tex:38 (window upper); .md:140 | MATCH |
| Window/ceiling `1\le A_K\le 4/(4-x)` | py:69-70 / wl:78-79 + monotonicity py:78/wl:87 | .tex:36-40 (eq:app-stage054-AK-window, boxed); .md:144,223 | MATCH |
| `dA_K/dy=-2xy/(\pi^2(1-x/4+xy^2/\pi^2)^2)` | py:77-78 / wl:86-87 (out py L24, wl L27) | .md:148 (states strictly decreasing; closed form is supporting) | MATCH |
| `y_{req}^2=(\pi^2/x)(1/\zeta_{req}-1+x/4)` | py:85-86 / wl:91,95 (out py L30, wl L30) | .md:191 | MATCH |
| `A_K,max=4/(4-x)` | py:89-90 / wl:92,96 (out py L31, wl L31) | .tex:39; .md:140,163 | MATCH |
| x-floor `=4-4/\zeta_{req}` | py:92-94 / wl:93,97-98 (out py L32-33, wl L32-34) | .md:173 | MATCH |
| Rescue criterion `\zeta_{req}\le4/(4-x)` | py:100 (final ledger) / paper text | .tex:41; .md:169,228 | MATCH |

INTERNAL (scaffolding, no prose expected): `B from Neumann bottom = A*tan(L*k)` / `bExpr` (intermediate BVP solution); `expect_zero`/`expectZero` PASS flags and residuals; banner strings; `K_X=K_W(1-x/4)`, `T_X=xL^2K_W/\pi^2` (x-parametrization substitutions, intermediate); `A_K,max(x_floor)-\zeta_{req}=0` (wl:99, internal inverse-consistency check); final-ledger summary prints.

reconciliation: complete; 14 values checked, 0 misaligned
