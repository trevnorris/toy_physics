---
unit_id: 207
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table.md]
  paper_appendix: present
---

# Audit unit 207 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_207.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 45, 236, 1337)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Exact primitive Hessian-envelope theorem, sign-adapted primitive drift table, certified primitive ray table, primitive elimination theorem, and primitive winner theorem." Stage 207 specializes the Stage 206 certified ray sieve to the five primitive log axes `(e_λ, e_c, e_γ, e_U, e_W)`. The notes enumerate six deliverables, of which the *certifiable* (script-checkable) ones are: (1) the canonical orientation rule `ε_i = -sgn(Γ_i)`, `K_i = -|Γ_i|`; (2) the primitive Hessian-envelope theorem `(ε_i e_i)^T H_h (ε_i e_i) = ∂_{ℓ_iℓ_i}` (diagonal-only, no off-diagonal); (3) the sign-adapted primitive drift table (Stage 204 exponents × orientation ε); (4) the certified monotone bracket root map `T(H0,k;c) = 2H0/(k+√(k²-2cH0))` solving `H0 - kτ + ½cτ² = 0`; (5) the certified turning bracket `√(-2H0/κ)`, κ<0; (6) the mixed-ray preview `H_s = a²∂_ii + 2ab∂_ij + b²∂_jj` showing off-diagonal entries first enter only for two-coordinate rays. The elimination/winner theorems are inequality-logic statements over those primitives, not residual identities. Card is checkpoint: False; status Mixed (exact closure + numerical). Card explicitly says "Mathematica audit: none yet" — the `.wl` is a pass-1 retrofit.

## What the script claims to verify

SymPy: §I diagonal Hessian reduction `e_i^T H e_i = H_ii` for both ±e_i; §II mixed quadratic form `(ae_i+be_j)^T H (ae_i+be_j) = a²H_ii+2abH_ij+b²H_jj` and isolates the cross term `2abH_lc`; §III orientation law `-sgn(Γ)·Γ + |Γ| = 0`; §IV substitutes `ε_i e_i` into generic Stage-204 drift forms and matches the five expected rows; §V posits the monotone root closed form and checks it solves the comparison quadratic; §VI posits the turning root and checks it; §VII prints the five-row certified table. Mathematica: M1 diagonal reduction, M2 mixed form + cross-coefficient over all 10 axis pairs, M3 orientation, M4 *Solves* the comparison quadratic and *selects* the branch then matches the closed form, M5 *Solves* the turning quadratic + selects the positive root then matches `√(-2H0/κ)`, M6 the drift table. These faithfully exercise the paper's load-bearing identities.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Canonical orientation `ε=-sgn(Γ)`, `K=-|Γ|` | py §III L94-99; wl M3 L123-125 | match |
| Primitive Hessian-envelope (diagonal-only) | py §I L61-67; wl M1 L96-102 | match |
| Off-diagonal first appears on mixed rays | py §II L77-85; wl M2 L105-120 | match |
| Sign-adapted primitive drift table (5 rows) | py §IV L162-169; wl M6 L223-231 | match |
| Certified monotone bracket root map | py §V/§VII L177-190,238-242; wl M4 L127-174 | match |
| Certified turning bracket | py §VI L201-213; wl M5 L176-198 | match |
| Primitive elimination theorem (inequality logic) | (not residual-checkable; structural) | n/a — paper carries as logical corollary of bracket admissibility, no identity to assert |
| Primitive winner theorem (pairwise ordering) | (not residual-checkable; structural) | n/a — same |

All residual-checkable deliverables are covered by both engines. The elimination/winner theorems are logical consequences of the certified-bracket admissibility test (`Δ≥0`, `τ_hi≤T_i`) plus Stage 206's pairwise ordering; they contain no new symbolic identity to assert, so their absence from the scripts is not a gap. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 66-67 | `expect_zero(e_i^T H e_i - H_ii)` ±e_i | Hessian-envelope (diag) | yes |
| A2 | sympy | 81 | `mixed_expr - expected_mixed == 0` | mixed-ray preview | yes |
| A3 | sympy | 82-85 | off-diag first-appearance == 0 | mixed-ray preview | yes |
| A4 | sympy | 99 | `K_oriented + Abs(Gam) == 0` | orientation law | yes |
| A5 | sympy | 169 | drift row - ε·expected == 0 (×5) | drift table | yes |
| A6 | sympy | 183-190 | `H0 - k·Tau + ½c·Tau² == 0` (lo,hi) | monotone bracket | yes |
| A7 | sympy | 207-213 | `H0 - ½a·Tau² == 0` (lo,hi) | turning bracket | yes |
| M1 | mathematica | 98-99 | `e.H.e - diag == 0` ±e_i | Hessian-envelope (diag) | yes |
| M2 | mathematica | 113,114-117 | form residual + `Coefficient[form,ab]-2H_ij == 0` (10 pairs) | mixed-ray preview | yes |
| M3 | mathematica | 125 | `orientedSlope + Abs[Gam] == 0` | orientation law | yes |
| M4 | mathematica | 171-174 | Solve+branch-select root, then `τ - closedForm == 0` + quadratic residual | monotone bracket | yes (derive) |
| M5 | mathematica | 195-198 | Solve+positive-root, then `τ - √(-2H0/κ) == 0` + turning residual | turning bracket | yes (derive) |
| M6 | mathematica | 225-228 | drift row - ε·expected == 0 (×5) | drift table | yes |

Every script-side check traces to a specific paper deliverable. No orphaned scaffolding. A1/A2/M1/M2/M6/A5 are elementary linear-algebra/substitution identities on shared monomial premises (non-tautological in the trivial sense — they would fail if H were mis-indexed or a row were mis-typed — but low discriminating power). A6/A7/M4/M5 are the load-bearing root-map checks.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.txt` (mtime 2026-05-11 12:49:16)
- vs `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.py` (mtime 2026-06-03 15:59:11)

**What's wrong:**
The committed SymPy output is older than the SymPy script by ~23 days. The captured transcript still prints the stale numbering banner `STAGE 190 — PRIMITIVE-RAY HESSIAN ENVELOPES AND CERTIFIED RAY TABLE` (output L11) and `STAGE 190 SYMPY AUDIT PASSED` (L242), whereas the current script prints `STAGE 207 …` (py L35) and `STAGE 207 SYMPY AUDIT PASSED` (py L244). Apart from the banner label, the structural content of the transcript (sections I–VII, all residuals = 0, `EXIT_CODE: 0`) matches what the current script would emit, so the staleness is a label/refresh artifact, not a math disagreement. The Mathematica output (mtime 2026-06-02 10:07:49) is fresh relative to its `.wl` (mtime 2026-06-02 10:06:44) and already carries the correct `STAGE 207` banner.

**Why this matters:**
The committed `.txt` is the authoritative record of what the script produces; a stale transcript with a wrong stage label could mislead later reconciliation and masks whether a banner/renumber fix has actually been re-captured.

**Required change:**
Re-run the SymPy script and re-capture the output so the transcript reflects the current `STAGE 207` banner. No script-math edit is required.

**Verification:**
Fresh `scripts/output/..._sympy_audit.txt` shows `STAGE 190` replaced by `STAGE 207` on the banner lines, mtime newer than the `.py`, `EXIT_CODE: 0`.

## Independent-derivation check (Mathematica)

Verdict: **INDEPENDENT.** The retrofit `.wl` is not a line-by-line port. The discriminator is the bracket root maps (the load-bearing objects of §V/§VI / M4/M5):

- SymPy §V (py L177-190) **posits** the closed form and verifies a residual:
  `TauL = sp.simplify(2*H0/(k + sp.sqrt(k**2 - 2*cL*H0)))` … `expect_zero("lower comparison quadratic", H0 - k*TauL + Rational(1,2)*cL*TauL**2)`.
- Mathematica M4 (wl L131-174) **derives** the root by solving and branch-selecting, then matches the closed form:
  `roots = stripConditional[tau /. Solve[comparisonPolynomial[cSym, tau] == 0, tau]]`; `positiveBranch = SelectFirst[roots, … root>0 && other>0 && root<other …]`; then `expectZero["M4 tau_lo selected root - closed form", tauLo - closedMonotoneRoot[kappaLo]]`. The `.wl` additionally enumerates the positive-root set via `Reduce[… && tau>0, tau, Reals]` and covers both `c>0` and `c<0` sign branches (Piecewise L160) — coverage the SymPy script does not attempt (it pins `cL,cU` positive).

Same contrast for the turning bracket: SymPy §VI posits `sqrt(2*H0/a_turn)` and checks the residual; M5 (wl L179-197) does `Solve[H0 - 1/2 (-curvature) tau^2 == 0, tau]`, `SelectFirst` the positive root, then `expectZero[…, tauTurnLo - Sqrt[-2 H0/kappaTurnLo]]`. This is derive-by-Solve vs posit-and-residual — genuinely independent extraction of the load-bearing root.

For the elementary blocks (M1/M2 quadratic-form expansion, M3 orientation, M6 drift substitution) both engines perform the same operation, but those share only the monomial/premise definitions (the 5×5 symmetric Hessian, the Stage-204 drift forms) and are the only operation those identities admit; sharing the premise is permitted. The independence test is met where it matters: the root-map (the object that actually carries the certified bracket) is posited in one engine and Solve-derived in the other. Rejecting the "each CAS runs its own simplifier" defense, this is still INDEPENDENT because the *method* differs (Solve+branch-selection vs posit+residual), not merely the simplifier.

## Engine cross-check

Both engines pass all checks (SymPy `EXIT_CODE: 0`; Mathematica `STAGE 207 MATHEMATICA AUDIT PASSED`, `Exit[0]`). Final symbolic forms agree: monotone bracket `2H0/(k+√(k²-2cH0))` (py output L166-177; the `.wl` `closedMonotoneRoot` L129 is identical); turning bracket `√2·√H0/√a_turn = √(-2H0/κ)` (py output L192-203; wl `Sqrt[-2H0/kappaTurnLo]` L195). Diagonal reductions all `0`, mixed cross-coefficient `2H_ij`, orientation `0`, all five drift rows `0` in both. No engine disagreement.

## Verdict justification

`findings`: one low-severity `stale_output` (SymPy `.txt` predates the `.py`, still showing the `STAGE 190` banner). Everything else holds up. Attacks tried that failed: (a) transliteration attack on the retrofit `.wl` — defeated, the root maps are Solve-derived in Mathematica vs posited+residual-checked in SymPy; (b) tautology attack — the diagonal/mixed/drift checks would fail under mis-indexing or a mis-typed row, and the bracket checks substitute the candidate back into the *independently stated* comparison quadratic (not into its own definition), so they are not algebraically guaranteed; (c) symbol-domain attack — SymPy pins `cL,cU` positive but the residual identities are sign-independent algebra, and the `.wl` covers both curvature signs, so no domain error is load-bearing; (d) parity/branch attack on the root selection — the `.wl` selects the smaller positive root for c>0 (matching the rationalized `2H0/(k+√…)`) and the unique positive root for c<0, both consistent with the closed form; (e) paper-misalignment / value attack — every emitted deliverable value reconciles to the card and notes (see below). I read the paper card, notes, and appendix; the script's claims match the paper's stated deliverables.

## Self-test notes

Checked: (1) variable independence — no `diff`/`D` of a non-dependent variable here; the only derivatives are quadratic-form expansions and `Solve`s, all over variables the expressions genuinely contain. (2) Parity/symmetry — n/a (no unbounded integrals); the relevant branch choice (smaller-positive root for c>0) was verified to match the rationalized closed form. (3) Trivial-case pre-check — substituting the candidate roots back into the comparison/turning quadratics gives `H0 - H0 = 0` analytically, confirming the residual checks are real and pass for the right reason. The only finding is a refresh (`stale_output`); no `paper_misalignment`, so no `## Resolve before fix_loop` block is needed.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 7 deliverable values checked, 0 misaligned.

All values the scripts emit are **symbolic closed forms / table rows** (no pinned floats, no benchmark numbers — this is a pure exact-closure compiler stage). Each reconciles to the paper card / notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Diagonal reduction `(ε_i e_i)^T H (ε_i e_i) = ∂_{ℓ_iℓ_i}` | py §I L61-67; wl M1; out py L17-36 / wl L9-28 | notes §3 L153-167 (`\mathcal H_{ŝ_i}=∂_{ℓ_iℓ_i}`); card Output "primitive Hessian-envelope theorem" | MATCH |
| Mixed form `a²H_ii+2abH_ij+b²H_jj`, cross coeff `2H_ij` | py §II L77-85; wl M2; out py L41-44 | notes §8 L397-405 | MATCH |
| Orientation `ε=-sgn(Γ)`, `K=-|Γ|` | py §III L93-99; wl M3 L123-125; out py L49-59 | notes §2.1 L117-128 (`ε_i:=-sgn(Γ_i)`, `K_i=-|Γ_i|`) | MATCH |
| Drift rows λ,c,γ,U,W (Stage-204 × ε) | py §IV L125-169; wl M6 L215-231 | notes §4 table L217-221; card Output "sign-adapted primitive drift table" | MATCH |
| Monotone root `T=2H0/(k+√(k²-2cH0))` | py §V L177-178, §VII L239-240; wl M4 L129; out py L166-177 | notes §5 L240 (root map `\mathcal T`); §5.1 L261-263 | MATCH |
| Monotone width `W_i=τ_hi-τ_lo` | py §V L191-193; out py L180-187 | notes §5.1 L268 (`W_i:=τ_{i,hi}-τ_{i,lo}`) | MATCH |
| Turning roots `√(-2H0/κ)`, width | py §VI L201-213; wl M5 L195,197; out py L192-203 | notes §5.2 L284-291 (`τ^{(tp)}:=√(-2H0/κ)`) | MATCH |

INTERNAL (scaffolding, no prose carrier expected): 5×5 symbolic Hessian symbols `Hll…HUW`; generic drift form symbols `s_λ…s_W`, `ε_λ…ε_W`; `a_star=(1+δU*)/(1+χ0*)` (intermediate carried from Stage 204, surfaced in notes §4 L209-211 as `\mathfrak a_*` — reconciles but used internally); residual `=0` print values; PASS flags; banner labels.

Note: the SymPy output `.txt` is stale (see F1) and prints the deliverables under a `STAGE 190` banner, but the *values themselves* (root maps, widths, drift rows, reductions) in that transcript are identical to what the current script defines and all reconcile; the staleness is a banner-label artifact only. Reconciliation is therefore based on script source confirmed against the (matching) transcript content.
