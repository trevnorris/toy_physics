---
unit_id: 214
batch: VI.1
auditor_model: Claude Opus 4.8 (1M context)
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md]
  paper_appendix: present
---

# Audit unit 214 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_214.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows: source-status map line 59; §subsec:app-part06-four-interior-and-splice lines 1058-1089)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.txt`

## What the paper claims

Stage 214 solves the three-parameter **interior** optimization of the certified
four-coordinate simplex objective. `\stagefield{Output}` reads:
"Finite quadruple-interior candidate set with preferred bound \(54\) per envelope,
projected fallback bound, and interior winner/non-improvement theorem against the
Stage~213 boundary ledger." The notes (§1-§6) and appendix
(§subsec:app-part06-four-interior-and-splice) enumerate the distinct deliverables:
(1) the certified interior objective `Φ_{Q,⋆}` / root `τ_{Q,⋆}` with discriminant
`Δ^♯`; (2) the three exact stationary-numerator derivative laws
`∂Φ/∂r = N_r/(2(1+r²+s²+t²)^{3/2}√Δ)` with `N_r = 2M_r√Δ + L_r`; (3) the lifted
polynomial system with degree ledger `(3,3,3,2)` and Bézout bound `3·3·3·2 = 54`
per envelope; (4) the square-root-free projected eliminants — three quintic
cross-consistency polynomials `C_rs,C_rt,C_st` and three sextic square conditions
`S_r,S_s,S_t`, giving the projected one-chart bound `5·5·6 = 150`; (5) the
diagonal-isotropic reduction (optimizer = gradient ray) and (6) the full-symmetry
reduction (equal-mix barycenter `(1,1,1,1)/2` is stationary); plus (7) the
interior winner / non-improvement order theorems against the imported Stage 213
boundary bracket `[β_lo, β_hi]`.

## What the script claims to verify

The SymPy script verifies all seven deliverables: (I) the three derivative-law
"compilers" `diff(Φ,·) − N_·/(2 den^{3/2}√Δ) = 0`; (II) the lifted polynomials
have total degrees `(3,3,3,2)` and product `54`; (III) the eliminants `C_··`,`S_·`
are non-vacuous, vanish at a concrete common stationary root
`(r,s,t,y)=(3/2,5/2,7/2,29/2)` with explicit coefficients, have degrees `(5,5,5)`
and `(6,6,6)`, product `5·5·6 = 150`; (IV) the diagonal-isotropic `Δ`/`τ`
reductions and gradient-ray unit-norm/slope; (V) symmetric `N_·(1,1,1)=0` and
equal-mix unit-norm; (VI) the winner and non-improvement order theorems by integer
brute-force enumeration (924 samples each). The Mathematica script (M1-M7) covers
the same algebraic claims I-V but does **not** include the winner/non-improvement
order theorems (VI).

## Paper ↔ script cross-check

| Paper deliverable | sympy | mathematica | status |
|---|---|---|---|
| (2) derivative laws `∂Φ/∂· = N_·/(2 den^{3/2}√Δ)` | A1-A3 | A8-A13 | match |
| (3) lifted degrees `(3,3,3,2)`, Bézout `54` | A4-A5 | A14-A18 | match |
| (4) eliminants `C_··` quintic, `S_·` sextic, bound `150` | A6-A7 | A19-A29 | match |
| (5) diagonal-isotropic + gradient ray | A8-A10 | A30-A33 | match |
| (6) full-symmetry equal-mix stationary | A11-A12 | A34-A37 | match |
| (7) winner/non-improvement order theorems | A13-A14 | — | partial (wl omits) |

`paper_alignment: aligned` — every deliverable is covered by at least the SymPy
engine; the order theorem (7) is single-engine (SymPy only), which is acceptable as
a trivial real-interval-ordering statement and noted, not a blocking finding.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 106 | `diff(Phi,r) − Nr/(2 den^{3/2}√Δ) == 0` | (2) | yes |
| A2 | sympy | 107 | `diff(Phi,s) − Ns/... == 0` | (2) | yes |
| A3 | sympy | 108 | `diff(Phi,t) − Nt/... == 0` | (2) | yes |
| A4 | sympy | 139-146 | `deg Fr,Fs,Ft==3, FDelta==2` | (3) | yes |
| A5 | sympy | 150 | `deg-product == 54` | (3) | yes |
| A6 | sympy | 187 | each eliminant non-vacuous | (4) | yes |
| A7 | sympy | 223 | all 10 polys vanish at concrete root | (3)(4) | yes |
| A8 | sympy | 239-247 | `deg C=5, S=6, product==150` | (4) | yes |
| A9 | sympy | 271-278 | iso `Δ`,`τ` reductions | (5) | yes |
| A10 | sympy | 284-285 | gradient ray norm/slope | (5) | yes |
| A11 | sympy | 316-318 | symmetric `N_·(1,1,1)==0` | (6) | yes |
| A12 | sympy | 321 | equal-mix norm | (6) | yes |
| A13 | sympy | 331-341 | winner order theorem (924 samples) | (7) | yes |
| A14 | sympy | 344-354 | non-improvement (924 samples) | (7) | yes |
| M8-M13 | wl | 107-112 | `directNum − stationaryNum==0`, `D[phi,·]−N/den==0` | (2) | yes |
| M14-M18 | wl | 121-130 | lifted degrees `(3,3,3,2)`, product `54` | (3) | yes |
| M19-M24 | wl | 150-155 | `Resultant[...,y]` eliminants == posited defs | (4) | yes |
| M25-M29 | wl | 161-172 | eliminant degrees `(5,5,5,6,6,6)`, product `150` | (4) | yes |
| M30-M33 | wl | 197-204 | iso `Δ`,`τ`, gradient ray | (5) | yes |
| M34-M37 | wl | 227-230 | symmetric `N_·(1,1,1)==0`, equal-mix norm | (6) | yes |

No tautological rows. A7 is a genuine non-trivial numeric anchor: the concrete
point with `A=-220/3, B=12, ... J=-85/3` at `(3/2,5/2,7/2,29/2)` is an actual
common root of all four lifted polynomials AND all six projected eliminants, so the
vanishing is not guaranteed by construction.

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt` (mtime 2026-06-02 11:59:40)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py` (mtime 2026-06-03 15:59:11)

**What's wrong:**
The committed SymPy `.txt` output predates the `.py` by ~28 h. The staleness is
confirmed by content, not just mtime: the saved output banner reads
`STAGE 197 — FULL INTERIOR FOUR-COORDINATE SIMPLEX OPTIMIZER` (line 3) and
`STAGE 197 SYMPY AUDIT COMPLETED SUCCESSFULLY` (line 1209), whereas the current
script's banner (`.py:41`) reads `STAGE 214 — FULL INTERIOR FOUR-COORDINATE SIMPLEX
OPTIMIZER`. The captured output is from the pre-renumber epoch and predates the
banner fix. (This is exactly the SCRIPT/OUTPUT-band stale-label class noted in the
project's numbering memory.) All *numerical/symbolic content* in the saved output
(degrees, `54`, `150`, the `=0` residuals, the `924`-sample sweep counts) matches
what the current script would emit — only the stage label differs. The Mathematica
output is internally consistent (banner reads STAGE 214) but shares the same
2026-06-02 mtime, predating no script edit of concern (the `.wl` mtime is
2026-06-02 11:41, older than its output, so the `.wl` output is fresh relative to
the `.wl`).

**Why this matters:**
Low. The content is correct; only the stale "197" label is wrong, and that is a
display artifact, not a math error. A fresh re-run will regenerate with the
correct "214" banner.

**Required change:**
Re-run `python3 scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py` and recommit the
captured output so the banner reads STAGE 214.

**Verification:**
After re-run, the saved `.txt` line 3 should read `STAGE 214 — ...` and the final
line `STAGE 214 SYMPY AUDIT COMPLETED SUCCESSFULLY`; all residuals remain `0` and
the sweep counts remain `924`.

## Independent-derivation check (Mathematica)

**Verdict: INDEPENDENT.**

The decisive discriminator is the projected eliminant family (`C_rs,C_rt,C_st` and
`S_r,S_s,S_t`), which is the load-bearing object behind the `150` bound. The two
engines extract these by genuinely different operations:

- **`.py` POSITS the closed-form elimination products** (`.py:158-164`):
  ```
  Crs = sp.expand(Ms * Lr - Mr * Ls)
  Sr  = sp.expand(Lr**2 - 4 * Mr**2 * Delta)
  ```
  i.e. it writes the pairwise-eliminated product and the squared condition by hand.

- **`.wl` DERIVES them by an actual resultant elimination of the auxiliary `y`**
  (`.wl:134-148`), then checks that derived resultant equals the posited closed
  form:
  ```
  crossEliminants = { -Cancel[Resultant[liftedPolys[[1]], liftedPolys[[2]], y]/2], ... };
  squareEliminants = Resultant[#, deltaLift, y] & /@ liftedPolys;
  ...
  expectZero["M4 C_rs resultant minus definition", crossEliminants[[1]] - crossDefinitions[[1]]];
  expectZero["M4 S_r resultant minus definition", squareEliminants[[1]] - squareDefinitions[[1]]];
  ```
  `Resultant[2 M_r y + L_r, 2 M_s y + L_s, y]` is a genuine elimination of `y`,
  not a re-typing of `M_s L_r − M_r L_s`; that the resultant reproduces the posited
  product is a non-trivial check the `.py` never performs. This is the
  derive-vs-posit split the policy requires.

A second independent route appears in the stationary numerators:

- **`.py` POSITS the monomial closed forms** (`.py:73-79`):
  `Mr = (1+r²+s²+t²)·kj − r·(ki+kj·r+kk·s+kl·t)`, `Lr = (1+r²+s²+t²)·∂_rΔ − 2r·Δ`.
- **`.wl` ASSEMBLES them from differentiated building blocks** (`.wl:91-98`):
  `metricNumerators = den·D[linearK,var] − D[den,var]·linearK/2`,
  `envelopeNumerators = den·D[quadEnvelope,var] − D[den,var]·quadEnvelope`.
  (`D[den,r]=2r`, so this equals `den·kj − r·linearK` — same object, different route:
  `.wl` differentiates `linearK`/`den`/`quadEnvelope` symbolically where `.py`
  expands the final monomial by hand.)

Both engines DO autodiff the same `Φ` for the derivative-law residual (`.py:106`
`sp.diff(Phi,r)` vs `.wl:101-104` `D[phi,var]`) and both count total degrees — but
those are the verification *comparison* and a trivial scaffolding count, not the
load-bearing extraction. For the load-bearing eliminants the `.wl` runs
`Resultant[...,y]` while `.py` posits the product; that single differing operation
is sufficient to clear the transliteration bar. Not a port.

## Engine cross-check

Both engines agree. SymPy: `54`, `150`, all derivative-law residuals `0`, all
eliminant degrees `(5,5,5,6,6,6)`, iso/symmetry residuals `0`. Mathematica:
identical `54`, `150`, `{3,3,3,2}`, `{5,5,5}`/`{6,6,6}`, all `M*` residuals `0` and
`PASS`. The single asymmetry is deliverable (7) (winner/non-improvement order
theorem), present only in SymPy (Section VI, 924 samples); Mathematica has no
counterpart M-section. Acceptable as single-engine for a trivial real-ordering
logical statement, but noted.

## Verdict justification

`findings: 1`. The math holds up under attack on both engines. Attacks tried that
failed: (a) checked A7 is not construction-guaranteed — the concrete coefficient
set is a real common root, so the vanishing is substantive; (b) checked the `.wl`
metric/envelope numerators are algebraically `Mr`/`Lr` but derived via `D[]` not
posited, and the eliminants are derived via `Resultant[...,y]` not retyped, so the
`.wl` is NOT a port of the `.py`; (c) confirmed degrees `(3,3,3,2)→54` and
`(5,5,6)→150` match paper + appendix + notes exactly; (d) confirmed the
diagonal-isotropic reduction `Δ_iso = linearK² − 2H₀ν·den` and gradient/equal-mix
rays match notes §5 and appendix eqs four-grad/four-eq. The only defect is the
stale SymPy `.txt` (stale "STAGE 197" banner from the pre-renumber epoch); content
is otherwise correct. Paper alignment confirmed: every Output deliverable maps to a
substantive script check.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 deliverable values checked, 0 misaligned.

The notes corrections flagged by the orchestrator HOLD:
- "218→150": the projected one-chart bound is `5·5·6 = 150` — notes line 435
  (`\boxed{5\cdot 5\cdot 6 = 150.}`), notes line 437 references it, and both scripts
  assert `150` (`.py:246` `bezout_projected != 150`, `.wl:172` `projectedProduct − 5*5*6`).
  No stray `218` survives in the notes (grep: no `218`).
- "230→162": no `162` (the Stage 217 five-coordinate bound) and no stray `230`
  appears anywhere in the Stage 214 notes (grep: no `162`, no `230`). The four-coord
  lifted bound is correctly `54`, not `162`. Confirmed corrected.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| lifted Bézout bound `54` | py:150 / wl:130; out `lifted Bezout bound = 54` (sympy out 122), `M3 = 54` (wl out 39) | `.tex:13,15` `(3,3,3)` /`54`; appendix `:1069`; notes `:317` `3·3·3·2=54` | MATCH |
| degree pattern `(3,3,3,2)` | py:139-146 / wl:120; wl out 26 `{3,3,3,2}` | `.tex:13` `(3,3,3,2)`; appendix `:1064`; notes `:307-313` | MATCH |
| projected one-chart bound `150` | py:246 / wl:172; sympy out 1160, wl out 76 | notes `:435` `5·5·6=150` | MATCH |
| eliminant degrees quintic `5` / sextic `6` | py:239-242 / wl:161-166; wl out 58-59 `{5,5,5}`/`{6,6,6}` | notes `:395` "quintic", `:415` "sextic" | MATCH |
| gradient-optimal ray `(k_i,k_j,k_k,k_l)/‖k‖` | py:281 / wl:201; wl out 87 | notes `:471-475`; appendix `:1042-1044` | MATCH |
| equal-mix barycenter `(1,1,1,1)/2` | py:320 / wl:225 | notes `:504-505`; appendix `:1048-1049` | MATCH |

INTERNAL (scaffolding, no finding): `Phi`,`tau`,`M_r/M_s/M_t`,`L_r/L_s/L_t`,
`N_·`, lifted polys `F_·`, the concrete stationary-point coefficients
(`A=-220/3,…,J=-85/3`, `(r,s,t,y)=(3/2,5/2,7/2,29/2)` — a verification probe, not a
reported result), the `924` integer-sweep sample counts, all `=0` residuals,
`zeroQ`/`PASS` flags, `derivativeDen`.

Note (documentation drift, not a value mismatch): the card's
`\stagefield{Verification}` (`.tex:11`) reads "Mathematica audit: none yet", yet a
`.wl` exists and passes (pass-1 dual-engine retrofit). This is a stale prose pointer
on the card, not a value reconciliation failure and not a script-side defect; flagged
here for the orchestrator's tracker-sync step but it does not change the verdict and
needs no Codex action (Codex does not edit paper/).

## Self-test notes

Variable-independence: checked `D[phi,var]`/`sp.diff(Phi,·)` genuinely depend on the
diff variable (Φ contains `r,s,t` in both numerator and `den`), so the derivative-law
residuals are non-trivially zero, not zero-because-constant. Trivial-case: A7's
concrete substitution reduces all 10 polynomials to literal `0` (verified the point
is a common root by the structure `y² = Δ` with `y=29/2 → 841/4` and the lifted
linear-in-`y` forms). Resultant route: confirmed `Resultant[2M_r y+L_r, 2M_s y+L_s, y]`
eliminates `y` to a multiple of `M_s L_r − M_r L_s` (the `/2` and sign in `.wl:135`
account for the leading-`y`-coefficient scaling), so the `.wl` check is correct, not a
port. No directive trap (only finding is stale_output + a paper-side prose note,
neither involving a new assertion).
