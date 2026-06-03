---
unit_id: 246
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 246 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_246.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (rows: line 90 status row; lines 224-238 block-contract source-family/moments equations; line 117 compiler-chain checklist)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

Stage 246 compiles the compensated source lane declared abstractly in Stage 243. It claims (card eqs. + notes Sec. 1-9, appendix row line 90, block-contract lines 224-238): (1) the mean-preserving two-mode family `σ_{a,b}(x)=1+a cos(πx)+b cos(2πx)` integrates to 1 on `[0,1]`, and rewrites as the exact quadratic `σ=1-b+ay+2by²` under `y=cos(πx)`; (2) the **exact minimum** `σ_min(a,b)` is the stated piecewise (boundary case `1+b-|a|` for `b≤0` or `|a|≥4b`; interior-vertex case `1-b-a²/(8b)` otherwise), giving the sign-change test `signchg ⇔ σ_min<0`; (3) the mouth-bias functional `g(a,b)=(2/π)(1+a/3-b/15)`; (4) the shell-loading functional `S(a,b)=(2 tanh(π/2)/π)(1+a/5+b/17)`; (5) the normalized two-moment map with `det M_src = 14/425>0` and the explicit inverse `a=(85/42)S̃+(25/14)g̃`, `b=(425/42)S̃-(85/14)g̃`; (6) the mixed-to-shell ratio `R=(g-r_F1)²/(1+r_F1²)` with `r_F1=√(12/π²·(37/20)²-1)`, the compensation line `b=5a+15-(15π/2)g_c`, and the **quarter-ratio theorem** `R=1/4` at `g_±=r_F1±½√(1+r_F1²)`; (7) the transported closure `a(r)=a0 s(r)`, `b(r)=b0 s(r)`, `s(r)=r_σ²/(r²+r_σ²)`; (8) on the Session-I orientation `a0>0,b0<0`, `σ_min(r)=1-(a0-b0)s(r)` with threshold `signchg ⇔ r<r_σ√(a0-b0-1)`; (9) the Session-I readback `g≈0.82823667`, `S≈0.67584771`, `R≈0.21677037`, `σ_min≈-0.08979545<0`, `g(0)≈1.1289>1`. `\stagefield{Verification}` reads: "SymPy audit: ... Mathematica audit: none yet." The appendix-row status is `\StatusExactClosure{}`.

## What the script claims to verify

The docstring enumerates 9 verified objects matching the 9 paper sections. The genuinely load-bearing checks are: the family integrates to 1 (L53); the vertex value formula (L55); the two moment integrals `g` and `S` against the closed forms (L101-102); `det M_src=14/425` (L125) and the inverse round-trip `M_src·[a_inv,b_inv]=[g̃,S̃]` (L126); the quarter-ratio `R(g_∓)=1/4` and the compensation line (L154-156); the transported `σ_min(r)` algebraic rearrangement (L196); and the numeric Session-I readback for `g, R, σ_min, r_thr, g(0)` (L230-234). Sections 6 (Π packet) and 7 (transported `g(r)/S(r)/R(r)` forms) are print-only.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1a) mean = 1 | L53 `sigma_mean==1` (real integral) | match |
| (1b) quadratic rewrite under y=cos(πx) | L54 asserts `expand(E)-E==0` on a self-defined `sigma_y`; substitution link never tested | mismatch (tautological) |
| (2) σ_min = exact minimum / sign-change law | L62-83 Piecewise self-tested vs identical Python branch logic; never tied to actual min of σ | partial/mismatch (tautological + insufficient) |
| (3) g(a,b) | L101 integral vs closed form | match |
| (4) S(a,b) | L102 integral vs closed form | match |
| (5) det=14/425 + inverse | L125 det; L126 forward·inverse round-trip | match |
| (6) R, compensation line, quarter-ratio | L154-156 | match |
| (7) transported closure | print-only (L177-184); exercised numerically in Sec 9 | partial |
| (8) Session-I σ_min(r) + threshold | L196 algebraic rearrangement of `1+b_r-a_r`; boundary-branch selection for a0>0,b0<0 not verified | partial (tautological) |
| (9) Session-I readback g/R/σ_min/g(0) | L230-234 numeric asserts; S readback (≈0.67584771) printed but not asserted | partial (S unasserted) |

Dominant pattern: the substantive functional/algebra claims (g, S, det, inverse, R, quarter-ratio, numeric readback) are faithfully and non-tautologically verified; the *minimum/sign-change* claims (1b, 2, 8) are encoded but not independently established. `paper_alignment: aligned` (no value/target disagreement with the published card — all formulas and constants match; the gaps are verification-strength gaps, not paper conflicts).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `sigma_mean == 1` | claim 1a | yes |
| A2 | sympy | 54 | `simplify(sigma_y - (1-b+ay+2by²))==0` | claim 1b | no (tautological) |
| A3 | sympy | 55 | `simplify(sigma_vertex - (1-b-a²/8b))==0` | claim 2 (vertex) | partial (on unproven quadratic) |
| A4 | sympy | 56 | `simplify(sigma(a=0,b=0)-1)==0` | claim 1 (recovery) | yes |
| A5 | sympy | 83 (×3) | `simplify(piece.subs - expected)==0` | claim 2 (σ_min) | no (tautological self-test) |
| A6 | sympy | 101 | `simplify(g_sigma - g_expected)==0` | claim 3 | yes |
| A7 | sympy | 102 | `simplify(S_sigma - S_expected)==0` | claim 4 | yes |
| A8 | sympy | 125 | `det_M_src == 14/425` | claim 5 | yes |
| A9 | sympy | 126 | `M_src*[a_inv,b_inv]-vec == 0` | claim 5 (inverse) | yes |
| A10 | sympy | 154 | `R_minus - 1/4 == 0` | claim 6 (quarter-ratio) | yes |
| A11 | sympy | 155 | `R_plus - 1/4 == 0` | claim 6 | yes |
| A12 | sympy | 156 | `b_from_gc - (5a+15-15πg_c/2)==0` | claim 6 (comp. line) | yes |
| A13 | sympy | 196 | `sigma_min_transport - (1-(a0-b0)s_r)==0` | claim 8 | no (tautological) |
| A14 | sympy | 230 | `abs(g_eval - 0.82823667) < 5e-9` | claim 9 | yes |
| A15 | sympy | 231 | `abs(R_eval - 0.21677037) < 5e-9` | claim 9 | yes |
| A16 | sympy | 232 | `abs(sigma_min_eval - (-0.08979545)) < 5e-9` | claim 9 | yes |
| A17 | sympy | 233 | `r_thr_eval > r_eval` | claim 9 (window) | yes |
| A18 | sympy | 234 | `g_zero_eval > 1.0` | claim 9 (g(0)>1) | yes |
| (none) | sympy | — | S_eval readback (0.6758...) | claim 9 (S) | missing assert |

## Findings

### F1 — missing_verification_script

**Severity:** high
**Files:**
- `(missing)` — no `.wl` exists for unit 246; manifest `mathematica.path: null`, `exists: false`.

**What's wrong:**
Stage 246 is SymPy-only. The unit is `is_checkpoint: false` but also `is_status_only_candidate: false`, so under the dual-engine rule both engines are required wherever Mathematica *can* independently verify. Every claim here is squarely verifiable in native Mathematica with a *different* decomposition than the SymPy script:
- `Integrate[sigma*Cos[Pi x/2], {x,0,1}]` and the `K_q` integral independently reproduce `g`/`S` (no reuse of the `.py`'s product-to-sum choreography).
- The σ_min/minimum claim can be verified the way the SymPy script does NOT: `Minimize[{sigma, 0<=x<=1}, x]` / `MinValue` directly against the piecewise, closing the tautology in F2/F3 from the second engine.
- `Det`/`Inverse` of `M_src`, the quarter-ratio algebra, the compensation line via `Solve`, the transported closure, the Session-I orientation branch via `Reduce`, and the numeric readback are all native primitives.

**Why this matters:**
No independent cross-check exists on a non-status, exact-closure stage. The σ_min and quarter-ratio claims are load-bearing for downstream Stage 247 (compiler chain, appendix line 117) and currently rest on a single engine — and, per F2/F3, on tautological encodings within that single engine.

**Required change:**
Add a NEW independent-route Mathematica audit using native primitives (claim manifest in the directive). Filename must be `moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_mathematica_audit.wl` in `mathematica/`.

**Verification:**
`redteam exec-mathematica 246` shows the new `.wl` present, exits 0, with `expectZero`/`expectTrue`/`expectApprox` covering M1-M9 in the directive.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py:40,54` (quadratic rewrite)
- `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py:62-83` (σ_min self-test)

**What's wrong:**
Two checks cannot fail by construction:

1. Line 40 defines `sigma_y = sp.expand(1 - b + a*y + 2*b*y**2)`, and line 54 asserts `simplify(sigma_y - (1 - b + a*y + 2*b*y**2)) == 0` — i.e. `expand(E) - E == 0`, guaranteed true. The paper's claim 1b is that *σ(x) equals this quadratic under the substitution `y=cos(πx)`* (notes lines 168-180: "Then `cos(2πx)=2y²-1`, so the compensated source becomes the exact quadratic"). That substitution link is never exercised; the SymPy `cos(πx)`/`cos(2πx)` form is never connected to the `y`-quadratic.

2. Lines 62-65 define `sigma_min_piece` as a SymPy `Piecewise` with conditions `Or(Le(b,0),Ge(Abs(a),4b))`. Lines 76-83 then recompute `sigma_min_expected` with a Python `if (bval<=0) or (abs(aval)>=4*bval)` using the *same* formulas, and assert the Piecewise substitution equals it. Both sides encode the identical piecewise definition; the assertion only confirms SymPy's `Piecewise` dispatch agrees with Python's `if` dispatch on the same conditions and formulas. It can never reveal whether the piecewise actually *is* the minimum of σ over the domain.

**Why this matters:**
The "exact sign-change test" (deliverable 1 in the notes' Status list, appendix line 90) is the headline claim of the stage, and it is currently unverified — the script asserts a definition against itself. A wrong `σ_min` formula (e.g. wrong vertex coefficient, wrong branch boundary `4b`) would pass both checks.

**Required change:**
- For the quadratic rewrite: substitute `y -> cos(pi*x)` into the asserted quadratic and assert it equals the original `sigma` symbolically, i.e. `assert sp.simplify(sigma - (1 - b + a*sp.cos(sp.pi*x) + 2*b*sp.cos(sp.pi*x)**2)) == 0` (uses `cos(2πx)=2cos²(πx)-1`). This is the non-tautological version of the rewrite. Keep line 54 or replace it.
- For σ_min: add a real minimization check tying the piecewise to the actual minimum of `sigma_y` on `y∈[-1,1]`. At each of the three existing test points, compute the true minimum as `min(boundary(+1), boundary(-1), vertex if -1<=y_star<=1 and 2b>0)` from the quadratic and assert `sigma_min_piece.subs == ` that true minimum. (See directive F2 for the concrete, trap-checked construction.)

**Verification:**
New asserts appear near lines 54 and 83; substituting a deliberately wrong vertex/boundary formula would now fail. Output Section 1/2 unchanged in printed values.

### F3 — tautological_check

**Severity:** medium
**Files:**
- `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py:190,196`

**What's wrong:**
Line 190 defines `sigma_min_transport = sp.simplify(1 + b_r - a_r)` with `a_r=a0*s_r`, `b_r=b0*s_r`. Line 196 asserts `simplify(sigma_min_transport - (1 - (a0-b0)*s_r)) == 0`, which is pure algebraic rearrangement of the definition (`1 + b0·s - a0·s ≡ 1 - (a0-b0)s`) and cannot fail. The paper's claim 8 (card lines 114-120, notes lines 443-464) is stronger: *on the Session-I orientation `a0>0,b0<0` the piecewise selects the boundary-minimum branch*, so `σ_min(r) = 1 + b(r) - |a(r)| = 1 + b_r - a_r`. The script hardcodes `1 + b_r - a_r` without verifying (i) that `a_r>0` makes `|a_r|=a_r`, and (ii) that the orientation `a0>0,b0<0` lands in the `b≤0` branch of `sigma_min_piece`. The connection to `sigma_min_piece` from Section 2 is never made.

**Why this matters:**
The transported sign-change threshold `r<r_σ√(a0-b0-1)` (asserted only indirectly via the numeric `r_thr_eval>r_eval` at line 233) depends on the boundary-branch claim. If the orientation actually fell in the interior-vertex branch, the threshold formula would be wrong, and the current check would not catch it.

**Required change:**
Tie `sigma_min_transport` to `sigma_min_piece`: declare `a0` positive and `b0` negative (or substitute the Session-I signs symbolically with `a0>0, b0<0`), evaluate `sigma_min_piece.subs({a:a_r, b:b_r})` under those assumptions, and assert it equals `1 - (a0-b0)*s_r`. This exercises the branch selection rather than restating the boundary formula. (See directive F3.)

**Verification:**
New assert near line 196 references `sigma_min_piece`; a wrong branch boundary would now surface. Output Section 8 printed values unchanged.

### F4 — insufficient_verification

**Severity:** low
**Files:**
- `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py:213,223` (S readback computed/printed, no assert)

**What's wrong:**
Deliverable 9 (card lines 510-514, notes lines 510-514) lists `S[σ](r_eval)≈0.675847711465632` as a recorded Session-I readback. The script computes `S_eval` (line 213) and prints it (line 223, output shows `0.6758477114656324`) but never asserts it against the recorded value, unlike `g`, `R`, and `σ_min` which all have numeric asserts (L230-232). The symbolic `S(a,b)` form is verified (A7), so this is only a missing numeric pin, not a paper conflict.

**Why this matters:**
A regression in the transported `S(r)` form (e.g. a sign error in the `b0/17` term) that left `g` and `R` intact would go uncaught at the readback gate.

**Required change:**
Add `assert abs(float(S_eval) - 0.67584771) < 5e-9` after line 232.

**Verification:**
New assert appears after line 232; output Section 9 `S[sigma](r_eval)` line unchanged.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration is not at issue. The required new `.wl` (F1) must NOT mirror the `.py`'s decomposition: it must use `Integrate` directly for the moments (not product-to-sum by hand), `Minimize`/`MinValue` for σ_min (closing the F2/F3 tautologies independently), and `Det`/`Inverse`/`Solve`/`Reduce` for the algebra. The directive specifies this explicitly.

## Engine cross-check

Only one engine present; `engines_agree: n/a`. The single engine's substantive checks (g, S, det, inverse, quarter-ratio, numeric readback) are internally consistent and reproduce the paper's published numerics to the stated precision (output lines 87-92).

## Verdict justification

`findings: 4`, `stop_cold: null`. The stage is well-aligned with the published card — every formula and constant in the script matches the paper (g, S, det=14/425, inverse coefficients, r_F1=1.77799…, quarter-ratio, transported closure, Session-I numerics all confirmed). No `paper_misalignment`. What does NOT hold up: the headline "exact sign-change test" and the quadratic rewrite are asserted against their own definitions (F2), and the transported σ_min is a definitional rearrangement that skips branch selection (F3); these are real tautologies that would pass a wrong minimum formula. Plus the stage is single-engine on a non-status exact-closure unit (F1) and the S readback is unpinned (F4). Attacks that failed: I tried to break the `g`/`S` integrals by hand (both reproduce `(2/π)(1+a/3-b/15)` and `(2tanh(π/2)/π)(1+a/5+b/17)` exactly), the det (14/425 confirmed), the inverse round-trip (coefficients confirmed), and the quarter-ratio (`(g_∓-r_F1)²=(1+r_F1²)/4` is a genuine identity) — all held. I confirmed I read the paper card, the notes, and the appendix block-contract before opening the script, and the script's verified formulas match the paper; the gaps are verification strength, not paper conflict.

## Self-test notes

Walked the prescribed fixes as SymPy/Mathematica: (1) Variable independence — the proposed σ_min minimization fix operates on `sigma_y(y)` which depends on `y,a,b`; no `diff` w.r.t. an absent variable is introduced (I avoided the independence-self-test trap by comparing min-values, not differentiating). (2) Trivial-case pre-check — at `(a,b)=(1/2,-1/5)` the true min is `min(1+a+b, 1-a+b)=min(1.3, 0.3)=0.3=3/10`, matching output; at `(5/2,1/4)` boundary `min(1+a+b,1-a+b)=min(3.75,-0.75)=-1.25=-5/4` and `|a|=2.5≥4b=1` so boundary branch, matching `-5/4`; at `(1/2,1/2)` interior `1-b-a²/8b=1-.5-.0625=.4375=7/16`, matching — so the fixed asserts reduce to the printed values. (3) Branch selection for F3 — with `a0>0,b0<0`, `b_r=b0 s<0` so the `Le(b,0)` condition holds and `Abs(a_r)=a_r` (since `a0>0,s>0`), giving `1+b_r-a_r=1-(a0-b0)s`, consistent. (4) Path — `.wl` target placed in `mathematica/` with the exact mandated stem.
