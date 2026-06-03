---
unit_id: 223
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 223 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_223.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows: line 58 status row; lines 574-597 detailed compatibility narrative)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 223 asks whether the same primitive finite-throat one-port branch can simultaneously satisfy the exact isotropic one-pole condition ($u_4=4u_2^2$) and the outgoing-normalization target ($P_0=P_{0,\rm target}$). The card's `\stagefield{Output}` states verbatim: "Exact primitive compatibility surface $N_0/P_{0,{\rm target}}=3(M+B_2+Z_2)^2/(B_4+Z_4)$ and finite dynamic survival windows on the compatible branch." The notes enumerate nine deliverables (notes §9): (1) the overlap constant $\kappa=2\sqrt2/\pi$; (2) the isotropic one-pole numerator identity $u_4-4u_2^2=[D_0(B_4+Z_4)-3(M+B_2+Z_2)^2]/D_0^2$; (3) the compatibility identity $N_0/P_{0,\rm target,compat}=3(M+B_2+Z_2)^2/(B_4+Z_4)$; (4) the primitive specialization of $P_{0,\rm target,compat}$; (5) the concrete sample values ($P_{0,\rm target,compat}\approx0.00207$, $K_{\rm compat}\approx24.47$, $D_{0,\rm compat}\approx24.24$); (6) the quartic pole census on the compatibility branch; (7) the four wall/internal residue-linewidth figures $\mathcal R_{Q,*}$; (8) the monotone $\lambda_W$ compatibility-family scan; (9) the four finite survival-window thresholds in $P_{0,\rm target,compat}$. The card carries a `Mixed` status (exact closure for formulas, numerical for the slice). The appendix (lines 574-597) restates the D-coefficient definitions and the boxed compatibility equation.

## What the script claims to verify

The SymPy script verifies the same nine items. It (a) computes $\kappa$ as an exact integral and asserts it equals $2\sqrt2/\pi$ (lines 92-93); (b) asserts the one-pole numerator identity vanishes via `cancel` (lines 149-152); (c) defines $P_{0,\rm target,compat}$ and $D_{0,\rm compat}$ then asserts the compatibility identity and a $D_0$ consistency identity vanish (lines 156-161), plus the primitive specialization (lines 163-166); (d) builds the exact quartic $F(y)$ and asserts degree 4 (lines 177-181); (e) evaluates the sample slice and asserts ~17 numeric literals via `assert_close` (lines 239-253); (f) computes the four-pole census and four $\mathcal R_{Q,*}$ values against literals (lines 263-286); (g) computes the two $\eta$-thresholds (lines 299-303); (h) runs the $\lambda_W$ compatibility-family scan against a 5-row `expected_scan` literal table with monotonicity asserts (lines 362-384); (i) bisects for the four survival-window thresholds (lines 389-417).

## Paper ↔ script cross-check

| Paper deliverable (notes §9) | Script-side check | Status |
|---|---|---|
| (1) $\kappa=2\sqrt2/\pi$ | lines 92-93 exact integral + assert | match |
| (2) one-pole numerator identity | lines 149-152 `cancel(...)==0` | match (but see F2: trivially true by construction) |
| (3) compatibility identity | lines 159-160 `cancel(...)==0` | match form, but tautological (F2) |
| (4) primitive specialization | lines 163-166 | match form, but tautological (F2) |
| (5) sample $P_{0,\rm compat},K_{\rm compat},D_{0,\rm compat}$ | lines 252-253 | match |
| (6) quartic pole census | lines 263-270 | match |
| (7) four $\mathcal R_{Q,*}$ | lines 279-286 | match |
| (8) monotone $\lambda_W$ scan table | lines 370-379 `expected_scan` | **mismatch at $\lambda_W=0.2$** (F1): notes give wall $\mathcal R_Q$ = 206.81 / 205.50; script asserts 138.81 / 137.50 |
| (9) four survival-window thresholds | lines 414-417 | match |

Set `paper_alignment: partial` — eight of nine deliverables match; deliverable (8) has a numeric mismatch at one scan node (F1).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 93 | `simplify(kappa - 2*sqrt2/pi) == 0` | claim 1 | yes |
| A2 | sympy | 152 | `cancel(u4-4u2^2 - num/D0^2) == 0` | claim 2 | partial (def-rewrite, see F2) |
| A3 | sympy | 160 | `cancel(N0/P0_compat - 3(...)^2/(B4+Z4)) == 0` | claim 3 | no (tautological, F2) |
| A4 | sympy | 161 | `cancel(D0_compat - N0/P0_compat) == 0` | claim 3 | no (tautological, F2) |
| A5 | sympy | 166 | `cancel(P0_compat - primitive_spec) == 0` | claim 4 | no (tautological def-inline, F2) |
| A6 | sympy | 181 | `Poly(F_y,y).degree() == 4` | claim 6 (setup) | yes |
| A7 | sympy | 239-253 | `assert_close` ×17 sample values | claim 5 | yes |
| A8 | sympy | 269-270 | `assert_close` 4 poles | claim 6 | yes |
| A9 | sympy | 285-286 | `assert_close` 4 $\mathcal R_{Q,*}$ | claim 7 | yes |
| A10 | sympy | 302-303 | `assert_close` 2 thresholds | claim 9 (setup) | yes |
| A11 | sympy | 377-379 | `assert_close` scan table | claim 8 | yes, but values conflict with notes (F1) |
| A12 | sympy | 381-384 | monotonicity asserts | claim 8 | yes |
| A13 | sympy | 414-417 | `assert_close` 4 windows | claim 9 | yes |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** value_mismatch
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md:351`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.py:371`

**What's wrong:**
The compatibility-family scan table disagrees at $\lambda_W=0.2$ in the two wall-like residue/linewidth columns.

Notes table (line 351, columns "lower wall $\mathcal R_Q$" / "upper wall $\mathcal R_Q$"):
> `| 0.2 | 0.000576970879843 | 29.3158464872314 | 206.814136942081 | 205.502546600713 |`

Script `expected_scan` (line 371) and saved output (line 63):
> `(0.000576970879843045, 29.3158464872314, 138.814136942081, 137.502546600713),`

The first two columns ($P_{0,\rm target,compat}=0.000576970879843$, $K_{\rm compat}=29.3158464872314$) agree exactly. Only the two $\mathcal R_Q$ columns differ: notes 206.814136942081 / 205.502546600713 vs. script 138.814136942081 / 137.502546600713. The fractional tails are byte-identical (`...814136942081`, `...502546600713`); only the integer parts differ by exactly 68 in both entries (206-68=138, 205-68=137). This is the signature of a hand-transcription typo on one side, not a genuine algebraic divergence. All other scan rows ($\lambda_W=0.4,0.6,0.8,1.0$) match between notes and script exactly, and the $\lambda_W=0.4$ row also matches the main-census $\mathcal R_{Q,*}$ values (30.20 / 36.17) computed independently earlier in the same script, so the script's scan formula is internally self-consistent.

**Why this matters:**
The notes table is the authoritative source for the paper card's deliverable (8) ("monotone compatibility-family scan"). If a reader transfers the 206.81/205.50 figures while the verified pipeline produces 138.81/137.50, the published narrative carries an unverified number. The monotonic-decrease claim holds under either value set (both 206 and 138 exceed the 0.4-row 30.20), so the qualitative conclusion is unaffected — but the specific figure is wrong on exactly one side and must be reconciled.

**Required change:**
Do NOT auto-resolve. Route to user. See `## Resolve before fix_loop` in the directive.

**Verification:**
After the user picks a direction: if the script is correct, the notes line 351 must read 138.814136942081 / 137.502546600713 (paper-side edit, Claude applies, not Codex); if the notes are correct, the script's `expected_scan[0]` at line 371 and the $\lambda_W=0.2$ scan output must read 206.814136942081 / 205.502546600713 and the script must reproduce them on re-run.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.py:156-161`

**What's wrong:**
The compatibility identity is asserted against a quantity defined to satisfy it. Line 156 defines
`P0_target_compat = N0*(B4+Z4)/(3*(M+B2+Z2)**2)`. Line 159 then forms
`compatibility_identity = cancel(N0/P0_target_compat - 3*(M+B2+Z2)**2/(B4+Z4))` and asserts it is zero (line 160). Substituting the line-156 definition, `N0/P0_target_compat = N0 / [N0*(B4+Z4)/(3*(M+B2+Z2)**2)] = 3*(M+B2+Z2)**2/(B4+Z4)` identically; the residual is zero by construction of the divisor. The assertion cannot fail for any choice of moments. The same applies to line 161: `D0_compat` is defined at line 157 as `K_pole - B0 - Z0`, and `K_pole` is defined at line 154 as `3*(M+B2+Z2)**2/(B4+Z4) + B0 + Z0`, so `D0_compat = 3*(M+B2+Z2)**2/(B4+Z4)`, which equals `N0/P0_target_compat` by the same definitional inversion — guaranteed zero. The primitive-specialization assert (line 166) is likewise a definition-inline: `P0_target_compat` and `primitive_specialization` are the identical expression with `N0=P**2/Delta**2`, `B4=C**2/varpi**6`, `B2=C**2/varpi**4` substituted, so `cancel(...)==0` only confirms that SymPy expands the same symbols the same way.

**Why this matters:**
These three asserts (A3, A4, A5) are the script's purported verification of paper deliverables (3) and (4) — the central compatibility surface that is the card's headline `\stagefield{Output}`. As written they verify nothing about the physics: they confirm that `a/(a*b/c) == c/b`, an algebra-of-fractions tautology. The genuine, non-tautological content for deliverable (3) is that $K_{\rm pole}=K_{\rm norm}$ implies $N_0/P_{0,\rm target}=3(M+B_2+Z_2)^2/(B_4+Z_4)$ — i.e. equating the two independently-defined wall stiffnesses ($K_{\rm pole}$ from the one-pole side, $K_{\rm norm}=N_0/P_{0,\rm target}+B_0+Z_0$ from the normalization side) and showing the result is the boxed equation. That equation IS independently exercised once via the one-pole numerator identity (line 152, which is a real algebra check on $u_4-4u_2^2$), so the claim is not entirely unverified — but the dedicated compatibility-surface asserts are circular.

**Required change:**
Replace the circular asserts with a check that equates the two independently-built wall stiffnesses. Concretely, the substantive identity is `cancel((K_pole - B0 - Z0) - N0/P0_target_compat) == 0` only AFTER `P0_target_compat` is obtained by *solving* `K_norm == K_pole` for `P0_target`, rather than by defining `P0_target_compat` as the target ratio up front. See directive F2 for the exact non-circular formulation (solve `K_norm - K_pole == 0` for `P0_target` with `sp.solve`, then assert the solution equals `N0*(B4+Z4)/(3*(M+B2+Z2)**2)`). This makes the assertion fail if the one-pole-side stiffness and the normalization-side stiffness did NOT actually coincide on the compatibility surface.

**Verification:**
After fix, a new `sp.solve(K_norm - K_pole, P0_target)` step should appear; the assert should compare its result to the boxed `P0_target_compat`. Removing the line-156 definition (or deriving it) and re-running must still exit 0. The one-pole numerator identity (line 152) remains untouched.

### F3 — missing_verification_script

**Severity:** medium
**Subtype:** missing_mathematica
**Files:**
- `(missing)` — no `.wl` exists for unit 223

**What's wrong:**
No Mathematica script exists for this unit (`find` returns only the `.py`, its `.txt`, and the notes `.md`; the `mathematica/` directory has no `stage223` file). The paper card line 11 explicitly states "Mathematica audit: none yet." This stage is `is_checkpoint: False` and `is_status_only_candidate: False`, so under the dual-engine rule a second engine is REQUIRED wherever Mathematica CAN independently verify — and Mathematica can: the κ integral, the one-pole numerator identity, the compatibility surface, the exact quartic, the rational moment algebra, and every numeric slice value are all native to Mathematica (`Integrate`, `Together`/`Cancel`, `Solve`/`NSolve`, `CoefficientList`). There is no genuine impossibility.

**Why this matters:**
Every load-bearing identity and the entire numeric census currently rest on a single engine (SymPy). The dual-engine policy exists precisely to catch transliteration-invisible algebra bugs and floating-point/root-selection artifacts (e.g. the `roots[-2:]` wall-pole selection at line 357, which is a hand-coded ordering convention that a native `NSolve` ordering could disagree with — exactly the kind of thing a second engine catches).

**Required change:**
Author an independent Mathematica audit at the exact target path below, deriving each claim from the physical premises with native primitives and a DIFFERENT decomposition (do not transliterate the SymPy choreography). See directive F3 for the claim manifest M1-M9 and the anti-transliteration guard.

**Verification:**
`redteam exec-mathematica 223` produces output; the new `.wl` exits 0 with all `expectZero`/`expect_close` checks passing, and its numeric slice values agree with the SymPy `.txt` (modulo the F1 resolution for the $\lambda_W=0.2$ row).

## Independent-derivation check (Mathematica)

No `.wl` exists, so no transliteration check applies yet. The future `.wl` (F3) must be an independent re-derivation, not a port.

## Engine cross-check

Only SymPy is present; no cross-check possible. This is itself finding F3.

## Verdict justification

Eight of nine paper deliverables are faithfully and (mostly) non-tautologically exercised by the SymPy script, and the saved output is fresh (output mtime 2026-05-11 12:50 > script mtime 11:56) and shows a clean PASS. The κ integral (A1), one-pole numerator identity (A2), quartic degree (A6), and the full numeric census including poles, residues, thresholds, the scan, and survival windows (A7-A13) hold up under attack — I checked the sample slice against the upstream Stage 222 slice (identical: $\lambda_B=1/2,\lambda_U=3/10,\lambda_W=2/5,\lambda_R=1/4,\Omega_U=1,\Omega_W=7/5,\varpi=2,M=1$), confirmed the scan formula is internally consistent (its $\lambda_W=0.4$ row reproduces the independently-computed main-census $\mathcal R_{Q,*}$), and verified the monotonicity claims survive both value sets. Three findings remain: (F1) a value_mismatch between the notes scan table and the script at $\lambda_W=0.2$ (206.81/205.50 vs 138.81/137.50, identical fractional tails — a transcription typo needing user-chosen direction); (F2) the compatibility-surface asserts (A3-A5) are circular definition-rewrites that cannot fail and so do not actually verify the card's headline Output; and (F3) the dual-engine rule requires a Mathematica audit that does not exist and is clearly possible. Verdict: `findings`, not stop_cold — F1 needs user resolution but does not mathematically propagate (the qualitative monotone-tradeoff conclusion is value-set-independent), and F2/F3 are script-side fixes Codex can apply after the user resolves F1's direction.

## Self-test notes

I checked: (1) variable-independence — the script contains no `sp.diff`; the residue $\mathcal R_{Q,*}$ has no derivative cancellation trap here. (2) Root selection — `roots[-2:]` (line 357) picks the two largest positive-$\omega^2$ roots as "wall-like"; verified this matches the main census labeling at $\lambda_W=0.4$ (1.998, 4.95 are indeed the largest two), so the scan's wall labels are consistent. (3) Tautology trace — confirmed A3/A4/A5 reduce to `a/(a*b/c)==c/b` and definition-inlining, guaranteed zero, hence F2; confirmed A2 (one-pole numerator) is NOT tautological because it rearranges $u_4-4u_2^2$ across distinct $D_0,D_2,D_4$ definitions and would expose a sign/factor error. (4) Numeric anchoring — all `assert_close` literals trace to notes §4/§5/§6 values; the only conflict is F1. (5) Paper round-trip for F2 fix — the proposed `sp.solve(K_norm-K_pole, P0_target)` reproduces exactly the boxed `P0_target_compat = N0(B4+Z4)/(3(M+B2+Z2)^2)`, introducing no new constant, so no new paper_misalignment.
