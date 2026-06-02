---
unit_id: 207
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table.md]
  paper_appendix: present
---

# Audit unit 207 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_207.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (subsection `Primitive certified table`, lines 783-831; status row line 45; narrative line 236)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

Stage 207 specializes the Stage 206 certified ray sieve to the five primitive axes `\(\{\lambda,c,\gamma,U,W\}\)` of the free-quintuple log space. The card's `\stagefield{Output}` reads verbatim: *"Exact primitive Hessian-envelope theorem, sign-adapted primitive drift table, certified primitive ray table, primitive elimination theorem, and primitive winner theorem."* The notes enumerate six deliverables: (1) the canonical orientation rule `\(\varepsilon_i=-\operatorname{sgn}(\Gamma_i)\)` giving oriented slope `\(K_i=-|\Gamma_i|<0\)`; (2) the diagonal-Hessian reduction `\(\mathcal H_{\widehat{\mathbf s}_i}=\partial_{\ell_i\ell_i}\)` (no off-diagonal entry needed for any primitive ray); (3) the certified monotone bracket via the Stage 240/206 root map `\(\mathcal T(H_0,k;c)=2H_0/(k+\sqrt{k^2-2cH_0})\)`, `\(k>0\)`; (4) the turning-row bracket `\(\tau^{\rm(tp)}=\sqrt{-2H_0/\kappa}\)` for `\(\overline\kappa_i<0\)`; (5) the sign-adapted primitive drift table (Section 4) giving the four dependent graph exponents `\(\sigma_\delta,\sigma_T,\sigma_{K_\eta},\sigma_\mu\)` per primitive axis, scaled by orientation sign `\(\varepsilon_i\)`; and (6) the statement that off-diagonal Hessian entries first appear only on genuine two-coordinate mixed rays `\(\mathcal H_{a\mathbf e_i+b\mathbf e_j}=a^2\partial_{ii}+2ab\partial_{ij}+b^2\partial_{jj}\)`. The elimination and winner theorems (Sections 6-7) are inequality/ordering statements built on the bracket data, not new closed-form identities.

## What the script claims to verify

The SymPy script (`STAGE 190` banner notwithstanding — see verdict note) verifies, with substantive `expect_zero` assertions, exactly the closed-form deliverables of the stage: Section I contracts a symbolic symmetric 5x5 Hessian against each `\(\pm\mathbf e_i\)` and asserts the result equals the diagonal entry (M1); Section II expands the two-coordinate quadratic form and asserts the off-diagonal coefficient is `\(2ab H_{ij}\)` (M2/M6 preview); Section III asserts `\(\varepsilon\Gamma+|\Gamma|=0\)` (M3); Section IV constructs generic `\(\sigma_\delta,\sigma_T,\sigma_{K_\eta},\sigma_\mu\)` expressions, evaluates them on each oriented primitive vector, and asserts equality with the paper Section-4 tabulated rows (M-rows); Section V builds the monotone root-map brackets `\(\tau_{\rm lo},\tau_{\rm hi}\)` and asserts each satisfies the comparison quadratic `\(H_0-k\tau+\tfrac12 c\tau^2=0\)` (M4 — the load-bearing check); Section VI does the analogous turning-row check `\(H_0-\tfrac12\kappa\tau^2=0\)` (M5); Section VII prints the five per-axis brackets (no assertion — summary only).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) canonical orientation `\(K_i=-|\Gamma_i|\)` | III: `expect_zero("canonical orientation law", eps*Gam + Abs(Gam))` | match |
| (2) diagonal Hessian reduction `\(\mathcal H_{\widehat e_i}=\partial_{ii}\)` | I: `\(\pm e_i^T H e_i - H_{ii}=0\)` for all 5 axes | match |
| (3) monotone certified bracket / root map | V: `tau_lo,tau_hi` satisfy `\(H_0-k\tau+\tfrac12 c\tau^2=0\)`; VII reprints per-axis | match |
| (4) turning certified bracket | VI: `tau^{(tp)}` satisfies `\(H_0-\tfrac12\kappa\tau^2=0\)` | match |
| (5) sign-adapted primitive drift table (Section 4) | IV: per-axis `\(\sigma_*\)` rows equal `\(\varepsilon_i\times\)` tabulated | match |
| (6) off-diagonal first appears on mixed rays | II: `\((a e_i+b e_j)^T H(\cdot)=a^2H_{ii}+2abH_{ij}+b^2H_{jj}\)` | match |
| primitive elimination theorem (Section 6) | none (inequality/discriminant logic, not a closed-form identity) | not-required |
| primitive winner theorem (Section 7) | none (ordering rule on bracket data) | not-required |

Convention note: the notes (lines 240, 261-263) and the script both use the `\(+k\)` form `\(\mathcal T=2H_0/(k+\sqrt{k^2-2cH_0})\)` with `\(k=|\Gamma_i|>0\)`. The appendix (lines 827-829) writes `\(\mathcal T(H_0,-k_i;\cdot)\)`, which is the identical map under the first equality `\(-2H_0/(-k-\sqrt{\cdot})=2H_0/(k+\sqrt{\cdot})\)` shown verbatim in notes line 240. No sign/branch misalignment; the three sources agree. The oriented slope in the comparison quadratic enters as `\(-k\)` (negative, matching `\(K_i<0\)`), and the script's Section V uses exactly `H0 - k*TauL + ...`, so the sign is correct.

The elimination/winner theorems are pure inequality and lexicographic-ordering rules over the already-certified `\((\tau_{\rm lo},\tau_{\rm hi},W_i,T_i)\)` data; they introduce no new closed-form identity to verify symbolically, so the absence of a dedicated assertion is not a `script_missing_paper_claim` defect. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 66-67 | `expect_zero(±e_i^T H e_i - H_ii)` (×5 axes) | (2) diagonal Hessian reduction | yes |
| A2 | sympy | 81 | `expect_zero(mixed_expr - (a²Hll+2abHlc+b²Hcc))` | (6) mixed-ray quadratic form | yes |
| A3 | sympy | 82-85 | `expect_zero(mixed - (a²Hll+b²Hcc) - 2abHlc)` | (6) off-diagonal first appearance | yes |
| A4 | sympy | 99 | `expect_zero(eps*Gam + Abs(Gam))` | (1) canonical orientation | yes |
| A5 | sympy | 169 | `expect_zero(row - eps_i*expected)` (×5 axes) | (5) sign-adapted drift table | yes |
| A6 | sympy | 183-186 | `expect_zero(H0 - k*TauL + ½cL*TauL²)` | (3) monotone bracket root map | yes |
| A7 | sympy | 187-190 | `expect_zero(H0 - k*TauU + ½cU*TauU²)` | (3) monotone bracket root map | yes |
| A8 | sympy | 207-210 | `expect_zero(H0 - ½a_turn*TauTurnLo²)` | (4) turning bracket | yes |
| A9 | sympy | 211-214 | `expect_zero(H0 - ½b_turn*TauTurnHi²)` | (4) turning bracket | yes |
| (V VII) | sympy | 238-242 | print only, no assertion | (3) per-axis summary | n/a (display) |

All assertions are non-tautological: in I/II the quadratic form is contracted from a fully symbolic symmetric matrix and compared to an independently-written target; in V/VI the closed-form root is substituted back into a separately-written defining quadratic (a genuine root-verification, not a definition-echo); in IV the `\(\sigma_*\)` rows are generic functions evaluated then matched against the paper table. No row is "no"-anchored. Section VII is display-only and carries no claim.

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` — no `.wl` for stage 207 (confirmed: directory holds 197 `.wl` files following `moving_throat_pde_stage<NNN>_<slug>_mathematica_audit.wl`; none matches stage 207).
- paper card line 11 acknowledges: `\stagefield{Verification}{... Mathematica audit: none yet.}`

**What's wrong:**
Subtype `missing_mathematica`. Stage 207 is non-status-only (`is_status_only_candidate: False`) and non-checkpoint (`is_checkpoint: False`). It computes multiple independently-verifiable closed-form results — the diagonal Hessian contraction, the mixed-ray quadratic form, the orientation-law identity, the sign-adapted exponent table, and the two certified root-map brackets. Every one of these is a symbolic identity that Mathematica can re-derive natively (matrix contraction, `Coefficient`/`Series`, `Solve`/`Reduce` on the comparison quadratic, `Sign`). Per the project dual-engine contract and the rendered prompt (line 118: "both scripts are required, and missing scripts are findings"), the absence of a second engine is a finding. The test is *can Mathematica independently verify*, and here it plainly can — the gap is not excused by the SymPy script passing.

**Why this matters:**
A single-engine verification has no cross-check against SymPy-specific simplification quirks or a transcription error in the closed-form brackets (Section V/VI), which are the load-bearing physics outputs. The two engines must agree via independent derivations.

**Required change:**
Add `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_mathematica_audit.wl` verifying the claim manifest M1–M6 (see directive). It must derive each result independently using native Mathematica primitives via a different decomposition than the `.py` — a line-by-line port is rejected as transliteration.

**Verification:**
`redteam exec-mathematica 207` produces a `.wl` output file; the new script exits 0 with all in-file checks (M1–M6) passing, and its residuals agree with the SymPy output for the shared identities.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot yet be assessed. The directive's anti-transliteration guard requires the new `.wl` to use a *different decomposition* (e.g. `Solve`/`Reduce` to obtain the bracket roots rather than constructing the closed form and substituting back, as the `.py` does; native symmetric-matrix contraction and `Coefficient` extraction for the Hessian forms).

## Engine cross-check

Only SymPy present; cross-check deferred until the `.wl` lands. The SymPy output (mtime 2026-05-11 12:49:16, newer than the script's 2026-05-11 11:56:51) is fresh, exit code 0, status PASS; all `expect_zero` lines report `= 0`.

## Verdict justification

Verdict `findings`. The SymPy script is substantive and faithfully aligned with the paper: all six closed-form deliverables map to non-tautological assertions, the root-map sign/branch convention agrees across notes, appendix, and script (the appendix's `\(-k_i\)` argument is the same map under the equality stated verbatim in the notes), and the Section-4 drift table is reproduced exactly. Attacks tried and failed: (a) hunted for a sign flip in the comparison quadratic — the oriented slope correctly enters as `\(-k\)` matching `\(K_i<0\)`; (b) checked whether Section V's root-substitution is a definitional tautology — it is not, the root is plugged into an independently-written quadratic; (c) checked the mixed-ray off-diagonal coefficient for a factor-of-2 error — it is `\(2ab\)`, correct; (d) checked the `\(\Gamma=0\)` branch of the orientation law — SymPy's `ConditionalExpression` handles it and the residual is 0. The only defect is the missing second engine (F1), which is the standing dual-engine requirement, not a math error. Not `UNFIXABLE` (math is sound) and not `CRITICAL_DOWNSTREAM` (adding a `.wl` changes no derived constant or forward). Cosmetic note (not a categorized finding, no math impact): the banners at lines 35/244 read `STAGE 190` rather than 207, and the Section-IV comment line 104 says "from Stage 204"; these are stale display labels carried from an earlier stage template and do not affect any verified identity.

## Self-test notes

Ran the required traps on the M1–M6 manifest. (1) Variable independence: no `D[]`/`diff` traps — M1/M2/M6 are pure matrix contractions, M3 is a `Sign` identity; M4/M5 verify a substituted root against a polynomial, no zero-derivative pitfall. (2) Symmetry/parity: no unbounded integrals in this stage; n/a. (3) Trivial-case pre-checks pass: M4 with `cL→0` gives `tau_lo=H0/k` and `H0-k(H0/k)+0=0`; M5 gives `H0-(1/2)a·(2H0/a)=0`; M1 gives `e_i^T H e_i - H_ii=0` identically. (4) Path: directive names the full target under `mathematica/`. (5) Paper round-trip: the manifest carries the notes/appendix-stated constants and the `\(+k\)` root-map form, introducing no new `paper_misalignment`.
