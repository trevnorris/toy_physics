---
unit_id: 021
batch: I.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T20:40:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage021_reduced_one_port_normal_form.md]
  paper_appendix: present
---

# Audit unit 021 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_021.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage021_reduced_one_port_normal_form.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row at line 64; group sentence line 9)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt`

## What the paper claims

The stage card (`\stagefield{Output}`, lines 77–81) states: "Stage~021 exports the reduced one-port self-energy \eqref{eq:app-stage021-self-energy}, the transfer factor \eqref{eq:app-stage021-transfer-factor}, and the wall-level outgoing quadrupole coefficient \eqref{eq:app-stage021-wall-odd}." Concretely the deliverables are: (1) the conservative mixed Maxwell self-energy `Σ = (g_A² W + 2 g_A g_W R + g_W² A)/(A W − R²)`; (2) the gauge invariance of the mixed fields `E_w`, `C_a`; (3) the low-frequency conservative coefficients `z0, z2, z4`; (4) the first-order outgoing transfer factor `N_l(ω) = [A g_W + R g_A]²/[A W − R²]²` with `N_l(0) ≥ 0`; (5) the compact outgoing `l=2` fingerprint `Ŷ₂ = 1 + a²ω²/(9c_s²) + 4a⁴ω⁴/(81c_s⁴) + i a⁵ω⁵/(27c_s⁵) + O(ω⁶)`; (6) the wall-level odd quadrupole coefficient `δD₂^odd = −i N₂(0) (a⁵/27c_s⁵) ω⁵ + O(ω⁷)`. The notes add a seventh deliverable absent from the card body but stated in the notes §8: (7) the scalar derivative-coupling compatibility check showing a derivative-coupled scalar outlet (`g_{W,0}=ηω`) pushes the dangerous odd scalar term from `iω` up to `iω³` at wall level. The appendix row (line 64) summarizes: "Retained reduced one-port self-energy, transfer factor, compact l=2 fingerprint, and scalar-compatibility criterion."

## What the script claims to verify

Both engines verify, in five sections: (I) `E_w` and `C_a` are gauge-invariant (variation = 0 under the stated gauge transformation); (II) the EL equations of the reduced Lagrangian, the closed-form conservative self-energy `Σ_cons` (sympy via direct A/W solve, Mathematica additionally via `LinearSolve` plus an explicit cross-check `sigmaCons from LinearSolve matches closed form`), and the three low-frequency coefficients `z0, z2, z4` both as a generic toy-rational series and as the substituted physical expressions; (III) the first-order outgoing correction giving the compact transfer factor `N(ω)=(A g_W + R g_A)²/Δ²`, its `N(0)` positive-square form, and the composed wall-odd correction `δD = −i Γ ω⁵ N(0)` with `Γ = a⁵/(27c_s⁵)`; (IV) the `l=2` Hankel fingerprint `Ŷ₂` and the extracted port coefficient `Γ₅ = a⁵/(27c_s⁵)`; (V) the derivative-coupled scalar lane giving `N_scalar ~ η²Ω_A⁴ω²/Δ₀²` and `Π₀·N_scalar = iγ₁η²Ω_A⁴ω³/Δ₀²`. Every section ends in `expect_zero` / `expectZero` residual checks; all are reported = 0 in both saved outputs.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) Σ_cons closed form (eq self-energy) | `Sigma z0/z2/z4` substitutions + `A/W exact solution residual` (py 148–149; wl 91,95–96) | match |
| (2) gauge invariance E_w, C_a | `E_w gauge variation`, `C_a gauge variation` (py 88–89; wl 55–56) | match |
| (3) z0,z2,z4 coeffs (eq z-coeffs) | `z0/z2/z4 formula` + `Sigma z0/z2/z4` (py 161–177; wl 104–119) | match |
| (4) transfer factor N(ω), N(0)≥0 (eq transfer-factor, n0-positive) | `N(omega) compact formula`, `N(0) positive-square form` (py 222–231; wl 146–147) | match |
| (5) Ŷ₂ fingerprint (eq outgoing-fingerprint) | `Y2_hat minimal branch`, `Gamma5_port` (py 290–302; wl 167–171) | match |
| (6) δD₂^odd wall-odd coeff (eq wall-odd) | `delta D_2^(odd) composed …` (py 241–251; wl 136–140) | match |
| (7) scalar derivative-coupling → iω³ (notes §8) | `N_scalar leading term`, `scalar odd order` (py 335–347; wl 191–192) | match |

`paper_alignment: aligned` — every paper-side deliverable, including the notes-only scalar criterion (which the appendix row also names "scalar-compatibility criterion"), has a faithful non-tautological script-side check, and the symbolic forms in the saved outputs match the card equations term-for-term.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 88 | `expect_zero(Ewp - Ew)` | claim 2 | yes |
| A2 | sympy | 89 | `expect_zero(Cap - Ca)` | claim 2 | yes |
| A3 | sympy | 128–130 | EL residuals = 0 | claim 1 (setup) | yes |
| A4 | sympy | 148–149 | A/W solution residual = 0 | claim 1 | yes |
| A5 | sympy | 161–163 | z0/z2/z4 toy series = formula | claim 3 | yes |
| A6 | sympy | 175–177 | Σ series z0/z2/z4 = substituted formula | claim 3 | yes |
| A7 | sympy | 222 | `N(ω) − (A g_W+R g_A)²/Δ²` = 0 | claim 4 | yes |
| A8 | sympy | 228 | `N(0) − (Ω_A² g_W+R g_A)²/(Ω_A²Ω_W²−R²)²` = 0 | claim 4 | yes |
| A9 | sympy | 241 | δD composed from N(0)·Γ₅ = 0 | claim 6 | yes |
| A10 | sympy | 290 | `Ŷ₂ − minimal branch` = 0 | claim 5 | yes |
| A11 | sympy | 302 | `Γ₅ − a⁵/(27c_s⁵)` = 0 | claim 5 | yes |
| A12 | sympy | 335 | `N_scalar − η²Ω_A⁴ω²/Δ₀²` = 0 | claim 7 | yes |
| A13 | sympy | 344 | scalar odd order = iγ₁η²Ω_A⁴ω³/Δ₀² | claim 7 | yes |
| B1 | mathematica | 55–56 | gauge variation = 0 | claim 2 | yes |
| B2 | mathematica | 77–79 | EL residuals = 0 | claim 1 (setup) | yes |
| B3 | mathematica | 91 | `LinearSolve` Σ = closed form | claim 1 (independent route) | yes |
| B4 | mathematica | 95–96 | A/W residuals = 0 | claim 1 | yes |
| B5 | mathematica | 104–119 | z0/z2/z4 toy + Σ substitution = 0 | claim 3 | yes |
| B6 | mathematica | 146–147 | N(ω), N(0) closed forms = 0 | claim 4 | yes |
| B7 | mathematica | 136–140 | δD composed = 0 | claim 6 | yes |
| B8 | mathematica | 167–171 | Ŷ₂, Γ₅ = 0 | claim 5 | yes |
| B9 | mathematica | 191–192 | scalar terms = 0 | claim 7 | yes |

No tautological rows: each `expect_zero` compares an independently-derived quantity (an EL output, a series expansion, a `LinearSolve` result, a Hankel-function series) against a separately-written closed form, so any algebra slip would surface as a nonzero residual. The `z0/z2/z4` checks are doubly anchored: the toy-rational series (A5/B5) verifies the generic expansion identity, and the substituted-Σ checks (A6) verify that the physical numerator/denominator data plug in correctly — neither is guaranteed by construction.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl` (mtime 1780523951 = 2026-06-03 15:59:11)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt` (mtime 1779769478 = 2026-05-25 22:24:38)

**What's wrong:**
The Mathematica saved output predates the current `.wl` script by ~9 days. Git (`git show e2a4780`) shows the only post-output change to the `.wl` was a single banner-string fix at line 35: `banner["STAGE 004 — …"]` → `banner["STAGE 021 — …"]` (the "numbering reconciliation Phase 1 (deterministic): 273 doc-only stage-label fixes" commit, byte-identical except the number token). The captured output therefore still shows the stale banner `STAGE 004 — MAXWELL + MIXED-SECTOR REDUCTION` (output line 3) and the closing line `Stage 004 Mathematica audit passed.` (output line 82). No equation, value, variable, or assertion changed; the script body at line 195 still prints `"Stage 004 Mathematica audit passed."` (this `Print` literal was NOT part of the label fix and is itself a leftover stale label inside the live script, but it is a non-load-bearing print string, not a result). Every math residual in the captured output is `= 0` and matches what the current script produces.

**Why this matters:**
Purely informational. The math content of the saved output is correct and current; only the banner/closing label text is stale. Per the audit policy `stale_output` is non-blocking unless the output's *content* disagrees with the current script, which it does not here (only a cosmetic header string differs).

**Required change:**
Re-run the Mathematica script to refresh `mathematica/output/…stage021…_mathematica_audit.txt` so the captured banner reads `STAGE 021`. Optionally, the live script's leftover closing label at `…stage021…_mathematica_audit.wl:195` `Print["Stage 004 Mathematica audit passed."];` can be corrected to `Stage 021` for consistency with the line-35 banner — cosmetic only, no math impact.

**Verification:**
After a fresh run, output line 3 reads `STAGE 021 — MAXWELL + MIXED-SECTOR REDUCTION` and the closing line reads `Stage 021 Mathematica audit passed.` (if the line-195 print is also fixed); all `PASS:` lines and `= 0` residuals remain unchanged.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration:
- **Self-energy elimination:** SymPy writes the solution closed forms `A_sol = (gA·Wker + gW·R)/Δ`, `W_sol = …` by hand and checks residuals (py 145–149). Mathematica instead forms the kernel matrix `matEAW = {{aKer,-r},{-r,wKer}}` and calls `solAW = LinearSolve[matEAW, {gA, gW}]` (wl 86–87), then adds an extra cross-check `sigmaCons from LinearSolve matches closed form` (wl 90–91) that has no SymPy counterpart. Genuinely different route to the same Σ.
- **Hankel fingerprint:** SymPy builds the `l=2` spherical functions by hand — `j2a = (3/za³ − 1/za)sin(za) − 3cos(za)/za²`, `y2a = …`, `h2a = j2a + I·y2a` (py 273–275). Mathematica uses the built-in `h2a = SphericalHankelH1[2, za]` (wl 155). Independent constructions of the same special function.
- **Euler–Lagrange:** SymPy `sympy.calculus.euler.euler_equations` (py 123–125) vs Mathematica `VariationalMethods`' `EulerEquations` (wl 76). Different library implementations.
These are three independent derivation choices, so the engines corroborate rather than echo each other. No `mathematica_transliteration` finding.

## Engine cross-check

The two saved outputs agree on every closed form (modulo cosmetic sign-grouping):
- Σ_cons: sympy `(−2Rg_Ag_W + g_A²(−Ω_W²+ω²) + g_W²(−Ω_A²+ω²))/(R² − (Ω_A²−ω²)(Ω_W²−ω²))` (sympy out 28–33) ≡ Mathematica `(gW²(−oA²+ω²) + gA²(ω−oW)(ω+oW) − 2gAgWr)/((oA−ω)(oA+ω)(ω−oW)(ω+oW) + r²)` (wl out 25). Same rational function.
- z0: both `(Ω_A²g_W² + Ω_W²g_A² + 2Rg_Ag_W)/(Ω_A²Ω_W² − R²)` (sympy out 58–63; wl out 37).
- N(0): both `(Ω_A²g_W + Rg_A)²/(Ω_A²Ω_W² − R²)²` (sympy out 123–127; wl out 54).
- Ŷ₂: both `1 + ω²a²/(9c_s²) + 4ω⁴a⁴/(81c_s⁴) + i ω⁵a⁵/(27c_s⁵)` (sympy out 152–157; wl out 65).
- Γ₅: both `a⁵/(27c_s⁵)` (sympy out 160–165; wl out 66).
- N_scalar / scalar odd order: identical (sympy out 171–183; wl out 75–76).
All `expect_zero`/`expectZero` residuals are `= 0` / `PASS` in both transcripts. `engines_agree: true`.

## Verdict justification

`verdict: findings` solely because of the single low-severity, informational `stale_output` (F1) — the Mathematica output's banner text is one numbering-fix commit behind the script, with zero math change. I attacked the math on every front and it held: (a) hand-checked the gauge-invariance cancellations (`∂_t∂_w χ = ∂_w∂_t χ`) — exact; (b) re-derived Σ as the Schur complement `[g_A g_W]·M⁻¹·[g_A;g_W]` — matches the card eq:app-stage021-self-energy; (c) confirmed `N(ω)=(A g_W+R g_A)²/Δ²` is the first `∂/∂Π` of Σ_full at Π=0 (matches eq:app-stage021-transfer-factor, and Mathematica computes it that very way at wl 133); (d) checked the `l=2` series term-by-term against eq:app-stage021-outgoing-fingerprint including the `4/81` and `i/27` coefficients; (e) verified the scalar `gA=0, gW=ηω` substitution genuinely pushes the leading order from ω to ω² in N and to ω³ after `Π₀=iγ₁ω`. The assertions are non-tautological (closed forms written independently of their derivations), both engines derive the results by genuinely different routes, every saved residual is zero, and every paper deliverable — including the notes-only scalar-compatibility criterion — is exercised. I read the card, the notes, and the appendix row; the script's claims match the paper exactly. No `paper_misalignment`, no symbol-domain error, no missing branch.

## Value Reconciliation (pass-2 augmentation)

`reconciliation: complete; 11 values checked, 0 misaligned`

Every RESULT/deliverable value the scripts emit is reflected, with the correct value, in the card and/or notes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Σ_cons = (g_A² W + 2 g_A g_W R + g_W² A)/(A W − R²)` | py 136 / wl 84; out: sympy L28–33, wl L25 | tex L29–34 (eq:app-stage021-self-energy); md L116–117 | MATCH |
| `E_w = −∂_t A_w − ∂_w A_0` gauge-invariant; `C_a = ∂_a A_w − ∂_w A_a` gauge-invariant | py 74–75,88–89 / wl 47–48,55–56; out: sympy L9–10, wl L9–11 | tex L20–22 (eq:app-stage021-mixed-fields); md L70–73 | MATCH |
| `z0 = N0/D0` (= `(Ω_A²g_W²+Ω_W²g_A²+2Rg_Ag_W)/(Ω_A²Ω_W²−R²)`) | py 175 / wl 117; out: sympy L58–63, wl L37 | tex L38 (eq:app-stage021-z-coeffs, z₀); md L147 | MATCH |
| `z2 = (N0 S2 − G2 D0)/D0²` | py 176 / wl 118; out: sympy L64–69, wl L38 | tex L39 (z₂); md L149 | MATCH |
| `z4 = (N0(S2²−D0) − S2 G2 D0)/D0³` | py 177 / wl 119; out: sympy L70–87, wl L39 | tex L40 (z₄); md L151 | MATCH |
| `N(ω) = (A g_W + R g_A)²/(A W − R²)²` | py 222 / wl 146; out: sympy L102–119, wl L53 | tex L48–50 (eq:app-stage021-transfer-factor); md L178 | MATCH |
| `N(0) = (Ω_A²g_W + R g_A)²/(Ω_A²Ω_W² − R²)² ≥ 0` | py 228 / wl 147; out: sympy L123–127, wl L54 | tex L55–57 (eq:app-stage021-n0-positive); md L182, L312 | MATCH |
| `Ŷ₂ = 1 + a²ω²/(9c_s²) + 4a⁴ω⁴/(81c_s⁴) + i a⁵ω⁵/(27c_s⁵) + O(ω⁶)` | py 290 / wl 169; out: sympy L152–157, wl L65 | tex L64–69 (eq:app-stage021-outgoing-fingerprint); md L218–219 | MATCH |
| `Γ₅(port) = a⁵/(27 c_s⁵)` | py 302 / wl 171; out: sympy L160–165, wl L66 | tex L74 (a⁵/27c_s⁵ inside eq:app-stage021-wall-odd); md L223 | MATCH |
| `δD₂^odd = −i N₂(0) (a⁵/27c_s⁵) ω⁵ + O(ω⁷)` | py 241 / wl 138; out: sympy L129, wl L50 | tex L73–75 (eq:app-stage021-wall-odd); md L234 | MATCH |
| scalar: `N_scalar ≈ η²Ω_A⁴ω²/Δ₀²` and `δD₀^odd = iγ₁η²Ω_A⁴ω³/Δ₀²` (→ iω³) | py 335,344 / wl 191–192; out: sympy L171–183, wl L75–76 | md L262, L270–272 (notes §8); appendix row L64 "scalar-compatibility criterion" | MATCH |

INTERNAL scaffolding (no finding): EL-equation residuals (`Q/A/W equation`); `A/W exact solution residual`; `sigmaCons from LinearSolve matches closed form` (Mathematica-only consistency check); the generic toy-rational `z0/z2/z4 formula` checks; `Sigma_full(ω)` intermediate (pre-expansion); `Lambda2(k)` series (intermediate to Ŷ₂); `D_cons(ω)` (printed-only, derived directly from Σ_cons, no separate paper-side number); all pass/fail flags and `= 0` residuals.

## Self-test notes

I checked the derivative-independence trap (Mathematica's `D[sigmaFull, piOut] /. piOut->0` at wl 133 genuinely depends on `piOut`, so the derivative is not identically zero — it reproduces `(A g_W + R g_A)²/Δ²`); the parity/series trap (the `Series[…, {omega,0,5}]`/`{omega,0,2}` expansions land on the claimed ω-powers, and the scalar lane's substitution `gW=ηω` correctly raises the leading order so `N_scalar` starts at ω² and `Π₀·N` at ω³, not ω); and the trivial-case trap (gauge variations cancel by the symmetry of mixed partials, exactly zero). No assertion is tautological — each closed form is written independently of its derivation route. The only finding is the cosmetic banner staleness; it carries no math consequence and is informational.
