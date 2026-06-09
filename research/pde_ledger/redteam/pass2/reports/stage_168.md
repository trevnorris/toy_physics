---
unit_id: 168
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage168_off_bundle_slippage.md]
  paper_appendix: present
---

# Audit unit 168 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_168.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage168_off_bundle_slippage.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 67 table row; §"Off-bundle slippages and grouped signature" lines 363-384; overview line 265)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage168_off_bundle_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage168_off_bundle_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.txt`

## What the paper claims

Stage 168 (`off_bundle_slippage`) reduces the first true first-order off-bundle defect to three scalar slippage variables. `\stagefield{Output}`: *"Reduces the first true off-bundle defect to axial-length, mixed-speed, and mouth-traction slippages."* The notes make this concrete: starting from the Stage-164 normal coordinate `δ_⊥` (with coefficients `A_*, B_*, C_*` built from `g_*` and `s=√(1+r_*²)`) and the Stage-165 exact lower-branch transport laws, substituting the laws gives `δ_⊥=0` (bundle tangency, §1). Defining the three off-bundle slippages `ε_L, ε_v, ε_T` (boxed, §2) and re-substituting yields the boxed exact identity `δ_⊥ = -[g_* ε_T + (g_*+1/(2s)) ε_v + (2g_*+3/(4s)) ε_L] = -ε_⊥` (§3). The notes further carry: mouth-bias transport `δΠ` (§5), the four outlet defects `δC, δE2, δE4, Δ_Q` re-expressed via `ε_⊥` (§6), the even-preservation theorem `ε_⊥=0, δκ_W=0` (§7), and Family-1 numeric coefficients (§4: `r_*≈1.77799353547498`, `g_*≈0.758035078944663`, and the three boxed `δ_⊥` coefficients; §5: the three `δΠ` slippage coefficients). The appendix (lines 363-384) confirms the three slippages as the stage's deliverable. The card itself is terse (status card, "not a second independent proof") — the notes are authoritative on the equations.

## What the script claims to verify

The SymPy docstring lists four checks: (1) the Stage-164 normal coordinate in microscopic log drifts; (2) exact cancellation of the Stage-165 transport laws (`expect_zero` on `delta_perp` evaluated on the exact lower branch → 0, i.e. bundle tangency); (3) exact reduction of `delta_perp` to the three slippages (`expect_zero(delta_perp_slip + eps_perp)` → 0); (4) exact transport of the mouth-bias variation (`deltaPi` identity → 0) and the four outlet defects (`dC, dE2, dE4, DeltaQ`, each `expect_zero` against the `ε_⊥`-form → 0). It then prints the Family-1 numeric coefficients (epsT/epsv/epsL in `delta_perp` and in `deltaPi`). The Mathematica script asserts the identical set via `expectZero`/`pass`/`fail` and prints the identical numeric block. Both scripts cover exactly the §1-§6 algebra; neither asserts the §7 even-preservation collapse (that is a downstream determinant argument inheriting Stage 159, and the appendix flags §7 as a re-reading, not a new identity).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Stage-164 `δ_⊥` in log drifts (§1) | `delta_perp` def (py 44-52 / wl 39-47), printed | match |
| Bundle tangency `δ_⊥=0` on exact lower branch (§1) | `expect_zero` py 65-68 / `expectZero` wl 58-61 → 0 | match |
| Three slippages `ε_L,ε_v,ε_T` (boxed, §2) | substitution `subs_slip` py 71-75 / wl 65-69 | match |
| Exact reduction `δ_⊥=-ε_⊥` (boxed, §3) | `expect_zero(delta_perp_slip + eps_perp)` py 80 / wl 73 → 0 | match |
| Mouth-bias `δΠ` transport (boxed, §5) | `deltaPi` identity py 84-90 / wl 76-78 → 0 | match |
| Outlet `δC` (boxed, §6) | py 99 / wl 85 → 0 | match |
| Outlet `δE2` (boxed, §6) | py 100-103 / wl 86-89 → 0 | match |
| Outlet `δE4` (boxed, §6) | py 104-107 / wl 90-93 → 0 | match |
| Outlet `Δ_Q` (boxed, §6) | py 108-111 / wl 94-97 → 0 | match |
| Family-1 `δ_⊥` coeffs (§4) | printed py 122-124 / wl 108-110 | match |
| Family-1 `δΠ` coeffs (§5) | printed py 130-132 / wl 111-113 | match |
| Even-preservation `ε_⊥=0,δκ_W=0` (§7) | (none — downstream determinant re-reading; not a new identity) | n/a |

`paper_alignment: aligned`. Every boxed identity and every numeric deliverable in the notes is exercised. §7 is a re-reading of Stage 159's determinant argument under the new variables; the appendix and notes both present it as inheritance, not as a fresh identity the audit card must independently re-prove, so its absence is not a gap.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 65-68 | `expect_zero(delta_perp on lower branch)` | bundle tangency (§1) | yes |
| A2 | sympy | 80 | `expect_zero(delta_perp_slip + eps_perp)` | exact reduction `δ_⊥=-ε_⊥` (§3) | yes |
| A3 | sympy | 90 | `expect_zero(deltaPi - deltaPi_expected)` | mouth-bias transport (§5) | yes |
| A4 | sympy | 99 | `expect_zero(dC + 16σ ε_⊥/s)` | outlet `δC` (§6) | yes |
| A5 | sympy | 100-103 | `expect_zero(dE2 - …)` | outlet `δE2` (§6) | yes |
| A6 | sympy | 104-107 | `expect_zero(dE4 - …)` | outlet `δE4` (§6) | yes |
| A7 | sympy | 108-111 | `expect_zero(DeltaQ - …)` | outlet `Δ_Q` (§6) | yes |
| B1 | math | 58-61 | `expectZero(delta_perp on lower branch)` | bundle tangency (§1) | yes |
| B2 | math | 73 | `expectZero(delta_perp_slip + eps_perp)` | exact reduction (§3) | yes |
| B3 | math | 78 | `expectZero(deltaPi - deltaPiExpected)` | mouth-bias transport (§5) | yes |
| B4 | math | 85 | `expectZero(dC + 16σ ε_⊥/s)` | outlet `δC` (§6) | yes |
| B5 | math | 86-89 | `expectZero(dE2 - …)` | outlet `δE2` (§6) | yes |
| B6 | math | 90-93 | `expectZero(dE4 - …)` | outlet `δE4` (§6) | yes |
| B7 | math | 94-97 | `expectZero(deltaQ - …)` | outlet `Δ_Q` (§6) | yes |

No tautological rows. Each `expect_zero` subtracts an **independently transcribed target** (`eps_perp`, the four outlet `ε_⊥`-forms, the `δΠ_expected` form) from a quantity assembled by **substitution** into `delta_perp`/`delta_perp_slip`. The targets are typed separately from the construction (e.g. `eps_perp = g*epsT + (g+B)*epsv + C*epsL` at py 77 is NOT used to build `delta_perp_slip`; if any boxed RHS coefficient in the notes were wrong, A2-A7 would produce a nonzero residual). A1 is the genuine cancellation of seven log-channel terms against the two transport laws — a real algebraic identity, not `X-X`. These are non-tautological and substantive.

## Findings

### F1 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl:35-121`

**What's wrong:**
The `.wl` is a verbatim, line-for-line port of the `.py`; it adds no independent derivation route. Every step corresponds 1:1, in the same order, with the same variable choreography and the same target literals. The second engine echoes the first engine's algebra rather than deriving the result independently. Three corresponding sections:

1. `delta_perp` construction — same seven-term sum, same coefficient names, same order:
   - py 40-52: `A = g + Rational(1,4)/s` … `delta_perp = A*Xi + 3*A*csw + B*cs - g*Tm - (g+B)*vw0 - 2*A*a - C*LW`
   - wl 35-47: `aCoeff = g + 1/(4*s)` … `deltaPerp = Expand[aCoeff*xi + 3*aCoeff*csw + bCoeff*cs - g*tM - (g+bCoeff)*vw0 - 2*aCoeff*a - cCoeff*lW]`

2. `eps_perp` target and reduction check — the SAME target literal is transcribed into both, then both subtract it from the same substituted expression:
   - py 77,80: `eps_perp = g*epsT + (g+B)*epsv + C*epsL` ; `expect_zero("delta_perp + eps_perp", delta_perp_slip + eps_perp)`
   - wl 71,73: `epsPerp = g*epsT + (g+bCoeff)*epsv + cCoeff*epsL` ; `expectZero["delta_perp + eps_perp", deltaPerpSlip + epsPerp]`

3. Outlet block — identical four checks, identical order, identical right-hand-side `ε_⊥`-forms:
   - py 94-111: `dC = 16*sigma*delta_perp_slip/s` … `DeltaQ - sigma*(-16*eps_perp/s - 27*dgammaW)/(3*(1-sigma))`
   - wl 81-97: `dC = 16*sigma*deltaPerpSlip/s` … `deltaQ - sigma*(-16*epsPerp/s - 27*dGammaW)/(3*(1-sigma))`

Even the trailing carry-forward `Print` strings are byte-for-byte identical (py 134-139 ↔ wl 116-121). This is the same situation flagged at Stage 165 (F1). The achievable independence is genuinely limited — the stage is pure algebraic substitution-and-cancellation of imported carry-forward formulas, with no `Solve`/`Series`/matrix route to take — so this is a low-severity policy-fidelity flag, not a correctness flag; the engines already agree exactly (residuals all 0).

**Why this matters:**
The second-engine policy requires both engines to reach the result independently so that a transcription slip in one is caught by the other. A verbatim port means a wrong coefficient typed into the shared target literal (e.g. `eps_perp`) would be copied identically into both engines and pass in both. The cross-check value is reduced to "two transcriptions of the same algebra agree," not "two independent derivations agree."

**Required change:**
Give the Mathematica script at least one genuinely independent route to the load-bearing reduction result, rather than transcribing the same `epsPerp` target. The achievable independence here is to DERIVE `epsPerp` instead of writing it as a literal: collect/solve `deltaPerpSlip` for the negated slippage combination and confirm the solved combination equals the boxed `ε_⊥`. Concretely (see directive): replace the literal `epsPerp` at wl:71 with a derived quantity obtained by `Coefficient`-extracting the `{epsT, epsv, epsL}` weights out of `-deltaPerpSlip`, then assert the derived weights match the boxed `(g_*, g_*+1/(2s), 2g_*+3/(4s))`. This makes the `.wl`'s path to `ε_⊥` distinct from the `.py`'s hand-typed target. If reconstruction is judged too underspecified to apply mechanically, append `## Blocked: F1`.

**Verification:**
After Codex applies, the verifier runs `redteam exec-mathematica 168` and confirms (a) wl:71 no longer holds the literal `epsPerp = g*epsT + (g+bCoeff)*epsv + cCoeff*epsL` verbatim but a derived form, (b) a new "epsPerp weights match boxed form" `expectZero` PASSES, (c) all seven original `expectZero` checks still PASS, and (d) the script exits 0.

## Independent-derivation check (Mathematica)

Not independent — verbatim transliteration (high confidence). See F1 for the three quoted corresponding sections. The `.wl` introduces no derivation route absent from the `.py` (no `Solve`, no `Series`, no matrix inverse, no implicit-function step — contrast Stages 163/164/166, which each added such a route and were cleared). The shared scaffolding (banner/expect helpers, same physical premises, same final target) is expected; the defect is that the derivation PATH is identical, including the same hand-typed `epsPerp` target.

## Engine cross-check

Both engines agree exactly. Side-by-side final outputs:

- `delta_perp` (printed): same expression up to ordering — py output line 5 (`-2*LW*g - 3*LW/(4*sqrt(r**2+1)) - Tm*g + Xi*g + …`) ↔ wl output line 5 (`-2*a*g + 3*csw*g - 2*g*lW - a/(2*s) + … + g*xi + xi/(4*s)`). Identical term set.
- `delta_perp with slippages`: py `(-3*epsL/4 - epsv/2 + g*sqrt(r**2+1)*(-2*epsL - epsT - epsv))/sqrt(r**2+1)` ↔ wl `-2*epsL*g - epsT*g - epsv*g - (3*epsL)/(4*s) - epsv/(2*s)` — same after clearing `s`.
- All seven `expectZero`/`expect_zero` residuals = 0 in both (py output lines 13,19,24,29-32; wl output lines 13-14,19-21,26-27,32-39 with PASS lines).
- Numeric coefficients agree to 20 digits: e.g. epsT in delta_perp `-0.75803507894466300000` (py) ↔ `-0.75803507894466304328…` (wl); epsL in deltaPi `-2.8791587799041576371` (py) ↔ `-2.87915877990415757043…` (wl). No `engine_disagreement`.

`stale_output`: not applicable. Both outputs are newer than their scripts (py: script mtime 1780005291 < output 1780006214; wl: script 1780005293 < output 1780006293). Outputs fresh.

## Verdict justification

`verdict: findings`, one finding (F1, `mathematica_transliteration`, low severity). The math is correct and the paper alignment is exact: every boxed identity (§1 tangency, §3 `δ_⊥=-ε_⊥`, §5 `δΠ`, §6 four outlet defects) is asserted with non-tautological `expect_zero`/`expectZero` checks that subtract an independently typed target, and both engines pass with zero residuals. Numeric Family-1 coefficients (§4, §5) reconcile to ≥12 digits against the notes. Attacks tried and failed: (a) checked A2-A7 for `X-X` tautology — the `eps_perp`/outlet targets are typed separately from the substituted construction, so a wrong boxed coefficient would fail; (b) checked symbol domains — `g>0, r>0` plus reals are the correct physical domain and do not over-constrain any `simplify`; (c) searched for the recurring stale `168π²/100π²` constant and the Family-1 radius `√(4107-100π²)/(10π)` — neither appears in this stage (the only `168` is the stage number), so no stale-constant contamination here; (d) confirmed `Sigma0=4.651033550168876` and `S_can=0.6703621156734617` trace to upstream Stages 156/158/163 and to this stage's §5, so they are sourced, not orphaned hardcodes. The single defect is a policy-fidelity one: the `.wl` is a verbatim port of the `.py` with no independent route, identical to the calibrated Stage-165 transliteration finding.

## Value Reconciliation (pass-2 augmentation)

Every RESULT/deliverable value the scripts emit, located in the docs:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `r_* = 1.77799353547498` | py:114 / wl:100 (input, echoed via coeffs) | notes:192 `\mathfrak r_*=\mathfrak r_{F1}\approx 1.77799353547498` | MATCH |
| `g_* = 0.758035078944663` | py:115 / wl:101 (input, echoed) | notes:194 `\mathfrak g_*\approx 0.758035078944663` | MATCH |
| coeff `ε_T` in `δ_⊥` = `-0.758035078944663…` | py out:37 / wl out:44 | notes:202 `-0.758035078944663\,\varepsilon_T` | MATCH |
| coeff `ε_v` in `δ_⊥` = `-1.00314310113848…` | py out:38 / wl out:45 | notes:203 `-1.00314310113848\,\varepsilon_v` | MATCH |
| coeff `ε_L` in `δ_⊥` = `-1.88373219118005…` | py out:39 / wl out:46 | notes:204 `-1.88373219118005\,\varepsilon_L` | MATCH |
| coeff `ε_T` in `δΠ` = `-1.15860596492310…` | py out:40 / wl out:47 | notes:273 `-1.15860596492310\,\varepsilon_T` | MATCH |
| coeff `ε_v` in `δΠ` = `-1.53323719829507…` | py out:41 / wl out:48 | notes:275 `-1.53323719829507\,\varepsilon_v` | MATCH |
| coeff `ε_L` in `δΠ` = `-2.87915877990416…` | py out:42 / wl out:49 | notes:277 `-2.87915877990416\,\varepsilon_L` | MATCH |
| symbolic `δ_⊥ = -ε_⊥` (`delta_perp_slip + eps_perp = 0`) | py:80 / wl:73 → 0 | notes:176-178 boxed `\delta_\perp=-\varepsilon_\perp` | MATCH |
| symbolic outlet `δC, δE2, δE4, Δ_Q` `ε_⊥`-forms | py:99-111 / wl:85-97 → 0 | notes:325-366 boxed §6 | MATCH |
| symbolic `δΠ` transport form | py:84-90 / wl:76-78 → 0 | notes:234-242 boxed §5 | MATCH |

INTERNAL (scaffolding/inputs not expected as new prose deliverables, no finding): `Sigma0_num = 4.651033550168876` and `S_num = 0.6703621156734617` (carry-forward inputs sourced from Stages 156/158/163 and used in this stage's §5 to produce the `δΠ` coefficients — they ARE present in §5's `δΠ_tan` line implicitly, e.g. `0.832409471081635 = 1 - S_can/4` and `1.16275838754222 = Sigma0/4` at notes:248-250, so even their derived combinations reconcile); `delta_perp` printed expression; `delta_perp_slip` printed expression; `Tm_br`/`vw0_br`/`LW_br` printed transport laws (notes §1, MATCH but intermediate); bundle-tangency residual (0); the carry-forward `Print` text block (definitional, mirrors notes §2-§3).

reconciliation: complete; 11 deliverable values checked, 0 misaligned.

## Self-test notes

Checked: (1) variable independence — no `sp.diff`/`D[]` in this stage, so the zero-derivative trap does not apply; all checks are substitution-and-simplify. (2) Symmetry/parity — no integrals; N/A. (3) Trivial-case pre-check — for the F1 directive, the proposed derived `epsPerp` via `Coefficient[-deltaPerpSlip, {epsT,epsv,epsL}]` reproduces exactly `{g, g+1/(2s), 2g+3/(4s)}` (read off the printed `delta_perp_slip = -2*epsL*g - epsT*g - epsv*g - 3*epsL/(4s) - epsv/(2s)`: `epsT`→`g`, `epsv`→`g+1/(2s)`, `epsL`→`2g+3/(4s)`), so the new `expectZero` reduces to 0 — verified mentally. (4) Path spec — F1 targets the existing `.wl` in `mathematica/`, no new file. (5) Paper round-trip — the F1 fix only re-derives the same `ε_⊥` weights the notes box, introducing no new constant, so no new `paper_misalignment`.
