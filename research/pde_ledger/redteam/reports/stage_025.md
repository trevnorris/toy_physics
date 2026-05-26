---
unit_id: 025
batch: II.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
  paper_appendix: present
---

# Audit unit 025 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_025.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (relevant rows around line 40 and the Part II target definition around lines 72-81)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.txt`

## What the paper claims

Per `\stagefield{Output}` (line 88): "Stage 025 outputs the explicit one-mode target [eq:app-stage025-normalization-target], the stability gate [eq:app-stage025-stability], and the first monotonicity signs [eq:app-stage025-monotonicity]." The body of the card lays out: definitions C, G_U, G_W, R (eq. couplings); the invariants Delta = OmegaU^2 OmegaW^2 - R^2, Q = G_U^2 OmegaW^2 + 2 G_U G_W R + G_W^2 OmegaU^2, P = OmegaU^2 G_W + R G_U; the BdG softening X = C^2/varpi^2; D_0 = K - X - Q/Delta and N_0 = P^2/Delta^2; the compact form P_0 = N_0/D_0 = P^2/(Delta (K Delta - Delta C^2/varpi^2 - Q)); the target equation mhat_rad^2 P_0 = N_Q^target; the stability gate Delta>0 and D_0>0; and the monotonicity signs partial_K P_0 = -N_0/D_0^2 < 0 and partial_X P_0 = +N_0/D_0^2 > 0. The Part II target N_Q^target := 54 G c_s^5 / (5 a^5 c^5) is defined in the appendix (lines 72-81). The card's `\stagefield{Checks}` enumerates four corollary checks (substitution to compact form; angular source map absence; sign of derivatives from positivity; P = 0 ⇒ N_0 = 0 ⇒ no positive target). Notes agree and add the operational reading.

## What the script claims to verify

Both scripts define the same invariants Delta, Q, P, B_0 = C^2/varpi^2, Z_0 = Q/Delta, N_0 = P^2/Delta^2, D_0 = K - B_0 - Z_0, and then assert: (i) raw P_0 = N_0/D_0 equals the compact form P^2/(Delta(K Delta - Delta C^2/varpi^2 - Q)); (ii) numerical sample point gives P_0 = 1/3 (Delta = 15, D_0 = 1/3); (iii) the target equation mhat^2 P_0 - 54 G c_s^5 / (5 a^5 c^5) is solvable with mhat^2 > 0 at the sample; (iv) algebraic identity Delta*D_0 = K Delta - Delta C^2/varpi^2 - Q; (v) algebraic identity N_0 = P^2/Delta^2 (this is tautological — see F1); (vi) Delta, D_0, P_0 all positive at the sample; (vii) closed-form derivatives dP_0/dK = -N_0/D_0^2 and dP_0/dX = +N_0/D_0^2 verified against sp.diff / D / Limit, with sample-point sign checks dP_0/dK < 0 and dP_0/dX > 0; (viii) the conservation identity dP_0/dK + dP_0/dX = 0. The Mathematica script additionally cross-checks the symbolic derivatives via Limit-versus-D and runs Reduce on the target equation for the mhat > 0 branch.

## Paper ↔ script cross-check

| Paper deliverable | Script side | Status |
|---|---|---|
| Definitions Delta, Q, P (eq:invariants) | sympy lines 74-76, math lines 41-43 | match |
| X = C^2/varpi^2 (eq:X) | sympy 78 (`B0 = C^2/varpi^2`), 60 (`X` as Section IV symbol with sample X=1), math 44, 35, 140 | match |
| D_0 = K - X - Q/Delta, N_0 = P^2/Delta^2 (eq:D0-N0) | sympy 80-81, math 46-47 | match (N_0 identity check is tautological — F1) |
| Compact P_0 form (eq:P0) | sympy 103, 108 (`P0 - P0_compact == 0`); math 67-68 | match |
| Target equation mhat_rad^2 P_0 = N_Q^target (eq:normalization-target) | sympy 123-138 (residual + sample mhat^2 > 0); math 84-100 (residual + sample + Reduce-based mhat > 0 solvability) | match |
| Stability gate Delta > 0, D_0 > 0 (eq:stability) | sympy 159-164, math 109-118 | partial (only sample-point spot check; algebraic positivity of P^2/Delta^2 is trivially true; no general positivity assertion) |
| Monotonicity signs partial_K P_0 = -N_0/D_0^2 < 0, partial_X P_0 = +N_0/D_0^2 > 0 (eq:monotonicity) | sympy 186-201, math 135-148 | match (closed-form identities + sample sign + sum-to-zero) |
| Substitution N_0,D_0 → P_0 compact form (Checks #1) | sympy 108, 150; math 68, 107 | match |
| Angular source map absent (Checks #2) | descriptive consistency — no algebraic check needed | n/a |
| Sign of partial P_0 follows from N_0>0, D_0>0 (Checks #3) | sympy 196-199, math 145-146 | partial (sign verified at sample; algebraic D_0^2 > 0 is automatic) |
| P = 0 ⇒ N_0 = 0 ⇒ no positive quadrupole target (Checks #4) | not exercised | missing (F2) |

`paper_alignment: aligned` — the script's bottom-line identities and the target constant 54/5 match the paper card and the Part II target definition exactly. One trivially-true assertion (F1) and one missing corollary (F2) are the only gaps.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 108 | `expect_zero("P0 - P0_compact", P0 - P0_compact)` | compact P_0 form (eq:P0) and Check #1 | yes |
| A2 | sympy | 114-115 | `P0_raw at sample == 1/3` | sanity on (i) | yes (numeric) |
| A3 | sympy | 116-117 | `P0_compact at sample == 1/3` | sanity on (i) | yes (numeric) |
| A4 | sympy | 137-138 | `mhat^2 at sample > 0` | target equation solvability | partial (sample only; no symbolic proof on the stability branch — Mathematica's Reduce supplies this) |
| A5 | sympy | 150 | `expect_zero("Delta*D0 - (K*Delta - Delta*C^2/varpi^2 - Q)", ...)` | distributive expansion of denominator | partial (mechanical algebra; already implied by A1) |
| A6 | sympy | 151 | `expect_zero("N0 - P^2/Delta^2", N0 - P^2/Delta^2)` | N_0 = P^2/Delta^2 (eq:D0-N0) | no — tautological (F1) |
| A7 | sympy | 159-160 | `Delta on sample > 0` | stability gate Delta > 0 | partial (sample only) |
| A8 | sympy | 161-162 | `D0 on sample > 0` | stability gate D_0 > 0 | partial (sample only) |
| A9 | sympy | 163-164 | `P0 on sample > 0` | stability ⇒ P_0 > 0 | partial (sample only) |
| A10 | sympy | 186 | `dP0/dK + N0/(K-X-Q/Delta)^2 == 0` | dP_0/dK = -N_0/D_0^2 | yes (algebraic; sp.diff vs hand form) |
| A11 | sympy | 187 | `dP0/dX - N0/(K-X-Q/Delta)^2 == 0` | dP_0/dX = +N_0/D_0^2 | yes |
| A12 | sympy | 188 | `dP0/dX + dP0/dK == 0` | corollary of A10+A11 | yes (genuine sum identity) |
| A13 | sympy | 196-197 | `dP0/dK at sample < 0` | sign of dP_0/dK | partial (sample only) |
| A14 | sympy | 198-199 | `dP0/dX at sample > 0` | sign of dP_0/dX | partial (sample only) |
| A15 | sympy | 200-201 | `dP0/dK + dP0/dX == 0 at sample` | corollary | yes |
| B1 | mathematica | 68 | `P0 - P0_compact` | compact P_0 form | yes |
| B2 | mathematica | 74-76 | `P0 at sample == 1/3` (raw and compact) | sanity on (i) | yes |
| B3 | mathematica | 92 | `mhat^2 at sample > 0` | solvability at sample | partial |
| B4 | mathematica | 94-99 | `Reduce[mhat^2 == target/p0Compact && mhat > 0 && delta > 0 && d0 > 0, mhat, Reals]` | target solvability on the stability branch | yes (genuinely independent of sympy) |
| B5 | mathematica | 107 | `Delta*D0 - (K*Delta - Delta*C^2/varpi^2 - Q)` | distributive expansion | partial (subset of B1) |
| B6 | mathematica | 108 | `N0 - P^2/Delta^2` | N_0 = P^2/Delta^2 | no — tautological (F1) |
| B7 | mathematica | 115-117 | `Delta`, `D0`, `P0` at sample > 0 | stability gate | partial (sample only) |
| B8 | mathematica | 135-136 | `Limit dP0/dK - D[p0,k] == 0`, same for X | sp.diff/D cross-check | yes (Limit-quotient vs D operator) |
| B9 | mathematica | 137-139 | `dP0/dK + N0/(...)^2 == 0`, `dP0/dX - N0/(...)^2 == 0`, `dP0/dK + dP0/dX == 0` | dP_0/dK and dP_0/dX closed forms | yes |
| B10 | mathematica | 145-147 | sample-point signs and sum | sign of derivatives | partial |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:151`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:108`

**What's wrong:**
In both engines `N_0` is *constructed* literally as `P**2 / Delta**2` (sympy line 80: `N0 = sp.simplify(P**2 / Delta**2)`; math line 46: `n0 = FullSimplify[p^2/delta^2, ...]`). The assertion `expect_zero("N0 - P^2/Delta^2", N0 - P**2/Delta**2)` therefore reduces to `(P^2/Delta^2) - P^2/Delta^2 == 0`, which holds by construction. The check cannot fail regardless of whether the physics is right; it adds zero adversarial signal. The paper's eq:app-stage025-D0-N0 *defines* N_0 = P^2/Delta^2, so what would be substantively verifiable here is the upstream Stage 023/024 derivation that gives this form — but stage 023 already owns that.

**Why this matters:**
A passing tautological check inflates the apparent verification surface. If a future edit ever decouples the construction of N_0 from this identity (e.g., introducing a corrected expression) the check would still silently pass.

**Required change:**
Replace the redundant identity with a substantive structural check that re-expands N_0 from primitive symbols and confirms it matches the constructed N_0. Concretely: assert `(P_raw)**2 / Delta**2 - N_0 == 0`, where `P_raw = OmegaU**2 * GW + R * GU` is built locally from the primitive symbols rather than referencing the cached `P` variable. The check then verifies that the cached `P` and `Delta` build N_0 consistently with their primitive definitions, instead of just comparing N_0 to its own definition. (Simplest alternative: drop the line outright.) See directive F1.

**Verification:**
After Codex applies, sympy line 151 and math line 108 either disappear or are replaced with the substantive check noted above. The next `redteam exec-sympy 025` and `redteam exec-mathematica 025` runs still exit 0; if replaced, the saved transcript contains a non-trivially-zero residual (still zero, but with primitive-symbol structure visible in the expression printed before the simplification step).

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py` (no test for P = 0 case)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl` (no test for P = 0 case)

**What's wrong:**
The paper card's `\stagefield{Checks}` (stage_025.tex line 94) enumerates four corollary checks. The fourth: "If P = 0, then N_0 = 0 and a positive quadrupole target cannot be reached." Neither engine substitutes P = 0 (equivalently Omega_U^2 G_W + R G_U = 0) and verifies N_0 = 0 nor that the target equation mhat^2 P_0 = N_Q^target has no positive-mhat solution. The corollary is a trivial substitution (N_0 = P^2/Delta^2 = 0 when P = 0, and N_Q^target > 0 by definition), but the paper explicitly enumerates it as a checked deliverable.

**Why this matters:**
The corollary anchors the physical reading "port-transfer combination P must not vanish" (notes section 3 line 137). A future change to the definition of P or to the target sign would silently keep the (currently absent) check passing. Treating the paper's enumerated checks as the audit contract requires the script to honor item #4.

**Required change:**
Add to Section II a short block that substitutes `P = 0` (e.g., set `GW_sub = -R*GU/OmegaU**2` so that `P = OmegaU**2*GW_sub + R*GU = 0` symbolically), confirms `N0.subs({GW: GW_sub}) == 0` symbolically (via expect_zero), and confirms `mhat**2 P0_compact.subs({GW: GW_sub}) - target` reduces to `-target` (negative, so no positive-mhat solution). Mirror in the Mathematica script. See directive F2.

**Verification:**
After Codex applies, the new block in both engines confirms `N0_at_Pzero == 0` and `(mhat^2 P0_compact - target) at P=0 simplifies to -target` (strictly negative when gConst, cs, cSpeed, a all positive). The `redteam exec-sympy 025` / `exec-mathematica 025` runs still exit 0, and the saved transcripts contain explicit `PASS: N0 vanishes when P=0` and `PASS: target unreachable when P=0` lines.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally parallel to the SymPy script: same section ordering (I → II.1 → II.2 → III → IV.1), same module/function decomposition (`zeroFrequencyCoefficients`, `normalizationFormula`, `stabilityAndPositivity`, `monotonicDerivatives`), same intermediate-variable names (`delta`, `q`, `p`, `b0`, `z0`, `n0`, `d0`, `p0`, `p0Compact`), and the same SAMPLE_POINT.

However, two assertions are genuinely independent of any algebraic choreography copied from SymPy:

1. Lines 94-99: `Reduce[mhat^2 == target/p0Compact && mhat > 0 && delta > 0 && d0 > 0, mhat, Reals]` — this is a global symbolic solvability statement on the stability branch, not present anywhere in the SymPy script (which only checks one sample point).
2. Lines 126-127, 135-136: derivatives via `Limit[((p0 /. k -> k+h) - p0)/h, h -> 0]` cross-checked against `D[p0, k]` — the SymPy script uses only `sp.diff`. The Limit-quotient route is an independent derivative computation.

Sympy: `Delta = OmegaU**2*OmegaW**2 - R**2`. Mathematica: `Delta = (omegaU*omegaW - r)*(omegaU*omegaW + r)` (factored). Same content, different canonical form chosen by each engine — this is mild evidence of engine independence.

Verdict on transliteration: borderline. The high-level choreography mirrors, but the two genuinely independent checks (Reduce and Limit-vs-D) keep it above the transliteration bar. Not flagged as a finding.

## Engine cross-check

The two engines compare identically on every shared computation:

| Quantity | SymPy | Mathematica | Match |
|---|---|---|---|
| `Delta` at sample | 15 | 15 | yes |
| `D_0` at sample | 1/3 | 1/3 | yes |
| `P_0` at sample (raw and compact) | 1/3 | 1/3 | yes |
| `mhat^2` at sample | 162/5 | 162/5 | yes |
| `dP_0/dK` at sample | -1 | -1 | yes |
| `dP_0/dX` at sample | +1 | +1 | yes |
| `P_0 - P0_compact` | 0 | 0 | yes |
| `Delta*D_0 - (K*Delta - Delta*C^2/varpi^2 - Q)` | 0 | 0 | yes |
| `dP_0/dK + N_0/D_0^2` | 0 | 0 | yes |
| `dP_0/dX - N_0/D_0^2` | 0 | 0 | yes |
| `dP_0/dX + dP_0/dK` | 0 | 0 | yes |

Engines agree at the level claimed.

## Verdict justification

The script's bottom-line identities — compact P_0, target equation, stability sample, closed-form derivatives, and conservation `partial_K + partial_X = 0` — match the paper card and Part II target definition exactly. Numerical sample agreement on Delta, D_0, P_0, mhat^2, and the two slopes is consistent across both engines. The target constant `54 G c_s^5 / (5 a^5 c^5)` matches `N_Q^target` in the appendix. Attacks attempted: (a) chase a sign-flip in `dP_0/dX` (would need d/dX of `K-X-Q/Delta` to be +1 instead of -1; checked: X is an independent sympy symbol, derivative is correct); (b) chase a tautology in N_0 vs P^2/Delta^2 (found one — F1); (c) check whether the `54/5` constant is paper-anchored (it is, in appendix line 78 and notes line 26); (d) check whether the Reduce-based mhat > 0 branch in Mathematica properly conditions on the stability gate (it does — line 95 includes `delta > 0 && d0 > 0`); (e) cross-check that the script's Section IV uses the same N_0 as Section I (it does — `N_0 = P^2/Delta^2`, independent of K and X); (f) look for a missing paper-card check (found one — F2, the P=0 corollary). Verdict is `findings` because of the two low-severity items. `paper_alignment: aligned` is preserved because both findings are script-internal (tautology, omitted corollary), not a paper-vs-script disagreement.

## Self-test notes

Checked: (1) variable independence — Section IV's `X` is an independent sympy symbol; `Q` and `Delta` do not depend on K or X, so `d(K-X-Q/Delta)/dK = 1` and `d/dX = -1`, giving the closed-form derivatives. (2) symmetry/parity — n/a, no integrals. (3) trivial-case pre-check for F2 — at P=0 (e.g., GW = -R*GU/OmegaU^2), N_0 = 0, so `mhat^2 P_0 = 0` and `mhat^2 P_0 - target = -target = -54 G c_s^5 / (5 a^5 c^5) < 0`, confirming the directive's proposed assertions are non-trivially true. (4) path specifications — F1 is an edit to existing lines, no new files; F2 adds checks to existing files at correct paths. (5) paper round-trip — neither fix introduces a new paper_misalignment: F1 removes/replaces a redundant identity (no constants involved); F2 adds the paper-card-enumerated P=0 corollary using the existing definitions of P, N_0, and N_Q^target as stated in the paper.
