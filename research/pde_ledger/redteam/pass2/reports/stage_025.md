---
unit_id: 025
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md]
  paper_appendix: present
---

# Audit unit 025 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_025.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 40; `\input` line 88)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.txt`

## What the paper claims

Stage 025 collapses the isotropic branch to the smallest radial/axial module (one BdG support mode `varpi`, one internal gauge coordinate `U`, one mixed coordinate `W`) and turns the formal ratio `N_0/D_0` into an explicit scalar. It defines the invariants `Δ = Ω_U²Ω_W² − R²`, `Q = G_U²Ω_W² + 2G_UG_WR + G_W²Ω_U²`, `P = Ω_U²G_W + RG_U` (eq:app-stage025-invariants), the BdG softening `X = C²/varpi²` (eq:app-stage025-X), and `D_0 = K − X − Q/Δ`, `N_0 = P²/Δ²` (eq:app-stage025-D0-N0). The bottom-line `\stagefield{Output}` is: "Stage 025 outputs the explicit one-mode target \eqref{eq:app-stage025-normalization-target}, the stability gate \eqref{eq:app-stage025-stability}, and the first monotonicity signs \eqref{eq:app-stage025-monotonicity}." The target equation (eq:app-stage025-normalization-target) is `mhat_rad² · P²/[Δ(KΔ − ΔC²/varpi² − Q)] = N_Q^target`; the notes pin `N_Q^target = 54 G c_s^5/(5 a^5 c^5)` (notes lines 101, 209, 236). The stability gate (eq:app-stage025-stability) is `Δ>0` and `D_0 = K − C²/varpi² − Q/Δ > 0`. The monotonicity (eq:app-stage025-monotonicity) is `∂P_0/∂K = −N_0/D_0² < 0` and `∂P_0/∂X = +N_0/D_0² > 0`. The card's `\stagefield{Checks}` adds: substitution of D0,N0 into P0=N0/D0 reproduces eq:app-stage025-P0; signs follow from N0>0,D0>0; and if `P=0` then `N0=0` and a positive quadrupole target is unreachable (Checks item 4).

## What the script claims to verify

Both engines build `Δ, Q, P, B0=C²/varpi², Z0=Q/Δ, N0=P²/Δ², D0=K−B0−Z0` from the symbolic definitions (sympy 74-81; wl 41-47). Section II asserts the raw `P0=N0/D0` equals the compact closed form `P²/[Δ(KΔ−ΔC²/varpi²−Q)]` symbolically (sympy 108; wl 68), pins both to `1/3` at the sample point (sympy 114-117; wl 74-75), forms the target `54 G c_s^5/(5 a^5 c^5)` carried forward from stage023, and checks the target equation `mhat²·P0_compact = target` is solvable for `mhat>0` (sympy 132-138 sample-point positivity; wl 89-99 sample positivity + a symbolic `Reduce` over `mhat>0 ∧ Δ>0 ∧ D0>0`). Section II.3 verifies the P=0 corollary: forcing `GW=−R·GU/ΩU²` zeroes `N0` and reduces the residual to `−target` (sympy 143-148; wl 102-111). Section III re-derives `Δ·D0` equals the compact denominator, reconstructs `N0` from raw symbols, and sample-checks `Δ>0, D0>0, P0>0` (sympy 160-176; wl 117-134). Section IV computes `∂P0/∂K` and `∂P0/∂X`, asserts each equals `∓N0/(K−X−Q/Δ)²`, that they sum to zero, and sign-checks them on the sample (sympy 191-213; wl 141-162).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Invariants Δ, Q, P (eq:app-stage025-invariants) | sympy 74-76 / wl 41-43, printed output lines 9-11 | match |
| X = C²/varpi², D0 = K−X−Q/Δ, N0 = P²/Δ² (eq:...-X, -D0-N0) | sympy 78-81 / wl 44-47, output 12-15 | match |
| P0 = N0/D0 = compact closed form (eq:app-stage025-P0) | sympy 108 expect_zero `P0 − P0_compact` / wl 68; out line 22 = 0 | match |
| One-mode target = N_Q^target = 54 G c_s^5/(5 a^5 c^5) (eq:...-normalization-target) | sympy 123-138 / wl 84-99; target solvable for mhat>0 (out: mhat²=162/5 at sample, Reduce ≠ False) | match |
| Stability gate Δ>0, D0>0 (eq:app-stage025-stability) | sympy 171-176 / wl 130-132 (sample) + sympy 160 / wl 118 (Δ·D0 ≡ compact denom) | match |
| Monotonicity ∂P0/∂K<0, ∂P0/∂X>0 (eq:...-monotonicity) | sympy 198-213 / wl 150-162 | match |
| Checks item 4: P=0 ⇒ N0=0, target unreachable | sympy 143-148 / wl 102-111 (II.3) | match |

`paper_alignment: aligned`. Every paper-side deliverable has a faithful, non-tautological script-side counterpart in both engines. No `extra` rows: every assertion traces to a paper claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 108 | `expect_zero(P0 − P0_compact)` | P0 closed form (eq:P0) | yes |
| A2 | sympy | 114-117 | `P0 raw/compact at sample == 1/3` | P0 numeric anchor | yes |
| A3 | sympy | 137-138 | `mhat_sq_at_sample > 0` | target solvability | yes |
| A4 | sympy | 145 | `expect_zero(N0 at P=0)` | Checks item 4 | yes |
| A5 | sympy | 148 | `expect_zero(residual + target)` at P=0 | Checks item 4 (target unreachable) | yes |
| A6 | sympy | 160 | `expect_zero(Δ·D0 − compact denom)` | stability denom identity | yes |
| A7 | sympy | 163 | `expect_zero(N0 − P_raw²/Δ_raw²)` | N0 definition | yes |
| A8 | sympy | 171-176 | sample `Δ,D0,P0 > 0` | stability gate | yes |
| A9 | sympy | 198-199 | `expect_zero(dP0/dK + N0/den²)`, `dP0/dX − N0/den²` | monotonicity forms | yes |
| A10 | sympy | 200 | `expect_zero(dP0/dX + dP0/dK)` | monotonicity antisymmetry | yes |
| A11 | sympy | 208-213 | sample `dP0/dK<0`, `dP0/dX>0`, sum=0 | monotonicity signs | yes |
| B1 | mathematica | 68 | `expectZero(p0 − p0Compact)` | P0 closed form | yes |
| B2 | mathematica | 74-75 | `p0Raw/p0Compact === 1/3` | P0 numeric anchor | yes |
| B3 | mathematica | 92, 94-99 | sample `mhat²>0` + `Reduce` solvable for mhat>0 | target solvability | yes |
| B4 | mathematica | 107 | `expectZero(N0 at P=0)` | Checks item 4 | yes |
| B5 | mathematica | 110 | `expectZero(residual + target)` at P=0 | Checks item 4 | yes |
| B6 | mathematica | 118 | `expectZero(Δ·D0 − compact denom)` | stability denom identity | yes |
| B7 | mathematica | 122 | `expectZero(N0 − pRaw²/ΔRaw²)` | N0 definition | yes |
| B8 | mathematica | 130-132 | sample `Δ,D0,P0 > 0` | stability gate | yes |
| B9 | mathematica | 150-151 | `expectZero(Limit dP0/dK − D[p0,k])`, X-analog | derivative method cross-check | yes |
| B10 | mathematica | 152-154 | `expectZero(dP0/dK + N0/den²)`, X-analog, sum=0 | monotonicity forms | yes |
| B11 | mathematica | 160-162 | sample `dP0/dK<0`, `dP0/dX>0`, sum=0 | monotonicity signs | yes |

## Findings

None. The attacks below all failed (see Verdict justification). One informational note on output mtime is recorded but does not rise to a finding.

## Independent-derivation check (Mathematica)

The `.wl` is not a transliteration of the `.py`. Three points of genuinely independent derivation:

1. **Derivatives.** SymPy uses symbolic differentiation `sp.diff(P0, K)` / `sp.diff(P0, X)` (sympy 191-192). Mathematica derives the same derivatives by an explicit *quotient limit* `dP0dK = Limit[((p0 /. k -> k + h) - p0)/h, h -> 0]` (wl 141-142) and then *additionally* cross-checks that the limit equals the direct `D[p0,k]` (wl 143-144, asserted at wl 150-151). This is a stronger, independent route, not an echo.
2. **Solvability.** SymPy establishes target solvability only by sample-point positivity of `mhat²` (sympy 132-138). Mathematica adds a symbolic `Reduce[mhat^2 == target/p0Compact && mhat > 0 && delta > 0 && d0 > 0, mhat, Reals]` (wl 94-99) over the full stability domain — a different and more general verification of the same claim.
3. **Algebraic canonicalization.** SymPy reaches the compact form via `sp.simplify` on `P²/(Δ(KΔ − ΔC²/varpi² − Q))` (sympy 103); Mathematica reaches it via `Together` then `Apart[..., k]` (wl 62-63). Different canonicalization machinery, same asserted equality (out line 22 / wl line 23 both = 0).

## Engine cross-check

Both engines produce identical bottom-line results. Δ, Q, P, B0, Z0, N0, D0 match (sympy out 9-15 vs wl out 9-15; Mathematica reports Δ factored as `(omegaU*omegaW − r)(omegaU*omegaW + r)`, algebraically identical to SymPy's `Ω_U²Ω_W² − R²`). `P0 − P0_compact = 0` in both (sympy out 22; wl out 22-23 PASS). Sample anchors identical: `P0 = 1/3`, `Δ = 15`, `D0 = 1/3`, `mhat² = 162/5` (note: SymPy out line 43 prints `162/5`; wl out line 32 prints `162/5` — agree). P=0 corollary: `N0 = 0` and residual `= −target` in both (sympy out 48-51; wl out 40-45). Monotonicity: `dP0/dK = −1`, `dP0/dX = 1` at sample in both (sympy out 97-98; wl out 101-102). All Mathematica checks emit `PASS:`; final banner "Stage 8 Mathematica audit passed." (wl out line 105 — the "Stage 8" label is the old stale self-banner, cosmetic only; see Self-test notes). No residual, sign, or factor disagreement.

## Verdict justification

`clean`. I confirmed I read the paper card, the notes, and the appendix row (part02 row 40), and the scripts' verified claims match the paper's `\stagefield{Output}` deliverables exactly. Attacks that failed:
- **Tautology hunt on the headline `P0 − P0_compact` check.** Not tautological: `P0` is built as `N0/D0` from independently-constructed `N0` and `D0` (sympy 80-81), while `P0_compact` is independently typed as `P²/(Δ(KΔ−ΔC²/varpi²−Q))` (sympy 103). The equality is a real algebraic identity that could fail under a typo in either form. Mathematica reaches the two sides by yet a third path (`Together`/`Apart`). Likewise A6/B6 (`Δ·D0 ≡ compact denom`) and A7/B7 (`N0` reconstructed from raw symbols) are real re-derivations, not `x==x`.
- **Sample-point-only suspicion.** The headline symbolic identities (A1, A6, A7, A9, A10) hold over all symbols, not just the sample; the sample checks (A2, A3, A8, A11) are positivity/sign witnesses, which is the right tool for the strict-inequality stability and monotonicity claims (those cannot be proven by `expect_zero`). Mathematica strengthens the solvability claim to a full-domain `Reduce` (B3). So the strict inequalities are not left at a single point in the `.wl`.
- **Symbol-domain attack.** SymPy declares `K, varpi, OmegaU, OmegaW` positive (sympy 53-58), `C, GU, GW, R` real-unrestricted (correct — couplings/overlaps can be either sign), `X` nonnegative (matches `X=C²/varpi²≥0`), `mhat` positive. These match the paper's stated physical setup (stability is the *claim being tested*, not an assumption baked into the symbol domains: D0's positivity is checked, not assumed). Mathematica's `$Assumptions` (wl 34-35) mirror the same domains. No over-strong assumption hides a branch.
- **Hardcoded-target attack.** The `54 G c_s^5/(5 a^5 c^5)` target is not an unanchored literal: the inline comment (sympy 119-122 / wl 78-83) cites stage023, and I verified upstream that stage023 line 357 anchors `Gamma5_port = a^5/(27 c_s^5)`, line 358 `gamma_GR = 2 G/(5 c^5)`, giving `gamma_GR/Gamma5_port = (2/5)/(1/27) = 54/5` — provenance is real and correctly cited. The notes carry the same value (notes 101/209/236).
- **Variable-independence trap (self-test #1).** `dP0/dX` (sympy 192) differentiates `P0 = N0/(K − X − Q/Δ)`, which genuinely depends on `X`; derivative is nonzero (out: `+N0/den²`), so the sign assertion is non-trivial. Same for `K`.

Informational note (not a finding): the SymPy saved output (mtime 2026-05-26) predates the SymPy script (mtime 2026-06-03). I checked the git diff: the only post-output change to the `.py` was commit e2a4780, a docstring label fix "Stage 8" → "Stage 25" (no logic/value change), and the II.3 corollary block that appears in the output was added in the earlier commit c6ee7d7. The captured output therefore still faithfully reflects what the current script produces; the stale mtime is benign and would not change content on a re-run.

## Value Reconciliation (pass-2 augmentation)

Enumerated every RESULT/deliverable value the two engines emit (from script source + committed `.txt` outputs; nothing executed) and located each in the `.tex` card and/or `.md` notes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Δ = Ω_U²Ω_W² − R²` | py 74 / wl 41; out line 9 | tex:31 (eq:...-invariants); md:61 | MATCH |
| `Q = G_U²Ω_W² + 2G_UG_WR + G_W²Ω_U²` | py 75 / wl 42; out 10 | tex:33; md:63 | MATCH |
| `P = Ω_U²G_W + RG_U` | py 76 / wl 43; out 11 | tex:35 (`\mathcal P`); md:65 | MATCH |
| `B0 = C²/varpi²` | py 78 / wl 44; out 12 | tex:40 (`X=C²/varpi²`); md:69 (`B_0`) | MATCH |
| `Z0 = Q/Δ` | py 79 / wl 45; out 13 | md:71 (`Z_0 = Q/Δ`); tex via D0 | MATCH |
| `N0 = P²/Δ²` | py 80 / wl 46; out 14 | tex:47; md:73 | MATCH |
| `D0 = K − C²/varpi² − Q/Δ` | py 81 / wl 47; out 15 | tex:45 (`D_0=K−X−Q/Δ`); md:77 | MATCH |
| `P0 = P²/[Δ(KΔ − ΔC²/varpi² − Q)]` (compact closed form) | py 103/108 / wl 63/68; out 21-22 | tex:59 (eq:...-P0); md:93 | MATCH |
| Target `N_Q^target = 54 G c_s^5/(5 a^5 c^5)` | py 123 / wl 84; out (target residual block) | md:101, 209, 236; tex:67 carries it as `N_Q^{\rm target}` symbol | MATCH (numeric value in notes; tex uses the symbol) |
| Normalization target eq: `mhat_rad²·P0 = N_Q^target` | py 124 / wl 85 | tex:65-67 (eq:...-normalization-target); md:101 | MATCH |
| Stability gate `Δ>0, D0>0` | py 171-174 / wl 130-131 | tex:75-77 (eq:...-stability); md:111-124 | MATCH |
| `∂P0/∂K = −N0/D0²` (<0) | py 198 / wl 152; out 92,97 | tex:82-83 (eq:...-monotonicity); md:159,165 | MATCH |
| `∂P0/∂X = +N0/D0²` (>0) | py 199 / wl 153; out 93,98 | tex:84 (sign +); md:161,168 | MATCH |
| P=0 corollary `N0=0`, target unreachable | py 145-148 / wl 107-110; out 48-51 | tex:94 (Checks item 4); md:136 (`P0>0` iff `P≠0`) | MATCH |

INTERNAL scaffolding (accounted for, no finding): sample-point substitutions and witnesses — `P0=1/3`, `Δ=15`, `D0=1/3`, `mhat²=162/5`, `dP0/dK=−1`, `dP0/dX=1` (all sample-point numerical witnesses of the symbolic claims, not deliverables); residual-zero PASS flags; `Reduce`/`Limit` intermediate forms; `P0 raw` un-simplified intermediate.

`reconciliation: complete; 14 deliverable values checked, 0 misaligned.` Every emitted deliverable value reconciles to the `.tex` card and/or `.md` notes. The notes are the natural carrier for the explicit numeric target `54 G c_s^5/(5 a^5 c^5)` (the `.tex` legitimately keeps it as the symbol `N_Q^{\rm target}`); that is a MATCH per the augmentation guard (value lives correctly in the `.md`).

## Self-test notes

Checked traps: (1) **Variable independence** — `sp.diff(P0, X)` and `sp.diff(P0, K)` (sympy 191-192) operate on `P0 = N0/(K−X−Q/Δ)`, which depends on both K and X; derivatives are non-vanishing (outputs `∓N0/den²`), so the sign asserts are substantive, not the zero-derivative trap. (2) **Strict inequalities** — stability (`Δ>0,D0>0,P0>0`) and monotonicity signs are correctly tested by sample-point witnesses (SymPy) and, in Mathematica, by a full-domain `Reduce` for solvability; `expect_zero` is reserved for the genuine equalities. (3) **Trivial-case** — P=0 substitution `GW=−R·GU/ΩU²` truly zeroes `P` and hence `N0` (algebraically `Ω_U²·(−R·GU/ΩU²)+R·GU = 0`), so the corollary residual reducing to `−target` is correct. (4) **Symbol domains** match the paper's setup and do not pre-assume the stability being tested. No directive is written (zero findings).
