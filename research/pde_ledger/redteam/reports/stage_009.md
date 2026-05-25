---
unit_id: 009
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 009 red-team report (v2 — paper-grounded re-audit)

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_009.tex`
- notes: (none) — em_projected per-stage notes for steps 004-020 not committed to repo, as noted in the audit prompt
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row for Stage 009 read at line 40)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.txt` (mtime May 21 11:26; script mtime May 21 11:11 — fresh)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.txt` (mtime May 21 11:51; script mtime May 21 11:11 — fresh)

## What the paper claims

The paper card (Stage 009: "Near-throat projection scale", anchor MTDC-T4) distinguishes symmetric interior averaging from one-sided mouth averaging. For a symmetric interior kernel, the first finite-width correction to a smooth profile is "even and begins at O(σ²)". For a mouth kernel of the form `W_ell(s) = (1/ell) w(s/ell)` on the half-line `s ≥ 0` with `∫_0^∞ w(u) du = 1`, the expansion is

```
<X>_ell = X(0) + ell·μ1·X'(0) + O(ell²),    μ1 = ∫_0^∞ u w(u) du.
```

Output verbatim: "Stage~009 explains why mouth-local projected EM data enter at first order in the mouth scale, while symmetric interior thickness effects enter later." The appendix row (line 40) adds: "Near-throat projected scale factors and the leakage/source normalization ledger." Distinct deliverables: (D1) symmetric-kernel expansion is O(σ²) leading correction; (D2) mouth-kernel expansion is `X(0) + ell·μ1·X'(0) + O(ell²)`, i.e. O(ell) leading correction whenever `X'(0) ≠ 0`; (D3) the scale-factor / normalization ledger (effective μ, ξ, source/gauge profile mismatch) referenced by the appendix row.

## What the script claims to verify

The SymPy script (a) writes the exact projected inhomogeneous Maxwell law `∂_μ <Z F^{μν}> + <∂_w(Z F^{wν})> + (1/ξ) ∂^ν <H B> = μ0 <J^ν>` (prose only, no assertion); (b) asserts that for an even kernel acting on a Taylor expansion `Q(w0+σu)` the moments satisfy `<Q>_σ = q0 + (m2 σ²/2) q2 + (m4 σ⁴/24) q4`, with a concrete Gaussian-kernel sanity check fixing `m2=1, m4=3`; (c) for the specific mouth kernel `W_ell(w) = exp(-w/ell)/ell` on `[0,∞)`, computes `<Q>_ell = q0 + ell q1 + ell² q2 + ell³ q3 + ell⁴ q4` and `<∂_w Q>_ell = q1 + ell q2 + ell² q3 + ell³ q4` by direct symbolic integration, and checks IBP recombination `[W Q]_0^∞ − <W' Q> = <∂_w Q>_ell` exactly; (d) derives O(σ²) and O(ell) effective parameter series `μ_eff = μ0 <S>/<Z>`, `ξ_eff = ξ <Z>/<H>` together with H = Z + ε Δh and S = C Z + ε Δs perturbative cancellation; (e) verifies concrete Gaussian-localizer fingerprints `<Z>_sym = 1 − σ²/λ² + 3σ⁴/(2λ⁴)` and `<Z>_mouth = 1 − 2ell²/λ² + 12ell⁴/λ⁴ − 120ell⁶/λ⁶`. The Mathematica script (M1-M5) re-checks the same five clusters.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check(s) | Status |
|---|---|---|
| D1: symmetric-kernel correction is O(σ²) | sympy lines 73-87 (general moment form + Gaussian sanity) + lines 148-156 (μ_eff, ξ_eff symmetric series begin at σ²); Mathematica M2, M4 | match |
| D2: mouth-kernel expansion `X(0)+ell·μ1·X'(0)+O(ell²)` | sympy lines 107-117 (direct integration verifying `<Q>_ell = q0 + ell q1 + ell² q2 + ...`, i.e. the linear coefficient is `q1 = X'(0)` because μ1=1 for `w(u)=e^{-u}`); Mathematica M1, M3 | match (verified for the specific exponential kernel where μ1=1; paper's generic μ1 reduces to 1 for this kernel) |
| D3: scale-factor / normalization ledger (effective μ, ξ, profile-alignment cancellation) | sympy lines 148-214 (μ_eff, ξ_eff series in both regimes + perturbative H=Z, S=CZ exact cancellation + first-order corrections); Mathematica M3, M4 | match |
| (extra) PDE reading of mixed-sector derivative `<∂_w(Z F^{wν})>` survival in near-throat limit | sympy section 1 prose + section 6 prose | extra (no assertion attached; the paper card mentions it implicitly via the kernel-expansion framing; not a finding because the bona fide assertions are anchored elsewhere) |

Dominant pattern: `match`. `paper_alignment: aligned`.

Note on D2: the paper's mouth-expansion identity is written with a generic kernel weight `w(u)` and a generic first moment `μ1`. The script verifies the expansion for one concrete kernel, `w(u) = e^{-u}` (chosen so `W_ell(w) = exp(-w/ell)/ell` is properly normalized on `[0,∞)`), for which `μ1 = ∫_0^∞ u e^{-u} du = 1`. The script's output `q0 + ell·q1 + ell²·q2 + ...` therefore equals `X(0) + ell·μ1·X'(0) + ell²·(higher moment)·X''(0)/...` evaluated at μ1=1. This is the paper's general formula specialized to the kernel actually used. The paper card does not require multi-kernel verification — its load-bearing claim is the *scaling* (O(ell) vs O(σ²)), not the precise numerical value of μ1 across a kernel family. The specialization is therefore a faithful exercise of the paper's claim, not a `partial` or `mismatch`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 86 | `assert_zero("Gaussian even-kernel Q moments", avg_Q_even_gauss − avg_Q_even[m2=1,m4=3])` | D1 (Gaussian instance of even kernel) | yes |
| A2 | sympy | 87 | `assert_zero("Gaussian even-kernel derivative moments", ...)` | D1 | yes |
| A3 | sympy | 111 | `assert_zero("half-line Q expansion", avg_Q_half − (q0 + ell q1 + ell² q2 + ell³ q3 + ell⁴ q4))` | D2 (mouth-kernel expansion) | yes (non-tautological: integrates symbolically, then asserts result equals literal polynomial) |
| A4 | sympy | 116 | `assert_zero("half-line boundary recombination", recombined_half − avg_dQ_half)` | D2 (IBP / boundary cancellation) | yes |
| A5 | sympy | 117 | `assert_zero("half-line derivative expansion", avg_dQ_half − (q1 + ell q2 + ell² q3 + ell³ q4))` | D2 | yes |
| A6 | sympy | 151-152 | `assert_zero("symmetric mu_eff series", ...)` | D1, D3 | yes |
| A7 | sympy | 154-155 | `assert_zero("symmetric xi_eff series", ...)` | D1, D3 | yes |
| A8 | sympy | 161-162 | `assert_zero("mouth mu_eff series", ...)` | D2, D3 | yes |
| A9 | sympy | 164-165 | `assert_zero("mouth xi_eff series", ...)` | D2, D3 | yes |
| A10 | sympy | 203 | `assert_zero("H=Z effective gauge (eps=0 cancellation)", xi_eff_HZ_zero − xi)` | D3 (exact cancellation when H=Z) | yes (non-tautological: integrates concrete profile, divides, simplifies) |
| A11 | sympy | 204 | `assert_zero("S=CZ effective coupling (eps=0 cancellation)", mu_eff_SCZ_zero − C mu0)` | D3 (exact cancellation when S=CZ) | yes |
| A12 | sympy | 211-212 | `assert_zero("xi_eff first-order correction in eps", ...)` | D3 (linear-in-ε correction) | yes |
| A13 | sympy | 213-214 | `assert_zero("mu_eff first-order correction in eps", ...)` | D3 | yes |
| A14 | sympy | 231-232 | `assert_zero("symmetric Gaussian asymptotic literal", series − (1 − σ²/λ² + 3σ⁴/(2λ⁴)))` | D1 (concrete Gaussian fingerprint) | yes |
| A15 | sympy | 245-246 | `assert_zero("mouth Gaussian integral equals erfc closed form", ...)` | D2 (closed-form integral guard) | yes |
| A16 | sympy | 251 | `assert_zero("mouth Gaussian asymptotic from erfc closed form", series − (1 − 2ell²/λ² + 12ell⁴/λ⁴ − 120ell⁶/λ⁶))` | D2 | yes |
| A17 | sympy | 252 | `assert_zero("mouth Gaussian asymptotic from Taylor integration", series − Taylor integral)` | D2 (cross-check series vs Taylor) | yes |
| M1a-c | math | 26-30 | `assertZero["M1a..c", ...]` | D2 | yes |
| M2a-b | math | 43-46 | `assertZero["M2a..b", ...]` | D1 | yes |
| M3a-b | math | 60-63 | `assertZero["M3a..b", ...]` | D2, D3 | yes |
| M4a-b | math | 77-80 | `assertZero["M4a..b", ...]` | D1, D3 | yes |
| M5a-b | math | 91-99 | `assertZero["M5a..b", ...]` | D1, D2 | yes |

All assertion rows are "yes" — each is non-tautological (the LHS is computed by symbolic integration / series expansion, then compared against a literal polynomial; the assertion would fail if the computed expansion produced any different coefficient).

## Findings

None.

### Adversarial attacks attempted (and what failed)

1. **Tautology hunt on the even-kernel moment form (A1, A2).** The script defines `avg_Q_even` as a literal expression rather than deriving it from integration. I checked whether this might be a tautological self-comparison. It is not: A1/A2 compare this literal against an actual SymPy integration of `Gaussian(u) · Q(σu)` and assert agreement when `m2=1, m4=3`. The "literal" form is the abstract moment series; the Gaussian integration is what's being tested. A real check.

2. **Specific-kernel attack on D2.** Could the half-line expansion `<Q>_ell = q0 + ell q1 + ell² q2 + ...` be specific to the exponential kernel and fail to represent the paper's generic claim? For `W_ell(w) = exp(-w/ell)/ell`, the moments are `<w^n> = n! ell^n`, so the literal coefficient of `ell^n` in `<Q>_ell` is `n!·qn/n! = qn`. That matches the paper's `μ1=1` instance. The paper does not require multi-kernel verification (its load-bearing point is the O(ell) scaling vs O(σ²)). No finding.

3. **Sign/symmetry attack on Gaussian asymptotic A16.** The asymptotic series in `ell` for `(1/ell) ∫_0^∞ exp(-w/ell) exp(-w²/λ²) dw` should expand as `1 − 2ell²/λ² + 12ell⁴/λ⁴ − 120ell⁶/λ⁶`. Cross-checked via the moments-of-`e^{-u}` formula: `<w^{2k}> = (2k)! ell^{2k}`, Z's Taylor expansion has coefficients `(−1)^k / (k! λ^{2k})`, so the ell^{2k} coefficient is `(−1)^k (2k)! / (k! λ^{2k})`. For k=0: 1. k=1: −2/λ². k=2: 4!/(2! λ⁴) = 12/λ⁴. k=3: 6!/(3! λ⁶) = 720/6 / λ⁶ = 120/λ⁶ with sign −. Matches A16.

4. **Series-window attack on A14, A16.** `sp.series(..., n).removeO()` truncates; could the literal RHS conceal a missing or wrong-sign coefficient inside the truncation? A14 uses order n=5 in sigma (which captures σ⁰, σ², σ⁴ — three terms). A16 uses asymptotic series at `r→∞` to order 8 in `r=1/ell` (captures ell⁰, ell², ell⁴, ell⁶ — four terms). Both compared term-by-term against literals; passes.

5. **`H=Z` exact-cancellation attack on A10.** I checked whether `xi * IZ_conc / IH_pert |_{eps=0}` actually reduces to `xi` rather than `xi · (some ratio that simplifies to 1)`. `IH_pert |_{eps=0} = ∫ W_conc · Z_conc dw = IZ_conc`, so the ratio is `xi · IZ_conc / IZ_conc = xi`. SymPy's `simplify` reliably handles this cancellation. Real check.

6. **Symbol-assumption attack.** All symbols positive where needed (`ell, sigma, lam, xi, mu0`), generic real otherwise. Mathematica uses parallel assumptions. No domain mismatch with the physical setup.

7. **Independent-derivation attack on Mathematica.** Sections M1, M3, M5 use the same definitions and integrals as the SymPy script. This is not strictly a "transliteration" — both engines independently perform the same symbolic integrations using built-in `Integrate`/`sp.integrate`. There is no algebraic identity being independently re-derived; the test is "does the integral come out to this polynomial in ell?", which is unambiguous and engine-agnostic. (The `mathematica_transliteration` category targets cases where the .wl script literally encodes the same intermediate-step algebra as the .py script. Here both scripts use the same physical premises and the same symbolic integrator — that's not transliteration, that's parallel verification of the same definite integral.) No finding.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally close to the SymPy script in that both use the same Taylor expansion `Qw = q0 + q1 w + q2 w²/2 + …` and the same kernel `W_ell = exp(-w/ell)/ell`, then compare `Integrate[Wel Qw, {w,0,∞}]` against the literal polynomial `q0 + ell q1 + …`. This is parallel symbolic integration in two engines, not algebraic transliteration: each engine computes the integral via its own integrator and returns the simplified result; neither engine echoes the other's algebraic intermediate steps. For Section M3/M4 the `Series[..., {x, 0, n}]` truncation is computed independently by each engine's series routine. The fact that both engines obtain the same simplified residual (zero) is meaningful agreement.

For comparison, two corresponding sections:

- SymPy lines 110-117: `avg_Q_half = sp.simplify(sp.integrate(Wexp*Qw, (w, 0, sp.oo)))` ⇒ `assert_zero("half-line Q expansion", avg_Q_half - (q0 + ell*q1 + ell**2*q2 + ell**3*q3 + ell**4*q4))`.
- Mathematica lines 20-30: `avgQ = Assuming[..., Integrate[Wel Qw, {w, 0, Infinity}]]` ⇒ `assertZero["M1c half-line Q expansion", avgQ - (q0 + ell q1 + ell^2 q2 + ell^3 q3 + ell^4 q4)]`.

These are direct symbolic-integration checks; the engines independently compute the integral and the result is engine-agnostic. No transliteration finding.

## Engine cross-check

| Check (SymPy) | Check (Mathematica) | Both pass? |
|---|---|---|
| half-line Q expansion (A3) | M1c | yes (both residuals = 0) |
| half-line derivative expansion (A5) | M1b | yes |
| half-line IBP recombination (A4) | M1a | yes |
| Gaussian even-kernel Q (A1) | M2a | yes |
| Gaussian even-kernel dQ (A2) | M2b | yes |
| mouth μ_eff series (A8) | M3a | yes |
| mouth ξ_eff series (A9) | M3b | yes |
| symmetric μ_eff series (A6) | M4a | yes |
| symmetric ξ_eff series (A7) | M4b | yes |
| symmetric Gaussian asymptotic (A14) | M5a | yes |
| mouth Gaussian asymptotic (A16) | M5b | yes |

Engines agree on every cross-checked assertion. The H=Z / S=CZ perturbative cancellation block (A10-A13) appears only in the SymPy script; this is a one-engine result but is a direct symbolic computation with no algebraic ambiguity and is consistent with the paper's "scale-factor / normalization ledger" deliverable.

## Verdict justification

Verdict: **clean**, paper_alignment **aligned**.

The paper card states two scaling deliverables (D1: symmetric O(σ²); D2: mouth O(ell) with leading coefficient μ1·X'(0)) and an appendix-level normalization-ledger deliverable (D3: effective μ_eff, ξ_eff with profile-alignment behaviour). Each deliverable maps to multiple non-tautological script-side assertions that pass under both SymPy and Mathematica. The script's choice of the specific exponential kernel `W_ell = exp(-w/ell)/ell` is a faithful concrete instance of the paper's general kernel-class statement (μ1 = 1 for this kernel) and does not understate the paper's claim, because the load-bearing physical content is the O(ell) vs O(σ²) scaling, which the script verifies. I tried tautology hunts on the moment series, sign-and-factor attacks on the Gaussian fingerprints, an exact-cancellation attack on the H=Z / S=CZ block, and a transliteration attack on the Mathematica script; none uncovered a defect. Outputs are fresh (txt mtimes are after the corresponding script mtimes). No findings; no directive required.

This v2 re-audit confirms and supersedes the v1 (script-only) audit. The prior v1 report flagged four findings — all of which were addressed in the May 21 script revisions (Mathematica mirror added; perturbative H=Z and S=CZ checks promoted from substitution to genuine integral-level cancellation; Gaussian asymptotic series now derived from the erfc closed form and cross-checked against Taylor integration). The current scripts hold up under the paper-grounded second pass.

## Self-test notes

Traps checked: (1) variable independence — no `diff` calls in the assertion blocks where the differentiated variable is absent from the expression; the only derivatives are `sp.diff(Qw, w)` and `D[Qw, w]`, both of which depend on `w` as expected. (2) Symmetry/parity — the symmetric kernel verification correctly uses an even kernel and even-power moments; the unit Gaussian gives `<u²>=1, <u⁴>=3` as expected, and odd moments vanish identically (the script's `q1` and `q3` and `q5` terms drop out under the integration in lines 81-87). (3) Trivial-case pre-check — for `q1=q2=q3=q4=0` the half-line integral reduces to `q0` (which matches `q0 + ell·0 + ell²·0 + ...`); for `lam → ∞` the Gaussian asymptotics give `<Z> = 1` (matches `1 − 0 − 0 − ...`). All consistent. No directive written because verdict is clean.
