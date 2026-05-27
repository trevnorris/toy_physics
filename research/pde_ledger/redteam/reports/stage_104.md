---
unit_id: 104
batch: IV.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
  paper_appendix: present
---

# Audit unit 104 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_104.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_104}` at L1242 — no per-stage row beyond the include)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.txt`

## What the paper claims

The stage card (label `stage:104`, titled "Exact Outgoing l=2 DtN Fingerprint") explicitly identifies its audit target as the quoted body equation:
> "Spherical Hankel DtN expansion gives \(\widehat Y_2^{\rm out}=1+z^2/9+4z^4/81+iz^5/27+\cdots\)."

The notes file is the authoritative source of detail. It defines `Λ_2^out(z) = z·h_2^(1)'(z)/h_2^(1)(z)`, gives its small-`z` expansion `-3 + z²/3 + z⁴/9 + iz⁵/9 − 2z⁶/27 − iz⁷/27 + O(z⁸)`, defines the normalized branch `Ŷ_2^out(z) := −3/Λ_2^out(z)` with expansion `1 + z²/9 + 4z⁴/81 + iz⁵/27 − 11z⁶/729 − iz⁷/243 + O(z⁸)`, and the restored-ω form `1 + a²ω²/(9c_s²) + 4a⁴ω⁴/(81c_s⁴) + ia⁵ω⁵/(27c_s⁵) + O(ω⁶)`, plus the canonical odd-coefficient consequence `Γ_{5,can}^{DtN} = a⁵/(27 c_s⁵)`. The card's `Purpose` line names the verification output (the quoted Y expansion) as the audit target; the `Derivation ledger` mention of `m̂_0²·χ_Q·N_Q=1` and `χ_Q=1` is contextual ledger metadata, not a numeric claim required to be tested in this stage.

## What the script claims to verify

Both scripts build the spherical Hankel mode (SymPy: `jn(2,z) + i*yn(2,z)`; Mathematica: `SphericalHankelH1[2,z]`), form `Λ = z·h'/h` and `Y = −3/Λ`, series-expand to O(z⁷), and assert (1) the static slot `Λ(0)+3 = 0`, (2) the five non-trivial coefficients of `Y(z)`: `z²→1/9`, `z⁴→4/81`, `z⁵→i/27`, `z⁶→−11/729`, `z⁷→−i/243`, and (3) after substituting `z=aω/c_s`, the three ω-frame coefficients `a²/(9c_s²)`, `4a⁴/(81c_s⁴)`, `i·a⁵/(27c_s⁵)`. Every assertion in both scripts maps to a coefficient in the notes' two boxed equations and to the body equation of the stage card.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Λ_2^out(0) = −3` (notes L40) | `expect_zero("static DtN slot", lim_{z→0} Λ + 3)` (sympy L37, math L42) | match |
| `Ŷ_2^out` z² coef = 1/9 (notes L57; card L16) | `Y z^2 coefficient` checks `Y.coeff(z,2) − 1/9` (sympy L38, math L43) | match |
| `Ŷ_2^out` z⁴ coef = 4/81 (notes L57; card L16) | `Y z^4 coefficient` checks `Y.coeff(z,4) − 4/81` (sympy L39, math L44) | match |
| `Ŷ_2^out` z⁵ coef = i/27 (notes L57; card L16) | `Y imag z^5 coefficient` checks `Y.coeff(z,5)/i − 1/27` (sympy L40, math L45) | match |
| `Ŷ_2^out` z⁶ coef = −11/729 (notes L57) | `Y z^6 coefficient` checks `Y.coeff(z,6) + 11/729` (sympy L41, math L46) | match |
| `Ŷ_2^out` z⁷ coef = −i/243 (notes L57) | `Y imag z^7 coefficient` checks `Y.coeff(z,7)/i + 1/243` (sympy L42, math L47) | match |
| ω-form coefs `a²/(9c_s²)`, `4a⁴/(81c_s⁴)`, `ia⁵/(27c_s⁵)` (notes L73–80) | three `expect_zero` calls (sympy L50–52, math L57–59) | match |
| `Γ_{5,can}^{DtN} = a⁵/(27 c_s⁵)` (notes L89) | identical to ω⁵ coefficient check above; also restated in script's printed RESULT | match |
| Card Check #1 "product m̂_0²·χ_Q·N_Q keeps source/conservative/outgoing factors separate" | no script check (structural, not numeric) | n/a (structural; addressed upstream) |
| Card Check #2 "higher odd terms begin beyond 2.5PN coefficient" | implicit: the series through z⁷ shows next odd term is iz⁷/243 | partial (no explicit ordering check, but the series exposes it) |
| Card Check #3 "outgoing l=2 DtN fingerprint vs. z=ωa/c_s expansion" | explicit ω-frame checks (sympy L50–52, math L57–59) | match |

Dominant pattern: `match`. Setting `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 37 | `lim(Λ,z,0) + 3 == 0` | static slot `Λ(0)=−3` | yes |
| A2 | sympy | 38 | `Y.coeff(z,2) − 1/9 == 0` | Y z² coef | yes |
| A3 | sympy | 39 | `Y.coeff(z,4) − 4/81 == 0` | Y z⁴ coef | yes |
| A4 | sympy | 40 | `Y.coeff(z,5)/i − 1/27 == 0` | Y z⁵ coef (imag) | yes |
| A5 | sympy | 41 | `Y.coeff(z,6) + 11/729 == 0` | Y z⁶ coef | yes |
| A6 | sympy | 42 | `Y.coeff(z,7)/i + 1/243 == 0` | Y z⁷ coef (imag) | yes |
| A7 | sympy | 50 | `Yw.coeff(ω,2) − a²/(9c_s²) == 0` | ω-frame ω² coef | yes |
| A8 | sympy | 51 | `Yw.coeff(ω,4) − 4a⁴/(81c_s⁴) == 0` | ω-frame ω⁴ coef | yes |
| A9 | sympy | 52 | `Yw.coeff(ω,5)/i − a⁵/(27c_s⁵) == 0` | ω-frame ω⁵ coef = Γ_{5,can} | yes |
| B1 | mathematica | 42 | `Coefficient[lamSeries,z,0] + 3 == 0` | static slot | yes |
| B2 | mathematica | 43 | `Coefficient[ySeries,z,2] − 1/9 == 0` | Y z² coef | yes |
| B3 | mathematica | 44 | `Coefficient[ySeries,z,4] − 4/81 == 0` | Y z⁴ coef | yes |
| B4 | mathematica | 45 | `Coefficient[ySeries,z,5]/I − 1/27 == 0` | Y z⁵ coef (imag) | yes |
| B5 | mathematica | 46 | `Coefficient[ySeries,z,6] + 11/729 == 0` | Y z⁶ coef | yes |
| B6 | mathematica | 47 | `Coefficient[ySeries,z,7]/I + 1/243 == 0` | Y z⁷ coef (imag) | yes |
| B7 | mathematica | 57 | `Coefficient[yOmega,ω,2] − a²/(9c²) == 0` | ω-frame ω² coef | yes |
| B8 | mathematica | 58 | `Coefficient[yOmega,ω,4] − 4a⁴/(81c⁴) == 0` | ω-frame ω⁴ coef | yes |
| B9 | mathematica | 59 | `Coefficient[yOmega,ω,5]/I − a⁵/(27c⁵) == 0` | ω-frame ω⁵ coef = Γ_{5,can} | yes |

Every assertion ties to a specific paper-side claim and is non-tautological: a wrong implementation of `h_2^(1)`, a wrong derivative, or a wrong rationalization would flip at least one of the literal rationals.

## Findings

None.

## Independent-derivation check (Mathematica)

The two scripts are not transliterations. SymPy constructs the Hankel mode via the explicit decomposition `h_2^(1) = j_2 + i·y_2` using `jn(2,z)+I*yn(2,z)` and then `sp.expand_func` to put it into closed elementary form (the printed `h_2^(1)(z)` involves the explicit `(−1/z + 3/z³)·sin(z) − (3/z²)·cos(z)` representation). Mathematica uses the built-in `SphericalHankelH1[2,z]` directly and leaves it abstract until `Series` resolves it. Both then independently compute `z·h'/h`, invert with the static slot, and series-expand. The series engines (`sp.series` vs. `Series`) are independent implementations. Final symbolic forms match exactly in both transcripts:
- Λ: `−3 + z²/3 + z⁴/9 + iz⁵/9 − 2z⁶/27 − iz⁷/27`
- Y: `1 + z²/9 + 4z⁴/81 + iz⁵/27 − 11z⁶/729 − iz⁷/243`

This is convergent agreement of two independent constructions, not transliteration.

## Engine cross-check

Both engines produce identical Λ and Y series and all 18 assertions (9 sympy + 9 mathematica) PASS with explicit `0` residuals. No engine disagreement.

## Verdict justification

`clean`. The stage card body equation, the boxed expansions in the notes, and the canonical odd-coefficient consequence (`a⁵/(27c_s⁵)`) are all faithfully tested by both engines. Attacks tried that failed: (1) checked whether the script defines the target rationals before asserting them — it doesn't; the rationals are literal targets and the series coefficients are computed from a Bessel-function construction, so any algebraic error would flip a sign or factor; (2) checked whether the `z⁵→i/27` check tolerates the wrong branch by dividing out `i` — no, dividing by `sp.I` is exact; if the imaginary part were not `i/27` the residual would be nonzero; (3) checked whether `assume(positive=True)` on `a,c_s` could hide an `Abs` collapse in the ω-frame check — the substitution `z→aω/c_s` is a real-positive monomial substitution, so no branch ambiguity; (4) checked whether the two engines are doing the same computation in disguise — they are not (SymPy decomposes into `jn`+`I·yn`; Mathematica uses `SphericalHankelH1` directly); (5) checked card's `Derivation ledger` references to `m̂_0²·χ_Q·N_Q=1` and `χ_Q=1` for missing-claim coverage — these are not standalone claims tested in this stage; the `Purpose` and quoted body equation explicitly name the Y expansion as the audit target, and the notes contain no such product. Paper, notes, and scripts are coherent and the assertions exercise every load-bearing identity.

## Self-test notes

I checked: (1) variable-independence (every series coefficient is genuinely a function of `z` via `h_2^(1)` derivatives — no zero-by-construction derivatives); (2) parity (the alternation of real/imag coefficients in Λ and Y is consistent with `h_2^(1)` having both real `j_2` and imaginary `y_2` parts — sign and parity of each checked coefficient match the notes); (3) trivial substitution (at `z=0` Λ→−3, so Y=−3/(−3)=1, matching `1` as the constant term); (4) the ω-frame substitution is a positive real monomial, no branch trap. No traps tripped.
