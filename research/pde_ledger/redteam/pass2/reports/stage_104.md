---
unit_id: 104
batch: IV.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md]
  paper_appendix: present
---

# Audit unit 104 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_104.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md` (one file)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read subsec "Canonical outgoing DtN branch", lines 297–343, plus rows 16–26 and 1176/1188/1242)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.txt`

## What the paper claims

The stage derives the exact outgoing spherical `l=2` Dirichlet-to-Neumann (DtN) operator directly from the spherical Hankel mode `h_2^(1)(z) = j_2(z) + i y_2(z)`, with `z := a*omega/c_s`. The operator is `Lambda_2^out(z) = z d/dz ln h_2^(1)(z)`, whose static slot is `Lambda_2^out(0) = -3`, and whose small-z series (notes box, lines 28–38) is `-3 + z^2/3 + z^4/9 + i z^5/9 - 2 z^6/27 - i z^7/27 + O(z^8)`. Normalizing by the static slot gives the canonical outgoing quadrupole admittance (notes box, lines 54–65; card line 16; appendix eq. `Yout-dtn`): `\widehat Y_2^out(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 - 11 z^6/729 - i z^7/243 + O(z^8)`. In ω-form (notes lines 72–81): `1 + a^2 omega^2/(9 c_s^2) + 4 a^4 omega^4/(81 c_s^4) + i a^5 omega^5/(27 c_s^5) + O(omega^6)`. The headline consequence (notes "Consequence", line 89): the canonical odd quadrupole coefficient is fixed, `Gamma_{5,can}^DtN = a^5/(27 c_s^5)`. The card's `Purpose` states "Its audit target is the verification output quoted below," and the quoted block is exactly the `\widehat Y_2^out` fingerprint. The `chi_Q=1` / `m̂_0^2 chi_Q N_Q = 1` framing in the card's "Derivation ledger" and in appendix eq. `chiQ-equals-one` is the downstream comparison that USES this fingerprint (against the retarded one-pole form of an earlier stage); it is narrative, not this script's own deliverable.

## What the script claims to verify

Both scripts build the spherical Hankel mode, form `Lam = z * d/dz(h2)/h2` and `Y = -3/Lam`, series-expand both to O(z^8)/O(z^7), and then assert, against external rational literals taken from the notes box: the static DtN slot (`Lambda(0) = -3`), the `Y` coefficients at `z^2 (1/9)`, `z^4 (4/81)`, `z^5 (i/27)`, `z^6 (-11/729)`, `z^7 (-i/243)`, and after the substitution `z -> a*omega/c_s` the ω-form coefficients at `omega^2 (a^2/9c_s^2)`, `omega^4 (4 a^4/81 c_s^4)`, `omega^5 (i a^5/27 c_s^5)`. The final printed RESULT states the model reproduces the canonical fingerprint including the odd coefficient `a^5/(27 c_s^5)` (= `Gamma_5`).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Static slot `Lambda_2^out(0) = -3` | sympy L37 `limit(Lam,z,0)+3`; wl L42 `Coefficient[lamSeries,z,0]+3` | match |
| `\widehat Y` coefficients z²,z⁴,z⁵,z⁶,z⁷ (notes box / card) | sympy L38–42; wl L43–47 (all five) | match |
| `\widehat Y(ω)` coefficients ω²,ω⁴,ω⁵ (notes box) | sympy L50–52; wl L57–59 | match |
| `Gamma_5 = a^5/(27 c_s^5)` (notes "Consequence") | covered by the ω⁵ assertion (sympy L52 / wl L59); printed RESULT | match |
| Full `Lambda_2^out` series z²/3,z⁴/9,z⁵/9,z⁶,z⁷ (notes box / appendix) | PRINTED (sympy out L11–15; wl out L6); equivalently captured via `Y=-3/Lam` and the asserted Y coeffs | match (printed, asserted indirectly) |
| `chi_Q = 1` / `m̂_0^2 chi_Q N_Q = 1` (card derivation-ledger narrative; appendix) | not tested here (cross-stage comparison object, not this stage's audit target) | not a deliverable of this script |

`paper_alignment: aligned`. The card's stated audit target (the `\widehat Y_2^out` fingerprint expansion) is fully and non-tautologically exercised by both engines. `chi_Q=1` is fixed by comparing this fingerprint to an earlier stage's retarded one-pole form in the appendix narrative; it is not an in-script deliverable of stage 104, so its absence from the script is not `script_missing_paper_claim`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 37 | `limit(Lam,z,0)+3 == 0` | static slot −3 | yes |
| A2 | sympy | 38 | `Y[z^2] − 1/9 == 0` | Y fingerprint | yes |
| A3 | sympy | 39 | `Y[z^4] − 4/81 == 0` | Y fingerprint | yes |
| A4 | sympy | 40 | `Y[z^5]/I − 1/27 == 0` | Y fingerprint (odd) | yes |
| A5 | sympy | 41 | `Y[z^6] + 11/729 == 0` | Y fingerprint | yes |
| A6 | sympy | 42 | `Y[z^7]/I + 1/243 == 0` | Y fingerprint (odd) | yes |
| A7 | sympy | 50 | `Yw[ω^2] − a²/9c_s² == 0` | ω-form / Gamma chain | yes |
| A8 | sympy | 51 | `Yw[ω^4] − 4a⁴/81c_s⁴ == 0` | ω-form | yes |
| A9 | sympy | 52 | `Yw[ω^5]/I − a⁵/27c_s⁵ == 0` | Gamma_5 | yes |
| B1 | mathematica | 42 | `Coefficient[lamSeries,z,0]+3 == 0` | static slot −3 | yes |
| B2 | mathematica | 43 | `Y[z^2] − 1/9 == 0` | Y fingerprint | yes |
| B3 | mathematica | 44 | `Y[z^4] − 4/81 == 0` | Y fingerprint | yes |
| B4 | mathematica | 45 | `Y[z^5]/I − 1/27 == 0` | Y fingerprint (odd) | yes |
| B5 | mathematica | 46 | `Y[z^6] + 11/729 == 0` | Y fingerprint | yes |
| B6 | mathematica | 47 | `Y[z^7]/I + 1/243 == 0` | Y fingerprint (odd) | yes |
| B7 | mathematica | 57 | `Yw[ω^2] − a²/9c² == 0` | ω-form | yes |
| B8 | mathematica | 58 | `Yw[ω^4] − 4a⁴/81c⁴ == 0` | ω-form | yes |
| B9 | mathematica | 59 | `Yw[ω^5]/I − a⁵/27c⁵ == 0` | Gamma_5 | yes |

Every assertion compares an engine-derived series coefficient against an EXTERNAL rational literal lifted from the notes/card box. No coefficient is defined as a symbol and then asserted equal to itself; there are no `x = expr; assert x == expr` round-trips. The imag-coefficient checks (A4/A6/B4/B6) divide the complex coefficient by `I` before subtracting the real literal — a spurious real part in those coefficients would survive as a `real/I` term and fail the zero-check, so they are genuinely falsifiable.

## Findings

None. No `tautological_check`, `hardcoded_result`, `symbol_assumption_error`, `missing_branch`, `engine_disagreement`, `mathematica_transliteration`, `missing_verification_script`, `insufficient_verification`, `stale_output`, or `paper_misalignment` found.

## Independent-derivation check (Mathematica)

The two engines construct the load-bearing primitive — the spherical Hankel function — by **different routes**:

- SymPy (L21): `h2 = sp.expand_func(jn(2, z) + sp.I * yn(2, z))`, i.e. explicitly forms `j_2 + i y_2` and expands each into elementary trig (the saved output L6–9 shows the resulting `sin/cos` rational-coefficient closed form).
- Mathematica (L31): `h2 = FullSimplify[SphericalHankelH1[2, z]]`, i.e. uses the built-in special-function object, which it carries symbolically (saved output L5 prints it unevaluated as `SphericalHankelH1[2, z]`).

The static-slot check is also done by **different mechanisms**: SymPy uses `sp.limit(Lam, z, 0) + 3` (L37) whereas Mathematica reads the `z^0` series coefficient `Coefficient[lamSeries, z, 0] + 3` (L42). The downstream choreography (`Lam = z*D[h2,z]/h2`, `Y = -3/Lam`, series, coefficient extraction, compare-to-external-literal) is necessarily parallel because that IS the single mathematical definition of the DtN operator stated in the notes — there is no second "route" to a DtN operator, and the parallelism here reflects the math, not a transliteration. Crucially, the two series engines (`sp.series` operating on an explicit trig closed form vs. `Series[]` operating on a built-in special function) and the differing static-slot mechanisms mean the `.wl` is NOT a line-by-line re-typing of the `.py`'s numeric algebra. Verdict: **independent** (the primitive and the static-slot mechanism differ; anchors are external literals, not cross-engine echoes). I considered `partial` because the mid-pipeline steps are textually similar, but the divergence at the two ends (primitive construction + slot mechanism) and the external anchoring place it on the independent side of the line.

## Engine cross-check

The engines agree exactly. SymPy output: `Lambda = -3 + z^2/3 + z^4/9 + i z^5/9 - 2 z^6/27 - i z^7/27`; `Y = 1 + z^2/9 + 4 z^4/81 + i z^5/27 - 11 z^6/729 - i z^7/243`; ω-form `1 + a²ω²/9c_s² + 4a⁴ω⁴/81c_s⁴ + i a⁵ω⁵/27c_s⁵`. Mathematica output L6–7,20: identical symbolic forms (with `aThroat/cSound` naming). All nine `expectZero`/`expect_zero` checks return `0` in both transcripts (sympy out L22–27,35–37; wl out L8–26). No residual sign/factor discrepancy.

## Verdict justification

`clean`. I read the card, the notes box, and the appendix "Canonical outgoing DtN branch" subsection before opening the scripts, and the script's verified identity (the exact `\widehat Y_2^out` fingerprint and its ω-form, including `Gamma_5 = a^5/(27 c_s^5)`) is precisely the card's stated audit target. Attacks tried and failed: (1) tautology — assertions anchor to external rational literals from the notes box, not to in-script re-substitutions, so each can fail; (2) imag-part extraction — dividing by `I` is correct and would expose any spurious real part; (3) symbol traps — `a,omega,c_s` declared `positive,real` and `z` real are physically justified (throat radius, sound speed, normalized frequency) and do not over-constrain the series; (4) transliteration — the Hankel primitive and the static-slot mechanism are constructed by genuinely different routes across engines; (5) staleness — both `.txt` outputs are newer than their scripts (sympy out 15:18 vs script 15:08; mma out 15:24 vs script 15:08); (6) paper misalignment — the `chi_Q=1` narrative is a cross-stage comparison, not this script's deliverable, so its absence is correct.

## Self-test notes

Checked: (1) variable-independence — only one `diff`/`D` appears (`d/dz h2`); `h2` genuinely depends on `z`, so the derivative is non-trivial (the printed `Lam` series confirms a nonzero z-dependence), no identically-zero-derivative trap. (2) Parity/symmetry — n/a (series-coefficient extraction on a bounded expansion, no unbounded-domain integral). (3) Trivial-case — the static slot `Lambda(0) = -3` and `Y(0)=1` are the simplest concrete probe and both reduce correctly; the asserted higher coefficients are nonzero literals matched against engine output, confirmed in both transcripts. No directive written (zero findings).

## Value Reconciliation (pass-2 augmentation)

Deliverable values the scripts emit (from source + saved outputs), located in `.tex`/`.md`:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Lambda_2^out(0) = -3` (static slot) | py L37 / wl L42; out: sympy L14 "−3", mma L6 | notes L40 "`Lambda_2^out(0)=-3`"; appendix L312 | MATCH |
| `Lambda_2^out(z) = -3 + z²/3 + z⁴/9 + i z⁵/9 − 2z⁶/27 − i z⁷/27` | py L28/out L11–15; wl L35/out L6 | notes box L28–38; appendix L306–313 (truncated O(z⁶)) | MATCH |
| `Y[z²] = 1/9` | py L38 / wl L43 | notes box L54–65; card L16; appendix L321 | MATCH |
| `Y[z⁴] = 4/81` | py L39 / wl L44 | notes box L54–65; card L16; appendix L321 | MATCH |
| `Y[z⁵] = i/27` | py L40 / wl L45 | notes box L54–65; card L16; appendix L321 | MATCH |
| `Y[z⁶] = −11/729` | py L41 / wl L46 | notes box L54–65 | MATCH |
| `Y[z⁷] = −i/243` | py L42 / wl L47 | notes box L54–65 | MATCH |
| `Y(ω)[ω²] = a²/(9c_s²)` | py L50 / wl L57 | notes box L72–81 | MATCH |
| `Y(ω)[ω⁴] = 4a⁴/(81c_s⁴)` | py L51 / wl L58 | notes box L72–81 | MATCH |
| `Y(ω)[ω⁵] = i a⁵/(27c_s⁵)` | py L52 / wl L59 | notes box L72–81 | MATCH |
| `Gamma_{5,can}^DtN = a⁵/(27 c_s⁵)` | py L55–56 RESULT; wl ω⁵ check | notes L89 `\Gamma_{5,\rm can}^{\rm DtN}=\frac{a^5}{27c_s^5}` | MATCH |

INTERNAL scaffolding (no finding): `h2` symbolic form, `lam`/`Lam` intermediate, `yOut`/`Y` intermediate, the `expectZero/expect_zero` residual `= 0` print lines, PASS/FAIL flags, the banner string.

reconciliation: complete; 11 values checked, 0 misaligned.
