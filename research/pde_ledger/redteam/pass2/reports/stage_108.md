---
unit_id: 108
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
  notes_stage_files: [moving_throat_pde_stage108_robustness_classes.md]
  paper_appendix: present
---

# Audit unit 108 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_108.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage108_robustness_classes.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read §"Outlet DtN and Robin outlet tests" rows 345–397 + the canonical `Yout-dtn` fingerprint rows 305–338; the row at 1250 inputs `stage_108`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage108_robustness_classes_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage108_robustness_classes_mathematica_audit.txt`

## What the paper claims

Stage 108 ("Robustness Classes for chi_Q") classifies which explicit isotropic DtN
deformations of the canonical outgoing l=2 response leave the normalized outgoing
quadrupole fingerprint `chi_Q = 1` invariant and which genuinely shift it. The
card's bottom line (the boxed quote) is: "Pure scale is harmless; pure argument
shift is harmless only if even moments remain canonical; additive throat-core data
can move chi_Q." The notes enumerate the deliverables: Class A (pure scale
`Λ_def = SΛ_out` ⇒ `chi_Q=1`), Class B (scale+argument `Λ_def = SΛ_out(βz)` ⇒
even-fingerprint preservation forces `β=1`, hence `chi_Q=1`; with `chi_arg=β⁵`),
Class C (additive core, `β=1`, even-match forces `Σ₂=−Σ₀/9, Σ₄=−Σ₀/27`, giving the
boxed `chi_Q = 3(S+9Σ₅)/(3S−Σ₀)`), and the exact preservation submanifold
`Σ₅ = S(1−β⁵)/9 − Σ₀/27`. The appendix (eq:app-part04-chiQ-general,
eq:app-part04-Ydef-expanded) states the general result `chi_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)`
after imposing the canonical even fingerprint, and confirms the canonical
fingerprint `Ŷ_out = 1 + z²/9 + 4z⁴/81 + i z⁵/27` (eq:app-part04-Yout-dtn).

## What the script claims to verify

Both engines build `Λ_out = -3 + z²/3 + z⁴/9 + i z⁵/9` and run four classes.
Class A: the normalized scaled response `-3S/(SΛ_out)` series-expanded must equal the
literal canonical fingerprint `1 + z²/9 + 4z⁴/81 + i z⁵/27` (a falsifiable external
anchor, not a bare S-cancellation). Class B: `Ŷ(βz)` has `m2=β²/9, m4=4β⁴/81,
chi_arg=β⁵`; the even-match `{β²=1, β⁴=1}` yields the two real roots `{−1,+1}`, and
`chi_arg(β=1)=1`. Class C: from the additive `Λ_add`, even-match solves
`Σ₂=−Σ₀/9, Σ₄=−Σ₀/27`, then `chi_add = 3(S+9Σ₅)/(3S−Σ₀)`, with the `chi_Q=1`
preservation locus `Σ₅=−Σ₀/27` anchored against the notes box (not re-substituted).
Class D (general): same with the β-argument retained, giving
`chi_gen = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)` and the submanifold `Σ₅ = S(1−β⁵)/9 − Σ₀/27`,
reducing to Class C at `β=1`. The `.wl` additionally cross-checks `chi_gen` against
a raw closed-form-coefficient route (no `Series`).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Canonical fingerprint `Ŷ_out = 1+z²/9+4z⁴/81+i z⁵/27` (eq:app-part04-Yout-dtn) | py L31–33 / wl L35–38 anchor scaled response to literal canonical | match |
| Class A: pure scale ⇒ `chi_Q=1` (notes box) | py L31–33 / wl L35–38 (S cancels, equals canonical) | match |
| Class B: `chi_arg=β⁵`, even-match ⇒ `β=±1` ⇒ `β=1` ⇒ `chi_Q=1` (notes) | py L36–49 / wl L40–56 (incl. coeff anchors m2,m4,chi in .wl) | match |
| Class C: `Σ₂=−Σ₀/9, Σ₄=−Σ₀/27` (notes) | py L59–64 / wl L66–73 | match |
| Class C box: `chi_Q = 3(S+9Σ₅)/(3S−Σ₀)` (notes) | py L65–67 / wl L75–77 | match |
| Class C preservation locus `Σ₅=−Σ₀/27` (notes special case) | py L69–73 / wl L79–81 (anchored to notes, non-tautological) | match |
| General `chi_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)` (eq:app-part04-chiQ-general) | py L77–98 / wl L85–110 (Class D) | match |
| Exact submanifold `Σ₅ = S(1−β⁵)/9 − Σ₀/27` (notes box) | py L99–110 / wl L111–121 | match |
| Linearized `Δ_Q = 5b + a₀/3 + 9a₅` (eq:app-part04-linear-DeltaQ-triple) | implied by the exact `chi_gen` the script verifies; not a separate deliverable | match (downstream consequence) |

`paper_alignment: aligned` — every paper-side and notes-side boxed/labeled result has a
faithful, non-tautological script-side counterpart in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 33 | `expect_zero(Y_scale − Y_can_literal)` | Class A / canonical fingerprint | yes |
| A2 | sympy | 47–48 | `{β roots} == {−1,1}` (raise on mismatch) | Class B even-match roots | yes |
| A3 | sympy | 49 | `expect_zero(chi_arg(β=1) − 1)` | Class B chi_Q=1 | yes |
| A4 | sympy | 60–61 | unique additive even-match solution | Class C Σ₂,Σ₄ | yes |
| A5 | sympy | 67 | `expect_zero(chi_add − 3(S+9Σ₅)/(3S−Σ₀))` | Class C box | yes |
| A6 | sympy | 73 | `expect_zero(Σ₅ locus − (−Σ₀/27))` | Class C preservation locus | yes (notes-anchored) |
| A7 | sympy | 101–104 | `expect_zero(Σ₅_gen − [S(1−β⁵)/9 − Σ₀/27])` | exact submanifold | yes |
| A8 | sympy | 107–110 | `expect_zero((Σ₅_gen + Σ₀/27)|β=1)` | Class C reduction of Class D | yes |
| A9 | math | 38 | `expectZero[yScale − yCanLiteral]` | Class A | yes |
| A10 | math | 49–51 | `expectZero[m2−β²/9], [m4−4β⁴/81], [chi−β⁵]` | Class B coeff closed forms | yes |
| A11 | math | 55–56 | β roots `{−1,1}` + `chi_arg(β=1)−1` | Class B | yes |
| A12 | math | 72–73 | `expectZero[Σ₂+σ₀/9], [Σ₄+σ₀/27]` | Class C even-match | yes |
| A13 | math | 77 | `expectZero[chiAdd − 3(sNorm+9σ₅)/(3sNorm−σ₀)]` | Class C box | yes |
| A14 | math | 81 | `expectZero[σ₅Pres − (−σ₀/27)]` | Class C locus | yes |
| A15 | math | 110 | `expectZero[chiGenAlt − chiGen]` (raw-coeff vs Series) | independent route cross-check | yes |
| A16 | math | 113–116 | `expectZero[σ₅PresGen − (sNorm(1−β⁵)/9 − σ₀/27)]` | exact submanifold | yes |
| A17 | math | 118–121 | `expectZero[(σ₅PresGen+σ₀/27)|β=1]` | Class C reduction | yes |

All load-bearing assertions trace to a specific paper/notes deliverable and are
non-tautological (each compares a script-derived quantity to an externally stated
closed form, not to itself).

## Findings

None. The script holds up against every attack tried (see Verdict justification).

## Independent-derivation check (Mathematica)

The `.wl` is **partially independent** — independent on the load-bearing general
result, parallel on the Class A/B Series extraction:

- **Class A/B (parallel / transliterated route).** Both engines do
  `Series[(-3*sNorm)/(sNorm*lambdaOut /. z->beta z), {z,0,5}]` (wl L40) mirroring
  `sp.series((-3*S)/(S*Lambda_out.subs(z, beta*z)), z, 0, 6)` (py L36), with the
  same coefficient extraction (`Coefficient[yArg,z,n]` vs `.coeff(z,n)`). However
  the `.wl` is NOT a bare echo here: it adds three substantive closed-form anchors
  the `.py` lacks — `expectZero["m2_arg - beta^2/9"]`, `["m4_arg - 4 beta^4/81"]`,
  `["chi_arg - beta^5"]` (wl L49–51) — so the second engine does strictly more
  verification than the first on this class.
- **Class C/D (genuinely independent route).** The `.wl` re-derives `chiGen` a
  *second* way, from hand-written closed-form L-coefficients with **no `Series`**:
  `L0raw = -3*sNorm + sigma0; L2raw = sNorm*beta^2/3 + sigma2; L4raw =
  sNorm*beta^4/9 + sigma4; L5raw = sNorm*beta^5/9 + sigma5;` (wl L103–106), solves
  the even-match, and asserts `chiGenAlt - chiGen == 0` (wl L110). I verified these
  raw coefficients by hand against `Λ_out(βz) = -3 + β²z²/3 + β⁴z⁴/9 + iβ⁵z⁵/9`
  scaled by S plus the additive slots — all four are correct. This route does not
  exist in the `.py` and is independent of the Series path.

Because the load-bearing general result (`chi_gen` / `chi_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)`)
is verified by a route in the `.wl` that the `.py` does not contain, this is NOT a
`mathematica_transliteration` finding under the second-engine policy. The parallel
Class A/B Series choreography is a low-information convenience step (extracting a
literal fixed truncation that both engines then anchor to the same external
appendix fingerprint), not the load-bearing derivation. Verdict: **partial**, with
the load-bearing result independently re-derived → no finding.

## Engine cross-check

Final outputs agree at the level claimed:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| pure scale invariance residual | `0` | `0` (PASS) |
| `Y_scale+arg` | `1 + β²z²/9 + 4β⁴z⁴/81 + i β⁵z⁵/27` | `1 + β²z²/9 + 4β⁴z⁴/81 + (I/27)β⁵z⁵` |
| `m2_arg`, `m4_arg`, `chi_arg` | `β²/9`, `4β⁴/81`, `β⁵` | `β²/9`, `4β⁴/81`, `β⁵` |
| β roots | `{−1, 1}` | `{−1, 1}` |
| `Σ₂(β=1)`, `Σ₄(β=1)` | `−Σ₀/9`, `−Σ₀/27` | residuals `+σ₀/9`, `+σ₀/27` → 0 (PASS) |
| `chi_add` | `3(S+9Σ₅)/(3S−Σ₀)` | `−3(9σ₅+sNorm)/(σ₀−3sNorm)` (≡ same) |
| `Σ₅` locus (Class C) | `−Σ₀/27` | `−σ₀/27` |
| `chi_gen(β)` | `3(Sβ⁵+9Σ₅)/(3S−Σ₀)` | `−3(9σ₅+β⁵sNorm)/(σ₀−3sNorm)` (≡ same) |
| `Σ₅` general locus | `S(1−β⁵)/9 − Σ₀/27` | `(−σ₀ − 3(−1+β⁵)sNorm)/27` (≡ same) |
| submanifold residual | `0` | `0` (PASS) |

Both exit 0. No `engine_disagreement`.

## Verdict justification

**clean.** I read the paper card, the notes (Class A/B/C + exact preservation
submanifold), and the part-04 appendix (general deformation algebra
eq:app-part04-chiQ-general, the canonical fingerprint eq:app-part04-Yout-dtn) before
opening the scripts. The pass-1 `F1 paper_misalignment` (the carried-forward
`chi_Q`/robustness result) is resolved and still holds: the script's literal
canonical fingerprint matches the appendix `Yout-dtn` byte-for-byte, and the boxed
`chi_Q = 3(S+9Σ₅)/(3S−Σ₀)` / general `3(Sβ⁵+9Σ₅)/(3S−Σ₀)` match the notes and
appendix exactly. Attacks tried and failed: (1) **tautology** — the Class C/D
preservation loci are solved from independently-built `chi` expressions then anchored
to the *notes'* closed form, and the round-trip self-check is explicitly demoted to a
print (py L106), so the load-bearing tests are non-tautological; (2) **bare
S-cancellation in Class A** — defeated by the literal external fingerprint anchor;
(3) **single-point / degenerate-parameter test** — the general `chi_gen(β)` retains
β symbolically (not fixed to a static value), and Class B exhibits the full β
dependence via the closed-form coefficient anchors; (4) **symbol-assumption trap** —
the only nontrivial assumption `3*sNorm - sigma0 != 0` (wl L31) is the physical
nondegeneracy of the `chi_Q` denominator and does not contradict the setup; the
`.py` needs no such assumption since it carries the rational function exactly; (5)
**transliteration** — the load-bearing general result is independently re-derived in
the `.wl` via a raw-coefficient route the `.py` lacks (wl L103–110). Outputs are
fresh (both `.txt` mtimes 00:21 > scripts 00:16–00:17) and their content matches what
the current scripts produce. No findings.

## Value Reconciliation (pass-2 augmentation)

Enumerating every RESULT/deliverable value the scripts emit and locating each in the
docs. Output `.txt` files present and fresh; reconciliation based on script source +
saved outputs.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Ŷ_scale` = `1 + z²/9 + 4z⁴/81 + i z⁵/27` (canonical fingerprint anchor) | py L32 / wl L37; py out L5, wl out L5 | appendix eq:app-part04-Yout-dtn (`stage_appendix_part04.tex:317–321`) | MATCH |
| `Ŷ_scale+arg` = `1 + β²z²/9 + 4β⁴z⁴/81 + i β⁵z⁵/27` | py out L6, wl out L7 | notes Class B (`...stage108...md:31–33`) | MATCH |
| `m2_arg = β²/9` | py out L7, wl out L8 | notes Class B (md:31–33, coeff of z²) | MATCH |
| `m4_arg = 4β⁴/81` | py out L8, wl out L9 | notes Class B (md:31–33, coeff of z⁴) | MATCH |
| `chi_arg = β⁵` | py out L9, wl out L10 | notes Class B (md:33); appendix `chi_Q` ∝ `β⁵` (`tex:373`) | MATCH |
| β roots `{−1, 1}` | py out L10, wl out L17 | notes Class B "β²=1, β⁴=1 ⇒ β=1" (md:38–44) | MATCH |
| `Σ₂(β=1) = −Σ₀/9` | py out L12, wl out L20/L22 | notes Class C (`md:62–65`) | MATCH |
| `Σ₄(β=1) = −Σ₀/27` | py out L13, wl out L24 | notes Class C (`md:62–65`) | MATCH |
| `chi_add = 3(S+9Σ₅)/(3S−Σ₀)` | py out L14, wl out L26 | notes Class C box (`md:67–70`) | MATCH |
| `Σ₅` locus (Class C) `= −Σ₀/27` | py out L16, wl out L29 | notes submanifold β=1 reduction (`md:90–92`) | MATCH |
| `Σ₂(β) = −Sβ²/3 + S/3 − Σ₀/9` | py out L18, wl out L32 | notes general algebra (consistent w/ md Class C at β=1); appendix general deformation (`tex:351–366`) | MATCH |
| `Σ₄(β) = −Sβ⁴/9 + S/9 − Σ₀/27` | py out L19, wl out L33 | as above | MATCH |
| `chi_gen(β) = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)` | py out L20, wl out L34 | appendix eq:app-part04-chiQ-general (`tex:370–373`) | MATCH |
| `Σ₅` general locus `= S(1−β⁵)/9 − Σ₀/27` | py out L21, wl out L37 | notes exact submanifold box (`md:89–92`) | MATCH |

Internal scaffolding (accounted for, no finding): `expect_zero`/`expectZero`
residuals (all `0`), `PASS:` flags, the demoted "general preservation locus check =
0" round-trip print, the "chiGen independent-route agreement = 0" cross-check
residual, and `# exit_code: 0`.

reconciliation: complete; 14 deliverable values checked, 0 misaligned

## Self-test notes

Checked: (1) **Variable independence** — `chi_gen` genuinely depends on β (output
`3(Sβ⁵+9Σ₅)/(3S−Σ₀)`), so the β=1 reduction and submanifold tests are not
degenerate-zero; the raw-coefficient route's L-coefficients each carry their stated β
power. (2) **Symbol assumptions** — `3*sNorm−sigma0 != 0` is the physical
nondegeneracy of the denominator, consistent with the setup; no positivity trap. (3)
**Trivial-case** — substituting the Class C locus `Σ₅=−Σ₀/27` into `chi_add` gives 1
(the demoted print confirms), and the anchors compare to externally-stated notes
forms rather than self, so no assertion passes by construction. No directive written
(zero findings).
