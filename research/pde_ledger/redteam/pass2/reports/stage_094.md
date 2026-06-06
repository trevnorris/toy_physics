---
unit_id: 094
batch: IV.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage094_isotropic_geometry_decoupling.md]
  paper_appendix: present
---

# Audit unit 094 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_094.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage094_isotropic_geometry_decoupling.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read intro + audit-path map + MTDC-T8 result-anchor rows; stage row at line 1222 `\input{stages/stage_094}`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.txt`

## What the paper claims

The stage proves an **isotropic geometry-decoupling theorem**: on an isotropic reference throat the scalar/geometry `l=0` lane and the grouped real `l=2` wall modes are exactly orthogonal in the quadratic wall theory, so the quadratic action is block-diagonal in `l` and there is no `l=0 <-> l=2` bilinear mixing (`M_{0<->2}=0`). The card's `\stagefield{Output}`-equivalent body line (stage_094.tex:16) states verbatim: "Orthogonality of \(l=0\) geometry and grouped \(l=2\) wall modes gives \(K_{g,2}=K_{g,4}=0\) on the isotropic branch." The derivation-ledger line (stage_094.tex:13) names the forced conservative carrier `\widehat Y_Q^{cons} = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)` and the obstruction variables `(eps_2, eps_4)`. The notes (Sections 1-3) make the mechanism explicit: every `l=0/l=2` bilinear reduces to one of `∫ Y00 Y2A`, `∫ grad Y00 · grad Y2A`, or `∫ Y00 (−Δ)Y2A`, all of which vanish because (i) `Y00 ⟂ l=2`, (ii) `grad Y00 = 0`, and (iii) `(−Δ)Y2A = 6 Y2A` reduces the third to the first. Consequently `K_{g,2}=K_{g,4}=0`, so `eps_2 = Omega_Q^2 K_{g,2}/K_pole = 0`, `eps_4 = Omega_Q^4 K_{g,4}/K_pole = 0`, and the obstruction collapses to `c_pole = 1/4`, `c_geom = 3/4`. The three named Checks are: (1) static limit `eps_2=eps_4=0` returns `c_pole=1/4`; (2) `l=0`/`l=2` orthogonality before the firewall; (3) success statements still carry the minimal-module hypothesis.

Distinct deliverables: (D1) the angular orthogonality integrals vanish; (D2) the `(−Δ)Y2A = 6 Y2A` eigenvalue identity; (D3) the generic isotropic cross coefficient between `l=0` and each `l=2` mode vanishes (block diagonality); (D4) `K_{g,2}=K_{g,4}=0 ⇒ eps_2=eps_4=0`; (D5) `c_pole=1/4`, `c_geom=3/4`.

## What the script claims to verify

Both scripts construct the five explicit real `l=2` harmonics plus `Y00`, integrate over `S^2` with the proper `sin(th)` measure, and assert for each `l=2` mode: normalization `=1`; `<Y00|Y2A> = 0` (overlap); `(−Δ)Y2A − 6 Y2A = 0` (eigenvalue); `<grad Y00 · grad Y2A> = 0` (gradient cross); `<Y00|(−Δ)Y2A> = 0` (Laplacian-weighted overlap); and a generic isotropic cross coefficient `C_{0,A} = mu·I_mass − Tw·I_mass − TOm·I_lap − K·I_pot = 0` built as a symbolic linear combination of those (already-vanishing) integrals. A final "static-limit" block sets `K_g2 = K_g4 = 0` literally, forms `eps_2 = Omega_Q^2·K_g2/K_pole`, `eps_4 = Omega_Q^4·K_g4/K_pole`, asserts both are zero, then sets `c_pole=1/4`, `c_geom=3/4` and asserts `c_pole + c_geom == 1`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 angular orthogonality `∫ Y00 Y2A = 0` | `domega(Y00*Y)`/`dOmega[y00*y]` per mode (py:40, wl:63) | match |
| D2 `(−Δ)Y2A = 6 Y2A` | `−lap_s2(Y) − 6Y` / `−lapS2[y] − 6 y` (py:36, wl:64) | match |
| D2' `∫ Y00 (−Δ)Y2A = 0` | `domega(Y00*(-lap_s2(Y)))` (py:42, wl:66) | match |
| D1' `∫ grad Y00 · grad Y2A = 0` | `grad_cross` (py:41, wl:65) — integrand identically 0 since grad Y00 = 0 | partial (mirrors paper's own argument; integral does no work) |
| D3 block diagonality `M_{0<->2}=0` | generic `C_{0,A}=0` (py:54, wl:67) | match (linear combo of D1/D2' integrals) |
| D4 `eps_2=eps_4=0` | py:62-65, wl:77-82 — `K_g2,K_g4` hardcoded to 0 | mismatch (tautological; not derived from D1) |
| D5 `c_pole=1/4`, `c_geom=3/4` | py:66-68, wl:79-83 — hardcoded, assert is only on the SUM `=1` | partial (values assigned, not asserted individually) |

Dominant pattern: the load-bearing main theorem (D1, D2, D2', D3) is faithfully and non-tautologically exercised, and every emitted value reconciles with the docs. The static-limit bookkeeping block (D4/D5, paper Check #1) is hollow. `paper_alignment: aligned` (no value or target mismatch; the lone defect is a script-internal tautology, not a paper↔script disagreement).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 40 | `domega(Y00*Y) == 0` | D1 | yes |
| A2 | sympy | 41 | `domega(grad_cross) == 0` | D1' | partial (integrand 0 by construction) |
| A3 | sympy | 42 | `domega(Y00*(-lap_s2(Y))) == 0` | D2' | yes |
| A4 | sympy | 36 (print) | `-lap_s2(Y) - 6*Y` printed 0 | D2 | yes (printed, not asserted — see note) |
| A5 | sympy | 54 | `simplify(Ccross) == 0` | D3 | yes |
| A6 | sympy | 64-65 | `eps_2 == 0`, `eps_4 == 0` | D4 | no (K_g2,K_g4 set to 0) |
| A7 | sympy | 68 | `c_pole + c_geom == 1` | D5 | no (1/4+3/4 trivially 1) |
| A8 | mathematica | 50,62 | `expectZero["...norm-1"]` | normalization | yes |
| A9 | mathematica | 63 | `expectZero["<Y00|Y..>"]` | D1 | yes |
| A10 | mathematica | 64 | `expectZero["(-Delta)Y.. - 6 Y.."]` | D2 | yes |
| A11 | mathematica | 65 | `expectZero["<grad Y00 . grad Y..>"]` | D1' | partial (integrand 0 by construction) |
| A12 | mathematica | 66 | `expectZero["<Y00|(-Delta)Y..>"]` | D2' | yes |
| A13 | mathematica | 67 | `expectZero["C_0,.."]` | D3 | yes |
| A14 | mathematica | 81-82 | `expectZero["eps_2/eps_4 (static)"]` | D4 | no (Kg2,Kg4 = 0) |
| A15 | mathematica | 83 | `expectZero["c_pole + c_geom - 1"]` | D5 | no (1/4+3/4 trivially 1) |

Note on A4: the SymPy eigenvalue residual is only printed, not `assert`ed; the Mathematica side (A10) does assert it via `expectZero`. The cross-engine assertion on A3/A12 (`<Y00|(−Δ)Y2A>`) implicitly exercises the eigenvalue identity in SymPy too, so D2 is covered overall — minor, not a finding.

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py:60-68`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl:73-83`

**What's wrong:**
The static-limit block — which the script's own comment ties to "paper Check #1" — is tautological. SymPy:
```
K_g2 = sp.Integer(0)   # line 60
K_g4 = sp.Integer(0)   # line 61
eps_2 = sp.simplify(Omega_Q**2 * K_g2 / K_pole)   # = 0 by construction
eps_4 = sp.simplify(Omega_Q**4 * K_g4 / K_pole)   # = 0 by construction
assert eps_2 == 0   # line 64 — cannot fail
assert eps_4 == 0   # line 65 — cannot fail
c_pole = sp.Rational(1, 4)
c_geom = sp.Rational(3, 4)
assert c_pole + c_geom == 1   # line 68 — checks 1/4 + 3/4 == 1, vacuous
```
Mathematica is identical (`Kg2 = 0; Kg4 = 0; ... expectZero["eps_2 (static limit)", eps2]; ... expectZero["c_pole + c_geom - 1", ...]`, wl:73-83). The paper's Check #1 is "Check the static limit `eps_2=eps_4=0` returns `c_pole=1/4`." The script proves none of this: `K_{g,2}`/`K_{g,4}` are *assigned* zero rather than *derived* from the orthogonality integrals computed in the block above, so `eps_2=eps_4=0` is true by definition and the assertions can never fail regardless of whether the orthogonality theorem holds. The `c_pole + c_geom == 1` assertion tests `1/4 + 3/4 == 1`, an arithmetic identity independent of the physics; it does not even pin `c_pole = 1/4` (it would pass for `c_pole=0.3, c_geom=0.7`), even though the named deliverable is specifically `c_pole = 1/4`.

The substantive theorem (the angular integrals A1, A3, A5 / A9, A12, A13) IS real and genuinely shows `K_{g,2}=K_{g,4}=0`; the defect is that the Check-#1 block re-hardcodes that conclusion instead of consuming it, so the named "static limit returns c_pole=1/4" check has no teeth.

**Why this matters:**
Paper Check #1 is a named verification obligation of this stage. As written, the script would still print "PASS"/"eps_2 = 0; c_pole = 1/4" even if the orthogonality block above were deleted or wrong — the two blocks are disconnected. A reviewer trusting the transcript would believe Check #1 was exercised when it was not.

**Required change:**
Make the static-limit block *consume* the orthogonality result rather than re-hardcode it. Concretely, derive `K_g2`/`K_g4` from the actually-computed angular cross integral so a nonzero overlap would propagate into a nonzero `eps`. For example, accumulate the worst-case (maximum-magnitude) `l=0/l=2` overlap integral across the five modes computed above into a symbol `K_g_overlap`, then set `K_g2 = Omega_Q**2 * K_g_overlap` style so that `eps_2` is `simplify`d FROM the integral, and assert `eps_2 == 0` only because that integral was proven 0 — i.e., `K_g2 = sum(domega(Y00*Y) for Y in Y2.values())` (this is the orthogonality sum, provably 0, but NOT typed as a literal 0). Then assert `c_pole == sp.Rational(1,4)` and `c_geom == sp.Rational(3,4)` individually (not just the sum), matching the named deliverable `c_pole=1/4`. Mirror identically in the `.wl` (`Kg2 = Sum of dOmega[y00*y]`, `expectZero` on the derived `eps2`, and an explicit `c_pole - 1/4` and `c_geom - 3/4` check). See directive for exact edits.

**Verification:**
After the fix, the SymPy block must contain `K_g2`/`K_g4` defined as a `domega(...)`-derived expression (no bare `sp.Integer(0)` assignment of `K_g2`), and the `.wl` block likewise via `dOmega[...]`. New/changed assertions: `eps_2 == 0`, `eps_4 == 0`, `c_pole == 1/4`, `c_geom == 3/4` (individually). Both scripts still exit 0; transcript still prints `eps_2 = 0; eps_4 = 0; c_pole = 1/4; c_geom = 3/4`.

## Independent-derivation check (Mathematica)

The `.wl` is structurally parallel to the `.py` (same `Y00`/`Y2A` definitions, same `dOmega`/`lapS2` operators, same per-mode loop). This parallelism is by necessity, not transliteration: both engines must use the *same* physical givens — the standard real `l=2` harmonics, the `S^2` measure `sin(th) dth dph`, and the angular Laplacian `(1/sin)∂_th(sin ∂_th) + (1/sin^2)∂_ph^2`. There is no intermediate *algebraic choreography* being echoed; each engine independently performs its own symbolic integration and `FullSimplify`/`simplify` of standard harmonics. Compare wl:36-39 (`lapS2`) vs py:20-21 (`lap_s2`) — same operator, independently coded in each language's syntax; wl:31-34 (`dOmega`) vs py:16-17 (`domega`) — same measure, independent integration calls. This does NOT rise to a `mathematica_transliteration` finding: the policy targets echoed derivation steps, not shared physical definitions. (The static-limit block is shared *and* tautological, but that is captured by F1 in both engines, not as a transliteration finding.)

## Engine cross-check

The two engines agree exactly. SymPy output (`...sympy_audit.txt`): every `<Y00|Y2A> = 0`, `(-Delta)Y2A - 6Y2A = 0`, `<grad Y00 . grad Y2A> = 0`, `<Y00|(-Delta)Y2A> = 0`, all five `C_0,A = 0`, and `eps_2 = 0; eps_4 = 0; c_pole = 1/4; c_geom = 3/4`. Mathematica output (`...mathematica_audit.txt`): identical residuals (every line `= 0` / `PASS`), normalizations `= 1`, `eps_2 = 0; eps_4 = 0; c_pole = 1/4; c_geom = 3/4`. No sign, factor, or residual disagreement. Both transcripts are fresh (output mtimes 1779913715 / 1779913764 are newer than script mtimes 1779902036 / 1779902045). `engines_agree: true`, `outputs_fresh: true`.

## Verdict justification

The main result — exact `l=0 <-> l=2` angular decoupling on the isotropic branch, hence `K_{g,2}=K_{g,4}=0` — is genuinely and non-tautologically verified by both engines via explicit `S^2` integrals of the standard harmonics, with the eigenvalue identity `(−Δ)Y2A = 6 Y2A` independently confirmed; the two engines agree exactly and all emitted values reconcile with the card and notes. Attacks that failed: I checked whether the orthogonality integrals were hardcoded (they are not — they are real `integrate` calls that could return nonzero), whether the `6` eigenvalue was wrong (it is `l(l+1)=6`, correct), whether `Y00` normalization or the harmonic forms were off (normalizations all return 1), and whether the engines secretly disagree (they do not). The one real defect is the static-limit block (paper Check #1): it sets `K_{g,2}=K_{g,4}=0` literally and then asserts the resulting `eps`'s are zero, plus a vacuous `1/4 + 3/4 == 1` — a `tautological_check` that disconnects Check #1 from the theorem it is supposed to consume. This is fixable script-side (route `K_g2/K_g4` through the computed overlap integrals; assert `c_pole=1/4`, `c_geom=3/4` individually). Hence `verdict: findings`, one finding, medium severity, no stop-cold.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `<Y00|Y2A> = 0` (all 5 modes) | py:40, wl:63; sympy out L3,8,13,18,23; wl out L11,24,37,50,63 | notes:55-66 (orthogonality argument); card:16 (`K_{g,2}=K_{g,4}=0`) | MATCH |
| `(−Δ)Y2A − 6 Y2A = 0` (all 5) | py:36, wl:64; sympy out L4,9,14,19,24 | notes:69 (`(−Δ)Y2A = 6 Y2A`) | MATCH |
| `<grad Y00 · grad Y2A> = 0` (all 5) | py:41, wl:65; sympy out L5,10,15,20,25 | notes:68 (`grad Y00 = 0`) | MATCH |
| `<Y00|(−Δ)Y2A> = 0` (all 5) | py:42, wl:66; sympy out L6,11,16,21,26 | notes:63-69 | MATCH |
| `C_0,A = 0` (all 5, block diagonality) | py:54, wl:67; sympy out L27-31 | notes:83-89 (`M_{0<->2}=0`) | MATCH |
| `eps_2 = 0` | py:64,69, wl:81; sympy out L32 | notes:31,114,157; card:22 | MATCH |
| `eps_4 = 0` | py:65,69, wl:82; sympy out L32 | notes:31,116,157; card:22 | MATCH |
| `c_pole = 1/4` | py:66,69, wl:79; sympy out L32 | notes:120; card:22 | MATCH |
| `c_geom = 3/4` | py:67,69, wl:80; sympy out L32 | notes:122; card:13 (3/4 in module) | MATCH |
| `K_{g,2}=K_{g,4}=0` | py:60-61, wl:73-74 | card:16; notes:27-28,108-110 | MATCH (but hardcoded — see F1) |

INTERNAL scaffolding (no finding): `Y00/Y2A normalization = 1`, symbolic wall coefficients `mu, Tw, TOm, K`, `Omega_Q`, `K_pole`, intermediate `I_mass`/`I_grad`/`I_lap`/`I_pot`, pass/fail flags.

reconciliation: complete; 10 deliverable values checked, 0 misaligned.

## Self-test notes

I checked the variable-independence trap (the `lap_s2`/`lapS2` derivatives act on `Y2A`, which genuinely depend on `th,ph`, so the eigenvalue residual is a real test — not an identically-zero derivative). Parity/vanishing: each `∫ Y00 Y2A dΩ` is `(constant l=0) × (l=2)` over the full sphere, which vanishes by `l`-orthogonality (the `cos(2ph)`/`sin(2ph)`/`cos(ph)`/`sin(ph)` factors integrate to 0 over `ph∈[0,2π]`, and `∫(3cos²−1)sin dth = 0` over `[0,π]`); I confirmed the `grad Y00` integrand is identically 0 (Y00 constant), so A2/A11 are vacuous-but-faithful to the paper's stated reason. Trivial-case: the eigenvalue `l(l+1)=6` is correct for `l=2`, and `Y00 normalization → 1` is the correct sphere integral; the only block where an `assert_zero` is structurally guaranteed independent of the physics is the static-limit block (F1).
