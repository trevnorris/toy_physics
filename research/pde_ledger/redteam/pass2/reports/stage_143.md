---
unit_id: 143
batch: IV.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-07T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage143_equal_normalized_singular_limit.md]
  paper_appendix: present
---

# Audit unit 143 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_143.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage143_equal_normalized_singular_limit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_143}` row at line 1320; no separate prose row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.txt`

## What the paper claims

The stage proves that the simple equal-normalized mouth-source branch `g_c = 1` is a *singular* limit, not a regular finite-bias branch. The card's load-bearing quote (`stage_143.tex:16`) is verbatim: "\(\mathfrak g_\Pi<1\) for all finite \(\Pi>0\); equal-normalized branch requires \(\Pi\to\infty\) and divergent traction." The notes (`...stage143...md`) carry three boxed deliverables: (1) `0 < g_Π < 1` for every finite `Π>0` (notes §1, lines 48-51), established by an exact three-term positive decomposition of the numerator of `1-g_Π` whose keystone is `e^Π - 1 - Π - Π²/2 > 0`; (2) `lim_{Π→∞} g_Π = 1`, so the equal-normalized branch is the `Π→+∞` delta-function/point-source limit (notes §2, lines 62-68); and (3) on that branch the normalized traction `T̂_m(Π) ~ sqrt(9/(20(1-R_∞))) Π^{1/2} ≈ 0.725669130700713 Π^{1/2}` diverges, where `R_∞ = (1-r_{F1})²/(1+r_{F1}²) ≈ 0.145454452260420` (notes §3, lines 86-119). The `g_Π`, `S_q`, `R_q`, `Σ₀`, `T̂_m`, and `r_{F1}` definitions are carried forward from Stage 142.

## What the script claims to verify

Both engines (a) re-derive `num := (1-g_Π)·(4Π²+π²)·(e^Π-1)` from the explicit `g_Π` family and assert it equals the claimed three-term decomposition `π²(e^Π-1-Π-Π²/2) + Π(π²-2π) + Π²(π²/2-4)`; (b) prove each term is positive for `Π>0` — the two constant coefficients by direct positivity, and the exponential remainder `R(Π)=e^Π-1-Π-Π²/2` by a Taylor-remainder monotonicity chain `R(0)=R'(0)=R''(0)=0, R'''(Π)=e^Π>0`; (c) compute the endpoint limits `lim_{Π→0+} g_Π = 2/π` and `lim_{Π→∞} g_Π = 1`; and (d) take the `Π→∞` limits of `R_q`, `S_q`, `Σ₀/Π`, `T̂_m/√Π`, asserting `R_∞=(1-r)²/(1+r²)`, `S_∞=1`, and `T̂_m/√Π → sqrt((9/20)/(1-R_∞))`, then print the numeric values. The Mathematica script additionally proves `num > 0` and `R(Π) > 0` directly via `Reduce[..., piM, Reals]`.

## Paper ↔ script cross-check

| paper/notes deliverable | script-side check | status |
|---|---|---|
| `0 < g_Π < 1` for all finite `Π>0` (numerator positivity via 3-term split) | num−decomp==0 (sympy 44 / wl 58); each term positive (sympy 50-61 / wl 71-88); wl also `Reduce[num>0]` (wl 61-65) | match |
| `lim_{Π→∞} g_Π = 1` (singular point-source limit) | `expect_equal lim==1` (sympy 69 / wl 97) | match |
| `T̂_m ~ sqrt(9/(20(1-R_∞)))·√Π`, diverges | `that_ratio == sqrt((9/20)/(1-Rinf))` (sympy 86 / wl 116); `Σ₀/Π → finite` printed | match |
| `R_∞ = (1-r_{F1})²/(1+r_{F1}²) ≈ 0.145454452260420` | `Rinf == (1-r)²/(1+r²)` (sympy 84 / wl 114); N[…,30] printed | match |
| `S_q → 1` (susceptibility closure) | `S_infty == 1` (sympy 85 / wl 115) | match |
| (none) — `lim_{Π→0+} g_Π = 2/π` | `expect_equal g0 == 2/pi` (sympy 68 / wl 96) | extra (supporting, not load-bearing; see Verdict) |

`paper_alignment: aligned`. Every boxed deliverable in the notes and the card's load-bearing quote has a faithful, non-tautological script-side check in both engines. The one script-side `extra` (the `Π→0+` limit) is a supporting computation, not a paper deliverable, and does not contradict the notes' weaker `0 < g_Π` bound.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 44 | `expect_zero(num - decomp)` | deliverable 1 (positivity split) | yes |
| A2 | sympy | 50 | `expect_positive(pi**2-2*pi)` | deliverable 1 (term 2) | yes |
| A3 | sympy | 51 | `expect_positive(pi**2/2-4)` | deliverable 1 (term 3) | yes |
| A4 | sympy | 57-59 | `expect_equal R(0),R'(0),R''(0)==0` | deliverable 1 (Taylor base) | yes |
| A5 | sympy | 60 | `expect_zero(R'''(Pi)-exp(Pi))` | deliverable 1 (Taylor keystone) | yes |
| A6 | sympy | 61 | `expect_positive(exp(Pi))` | deliverable 1 (R'''>0 ⟹ R>0) | yes |
| A7 | sympy | 68 | `expect_equal g0 == 2/pi` | none (supporting extra) | n/a |
| A8 | sympy | 69 | `expect_equal ginf == 1` | deliverable 2 | yes |
| A9 | sympy | 84 | `expect_equal Rinf == (1-r)**2/(1+r**2)` | deliverable 4 | yes |
| A10 | sympy | 85 | `expect_equal Sinf == 1` | deliverable 5 | yes |
| A11 | sympy | 86 | `expect_equal that_ratio == sqrt((9/20)/(1-Rinf))` | deliverable 3 | yes |
| B1 | mathematica | 58 | `expectZero(num - decomp)` | deliverable 1 | yes |
| B2 | mathematica | 61-65 | `Reduce[num>0,piM,Reals]==True` | deliverable 1 (independent global proof) | yes |
| B3 | mathematica | 75-80 | `Reduce[R>0,piM,Reals]==True` | deliverable 1 (independent global proof) | yes |
| B4 | mathematica | 71-72 | `expectPositive(Pi^2-2*Pi),(Pi^2/2-4)` | deliverable 1 (terms 2,3) | yes |
| B5 | mathematica | 84-88 | Taylor chain R(0..2)=0, R'''=Exp, Exp>0 | deliverable 1 (backing route) | yes |
| B6 | mathematica | 97 | `expectEqual gInf == 1` | deliverable 2 | yes |
| B7 | mathematica | 114-116 | `rInf/sInf/tHatRatio` asserts | deliverables 3,4,5 | yes |

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an **INDEPENDENT** derivation of the load-bearing claim, not a transliteration. Three lines of evidence:

1. **Variable-naming inversion.** In sympy, `Pi` is the bias variable and `pi = sp.pi` (sympy 31-32). In Mathematica, `Pi` is the built-in π and `piM` is the bias variable (`$Assumptions = ... piM > 0`, wl 48; `r = Sqrt[4107 - 100*Pi^2]/(10*Pi)`, wl 50). A mechanical port would keep the role of each symbol aligned; here the two names are swapped between engines, indicating the `.wl` was written from the physics, not echoed from the `.py`.

2. **Distinct proof method for the central positivity claim.** The `.py` proves `num > 0` *only* through the Taylor-remainder monotonicity chain (sympy 52-61). The `.wl` proves it through a genuinely different decision procedure — `numPositiveCheck = Reduce[num > 0, piM, Reals]` returning `True` (wl 61-65, output line 12-13) — and *separately* via `Reduce[exp remainder > 0, piM, Reals]` (wl 75-80, output line 22-23). `Reduce` over the reals is cylindrical-algebraic-style quantifier elimination, an independent route to global positivity that the sympy side does not use. The Taylor chain is retained in the `.wl` only as a "backing route" (wl 81-88).

3. **Shared scaffolding is unavoidable common premise.** The identical `g_Π`, `S_q`, `R_q`, `Σ₀`, `T̂_m`, and `r_{F1}` definitions are the *same physical system* both engines must model (carried forward from Stage 142); reusing them is required, not transliteration. The limit checks (`lim g_Π`, `R_∞`, `S_∞`) are computed by each engine's own `limit`/`Limit` independently.

The positivity proof is genuinely **global**, not point-sampled: sympy establishes `R(Π)>0` ∀Π>0 by the exact derivative chain (every term symbolic, no numeric substitution beyond the `Π=0` boundary values), and Mathematica confirms `num>0` and `R>0` over all reals via `Reduce`. No stale `168π²` appears anywhere — the script uses `4107 - 100π²` (sympy 36 / wl 50), which is the correct Stage-142 `r_{F1}` constant (confirmed against `stage142...sympy_audit.py:34`, identical form).

## Engine cross-check

Both engines pass every assertion and print matching results:

| quantity | sympy output | mathematica output | agree |
|---|---|---|---|
| num − decomp | `0` (line 10) | `0` + PASS (lines 10-11) | yes |
| lim_{Π→0+} g_Π | `2/pi` (line 26) | `2/Pi` (line 38, `Pi`=π) | yes |
| lim_{Π→∞} g_Π | `1` (line 27) | `1` (line 39) | yes |
| R_∞ (closed form) | `(-sqrt(4107-100π²)+10π)²/4107` (line 30) | `1 - 20π√(4107-100π²)/4107` (line 44) | yes (algebraically equal; both ⇒ residual 0 vs `(1-r)²/(1+r²)`) |
| R_∞ (numeric) | `0.145454452260420126...` (line 50) | `0.1454544522604201261014...` (line 67) | yes |
| lim T̂_m/√Π (numeric) | `0.725669130700713219...` (line 51) | `0.7256691307007132197...` (line 68) | yes |

The closed forms for `R_∞` are printed in superficially different but algebraically identical forms; both engines independently verify `R_∞ == (1-r)²/(1+r²)` with residual 0, so the apparent difference is cosmetic, not an `engine_disagreement`.

## Verdict justification

`clean`. I attacked: (1) the decomposition assertion — it is non-tautological because `num` is computed from the `g_Π` *definition* while `decomp` is an independently written claimed form; a wrong split would fail. (2) The positivity argument — I checked the Taylor chain closes globally (R'''=e^Π>0 ∀Π, with R(0)=R'(0)=R''(0)=0 ⟹ R>0 ∀Π>0) and confirmed Mathematica's `Reduce` returns `True` (a genuinely independent global proof, not point-sampling). (3) The carried-forward `r_{F1}` constant — `4107-100π²` matches Stage 142 verbatim; no stale `168π²`. (4) The numeric reconciliation — `R_∞` and `T̂_m/√Π` match the notes' boxed values to all printed digits. (5) Transliteration — the `.wl` uses inverted symbol roles and an independent `Reduce` proof method, so it is not a port. The one script-side extra (`lim_{Π→0+} g_Π = 2/π`) is a correct supporting computation absent from the notes/card; it is not a deliverable and does not contradict the notes' weaker `0 < g_Π` bound, so it is not a finding. The card's `Checks` items 1 (gain pair `(M_s,M_q)`) are generic ledger-block checklist boilerplate, not paper Outputs; the load-bearing claim (the boxed `g_Π<1`/singular-limit quote) is fully and independently exercised by both engines. I read the card, the notes, and the appendix row, and the script's verified claim matches the paper's claim.

## Self-test notes

Variable-independence trap: the only derivatives are of `R(Π)=e^Π-1-Π-Π²/2` w.r.t. the bias variable (`Pi` in sympy, `piM` in wl), which `R` genuinely depends on, so none are identically-zero artifacts; `R'''=e^Π` is verified `≠0`. Trivial-case pre-check: at `Π→0+`, `num` and `R` both → 0 with the boxed limits matching; for `Π>0`, `R'''=e^Π>0` makes `R` strictly increasing from 0, confirming the `assert_positive`s give nonzero/positive literals. No unbounded-domain integrals appear, so parity is N/A. Paper round-trip: the audit prescribes no fix, so no new `paper_misalignment` can be introduced.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `0 < g_Π < 1` ∀ finite Π>0 (3-term positive split of num of `1-g_Π`) | py 43-44, wl 57-58; sympy out line 10, wl out lines 10-13 | notes:48-51 (boxed `0<g_Π<1`); card quote stage_143.tex:16 | MATCH |
| `lim_{Π→∞} g_Π = 1` (singular point-source limit) | py 69, wl 97; sympy out line 27/29, wl out line 39/43 | notes:62-68 (boxed `Π→+∞`); card stage_143.tex:16 | MATCH |
| `R_∞ = (1-r_{F1})²/(1+r_{F1}²)` (closed form) | py 84, wl 114; sympy out line 30/34, wl out line 44/49 | notes:92-94 | MATCH |
| `R_∞ ≈ 0.145454452260420` (numeric) | py 98, wl 128; sympy out line 50, wl out line 67 | notes:95 (`≈ 0.145454452260420`) | MATCH |
| `T̂_m/√Π → sqrt((9/20)/(1-R_∞)) ≈ 0.725669130700713` | py 86/99, wl 116/129; sympy out line 33/36/51, wl out line 47/53/68 | notes:108-111 (boxed, `≈ 0.725669130700713`) | MATCH |
| `S_q → 1` (susceptibility closure) | py 85, wl 115; sympy out line 31/35, wl out line 45/51 | notes:90 (`S_q(Π)→1`) | MATCH |

INTERNAL (accounted, no finding): `r_{F1} = sqrt(4107-100π²)/(10π)` (carried-forward input from Stage 142, not a Stage-143 deliverable; documented at stage142 sympy:30-34 / notes:38); `num` raw expanded form, `decomp`, `Σ₀/Π` limit (`4107/(20π√(4107-100π²))`, an intermediate feeding `T̂_m`); `lim_{Π→0+} g_Π = 2/π` (supporting endpoint limit, not a notes/card deliverable; absent from notes but a correct value that does not contradict the weaker `0<g_Π` box, classified INTERNAL/supporting per the false-positive guard); Taylor-base values `R(0)=R'(0)=R''(0)=0`, `R'''=e^Π` (proof scaffolding); the `2/π < g_Π < 1` ledger printout (a non-asserted print annotation, stronger than the proven/boxed `0<g_Π`, but consistent with it since `2/π>0`).

Every stated deliverable value reconciles between the scripts, the saved outputs, the `.tex` card quote, and the `.md` notes.
