---
unit_id: 130
batch: IV.4
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage130_mouth_bias_map.md"]
  paper_appendix: present
---

# Audit unit 130 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_130.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage130_mouth_bias_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (present; for stage 130 it contains only `\input{stages/stage_130}` at line 1294 — no separate row/narrative, the card is the sole paper-side authority)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.txt`

## What the paper claims

Stage 130 builds an exact "mouth-bias map" for the explicit boundary-layer profile `sigma_Pi(z) = Pi*e^{-Pi z/L}/(L(1-e^{-Pi}))` measured against the first D/N derivative shape `cos(pi z/(2L))`. The card's bottom line (the `Derivation ledger` block-quote, the de-facto Output) is: "Exact `g_Pi`, monotonicity, and unique Family--1 point `Pi_* ≈ 1.5088295`." The notes make this precise with four boxed deliverables: (1) the closed form `g_Pi = 2*Pi*(2*Pi*e^Pi + pi)/((4*Pi^2+pi^2)(e^Pi-1))` (notes §1, lines 33-40); (2) strict monotonicity `dg_Pi/dPi = -(1/L)Cov_Pi(cos(pi z/2L), z) > 0`, with the range `g_Pi: 2/pi → 1` as `Pi: 0^+ → ∞` (notes §2, lines 64-89); (3) the unique compensation point `Pi_* ≈ 1.50882951349316` solving `g_Pi = g_-^{F1}` against the imported lower-branch target `g_-^{F1} ≈ 0.758035078944663` (notes §3, lines 94-110); and (4) the equivalent penetration value `x_* = 1/Pi_* ≈ 0.662765402623160` (notes lines 113-119). The card marks the result `\StatusExactClosure / \StatusNumerical` and lists `g_-^{F1}` as an imported Input (from the stage 122 lower compensated core branch), not derived here.

## What the script claims to verify

Both scripts (a) independently evaluate the integral `∫_0^L sigma*cos(pi z/2L) dz` and assert it equals the paper boxed closed form for `g_Pi`; (b) verify the covariance identity `dg/dPi + (E[fz] - g·E[z])/L = 0`; (c) verify the endpoint limits `g(0+) = 2/pi` and `g(∞) = 1`; (d) certify global strict monotonicity `dg/dPi > 0` on `Pi>0` via a normalized-density check, a symmetrized FKG/Chebyshev double-integral identity `Cov = (1/2)∬(f(z1)-f(z2))(z1-z2)p(z1)p(z2)`, the closed form and sign of `f'(z)` (strictly decreasing `f`), and a `dg/dPi = -Cov/L` consistency tie; (e) bracket `2/pi < g_- < 1` to guarantee a unique root; and (f) solve for `Pi_*`, reporting `Pi_*`, `x_*`, `g(Pi_*)`, `g'(Pi_*)`. The Mathematica script additionally asserts the root with `expectApprox[..., gMinus, 10^-20]` (residual 8.9e-58). The SymPy `g_-` is the decimal literal `0.758035078944663`; the Mathematica `g_-` is the exact stage-122 radical `(2*Sqrt[4107-100*Pi^2]-37*Sqrt[3])/(20*Pi)`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) `g_Pi` boxed closed form (notes §1) | sympy L15-18 `integrate` vs `gPi_boxed`; wl L42-44 `Integrate` vs `gPiBoxed` (`expectZero`) | match |
| (2a) `dg/dPi = -(1/L)Cov` (notes §2 boxed) | sympy L21,31 `cov_id`; wl L47,51 `covId` | match |
| (2b) `dg/dPi > 0` strict monotonicity (notes §2 boxed) | sympy L43-87 (density norm + symmetrized cov identity + f'(z) sign); wl L65-100 (same, sin-positivity via `Reduce`) | match |
| (2c) range `g: 2/pi → 1` (notes §2) | sympy L26,33,36; wl L54-59 limit checks | match |
| (3) unique `Pi_* ≈ 1.50882951349316` (notes §3 boxed, card) | sympy L93-103 bracket + `nsolve`; wl L104-115 bracket + `FindRoot` + `expectApprox` | match |
| (4) `x_* ≈ 0.662765402623160` (notes) | sympy L104; wl L112 | match |
| Input `g_-^{F1} ≈ 0.758035078944663` (card Input, notes §3) | sympy L10 literal; wl L37 exact radical (= stage 122 derivation) | match (imported, anchored upstream) |

`paper_alignment: aligned` — every paper deliverable has a faithful, non-tautological script-side check; nothing is verified that the paper does not claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 17 | `simplify(gPi - gPi_boxed) != 0` raise | claim 1 (g_Pi closed form) | yes |
| A2 | sympy | 31 | `cov_id != 0` raise | claim 2a (covariance form) | yes |
| A3 | sympy | 33 | `simplify(g0 - 2/pi) != 0` raise | claim 2c (lower limit) | yes |
| A4 | sympy | 36 | `simplify(ginf - 1) != 0` raise | claim 2c (upper limit) | yes |
| A5 | sympy | 47 | `simplify(norm_p - 1) != 0` raise | claim 2b (density normalization, support) | yes |
| A6 | sympy | 64 | `simplify(Cov_double - Cov) != 0` raise | claim 2b (symmetrized cov identity) | yes |
| A7 | sympy | 76 | `simplify(gprime - (pi/2L)sin) != 0` raise | claim 2b (f' sign) | yes |
| A8 | sympy | 85 | `simplify(dgPi + Cov/L) != 0` raise | claim 2a/2b (consistency tie) | yes |
| A9 | sympy | 95 | `not (g_lo < g_minus < g_hi)` raise | claim 3 (uniqueness bracket) | yes |
| A10 | sympy | 102 | `nsolve(gPi - g_minus, tol=1e-30)` (raises on non-convergence) | claim 3 (Pi_*) | yes (solver-tol gate; printed `g(Pi_*)` confirms) |
| B1 | mathematica | 44 | `expectZero[gPi - gPiBoxed]` | claim 1 | yes |
| B2 | mathematica | 51 | `expectZero[covId]` | claim 2a | yes |
| B3 | mathematica | 58 | `expectZero[g0 - 2/Pi]` | claim 2c | yes |
| B4 | mathematica | 59 | `expectZero[gInf - 1]` | claim 2c | yes |
| B5 | mathematica | 68 | `expectZero[normP - 1]` | claim 2b | yes |
| B6 | mathematica | 83 | `expectZero[covDouble - cov]` | claim 2b | yes |
| B7 | mathematica | 90 | `expectZero[fPrime + (Pi/2lM)Sin]` | claim 2b | yes |
| B8 | mathematica | 91-95 | `Reduce[Sin>0 && 0<z<lM]` → not False | claim 2b (f' sign, decidable) | yes |
| B9 | mathematica | 99 | `expectZero[dgPi + cov/lM]` | claim 2a/2b | yes |
| B10 | mathematica | 105 | `If[gLo < gMinus < 1]` pass/fail | claim 3 (bracket) | yes |
| B11 | mathematica | 115 | `expectApprox[gPi@piStar, gMinus, 1e-20]` (res 8.9e-58) | claim 3 (Pi_* root) | yes |

No tautological rows. The only redundancy is A8/B9 (`dg/dPi = -Cov/L`) overlapping with the covariance identity A2/B2; but the redundancy is intentional ("non-tautological: a wrong gPi or wrong Cov breaks it") and serves to tie the certified sign back to the verified derivative — it is not a defect.

## Findings

None.

### Adversarial attacks attempted (and why they failed)

- **Tautology on the g_Pi check?** No — `gPi` is produced by an independent `integrate`/`Integrate` of `sigma*f`, and `gPi_boxed`/`gPiBoxed` is a separately typed closed form; the assertion compares two independently-produced expressions. If the integral evaluated to a different rational function the check would fail.
- **`nsolve` with no assertion in SymPy?** SymPy reports `Pi_*` via `nsolve(... , tol=1e-30)`, which raises if it cannot converge to that tolerance, and prints `g(Pi_*) = 0.758035078944663043...` matching `g_minus`. The deliverable is genuinely exercised. Mathematica additionally has an explicit `expectApprox` (residual 8.9e-58). Not `insufficient_verification`.
- **Symbol-domain error from the `Pi` overload in Mathematica?** Checked: `.wl` uses built-in `Pi` (=π) inside `Cos[Pi z/(2 lM)]` and `piM` for the bias parameter Π. No collision — output line 7 shows `(2*piM*(Pi + 2*E^piM*piM))/...`, correctly distinguishing the two. SymPy uses `Pi` for the bias and `sp.pi` for π. Both consistent.
- **Monotonicity only proved at a sample?** No — the comment explicitly disclaims a finite sweep; the proof is a symbolic certificate (normalized density + symmetrized covariance identity + closed-form sign of `f'` decidable on the bounded domain), valid for all `Pi>0`. The Mathematica `Reduce` step independently decides `Sin>0` on `(0,lM)`.
- **Missing-branch on monotonicity range?** Both endpoint limits `2/pi` and `1` are asserted, matching notes §2.
- **`g_minus` hardcoded with no provenance?** It is an Input per the card, derived at stage 122: notes §3 line 98 and stage 122 notes line 63 give `0.758035078944663`; the exact form `(2√(4107-100π²)-37√3)/(20π)` is derived at stage 122 (notes line 56, stage 122 sympy line 60, stage 122 wl line 82). Properly anchored upstream; not an orphaned literal.

## Independent-derivation check (Mathematica)

The two scripts follow the same paper-prescribed derivation outline (the notes themselves lay out: §1 g_Pi, §2 covariance/monotonicity, §3 compensation point), so structural parallelism is expected and appropriate. Within that outline each engine derives independently: every integral (`Integrate` vs `integrate`), derivative (`D` vs `diff`), and limit is recomputed by the respective CAS, and the two engines use genuinely different methods at the load-bearing points:
- `g_-^{F1}` representation: Mathematica uses the exact radical `(2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi)` (wl line 37); SymPy uses the decimal literal `0.758035078944663` (py line 10). These are not echoes of each other — Mathematica re-derives the algebraic value, SymPy carries the decimal.
- `f'` positivity: Mathematica decides it with `Reduce[Sin[...] > 0 && 0 < z < lM, z, Reals]` (wl line 91, output line 23: `lM > 0 && 0 < z < lM`); SymPy instead certifies via closed-form equality `gprime - (pi/2L)sin(...) == 0` plus a prose argument (py lines 73-80). Genuinely different machinery for the same lemma.
This is an independent re-derivation, not a `mathematica_transliteration`.

## Engine cross-check

Both engines produce matching results at the precision they claim:

| quantity | SymPy | Mathematica | agree? |
|---|---|---|---|
| `g_Pi` closed form | `2*Pi*(2*Pi*exp(Pi)+pi)/((4*Pi**2+pi**2)*(exp(Pi)-1))` | `(2*piM*(Pi+2*E^piM*piM))/((-1+E^piM)*(Pi^2+4*piM^2))` | yes (same function, symbol roles swapped) |
| `cov_id` residual | 0 | 0 | yes |
| limits | `2/pi`, `1` | `2/Pi`, `1` | yes |
| `Cov_double` | rational fn (out L5) | rational fn (out L25) | yes (same, π/Π roles swapped) |
| `Pi_*` | 1.5088295134931586… | 1.5088295134931556… | yes — agree to ~14 sig figs; last-digit divergence is from SymPy's rounded `g_-` literal vs Mathematica's exact radical, well inside the 15-digit value the notes/card quote |
| `x_*` | 0.6627654026231601… | 0.6627654026231614… | yes (same, ~14 sig figs) |
| `g'(Pi_*)` | 0.07144535580831948… | 0.07144535580831952… | yes (~14 sig figs) |

No `engine_disagreement`: the only differences are at the precision floor set by SymPy's decimal `g_-` input, and both round to the values stated in the notes and card.

## Verdict justification

`clean`. I read the card, the notes, and the appendix `\input` line first, built the four-deliverable model (g_Pi closed form, covariance-based strict monotonicity with range 2/pi→1, unique Pi_*≈1.50882951349316, and x_*≈0.662765402623160), then attacked the scripts against it. Every deliverable maps to a non-tautological, well-anchored assertion in both engines; the monotonicity claim is a genuine global symbolic certificate (not a sweep); the imported `g_-^{F1}` is anchored to its stage-122 derivation in both decimal and exact-radical form; the two engines agree to the precision the notes/card state, with the only divergence explained by SymPy's rounded input literal. Attacks on tautology, missing branch, symbol-overload (the `Pi`/`piM` distinction), print-only verification of `Pi_*`, and transliteration all failed. The Mathematica is an independent re-derivation using different methods at the load-bearing points (`Reduce` for sign, exact radical for `g_-`).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 7 deliverable values checked, 0 misaligned

Outputs are fresh: sympy `.txt` mtime 2026-05-29 15:39:39 > sympy `.py` 15:32:01; mathematica `.txt` 15:39:39 > `.wl` 15:31:44. Reconciliation is based on script source plus the committed saved outputs.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `g_Pi = 2*Pi*(2*Pi*e^Pi+pi)/((4*Pi^2+pi^2)(e^Pi-1))` | sympy py:16, out L1; wl:43, out L7 | notes L33-40 (boxed §1); card L16 ("Exact `g_Pi`") | MATCH |
| `g_Pi(0^+) = 2/pi` | sympy py:26/33, out L3; wl:54/58, out L11 | notes L87 (range `2/π → 1`) | MATCH |
| `g_Pi(∞) = 1` | sympy py:29/36, out L4; wl:55/59, out L12 | notes L87 | MATCH |
| `dg_Pi/dPi > 0` (strict monotonicity) | sympy py:80-87, out L7; wl:99-100, out L28 | notes L65-82 (boxed §2); card L16 ("monotonicity") | MATCH |
| `Pi_* ≈ 1.50882951349316` | sympy py:103, out L9 (1.50882951349315861…); wl:111, out L31 (1.5088295134931556…) | notes L108 (boxed §3), L144; card L16 (`1.5088295`) | MATCH |
| `x_* ≈ 0.662765402623160` | sympy py:104, out L10 (0.662765402623160072…); wl:112, out L32 (0.6627654026231614…) | notes L115 (boxed §3) | MATCH |
| `g_-^{F1} ≈ 0.758035078944663` (imported Input; exact `(2√(4107-100π²)-37√3)/(20π)`) | sympy py:10, out L8; wl:37, out L30 | notes L98 (§3); card L9 (Input); stage 122 notes L56/L63 | MATCH (anchored upstream at stage 122) |

INTERNAL scaffolding (accounted for, no finding): covariance identity residual (0), `sigma` normalization residual (0), symmetrized covariance identity residual (0), `f'(z)` closed-form residual (0), `Cov_double` intermediate rational form (drives the sign certificate), `dg/dPi = -Cov/L` consistency residual (0), bracket inequality `2/pi < g_- < 1`, `g(Pi_*)` root residual (sympy printed 0.7580350789446630…; wl `expectApprox` 8.9e-58), `g'(Pi_*)` value (≈0.0714453558…; reported as confirmation the root is non-degenerate, not a stated deliverable).

All 7 stage deliverable values reconcile against the card and/or notes. No MISMATCH and no MISSING-DELIVERABLE.

## Self-test notes

I checked: (1) variable independence — `diff(gPi, Pi)`/`D[gPi, piM]` is taken w.r.t. the bias parameter that `gPi` genuinely depends on (the integral result is a non-trivial function of Pi/piM), so the derivative is not identically zero; (2) symmetry/parity — the symmetrized covariance integrand `(1/2)(f(z1)-f(z2))(z1-z2)p1p2` over `[0,L]^2` is verified equal to `Cov` by the engines and its sign is established by `f` strictly decreasing, consistent with `Cov<0`; (3) trivial-case — the endpoint limits `2/pi` and `1` are the correct uniform-source and point-source values and are asserted, and the root `g(Pi_*)=g_-` reduces to the imported target to 1e-58 in Mathematica. No directive is written (zero findings).
