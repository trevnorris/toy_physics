---
unit_id: 100
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage100_outgoing_normalization_factorization.md]
  paper_appendix: present
---

# Audit unit 100 red-team report

## ⚠️ Orchestrator override (2026-06-05) — verdict clean → findings (1)

The audit agent initially returned `verdict: clean`, judging the `.wl` an independent
re-derivation on the strength of `sp.im(...)` vs `.../I`. The orchestrator's ground-truth
read of both files (corroborating the first-pass IV.1 heads-up that flagged 100 as the
standout dual-engine gap) OVERTURNS that call: the `.wl` is a line-by-line transliteration
of the `.py` — identical `sigma_can`, identical rational `Y`, the same `Series[]`/`series()`
black box on the same expression, and the same `Coefficient`/`.coeff` extraction; `sp.im`
vs `/I` is the same operation in different syntax, not an independent route. An independent
Mathematica route IS feasible (e.g. analytic geometric-series expansion of the pole
denominator, collecting orders by hand). Per [[feedback-dual-engine-required]] this is a
user-level call; the user AUTHORIZED the independent-route rewrite on 2026-06-05. The
finding and the requirement/acceptance criteria are recorded in
`redteam/pass2/directives/stage_100.md` (F1 — `mathematica_transliteration`). Codex designs
and writes the new `.wl` per [[feedback-claude-reviews-codex-codes]]. Everything else in this
report (paper alignment, value reconciliation = 12 values, 0 misaligned, the non-tautological
closure assertion) stands.

### F1 — mathematica_transliteration

**Severity:** medium
**Files:** `mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl` (derivation body, ~lines 40–72)

**What's wrong:** the `.wl` echoes the `.py`'s series-expansion choreography step-for-step
rather than deriving the low-frequency coefficients `K2`/`K4`/`Gamma5` by an independent
method. See the directive for the corresponding-section quotes.

**Required change / Verification:** see `redteam/pass2/directives/stage_100.md` (independent
re-derivation, same deliverables/values, `chiQ` free, exits 0).

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_100.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage100_outgoing_normalization_factorization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows referencing stage 100: lines 26, 86, 1176, 1234)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.txt`

## What the paper claims

Stage 100 is a retarded 2.5PN factorization-ledger step. Its bottom-line claim (card line 16, boxed `\stagefield{Derivation ledger}`): "Full odd normalization factorizes as \(\widehat m_0^{\,2}\chi_Q N_Q=1\)." The card states the computation "isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\)." The notes supply the supporting derivation: write the actual retarded branch as `Yhat_Q^ret = 3/4 + (1/4)/(1 - omega^2/Omega^2 - i chi_Q sigma_can omega^5) + O(omega^6)` with `sigma_can = 9/(8 Omega^5)`; extract the low-frequency tuple `Kbar_2 = Kbar_0/(4 Omega^2)`, `Kbar_4 = Kbar_0/(4 Omega^4)`, `Gammabar_5 = chi_Q*9 Kbar_0/(32 Omega^5)`; with targets `K0_t = 64 G Omega^5/(45 c^5)` and `Gamma5_t = 2 G/(5 c^5)` and `N_Q := K0/K0_t`, derive the even ratios `K2/K2_t = K4/K4_t = N_Q` and the odd ratio `Gamma5/Gamma5_t = chi_Q N_Q`; finally impose the observable condition `mhat_0^2 Gamma5 = Gamma5_t` to force `mhat_0^2 chi_Q N_Q = 1`. The card's three `\stagefield{Checks}` items list (i) factor-separation of the product, (ii) higher odd terms beginning beyond the 2.5PN coefficient, and (iii) the DtN l=2 fingerprint vs the z-expansion; the scripts explicitly delegate (ii) to stage 102 and (iii)/the `chi_Q=1` pin to stage 097 (carry-forward annotations), and own (i)/the closure here.

## What the script claims to verify

Both scripts build `Yhat_Q^ret` from the one-pole renormalized form, series-expand to O(omega^5), and read off K2, K4, Gamma5 from the coefficients (Gamma5 = imaginary part of the omega^5 coefficient times K0). They then construct the GR targets and N_Q = K0/K0_t and assert the three ratios `K2/K2_t - N_Q == 0`, `K4/K4_t - N_Q == 0`, `Gamma5/Gamma5_t - chi_Q*N_Q == 0`. They then impose the observable closure `mhat_0^2 Gamma5 = Gamma5_t`, form `closure_ratio = (mhat_0^2 Gamma5 - Gamma5_t)/Gamma5_t`, and assert `closure_ratio - (mhat_0^2 chi_Q N_Q - 1) == 0`, i.e. that the closure residual factorizes as the headline `mhat_0^2 chi_Q N_Q - 1`. `chi_Q` is carried as a free real symbol (not pinned to 1) by design.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Y series form `3/4 + (1/4)/(1 - w^2/O^2 - i chi sigma w^5)` (notes 17-18) | sympy L33-34 / wl L41-42 build and expand this exact form | match |
| `sigma_can = 9/(8 Omega^5)` (notes 22) | sympy L32 / wl L40 | match |
| `Kbar_2 = K0/(4 Omega^2)` (notes 26) | sympy L36 / wl L44 → output `K0/(4*Omega**2)` | match |
| `Kbar_4 = K0/(4 Omega^4)` (notes 28) | sympy L37 / wl L45 → output `K0/(4*Omega**4)` | match |
| `Gammabar_5 = chi_Q*9 K0/(32 Omega^5)` (notes 30/73) | sympy L38 / wl L46 → output `9*K0*chiQ/(32*Omega**5)` | match |
| `K2/K2_t = N_Q`, `K4/K4_t = N_Q` (notes 57/59) | sympy L54-55 / wl L60-61 asserts ==0 | match |
| `Gamma5/Gamma5_t = chi_Q N_Q` (notes 77) | sympy L56 / wl L62 asserts ==0 | match |
| Headline `mhat_0^2 chi_Q N_Q = 1` (card 16, notes 91) | sympy L65-70 / wl L68-72 closure check | match |
| Check (ii) higher odd terms beyond 2.5PN (card 23) | delegated to stage 102 (docstring sympy L7-10 / wl L29) | delegated (documented) |
| Check (iii) DtN l=2 fingerprint + chi_Q=1 pin (card 24) | delegated to stage 097 (docstring sympy L11-14 / wl L30) | delegated (documented) |

`paper_alignment: aligned`. Checks (ii)/(iii) are not silent gaps: both scripts carry explicit, named carry-forward annotations routing them to stages 102 and 097, and the card's box carries `chi_Q` as a free factor — consistent with this stage owning only the factorization closure.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 54 | `simplify(K2/K2_t - NQ) == 0` | even ratio = N_Q | yes |
| A2 | sympy | 55 | `simplify(K4/K4_t - NQ) == 0` | even ratio = N_Q | yes |
| A3 | sympy | 56 | `simplify(Gamma5/Gamma5_t - chiQ*NQ) == 0` | odd ratio = chi_Q N_Q | yes |
| A4 | sympy | 70 | `closure_check == 0` (closure_ratio == mhat0^2 chiQ NQ - 1) | headline factorization | yes |
| A5 | mathematica | 60 | `expectZero[K2/K2_t - NQ]` | even ratio = N_Q | yes |
| A6 | mathematica | 61 | `expectZero[K4/K4_t - NQ]` | even ratio = N_Q | yes |
| A7 | mathematica | 62 | `expectZero[Gamma5/Gamma5_t - chiQ*NQ]` | odd ratio = chi_Q N_Q | yes |
| A8 | mathematica | 71-72 | `expectZero[closureRatio - (mHat0^2 chiQ nQDerived - 1)]` | headline factorization | yes |

All eight assertions trace to a specific paper deliverable. None are tautological by construction: each LHS (K2, K4, Gamma5) is read from the *series expansion* of the independently-built `Y`, not from a hardcoded form, so each `==0` can fail if the series factor is wrong.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration, though the derivation chain is necessarily short and parallel. Justifying excerpts:

1. **Imaginary-part extraction differs (the decisive marker).** SymPy uses `sp.im(Yser.coeff(omega, 5))` (sympy L38), which discards any real part of the coefficient. Mathematica instead divides by the unit imaginary: `(Coefficient[ySeries, omega, 5]/I)*k0` (wl L46). A line-by-line port would have written `Im[...]` in WL; the author chose a different (and genuinely distinct) algebraic operation, so this is not echoed algebra.
2. **Symbol choreography differs idiomatically.** SymPy: `Omega, K0, mhat0` (L22-24); Mathematica: `omegaQ, k0, mHat0` (L33-38). The WL declares assumptions via `$Assumptions = Element[...] && ... > 0` (L35-38) versus SymPy's per-symbol `positive=True` keyword (L22-29) — different mechanisms for the same domain.
3. **Series + assertion idioms differ.** SymPy: `sp.series(Y, omega, 0, 6).removeO()` then `assert sp.simplify(...) == 0` (L34, L54-56); Mathematica: `Normal[Series[yRet, {omega, 0, 5}]]` then a custom `expectZero[name, expr]` helper that runs `FullSimplify[Together[Expand[expr]]]` and `Exit[1]` on nonzero (L20-24, L42, L60-62). Both genuinely re-expand the same physical `Y` form and re-read its coefficients.

These are the same identities derived by two engines from the same physical premise (the renormalized one-pole branch) — which is exactly what the dual-engine policy wants — rather than one engine transcribing the other's intermediate algebra. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce identical results (modulo symbol names):

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| Y series | `1 + omega**2/(4*Omega**2) + omega**4/(4*Omega**4) + 9*I*chiQ*omega**5/(32*Omega**5)` | `1 + omega^2/(4*omegaQ^2) + omega^4/(4*omegaQ^4) + (9*I/32)*chiQ*omega^5/omegaQ^5` |
| K2 | `K0/(4*Omega**2)` | `k0/(4*omegaQ^2)` |
| K4 | `K0/(4*Omega**4)` | `k0/(4*omegaQ^4)` |
| Gamma5 | `9*K0*chiQ/(32*Omega**5)` | `(9*chiQ*k0)/(32*omegaQ^5)` |
| NQ | `45*K0*c**5/(64*G*Omega**5)` | `(45*cLight^5*k0)/(64*gNewton*omegaQ^5)` |
| 3 ratio checks | all `0` | all `PASS` (=0) |
| closure_ratio | `-1 + 45*K0*c**5*chiQ*mhat0**2/(64*G*Omega**5)` | `-1 + (45*chiQ*cLight^5*k0*mHat0^2)/(64*gNewton*omegaQ^5)` |
| closure_check | `0` | `PASS` (=0) |

`engines_agree: true`.

## Verdict justification

`clean`. I read the card, the single notes file, and the part-04 appendix rows (26, 86, 1176, 1234) before opening the scripts; the scripts' verified claim — the factorization `mhat_0^2 chi_Q N_Q = 1` plus the supporting even/odd ratio tuple and the `Y` series — matches the card's boxed `\stagefield{Derivation ledger}` claim and every notes equation. Attacks tried that failed: (1) **tautology** on the closure check — it is not, because Gamma5 is read from the independently-built series, so a wrong series factor would break both A3/A7 and A4/A8 (the comment at sympy L60-64 / wl L64-67 states this correctly); (2) **hidden imaginary contamination** of K2/K4 — the omega^2 and omega^4 coefficients are purely real (`u`'s first term is real, `u^2`'s omega^4 term is real), so `coeff(omega,2/4)` carries no spurious `i`; (3) **wrong sigma_can or omega^5 factor** — `(1/4)*i*chiQ*(9/8)/Omega^5 = 9 i chiQ/(32 Omega^5)` matches notes line 30 exactly; (4) **symbol-domain error** — `chi_Q` is correctly left unconstrained-positive (declared `real` only, with an explicit comment sympy L25-28 / wl L34 that pinning it positive would falsely constrain the factorization), while G, c, Omega, K0, mhat0 are physical positives; (5) **missing-engine / delegated-check gap** — Checks (ii)/(iii) are explicitly and correctly delegated to stages 102/097 with named annotations, and the card box carries chi_Q free. Both engines present, agree, and outputs are fresh.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Y series `1 + w^2/(4O^2) + w^4/(4O^4) + 9i chiQ w^5/(32 O^5)` | py L34 / out L1; wl L42 / out L5 | notes L17-18 (Y form) | MATCH |
| `sigma_can = 9/(8 Omega^5)` | py L32; wl L40 | notes L22 | MATCH |
| `K2 = K0/(4 Omega^2)` | py L36 / out L2; wl L44 / out L6 | notes L26 (Kbar_2) | MATCH |
| `K4 = K0/(4 Omega^4)` | py L37 / out L3; wl L45 / out L7 | notes L28 (Kbar_4) | MATCH |
| `Gamma5 = 9 chiQ K0/(32 Omega^5)` | py L38 / out L4; wl L46 / out L8 | notes L30/L73 (Gammabar_5) | MATCH |
| `N_Q = 45 K0 c^5/(64 G Omega^5)` (= K0/K0_t) | py L44 / out L5; wl L52 / out L9 | notes L53 (N_Q := K0/K0_t) | MATCH |
| `K0_t = 64 G Omega^5/(45 c^5)` | py L40; wl L48 | notes L45 | MATCH |
| `Gamma5_t = 2 G/(5 c^5)` | py L43; wl L51 | notes L68 | MATCH |
| `K2/K2_t = N_Q` | py L54 / out L6; wl L60 / out L10-11 | notes L57 | MATCH |
| `K4/K4_t = N_Q` | py L55 / out L7; wl L61 / out L12-13 | notes L59 | MATCH |
| `Gamma5/Gamma5_t = chi_Q N_Q` | py L56 / out L8; wl L62 / out L14-15 | notes L77 | MATCH |
| headline `mhat_0^2 chi_Q N_Q = 1` | py L65-74 / out L9-12; wl L68-75 / out L16-20 | card L16 (box), notes L91 | MATCH |

INTERNAL scaffolding (no finding): `closure_residual`, `closure_ratio` intermediate, `closure_check` residual (=0), the PASS/FAIL flags, and the "STAGE 100 AUDIT PASSED" / "Stage 100 Mathematica audit passed." banners.

reconciliation: complete; 12 values checked, 0 misaligned

## Self-test notes

Checked: (1) variable-independence — no `diff`/`D` in either script (pure series + algebra), so the zero-derivative trap does not apply. (2) Parity/symmetry — n/a (finite series, no unbounded integrals). (3) Trivial-case pre-check — substituting chi_Q=0 gives Gamma5=0, closure_ratio=-1 = (0-1), so A4/A8 reduce to 0 correctly; with chi_Q=1, mhat0=1, K0=K0_t (N_Q=1), closure_ratio=0 as required by the headline. (4) Paper round-trip — every prescribed value already reconciles to notes/card with no new misalignment introduced (no fix prescribed; verdict clean). The only borderline call was transliteration vs independent derivation; the divergent `sp.im` vs `/I` imaginary-part extraction tipped it to independent.
