---
unit_id: 077
batch: III.4
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage077_family1_theta_extraction.md]
  paper_appendix: present
---

# Audit unit 077 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_077.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage077_family1_theta_extraction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 132; coefficient `25` anchored upstream at row 130 / stage 076)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.txt`

## What the paper claims

Stage 077 extracts the concrete Family-1 wall-depth datum `Theta_w` on the balanced reference branch (`alpha_r = 10`, `epsilon_r = 0.05`, `p_r = 2`). Carrying forward stage 076's lock `Theta_w = 25 lambda_mu^2 rho_w^2` (`eq:app-stage076-theta`), it weights the explicit radial Thomas–Fermi wall profile `rho_r(xi) = [1 - alpha_r S(xi)^2]_+^(1/4)` (with `S(xi)=(1+tanh xi)/2`) by the canonical support weight `chi_phi(xi)=S'(xi)=(1/2)sech^2 xi`. The `\stagefield{Output}` is the "Natural and conservative Family--1 wall-depth data" given by the two boxed equations: the natural shell-weighted datum `Theta_w^{(\chi)} = 25 lambda_mu^2 <rho_r^2>_chi ≈ 4.06863235008162 lambda_mu^2` and the conservative Jensen lower envelope `Theta_w^{(J)} = 25 lambda_mu^2 <rho_r>_chi^2 ≈ 0.927552032539308 lambda_mu^2`. The card also states the intermediate moments `<rho_r>_chi ≈ 0.192619005556493` and `<rho_r^2>_chi ≈ 0.162745294003265` (eq:app-stage077-rho-weights). The notes additionally pin the exact support cut point `xi_* = artanh(2/sqrt(alpha_r)-1) ≈ -0.3855810692` and the normalization `I_f = int chi^2 = 1/3`, and assert `Theta_w^{(\chi)} >= Theta_w^{(J)}` (Jensen).

## What the script claims to verify

Both scripts (a) symbolically derive `chi_phi=(1/2)sech^2 xi` and prove `I_f = int_{-oo}^{oo} chi^2 = 1/3`; (b) form `xi_* = atanh(2/sqrt(alpha_r)-1)` and prove it is the exact cut point by substituting back and confirming `1 - alpha_r S(xi_*)^2 = 0`; (c) numerically integrate the explicit `alpha_r=10` profile against the support weight to produce `<rho>_chi` and `<rho^2>_chi`, then `Theta_w^{(\chi)} = 25 <rho^2>_chi` and `Theta_w^{(J)} = 25 <rho>_chi^2`, each compared (~1e-26..1e-28 tol) against the same closed numeric values the card/notes report; and (d) assert the ordering `Theta_w^{(\chi)} >= Theta_w^{(J)} > 0`. This is exactly the stage's deliverable set.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `Theta_w = 25 lambda^2 rho^2` carry-in (stage 076) | coefficient `25` hardcoded in both `Theta_chi`/`Theta_J`; matches stage 076 card (row 130) and notes | match |
| `I_f = 1/3` (notes §2) | sympy L36-39 / wl L45-50 symbolic integral + `expect_zero(If-1/3)` | match |
| `xi_* = artanh(2/sqrt(alpha)-1) ≈ -0.38558106921542562404` (notes §1) | sympy L41-47 / wl L46-55: derive + back-substitute, `expect_zero(1 - alpha S(xi_*)^2)`; numeric `xi_*` printed | match |
| `<rho>_chi ≈ 0.192619005556493` | sympy L67/72 + `expect_close` L88-93; wl L71/86/92 + `expectApprox` | match |
| `<rho^2>_chi ≈ 0.162745294003265` | sympy L68/73 + `expect_close` L94-99; wl L75/87/93 + `expectApprox` | match |
| `Theta_w^{(\chi)} ≈ 4.06863235008162` | sympy L76/79 + `expect_close` L100-105; wl L82/89/94 | match |
| `Theta_w^{(J)} ≈ 0.927552032539308` | sympy L77/80 + `expect_close` L106-111; wl L83/90/95 | match |
| `Theta_w^{(\chi)} >= Theta_w^{(J)}` (Jensen, notes §5) | sympy L112-113; wl L96 | match |

`paper_alignment: aligned` — every paper deliverable has a faithful, non-tautological script-side check, and every script-emitted value reconciles (see Value Reconciliation section).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `expect_zero(If - 1/3)` (symbolic `int chi^2`) | `I_f = 1/3` | yes |
| A2 | sympy | 47 | `expect_zero(1 - alpha_r S(xi_*)^2)` (back-sub) | `xi_*` cut point | yes |
| A3 | sympy | 88-93 | `expect_close(R1, 0.1926190055…, 1e-28)` | `<rho>_chi` | yes |
| A4 | sympy | 94-99 | `expect_close(R2, 0.1627452940…, 1e-28)` | `<rho^2>_chi` | yes |
| A5 | sympy | 100-105 | `expect_close(25*R2, 4.0686323500…, 1e-26)` | `Theta_w^{(\chi)}` | yes |
| A6 | sympy | 106-111 | `expect_close(25*R1^2, 0.9275520325…, 1e-27)` | `Theta_w^{(J)}` | yes |
| A7 | sympy | 112-113 | `assert Theta_chi >= Theta_J > 0` | Jensen ordering | yes |
| B1 | mathematica | 50 | `expectZero(ifMom - 1/3)` (symbolic `Integrate`) | `I_f = 1/3` | yes |
| B2 | mathematica | 55 | `expectZero(1 - alphaR S(xi_*)^2)` (back-sub) | `xi_*` cut point | yes |
| B3 | mathematica | 92 | `expectApprox(r1, 0.1926190055…, 1e-28)` | `<rho>_chi` | yes |
| B4 | mathematica | 93 | `expectApprox(r2, 0.1627452940…, 1e-28)` | `<rho^2>_chi` | yes |
| B5 | mathematica | 94 | `expectApprox(25*r2, 4.0686323500…, 1e-26)` | `Theta_w^{(\chi)}` | yes |
| B6 | mathematica | 95 | `expectApprox(25*r1^2, 0.9275520325…, 1e-27)` | `Theta_w^{(J)}` | yes |
| B7 | mathematica | 96 | `expectTrue(thetaChi >= thetaJ && thetaJ>0)` | Jensen ordering | yes |

None of the rows are tautological. A1/B1 (symbolic `int chi^2 = 1/3`) and A2/B2 (back-substitution into `1 - alpha S^2`) are independent derivations, not assert-what-you-defined. A3–A6 compare a live-computed `mp.quad`/`NIntegrate` integral against a pinned target: if the profile, weight, cut point, or coefficient `25` were wrong the integral would diverge from the pinned value and the assertion would fail (not tautological — the target is the paper's reported value and the LHS is recomputed from first principles).

## Findings

### F1 — symbol_assumption_error

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py:33`

**What's wrong:**
The symbolic integration variable is declared positive:

```python
xi, alpha_r, lambda_mu = sp.symbols("xi alpha_r lambda_mu", positive=True, real=True)
```

But `xi` is the variable of the symmetric, full-line integral `If = sp.integrate(sp.simplify(chi**2), (xi, -sp.oo, sp.oo))` (line 36), and the physically relevant support cut point is `xi_* ≈ -0.3856 < 0` (the stage's own notes §1: "the active support weight lies mostly on the inner edge ... `xi_* ≈ -0.3855810692`"). Declaring `xi` `positive=True` contradicts the setup: `xi` ranges over all reals and the load-bearing point is negative. The Mathematica mirror correctly declares `Element[{xi, alphaR}, Reals]` (`.wl:41`), not positive — so the two engines disagree on the domain of `xi`.

This is non-load-bearing in the current script: the only symbolic use of `xi` is the definite integral with explicit `(-oo, oo)` bounds (the explicit bounds, not the symbol assumption, set the domain) and the result `1/3` is produced and independently confirmed by the correctly-assumed Mathematica engine (both outputs show `I_f = 1/3`). The numeric extraction uses mpmath floats, not the symbol, so it is unaffected. Hence no result is corrupted today.

**Why this matters:**
The declaration is a latent trap. Any future symbolic manipulation added to this block (e.g. an antiderivative, an `assuming`/`refine` simplification, a parity argument, or replacing the explicit-bounds integral with an indefinite one) would silently inherit `xi > 0` and could produce a wrong simplification or mask a sign/parity error — precisely the class of bug the second-engine policy exists to catch. It also makes the two engines model different domains for the same variable, weakening the independence claim.

**Required change:**
On `scripts/..._sympy_audit.py:33`, drop `positive=True` for `xi` (keep `real=True`); keep `positive=True` for `alpha_r` and `lambda_mu`, which are genuinely positive. Concretely, declare `xi` as a separate real symbol:

```python
xi = sp.symbols("xi", real=True)
alpha_r, lambda_mu = sp.symbols("alpha_r lambda_mu", positive=True, real=True)
```

**Verification:**
After the edit, `python3 scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py` still exits 0 and the saved output still prints `I_f = 1/3`, `I_f - 1/3 = 0`, and the unchanged numeric block (`<rho>_chi = 0.19261900555649309…`, `Theta_w^(chi) = 4.0686323500816155…`). The symbolic results must be identical because the integral bounds were already explicit; the change only removes the spurious positivity assumption.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. Evidence: (1) it derives a *different closed form* for the same cut point — `xi_* = -1/2*Log[-1+Sqrt[alphaR]]` (output L9) versus SymPy's `-atanh(1 - 2/sqrt(alpha_r))` (output L8) — both correct and equal, but reached via Mathematica's own `FullSimplify` of `ArcTanh`, not a copy of the SymPy expression. (2) It computes `<rho^2>_chi` via `rhoSqNum[x] = Sqrt[1 - alphaNum*sNum[x]^2]` (`.wl:65`, a direct `val^(1/2)`), whereas SymPy squares the quartic root `rho_num(x)**2 = (val**0.25)**2` (`.py:68`) — different code paths to the same quantity. (3) It uses `NIntegrate[..., WorkingPrecision->60, AccuracyGoal->28]` with `(-Infinity, xiCut)` bounds while SymPy uses `mp.quad(..., [-mp.inf, -4, xi_cut])` with an explicit interior split node at `-4`; the quadrature strategies differ. (4) The symbolic-integral engines differ (`Integrate` with `Element[...,Reals]` vs `sp.integrate` over `(-oo,oo)`). This is a genuine second engine.

## Engine cross-check

The two engines agree to full precision:

| quantity | sympy output | mathematica output |
|---|---|---|
| `xi_*(alpha=10)` | `-0.38558106921542562403635498846713378847348301441599` | `-0.385581069215425624036354988467133788473483014415991…` |
| `<rho>_chi` | `0.19261900555649309777068139356018510792903510747507` | `0.19261900555649309777068139356018510792903510747506717…` |
| `<rho^2>_chi` | `0.16274529400326462037087418498629868328210821103971` | `0.16274529400326462037087418498629868328210821103971427…` |
| `Theta_w^(chi)` | `4.0686323500816155092718546246574670820527052759928` | `4.06863235008161550927185462465746708205270527599285687…` |
| `Theta_w^(J)` | `0.92755203253930797183993260663904217023332624032789` | `0.92755203253930797183993260663904217023332624032789843…` |
| `denominator (I_f)` | `0.33333333333333333333333333333333333333333333333333` | `0.33333333333333333333333333333333333333333333333333268…` |

Both pass all checks (`I_f - 1/3 = 0`; cut-point residual `= 0`; all numeric diffs ~1e-50 << tolerances; ordering `True`). No engine disagreement.

## Verdict justification

The scripts faithfully verify exactly what the paper claims: the two symbolic identities (`I_f = 1/3`, the exact cut point `xi_*`) are non-tautological independent derivations, and the four numeric deliverables (`<rho>_chi`, `<rho^2>_chi`, `Theta_w^{(\chi)}`, `Theta_w^{(J)}`) are recomputed from the explicit profile and compared to the card/notes values to ~1e-50 by two independent engines, plus the Jensen ordering is asserted. The coefficient `25` matches the upstream stage 076 lock and the appendix. Every emitted value reconciles to the `.tex` and/or `.md` (Value Reconciliation: 9 values, 0 misaligned). Attacks that failed: (i) the numeric-vs-pinned comparisons are not tautological — the LHS is a live integral; (ii) the back-substitution `1 - alpha S(xi_*)^2 = 0` genuinely tests the `atanh` form, it does not assert a defined quantity; (iii) cutting the numerator integral at `xi_cut` while the denominator spans the full line is correct because `rho` is clamped to 0 beyond the cut (`[...]_+` in notes), so it equals the full-line integral of the clamped profile; (iv) the `25` coefficient and all reported decimals match the paper. The single finding (F1) is a real but non-load-bearing `symbol_assumption_error`: declaring the symmetric integration variable `xi` `positive=True` contradicts the setup (`xi_* < 0`, integral over `(-oo,oo)`) and disagrees with the Mathematica domain, but the explicit integral bounds neutralize it today. Verdict is `findings` (one low-severity script-side fix), not clean, and not stop_cold — fixing it cannot change any result.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 9 values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `chi_phi(xi) = (1/2)sech^2 xi` | py L37 / wl L48; out L5 | notes:56 `chi_phi(xi) = S'(xi) = (1/2) sech^2 xi` | MATCH |
| `I_f = 1/3` | py L38-39 / wl L49-50; out L6-7 | notes:59 `I_f = int dxi chi_phi(xi)^2 = 1/3` | MATCH |
| `xi_*` symbolic `= artanh(2/sqrt(alpha)-1)` | py L41 / wl L46; out L8-9 | notes:41 `xi_* = artanh( 2/sqrt(alpha_r) - 1 )` | MATCH |
| `xi_*(alpha=10) ≈ -0.38558106921542562404` | py L43 / wl L52; out L9/L11 | notes:45 `xi_* ≈ -0.3855810692` | MATCH |
| `<rho>_chi ≈ 0.192619005556493` | py L72 / wl L86; out L12 | tex:17 `\langle\rho_r\rangle_\chi\simeq0.192619005556493`; notes:87 | MATCH |
| `<rho^2>_chi ≈ 0.162745294003265` | py L73 / wl L87; out L13 | tex:18 `\langle\rho_r^2\rangle_\chi\simeq0.162745294003265`; notes:89 | MATCH |
| `Theta_w^{(\chi)} ≈ 4.06863235008162` | py L79 / wl L89; out L15 | tex:24-25 `\Theta_w^{(\chi)}=25\lambda_\mu^2\langle\rho_r^2\rangle_\chi\simeq4.06863235008162`; notes:93-94 | MATCH |
| `Theta_w^{(J)} ≈ 0.927552032539308` | py L80 / wl L90; out L16 | tex:31-32 `\Theta_w^{(J)}=25\lambda_\mu^2\langle\rho_r\rangle_\chi^2\simeq0.927552032539308`; notes:104-105 | MATCH |
| coefficient `25` (in `Theta_chi`, `Theta_J`) | py L76-77 / wl L82-83 | tex:24/31 `25\lambda_\mu^2…`; anchored upstream stage_076.tex:17 / appendix row 130 | MATCH |

INTERNAL (scaffolding, no finding): `denominator = 1/3` (numeric recomputation of `I_f`, the normalization; py L74/out L14, wl L88); the four `*_diff` residual values (~1e-50, tolerance scaffolding); the pass/fail flags. The numeric `<rho>`/`<rho^2>` are reported to ~15 digits in the card and to 50 digits in the script — the card's truncation is a legitimate terse rounding of the same value (MATCH, not MISMATCH).

## Self-test notes

(1) Variable independence: the only symbolic `sp.diff`/`D` is `chi = d/dxi S(xi)` and `S` does depend on `xi`, so the derivative is non-trivial (`= (1/2)sech^2 xi`) — no identically-zero-derivative trap. (2) Symmetry/parity: the denominator integrand `chi^2` is even on `(-oo,oo)` → finite `1/3` (correct, not a spurious zero); the numerator integrals are over `[-oo, xi_*]` only (asymmetric, by design, because `rho` is clamped to 0 beyond `xi_*`) so no parity-cancellation is claimed or needed. (3) Trivial-case: the F1 fix only removes `positive=True` from `xi`; the integral `int_{-oo}^{oo}(1/2 sech^2 xi)^2 = 1/3` is independent of any `xi` assumption (explicit bounds), so the assert still passes — confirmed against the Mathematica engine which already uses the correct (non-positive) domain and gets `1/3`. (4) Paper round-trip: the fix introduces no new constant and changes no value, so it cannot create a new paper_misalignment.
