---
unit_id: 075
batch: III.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: misaligned
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage075_family1_threshold_window.md
  paper_appendix: present
---

# Audit unit 075 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_075.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage075_family1_threshold_window.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (only the row mentioning stage 075 was inspected for status — the rest is a longtable of unrelated stages)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.txt`

## What the paper claims

The paper card declares Stage 075 as the explicit numerical evaluation of the Family-1 threshold window. The card's stated `\stagefield{Inputs}` (line 7) are `kappa = 12321/5`, `eta = 37`, and `Upsilon_w = 117 Theta_w`. The body equations give `Delta_0(12321/5, 37) ~ 1.73302079e-4` and `Delta_inf(12321/5, 37) ~ 2.01447566e-2` (lines 14-23). The `\stagefield{Output}` is the boxed pair of inequalities (lines 27-35):

> `Theta_w <= 3.62605617972939e-4 Pe_req => fail`
> `Theta_w >= 4.21495341569977e-2 Pe_req => succeed`

The notes file enumerates the same chain (Delta_0, Delta_inf, Upsilon_fail, Upsilon_suff, Xi_fail, Xi_suff) but states the wall-depth reduction differently: notes section 3 defines `V0 = alpha_r mu_*` with `alpha_r = 10` and concludes `Upsilon_w = alpha_r^2 Theta_w`, then states "On the balanced Family-1 reference branch, `Upsilon_w = 168 Theta_w`." Notes section 4 then computes `Theta_fail = Upsilon_fail / 168 ~ 3.62605617972939e-4 Pe_req` — but `0.0362605617972939 / 168 = 2.158e-4`, not `3.626e-4`. The notes' final `Theta` numbers are arithmetically consistent with dividing by 100 (i.e., `alpha_r^2 = 100`), not 168.

## What the script claims to verify

The SymPy script's docstring states four tasks: evaluate `Delta_0` and `Delta_inf` on the explicit branch (`alpha = sqrt(12321/5)`, `eta = 37`); compute the `Upsilon` and `Xi` thresholds; reduce the remaining amplitude via `Upsilon_w = alpha_r^2 Theta_w` with `alpha_r = 10`; and compute the `Theta_w` threshold window. The current assertions are: (i) two free-symbol "algebraic identities" for `Delta_0` and `Delta_inf` (lines 50-61 SymPy; lines 83-93 Mathematica), and (ii) two "round-trip" checks `Upsilon_fail - alpha_r^2 * Theta_fail == 0` and `Upsilon_suff - alpha_r^2 * Theta_suff == 0` (lines 99-104 SymPy; lines 95-96 Mathematica). The Mathematica script additionally compares its numeric outputs against eight literal targets (`expectApprox`, lines 98-105).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Delta_0(12321/5, 37) ~ 1.73302079e-4` (eq. app-stage075-Delta0) | sympy L31-34 / wl L43-46 plus numeric back-checks | match |
| `Delta_inf(12321/5, 37) ~ 2.01447566e-2` (eq. app-stage075-Deltainf) | sympy L35-38 / wl L47-50 plus numeric back-checks | match |
| `Theta_w <= 3.62606e-4 Pe_req => fail` (boxed) | sympy L78 / wl L65 numeric for `Theta_fail/Pe_req` | match (numeric) |
| `Theta_w >= 4.21495e-2 Pe_req => succeed` (boxed) | sympy L79 / wl L66 numeric for `Theta_suff/Pe_req` | match (numeric) |
| `Upsilon_w = 117 Theta_w` (paper Inputs line) | script uses `Upsilon_w = alpha_r^2 Theta_w` with `alpha_r = 10` → `100 Theta_w` | mismatch |
| `Upsilon_w = 168 Theta_w` (notes section 3) | script uses `Upsilon_w = 100 Theta_w` | mismatch |
| (none) `Upsilon_w = 100 Theta_w` consistent with boxed numerics | script uses 100 | extra (script-side numerics agree with paper's boxed Theta values, but no paper text states 100) |

Paper alignment: `misaligned` — the numeric `Theta_fail` and `Theta_suff` values in the script match the paper's boxed output exactly, but the **conversion factor stated in the paper's Inputs line (117)** and the **conversion factor stated in the notes (168)** both contradict the script's actual use of `alpha_r^2 = 100`, and neither external value is consistent with the paper's own boxed final numbers.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 50-60 | `simplify((alpha_sym*sinh + eta_sym*cosh) * Delta0_sym - eta_sym*(cosh-1)/alpha_sym^2) == 0` | Delta_0 closed form | no — tautological by definition (see F1) |
| A2 | sympy | 54-61 | `simplify((alpha_sym*sinh + eta_sym*cosh) * Deltainf_sym - (cosh + (eta_sym/alpha_sym)*sinh - 1)) == 0` | Delta_inf closed form | no — tautological by definition |
| A3 | sympy | 99-103 | `simplify(Upsilon_fail - alpha_r**2 * Theta_fail) == 0` where `Theta_fail := Upsilon_fail / alpha_r**2` | Upsilon_w = alpha_r^2 Theta_w | no — tautological: `Theta_fail` is defined as `Upsilon_fail/100`, so `100*Theta_fail = Upsilon_fail` by construction |
| A4 | sympy | 100-104 | same shape for `Upsilon_suff`/`Theta_suff` | same | no — same tautology |
| M1 | mathematica | 87-89 | `expectZero[(aSym*Sinh+eSym*Cosh)*delta0Sym - eSym*(Cosh-1)/aSym^2]` | Delta_0 closed form | no — same tautology as A1 |
| M2 | mathematica | 90-92 | `expectZero[(aSym*Sinh+eSym*Cosh)*deltaInfSym - (Cosh + (eSym/aSym)*Sinh - 1)]` | Delta_inf closed form | no — same tautology as A2 |
| M3 | mathematica | 95 | `expectZero["Upsilon_fail - alphaR^2 * Theta_fail", upsilonFail - alphaR^2*thetaFail]` | rescaling | no — `thetaFail := upsilonFail/alphaR^2`, so check is `upsilonFail - upsilonFail = 0` |
| M4 | mathematica | 96 | same for suff | rescaling | no — same tautology |
| M5-M12 | mathematica | 98-105 | eight `expectApprox` calls vs literal SymPy targets | Delta_0, Delta_inf, Upsilon/Xi/Theta numerics | no — literal targets copied from SymPy output (the v1 audit's F2 left these in place as "informational"; they remain a hardcoded-result issue) |

None of the substantive assertions are anchored. The two transcripts both pass, but the pass is content-free.

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:77-104`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:65-96`

**What's wrong:**
The v1 audit's F1 "round-trip" fix is structurally still tautological. The SymPy script defines (lines 77-79):
```
Theta = sp.symbols("Theta_w", positive=True, real=True)
Theta_fail = sp.simplify(Upsilon_fail / alpha_r**2)
Theta_suff = sp.simplify(Upsilon_suff / alpha_r**2)
```
then asserts (lines 99-104):
```
Upsilon_fail_from_Theta = sp.simplify(alpha_r**2 * Theta_fail)
Upsilon_suff_from_Theta = sp.simplify(alpha_r**2 * Theta_suff)
...
assert sp.simplify(Upsilon_fail - Upsilon_fail_from_Theta) == 0
assert sp.simplify(Upsilon_suff - Upsilon_suff_from_Theta) == 0
```
Since `Theta_fail := Upsilon_fail / alpha_r**2` by construction, `alpha_r**2 * Theta_fail = Upsilon_fail` algebraically, so the residual is identically zero no matter what `alpha_r` is or what `Upsilon_fail` is. The same defect appears in the Mathematica script at lines 65-66 (definition) and lines 95-96 (assertion).

Additionally, the v1 "free-symbol identity" check (A1/A2, M1/M2) is **also tautological**:
```
Delta0_sym = eta_sym * (sp.cosh(alpha_sym) - 1) / (
    alpha_sym**2 * (alpha_sym * sp.sinh(alpha_sym) + eta_sym * sp.cosh(alpha_sym))
)
...
delta0_identity = sp.simplify(
    (alpha_sym * sp.sinh(alpha_sym) + eta_sym * sp.cosh(alpha_sym)) * Delta0_sym
    - eta_sym * (sp.cosh(alpha_sym) - 1) / alpha_sym**2
)
assert delta0_identity == 0
```
The identity `(A) * (X/A) - X = 0` is an algebraic-fraction cancellation, not a non-trivial property of `Delta_0`. The relation tested is the *definition* of `Delta_0` itself (you derive the closed form by rearranging `(alpha*sinh+eta*cosh) * Delta_0 = eta*(cosh-1)/alpha^2`). The check confirms that `sympy.simplify` and `Mathematica.FullSimplify` can cancel a common denominator — a feature of the CAS, not a verification of the physics.

A real verification of `Delta_0` and `Delta_inf` would derive them from their physical origin (the linearized profile ODE Green's function evaluated at the Stage-41/42 thresholds, or a series-expansion check at small alpha, or an independent closed-form derivation from the upstream `Upsilon_fail = Pe_req / (Lambda_ell^2 Delta_inf)` definition). None of these is present.

**Why this matters:**
Stages 076-078 consume `Theta_fail` and `Theta_suff` from this stage. If the closed forms for `Delta_0` or `Delta_inf` have a wrong factor (numerator or denominator), neither engine catches it: both engines evaluate the same closed forms, the free-symbol "identity" trivially cancels, the round-trip is by definition, and the numeric back-checks (M5-M12) all use literals lifted from the SymPy transcript. The PASS verdict tells you only that the CAS can simplify a fraction; it tells you nothing about whether `Delta_0 = eta*(cosh-1)/(alpha^2*(alpha*sinh+eta*cosh))` is the right answer.

**Required change:**
Add a substantive check that exercises `Delta_0` and `Delta_inf` against an *independent* construction, not a rearrangement of their own definitions. Two acceptable paths:

(a) **Limit/asymptotic check.** From the notes section 2: `Delta_inf ~ 1/alpha` and `Delta_0 ~ eta / [alpha^2 (alpha + eta)]` in the large-alpha regime. Verify that `lim_{alpha → infinity} alpha * Delta_inf = 1` and `lim_{alpha → infinity} alpha^2 * (alpha + eta) * Delta_0 / eta = 1` (both as exact limits in SymPy with `sp.limit(..., alpha_sym, sp.oo)` and as `Limit[..., aSym -> Infinity]` in Mathematica with the appropriate assumption). These limits are non-trivial consequences of the closed form — if a factor in the numerator or denominator were wrong by a constant, the limit would not be 1.

(b) **Small-alpha series check.** Verify the leading-order term: `Delta_0` should behave as `eta / (alpha^2 * eta * (1) ) * (alpha^2/2) = 1/2` at leading order in alpha → 0? Compute `sp.series(Delta_0_sym, alpha_sym, 0, 3)` and compare against a hand-computed leading term derived from `cosh(x) ~ 1 + x^2/2` and `sinh(x) ~ x`. If both `(alpha*sinh+eta*cosh) ~ alpha^2 + eta` and `cosh - 1 ~ alpha^2/2`, then `Delta_0 ~ eta * (alpha^2/2) / (alpha^2 * (alpha^2 + eta)) → 1/2 * eta/(alpha^2 + eta) → 1/2` as alpha → 0. Assert this with `sp.limit(Delta_0_sym, alpha_sym, 0) == sp.Rational(1, 2)`.

Either (a) or (b) would catch a wrong factor. Apply the analogous independent construction in both engines.

Also, replace the `Upsilon_w = alpha_r^2 Theta_w` round-trip check with a numeric *independence* check: instead of defining `Theta_fail` from `Upsilon_fail` and then verifying `alpha_r^2 * Theta_fail == Upsilon_fail`, define `Theta_fail` *independently* (e.g., via the notes' Theta_w definition `Theta_w = 4 rho_w^2 mu_*^2 / (hbar^2 c_{s,w}^2)` evaluated against a constructed reference, or just by the numerical value `0.00036260561797293886969 Pe_req`) and verify that `alpha_r^2 * Theta_fail` recovers `Upsilon_fail`. Note this would also force a paper-side decision about which `alpha_r^2` value is intended.

**Verification:**
After Codex patches, the script transcripts must show: (i) the large-alpha limit checks (or small-alpha leading-term checks) for `Delta_0` and `Delta_inf` printing PASS; (ii) a `Upsilon_w` reduction check where `Theta_fail` is constructed independently of `Upsilon_fail` (e.g., from the numeric target or from the microscopic definition).

### F2 — paper_misalignment

**Severity:** high
**Subtype:** `value_mismatch`
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_075.tex:7` — quote: `\stagefield{Inputs}{\(\kappa=12321/5\), \(\eta=37\), and \(\Upsilon_w=117\Theta_w\).}`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:108` — quote: `Upsilon_w = 168 Theta_w.`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:124-128` — quote: `Theta_fail = Upsilon_fail / 168 ≈ 3.62605617972939e-4 * Pe_req,` (but `0.0362605617972939 / 168 = 2.158e-4`, NOT `3.626e-4`; the notes' arithmetic is internally inconsistent)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:24` — quote: `alpha_r = sp.Integer(10)`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:37` — quote: `alphaR = 10;`

**What's wrong:**
Three different conversion factors are in play for the same physical relation `Upsilon_w = alpha_r^2 Theta_w`:
- Paper card Inputs line (line 7): `Upsilon_w = 117 Theta_w`.
- Notes section 3 (line 108): `Upsilon_w = 168 Theta_w` (with `alpha_r = 10` stated explicitly on line 92; the conclusion "168" is unsourced inside the notes and inconsistent with `alpha_r = 10` → `alpha_r^2 = 100`).
- Script: `alpha_r = 10` → `Upsilon_w = 100 Theta_w` (line 24 SymPy, line 37 Mathematica).

The script's numeric final `Theta_fail = 3.62605617972939e-4 Pe_req` and `Theta_suff = 4.21495341569977e-2 Pe_req` *match* the paper's boxed final equations (lines 28, 34). Those final numbers can be reproduced *only* by dividing the Upsilon thresholds by 100, not by 117 and not by 168.

So the paper's Inputs line is internally inconsistent with the paper's own final boxed Theta values; the notes' "168" likewise contradicts both the script and the notes' own final Theta numbers. The script's `alpha_r = 10` is consistent with everyone's *final numbers* but inconsistent with the *stated conversion factors* in both the paper Inputs line and the notes.

This is a classic stale-text vs. live-script disagreement. Three resolutions are possible:
- Paper Inputs line should be `Upsilon_w = 100 Theta_w` (consistent with paper's boxed final numbers and the script).
- Paper's boxed Theta values should be recomputed using `Upsilon_w = 117 Theta_w` (i.e., `Theta_fail = 0.036260562 / 117 = 3.099e-4 Pe_req` and `Theta_suff = 4.215 / 117 = 3.602e-2 Pe_req`).
- Script's `alpha_r` value is wrong and should be changed to match whichever paper choice the user prefers.

The notes' "168" is almost certainly a stale artifact from an earlier `alpha_r` value; the notes' Theta-numeric values already match `alpha_r^2 = 100`, so the "168" is an isolated text error.

**Why this matters:**
Stages 076-077 derive `Theta_w` from a wall profile (paper appendix row 076: `Theta_w = 25 lambda_mu^2 rho_w^2`), then Stage 078 compares it to *this stage's* threshold window. If the conversion factor `alpha_r^2` is uncertain, the comparison gate is uncertain. Even though the script's numeric output happens to match the paper's boxed Theta numbers, a reader cannot tell which of 100, 117, or 168 is the load-bearing constant — and any downstream consumer that uses the paper Inputs line literally will be off by 17%.

**Required change:**
This is a `paper_misalignment` finding. Codex must not silently edit either side. The user must choose the resolution direction (see `## Resolve before fix_loop` in the directive).

**Verification:**
After user resolution and follow-up directive, either the paper card and notes are updated to state `Upsilon_w = 100 Theta_w` (no script change), or the script's `alpha_r` is changed and downstream numerics propagate.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:34-66`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:23-79`

**What's wrong:**
The v1 audit's F3 fix declared the free-symbol identity check (M1/M2) to be "the independent-derivation leg". But as documented in F1, that identity is itself tautological — `(A)*(X/A) - X = 0` simplifies by algebraic cancellation, not by any non-trivial property of the closed form. So F3 remains: the Mathematica script is structurally a line-by-line transliteration of the SymPy script (same symbol declarations, same closed-form expressions for `delta0`/`deltaInf`, same rescaling cascade for `upsilonFail`/`upsilonSuff`/`xiFail`/`xiSuff`/`thetaFail`/`thetaSuff`, same hardcoded constants), and the only "independent" check (M1/M2) doesn't actually exercise the physics.

The comment inserted by v1's F3 fix ("This identity check is the independent-derivation leg…") is misleading. Both engines compute identical closed forms and then back-check against numerics produced by the first one.

**Why this matters:**
Two engines computing the same closed form by the same syntactic recipe cannot disagree on the answer. The "engine cross-check" claim collapses. A wrong factor or sign in the stated closed form of `Delta_0` or `Delta_inf` would propagate identically through both engines.

**Required change:**
F3's fix is structurally tied to F1's fix. If F1 is corrected by adding a genuine independent check (e.g., the large-alpha limit `lim alpha * Delta_inf = 1`), and both engines compute that limit independently, F3 is discharged. Remove or rewrite the misleading "independent-derivation leg" comment in the Mathematica script (lines 78-82) so it accurately describes the new substantive check from F1.

**Verification:**
After F1's required change lands, the Mathematica script must independently compute the large-alpha (or small-alpha) limit of `Delta_0` / `Delta_inf` via `Limit[..., aSym -> Infinity]` or `Series[..., {aSym, 0, n}]` — Mathematica's symbolic engine working independently from SymPy's `sp.limit` / `sp.series`. The two engines now compute genuinely independent quantities (limits are evaluated by separate algorithms), and the comment should accurately describe this.

## Independent-derivation check (Mathematica)

After the v1 fix, the Mathematica script still uses the same closed-form recipe as SymPy. The only check that touches free symbols (M1/M2) is a definitional rearrangement that any CAS will simplify to 0 by algebraic cancellation. The Mathematica script is therefore still a transliteration; the "second engine" only confirms that two CAS implementations of the same syntactic expression produce the same floating-point output. See F1, F3.

## Engine cross-check

Both engines produce identical numeric values (within `1e-13` to `1e-18`). Mathematica transcript:
```
Delta_0 = 0.00017330207902152514905715619654992403
Delta_inf = 0.02014475655405215942710329560991777563
```
SymPy transcript:
```
Delta_0 (numeric) = 0.00017330207902152514906
Delta_inf (numeric) = 0.020144756554052159427
```
The agreement is by construction (identical closed forms, identical numeric substitutions). It is not evidence of independent derivation. No `engine_disagreement` finding.

## Verdict justification

The v1 audit identified the right three findings, but Codex's v1 patches did not actually discharge any of them: the "round-trip" check substitutes one tautology for another (Theta_fail is *defined* as Upsilon_fail/alpha_r^2, so the round-trip is trivial); the "free-symbol identity" check is the algebraic rearrangement that defines the closed form (so the identity is a CAS-level fraction cancellation, not a physics check); the "independent-derivation leg" comment in the Mathematica script labels a tautology as substantive. None of the assertions in either script can fail in a way that would reveal a wrong factor in `Delta_0` or `Delta_inf`. Layered on top is a clear `paper_misalignment` (F2): paper Inputs line says `Upsilon_w = 117 Theta_w`, notes say `168`, script uses `100`, and the script's final numerics match the paper's boxed Theta values — so the script and the boxed paper output agree on `100` while two text descriptions disagree with both. Verdict: `findings` (3). Attacks I tried: (i) verify whether the F1-replacement round-trip assertion can fail — no, it's `Upsilon_fail - 100*(Upsilon_fail/100)` by construction; (ii) verify whether the free-symbol identity check can fail — no, it's `(A)(X/A) - X = 0` purely algebraically; (iii) cross-check paper Inputs `117` against the script's `Theta_fail` — `0.0362605617972939 / 117 = 3.099e-4`, while the paper's box says `3.626e-4`, so the paper Inputs line is inconsistent with the paper's own final output; (iv) cross-check notes `168` — `0.0362605617972939 / 168 = 2.158e-4`, also inconsistent.

## Self-test notes

I checked: (a) the SymPy round-trip "fix" — substituted `Theta_fail = Upsilon_fail/100` into `Upsilon_fail - 100*Theta_fail` and got `Upsilon_fail - Upsilon_fail = 0` regardless of physics, confirming F1 still active. (b) The free-symbol identity check — `(alpha*sinh+eta*cosh) * Delta0_sym - eta*(cosh-1)/alpha^2` with `Delta0_sym := eta*(cosh-1)/(alpha^2*(alpha*sinh+eta*cosh))` reduces to `eta*(cosh-1)/alpha^2 - eta*(cosh-1)/alpha^2 = 0` after canceling the `(alpha*sinh+eta*cosh)` factor in numerator/denominator; confirmed tautology. (c) Arithmetic of paper Inputs vs. script: `0.0362605617972939 / 100 = 3.62605617972939e-4` (matches paper box); `/117 = 3.099e-4` (does not match); `/168 = 2.158e-4` (does not match) — confirming paper Inputs line `117` and notes `168` are both inconsistent with the paper's boxed final result. (d) Output mtimes: sympy script `1779513219`, sympy output `1779513365` (output fresher by 146s); mathematica script `1779513219`, mathematica output `1779513400` (output fresher by 181s) — no `stale_output` finding. (e) Proposed F1 fix (large-alpha limit `lim alpha * Delta_inf = 1`): with `Delta_inf = (cosh + (eta/alpha) sinh - 1)/(alpha*sinh + eta*cosh)`, large `alpha` gives `(e^alpha/2)(1 + eta/alpha)/((alpha + eta)(e^alpha/2)) = (1 + eta/alpha)/(alpha + eta) → 1/alpha` so `alpha * Delta_inf → 1`. Non-trivial. (f) For F2, deliberately did not propose a Codex-side fix; routed to user-resolution block in directive.
