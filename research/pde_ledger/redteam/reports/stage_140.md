---
unit_id: 140
batch: IV.5
auditor_model: claude-opus-4-7
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage140_selfmatched_mouth_susceptibility.md]
  paper_appendix: present
---

# Audit unit 140 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_140.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage140_selfmatched_mouth_susceptibility.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only includes `\input{stages/stage_140}` at line 1314 — no separate row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.txt`

## What the paper claims

The stage card declares the audit target verbatim: "Same-layer susceptibility gives `Sigma_0 = M_s = (20/9) That_m^2`". The card's inputs identify the construction as a coupled mouth fixed point and gain selection ledger step inside the reduced branch, importing the shell/mixed core, the mouth source law, outlet consistency, core-to-mouth gain maps, and self-matched susceptibility closure. The notes file makes the construction explicit: with `Theta_sigma = H_w J_s`, `g_s = T_m J_s`, `K_s = 3 pi a^2 hbar^2 / (5 m_psi rho_w ell)`, and `Sigma_0 = L g_s^2 / (K_s Theta_sigma)`, direct simplification yields `Sigma_0 = 20 L ell^2 rho_w^2 T_m^2 / (9 hbar^2 c_{s,w}^2)`, which under the normalized traction definition `That_m = rho_w ell sqrt(L) T_m / (hbar c_{s,w})` reduces to `Sigma_0 = (20/9) That_m^2`. The notes also list (Section 3) the numerically located fixed points `M_s^nat ≈ 1.66854252965624`, `M_s^comp ≈ 1.80594111095636`, the derived `That_m^nat ≈ 0.866512630228382`, `That_m^comp ≈ 0.901484054174206`, and the fractional enhancement ≈ 0.0403588161624. The card's checks list also flags M_q on the family-1 branch and outlet consistency as items to track, though only the `Sigma_0 = (20/9) That^2` identity is quoted as the audit target.

## What the script claims to verify

The sympy script (25 lines) builds `Sigma_0 = L * gs^2 / (Ks * Theta)` from independent definitions of `Js`, `Hw`, `Theta`, `Ks`, and `gs`, simplifies it, substitutes `T_m → That * hbar * cs / (rho * ell * sqrt(L))`, and asserts `simplify(Sigma0_hat - (20/9) * That^2) == 0` (line 16). It then prints the numerically located `That_nat`, `That_comp`, and the fractional enhancement, using hardcoded `Ms_nat` and `Ms_comp` from upstream (these are not asserted; only printed). The Mathematica script (59 lines) mirrors this exactly: same symbol choreography, same `sigma0 = lM*gS^2/(kS*thetaSigma)`, same substitution, and the same `expectZero["Sigma_0_hat - 20 That^2/9", sigma0Hat - (20/9)*tHat^2]` (line 44). The audit target stated by the card is therefore exercised by a non-trivial multi-step algebraic identity, not a tautology.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Sigma_0 = M_s = (20/9) That_m^2` (Output, body equation) | sympy line 16, mathematica line 44 | match |
| `That_m^nat ≈ 0.866512...`, `That_m^comp ≈ 0.901484...`, fractional enhancement ≈ 0.04036 (notes §3) | sympy lines 18-24, mathematica lines 46-53 (printed, not asserted) | partial (numerics printed but not asserted; values agree with notes to many digits) |
| Checks list "gain pair (M_s, M_q) against outlet consistency" | not exercised in this stage's scripts | missing (but card frames this as a downstream check, not as the audit target) |
| Checks list "numerical fixed points recorded as numerically located, not closed-form" | printed numerics in both scripts | match |

Dominant pattern: aligned. The Output box is what the card lists as the audit target ("Stage~140 is a coupled mouth fixed point and gain selection ledger step. Its audit target is the verification output quoted below."), and the audit target is the (20/9) That^2 identity, which both scripts exercise as a non-tautological assertion.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 16 | `assert sp.simplify(Sigma0_hat - sp.Rational(20,9)*That**2) == 0` | `Sigma_0 = (20/9) That^2` (Output) | yes |
| A2 | mathematica | 44 | `expectZero["Sigma_0_hat - 20 That^2/9", sigma0Hat - (20/9)*tHat^2]` | `Sigma_0 = (20/9) That^2` (Output) | yes |

The numerics on lines 18-24 (sympy) and 46-53 (mathematica) are reported via `print` only — no assertions. That matches the notes (Section 3) and the card's check item that these be reported as "numerically located, not closed-form constants", so the absence of an assertion is acceptable.

## Findings

### F1 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl:26-44`

**What's wrong:**
The Mathematica script reproduces the sympy script's algebraic choreography line-for-line, with identical variable order and identical algebraic steps. Compare:

sympy lines 6-16:
```
Js = 4*sp.pi*a**2*ell/3
Hw = m*cs**2/rho
Theta = Hw*Js
Ks = 3*sp.pi*a**2*hbar**2/(5*m*rho*ell)
gs = Tm*Js
Sigma0 = sp.simplify(L*gs**2/(Ks*Theta))
...
Sigma0_hat = sp.simplify(Sigma0.subs(Tm, That*hbar*cs/(rho*ell*sp.sqrt(L))))
assert sp.simplify(Sigma0_hat - sp.Rational(20,9)*That**2) == 0
```

mathematica lines 33-44:
```
jS = 4*Pi*a^2*ell/3;
hWall = mPsi*cSound^2/rhoW;
thetaSigma = hWall*jS;
kS = 3*Pi*a^2*hbar^2/(5*mPsi*rhoW*ell);
gS = tM*jS;
sigma0 = FullSimplify[lM*gS^2/(kS*thetaSigma), Assumptions -> $Assumptions];
...
sigma0Hat = FullSimplify[sigma0 /. tM -> tHat*hbar*cSound/(rhoW*ell*Sqrt[lM]), Assumptions -> $Assumptions];
expectZero["Sigma_0_hat - 20 That^2/9", sigma0Hat - (20/9)*tHat^2];
```

Additionally, the banner on `.wl` line 26 reads `"STAGE 123 — SELFMATCHED MOUTH SUSCEPTIBILITY CLOSURE"`, a copy-paste artifact whose presence further supports that the file was ported from another stage's `.wl` rather than independently composed.

**Why this matters:**
The two-engine policy expects each engine to derive the claim from physical premises independently, so that algebra bugs do not propagate identically. Here, the claim is a single closed-form algebraic identity (substitute four explicit definitions into `L*g^2/(K*Theta)` and simplify); there is genuinely only one derivation pathway, so the structural similarity is partly intrinsic to the problem. Severity is low because (a) the underlying CAS engines are independent and run their own `Simplify`/`FullSimplify` algorithms, and (b) any algebra slip would surface as a non-zero residual in either engine. The stale "STAGE 123" banner is purely cosmetic but should be corrected to "STAGE 140".

**Required change:**
- Correct the Mathematica banner at line 26 to read `"STAGE 140 — SELFMATCHED MOUTH SUSCEPTIBILITY CLOSURE"`.
- Optionally, give the Mathematica derivation an alternate ordering (e.g., precompute `J_s^2 / K_s` and `Theta_sigma` separately, or factor through dimensionless `lambda = rho ell / (hbar c_s) sqrt(L)` directly) so the algebraic path differs from the sympy script. This is a discretionary structural change; the math is already correct.

**Verification:**
After the banner fix, the printed banner in the Mathematica output transcript should read "STAGE 140 — ..." instead of "STAGE 123 — ...". The `expectZero` residual must remain 0.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.py:18-24`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl:46-53`

**What's wrong:**
The Section-3 numerics in the notes (`That_nat ≈ 0.866512630228382`, `That_comp ≈ 0.901484054174206`, fractional enhancement ≈ 0.0403588161624) are reproduced by the scripts via `print` only, with no `assert` / `expectZero` check. If a future edit perturbs the input `Ms_nat`, `Ms_comp` literals or the `sqrt(9*Ms/20)` formula, the script will still exit 0 and the regression will go unnoticed until the transcript is read by a human.

The notes explicitly box these values:
> `That_m^{nat,*} = sqrt(9 M_s^{nat,*} / 20) ≈ 0.866512630228382`
> `That_m^{comp,*} = sqrt(9 M_s^{comp,*} / 20) ≈ 0.901484054174206`
> `That_m^{comp,*}/That_m^{nat,*} - 1 ≈ 0.0403588161624`

so they are paper-side claims, not informal asides.

**Why this matters:**
The card's check list explicitly includes "Check numerical fixed points are recorded as numerically located, not closed-form constants" — the script does record them, but it does not check that the recorded values match the notes' boxed numbers. A cheap `assert abs(That_nat_value - sp.Rational('866512630228382',10**15)) < 1e-12` (or `expectClose`-style check) would close this gap.

**Required change:**
- sympy: after computing `That_nat`, `That_comp`, and the fractional enhancement, add three `assert sp.Abs(value - expected) < tol` checks where `expected` is taken from the notes' boxed numeric values and `tol ≈ 1e-12`.
- mathematica: add the corresponding `expectZero["...", N[value - expected, 30]]` or `If[Abs[...] > tol, fail[...], pass[...]]` checks.

**Verification:**
After the change, both scripts must still exit 0; the transcript should show three new `PASS:` (or assertion-pass) lines for the three numerics, and a deliberate perturbation of any input literal should make the corresponding assertion fail.

## Independent-derivation check (Mathematica)

See F1. The `.wl` script is a structural transliteration of the `.py` script (same variable assignment order, same intermediate quantities, same substitution path). The two engines do call their own simplifiers (`sp.simplify` vs. `FullSimplify[Together[Expand[...]], Assumptions -> ...]`), which provides genuine cross-engine confirmation that the algebraic identity holds. Severity is therefore kept low.

## Engine cross-check

Both engines arrive at identical results:

sympy transcript:
```
Sigma_0 = 20*L*T_m**2*ell**2*rho**2/(9*c_s**2*hbar**2)
Sigma_0 in terms of That = 20*That**2/9
That_nat = 0.866512630228381532619923771658
That_comp = 0.901484054174205568349216979026
fractional traction enhancement = 0.040358816162445119166
```

mathematica transcript:
```
Sigma_0 = (20*ell^2*lM*rhoW^2*tM^2)/(9*cSound^2*hbar^2)
Sigma_0 in terms of That = (20*tHat^2)/9
Sigma_0_hat - 20 That^2/9 = 0
PASS: Sigma_0_hat - 20 That^2/9
That_nat = 0.8665126302283815392959569642047049776046752355560060900671`30.
That_comp = 0.9014840541742055424613909744622661805033974943213822507436`30.
fractional traction enhancement = 0.04035881616244508127491719985972103265`20.
```

Closed-form `Sigma_0` matches; `That_nat`/`That_comp` match to ≈ 15 digits (the limit of the sympy 30-digit float used in the script, which is initialized from a string literal carrying only 17 digits of `Ms_nat`); fractional enhancement matches. Engines agree.

## Verdict justification

`findings` (not `clean`, not `stop_cold`). The audit target stated in the card — `Sigma_0 = M_s = (20/9) That_m^2` — is exercised by a non-tautological, multi-step algebraic identity in both engines and the identity holds. Both findings are low severity: F1 flags the Mathematica banner copy-paste and the line-for-line transliteration (the math is correct; this is a structural and cosmetic gap), and F2 flags the absence of assertions on the Section-3 numerics that the notes explicitly box. Attacks tried that failed: (a) re-derived `L*gs^2/(Ks*Theta)` by hand and reproduced `20 L ell^2 rho^2 T_m^2 / (9 hbar^2 cs^2)`, then the dimensionless `20/9 That^2`; (b) checked `That_nat = sqrt(9 * 1.66854252965624 / 20) ≈ 0.866512630228382` and `That_comp = sqrt(9 * 1.80594111095636 / 20) ≈ 0.901484054174206` matched the transcripts; (c) verified all positivity assumptions are physically justified (lengths, densities, sound speeds, moduli, normalized traction); (d) verified the substitution `T_m → That * hbar * cs / (rho * ell * sqrt(L))` is the inverse of the notes' boxed `That = rho ell sqrt(L) T_m / (hbar cs)`. Paper alignment is `aligned`; no `paper_misalignment` finding is warranted.

## Self-test notes

- Variable-independence trap: no `sp.diff` / `D[...]` calls in either script, so the "VAR not in EXPR" failure mode does not apply.
- Symmetry/parity trap: no integrals, so the odd-function vanishing failure mode does not apply.
- Trivial-case pre-check: substituting any specific positive values for `(L, a, ell, rho, T_m, hbar, c_s, m)` into `Sigma_0 = L*gs^2/(Ks*Theta)` cancels `a` and `m` and reduces to `20 L ell^2 rho^2 T_m^2 / (9 hbar^2 c_s^2)`; the boxed identity holds independently of those four variables, which the assertion correctly captures.
- Paper round-trip: F2's proposed assertions use the notes' boxed values directly, so they do not introduce a new constant the paper does not state.
