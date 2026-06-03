---
unit_id: 251
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 251 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_251.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (rows at lines 100, 304-319)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 251 replaces the Session-IV phenomenological `\gamma_{\rm tot}\dot V` envelope law with the first *microscopic* odd export kernel seen by the active dressing leg `V`. The card (Verification line: "SymPy audit: …; Mathematica audit: none yet") states the deliverables: (1) the cubic coefficient `\Gamma_{3,0}=\gamma_1\eta_0^2\Omega_{U,0}^4/\Delta_0^2` derived from the derivative-coupled scalar outlet via the Stage-021 transfer factor `N_0(\omega)`, projected as `\Gamma_3=\Pi_{V0}^2\Gamma_{3,0}`; (2) the quintic coefficient `\Gamma_{5,-}=\frac{a^5}{27c_s^5}P_{0,-}`, `P_{0,-}=\beta_0 s_-/\lambda_-`, imported from the selected mixed quadrupole lane, projected as `\Gamma_5=\Pi_{V-}^2\Gamma_{5,-}`; (3) the super-Ohmic kernel `K_{\rm exp}(s)=\Gamma_3 s^3+\Gamma_5 s^5`; (4) the exact Schott power identity yielding positive-definite exported power `\mathcal P_{\rm exp}=\Gamma_3\ddot V^2+\Gamma_5\dddot V^2\ge0`; (5) the characteristic polynomial `F(s)=\Gamma_5 s^5+\Gamma_3 s^3+\mu_\eta s^2-\kappa_V` with `F(0)<0`, `F\to+\infty`, `F'(s)>0` for `s>0`, hence exactly one positive growth root (no finite microscopic `\gamma_{\rm crit}`); (6) the small-kernel slowdown `s_+=s_0-(\Gamma_3 s_0^2+\Gamma_5 s_0^4)/(2\mu_\eta)+O(\Gamma^2)`; (7) the event-safe half-plane `\widehat\Gamma_3+s_c^2\widehat\Gamma_5\ge(s_0^2-s_c^2)/s_c^3`; and (8) the Session-IV cold-event benchmark `\widehat\Gamma_3+0.3013336471\,\widehat\Gamma_5\ge289.61004918`, with `\widehat\Gamma_{3,\rm safe}\approx289.61004918`, `\widehat\Gamma_{5,\rm safe}\approx961.09429528`. The notes (Section 9) explicitly flag `\Pi_{V0}`, `\Pi_{V-}`, `\gamma_1`, and `P_{0,-}` (hence `\Gamma_{5,-}`) as still-open, branch-dependent objects this stage does not compute — they are carried as symbolic placeholders.

## What the script claims to verify

The docstring lists seven verified objects (cubic coeff, quintic coeff, projected kernel, Schott identities, characteristic polynomial + monotonicity data, event-safe surface, benchmark numbers). The assertions actually exercise: the omega^2 coefficient of `N_0` matches `\eta_0^2\Omega_{U,0}^4/\Delta_0^2` (line 59, a genuine limit extraction); the projected quintic equals its own factored form (line 81); the Schott power balance reduces to `-\Gamma_3\ddot V^2-\Gamma_5\dddot V^2` (line 109, genuine differential-algebra identity); `F'(s)` form (line 137); the solve-for-`\Gamma_3` safe-surface forms (lines 138-139); six numerical benchmark identities (lines 178-183); and that the normalized characteristic polynomial at the pure-cubic and pure-quintic safe coefficients has exactly one positive root equal to `s_c` (lines 204-207, the strongest checks). The small-kernel slowdown expression and the dimensionless/half-plane restatements are printed but not asserted.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) cubic `\Gamma_{3,0}=\gamma_1\eta_0^2\Omega_{U,0}^4/\Delta_0^2`, `\Gamma_3=\Pi_{V0}^2\Gamma_{3,0}` | line 49-59 limit extraction from `N_0` + boxed-form assert | match |
| (2) quintic `\Gamma_{5,-}=a^5/(27c_s^5)P_{0,-}`, `\Gamma_5=\Pi_{V-}^2\Gamma_{5,-}` | line 69-81 defines then re-asserts product | partial (tautological; see F2) |
| (3) `K_{\rm exp}(s)=\Gamma_3 s^3+\Gamma_5 s^5` / `\Sigma_{\rm exp}(\omega)` | line 87 builds `Sigma_exp`, printed, no assert | partial (display-only, low concern) |
| (4) Schott `\mathcal P_{\rm exp}=\Gamma_3\ddot V^2+\Gamma_5\dddot V^2\ge0` | line 98-109 power-balance assert | match |
| (5) `F(s)`, `F'(s)>0`, one positive root | line 119-137 `F`/`F'` form; line 204-207 root-count at two safe pts | partial (monotonicity *form* + 2 sample root-counts; `F(0)<0`, `F\to+\infty`, general uniqueness not asserted) |
| (6) slowdown `s_+=s_0-(\Gamma_3 s_0^2+\Gamma_5 s_0^4)/(2\mu_\eta)+O(\Gamma^2)` | line 126 `root_shift` printed, **no assert, not derived from F** | missing (see F3) |
| (7) safe half-plane `\widehat\Gamma_3+s_c^2\widehat\Gamma_5\ge(s_0^2-s_c^2)/s_c^3` | line 121-139 solve + asserts; line 204-207 root-at-`s_c` | match |
| (8) benchmark 289.61.../961.09... | line 156-183 numerical asserts | match |
| (—) `K_{\rm exp}` derived independently of `\gamma_{\rm tot}` (no circular check against the law it replaces) | script never references `\gamma_{\rm tot}` | match (the warned circularity is absent) |

`paper_alignment: partial` — no value disagreements (no `paper_misalignment`); the gaps are missing/tautological script-side coverage of deliverables (2), (5)-partial, and (6), all script-side fixable.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `simplify(N0_coeff2 - eta0^2 OmU0^4/Delta0^2)==0` | claim 1 (cubic) | yes (limit extracted from independent `N_0`) |
| A2 | sympy | 81 | `simplify(Gamma5 - PiVm^2 a^5 beta0 sminus/(27 cs^5 lamminus))==0` | claim 2 (quintic) | no (tautological re-expansion of defined product) |
| A3 | sympy | 109 | `power_balance == -G3 V''^2 - G5 V'''^2` | claim 4 (Schott) | yes (genuine identity) |
| A4 | sympy | 137 | `simplify(Fprime - (5 G5 s^4+3 G3 s^2+2 mu s))==0` | claim 5 (monotonicity form) | yes (form only; positivity not asserted) |
| A5 | sympy | 138 | `simplify(safe_eq - mu(s0^2-sc^2-G5 sc^5/mu)/sc^3)==0` | claim 7 (safe surface) | yes (solve verified) |
| A6 | sympy | 139 | `simplify(safe_eq_hat - ((s0^2-sc^2)/sc^3 - (G5/mu) sc^2))==0` | claim 7 (safe surface) | yes |
| A7 | sympy | 178 | `abs(s0_from_t - s0) < 1e-6` | claim 8 (benchmark consistency) | yes |
| A8 | sympy | 179-183 | five `abs(... - <literal>) < tol` | claim 8 (benchmark numbers) | yes (literals match paper card) |
| A9 | sympy | 204-207 | `len(pos)==1`, `abs(root - sc)<1e-9` (cubic & quintic) | claims 5 (uniqueness) + 7 (safe surface marginality) | yes (genuine `nroots`) |
| — | sympy | 126 | `root_shift = ...` (no assert) | claim 6 (slowdown) | NO ASSERT |
| — | sympy | 87,146-147 | `Sigma_exp`, `half_plane`, `dimless_safe` (no assert) | claims 3, 7 (display) | NO ASSERT |

16 asserts total vs. 7 docstring items — count is adequate; the defect is coverage, not quantity.

## Findings

### F1 — missing_verification_script (subtype: missing_mathematica)

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_mathematica_audit.wl` (does not exist)

**What's wrong:**
Stage 251 is SymPy-only. The card's Verification line itself states "Mathematica audit: none yet." Every claim this stage makes is fully reachable with native Mathematica primitives: the cubic coefficient via `Series[N0,{omega,0,4}]` (a *different* decomposition than the SymPy `Limit[N0/omega^2,...]` route), the Schott identity via `D[...]`, the characteristic-polynomial monotonicity/sign data via `Reduce`/`Resolve` over the positive orthant, the safe surface via `Solve[F[sc]==0,Gamma3]`, the small-kernel slowdown via `Series[F[s0+eps],{Gamma...,0,1}]`, and the benchmark roots via `NSolve`. None of this is impossible; therefore the dual-engine rule requires an independent `.wl`.

**Why this matters:**
With a single engine, every coefficient and identity rests on one CAS's algebra; the second engine is the cross-check that catches a SymPy `series`/`limit` artifact or a sign convention. It is also the natural place to verify the parts the SymPy script leaves only display-level (uniqueness via `Reduce`, `F(0)<0`, `F\to+\infty`) and to re-derive the cubic by a genuinely different route.

**Required change:**
Add a new independent-route Mathematica audit at the path above. See the directive's claim manifest M1-M7. Must use native Mathematica primitives and a different decomposition than the `.py` (e.g., `Series` not `Limit` for the cubic; `Reduce` for monotonicity/uniqueness), not a transliteration.

**Verification:**
`redteam exec-mathematica 251` produces the `.wl`, it exits 0, and the saved output shows `expectZero`/`expectTrue` lines covering M1-M7.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.py:69-81`

**What's wrong:**
The quintic coefficient is *defined* from its final factors and then *asserted equal* to those same factors:
```
69  P0_minus    = beta0 * sminus / lamminus
70  Gamma2_port = a**5 / (27 * cs**5)
71  Gamma5_minus = Gamma2_port * P0_minus
73  Gamma5      = PiVm**2 * Gamma5_minus
81  assert sp.simplify(Gamma5 - PiVm**2 * a**5 * beta0 * sminus / (27 * cs**5 * lamminus)) == 0
```
The assertion is `(PiVm^2 * a^5/(27 cs^5) * beta0 sminus/lamminus) - (PiVm^2 a^5 beta0 sminus/(27 cs^5 lamminus)) == 0`, which is identically zero by construction and cannot fail for any physics. Contrast with the cubic (line 49-59), which is genuine: `N0_coeff2` is *extracted* (`sp.limit(N0_scalar/omega^2,...)`) from an independently built rational function and then compared to the boxed form, so that assertion *can* fail if the transfer factor or `\Delta_0` is wrong.

This is NOT a `paper_misalignment`: the paper card and notes (Section 9) explicitly carry `\Gamma_{5,-}=a^5/(27c_s^5)P_{0,-}` and `P_{0,-}=\beta_0 s_-/\lambda_-` as *imported placeholders* from the selected mixed quadrupole lane, "still missing" on the realized branch. So the script faithfully reproduces the paper's symbolic form; the defect is that the check exercises no derivation — it is bookkeeping dressed as verification.

**Why this matters:**
The script's docstring item 2 claims it verifies "the exact quintic coefficient imported from the selected outgoing quadrupole lane." A reader trusts that the `a^5/(27c_s^5)` prefactor and the `\Gamma_5=\Pi_{V-}^2\Gamma_{5,-}` projection structure were checked, but nothing was. If a later revision mistyped the prefactor (say `a^5/(27c_s^4)`), the assertion would silently still pass because it re-derives from whatever was typed.

**Required change:**
Make the assertion non-tautological by checking the *structure the paper asserts*, not the literal you just built. Concretely, verify the projection homogeneity and the dim-consistency the paper requires without reusing the assembled `Gamma5` symbol on both sides — e.g., assert the projection scales as `\Pi_{V-}^2` (`Gamma5.subs(PiVm, k*PiVm) - k**2*Gamma5 == 0` for a free symbol `k`), assert `Gamma5_minus` factors as `Gamma2_port*P0_minus` where `Gamma2_port` is the *independently re-parsed* `a**5/(27*cs**5)` from the kernel `K_exp = Gamma3*s**3 + Gamma5*s**5` coefficient of `s**5`, and assert `P0_minus` is first-order in each of `beta0`, `sminus`, `1/lamminus`. This turns a self-equality into checks that can actually fail if the prefactor power or projection power is wrong. (Keep the SymPy symbolic form identical to the paper — do not change any value.)

**Verification:**
The diff shows line 81 replaced by ≥2 assertions that reference distinct constructions (e.g., `Poly(K_exp, s).coeff_monomial(s**5)` vs. the factored form, plus the `\Pi_{V-}^2` homogeneity scaling). Output exits 0.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.py:126,135`

**What's wrong:**
The small-kernel slowdown is a *boxed paper deliverable* (card eq `app-part08-stage251-small-kernel`, notes Section 4.2):
`s_+ = s_0 - (\Gamma_3 s_0^2 + \Gamma_5 s_0^4)/(2\mu_\eta) + O(\Gamma^2)`. The script only *defines and prints* the shift expression:
```
126  root_shift = sp.simplify(-(G3 * s0**2 + G5 * s0**4) / (2 * mu_eta))
135  print("delta s (weak export) =", root_shift)
```
There is no assertion, and the expression is hand-written rather than derived from `F(s)`. It "verifies" nothing — it merely restates the paper's answer. (This is not a value disagreement, so it is `insufficient_verification`, not `paper_misalignment`; the fix is a pure script-side addition.)

**Why this matters:**
The slowdown formula is the quantitative content of Section 4.2 ("both channels slow the collapse, quintic weighted by extra `s_0^2`"). If the formula were wrong (e.g., a missing factor of 2 or the wrong power of `s_0`), the script as written would not catch it.

**Required change:**
Derive and assert the slowdown from the characteristic polynomial. Substitute `\kappa_V = \mu_\eta s_0^2` and `s = s_0 + \varepsilon` into `F`, treat `\Gamma_3,\Gamma_5` as the small parameter (introduce a bookkeeping scalar `g` multiplying both, i.e. `G3 -> g*G3`, `G5 -> g*G5`), solve `F=0` for `\varepsilon` to first order in `g`, and assert the `O(g^1)` part equals `-(G3*s0**2 + G5*s0**4)/(2*mu_eta)`. The audit self-test below confirms this reduces to exactly the boxed coefficient.

**Verification:**
A new `assert sp.simplify(eps_first_order - (-(G3*s0**2 + G5*s0**4)/(2*mu_eta))) == 0` appears near line 126 and exits 0; output prints the derived (not hand-written) coefficient.

## Independent-derivation check (Mathematica)

No `.wl` exists, so `mathematica_transliteration` does not apply. The directive's new `.wl` (F1) must be an independent derivation: cubic via `Series` (not `Limit`), monotonicity/uniqueness via `Reduce`/`Resolve` (not just the algebraic `F'` form), and the slowdown via `Series` in the kernel strength. Filename: `moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_mathematica_audit.wl` in `mathematica/`.

## Engine cross-check

Only one engine present; `engines_agree: n/a`. The SymPy output (exit 0, all checks pass) is internally consistent and the benchmark numbers reproduce the paper card (`289.61004917557426` vs card `289.61004918`; `961.094295282802` vs card `961.09429528`; weight `0.3013336470698294` vs card `0.3013336471`).

## Verdict justification

`findings`. No `paper_misalignment`: every numeric constant and symbolic form in the script matches the card and notes, and the warned circularity (checking the microscopic kernel against the `\gamma_{\rm tot}` law it replaces) is absent. The genuine, load-bearing checks hold under attack: the cubic coefficient is a real limit extraction (A1), the Schott power identity is a real differential-algebra identity (A3), the safe surface is confirmed by an independent `nroots` root-at-`s_c` computation (A9), and the benchmark arithmetic matches the card. Three real defects remain: (F1) the dual-engine rule requires an independent Mathematica route, which is fully possible here; (F2) the quintic-coefficient assertion is a tautological re-expansion of a self-defined product and exercises no physics; (F3) the boxed small-kernel slowdown is printed but neither derived nor asserted. None is `UNFIXABLE` or `CRITICAL_DOWNSTREAM` — F2/F3 are script-side additions that do not change any quoted value, and F1 adds an engine; no downstream constant changes.

## Self-test notes

Variable-independence trap: F3's proposed check differentiates/series-expands `F(s_0+eps)` in the *kernel strength* `g`, which `F` genuinely depends on, and in `eps`; no derivative is taken w.r.t. an absent variable, so no vacuous identically-zero trap. Trivial-case pre-check: I hand-verified F3 — `F(s_0)=Gamma_3 s_0^3+Gamma_5 s_0^5`, `F'(s_0)|_{g=0}=2 mu_eta s_0`, so `eps = -F(s_0)/F'(s_0) = -(Gamma_3 s_0^2+Gamma_5 s_0^4)/(2 mu_eta)`, matching the box exactly. Symmetry/parity: no unbounded integrals introduced. Paper round-trip: F2 and F3 keep the paper's symbolic forms and constants verbatim (only the *checks* change), so no new paper_misalignment is introduced; the new `.wl` (F1) is instructed to use the same constants the card states.
