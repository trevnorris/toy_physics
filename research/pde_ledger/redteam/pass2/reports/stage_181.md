---
unit_id: 181
batch: V.2
auditor_model: Claude Opus 4.8 (1M context)
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage181_coherent_tracking_defect.md]
  paper_appendix: present
---

# Audit unit 181 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_181.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage181_coherent_tracking_defect.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 93, 265, 721–817 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.txt`

## What the paper claims

Stage 181 substitutes the *actual* coherent local D/N tracking branch into the
Stage 180 selected-branch transfer-shape law and proves the support-blindness
theorem. The card's `\stagefield{Output}` states verbatim: "On the coherent local
D/N tracking branch, \(\Xi_1\) is carried only by mixed/outgoing placement
variables; coherent support ratio \(\zeta\) drops out." The notes + appendix
(eqs `coherent-T`, `epsilon-split`, `Xi1-coherent-physical`, `epsilon1-coherent`,
`support-blindness`) enumerate the distinct deliverables:
(D1) the transfer-shape identity `T² = Z_W(1+χ₀)²/(Ω_W²(1-ε)²) = (27π²Gc_s⁵/20a⁵c⁵)(1-ε_η)/R_target`;
(D2) the split-blocking ratio `ε = ε_W(1 - (2/11)δ_U/(1+δ_U))`;
(D3) the split-blocking drift `ε₁ = (1-(2/11)δ_U/(1+δ_U))ε_W - (2ε_W/(11(1+δ_U)²))δ_{U,1}`;
(D4) the defect law `Ξ₁ = ζ_Z - ω_W + 2χ₁/(1+χ₀) + 2ε₁/(1-ε)`;
(D5) support-blindness `∂_ζ ln T² = ∂_ζ ln R_target = ∂_ζ Ξ₁ = 0`;
(D6) the selected-branch identity `R₁ = ω_W - η₁/(1-ε_η) - ζ_Z - 2χ₁/(1+χ₀) - 2ε₁/(1-ε)` with `Ξ₁ = -η₁/(1-ε_η) - R₁` (notes §4);
(D7) the tracking-factor drift `Θ₁ = -(χ₀(1+χ₀)δ_{U,1}+δ_U(1+δ_U)χ₁)/((1+χ₀)(1+δ_U)(1+χ₀+δ_U))` plus the claim that `Θ₁=0` does not force `Ξ₁=0` (notes §5).

## What the script claims to verify

Both engines build the coherent-branch objects (`Lam`, `eps`, `R_target`,
`T2_direct`, the support-loaded mass route) and assert: (i) the direct↔selected
transfer-shape identity (D1); (ii) a support-loaded reconstruction of `R_target`
and `T²` from a `ζ`-bearing route (D5 substrate); (iii) `∂_ζ ln T² = ∂_ζ ln
R_target = 0`, with a deliberately *spoiled* support packet confirming the check
can fail (D5); (iv) the split-blocking drift `eps_1` against a hand-written
closed form (D3); (v) `R_1` derived by `s`-perturbing `R_target` matches the
closed form (D6); (vi) `Xi_1` derived by `s`-perturbing `T²` matches the defect
law, plus `∂_ζ Xi_1 = 0` through the `ζ`-bearing loaded route (D4, D5);
(vii) the tracking-factor drift `Theta_1` against a hand-written closed form (D7);
and (viii) a support-rigid specialization `χ₁=δ_{U,1}=0` leaving a nonzero `Xi_1`
while `Theta_1=0` (D7 non-sufficiency). The docstring enumerates exactly these.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 T² identity | `T2_direct - T2_selected` (py:54 / wl:46) | match |
| D2 ε split form | used in `eps` def (py:41 / wl:33); exercised via D3 | match |
| D3 ε₁ drift | `split-blocking drift eps_1` vs closed form (py:81 / wl:79) | match |
| D4 Ξ₁ defect law | `Xi_1 derived from T^2 matches defect law` (py:114 / wl:99) | match |
| D5 support-blindness | `d/dzeta ln T²`, `ln R_target`, `Xi_1` (py:58,59,123 / wl:50,51,109) + spoiler | match |
| D6 selected-branch R₁ / identity | `R_1 derived...` (py:103/wl:90) + `selected-branch identity` (py:105/wl:91) | match |
| D7 Θ₁ + non-sufficiency | `tracking-factor drift` (py:133/wl:119) + support-rigid residual (py:143/wl:128) | match |

Every paper deliverable maps to a substantive, non-tautological script check.
`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 54 | `simplify(T2_direct - T2_selected)==0` | D1 | yes |
| A2 | sympy | 55 | `simplify(R_target_loaded - R_target)==0` | D5 substrate | yes |
| A3 | sympy | 56 | `simplify(T2_loaded - T2_direct)==0` | D5 substrate | yes |
| A4 | sympy | 57 | `R_target_loaded*Mtr - product_loaded==0` | none (round-trip) | NO — tautological |
| A5 | sympy | 58 | `diff(log(T2_loaded),zeta)==0` | D5 | yes |
| A6 | sympy | 59 | `diff(log(R_target_loaded),zeta)==0` | D5 | yes |
| A7 | sympy | 66 | spoiled `dlnR != 0` (counter-check) | D5 discriminating power | yes |
| A8 | sympy | 81 | `eps1 - eps1_expected==0` | D3 | yes |
| A9 | sympy | 103 | `R1_derived - R1==0` | D6 | yes |
| A10 | sympy | 105 | `Xi1 + eta1/(1-eps_eta) + R1==0` | D6 | yes |
| A11 | sympy | 114 | `Xi1_derived - Xi1==0` | D4 | yes |
| A12 | sympy | 123 | `diff(Xi1_loaded, zeta)==0` | D5 (Ξ₁) | yes |
| A13 | sympy | 133 | `Theta1 - Theta1_expected==0` | D7 | yes |
| A14 | sympy | 143 | `Xi_support_rigid != 0` (counter-check) | D7 non-sufficiency | yes |
| B1–B14 | mathematica | 46–128 | same forms (`FullSimplify[...]===0`, `D[...]`, `fail` guards) | same claims | yes (except B-analog of A4 at wl:49 — tautological) |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl:1-141`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py:1-151`

**What's wrong:**
The `.wl` is a line-for-line port of the `.py`, not an independent re-derivation.
The two scripts share: identical symbol names for the drift coefficients
(`zetaZ, omegaW, chi1, epsW1, deltaU1, eta1, s`); identical intermediate
definitions in the same order (`Lam`/`lamNorm`, `eps`, `R_target`/`rTarget`,
`T2_direct`/`t2Direct`, `Mmix`/`mMix`, `Msupp`/`mSupp`, `Ssupport`/`sSupport`,
`product_loaded`/`productLoaded`); the identical assertion list in the identical
order with identical printed names (e.g. "split-blocking drift eps_1",
"R_1 derived from R_target matches closed form", "Xi_1 derived from T^2 matches
defect law", "tracking-factor drift"); the identical `s`-perturbation route
(`epsPert`, `rTargetPert`, `t2DirectPert`, `t2LoadedPert`); the identical
`bad`-spoiler counter-check; and — most tellingly — the *identical hand-written
closed forms* fed to the comparison checks:
- py:77–80 `eps1_expected = (1 - Rational(2,11)*deltaU/(1+deltaU))*epsW1 - (Rational(2,11)*epsW/(1+deltaU)**2)*deltaU1`
- wl:75–78 `eps1Expected = (1 - (2/11)*deltaU/(1 + deltaU))*epsW1 - (2*epsW*deltaU1)/(11*(1 + deltaU)^2)`
and
- py:129–132 `Theta1_expected = -(chi0*(1+chi0)*deltaU1 + deltaU*(1+deltaU)*chi1)/((1+chi0)*(1+deltaU)*(1+chi0+deltaU))`
- wl:114–118 same.
Both engines compute `eps1`/`theta1` by the same `D[..,epsW]*epsW1 + ...`
construction and compare to the same transcribed target. There is no independent
route (no series expansion vs. direct, no alternate parameterization, no closed
form re-derived a second way). The only place the algebra diverges is the
support-loaded mass: SymPy uses `Mtr = Mmix + Msupp` (py:49) while Mathematica
uses `loadMassFromSupport = mMix*sSupport` (wl:42) — a genuine additive-vs-
multiplicative-via-S difference — but this is a single local divergence inside an
otherwise verbatim transliteration.

**Why this matters:**
The second-engine policy requires Mathematica to derive the result independently
from the physical premises, not echo SymPy's algebra. As written, an error
transcribed identically into both `eps1_expected`/`Theta1_expected` closed forms
would pass both engines. The cross-engine agreement is therefore weaker evidence
than it appears.

**Required change:**
Re-author the `.wl` so at least the load-bearing closed forms (D3 `ε₁`, D7 `Θ₁`,
D4 `Ξ₁`) are obtained by a route the SymPy script does not use — e.g. derive `ε₁`
and `Θ₁` via `Series[..., {s,0,1}]` Coefficient extraction of the perturbed
`eps`/`Rtr` (an `s`-expansion the SymPy `eps1` check does NOT use, which currently
uses direct `D[eps,epsW]*epsW1`), and compare that series-extracted coefficient to
the perturbed-`T²`/`R_target` route rather than to a hand-transcribed literal.
Keep the same physical premises; change the derivation path so the two engines are
genuinely independent. Do NOT merely rename variables.

**Verification:**
The verifier confirms the `.wl` no longer contains `eps1Expected`/`theta1Expected`
hand-written literals identical to the `.py`, and that the `ε₁`/`Θ₁`/`Ξ₁` checks
are anchored to an independently-derived (series-extracted) object; output still
exits 0 with all checks at 0.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py:57`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl:49`

**What's wrong:**
The "support-loaded branch product law" check is `R_target_loaded*Mtr -
product_loaded == 0` (py:57) / `rTargetLoaded*loadMassFromSupport - productLoaded
== 0` (wl:49). But `R_target_loaded` is *defined* as `product_loaded / Mtr`
(py:51) / `productLoaded / loadMassFromSupport` (wl:43). Substituting the
definition, the check reduces to `product_loaded - product_loaded == 0`, which is
algebraically guaranteed by construction and cannot fail for any physics. It is a
round-trip of the immediately-preceding definition.

**Why this matters:**
The check appears in the transcript as a PASS line ("support-loaded branch product
law = 0") that adds no verification content; it gives a false impression of an
extra independent confirmation. It does not threaten correctness (the real
support-blindness content is carried by A2/A3/A5/A6/A12), but it is dead
scaffolding.

**Required change:**
Either remove the "support-loaded branch product law" check (py:57 / wl:49), or
replace it with a check that exercises the product law non-tautologically — e.g.
construct `product_loaded` and `Mtr` independently and assert
`R_target_loaded - product_loaded/Mtr == 0` only if `R_target_loaded` is instead
defined by a different route (e.g. from `T2_loaded`). Simplest correct action:
delete the redundant check in both engines.

**Verification:**
The verifier confirms the round-trip check is gone (or rewritten to use an
independently-defined `R_target_loaded`); both scripts still exit 0.

## Independent-derivation check (Mathematica)

Not independent — see F1. The `.wl` mirrors the `.py` choreography, symbol names,
assertion order, and hand-written closed-form targets. The single non-mirrored
construct is the support-loaded mass (`mMix*sSupport` vs `Mmix+Msupp`); everything
else (the `s`-perturbation derivations, the `D[..]*coeff` drift constructions, the
literal `eps1Expected`/`theta1Expected` targets) is a transcription.

## Engine cross-check

Both transcripts agree at the claimed level. Every `expectZero` lands at 0 in both
(`split-blocking drift eps_1 = 0`, `Xi_1 derived from T^2 matches defect law = 0`,
`tracking-factor drift = 0`, support-blindness lines = 0). The emitted closed forms
agree up to formatting:
- SymPy `Theta_1 = (-chi0*deltaU_1*(chi0 + 1) - chi_1*deltaU*(deltaU + 1))/((chi0 + 1)*(deltaU + 1)*(chi0 + deltaU + 1))`
- Mathematica `Theta_1 = -((chi1*deltaU*(1 + deltaU) + chi0*(1 + chi0)*deltaU1)/((1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU)))`
are identical. The support-rigid residuals and `Xi_1`/`R_1` forms likewise agree.
The spoiled-drift residuals differ only in surface form (both nonzero rational
functions) — expected, since `FullSimplify`/`simplify` pick different
representations. No `engine_disagreement`.

## Verdict justification

The scripts faithfully verify every Stage 181 paper deliverable (D1–D7) with
substantive, non-tautological, well-anchored checks, including a discriminating
spoiler counter-check for support-blindness and a support-rigid counter-check for
tracking non-sufficiency. The math holds up: I attacked the support-blindness
checks (could `∂_ζ ln T² = 0` be trivially true? No — `T2_loaded` genuinely carries
`ζ` through `product_loaded` and `Mtr`, and the spoiled packet confirms the check
fails when blindness is broken) and the drift checks (are the perturbation
derivatives identically zero? No — `eps`, `R_target`, `T²` all depend on the
perturbed symbols). The two findings are: (F1) the Mathematica script is a
transliteration of the SymPy one, weakening the cross-engine guarantee — medium;
and (F2) one redundant round-trip "product law" check that cannot fail — low.
Neither is a `paper_misalignment`; the value reconciliation (below) is complete
with 0 misaligned. Verdict: `findings` (script-side only), no stop-cold.

## Value Reconciliation (pass-2 augmentation)

Outputs are FRESH (sympy `.txt` 02:06 > `.py` 01:15; mathematica `.txt` 02:07 >
`.wl` 02:02), so reconciliation uses the committed transcripts directly.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| T² identity `Z_W(1+χ₀)²/(Ω_W²(1-ε)²) = (27π²Gc_s⁵/20a⁵c⁵)(1-ε_η)/R_target` | py:44–45 / wl:36–37; out line 5 | `.tex` appendix eq:app-part05-coherent-T (part05:722–729); notes eqs L54–62/120–127 | MATCH |
| ε split `ε_W(1-(2/11)δ_U/(1+δ_U))` | py:41 / wl:33 | `.tex` eq:app-part05-epsilon-split (part05:731–736); notes L108–109 | MATCH |
| ε₁ drift `(1-(2/11)δ_U/(1+δ_U))ε_W - (2ε_W/(11(1+δ_U)²))δ_{U,1}` | py:77–80 / wl:75–78; out "split-blocking drift eps_1 = 0" | `.tex` eq:app-part05-epsilon1-coherent (part05:749–756); notes L74–80/238–246 | MATCH |
| Ξ₁ defect law `ζ_Z-ω_W+2χ₁/(1+χ₀)+2ε₁/(1-ε)` | py:83–88 / out line 28 (wl out) | `.tex` eq:app-part05-Xi1-coherent-physical (part05:738–747); notes L66–69/224–229 | MATCH |
| support-blindness `∂_ζ ln T² = ∂_ζ ln R_target = ∂_ζ Ξ₁ = 0` | py:58,59,123 / wl:50,51,109; out lines 9,10,22 | `.tex` eq:app-part05-support-blindness (part05:758–765); notes L84–92/163–190 | MATCH |
| R₁ selected-branch `ω_W-η₁/(1-ε_η)-ζ_Z-2χ₁/(1+χ₀)-2ε₁/(1-ε)` | py:90–96 / out line 29 | notes §4 L299–308 (boxed); `.tex` card body implicit via Output | MATCH (notes) |
| Θ₁ tracking drift `-(χ₀(1+χ₀)δ_{U,1}+δ_U(1+δ_U)χ₁)/((1+χ₀)(1+δ_U)(1+χ₀+δ_U))` | py:129–132 / wl:114–118; out line 40 | `.tex` eq:app-part05-Theta-Sigma-tr (part05:807–814, σ-form); notes §5 L331–338 (boxed, χ₁/δ_{U,1} form) | MATCH |
| support-rigid Ξ₁ residual (χ₁=δ_{U,1}=0, nonzero) | py:138/143 / wl:123/128; out line 45 | notes §5 L351–368 (example: Θ₁=0 but Ξ₁≠0) | MATCH |

INTERNAL (scaffolding, no prose-carrier expected, raise no finding):
`Lam`/`lamNorm`, `Mmix`/`Msupp`/`Ssupport`/`product_loaded`/`R_target_loaded`/
`T2_loaded` (support-loaded reconstruction objects; the explicit `8/π²`,`(1-ε_η)`
mass prefactors are intermediate-only — notes carry only `S(ζ;ε)` and
`M_tr=M_mix S` symbolically, L150–153), `T2_selected`, the spoiled-drift residual,
`R1`/`Xi1` round-trip identity residuals, all pass/fail flags.

reconciliation: complete; 8 deliverable values checked, 0 misaligned

## Self-test notes

Checked (1) variable independence: every `diff`/`D` perturbation target genuinely
depends on the differentiation variable — `eps` depends on `epsW`,`deltaU`;
`R_target`/`T²` on all five drift symbols; `Rtr` on `chi0`,`deltaU` — so no
identically-zero-derivative trap, and the support-blindness `∂_ζ` checks act on
objects that genuinely carry `ζ` (confirmed by the wl:107 `FreeQ` guard and the
spoiler). (2) Trivial-case: the support-blindness zeros are non-trivial (the
spoiled packet at py:66/wl:64 gives a nonzero rational, proving discrimination);
the support-rigid `Xi_1` counter-check (py:143/wl:128) correctly remains nonzero.
(3) The proposed F1 re-derivation (series-extraction route for `ε₁`/`Θ₁`) does not
introduce a new constant — it reproduces the same `2/11` and `(1+δ_U)²` the paper
states, so no new paper_misalignment. F2 removal/rewrite touches only redundant
scaffolding. No UNFIXABLE/CRITICAL_DOWNSTREAM: the closed forms feed Stage 182's
slippage decomposition unchanged, and neither fix alters a derived constant or
sign.
