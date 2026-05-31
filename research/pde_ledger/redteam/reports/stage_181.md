---
unit_id: 181
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
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
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage181_coherent_tracking_defect.md` (present; internally titled "Stage 249" but is the source for stage 181)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 93, 721–766, 817 cover this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.txt`

## What the paper claims

Stage 181 inserts the coherent local D/N tracking branch into the one-port continuum transfer-shape law and proves a "support-blindness" theorem. The card's `Output` states verbatim: *"On the coherent local D/N tracking branch, \(\Xi_1\) is carried only by mixed/outgoing placement variables; coherent support ratio \(\zeta\) drops out."* The appendix (lines 721–766) and notes enumerate the concrete deliverables: (1) the coherent-branch transfer-shape identity \(\mathcal T^2 = Z_W(1+\chi_0)^2/(\Omega_W^2(1-\epsilon)^2) = (27\pi^2 G c_s^5/20a^5c^5)(1-\epsilon_\eta)/R_{\rm target}\) (eq:app-part05-coherent-T); (2) the split-blocking ratio \(\epsilon = \epsilon_W(1-\tfrac{2}{11}\delta_U/(1+\delta_U))\) (eq:epsilon-split); (3) the weak-axisymmetric defect law \(\Xi_1 = \zeta_Z-\omega_W+2\chi_1/(1+\chi_0)+2\epsilon_1/(1-\epsilon)\) (eq:Xi1-coherent-physical), itself *defined* in the notes (§3, eqs 224–229) as \(\Xi_1 = \delta\ln\mathcal T_A^2/(\epsilon\lambda_A)\); (4) the \(\epsilon_1\) drift formula (eq:epsilon1-coherent); (5) the support-blindness triple \(\partial_\zeta\ln\mathcal T^2 = \partial_\zeta\ln R_{\rm target} = \partial_\zeta\Xi_1 = 0\) (eq:support-blindness — the headline). The notes additionally give (6) the selected-branch identity \(\Xi_1 = -\eta_1/(1-\epsilon_\eta) - \mathcal R_1\) (§4) and (7) the tracking-factor drift \(\Theta_1\) (§5), with the physical conclusion that exact tracking rigidity \(\Theta_1=0\) does not force \(\Xi_1=0\).

## What the script claims to verify

The docstring lists six checks: the coherent-branch transfer-shape identity; support-blindness (d/dζ of T² and R_target vanish); the split-blocking drift formula; the coherent-branch defect law \(\Xi_1\); the selected-branch identity; and the tracking-factor drift. In practice the load-bearing, genuinely-failable assertions are: the `direct-selected transfer-shape identity` (line 54), the four `support-loaded` reconstructions/derivatives (lines 55–59) plus a `spoiled` adversarial control (lines 61–67), the `split-blocking drift eps_1` (line 80), the `selected-branch identity` (line 97), and the `tracking-factor drift` (line 110). \(\Xi_1\) and \(\mathcal R_1\) are *defined* by typing the paper formulas (lines 82–95); the support-rigid specialization (lines 115–121) is a print + nontriviality guard. The Mathematica script mirrors all of this except it omits the `spoiled` control.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Verdict |
|---|---|---|
| (1) Transfer-shape identity (eq:coherent-T) | line 54 `T2_direct - T2_selected == 0` (hinges on Lam carrying Ω_W²) | match |
| (2) ε split-blocking (eq:epsilon-split) | line 41 definition; exercised via eps1 | match |
| (3) Ξ₁ defect law (eq:Xi1-coherent-physical) | line 82 — **definition only**, never derived from δ ln T²/(ελ) | **partial / missing** |
| (4) ε₁ drift (eq:epsilon1-coherent) | line 80 `eps1 - eps1_expected == 0` (SymPy diff vs paper closed form) | match |
| (5a) ∂_ζ ln T² = 0 | line 58 (support-loaded route — genuine ζ-cancellation) | match |
| (5b) ∂_ζ ln R_target = 0 | line 59 (support-loaded route) + line 55 reconstruction | match |
| (5c) **∂_ζ Ξ₁ = 0** (headline Output) | **no check** — Ξ₁ defined ζ-free at line 82, so trivially true and not even tested | **missing** |
| (6) Selected-branch identity (notes §4) | line 97 — **tautological** (R1 defined as −Xi1 − eta1/(1−eps_eta)) | mismatch (vacuous) |
| (7) Tracking-factor drift Θ₁ (notes §5) | line 110 `Theta1 - Theta1_expected == 0` (SymPy log-diff vs paper) | match |

`paper_alignment: partial` — the script and paper never disagree on a value or identity (no value/target mismatch), but two headline deliverables (the Ξ₁ defect law derivation and ∂_ζ Ξ₁ = 0) are under-verified, and the selected-branch identity is checked tautologically.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 54 | `simplify(T2_direct - T2_selected) == 0` | claim 1 | yes |
| A2 | sympy | 55 | `R_target_loaded - R_target == 0` | claim 5b | yes |
| A3 | sympy | 56 | `T2_loaded - T2_direct == 0` | claim 1/5a (support route) | yes |
| A4 | sympy | 57 | `R_target_loaded*Mtr - product_loaded == 0` | claim 5b (route consistency) | yes |
| A5 | sympy | 58 | `diff(log(T2_loaded), zeta) == 0` | claim 5a | yes |
| A6 | sympy | 59 | `diff(log(R_target_loaded), zeta) == 0` | claim 5b | yes |
| A7 | sympy | 66 | `dlnR_spoiled != 0` (adversarial control) | claim 5b (negative control) | yes |
| A8 | sympy | 80 | `eps1 - eps1_expected == 0` | claim 4 | yes |
| A9 | sympy | 97 | `Xi1 + eta1/(1-eps_eta) + R1 == 0` | claim 6 | **no — tautological** |
| A10 | sympy | 110 | `Theta1 - Theta1_expected == 0` | claim 7 | yes |
| A11 | sympy | 120 | `Xi_support_rigid != 0` (nontriviality) | claim 7 consequence | yes (guard) |
| A12 | math | 46–51 | `expectZero[...]` ×6 | claims 1, 5a, 5b | yes (but transliterated) |
| A13 | math | 62 | `expectZero["split-blocking drift eps_1"]` | claim 4 | yes |
| A14 | math | 69 | `expectZero["selected-branch identity"]` | claim 6 | **no — tautological** |
| A15 | math | 82 | `expectZero["tracking-factor drift"]` | claim 7 | yes |
| — | — | — | **(no assertion)** | **claim 3 (Ξ₁ law), claim 5c (∂_ζ Ξ₁)** | **missing** |

## Findings

### F1 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py:82-87`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl:64`

**What's wrong:**
The paper's central deliverable is the coherent-branch defect law (eq:app-part05-Xi1-coherent-physical), and the notes (§3, eqs 224–229) are explicit that it is *derived*: \(\Xi_1 = \delta\ln\mathcal T_A^2/(\epsilon\lambda_A)\), obtained by perturbing the transfer shape \(\mathcal T^2 = Z_W(1+\chi_0)^2/(\Omega_W^2(1-\epsilon)^2)\) in the five primitive drifts \((\zeta_Z,\omega_W,\chi_1,\varepsilon_W,\delta_{U,1})\). The script never performs this perturbation. Instead it *types the answer in as a definition*:
```python
Xi1 = sp.simplify(zetaZ - omegaW + 2 * chi1 / (1 + chi0) + 2 * eps1 / (1 - eps))
```
There is no `expect_zero` comparing this `Xi1` to the log-derivative of `T2_direct`. The script verifies `T2_direct` (A1) and verifies `eps1` (A8), but never closes the loop that \(\Xi_1 = \delta\ln\mathcal T^2/(\epsilon\lambda)\). Consequently the headline `Output` claim "coherent support ratio \(\zeta\) drops out" — formalized as \(\partial_\zeta\Xi_1 = 0\) (eq:support-blindness, third equality) — is verified only for \(\mathcal T^2\) and \(R_{\rm target}\) (A5, A6), and is true for `Xi1` *only because line 82 was written in ζ-free variables*, i.e. by construction. The single most important conclusion of the stage is asserted, not proven.

**Why this matters:**
A reviewer reading the transcript sees `Xi_1 = ...` printed and the support-blindness conclusion stated, and would reasonably believe the defect law was checked. It was not. If the paper's \(\Xi_1\) formula (or the sign of any of its four terms) were wrong, this script would still pass, because `Xi1` is whatever line 82 says. The derivation from \(\mathcal T^2\) is exactly the load-bearing algebra of the stage.

**Required change:**
Add an `expect_zero` that derives \(\Xi_1\) from the transfer shape and compares it to the line-82 formula, and add the missing \(\partial_\zeta\Xi_1 = 0\) check via the support-loaded route. Concretely (SymPy), introduce a perturbation parameter `s` (standing for \(\epsilon\lambda\)) and perturb the primitive variables of `T2_direct`:
```python
s = sp.symbols("s", real=True)
T2_pert = (ZW * (1 + s * zetaZ)) * (1 + (chi0 + s * chi1))**2 \
          / ((OmegaW2 * (1 + s * omegaW)) * (1 - epsW_of(epsW + s * epsW1, deltaU + s * deltaU1))**2)
```
where `epsW_of(eW, dU) = eW * (1 - sp.Rational(2,11) * dU / (1 + dU))` is the split-blocking map (matching line 41). Then
```python
Xi1_derived = sp.simplify(sp.diff(sp.log(T2_pert), s).subs(s, 0))
expect_zero("Xi_1 derived from T^2 matches defect law", Xi1_derived - Xi1)
```
For the missing ζ-blindness of the defect itself, derive \(\Xi_1\) from the ζ-bearing `T2_loaded` (which line 56 already proves equals `T2_direct`) using the *same* primitive perturbation, then assert it carries no ζ:
```python
# perturb the same primitives inside the support-loaded route, call the result T2_loaded_pert
Xi1_loaded = sp.simplify(sp.diff(sp.log(T2_loaded_pert), s).subs(s, 0))
expect_zero("d/dzeta Xi_1 (support-loaded route)", sp.diff(Xi1_loaded, zeta))
```
Mirror both new checks in the Mathematica script (`D[Log[t2Pert], s] /. s -> 0`, etc.). See directive for the exact perturbation map and self-test confirmation.

**Verification:**
New transcript lines `Xi_1 derived from T^2 matches defect law = 0` and `d/dzeta Xi_1 (support-loaded route) = 0` appear in both engines' outputs; scripts exit 0.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py:89-97`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl:65-69`

**What's wrong:**
The "selected-branch identity" check (line 97) is `expect_zero("selected-branch identity", Xi1 + eta1/(1 - eps_eta) + R1)`. But `R1` is *defined* (lines 89–95) as
```python
R1 = omegaW - eta1/(1 - eps_eta) - zetaZ - 2*chi1/(1 + chi0) - 2*eps1/(1 - eps)
```
i.e. exactly `-(Xi1) - eta1/(1 - eps_eta)` written out term-by-term (since `Xi1 = zetaZ - omegaW + 2*chi1/(1+chi0) + 2*eps1/(1-eps)`). Therefore `Xi1 + eta1/(1-eps_eta) + R1` is identically zero by construction — the assertion cannot fail for any physics. It verifies nothing about the paper's selected-branch identity (notes §4, eq:R1-coherent), which is the *claim* that the directly-computed \(\mathcal R_1 := \delta\ln R_{{\rm target},A}/(\epsilon\lambda_A)\) equals \(\omega_W - \eta_1/(1-\epsilon_\eta) - \zeta_Z - 2\chi_1/(1+\chi_0) - 2\epsilon_1/(1-\epsilon)\).

**Why this matters:**
This is the only check that the direct-port and selected-branch descriptions "agree exactly on the coherent tracking branch" (notes §4). As written, the agreement is imposed by definition, so a sign or factor error in the \(R_{\rm target}\) drift would never be caught here.

**Required change:**
Compute `R1` independently as the log-derivative of `R_target` w.r.t. the perturbation parameter `s` (using the same primitive perturbation map as F1, plus the \(\eta\) drift \(\delta\epsilon_\eta = \epsilon\lambda\,\eta_1\)), then compare to the paper's closed form. Concretely:
```python
R_target_pert = (Lam_pert) * (1 - (eps_eta + s*eta1)) * (1 - eps_pert)**2 \
                / ((ZW*(1 + s*zetaZ)) * (1 + (chi0 + s*chi1))**2)
R1_derived = sp.simplify(sp.diff(sp.log(R_target_pert), s).subs(s, 0))
expect_zero("R_1 derived from R_target matches closed form", R1_derived - R1)
```
where `Lam_pert` perturbs `OmegaW2 -> OmegaW2*(1 + s*omegaW)` inside `Lam` and `eps_pert` is the perturbed split-blocking ratio. Keep `R1` as the paper closed form (line 89), but make it a *verified* target. Then the existing line-97 identity becomes a genuine consistency check between the two derived drifts. Mirror in Mathematica.

**Verification:**
New line `R_1 derived from R_target matches closed form = 0` in both transcripts; line-97 `selected-branch identity = 0` retained but now non-vacuous because `R1` is independently anchored; scripts exit 0.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl:32-89`
- (compared against) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py:40-118`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py` rather than an independent re-derivation. Corresponding sections:

SymPy 44–52:
```python
T2_direct = ... ; T2_selected = ... ; Mmix = ... ; Msupp = ... ;
Ssupport = ... ; Mtr = Mmix + Msupp ; product_loaded = ... ;
R_target_loaded = product_loaded / Mtr ; T2_loaded = (Lam/OmegaW2)*(1-eps_eta)/R_target_loaded
```
Mathematica 36–44:
```wl
t2Direct = ... ; t2Selected = ... ; mMix = ... ; mSupp = ... ;
sSupport = ... ; mTr = mMix + mSupp ; productLoaded = ... ;
rTargetLoaded = productLoaded/mTr ; t2Loaded = (lamNorm/omegaW2)*(1-epsEta)/rTargetLoaded
```
Identical intermediate-variable choreography, same names transliterated, same order. Likewise `eps1 = D[eps,epsW]*epsW1 + D[eps,deltaU]*deltaU1` (math 57) is the exact mirror of SymPy line 75, and the `Theta1` log-derivative (math 76) mirrors SymPy 105. The same six `expectZero` assertions appear in the same sequence. The two engines do not constitute independent derivations; the Mathematica run re-runs the SymPy algebra in WL syntax. (It is also not even a complete port: the SymPy `spoiled` adversarial control at lines 61–67 has no Mathematica counterpart.)

**Why this matters:**
The second-engine policy exists to catch algebra errors that one CAS would make consistently (e.g. a mis-simplification both engines inherit if they follow identical steps). A transliteration provides little of that protection.

**Required change:**
Restructure the Mathematica derivation so it reaches the same load-bearing results by an independent route, rather than echoing the SymPy intermediate-variable sequence. Minimal acceptable independence for this unit: (a) derive the support-blindness of \(\mathcal T^2\) by directly showing \(\mathcal T^2\) (as the closed form `zW*(1+chi0)^2/(omegaW2*(1-eps)^2)`) contains no `zeta` and that the support-loaded reconstruction returns to it — but compute the loaded route by a different grouping (e.g. solve `rTargetLoaded == rTarget` for the cancellation rather than constructing `mTr` term-by-term); and (b) add the F1/F2 derivation checks (`D[Log[t2Pert], s]`, `D[Log[rTargetPert], s]`) which are independent of the SymPy intermediate choreography. Also port the `spoiled` negative control. The goal is that the `.wl` no longer mirrors the `.py` line for line.

**Verification:**
The `.wl` no longer reproduces the SymPy `mMix/mSupp/mTr/productLoaded` choreography verbatim; it includes the new F1/F2 derivation checks and a negative control; output still shows all `PASS` and exits 0.

## Independent-derivation check (Mathematica)

Not independent. The `.wl` (lines 32–89) is a direct transliteration of the SymPy script's variable choreography and assertion sequence (see F3 for quoted correspondences). It also omits the SymPy `spoiled` adversarial control (sympy 61–67). This is the basis for finding F3.

## Engine cross-check

Both engines agree at the level they claim. The printed final forms match after normalization:
- SymPy `Xi_1` (output line 26) and Mathematica `Xi_1` (output line 33) are the same rational function (denominators `11*(deltaU+1)^2*(epsW*(9deltaU+11)/(11(deltaU+1)) - 1)` vs `(1+deltaU)*(-11(1+deltaU)+(11+9deltaU)epsW)` are equal up to the common factor; numerators agree).
- `Theta_1` (sympy output 33 vs math output 41) agree.
- All shared `expect_zero`/`expectZero` checks print `0` / `PASS` in both transcripts; both exit 0.

No `engine_disagreement` finding. Both engines share the same blind spots (the missing Ξ₁ derivation and the tautological selected-branch identity), which is expected given the transliteration.

## Verdict justification

`verdict: findings`, `stop_cold: null`. The arithmetic the script *does* exercise holds up under attack: the transfer-shape identity (A1) genuinely hinges on Lam carrying the Ω_W² factor; the support-loaded ζ-cancellation (A2–A6) is real algebra (I verified by hand that Ssupport cancels in R_target_loaded → R_target), and the `spoiled` control (A7) confirms the check can fail; the ε₁ drift (A8) and Θ₁ drift (A10) are genuine SymPy-diff-vs-paper-closed-form comparisons. But three deliverables fall short of the paper: the headline Ξ₁ defect law is typed in as a definition rather than derived from δ ln T²/(ελ) (F1, high), the support-blindness of Ξ₁ itself (the `Output` line "ζ drops out") is consequently true only by construction and never tested, the selected-branch identity is tautological because R1 is defined as −Xi1−η₁/(1−ε_η) (F2, medium), and the Mathematica engine is a transliteration of the SymPy one (F3, medium). None is a paper↔script value/target disagreement, so no `paper_misalignment` and no user resolution is needed; all three are fixable by Codex adding genuine derivation checks. No downstream constant changes, so not `CRITICAL_DOWNSTREAM`. Minor note: both scripts' banner mislabels the stage as "STAGE 164" (sympy line 33, math line 26) though the Mathematica closing print correctly says "Stage 181" (line 101) — cosmetic, folded into the directive.

## Self-test notes

I checked variable independence for the proposed F1/F2 perturbation: `T2_pert` depends on all of ZW, chi0, OmegaW2, epsW, deltaU, so `diff(log(T2_pert), s)|_{s=0}` is nonzero and depends on every one of (ζ_Z, ω_W, χ₁, ε_W, δ_U,1) — no derivative is identically zero, and the result matches the line-82 formula (hand-verified: δ ln T²/(ελ) = ζ_Z − ω_W + 2χ₁/(1+χ₀) + 2ε₁/(1−ε)). No integrals/parity traps in this unit. Trivial-case check: at s=0 the perturbed forms reduce to the unperturbed closed forms, so the `expect_zero(... - Xi1)` residual collapses to 0 as required. Paper round-trip: the proposed derivation uses only the constants already in the paper (split-blocking 2/11, the Lam normalization), introduces no new constant, so it cannot create a new paper_misalignment.
