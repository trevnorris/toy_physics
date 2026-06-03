---
unit_id: 252
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T15:32:13Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.md
  paper_appendix: present
---

# Audit unit 252 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_252.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (row at line 102 references stage 252)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 252 is a channel-resolved energy compiler built on the Stage-251 microscopic export kernel `K_exp(s)=Γ₃s³+Γ₅s⁵`. It splits each coefficient into vacuum/lattice pieces (`Γ₃=Γ₃ᵛ+Γ₃ˡ`, `Γ₅=Γ₅ᵛ+Γ₅ˡ`, all ≥0) and, for an event window `[0,T]` with shape integrals `I₂=∫V̈²dt`, `I₃=∫V⃛²dt`, gives exact exported energies `E_vac=Γ₃ᵛI₂+Γ₅ᵛI₃`, `E_lat=Γ₃ˡI₂+Γ₅ˡI₃`, `E_exp=Γ₃I₂+Γ₅I₃`. The appendix row (line 102) summarizes the deliverables: "Exact vacuum/lattice exported energies, one-scalar event-shape quotient, single-growth cold-event specialization, safe-edge exported-energy theorem, and event-equivalent export rates." The distinct deliverables are: (1) the channel-resolved energy ledger; (2) the one-scalar quotient `r_V=I₃/I₂` with fractions `f_vac, f_lat` and `f_vac+f_lat=1`; (3) the speed-drift law `df_lat/dr_V=(Γ₅ˡΓ₃ᵛ−Γ₃ˡΓ₅ᵛ)/(Γ₃+Γ₅r_V)²` with endpoint limits `f_lat(0)=Γ₃ˡ/Γ₃`, `f_lat(∞)=Γ₅ˡ/Γ₅`; (4) the exponential-event specialization `I₂=Vin²s³(e^{2sT}−1)/2`, `I₃=s²I₂`, `r_V=s²`, with event-equivalent rates `γ_vac^eq=Γ₃ᵛs²+Γ₅ᵛs⁴` etc.; (5) the safe-edge theorem `E_exp,min^safe=(Vin²/2)(e²−1)μ_η(s₀²−s_c²)` and rate `γ_eff,safe^eq=μ_η(s₀²−s_c²)/s_c`; (6) the 3:1 split as the coefficient surface `f_lat(s_c)=3/4 ⟺ Γ₃ˡ+Γ₅ˡs_c²=3(Γ₃ᵛ+Γ₅ᵛs_c²)` plus the speed-independent family `f_lat=φ` when both channels share the split; (7) a Session-IV numeric benchmark calibration.

## What the script claims to verify

The docstring enumerates seven verified objects matching the paper's deliverables. Section 1 builds `E_vac/E_lat/E_tot` from the symbols and derives `f_vac, f_lat` from `Ev/Et`, asserting they reduce to the paper's closed forms and sum to 1. Section 2 differentiates `f_lat` w.r.t. `r_V` and takes its endpoint limits, comparing to the paper's drift law and limits. Section 3 integrates the actual derivatives of `V=Vin e^{st}` and checks `I₂`, `I₃=s²I₂`, and `E_vac/E_lat = γ^eq·I₁`. Section 4 substitutes the exponential integrals at the safe edge and checks the energy theorem and the safe-edge rate. Section 5 attempts the 3:1 coefficient surface and the φ-family reduction. Section 6 checks Session-IV benchmark numbers against hardcoded literals.

## Paper ↔ script cross-check

| Deliverable | Script coverage | Status |
|---|---|---|
| (1) channel-resolved energy ledger | L45-60: Ev/El/Et built, fractions derived from Ev/Et, asserted vs paper closed forms | match |
| (2) one-scalar quotient + f_vac+f_lat=1 | L48-49, L60 | match |
| (3) speed-drift law + endpoint limits | L66-77 | match |
| (4) exponential specialization + γ^eq rates | L89-107 | match |
| (5) safe-edge energy theorem + rate | L116-142; energy theorem fully carried (L139), but the rate-bridge identity `Γ₃s_c²+Γ₅s_c⁴ = μ_η(s₀²−s_c²)/s_c` is split into one good check (L140) and one tautology (L141) | partial |
| (6) 3:1 coefficient surface (the iff) | L148/L164 only re-expands an algebraic expression; the iff connecting `f_lat(s_c)=3/4` to the surface is never exercised | mismatch (tautological) |
| (6b) speed-independent φ-family | L150-165 | match |
| (7) Session-IV numeric benchmark | L171-206, all against hardcoded literals; partition 0.25/0.75 hardcoded, not derived from f_vac(s_c) on φ=3/4 | partial |
| Mathematica second engine | absent | missing |

`paper_alignment: partial` — the symbolic core (deliverables 1-4, 6b) is faithfully exercised; deliverable 6's iff is verified by a tautology, and the second engine is missing.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 58 | `simplify(fv − (G3v+G5v·rV)/(G3+G5·rV))==0` | (1)/(2) f_vac form | yes |
| A2 | sympy | 59 | `simplify(fl − (G3l+G5l·rV)/(G3+G5·rV))==0` | (2) f_lat form | yes |
| A3 | sympy | 60 | `simplify(fv+fl−1)==0` | (2) partition sums to 1 | yes |
| A4 | sympy | 75 | `simplify(dfl − dfl_expected)==0` | (3) drift law | yes |
| A5 | sympy | 76 | `simplify(fl_r0 − G3l/G3)==0` | (3) limit r_V→0 | yes |
| A6 | sympy | 77 | `simplify(fl_rinf − G5l/G5)==0` | (3) limit r_V→∞ | yes |
| A7 | sympy | 104 | `simplify(I2_exp − Vin²s³(e^{2sT}−1)/2)==0` | (4) I₂ form | yes |
| A8 | sympy | 105 | `simplify(I3_exp − s²·I2_exp)==0` | (4) I₃=s²I₂ | yes |
| A9 | sympy | 106 | `simplify(Ev.subs(I2,I3) − γ_v^eq·I1_exp)==0` | (4) E_vac=γ_vac^eq·I₁ | yes |
| A10 | sympy | 107 | `simplify(El.subs(I2,I3) − γ_l^eq·I1_exp)==0` | (4) E_lat=γ_lat^eq·I₁ | yes |
| A11 | sympy | 137 | `simplify(I1_safe − Vin²s_c(e²−1)/2)==0` | (5) I₁ at safe edge | yes |
| A12 | sympy | 138 | `simplify(Et_safe_expected − Et_safe)==0` | (5) total = γ_eff^eq·I₁ form | yes |
| A13 | sympy | 139 | `simplify(Et_safe_reduced − Vin²(e²−1)μ_η(s₀²−s_c²)/2)==0` | (5) safe-edge energy theorem | yes |
| A14 | sympy | 140 | `simplify(gamma_safe_eq_expected − gamma_eq.subs(s,s_c))==0` | (5) γ_eff,safe^eq = Γ₃s_c²+Γ₅s_c⁴ | yes |
| A15 | sympy | 141 | `simplify(gamma_safe_eq − μ_η(s₀²−s_c²)/s_c)==0` | (5) γ_eff,safe^eq = μ_η(...)/s_c | **no — X−X tautology** |
| A16 | sympy | 142 | `simplify(Ev_safe+El_safe − Et_safe_reduced)==0` | (5)/(2) channel split sums to total | partial (redundant w/ A3) |
| A17 | sympy | 164 | `simplify(split_surface) == G3l+G5l·s_c²−3·G3v−3·G5v·s_c²` | (6) 3:1 surface iff | **no — expand(X)==expand(X)** |
| A18 | sympy | 165 | `simplify(fl_phi − phi)==0` | (6b) φ-family | yes |
| A19 | sympy | 197-206 | `abs(num − literal) < tol` (×10) | (7) Session-IV benchmark | partial (hardcoded literals; 0.25/0.75 not derived) |

## Findings

### F1 — missing_verification_script (subtype missing_mathematica)

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/redteam/MANIFEST.yaml:8927-8929` (mathematica path null, exists false)

**What's wrong:**
Stage 252 is `is_checkpoint: false`, `is_status_only_candidate: false`. Per the dual-engine rule, a `.wl` is required wherever Mathematica CAN independently verify the stage. Every claim here is squarely within Mathematica's native algebra primitives: the rational-fraction partition (`Together`/`Cancel`/`Simplify`), the drift derivative (`D`), endpoint limits (`Limit`), the exponential shape integrals (`Integrate` of `(D[Vin Exp[s t],{t,2}])^2` over `{t,0,T}`), the safe-edge substitution (`/.` plus `Simplify`), the φ-family reduction, and the numeric benchmark (`N`). There is no obstruction that would make this single-engine-impossible — Mathematica can re-derive each result independently from the physical premises.

**Why this matters:**
A checkpoint-adjacent symbolic stage with no second engine has zero cross-engine corroboration. The SymPy script also contains a tautology (F2) and a re-expansion sham (F3) that a genuinely independent Mathematica derivation would not reproduce, so the second engine is exactly the catch the policy is designed for.

**Required change:**
Author a new independent Mathematica audit script (see directive F1) at `mathematica/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_mathematica_audit.wl`, deriving the claims via native Mathematica primitives and a different decomposition than the .py (NOT a transliteration).

**Verification:**
`redteam exec-mathematica 252` produces a `.wl` that exits 0 with `expectZero`/`expectApprox` checks covering claims M1–M9 in the directive.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py:122,141`

**What's wrong:**
L122 defines `gamma_safe_eq = sp.simplify(mu_eta * (s0**2 - sc**2) / sc)`. L141 then asserts `sp.simplify(gamma_safe_eq - mu_eta * (s0**2 - sc**2) / sc) == 0`. This is `X − X == 0`; it cannot fail regardless of the physics. The paper's deliverable (5) safe-edge rate is the *bridge* identity `Γ₃s_c²+Γ₅s_c⁴ = μ_η(s₀²−s_c²)/s_c` (card line 359-365). That bridge is the substance: `gamma_safe_eq_expected = safe_combo/s_c = Γ₃s_c²+Γ₅s_c⁴` (L121) should be shown equal to `μ_η(s₀²−s_c²)/s_c` *under the safe equality* `Γ₃s_c³+Γ₅s_c⁵=μ_η(s₀²−s_c²)` (defined at L114 but never used in an assertion). The script verifies the two halves separately (L140 gives expected == gamma_eq; L141 is a self-identity) but never connects them through the safe equality, so the rate-version of the safe-edge theorem is not actually exercised.

**Why this matters:**
The whole point of Section 4 is that the cold-survival edge *forces* the event-equivalent export rate to equal `μ_η(s₀²−s_c²)/s_c`. As written, the assertion that delivers this conclusion is a vacuous self-comparison; if the paper's bridge were wrong, the script would still pass.

**Required change:**
Replace the tautological L141 with a substitution-based check that carries the safe equality. Concretely (mirroring the energy-side pattern at L120/L139): assert that `gamma_safe_eq_expected` (= `safe_combo/sc`) with `safe_combo → mu_eta*(s0**2 - sc**2)` substituted equals `mu_eta*(s0**2 - sc**2)/sc`:
`assert sp.simplify(gamma_safe_eq_expected.subs(safe_combo, mu_eta*(s0**2 - sc**2)) - mu_eta*(s0**2 - sc**2)/sc) == 0`.
(If `simplify` has expanded `safe_combo` out of `gamma_safe_eq_expected`, substitute via `Gamma3*sc**3+Gamma5*sc**5 → mu_eta*(s0**2-sc**2)` on the unsimplified `safe_combo/sc` instead — see directive for the robust form.) Keep L140 as-is.

**Verification:**
L141 now references the safe equality (`mu_eta`, `s0`, AND `safe_combo`/Γ symbols all appear), so it can fail if the bridge is wrong; script still exits 0.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py:148,164`

**What's wrong:**
The paper's deliverable (6) is the *iff* (card lines 388-395, notes lines 388-396):
`f_lat(s_c)=3/4 ⟺ Γ₃ˡ+Γ₅ˡs_c² = 3(Γ₃ᵛ+Γ₅ᵛs_c²)`.
The script only computes `split_surface = simplify(expand((G3l+G5l*sc**2) - 3*(G3v+G5v*sc**2)))` (L148) and then asserts (L164) that `simplify(split_surface) == G3l + G5l*sc**2 - 3*G3v - 3*G5v*sc**2`. The RHS is just the expanded form of the LHS definition, so this is `expand(X) == expand(X)` — it verifies nothing about the partition fraction. The connection to `f_lat(s_c)` (the actual deliverable) is never made: the script never sets `f_lat(s_c) = 3/4` and shows it is equivalent to the surface vanishing.

**Why this matters:**
Deliverable (6) — the headline "the old 3:1 split becomes a microscopic coefficient surface" — is the new physics of the stage relative to Session IV. The script claims to verify it (docstring item 5) but actually only re-expands a polynomial. The iff could be false and the script would still pass.

**Required change:**
Replace the vacuous L164 with a check that exercises the iff. The clean route: form `f_lat(s_c) − 3/4` over a common denominator and show its numerator equals the surface expression (up to a positive factor). Concretely:
`surface = G3l + G5l*sc**2 - 3*(G3v + G5v*sc**2)`
`resid = sp.together(fl_sc - sp.Rational(3,4))`  (fl_sc is already defined at L124 as f_lat(s_c))
`num = sp.numer(resid)`
`assert sp.simplify(4*num - surface) == 0`  (the factor 4 clears the 3/4; verify the exact factor by self-test)
This makes the assertion fail unless `f_lat(s_c)=3/4` is genuinely equivalent to the surface vanishing. See directive for the self-tested factor.

**Verification:**
New assertion ties `fl_sc` (the partition fraction) to the surface; removing the `−3*(...)` term would break it. Script exits 0.

### F4 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py:174-175,180-181,184-185,204-206`

**What's wrong:**
Section 6 hardcodes the partition fractions `frac_v_num = 0.25`, `frac_l_num = 0.75` (L174-175) and then both the rates (`gamma_v_eq_num = 0.25*gamma_eq_safe_num`, L180) and the channel energies (`E_v_num = 0.25*E_diss_num`, L184) are computed by multiplying by these literals. The benchmark closes (L204: `E_v_num ≈ 0.00258365`) but only by confirming `0.25 × 0.01033460 = 0.00258365` — i.e. literal arithmetic against a literal. The paper (card §5, §5.1; lines 458-465, 482-492) obtains 0.25/0.75 specifically from the speed-independent φ=3/4 family, i.e. `f_vac(s_c)=1/4` on that family. The script never evaluates the symbolic `f_vac(s_c)`/`f_lat(s_c)` on the φ=3/4 family to confirm they equal 0.25/0.75; it just asserts the products of hardcoded fractions.

**Why this matters:**
This is the calibration step that ties the symbolic compiler to the Session-IV numbers. As written it never touches the symbolic partition — it is decoupled arithmetic on benchmark constants. The benchmark is a "calibration consistency check, not a theorem" (notes line 496), so severity is low, but the check should at least connect the hardcoded 0.25/0.75 to the φ=3/4 family it claims to instantiate.

**Required change:**
Before the numeric block, derive the fraction from the symbolic φ-family: confirm `f_vac` evaluated on `{G3l:phi*G3T, G3v:(1-phi)*G3T, G5l:phi*G5T, G5v:(1-phi)*G5T}` equals `1-phi`, then set `frac_v_num = float((1 - sp.Rational(3,4)))` and `frac_l_num = float(sp.Rational(3,4))` so the literals are derived from φ=3/4 rather than typed in. (Do not change the paper-stated benchmark inputs `t_cross_num`, `s0_num`, `E_diss_num`, or the literal targets in the asserts — those are anchored to card §5.) See directive.

**Verification:**
`frac_v_num`/`frac_l_num` are now derived from `f_vac` on φ=3/4; the numeric asserts (L197-206) still pass with the same literals.

## Independent-derivation check (Mathematica)

No `.wl` exists for stage 252, so transliteration is not applicable. F1 prescribes a NEW independent-route `.wl` and the directive's claim manifest specifies a *different decomposition* (e.g. derive fractions via `Together`/`Cancel` on `E_lat/E_exp` rather than re-substituting `I3 → rV*I2`; obtain the drift law via `D` of the simplified `f_lat` then `Together`; integrate the exponential shape integrals directly with `Integrate`; verify the 3:1 surface via the iff route F3 prescribes, not a re-expansion).

## Engine cross-check

Only one engine present; not applicable. `engines_agree: n/a`.

## Verdict justification

The symbolic core of stage 252 (deliverables 1-4 and the φ-family of 6) holds up under attack: the partition fractions are genuinely *derived* from the energy ledger and checked against the paper's closed forms; the drift derivative is taken w.r.t. a variable `f_lat` truly depends on (no variable-independence trap — `f_lat` contains `r_V`); the endpoint limits, the exponential shape integrals (`I₂=Vin²s³(e^{2sT}−1)/2`, factor verified by hand, `I₃=s²I₂`), and the energy=rate×I₁ relations are all non-tautological. The safe-edge *energy* theorem (A13) correctly carries the safe equality through `safe_combo`. The verdict is `findings`, not `clean`, because: (F1) the required second engine is absent and Mathematica can independently verify every claim; (F2) the safe-edge *rate* identity is asserted as a vacuous `X−X`, leaving the rate-version of the theorem unexercised; (F3) the headline 3:1-coefficient-surface iff is "verified" by re-expanding a polynomial (`expand(X)==expand(X)`), never touching `f_lat(s_c)`; and (F4) the benchmark partition is hardcoded 0.25/0.75 rather than derived from the φ=3/4 family. No `paper_misalignment` (the paper and script agree on all values; the issues are weak/vacuous *checks*, not wrong *targets*). No stop-cold: all findings are mechanically fixable and none changes a downstream-consumed constant (stage 253 consumes `γ_lat^eq(s_c)=Γ₃ˡs_c²+Γ₅ˡs_c⁴`, which is correctly defined here and untouched by the fixes).

## Self-test notes

Variable-independence trap: checked F2/F3 fixes — F2's substitution references `mu_eta, s0, sc` AND the Γ symbols (via safe_combo), so it is not a same-variable derivative or self-identity; F3 ties `fl_sc` (which depends on `G3l,G5l,G3v,G5v,sc`) to `surface`, not a w.r.t.-derivative, so no identically-zero trap. Trivial-case pre-check: for F3, with `G3l=3*G3v, G5l=3*G5v` (the surface = 0 case) `f_lat(s_c)` reduces to `(3G3v+3G5v sc²)/(4G3v+4G5v sc²)=3/4`, confirming the iff direction and that `4*numer(fl_sc−3/4)` equals `surface` up to sign/factor (factor to be pinned by Codex's run). Symmetry/parity: only definite integrals over `[0,T]` of strictly positive `exp` integrands appear — no symmetric-domain cancellation risk. Paper round-trip: F2/F3/F4 fixes introduce no new constants and leave all paper-stated literals (`t_cross`, `s0`, `E_diss`, `0.25/0.75`, the safe-edge prefactor) intact, so no new paper_misalignment.
