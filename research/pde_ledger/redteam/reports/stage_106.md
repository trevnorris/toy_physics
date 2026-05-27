---
unit_id: 106
batch: IV.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md"]
  paper_appendix: present
---

# Audit unit 106 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_106.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read line 1246 plus all stage-106-relevant rows 17, 19, 220-295, 300-345, 1175-1200)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.txt`

## What the paper claims

The card defines Stage 106 as the "reduced 2.5PN closure on the canonical outgoing DtN branch". Its `Derivation ledger` says the computation "isolates the reduced product `\widehat m_0^{\,2}\chi_Q N_Q=1` and the canonical condition `\chi_Q=1`". Its `Checks` block (verbatim) lists THREE distinct verification items: (i) the product `\widehat m_0^{\,2}\chi_Q N_Q` keeps source, conservative, and outgoing factors separate; (ii) higher odd terms begin beyond the point-particle 2.5PN coefficient; (iii) the outgoing `l=2` DtN fingerprint matches the normalized `z=\omega a/c_s` expansion. The notes additionally specify the canonical invariant coefficients `\overline K_0=54Gc_s^5/(5a^5c^5)`, `\overline K_2=6Gc_s^3/(5a^3c^5)`, `\overline K_4=8Gc_s/(15ac^5)`, `\overline\Gamma_5=2G/(5c^5)`, and the conditional outcome `\gamma_{\rm quad}^{\rm eff}=\widehat m_0^{\,2}\Gamma_5=2G/(5c^5)` on the canonical branch. Appendix eq. `app-part04-Lambda-out-dtn` (line 308) gives the explicit DtN expansion `-3 + z^2/3 + z^4/9 + i z^5/9 + O(z^6)`, and eq. `app-part04-Yout-dtn` (line 317) gives the normalized response `1 + z^2/9 + 4z^4/81 + i z^5/27 + O(z^6)`, from which the appendix derives `\chi_Q=1`.

## What the script claims to verify

The script (banner: "STAGE 89 — REDUCED 2.5PN CLOSURE ON CANONICAL OUTGOING DtN BRANCH" — note the stale stage label) solves the algebraic constraint `m0hat^2 * chi_Q * N_Q == 1` for `N_Q` to get `N_Q = 1/(chi_Q*m0hat^2)`, substitutes `chi_Q=1` to get `N_Q = 1/m0hat^2`, then substitutes `m0hat=1` to check `N_Q=1`. It hardcodes the four target invariants `K0_target = 54*G*c_s**5/(5*a**5*c**5)`, `K2_target = 6*G*c_s**3/(5*a**3*c**5)`, `K4_target = 8*G*c_s/(15*a*c**5)`, `Gamma5_target = 2*G/(5*c**5)`, then forms `K0 = NQ_general * K0_target` and derives `K2 = K0/(4*(3*c_s/(2*a))**2)`, `K4 = K0/(4*(3*c_s/(2*a))**4)` (note `K2_target` and `K4_target` are computed but never used as anchors). It then checks `K4 - 4*K2^2/K0 == 0` and `Gamma5 - 9*sqrt(K2^5/K0^3) == 0`, plus `gamma_eff_canonical - Gamma5_target == 0` where `gamma_eff_canonical = (m0hat^2 * Gamma5).subs(chi_Q,1)`. The Mathematica script is a line-by-line transliteration of the same algebra.

## Paper ↔ script cross-check

| Paper deliverable | Script-side coverage | Status |
|---|---|---|
| Reduced product `\widehat m_0^2\chi_Q N_Q=1` displays three independent factors (Check 1) | Constraint is declared and solved for N_Q; the three factors appear symbolically but no orthogonality/independence check is performed | partial |
| Higher odd terms begin beyond point-particle 2.5PN (Check 2) | No `\omega`-expansion is constructed; no check that the next nonzero odd contribution sits at `\omega^7` or higher | missing |
| Outgoing `l=2` DtN fingerprint matches `z=\omega a/c_s` expansion (Check 3) — *this is the derivation that fixes `\chi_Q=1`* | The script imports `\chi_Q=1` as a substitution; never constructs `\Lambda_2^{\rm out}(z)` or `\widehat Y_2^{\rm out}(z)`; never verifies the expansion `1+z^2/9+4z^4/81+i z^5/27+O(z^6)` | missing |
| `N_Q=1` in strict point-particle limit (notes line 35) | Verified at `m0hat=1, chi_Q=1` (script line 38) | match |
| Canonical invariants `\overline K_0=54Gc_s^5/(5a^5c^5)`, `\overline K_2=6Gc_s^3/(5a^3c^5)`, `\overline K_4=8Gc_s/(15ac^5)`, `\overline\Gamma_5=2G/(5c^5)` (notes lines 42-51) | Hardcoded literally as `K0_target` etc.; no derivation; `K2_target`, `K4_target` are defined but never anchored to a check (the script re-derives K2, K4 from K0 via `K0/(4Ω^2)^n`) | partial |
| `\gamma_{\rm quad}^{\rm eff}=\widehat m_0^2\Gamma_5=2G/(5c^5)` on canonical branch (notes line 60) | Verified algebraically at script line 64 | match |

Dominant pattern: paper-side Checks (ii) and (iii) are entirely absent from the script. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 38 | `expect_zero(NQ_canonical.subs(m0hat,1) - 1)` | point-particle `N_Q=1` (notes line 35) | yes (non-trivially solves the constraint, but tautologically true once the constraint is given) |
| A2 | sympy | 55 | `expect_zero(K4 - 4*K2**2/K0)` | none — algebraic identity true by construction (K2≡K0/(4Ω^2), K4≡K0/(4Ω^4)) | no (tautological) |
| A3 | sympy | 56-59 | `expect_zero(Gamma5 - 9*sqrt(K2^5/K0^3))` | weak proxy for K0_target/Γ5_target consistency | partial (passes only because hardcoded literals happen to be mutually consistent; does not verify either against an upstream derivation) |
| A4 | sympy | 61-65 | `expect_zero(gamma_eff_canonical - Gamma5_target)` | `\gamma_{\rm quad}^{\rm eff}=2G/(5c^5)` on canonical branch | partial (follows directly from N_Q=1/m0hat^2 substitution; equivalent to A1 up to a multiplicative target) |
| A5 | mathematica | 40-43 | `expectZero[(nqCanonical /. m0hat -> 1) - 1]` | mirror of A1 | yes (same content as A1) |
| A6 | mathematica | 60 | `expectZero[k4 - 4*k2^2/k0]` | mirror of A2 | no (tautological, mirrors A2) |
| A7 | mathematica | 61-64 | `expectZero[gamma5 - 9*Sqrt[k2^5/k0^3]]` | mirror of A3 | partial (mirrors A3) |
| A8 | mathematica | 66-67 | `expectZero[gammaEffCanonical - gamma5Target]` | mirror of A4 | partial (mirrors A4) |
| --- | --- | --- | --- | Check (ii) higher odd > 2.5PN | **no assertion** |
| --- | --- | --- | --- | Check (iii) DtN `z`-expansion fingerprint | **no assertion** |

## Findings

### F1 — paper_misalignment

**Severity:** high
**Subtype:** script_missing_paper_claim
**Files:**
- paper side: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_106.tex:21-25` (Checks block items ii and iii)
- paper side: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:300-327` (eq:app-part04-Lambda-out-dtn, eq:app-part04-Yout-dtn, eq:app-part04-chiQ-equals-one)
- script side: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:25-71` (no z-expansion code anywhere)
- script side: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:27-75` (no z-expansion code anywhere)

**What's wrong:**
The paper card lists three distinct Checks. The scripts cover only the first (and weakly). Check (ii) "higher odd terms begin beyond the point-particle 2.5PN coefficient" is never tested — there is no expansion in ω, no demonstration that the next-to-leading odd coefficient sits at ω^7 or beyond, and no quote of the `O(ω^6)` remainder in the reduced response. Check (iii) "outgoing l=2 DtN fingerprint against the normalized z=ωa/c_s expansion" is the very calculation that *derives* `χ_Q = 1` in the appendix (eqs. 308-327). The scripts never construct `Λ_2^out(z) = z d/dz ln h_2^(1)(z)` and never verify its expansion `-3 + z^2/3 + z^4/9 + i z^5/9 + O(z^6)` nor `Ŷ_2^out(z) = 1 + z^2/9 + 4z^4/81 + i z^5/27 + O(z^6)`. Instead, `χ_Q = 1` is imported as a substitution from "Stage 88" without script-side verification.

Paper quote (stage_106.tex:21-25):
> Check the product `\widehat m_0^{\,2}\chi_Q N_Q` keeps source, conservative, and outgoing factors separate. Check that higher odd terms begin beyond the point-particle 2.5PN coefficient. Check the outgoing l=2 DtN fingerprint against the normalized `z=\omega a/c_s` expansion.

Script behavior: the entire algebraic content is contained between sympy lines 30-65; nothing constructs a z-expansion of any DtN operator.

**Why this matters:**
A reader of the paper card reasonably expects all three Checks to have script-side coverage. Two of them are entirely absent. The most important missing check is the DtN expansion (iii), because it is the upstream derivation that fixes `χ_Q=1` (everything else in the script depends on `χ_Q=1` as a substitution). With Check (iii) missing, the script's apparent "closure" is purely algebraic: solving `m0hat^2·χ_Q·N_Q=1` for N_Q, then setting both m0hat=1 and χ_Q=1 to recover 1. That is not a verification of the canonical branch realization — it is a tautology dressed as one.

**Required change:**
See directive F1 — pending user resolution. The user must choose whether (a) the script is expected to verify Checks (ii) and (iii) directly and must be extended, or (b) Checks (ii) and (iii) are carried forward from upstream Stage 88 (DtN derivation) and the paper card should explicitly say so rather than list them as Stage-106-side checks.

**Verification:**
After user resolution: either (a) the scripts contain new asserts on the z-expansion of `Λ_2^out` and the `O(ω^7)` cutoff, or (b) the paper card's Checks block is updated to reference the upstream stage that owns those checks.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:46-47, 55`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:51-52, 60`

**What's wrong:**
The check `expect_zero("branch identity K4 - 4 K2^2 / K0", K4 - 4 * K2**2 / K0)` (sympy line 55) is algebraically tautological by the script's own definitions. The script defines `K2 = K0/(4*(3*c_s/(2*a))**2) = K0/(4Ω^2)` and `K4 = K0/(4*(3*c_s/(2*a))**4) = K0/(4Ω^4)`. Then `4*K2^2/K0 = 4*(K0/(4Ω^2))^2/K0 = K0^2/(4Ω^4 K0) = K0/(4Ω^4) ≡ K4`. The identity holds regardless of the value of K0, regardless of N_Q, regardless of the physics. It cannot fail no matter what numerical constants the script picks. The Mathematica mirror at line 60 has the same defect.

**Why this matters:**
This is presented in the script's print output as a verification of a "branch identity", but it verifies nothing about the branch — it only verifies that the script's own definitions are internally consistent at the level of polynomial algebra. A reader (or auditor) is misled into thinking the canonical even-fingerprint relation `K_0 K_4 = 4 K_2^2` (appendix line 114) has been independently confirmed; in fact the script enforces it by construction.

**Required change:**
Either (a) remove the assertion since it carries no information, or (b) re-derive K2 and K4 from the explicit per-coefficient targets `K2_target = 6Gc_s^3/(5a^3c^5)` and `K4_target = 8Gc_s/(15ac^5)` (which are computed at lines 41-42 but never referenced in any check), then assert `K4_target - 4*K2_target^2/K0_target == 0`. Option (b) actually exercises a non-trivial numerical-consistency claim about the four upstream constants.

**Verification:**
After fix: the assertion residual still simplifies to zero, but now it depends on the values of K0_target, K2_target, K4_target rather than on the script-defined K2=K0/(4Ω^2) etc. Substituting any one of those four numerical constants with a wrong value (e.g. K2_target → 7Gc_s^3/(5a^3c^5)) should cause the assertion to fail.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py` (entire)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl` (entire)

**What's wrong:**
The Mathematica script is a line-by-line transliteration of the SymPy script: identical variable names (`m0hat`, `chi_Q`/`chiQ`, `N_Q`/`NQ`, `K0`, `K2`, `K4`, `Gamma5`/`gamma5`), identical intermediate quantities defined in the identical order, identical RESULT prose at the bottom. Compare:

SymPy lines 30-37:
```
constraint = sp.Eq(m0hat**2 * chi_Q * N_Q, 1)
print("Normalization constraint =", constraint)
NQ_general = sp.simplify(sp.solve(constraint, N_Q)[0])
print("N_Q on the general outgoing branch =", NQ_general)
NQ_canonical = sp.simplify(NQ_general.subs(chi_Q, 1))
print("N_Q on the canonical outgoing branch =", NQ_canonical)
```

Mathematica lines 32-39:
```
constraint = m0hat^2*chiQ*NQ == 1;
Print["Normalization constraint = ", fmt[constraint]];
nqGeneral = FullSimplify[NQ /. First[Solve[constraint, NQ, Reals]], Assumptions -> $Assumptions];
Print["N_Q on the general outgoing branch = ", fmt[nqGeneral]];
nqCanonical = FullSimplify[nqGeneral /. chiQ -> 1, Assumptions -> $Assumptions];
Print["N_Q on the canonical outgoing branch = ", fmt[nqCanonical]];
```

SymPy lines 40-48 define `K0_target, K2_target, K4_target, Gamma5_target` then `K0 = NQ_general * K0_target` then `K2 = K0/(4*(3*c_s/(2*a))**2)`. Mathematica lines 45-53 do the exact same five definitions in the exact same order with the exact same intermediate-variable choreography. Both end with the same three assertions (`K4 - 4K2^2/K0`, `Gamma5 - 9*sqrt(K2^5/K0^3)`, `gamma_eff_canonical - target`) in the exact same sequence.

**Why this matters:**
The second-engine policy requires Mathematica to derive the result independently, so that any algebra bug in one engine has a chance to be caught by the other. Transliteration defeats this purpose: a bug in the SymPy algorithm (e.g., picking the wrong branch of Solve, or a sign error in N_Q solution) will be reproduced verbatim in Mathematica.

**Required change:**
Re-author the Mathematica script so that it derives the canonical-branch closure from a structurally different starting point. Suggested independent path: (a) start from the appendix's retarded one-pole form `\widehat Y_Q^{\rm ret}(\omega) = 3/4 + 1/4 · 1/(1 - \omega^2/\Omega_Q^2 - i\chi_Q\sigma_Q^{\rm can}\omega^5)` (eq. app-part04-retarded-one-pole) and series-expand in ω to extract `\overline\Gamma_5 = \chi_Q · 9\overline K_0/(32\Omega_Q^5)` (eq. app-part04-gamma5-chiN), then derive `\widehat m_0^2\chi_Q N_Q = 1` from `\overline\Gamma_5/(2G/(5c^5)) = \chi_Q N_Q` and the source-map relation `\widehat m_0 → 1`. This path uses entirely different intermediate variables (`Y_Q^{\rm ret}`, `\sigma_Q^{\rm can}`, series expansion) from the SymPy path.

**Verification:**
After re-authoring: variable names, intermediate-quantity choreography, and order of operations should be visibly different from the SymPy script while reaching the same bottom-line `N_Q = 1` on the canonical point-particle branch.

### F4 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:38, 61-65`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:40-43, 66-67`

**What's wrong:**
The two "real" assertions are A1 (`NQ_canonical.subs(m0hat,1) - 1`) and A4 (`gamma_eff_canonical - Gamma5_target`). Both follow trivially from the algebra: A1 reduces to `1/1^2 - 1 = 0` and A4 reduces to `m0hat^2 · (1/(chi_Q m0hat^2)) · Gamma5_target |_{chi_Q=1} - Gamma5_target = Gamma5_target - Gamma5_target = 0`. Once the constraint `m0hat^2·chi_Q·N_Q=1` is solved for N_Q, these are forced identities at the substituted values, not independent checks. The script never exercises the more interesting question — what happens when chi_Q ≠ 1, i.e., what `\Delta_Q := chi_Q - 1` (appendix eq. app-part04-DeltaQ-def, line 293) does to N_Q, Γ5, K0,…, which the appendix says is "the entire reduced failure" at leading order (line 342). A more substantive check would series-expand the deviation `\gamma_{\rm quad}^{\rm eff} - 2G/(5c^5)` to first order in `\Delta_Q` and assert the coefficient equals `2G/(5c^5)` (i.e., `gamma_eff = 2G/(5c^5) · (1 + \Delta_Q) + O(\Delta_Q^2)` on the natural source-map branch).

**Why this matters:**
The current script "passes" trivially: the assertions cannot fail. They exercise neither the canonical-branch realization (which is what Checks ii and iii address; see F1) nor the off-canonical sensitivity (which is the conditional theorem's content per app-part04 line 342). Without one of these, the audit transcript provides no evidence of physical content — only of substitution algebra.

**Required change:**
Add a sensitivity assertion: define `gamma_eff_off = (m0hat^2 * Gamma5)` (no chi_Q→1 substitution), expand to first order in `Delta_Q = chi_Q - 1` about `chi_Q=1` on the natural source-map branch `m0hat=1`, and assert the linear coefficient equals `-2G/(5c^5)` (i.e., `gamma_eff = 2G/(5c^5) - 2G/(5c^5)·\Delta_Q + O(\Delta_Q^2)`). The minus sign comes from `N_Q = 1/chi_Q = 1 - \Delta_Q + \Delta_Q^2 - ...` at m0hat=1, so `gamma_eff = N_Q · Gamma5_target = Gamma5_target·(1 - \Delta_Q + ...)`, so the first-order coefficient is `-Gamma5_target = -2G/(5c^5)`. This exercises real algebra: a sign or factor mistake in N_Q would change this coefficient.

**Verification:**
After fix: a new `expect_zero` of the form `expect_zero("Delta_Q first-order sensitivity", series_coeff(gamma_eff at m0hat=1, Delta_Q, 1) - (-2*G/(5*c**5)))` should appear (or analogous expansion) and pass. The Mathematica mirror should be added independently.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script — see F3 for quoted side-by-side excerpts. Both scripts use the identical variable choreography (`m0hat`, `chi_Q`/`chiQ`, `N_Q`/`NQ`, then four `*_target` constants, then `K0/K2/K4/Gamma5`), the identical sequence of three `expect_zero`/`expectZero` calls, and identical RESULT prose. There is no second-engine independence.

## Engine cross-check

The two engines do agree on the algebra they share. Both produce the identical residuals:

SymPy output (`scripts/output/...:14-23`):
```
N_Q on the general outgoing branch = 1/(chi_Q*m0hat**2)
N_Q on the canonical outgoing branch = m0hat**(-2)
point-particle canonical branch gives N_Q = 1 = 0
K0 = 54*G*c_s**5/(5*a**5*c**5*chi_Q*m0hat**2)
K2 = 6*G*c_s**3/(5*a**3*c**5*chi_Q*m0hat**2)
K4 = 8*G*c_s/(15*a*c**5*chi_Q*m0hat**2)
Gamma5 = 2*G/(5*c**5*chi_Q*m0hat**2)
branch identity K4 - 4 K2^2 / K0 = 0
branch identity Gamma5 - 9 K2^(5/2)/K0^(3/2) = 0
canonical gamma_eff - target = 0
```

Mathematica output (`mathematica/output/...:14-27`):
```
N_Q on the general outgoing branch = 1/(chiQ*m0hat^2)
N_Q on the canonical outgoing branch = m0hat^(-2)
point-particle canonical branch gives N_Q = 1 = 0
K0 = (54*cS^5*G)/(5*a^5*c^5*chiQ*m0hat^2)
K2 = (6*cS^3*G)/(5*a^3*c^5*chiQ*m0hat^2)
K4 = (8*cS*G)/(15*a*c^5*chiQ*m0hat^2)
Gamma5 = (2*G)/(5*c^5*chiQ*m0hat^2)
branch identity K4 - 4 K2^2 / K0 = 0
branch identity Gamma5 - 9 K2^(5/2)/K0^(3/2) = 0
canonical gamma_eff - target = 0
```

Engines agree, but since the Mathematica script is a transliteration this agreement is weak evidence.

## Verdict justification

`verdict: findings`. The two assertions that are not tautological are also not substantive: they verify forced consequences of the constraint `m0hat^2·chi_Q·N_Q=1` at the substituted values, but they do not exercise the canonical-branch derivation (which lives in the DtN z-expansion the script never constructs), nor the off-canonical sensitivity (the `\Delta_Q` linearization that measures the reduced theorem's failure mode). The K4-vs-K2 "branch identity" check is a definitional identity, not a physics test. The Mathematica script duplicates the SymPy algorithm step for step, so the engine agreement is structurally guaranteed. Three of the four findings are script-side fixable (F2 tautological-check rewrite, F3 Mathematica re-derivation, F4 sensitivity assertion). The fourth (F1) is a paper_misalignment requiring user resolution: paper Checks (ii) and (iii) are not exercised by the scripts, and either the scripts must be extended or the card must be corrected to reference an upstream stage that owns those checks. `stop_cold` is null because no finding's correction propagates downstream irreversibly — Stage 106 imports `chi_Q=1` from upstream Stage 88, so fixing the script-side coverage here does not change any output downstream stages depend on; it only adds verification depth.

Secondary observation (not a separate finding): both scripts carry stale banner labels ("STAGE 89" / "STAGE 089" in sympy line 25 and mathematica line 27), and the paper card section heading reads "Stage~123" (stage_106.tex:1) — three different stage numbers for the same unit. The math content matches stage_106 correctly; only the labels are stale. Not blocking.

## Self-test notes

I checked: (1) the K4 - 4K2^2/K0 "branch identity" simplifies to 0 by pure substitution of the script's K2,K4 definitions — confirmed tautological; (2) the Gamma5 - 9·sqrt(K2^5/K0^3) check evaluates to 0 by substituting K0_target = 54Gc_s^5/(5a^5c^5) and Gamma5_target = 2G/(5c^5), confirming the four hardcoded literals are mutually consistent (which is the only physics content in any assertion); (3) the proposed F4 first-order Delta_Q sensitivity coefficient of gamma_eff at m0hat=1 — expanding N_Q = 1/(1+Delta_Q) = 1 - Delta_Q + Delta_Q^2 - ... and multiplying by Gamma5_target gives leading coefficient -2G/(5c^5), confirmed non-zero (so an assert_nonzero on the linear coefficient or assert_zero on the sum coefficient + 2G/(5c^5) would both behave correctly); (4) the proposed F2 fix re-derives the K0·K4 = 4K2^2 identity from the four `*_target` literals — substituting numerically: 4*(6/5)^2 = 144/25, and (54/5)*(8/15) = 432/75 = 144/25 ✓, so the rewritten assertion is non-tautological and passes only because the four literal coefficients happen to satisfy the identity. (5) Paper round-trip on F4: the deviation coefficient `-2G/(5c^5)` matches the appendix's statement that `\Delta_Q = \chi_Q - 1` measures "the entire reduced failure" at leading order (line 342), so no new paper misalignment is introduced.
