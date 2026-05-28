---
unit_id: 164
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage164_microscopic_log_channels.md]
  paper_appendix: present
---

# Audit unit 164 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_164.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage164_microscopic_log_channels.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at L59, L283-360, L411)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage164_microscopic_log_channels_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage164_microscopic_log_channels_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Rewrites the Stage 163 off-family scalar in explicit branch-variable drifts of \((\mathcal Z_q, \rho_w, c_{s,w}, c_s, \mathcal T_m, v_{w0}, a, L_W)\)." The card and appendix (eq:app-part05-log-channels-rg, eq:app-part05-delta-perp-explicit, eq:app-part05-ABC-defs) state the stage proves: (1) the two Stage-163 logarithmic channels equal `dln(g/r)` and `-dln r_c` exactly via the parent-ratio definitions `r=lam/sqrt(Ks Kq)`, `g=gq sqrt(Ks)/(gs sqrt(Kq))`, `rc=lam^2/(Ks Kq)`; (2) on the general parent-action formulas, then under uniform-overlap + D/N closure, then under the healing-locked shell branch, the two products collapse to the explicit monomials `g_q K_s/(g_s lam) = -(27 pi m_psi^2/(40 hbar mu_0 q_*)) Z_q c_sw^3/(rho_w T_m v_w0 a^2 L_W^2)` and `K_s K_q/lam^2 = (27 pi^3 m_psi^2/(320 hbar mu_0 q_*^2)) Z_q c_s^2 c_sw^3/(rho_w v_w0^2 a^2 L_W^3)`; (3) `delta_perp = g_* dln(...) + (1/(4 sqrt(1+r_*^2))) dln(...)` reduces to the compressed `A_*,B_*,C_*` form; (4) the numeric Family-1 coefficients at `g_* ≈ 0.758035…`, `r_* ≈ 1.77799…` (e.g. `A_*≈0.880589`, `B_*≈0.245108`, L_W coeff `≈ -1.88373`). The notes (Secs. 1-6) supply each derivation step and the tangency-law solve for `dln T_m`.

## What the script claims to verify

The SymPy script asserts (`expect_zero` = `simplify(expand(.)) == 0`): the two parent-ratio identities `g/r = gq Ks/(gs lam)` and `1/rc = Ks Kq/lam^2` (L41-42); the healing-locked product monomials against the literal expected forms (L101-102); and the compressed `delta_perp` form against the hand-written `A_*/B_*/C_*` combination (L130). It then prints (no assert) the general and uniform products, the tangency-law solution, and the numeric Family-1 coefficients. The Mathematica script asserts the same parent-ratio identities (L57-58), adds asserts at the uniform stage (L93-94) the SymPy lacks, asserts the healing products (L115-116) and the compressed `delta_perp` (L133), and adds five `expectApprox` numeric cross-checks (L162-166) on `A_*`, `B_*`, `3A_*`, the `v_w0` coeff, and the `L_W` coeff against the paper's published values.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) channels = `dln(g/r)`, `-dln rc` via parent ratios | sympy L41-42, math L57-58 `expect_zero` | match |
| (2a) general products | printed only (sympy L64-65, math L70-71); math also asserts at uniform L93-94 | partial (general printed, downstream stages asserted) |
| (2b) healing-locked product monomials | sympy L101-102, math L115-116 `expect_zero` vs literal expected | match |
| (3) compressed `delta_perp` `A_*/B_*/C_*` form | sympy L130, math L133 `expect_zero` | match |
| (4) numeric Family-1 coefficients | sympy printed only; math `expectApprox` L162-166 | match (asserted on math side) |
| tangency-law `dln T_m` solve | printed only (sympy L134, math L140) | partial (printed, follows from delta_perp=0) |

`paper_alignment: aligned` — every load-bearing identity in the card/appendix is asserted, with correct constants and parameter forms. The general-product and tangency rows are print-only but each is a strict algebraic consequence of an asserted neighbor (the healing products are asserted; the tangency law is the linear solve of `delta_perp=0`, whose form is asserted), so no deliverable is left unverified.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41 | `expect_zero(g/r - gq Ks/(gs lam))` | claim 1 | yes |
| A2 | sympy | 42 | `expect_zero(1/rc - Ks Kq/lam^2)` | claim 1 | yes |
| A3 | sympy | 101 | `expect_zero(first_prod_heal - first_expected)` | claim 2b | yes |
| A4 | sympy | 102 | `expect_zero(second_prod_heal - second_expected)` | claim 2b | yes |
| A5 | sympy | 130 | `expect_zero(delta_perp - delta_perp_expected)` | claim 3 | yes |
| M1 | mathematica | 57 | `expectZero(g/r - gq ks/(gs lam))` | claim 1 | yes |
| M2 | mathematica | 58 | `expectZero(1/rc - ks kq/lam^2)` | claim 1 | yes |
| M3 | mathematica | 93 | `expectZero(firstProdUniform - firstUniformExpected)` | claim 2a→2b bridge | yes |
| M4 | mathematica | 94 | `expectZero(secondProdUniform - secondUniformExpected)` | claim 2a→2b bridge | yes |
| M5 | mathematica | 115 | `expectZero(firstProdHeal - firstHealExpected)` | claim 2b | yes |
| M6 | mathematica | 116 | `expectZero(secondProdHeal - secondHealExpected)` | claim 2b | yes |
| M7 | mathematica | 133 | `expectZero(deltaPerp - deltaPerpExpected)` | claim 3 | yes |
| M8 | mathematica | 162-166 | `expectApprox(...)` x5 | claim 4 | yes |

All assertions are non-tautological: A3/A4/M5/M6 derive the product by an independent substitution chain (`g_q,K_q,lam,K_s,ell` → explicit formulas → uniform overlap → healing lock) and then compare against a *separately written* closed-form literal — if any substitution or simplification were wrong, the residual would be nonzero. A5/M7 build `delta_perp` two ways (weighted-sum of channel formulas vs. compressed `A_*/B_*/C_*` form) and compare; the equality is the physics content, not a definitional identity. M8 compares script-derived `aNum/bbNum` (built from `gNum,rNum`) against the paper's published constants.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage164_microscopic_log_channels_sympy_audit.py:30-166`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.wl:32-173`

**What's wrong:**

The Mathematica `.wl` is a near-line-for-line transliteration of the SymPy `.py`: identical variable choreography, identical intermediate-step ordering, identical hand-written target forms, and the same copy-paste banner typo. The math here is not "short" (unlike Stage 162) — it is a multi-stage algebraic reduction chain — so the structural parallelism cannot be dismissed as the irreducible minimum of a second engine. Three corresponding sections:

(a) Both open with the identical mislabelled banner `banner("STAGE 147 — MICROSCOPIC LOG-IMBALANCE CHANNELS")` (sympy L30, wl L32) for a Stage-164 unit. A preserved wrong stage number is the classic transliteration tell.

(b) The healing-locked products are built by the identical substitution choreography. SymPy L89-93:
```
ell_lock = hbar / (2 * mpsi * csw)
Ks_lock = sp.simplify(3 * sp.pi * a**2 * hbar**2 / (5 * mpsi * rho_w * ell_lock))
first_prod_heal = sp.simplify(first_prod_uniform.subs({Ks: Ks_lock, ell: ell_lock}))
second_prod_heal = sp.simplify(second_prod_uniform.subs({Ks: Ks_lock, ell: ell_lock}))
```
Mathematica L98-108 mirrors this exactly:
```
ellLock = hbar/(2*mpsi*csw);
ksLock = FullSimplify[3*Pi*a^2*hbar^2/(5*mpsi*rhoW*ellLock), ...];
firstProdHeal = FullSimplify[firstProdUniform /. {ks -> ksLock, ell -> ellLock}, ...];
secondProdHeal = FullSimplify[secondProdUniform /. {ks -> ksLock, ell -> ellLock}, ...];
```

(c) The compressed `delta_perp` target is hand-written identically in both. SymPy L120-129 and Mathematica L123-132 construct `firstHeal`, `secondHeal`, `deltaPerp = gstar*firstHeal + b*secondHeal`, and `deltaPerpExpected = (gstar+b)*(dlnZ-dlnrho) + 3*(gstar+b)*dlncsw + 2*b*dlncs - gstar*dlnTm - (gstar+2*b)*dlnv - 2*(gstar+b)*dlna - (2*gstar+3*b)*dlnLw` term-for-term in the same order with the same symbol names. The two scripts therefore verify the same algebra through the same path; a transcription or simplification slip in one engine's choreography would be reproduced, not caught, by the other.

**Why this matters:**

The two-engine policy exists so that an algebra slip or CAS-specific simplification artifact in one engine is caught by an independent re-derivation in the other. A line-for-line port makes cross-engine agreement structurally guaranteed rather than an independent confirmation, so the Mathematica run provides no real second-source assurance for the load-bearing healing-locked monomials and the `delta_perp` compression.

**Required change:**

Edit only the Mathematica script. Keep all existing assertions (they may stay and continue to pass) and add independent-derivation assertions that obtain the two healing-locked log-channel coefficient vectors and the compressed `delta_perp` by a route with no SymPy counterpart (a `Series`/`Coefficient` perturbation of the explicit healing-locked product monomials), then reconcile against the hand-written `firstHeal`/`secondHeal`/`deltaPerpExpected`. See the directive for the exact block. Also correct the banner stage number at wl L32.

**Verification:**

After Codex applies, `redteam exec-mathematica 164` should still exit 0, and new `PASS` lines should appear: `first channel via series route`, `second channel via series route`, and `delta_perp via series route`. All pre-existing `PASS` lines must remain.

## Independent-derivation check (Mathematica)

Not independent. The `.wl` reproduces the `.py` step-for-step: same symbol set (`ks,kq,lam,gs,gq,tm,js,zq,mu0,lw,qstar,vw0,isq,cs,a,ell,mpsi,hbar,rhoW,csw,...`), same six-section structure in the same order, the same `gsGen/gqGen/kqGen/lamGen` definitions, the same `jsClosure/iq/isqClosure` closure, the same `ellLock/ksLock` healing lock, the same hand-typed `deltaPerpExpected`, and the same preserved "STAGE 147" banner typo. The only Mathematica-only additions (`firstUniformExpected`/`secondUniformExpected` asserts at L88-94 and the `expectApprox` numeric checks at L162-166) are extra asserts layered on the identical algebraic path, not an independent derivation of the core monomials. → `mathematica_transliteration` (F1).

## Engine cross-check

Both engines agree at the level they claim. The parent-ratio residuals are `0` in both (sympy L13-14, math L13-16). The healing-locked products print identically (`-27 pi Z_q c_sw^3 mpsi^2/(40 L_W^2 T_m a^2 hbar mu_0 qstar rho_w v_w0)` and `27 pi^3 Z_q c_s^2 c_sw^3 mpsi^2/(320 ...)`) and both residuals are `0`. The compressed `delta_perp` residual is `0` in both. Numeric coefficients agree to <5e-16: `A_*=0.880589090041570`, `B_*=0.245108022193813`, `coeff[ln c_sw]=2.641767270124709`, `coeff[ln a]=-1.761178180083139`, `coeff[ln L_W]=-1.883732191180046`, matching the notes (Sec. 5) exactly. No `engine_disagreement`.

## Verdict justification

The math holds against the paper. I independently re-derived the uniform-overlap and healing-locked monomials (`I_sq = 8 a^2 ell sqrt(2 L_W)/3`; `K_s_lock = 6 pi a^2 hbar c_sw/(5 rho_w)`; first product coeff `-27 pi/40`, second `27 pi^3/320`) and confirmed they match `first_expected`/`second_expected` and the notes' boxed forms; the log-channel coefficient vectors, the `A_*/B_*/C_*` compression, and the numeric Family-1 constants all match. Attacks that failed: factor/sign of `pi`, `sqrt2`, and the `27/40` vs `27/320` coefficients (correct); the `b=1/(4 sqrt(1+r_*^2))` weight vs `B_*=2b` (consistent); positivity assumptions (all symbols are physical magnitudes; `lam_gen=-q_* v_w0 I_sq` is the only signed quantity and its sign is carried correctly into the `-27/40` product). The `g_*`/`r_*` literals are upstream Family-1 carry-forwards stated verbatim in the notes/appendix, not in-stage results, so no `hardcoded_result`. The only real finding is structural: the Mathematica file is a transliteration of the SymPy file (F1), reproducing the same algebraic path and the same "STAGE 147" banner typo, so it is not the independent second engine the policy requires. Not `clean`. Not `UNFIXABLE` (the identities are correct; only an independent route must be added) and not `CRITICAL_DOWNSTREAM` (the verified monomials/coefficients are correct, so adding an independent check changes no downstream value). The banner mislabel itself is a codebase-wide renumbering artifact (offset ~17, matching 160→143, 161→144, 162→145, 163→146) consistent with prior reports; it is folded into F1 as a transliteration tell + a one-line label fix rather than a standalone finding.

## Self-test notes

I checked (1) variable independence: no `D[]`/`diff` over an absent variable — the new check uses `Coefficient[Series[ratio,{eps,0,1}],eps]` on monomials, which is well-posed. (3) trivial-case pre-check: with all `dln*=0` the series ratio is `1`, so its O(eps) coeff is `0` (matches `first_heal=0`); setting only `dlnZ=1` gives O(eps) coeff `1` for the first channel and `1` for the second, matching both hand forms' `+dlnZ` coefficient; the `c_sw^3` factor yields `+3 dlncsw`, the `L_W^2`/`L_W^3` factors yield `-2 dlnLw`/`-3 dlnLw` — all consistent with `first_heal`/`second_heal`. (5) paper round-trip: the proposed independent route reuses the already-asserted `firstHealExpected`/`secondHealExpected` monomials as its input and reconciles to the same `deltaPerpExpected`, introducing no new constant, so it cannot create a new `paper_misalignment`. Using the product *ratio* (not `Log`) sidesteps the negative-sign branch issue on the first (negative) product.
