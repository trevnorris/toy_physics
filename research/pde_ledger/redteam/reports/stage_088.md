---
unit_id: 088
batch: III.5
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage088_loading_ratio_from_minimal_module.md]
  paper_appendix: present
---

# Audit unit 088 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_088.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage088_loading_ratio_from_minimal_module.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 154 and 294 referencing stage 088)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.txt`

## What the paper claims

Stage 088 reads the minimal isotropic conservative quadrupole precursor `Y_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` as a natural contact-plus-pole expansion `Y_Q^cons = 1/rho_alpha + (rho_alpha-1)/rho_alpha * 1/(1 - omega^2/Omega_Q^2)` and matches coefficients to extract the loading ratio. `\stagefield{Output}{The minimal-module ratio \eqref{eq:app-stage088-rho}.}` — namely `rho_alpha = 1/c0 = 4/3`, `zeta_req = c1/c0 = 1/3`, and `Pi_tr = (4/3) C_mix`. The card adds the regime classification `C_mix < Pi_tr < 2 C_mix` (symmetric-lowest-twin regime), `0 < zeta_req = 1/3 < 1`. The notes additionally specify the Stage 68 identity `Pi_tr / C_mix = alpha_req / alpha_mix = rho_alpha` as the input that converts the rho_alpha result into the Pi_tr product-language statement.

## What the script claims to verify

Both scripts:
1. Write down `Y_loading` (in alpha_mix/alpha_req variables) and `Y_rho` (in rho_alpha) and verify they coincide under `alpha_req -> rho_alpha * alpha_mix`.
2. Verify the contact-plus-pole identity `Y_rho == 1/rho_alpha + (rho_alpha-1)/rho_alpha * 1/(1 - omega^2/Omega_Q^2)` and `c0+c1=1`.
3. Verify the round-trip `rho(c0(rho)) = rho`, `rho(c1(rho)) = rho`, `zeta(c(rho)) = rho - 1`.
4. Substitute c0=3/4 and c1=1/4 and assert `rho_min = 4/3`, `zeta_min = 1/3`.
5. Assert `Pi_tr / C_mix = 4/3` after defining `Pi_tr := (4/3) C_mix` directly (in both engines).
6. Mathematica only: `expectTrue["C_mix < Pi_tr < 2 C_mix", cMix < (4/3)*cMix < 2*cMix]`.

## Paper <-> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `rho_alpha = 4/3` | rho_min = 1/(3/4) checked = 4/3 (sympy L78, math L64) | match |
| `zeta_req = 1/3` | zeta_min = (1/4)/(3/4) checked = 1/3 (sympy L79, math L65) | match |
| `Pi_tr = (4/3) C_mix` | `(4/3)*C_mix / C_mix - 4/3 == 0` (sympy L83, math L66) | mismatch — tautological; the link via Stage 68 `Pi_tr/C_mix = rho_alpha` is not exercised; the script just hardcodes Pi_tr = (4/3)C_mix |
| Regime `C_mix < Pi_tr < 2 C_mix` | sympy: print only (L86-87); math: expectTrue (L67) | partial — sympy emits prose; math checks 1 < 4/3 < 2 (numeric tautology) |
| Contact-plus-pole identification (eq. \ref{eq:app-stage088-contact-pole}) | algebraic identity asserted = 0 (sympy L41-44, math L44) | match (algebraic identity, holds by construction; modest content) |

`paper_alignment: aligned` — every paper deliverable maps to *some* script-side line, the numeric answers agree, and no value disagreement exists. Two of the entries are weakened by tautology, captured below as findings.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41-44 | `simplify(Y_loading.subs(alpha_req,rho_alpha*alpha_mix) - Y_rho) == 0` | contact-plus-pole identification | partial (identity by construction) |
| A2 | sympy | 49-52 | `Y_rho - (c0_from_rho + c1_from_rho/(1-omega^2/Omega_Q^2)) == 0` | contact-plus-pole identification | no (tautological — c0_from_rho, c1_from_rho are defined as those very pieces of Y_rho) |
| A3 | sympy | 57 | `c0_from_rho + c1_from_rho - 1 == 0` | normalized static limit | no (tautological — 1/r + (r-1)/r = 1 algebraically by construction) |
| A4 | sympy | 67 | `rho_from_c0.subs(c0,c0_from_rho) - rho_alpha == 0` | inverse map | no (tautological — substituting 1/r into 1/c gives r) |
| A5 | sympy | 68 | `rho_from_c1.subs(c1,c1_from_rho) - rho_alpha == 0` | inverse map | no (tautological — by construction) |
| A6 | sympy | 69 | `zeta_from_c.subs(...) - (rho-1) == 0` | inverse map | no (tautological — by construction) |
| A7 | sympy | 78 | `rho_min - 4/3 == 0` after `c0 -> 3/4` | `rho_alpha = 4/3` | yes (arithmetic: 1/(3/4)=4/3) |
| A8 | sympy | 79 | `zeta_min - 1/3 == 0` after `c0->3/4, c1->1/4` | `zeta_req = 1/3` | yes (arithmetic: (1/4)/(3/4)=1/3) |
| A9 | sympy | 83 | `(4/3)*C_mix / C_mix - 4/3 == 0` | `Pi_tr = (4/3) C_mix` | no (tautological — Pi_tr literal substituted on both sides) |
| B1 | mathematica | 44 | `(yLoading /. alphaReq->rhoAlpha*alphaMix) - yRho == 0` | contact-plus-pole identification | partial |
| B2 | mathematica | 48 | `yRho - (c0FromRho + c1FromRho/(1-omega^2/omegaQ^2)) == 0` | contact-plus-pole identification | no (tautological, mirrors A2) |
| B3 | mathematica | 49 | `c0FromRho + c1FromRho - 1 == 0` | normalized static limit | no (mirrors A3) |
| B4 | mathematica | 54-56 | round-trip identities | inverse map | no (mirror A4-A6) |
| B5 | mathematica | 64 | `rhoMin - 4/3 == 0` | `rho_alpha = 4/3` | yes |
| B6 | mathematica | 65 | `zetaMin - 1/3 == 0` | `zeta_req = 1/3` | yes |
| B7 | mathematica | 66 | `piMin/cMix - 4/3 == 0` with piMin = (4/3)*cMix | `Pi_tr = (4/3) C_mix` | no (mirrors A9) |
| B8 | mathematica | 67 | `cMix < (4/3)*cMix < 2*cMix` | regime classification | no (literal numeric trichotomy 1<4/3<2) |

The only assertions that exercise something the math could "get wrong" are A7/B5 (`1/(3/4) == 4/3`) and A8/B6 (`(1/4)/(3/4) == 1/3`). Everything else either restates a definition or substitutes a hardcoded constant.

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py:83`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl:60,66`

**What's wrong:**
The paper claims `Pi_tr = (4/3) C_mix` follows from `rho_alpha = 4/3` via the Stage-68 identity `Pi_tr / C_mix = alpha_req / alpha_mix = rho_alpha` (notes section 2: "since Stage 68 proved `Pi_tr / C_mix = alpha_req / alpha_mix`, the same result is `Pi_tr = (4/3) C_mix`"). The SymPy script instead writes `expect_zero("Pi_tr/C_mix - 4/3", sp.simplify((sp.Rational(4,3)*C_mix)/C_mix - sp.Rational(4,3)))` — i.e. it labels the check "Pi_tr/C_mix - 4/3" while in fact substituting `Pi_tr = (4/3) C_mix` literally on the left. The expression `(4/3*C)/C - 4/3` is identically zero with no physics in play. Mathematica mirrors this at L60 (`piMin = FullSimplify[(4/3)*cMix]`) and L66 (`piMin/cMix - 4/3`).

**Why this matters:**
The script presents this as verification of the Pi_tr identity, but the Stage-68 identity (`Pi_tr/C_mix = rho_alpha`) is never invoked. If the Stage-68 identity were ever wrong, this assertion would still pass because Pi_tr is hardcoded.

**Required change:**
Replace the literal hardcoding with an explicit symbolic chain that exercises rho_alpha: define `Pi_tr_from_rho = rho_alpha * C_mix` (this is the Stage-68 identity, the only thing being assumed here), substitute `rho_alpha -> rho_min` (where rho_min was just derived from c0=3/4 above), and then assert `Pi_tr_from_rho - (4/3)*C_mix == 0`. That way the check actually exercises the chain `c0=3/4 -> rho_min=4/3 -> Pi_tr=(4/3)C_mix` rather than restating the answer.

**Verification:**
The new assertion residual still reduces to 0 (since 4/3 * C_mix - (4/3)*C_mix = 0), but now mutating `rho_min` (e.g. swapping c0=1/2) would propagate to a nonzero residual.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py:49-52,57,67-69`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl:46-49,54-56`

**What's wrong:**
The script defines `c0_from_rho = 1/rho_alpha`, `c1_from_rho = (rho_alpha-1)/rho_alpha`, then asserts `Y_rho - (c0_from_rho + c1_from_rho/(1-omega^2/Omega_Q^2)) == 0`. But `Y_rho` is *defined* (L36) as exactly `1/rho_alpha + (rho_alpha-1)/rho_alpha/(1-omega^2/Omega_Q^2)`, so the residual is identically zero by construction. Similarly `c0_from_rho + c1_from_rho - 1` is `1/r + (r-1)/r - 1 == 0` algebraically, and the round-trip identities `rho_from_c0.subs(c0, c0_from_rho) - rho_alpha`, etc. are forced to zero by the same definitions. None of these can fail no matter what the physics is.

**Why this matters:**
The script appears to verify the contact-plus-pole identification but actually only verifies that algebraic simplification is consistent. The only substantive content of the stage — that the natural reading of the precursor maps c0 to 1/rho_alpha — is asserted by definition, not derived.

**Required change:**
Either (a) keep these as labeled prints (not `expect_zero` assertions) since they are book-keeping identities, or (b) re-derive `c0` and `c1` by an independent path and then assert they equal `1/rho_alpha` and `(rho_alpha-1)/rho_alpha`. The simplest substantive check: take the partial-fraction decomposition `apart(Y_rho, omega)` and read off the constant term and pole residue, then assert those equal `c0_from_rho` and `c1_from_rho`. That exercises the *coefficient-matching* claim of the stage rather than restating the parameterization.

**Verification:**
After the fix, the assertion should still pass but should exercise sympy's `apart`/`series` machinery; the residue extraction must agree with `(rho_alpha-1)/rho_alpha`, and the static-limit constant with `1/rho_alpha`.

### F3 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl:39-67`

**What's wrong:**
The Mathematica script is structurally a line-by-line transliteration of the SymPy script:
- sympy L35 `Y_loading = sp.simplify(alpha_mix/alpha_req + (alpha_req-alpha_mix)/alpha_req/(1 - omega**2/Omega_Q**2))` ↔ math L39 `yLoading = FullSimplify[alphaMix/alphaReq + ((alphaReq - alphaMix)/alphaReq)/(1 - omega^2/omegaQ^2), ...]` — same expression in same order.
- sympy L46-47 `c0_from_rho = sp.simplify(1/rho_alpha); c1_from_rho = sp.simplify((rho_alpha-1)/rho_alpha)` ↔ math L46-47 `c0FromRho = FullSimplify[1/rhoAlpha,...]; c1FromRho = FullSimplify[(rhoAlpha - 1)/rhoAlpha, ...]` — same construction.
- sympy L59-69 round-trip checks (rho_from_c0, rho_from_c1, zeta_from_c, three asserts) ↔ math L51-56 — same five checks in the same order, with the same intermediate symbol names rewritten in camelCase.
- Assertion list 1-9 in both engines is identical in order, identical in label, identical in form.

**Why this matters:**
The second engine is supposed to derive the claim independently from the physical premises. Here it follows the SymPy script's algebraic choreography step-for-step, so a derivation error in the SymPy script (e.g., wrong parameterization of the contact-plus-pole form) would be silently reproduced by Mathematica. (This stage's central claim is admittedly thin enough that "independent derivation" is borderline, but the verbatim mirroring is the textbook case.)

**Required change:**
Restructure the Mathematica script to attack the claim from a different direction. Suggested alternative path:
1. Start from `YQ = 3/4 + (1/4)/(1 - omega^2/omegaQ^2)` (the literal paper-Input form).
2. Use `Apart[YQ, omega]` or extract `Limit[YQ, omega -> 0]` for c0 and `Residue` at omega = omegaQ scaled appropriately for c1, to recover c0=3/4 and c1=1/4 independently.
3. Compute `rhoAlpha = 1/c0`, `zetaReq = c1/c0`, and assert these equal 4/3 and 1/3 respectively.
4. Verify the Stage-68 product identity `Pi_tr = rhoAlpha * cMix` and substitute to get `(4/3) cMix`.

The two engines then converge on the same numeric answer through structurally different machinery (Apart/Residue/Limit vs Eq-equation substitution).

**Verification:**
After the fix, the Mathematica file's intermediate symbols and call order should not align line-for-line with the SymPy file; the final printed numeric assertions (rho_alpha=4/3, zeta_req=1/3) must still pass.

## Independent-derivation check (Mathematica)

The Mathematica script is a line-by-line port of the SymPy script (see F3 for the side-by-side correspondence). The same variable choreography (`Y_loading`, `Y_rho`, `c{0,1}_from_rho`, `rho_from_c{0,1}`, `zeta_from_c`), the same assertions in the same order, the same labels. No independent algebraic path — `Apart`, `Residue`, `Series`, or `Limit` are nowhere invoked. This is `mathematica_transliteration`.

## Engine cross-check

Both engines exit 0 and print the same numerical answers: `rho_alpha = 4/3`, `zeta_req = 1/3`, `Pi_tr/C_mix - 4/3 = 0`. Agreement is real but trivial — both engines hardcode the same constants and run the same algebra. No engine disagreement.

## Verdict justification

The math the paper claims (rho_alpha = 4/3, zeta_req = 1/3, Pi_tr = (4/3) C_mix, symmetric-lowest-twin regime) is correct. Both scripts exit 0 with answers consistent with the paper. **However**: most of the script's "checks" are tautological by construction — `c0+c1=1` is `1/r + (r-1)/r = 1`, the round-trip checks are forced by their definitions, and the load-bearing `Pi_tr/C_mix = 4/3` check substitutes `Pi_tr = (4/3) C_mix` literally instead of invoking the Stage-68 identity that the paper says supplies the link. The only assertions that exercise real arithmetic are `1/(3/4) = 4/3` and `(1/4)/(3/4) = 1/3`. The Mathematica script is a structural transliteration of the SymPy script. Findings: F1 (Pi_tr tautology), F2 (contact-plus-pole identity tautology), F3 (Mathematica transliteration). Paper alignment is `aligned`; the verdict is `findings` because the verification is thinner than the paper's claim warrants but no value disagrees with the paper.

## Self-test notes

I checked (1) variable-independence traps: there are no `sp.diff` or `D[expr, var]` calls in this script — no derivative dependence to mis-thread; (2) symmetry/parity: no unbounded integrals in this stage; (3) the trivial-case pre-check on the proposed F1 fix: with `rho_min = 4/3`, `rho_min * C_mix - (4/3)*C_mix = (4/3)C_mix - (4/3)C_mix = 0` ✓, and mutating rho_min to e.g. 5/4 would yield `(5/4 - 4/3)*C_mix = -(1/12)*C_mix ≠ 0` ✓; (4) the proposed F2 fix uses `sp.apart`, which on `Y_rho = (Omega_Q^2 rho_alpha - omega^2)/(rho_alpha(Omega_Q^2 - omega^2))` against `omega` should yield `1/rho_alpha + ((rho_alpha-1)/rho_alpha)/(1-omega^2/Omega_Q^2)` — confirming c0=1/rho_alpha, c1=(rho_alpha-1)/rho_alpha; (5) paper round-trip: none of F1-F3's proposed changes touch the numeric values (4/3, 1/3, c0=3/4, c1=1/4) — only the algebraic path that arrives at them — so no new paper_misalignment is introduced.
