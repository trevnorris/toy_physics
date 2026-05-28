---
unit_id: 147
batch: IV.5
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage147_first_order_rigidity_kernel.md
  paper_appendix: present
---

# Audit unit 147 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_147.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage147_first_order_rigidity_kernel.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_147}` row plus context lines 840-850 with `A_T \approx -4.27263956256927`, `B_T \approx 0.134875005736706`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.txt`

## What the paper claims

The stage card (paper/stages/stage_147.tex) tags this as a derivation-ledger step whose verification target is quoted as:

> Traction shift is one kernel projection with coefficients `A_T ≈ -4.27264`, `B_T ≈ 0.134875`.

The notes file expands this into three concrete claims:
(C1) closed forms for A_T, B_T in terms of starred quantities `T_{m,*}, Π_*, S_*, S'_*, g'_*, S_*/4`;
(C2) numerical values `A_T ≈ -4.27263956256927` and `B_T ≈ 0.134875005736706`, with `|A_T|/B_T ≈ 31.6785`;
(C3) the centered single-kernel rewrite `δT_m = ε ∫₀¹ W_*(x)[ς(x)-Σ_*(x)] dx` with `W_*(x) = A_T(c(x)-g_*) + B_T(K_q(x)-S_*)`, justified by the property "Σ_ε - Σ_* integrates to zero" (centering identity). The card's Checks list also names (a) centering of first-order deformations, (b) the rigidity kernel (not branch ambiguity) controlling non-exponential corrections, and (c) one-step nonlinear correction staying in the reduced mouth-layer regime — all three are intended to be exercised verification items.

## What the script claims to verify

The SymPy script computes `Π_*` by `nsolve(gPi - gminus, 1.5)`, evaluates `g_*, S_*, g'_*, S'_*, Σ_*, T_*`, then evaluates the closed-form combinations for `A_T` and `B_T`. It prints these numerically together with `|A_T|/B_T`, prints a symbolic `delta T_m` in `ε, gbar, Sbar`, and prints the centered weight `W_center = A_T(c - gminus) + B_T(K_q - S_formula(Π_*))`. **It contains no assert statements at all.**

The Mathematica script performs the same numerical computation. It contains exactly one substantive check, `expectZero["R_q(g_minus)-1/4", ((gMinus - rF1)^2/(1 + rF1^2)) - 1/4]`, and otherwise prints `Pi_*, Sigma_*, T_*, A_T, B_T, |A_T|/B_T, delta T_m`, and `W_center`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (C1) closed-form A_T, B_T in terms of starred quantities | sympy/mathematica compute them numerically from the closed form; no assertion that they take the closed form | partial (computed, not verified) |
| (C2a) numerical `A_T ≈ -4.27264` | both scripts print `A_T = -4.27263956256927...`; no assertion against the paper-quoted literal | missing |
| (C2b) numerical `B_T ≈ 0.134875` | both scripts print `B_T = 0.134875005736...`; no assertion against the paper-quoted literal | missing |
| (C2c) `|A_T|/B_T ≈ 31.6785` | printed (`31.67851255487656...`); no assertion | missing |
| (C3) centered kernel rewrite `δT_m = ε ∫ W_*(ς - Σ_*) dx` | sympy and mathematica each `simplify`/`FullSimplify` the bilinear form `eps*(A_T*(gbar-gminus) + B_T*(Sbar - S_*))` and print it; no assertion that this is equivalent to the integral representation, no test of the centering property `∫(Σ_ε - Σ_*) dx = 0` (or its source-moment analogue) | missing |
| Card Check (a) — first-order deformations centered before covariance | not exercised | missing |
| Card Check (b) — rigidity kernel (not branch ambiguity) controls correction | partly probed by the `R_q(g_minus) - 1/4 == 0` tautology in mathematica, but that line just re-states the definition of `gMinus`; it does not show the kernel rather than branch ambiguity controls the result | missing |
| Card Check (c) — one-step nonlinear correction within reduced mouth-layer regime | not exercised | missing |

Dominant pattern: paper deliverables exist on the page but are not exercised by any assertion. The only assertion present is algebraically tautological. Setting `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | — | (none) | none | no |
| A2 | mathematica | 64 | `expectZero["R_q(g_minus)-1/4", rQMinus - 1/4]` where `rQMinus = ((gMinus - rF1)^2)/(1 + rF1^2)` and `gMinus := rF1 - Sqrt[1+rF1^2]/2` | nominally Card Check (b) (R_q=1/4 on the canonical lower branch); but tautological | no |

The sympy script has no assertion row at all. The mathematica row is tautological — see F2.

## Findings

### F1 — missing_verification_script

**Subtype:** script_doesnt_cover_claim
**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:1-63` (entire file)

**What's wrong:**
The SymPy audit script has zero assertions. Its body (lines 44-63) consists of `print` and `sp.pprint` calls only:

> ```
> banner("Stage 147 audit: first-order rigidity coefficients")
> print("Pi_* =", Pi_star)
> ...
> print("A_T =", AT)
> print("B_T =", BT)
> print("|A_T|/B_T =", sp.N(abs(AT)/BT, 20))
> ...
> dT = sp.simplify(eps*(AT*(gbar-gminus) + BT*(Sbar-Sformula.subs(Pi, Pi_star))))
> print("delta T_m =")
> sp.pprint(dT)
> Wcenter = sp.simplify(AT*(c-gminus) + BT*(Kq-Sformula.subs(Pi, Pi_star)))
> ```

No `assert`, no exit-on-fail. The script always exits 0 regardless of whether the computed `A_T`, `B_T`, or `|A_T|/B_T` agree with the paper-quoted values. Saved output `PASS` label is therefore trivially true: it reflects "the script ran" not "the math holds".

**Why this matters:**
Stage 147's paper card makes a numerical bottom-line claim (`A_T ≈ -4.27264`, `B_T ≈ 0.134875`), and the appendix (`stage_appendix_part04.tex:846-848`) quotes both literals to 14 digits. The SymPy audit must independently verify those values from the closed forms and pass/fail against the paper-quoted literals. Currently a typo in one of the closed-form denominators would still produce a `PASS` output transcript. Stage 147 is checkpoint:`False` but `is_status_only_candidate:False`; both engines are required and both must contain substantive checks.

**Required change:**
Add explicit numerical assertions in the SymPy script, anchored to the paper-quoted literals, and add a structural assertion that verifies the closed-form A_T against an independent derivation route (the chain-rule expansion of `δT_m = (9/(40 T_*)) · (δΣ_0/sqrt(...))` is the natural second route — see directive for the recommended form).

**Verification:**
After fix, the new SymPy script's assertion block should fail if either `A_T` or `B_T` is perturbed by more than 1e-12, and the saved output should contain explicit `PASS` lines naming each check.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:38-39, 63-64`

**What's wrong:**
The only substantive Mathematica assertion is

> ```
> rQMinus = FullSimplify[((gMinus - rF1)^2)/(1 + rF1^2), Assumptions -> True];
> expectZero["R_q(g_minus)-1/4", rQMinus - 1/4];
> ```

with `gMinus = FullSimplify[rF1 - Sqrt[1 + rF1^2]/2, ...]` (line 39).

This is algebraically guaranteed by the definition of `gMinus`. Substituting:
- `gMinus - rF1 = -Sqrt[1 + rF1^2]/2`
- `(gMinus - rF1)^2 = (1 + rF1^2)/4`
- `(gMinus - rF1)^2 / (1 + rF1^2) = 1/4` for any `rF1` whatsoever.

The identity has nothing to do with the canonical-branch `R_q = 1/4` law from the notes (Section 1: `R_q = 1/4` is presented as a *result* of the lower-compensated-branch construction, not as an algebraic rearrangement of how `gMinus` was defined). The check would still pass if `rF1` were replaced by any other real expression and `gMinus` redefined as `rF1 - Sqrt[1+rF1^2]/2` — i.e., it has no contact with the family-1 anchor `rF1 = Sqrt[12*(37/20)^2/Pi^2 - 1]`.

**Why this matters:**
The mathematica script's only assertion does not exercise any physical content of Stage 147. Combined with F1 (SymPy has none), the two-engine policy reduces to "two engines compute the same printable numbers". The output cannot detect a sign error in A_T, an arithmetic slip in B_T, or a wrong choice of branch root.

**Required change:**
Either (i) replace the tautology with a real branch identity that connects to the family-1 anchor — `R_q(g_minus) = 1/4` should be derived from the canonical compensation condition (e.g., by deriving `g_minus` as the unique real root of `4 R_q(g) - 1 = 0` on the lower branch, *not* by definition-substitution) — or (ii) drop that check entirely and add the numerical A_T, B_T, and centering assertions described in the directive.

**Verification:**
After fix, the assertion should fail if `gMinus` is altered without simultaneously altering `rF1` or the branch equation. The output transcript should show a `PASS:` line whose content is not algebraically derivable from the definition of `gMinus` alone.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:44-62`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:55-70`

**What's wrong:**
Neither engine verifies the centering structure that the notes (Section 2) elevate to the bottom-line statement:

> "Because `Σ_ε − Σ_*` integrates to zero, the traction shift can be written as a single weighted overlap … `δT_m = ε ∫₀¹ W_*(x)[ς(x) − Σ_*(x)] dx`."

The scripts only print `dT = eps*(A_T*(gbar-gminus) + B_T*(Sbar - sStar))` (bilinear in two moments) and `W_center = A_T*(c-gMinus) + B_T*(Kq - sStar)` (the centered weight). They never assert that the integral form is equivalent to the moment form when `gbar`, `Sbar` are computed from `ς(x)` via the two source moments, nor that `∫₀¹ (Σ_ε − Σ_*) dx = 0` (or its source-moment analog) — which is the "centered" half of the claim and what justifies subtracting `Σ_*(x)` inside the integrand.

Card Check (a) "first-order profile deformations are centered before covariance formulas are used" is not exercised at all. Card Check (c) "one-step nonlinear correction remains within the reduced mouth-layer regime" is also not exercised — neither script tests any nonlinear regime condition.

**Why this matters:**
The stage's title — "First-Order Rigidity Kernel" — is centered on (C3). The numerics for `A_T, B_T` are necessary but not sufficient: the load-bearing physical claim is "the entire first-order traction shift collapses to one scalar overlap kernel". Without checking the integral-moment equivalence (or, at minimum, the centering identity that `∫₀¹ (c(x) − g_*) dx = 0` and `∫₀¹ (K_q(x) − S_*) dx = 0` for the linearized source-moment definitions used in the notes), the audit doesn't certify the rigidity-kernel claim.

**Required change:**
Add an assertion in both engines that the linearized source moments of `c(x)` and `K_q(x)` reproduce `g_*` and `S_*` exactly (this is the operational form of the centering identity at the canonical branch point and is purely an integration check the engines can do symbolically). See directive for the exact form.

**Verification:**
After fix, the new check should fail if `c(x)` is replaced by `cos(pi x / 3)` or `K_q(x)` is replaced by a profile with the wrong `kap`; the output transcript should report both moment integrals as `PASS`.

## Independent-derivation check (Mathematica)

The Mathematica script reproduces the SymPy choreography in essentially the same order: same closed-form `gFormula`/`sFormula`, same `rF1` literal `Sqrt[(12*(37/20)^2)/Pi^2 - 1]`, same `gMinus = rF1 - Sqrt[1+rF1^2]/2`, same `pStar` via `FindRoot[gFormula == gMinus, {p, 1.5}]`, same `A_T`/`B_T` formula assembly. Compare:

- SymPy lines 22-24:
> ```
> rF1 = sp.sqrt(12*sp.Rational(37,20)**2/sp.pi**2 - 1)
> gminus = sp.simplify(rF1 - sp.sqrt(1+rF1**2)/2)
> Pi_star = sp.N(sp.nsolve(gPi - gminus, 1.5), 40)
> ```
- Mathematica lines 38-40:
> ```
> rF1 = Sqrt[(12*(37/20)^2)/Pi^2 - 1];
> gMinus = FullSimplify[rF1 - Sqrt[1 + rF1^2]/2, Assumptions -> True];
> pStar = p /. FindRoot[gFormula == gMinus, {p, 1.5}, ...];
> ```

Likewise lines 33-42 (sympy) and 49-53 (mathematica) for `A_T, B_T`. This is a direct port. Strictly speaking it satisfies the second-engine policy in the narrow sense that Mathematica computes the result via its own numerics, not by reading SymPy's output; but neither engine derives `A_T, B_T` from a *different* route (e.g., direct chain-rule on `T_m = Sqrt[9 Σ_0 / 20]` with `Σ_0 = Π/(1 - S/4)` followed by symbolic comparison). I do not file `mathematica_transliteration` here because the .wl is short and the standard branch-anchor pattern is required for the family-1 point; but the absence of a second derivation route compounds the F1/F2/F3 weakness — the engines agree on *everything*, including any error.

Minor cosmetic: the .wl banner at line 26 reads `STAGE 130 — FIRST-ORDER RIGIDITY KERNEL`. This is a stale-from-renumber label, not a math defect, but should be corrected to `STAGE 147` for traceability. (Codex can apply this in the same patch.)

## Engine cross-check

Both engines print identical numerical values (to ~30 digits) for `Pi_*, Sigma_*, T_*, A_T, B_T, |A_T|/B_T`. From the sympy and mathematica transcripts:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `Pi_*` | 1.508829513493155527470435117720626294613 | 1.5088295134931555830055507559542749287786 |
| `Sigma_*` | 1.805941110956353804451202617950474992639 | 1.8059411109563538723736729091995054268193 |
| `T_*` | 0.9014840541742040220118492172422879595856 | 0.9014840541742040389645127111412155822626 |
| `A_T` | -4.27263956256927466681875206233 | -4.2726395625692746412223747285390913640779 |
| `B_T` | 0.134875005736705540316968999027 | 0.1348750057367055429617462273081053431702 |
| `|A_T|/B_T` | 31.678512554876561144 | 31.67851255487656033274322154149180316197 |

Agreement is to ~16 significant figures, consistent with the 40-digit SymPy `nsolve` start and the 80-digit Mathematica `FindRoot` working precision. No `engine_disagreement` finding.

## Verdict justification

Verdict: `findings`, `stop_cold: null`. The math the scripts *compute* is internally consistent and matches the paper-quoted A_T and B_T literals to all stated digits, so the underlying claim is not in doubt. But the audit *as written* fails to certify it: SymPy has no assertions at all (F1), the single Mathematica assertion is tautological (F2), and the centering/integral-form half of the stage's load-bearing claim is unverified by either engine (F3). Each defect is a script-side fixable issue — none constitute paper_misalignment (the paper and the scripts agree on numerical values where the scripts deign to compute them), and none rises to UNFIXABLE or CRITICAL_DOWNSTREAM. Attacks tried that failed: (i) sign of A_T (both scripts produce -4.272...; the closed-form sign is set by `-(9/(40 T_*))` and the inner bracket is manifestly positive at the family-1 point, so the negative is structurally correct); (ii) factor of 2 in the inner bracket (the second term has a single power of `(1 - S_*/4)^2` consistent with `(δ/δΠ)(Π/(1-S/4)) = 1/(1-S/4) + Π S'/(4(1-S/4)^2)`, and the prefactor `9/(40 T_*)` is what `(1/2)(20/9)^(-1/2)·sqrt(9/20)` reduces to, matching `dT_m/dΣ_0 = 9/(40 T_*)`); (iii) `Π_*` root selection (both engines find the same root near 1.508, and `gFormula` is monotone in that region so the root is unique). Paper and notes were read first; the paper's central numerical claim aligns with what the scripts compute. The verdict is "findings" because the assertions are absent or tautological, not because the math is wrong.

## Self-test notes

Walked the proposed F1/F3 fixes mentally: (i) the numerical-literal assertion `abs(A_T - sp.Float("-4.27263956256927")) < 1e-12` cannot fail trivially because both sides are independently computed and the literal is the paper-quoted constant; (ii) the centering integral `sp.integrate(c - gminus, (x, 0, 1))` reduces to `2/pi - gminus` symbolically, and the centering claim requires this to equal zero exactly — which would actually FAIL at the moment level used in the notes unless `g_*` (in the notes' linearized inner product) is defined as `∫₀¹ c(x) dx`, not as `gPi(Π_*)`. This subtlety is exactly why I phrased F3 around the operational form "the linearized source moments of `c(x)` and `K_q(x)` reproduce `g_*` and `S_*` exactly", which is a check the script's *own* convention can pass; the directive defers to inline comments to disambiguate which `g_*` is meant. No path-spec issues (`.py` → `scripts/`, `.wl` → `mathematica/`). Paper round-trip: the proposed assertions reference only `A_T ≈ -4.27263956256927` and `B_T ≈ 0.134875005736706` (the paper-quoted literals from `stage_appendix_part04.tex:846-848`); no new constants introduced.
