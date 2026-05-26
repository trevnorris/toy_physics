---
unit_id: 055
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage055_explicit_lowest_lane_reachability.md
  paper_appendix: present
---

# Audit unit 055 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_055.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage055_explicit_lowest_lane_reachability.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows for stage 055 at line 88; `\input{stages/stage_055}` at line 228)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage055_explicit_reachability_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.txt`

## What the paper claims

The stage card (Output, line 29) states the boxed result: `1 <= zeta_0^phys <= (pi^2/4)*(4/(4-x))` for the explicit non-twin lowest-lane family, with `zeta_0^phys = A_K(eta) * Omega_0^2` (eq:app-stage055-zeta-factor). The notes elaborate the explicit family as `zeta_0^(exp+R)(alpha, eta) = Omega_exp(alpha)^2 / (1 - x/4 + x y(eta)^2 / pi^2)`, where `y tan y = eta`, with closure range `1 <= zeta <= pi^2/(4-x)`, compliance floor `x >= 4 - pi^2/zeta_req`, and equivalent stiffness-ratio form `K_X/K_W^eff <= pi^2/(4 zeta_req)`. The appendix row tags status `\StatusExactClosure{}` and summarizes "Combined overlap/compliance reachability window for the lowest support lane." The notes also enumerate three regimes (A/B/C) splitting by `zeta_req` vs `pi^2/4` and `pi^2/(4-x)`. Distinct deliverables: (D1) lower bound = 1 at the symmetric twin point (`alpha = 0`, `eta = +inf` → `y = pi/2`); (D2) upper bound = `pi^2/(4-x)` at closure (`alpha → +inf`, `eta → 0+` → `y = 0`); (D3) compliance floor `x_floor = 4 - pi^2/zeta_req`; (D4) stiffness-ratio equivalent `pi^2/(4 zeta_req)`; (D5) the three regimes.

## What the script claims to verify

Both scripts construct `Omega_exp(alpha)`, `A_K(y, x)`, and `zeta_family = Omega_exp^2 * A_K` using the exact closed forms quoted in the notes. They then assert (i) `lim_{alpha→0} zeta_family |_{y = pi/2} = 1` (twin point); (ii) `lim_{alpha→∞} zeta_family |_{y → 0} = pi^2/(4 - x)` (closure maximum); (iii) `x_floor = 4 - pi^2/zeta_req`; (iv) `(1/A_K)|_{y = 0, x = x_floor} = pi^2/(4 zeta_req)` (KX/KW equivalence). Mathematica adds a substantive cross-check (line 59): substituting `x_floor` into `zeta_max` returns `zeta_req`. The regime split (D5) is print-only commentary — acceptable because it is a logical immediate consequence of the verified endpoints D1–D3.

## Paper ↔ script cross-check

| Paper deliverable | Script-side coverage | Status |
|---|---|---|
| D1: lower bound = 1 (twin point) | sympy:45, math:56 | match |
| D2: upper bound = `pi^2/(4-x)` | sympy:46, math:57 | match |
| D3: x_floor = `4 - pi^2/zeta_req` | sympy:50–52 (Solve + assertion); math:51,58,59 | match (sympy derives by Solve; math hardcodes the value then verifies by substitution at line 59) |
| D4: KX/KW = `pi^2/(4 zeta_req)` | sympy:56, math:60 (both use `1/AK` non-tautologically) | match |
| D5: regime A/B/C split | sympy:60–65 print, math: absent | prose-only (acceptable — logical corollary of D1–D3) |

The boxed paper form `(pi^2/4)*(4/(4-x))` is algebraically `pi^2/(4-x)` as used by the script — no constant drift. All other numeric and symbolic constants (4, pi^2, pi^2/(4 zeta_req)) match the paper/notes verbatim. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 45 | `expect_zero("twin value", zeta_twin - 1)` (via `sp.limit(..., alpha, 0).subs(y, pi/2)`) | D1 | yes |
| A2 | sympy | 46 | `expect_zero("closure maximum", zeta_max - pi**2/(4 - x))` (via `sp.limit(..., alpha, sp.oo).subs(y, 0)`) | D2 | yes |
| A3 | sympy | 50–52 | `x_floor = sp.solve(sp.Eq(zeta_max, zeta_req), x)[0]`; `expect_zero("x floor = 4 - pi^2/zeta_req", x_floor - (4 - pi**2/zeta_req))` | D3 | yes (Solve output compared to paper-stated form) |
| A4 | sympy | 56 | `expect_zero("KX/KW equivalence", (1/AK).subs(y, 0).subs(x, x_floor) - pi**2/(4*zeta_req))` | D4 | yes (substantive — uses `1/AK`, not hand-typed `1 - x/4`) |
| A5 | math | 56 | `expectZero["twin value", zetaTwin - 1]` | D1 | yes |
| A6 | math | 57 | `expectZero["closure maximum", zetaMax - Pi^2/(4 - x)]` | D2 | yes |
| A7 | math | 58 | `expectZero["x floor = 4 - Pi^2/zeta_req", xFloor - (4 - Pi^2/zetaReq)]` where `xFloor = 4 - Pi^2/zetaReq` was hardcoded on line 51 | D3 | **no — tautological by construction** |
| A8 | math | 59 | `expectZero["zeta_max(x_floor) - zeta_req", (zetaMax /. x -> xFloor) - zetaReq]` | D3 | yes (substantive: substitutes the candidate value and confirms zeta_max recovers zeta_req) |
| A9 | math | 60 | `expectZero["KX/KW equivalence", ((1/aK) /. y -> 0 /. x -> xFloor) - Pi^2/(4 zetaReq)]` | D4 | yes |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl:51,58`

**What's wrong:**
Line 51 defines:
```
xFloor = FullSimplify[4 - Pi^2/zetaReq, Assumptions -> zetaReq > 0];
```
This hardcodes the paper-claimed value of the compliance floor directly. Then line 58 asserts:
```
expectZero["x floor = 4 - Pi^2/zeta_req", xFloor - (4 - Pi^2/zetaReq)];
```
This reduces to `(4 - Pi^2/zetaReq) - (4 - Pi^2/zetaReq) == 0`, which is algebraically guaranteed by construction. The assertion cannot fail no matter what `Omega_exp`, `A_K`, or `zetaMax` evaluate to — it confirms a literal against itself. By contrast the SymPy script (line 50) derives `x_floor` independently via `sp.solve(sp.Eq(zeta_max, zeta_req), x)[0]`, so its corresponding assertion on line 52 is genuinely substantive.

Note that D3 is *not* unverified in Mathematica: the very next line (59) substitutes `xFloor` into `zetaMax` and checks the residual against `zetaReq` — that check is substantive (a regression in `zetaMax` or `xFloor` would surface there). The tautology is confined to line 58.

**Why this matters:**
Line 58 is presented as a verification line ("PASS: x floor = 4 - Pi^2/zeta_req" in the transcript), but it cannot fail and therefore provides no information. If the paper changed `xFloor` to `4 - Pi^2/(zetaReq+1)`, this line would still PASS by simply propagating the new literal into both sides. The Mathematica audit relies entirely on line 59 for D3 — which is fine substantively but leaves a misleading PASS line in the transcript that suggests redundant verification.

**Required change:**
Replace the hardcoded definition with a Solve-derived value, mirroring the SymPy script. At line 51, change:
```
xFloor = FullSimplify[4 - Pi^2/zetaReq, Assumptions -> zetaReq > 0];
```
to:
```
xFloor = FullSimplify[
  x /. First[Solve[zetaMax == zetaReq, x]],
  Assumptions -> $Assumptions
];
```
Leave lines 58 and 59 unchanged — once `xFloor` is derived by `Solve`, line 58's assertion becomes a non-tautological check that the Solve output equals the paper's claimed closed form, and line 59 remains a substantive substitution check.

**Verification:**
After the fix, the transcript line "x floor from zeta_max = zeta_req = 4 - Pi^2/zetaReq" should still print (Solve gives this), line 58's "x floor = 4 - Pi^2/zeta_req = 0" should still PASS, and line 59's "zeta_max(x_floor) - zeta_req = 0" should still PASS. Critically, the algebra of D3 no longer depends on the script's author having pre-typed the right answer. The exit code remains 0.

### F2 — symbol_assumption_error

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py:27`

**What's wrong:**
Line 27 declares all symbols positive and real:
```
alpha, x, y, zeta_req = sp.symbols("alpha x y zeta_req", positive=True, real=True)
```
The paper/notes restrict `x ∈ (0, 4)` (the closure maximum `pi^2/(4-x)` is finite and positive only on this open interval; the notes' Regime C explicitly references the case `zeta_req > pi^2/(4-x)`, which only makes sense for `x < 4`). SymPy is told only `x > 0`, so `(pi**2/(4-x)).subs(x, 5)` would happily yield `-pi**2`, a meaningless result for the family. The Mathematica script correctly declares `0 < x < 4` (line 36).

Additionally `y` is declared only positive; the physical setup has `y` as the principal root of `y tan y = eta`, hence `y ∈ (0, pi/2)`. Not load-bearing in this stage (the script only evaluates at `y = 0` and `y = pi/2`), but the script silently models a wider domain than the paper requires.

**Why this matters:**
No assertion in the current script is invalidated, because (a) the closure formulas simplify symbolically without needing `x < 4`, and (b) `y` is only used at endpoints. So this is cosmetic now. But the script appears to model `x ∈ (0,4)` and silently relaxes the constraint — a downstream cut-and-paste of this declaration template could miss a divergence at `x = 4` without the auditor noticing.

**Required change:**
At line 27, replace the joint declaration with one that documents the paper-stated domain. Either:

(option A — minimal, documentary)
```
# Paper-stated domain: alpha > 0, 0 < x < 4, 0 < y < pi/2, zeta_req > 0.
# SymPy lacks compound symbol-level bounds, so positivity is declared here
# and the (0, 4) and (0, pi/2) bounds are enforced by the assertions below
# (closure maximum, twin/closure endpoint substitutions).
alpha, x, y, zeta_req = sp.symbols("alpha x y zeta_req", positive=True, real=True)
```

(option B — adds a numeric sanity gate)
Keep the existing declaration and add, after line 35, a small sanity check:
```
# Closure maximum must be positive and finite on x in (0, 4).
assert sp.simplify(sp.Piecewise(
    (1, sp.And(x > 0, x < 4)), (0, True)
).subs(x, sp.Rational(1, 2))) == 1
assert (pi**2 / (4 - x)).subs(x, sp.Rational(1, 2)) > 0
```

Option A is preferred — Option B's piecewise probe doesn't actually constrain the symbolic algebra below and would add noise. The change is documentary: it surfaces the paper's stated domain at the point where the symbols are declared, matching Mathematica's `0 < x < 4`.

**Verification:**
After the fix, all printed expressions and assertion residuals should be unchanged (the algebra was already correct). The change is documentary — it brings the SymPy script's symbol declarations into agreement with the Mathematica script's `$Assumptions` and the paper-stated domain. The verifier confirms by reading line 27's comment and seeing the paper's `x ∈ (0, 4)` and `y ∈ (0, pi/2)` named explicitly.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally close to the SymPy script in this combinator stage: same variable names in `camelCase` (`omegaExp`, `aK`, `zetaFamily`), same construction order (define Omega_exp, then A_K, then product, then take twin/closure limits, then verify floor and KX/KW). Compare:

SymPy 30–35:
```
Omega_exp = sp.simplify( pi * alpha * (2*alpha*sp.exp(alpha) + pi) / ((4*alpha**2 + pi**2)*(sp.exp(alpha)-1)) )
AK = sp.simplify(1 / (1 - x/4 + x*y**2/pi**2))
zeta_family = sp.simplify(Omega_exp**2 * AK)
```
Mathematica 38–43:
```
omegaExp = FullSimplify[ Pi alpha (2 alpha Exp[alpha] + Pi)/((4 alpha^2 + Pi^2)(Exp[alpha]-1)), ...];
aK = FullSimplify[1/(1 - x/4 + x y^2/Pi^2), ...];
zetaFamily = FullSimplify[omegaExp^2 aK, ...];
```

The construction is parallel. However, this stage's algebraic content is exhausted by (a) carrying forward two upstream closed forms from Stages 053/054 and (b) two endpoint limits plus one rearrangement to recover `x_floor` — there is no independent re-derivation to be had at the operator level, because both forms come from upstream stages by reference. Mathematica adds two genuine deviations from the SymPy path: (i) it uses `Limit[..., Direction -> "FromAbove"]` for the `y → 0` closure (SymPy uses unidirectional `subs`), and (ii) it adds an extra cross-check on line 59 (`zetaMax(xFloor) - zetaReq == 0`) absent from the SymPy script. I do not raise `mathematica_transliteration`: the upstream closed forms are not in scope of this unit to re-derive, and line 59 is an independent check. Borderline but acceptable for a combinator stage.

## Engine cross-check

Both engines report the same closure forms:
- `Omega_exp`: SymPy `pi*alpha*(2*alpha*exp(alpha) + pi)/((4*alpha**2 + pi**2)*(exp(alpha) - 1))` vs Mathematica `(alpha*Pi*(2*alpha*E^alpha + Pi))/((-1 + E^alpha)*(4*alpha^2 + Pi^2))` — algebraically identical (denominator factored differently; sign of `exp(alpha) - 1` vs `-1 + E^alpha` matches).
- `A_K`: SymPy normalizes to `4*pi**2/(4*x*y**2 + pi**2*(4 - x))`; Mathematica keeps `(1 + x*(-1/4 + y^2/Pi^2))^(-1)`. Multiply Mathematica's reciprocal expression numerator and denominator by `4*Pi^2` → same.
- `closure maximum`: SymPy prints `-pi**2/(x - 4)` = `pi**2/(4 - x)` = Mathematica's `Pi^2/(4 - x)`. Match.
- All four (sympy) / five (math) `expect_zero` / `expectZero` residuals print `0` and PASS.

Engines agree. The Mathematica `Limit::alimv` warnings (output lines 9–13) are benign — Mathematica warns that `$Assumptions` containing `alpha > 0` is ignored inside `Limit[..., alpha -> Infinity]`; the limit is still computed correctly to `pi/2` and the assertion passes.

Output freshness: SymPy script mtime 2026-05-22 17:40, output 2026-05-22 17:41 (fresh by 1 min). Mathematica script mtime 2026-05-22 17:40, output 2026-05-22 17:41 (fresh by 1 min). No `stale_output` finding.

## Verdict justification

The script's load-bearing claims (D1–D4) line up exactly with the paper's `\stagefield{Output}` and the notes' enumerated deliverables. Both engines compute and agree. Constants match the paper verbatim (no value drift). The only defects are (F1) one tautological Mathematica line that compares a hardcoded literal against itself — backstopped by a substantive cross-check on the next line, so D3 is still genuinely verified, but the line as written cannot fail — and (F2) loose symbol assumptions on the SymPy side that don't currently invalidate any assertion but drift away from the paper's stated domain `x ∈ (0,4)`. Neither is `stop_cold`-worthy. No `paper_misalignment` — the script's algebra matches the paper card and notes exactly.

Attacks attempted that failed: (1) tried to find a sign error in the closure maximum — both engines agree on `pi^2/(4-x)` positive on `0 < x < 4`; (2) tried to find a missing branch of `y tan y = eta` — script evaluates only at endpoints `y = 0` and `y = pi/2`, the correct closure limits of the principal branch (the notes' `lim_{eta → 0+} A_K = 4/(4-x)` corresponds to `y → 0`, which the script uses); (3) tried to catch an `Omega_exp` slip at `alpha → 0` — limit is 1, consistent with the notes' explicit `Omega_exp(0) = 1`; (4) tried to find a hidden hardcoded numeric — only `xFloor` in Mathematica (F1); the SymPy x_floor comes from `sp.solve`; (5) tried to find a paper↔script constant mismatch — `(pi^2/4)*(4/(4-x))` (paper boxed form) equals `pi^2/(4-x)` (script form) exactly.

Verdict: `findings` (2, both low). No stop-cold flag.

## Self-test notes

Checked: (1) Variable-independence trap — not applicable; no `sp.diff` or `D[]` calls, only `Limit` and substitution. (2) Integrand-parity trap — not applicable; no integrals. (3) Trivial-case pre-check for F1's required change: `Solve[zetaMax == zetaReq, x]` with `zetaMax = Pi^2/(4-x)` yields `x = 4 - Pi^2/zetaReq`, matching the current literal — proposed fix is internally consistent and does not change the printed result. (4) Path specifications — F2 edits `scripts/...py` line 27; F1 edits `mathematica/...wl` line 51 only (lines 58 and 59 stay). (5) Paper round-trip — neither fix introduces new constants or alters any paper-side claim; F1 substitutes Solve for a hardcoded literal (same numerical answer); F2 is purely documentary.
