---
unit_id: 134
batch: IV.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md
  paper_appendix: present
---

# Audit unit 134 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_134.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows around lines 660-725 and the `\input{stages/stage_134}` at line 1302)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.txt`

## What the paper claims

Stage 134 is the Family-1 shell + first mixed D/N tube reduction of the coupled mouth-layer fixed-point law. The stage card's purpose line says it is a "coupled mouth fixed point and gain selection ledger step" and quotes the bottom-line equation: "First explicit reduction gives `Pi = M_s + M_q * S_q(Pi)`." The notes give the explicit closed form `S_q(Pi) = Pi*[(pi/2)*tanh(pi/2) + Pi*(e^-Pi * sech(pi/2) - 1)] / [(1 - e^-Pi)*(pi^2/4 - Pi^2)]`, the shell-channel limit `kappa_s = 0 => S = 1`, the mixed lane `kappa_q = pi/2`, the numerical evaluation `S_q(Pi_*) ≈ 0.658075937605428` at `Pi_* ≈ 1.50882951349316`, and the canonical gain line `M_s ≈ 1.50882951349316 - 0.658075937605428 * M_q`. The card's `Checks` list also asks for checks against (i) outlet consistency, (ii) self-matched susceptibility closure, and (iii) numerical-only fixed points.

## What the script claims to verify

The SymPy script defines `S(Pi, kappa)` from the kernel formula, prints `S_shell = limit(S, kappa->0) = 1`, prints `S_q(Pi) = S(Pi, pi/2)`, prints the fixed-point law `Pi_eq = M_s + M_q * S_q`, and prints the numerical `S_q(Pi_*)` and gain line at the hardcoded `Pi_* = 1.50882951349316`. The SymPy script contains **zero assertions** — it only prints. The Mathematica script does the analogous flow but adds two assertions: `expectZero["static shell channel", sShell - 1]` and `expectZero["specialized D/N kernel", sQ - sQExpected]`. Both scripts print a banner "STAGE 117 — FAMILY-1 FIXED-POINT REDUCTION" (wrong stage number).

## Paper ↔ script cross-check

| paper deliverable | script-side coverage | status |
|---|---|---|
| `Pi = M_s + M_q*S_q(Pi)` (fixed-point law) | both scripts print the symbolic `Pi_eq = Ms + Mq*S_q`; no assertion comparing to a target form | partial (printed, not asserted) |
| Closed form `S_q(Pi) = Pi*[(pi/2)tanh(pi/2) + Pi*(e^-Pi*sech(pi/2)-1)] / [(1-e^-Pi)*(pi^2/4-Pi^2)]` | sympy prints `S(Pi, pi/2)` (no comparison); mathematica asserts `sQ - sQExpected == 0` where `sQExpected` is the same kernel formula re-typed at k=pi/2 | mismatch (sympy missing); tautological in mathematica |
| Static shell limit `S(Pi,0)=1` | sympy prints `S_shell = 1` only; mathematica asserts `sShell - 1 == 0` | partial (mathematica only) |
| `S_q(Pi_*) ≈ 0.658075937605428` | both scripts compute and print, but no numerical-tolerance assertion comparing to the expected literal | partial (printed, not asserted) |
| Gain line `M_s ≈ 1.50882951349316 - 0.658075937605428*M_q` | both scripts print the line `Pi_star - sStar*M_q`; no assertion comparing to expected coefficients | partial (printed, not asserted) |
| Check gain pair (Ms, Mq) against outlet consistency (per `\stagefield{Checks}`) | not exercised by either script | missing |
| Check self-matched susceptibility closure (per `\stagefield{Checks}`) | not exercised by either script | missing |
| Numerical fixed points recorded numerically (per `\stagefield{Checks}`) | `Pi_*` is hardcoded as a float literal `1.50882951349316`, sourced from stage 233 per notes; no re-derivation here, which is consistent with "numerically located, not closed-form" | partial (uses hardcoded `Pi_*` literal without provenance comment) |

Dominant pattern: most deliverables are printed but not asserted; the only non-trivial assertion (mathematica's `sShell - 1`) covers one auxiliary claim. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | (none) | none | n/a | n/a |
| A2 | mathematica | 48 | `expectZero["static shell channel", sShell - 1]` | shell limit `S(Pi,0)=1` | yes |
| A3 | mathematica | 51 | `expectZero["specialized D/N kernel", sQ - sQExpected]` | closed form `S_q(Pi)` | no — tautological; `sQExpected` is the same formula re-typed at `k=pi/2` |

## Findings

### F1 — missing_verification_script

**Severity:** high
**Subtype:** `script_doesnt_cover_claim` (sympy side)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:27-43`

**What's wrong:**
The SymPy audit script for stage 134 contains zero assertions. It only prints values:

```
print("S_shell =", S_shell)
print("S_q(Pi) =")
sp.pprint(S_q)
print("Fixed-point law Pi =")
sp.pprint(Pi_eq)
...
print("S_q(Pi_star) =", S_star)
print("Canonical gain line: M_s = Pi_star - S_q(Pi_star) M_q")
sp.pprint(sp.N(Pi_star - S_star*M_q, 30))
```

There is no `assert`, no `==` check, no `expect_close`. The script's exit code 0 only proves it ran without raising; it does not certify any of the printed expressions. The stage card's `\stagefield{Verification}` cites this file as an audit, but it does not audit anything.

**Why this matters:**
The bottom-line claim of the stage (the fixed-point law `Pi = M_s + M_q*S_q(Pi)` with the specific closed-form `S_q`, its static-shell limit `=1`, its value at `Pi_*`, and the resulting canonical gain line) is not protected by any executable check on the sympy side. A future edit could silently break any of these and the script would still "PASS."

**Required change:**
Add concrete assertions in the SymPy script covering every paper-side deliverable. Specifically:
- `assert sp.simplify(S_shell - 1) == 0` (verify the kappa→0 limit equals 1).
- `assert sp.simplify(S_q - S_q_paper_form) == 0` where `S_q_paper_form` is constructed **as an independent re-derivation** (e.g., by solving the boundary-value problem for `kappa=pi/2`, or by deriving `S_q` from the diagonalized 1-D D/N problem), not by re-typing the same kernel literal at `k=pi/2`. If an independent derivation is out of scope for this stage, instead assert a non-tautological invariant of the closed form: e.g., evaluate `S_q` at three independent numeric values of `Pi` (say `Pi=0.5, 1.0, 2.0`) against high-precision numeric values stored as literals in the assertion.
- `assert abs(float(S_q.subs(Pi, Pi_star)) - 0.658075937605428) < 1e-12` to lock the numerical `S_q(Pi_*)` against the paper's quoted value.
- `assert abs(float(Pi_star - 0.658075937605428*M_q - (1.50882951349316 - 0.658075937605428*M_q)).subs(M_q, 0)) < 1e-12` (or simpler: assert the two coefficients of the gain line match the paper's quoted constants to 1e-12).

**Verification:**
After fix, the SymPy script contains at least three `assert` statements with non-tautological content, and `redteam exec-sympy 134` produces output showing each assertion's check explicitly (e.g., printed `OK: S_shell limit = 1`, `OK: S_q(Pi_*) ≈ 0.658075937605428`, etc.).

### F2 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:40-51`

**What's wrong:**
The Mathematica script asserts `expectZero["specialized D/N kernel", sQ - sQExpected]` where

```
sQ = FullSimplify[sKernel[p, Pi/2], Assumptions -> p > 0];
sQExpected = FullSimplify[
  p*((Pi/2)*Tanh[Pi/2] + p*(Exp[-p]/Cosh[Pi/2] - 1))/((1 - Exp[-p])*(Pi^2/4 - p^2)),
  Assumptions -> p > 0];
```

But `sKernel[p, k] := FullSimplify[p*(k*Tanh[k] + p*(Exp[-p]/Cosh[k] - 1))/((1 - Exp[-p])*(k^2 - p^2)), …]` is literally the same closed form. Evaluating `sKernel[p, Pi/2]` substitutes `k -> Pi/2` everywhere, producing an expression identical (after simplification) to `sQExpected`. The assertion is algebraically guaranteed to pass and tests only that `FullSimplify` is deterministic; it does not test that either form actually matches the paper's quoted closed form for `S_q`.

**Why this matters:**
The paper's central closed form for `S_q(Pi)` (notes Section 1, appendix eq. `app-part04-S-kernel` and `app-part04-F1-mouth-fixedpoint`) is the load-bearing identity of this stage. The mathematica script currently appears to "verify" it but the check is tautological.

**Required change:**
Replace the `sQExpected` definition with an **independent** construction of `S_q` — e.g., by solving the 1-D mouth boundary-value problem `[-d^2/dx^2 + kappa^2] u(x) = Sigma_Pi(x)` with `u(0)=0`, `u'(1)=0`, at `kappa=pi/2`, and reading off the D/N response. If that is out of scope here (this stage carries forward the kernel from earlier stages), then at minimum exercise the closed form numerically against three independent target values, e.g.:

```
expectClose["S_q at p=0.5", sQ /. p -> N[1/2, 30], <high-precision literal target>, 1*^-15];
expectClose["S_q at p=1.0", sQ /. p -> N[1, 30], <literal>, 1*^-15];
expectClose["S_q at p=2.0", sQ /. p -> N[2, 30], <literal>, 1*^-15];
```

The literal targets should be computed once by hand (or with arbitrary-precision arithmetic in a separate session) and pasted in, not derived from `sKernel` itself.

**Verification:**
After fix, the assertion `sQ - sQExpected == 0` is replaced by checks whose target values come from outside `sKernel`, so the check could in principle fail if `sKernel` had a typo. The script transcript must show each numeric check passing against an independently sourced literal.

### F3 — paper_misalignment

**Severity:** medium
**Subtype:** `script_missing_paper_claim`
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_134.tex:21-25` (the `\stagefield{Checks}` checklist)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:27-43` (entire body)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:26-67` (entire body)

**What's wrong:**
The paper card's `\stagefield{Checks}` block enumerates three checks the stage is supposed to perform:

```
\item Check the gain pair (M_s, M_q) against outlet consistency.
\item Check the self-matched susceptibility closure before using the one-scalar branch law.
\item Check numerical fixed points are recorded as numerically located, not closed-form constants.
```

Neither script exercises items 1 or 2. Item 1 (outlet consistency) corresponds to appendix eq. `app-part04-outlet-consistent-gains`: `M_s = 4*Sigma_m`, `M_q = -Sigma_m`, with `Sigma_m^* ≈ 0.451485277739090` (appendix eq. `app-part04-Sigmam-star`). The script reports the gain line `M_s ≈ 1.50882951349316 - 0.658075937605428*M_q` but never verifies that this is consistent with `(M_s, M_q) = (4*Sigma_m, -Sigma_m)` for the canonical `Sigma_m^*`. Item 2 (self-matched susceptibility closure) is not even printed.

Item 3 is implicitly honored (the scripts hardcode `Pi_* = 1.50882951349316` as a float literal sourced from stage 233, rather than re-deriving it as a closed form), but no comment in the script explicitly marks this as a numerical carry-forward.

**Why this matters:**
The paper card has a verification checklist that the audit scripts do not honor. Either the paper card overstates what stage 134 is supposed to verify (in which case the card needs trimming), or the scripts are incomplete (in which case checks for outlet consistency and self-matched susceptibility closure must be added). Direction of fix is a user call.

## Resolve before fix_loop

The paper card's checklist items "Check the gain pair (M_s, M_q) against outlet consistency" and "Check the self-matched susceptibility closure before using the one-scalar branch law" are not exercised by either script. Should these checks be added to the stage 134 scripts (with the appendix-quoted gain relations `M_s = 4*Sigma_m`, `M_q = -Sigma_m` and `Sigma_m^* ≈ 0.451485277739090`), or should they be moved to the appropriate downstream stage cards (e.g., 135 or later) where the outlet-consistent reduction and susceptibility closure actually live?

Possible directions:
- (a) Both checks belong to stage 134 → add `assert M_s_canonical == 4*Sigma_m_star`, `assert M_q_canonical == -Sigma_m_star`, and a check that `Sigma_m^* * (4 - S_q(Pi_*)) == Pi_*` to both scripts; cite source equations in comments.
- (b) These checks belong to a later stage → remove items 1 and 2 from `paper/stages/stage_134.tex:22-24`, leaving only the numerical-fixed-point note; no script change.
- (c) Some checks belong here and some elsewhere → user specifies the split.

The orchestrator will not invoke Codex on this misalignment until the user has chosen a direction.

### F4 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:21-43`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:31-61`

**What's wrong:**
The two scripts follow the same line-by-line choreography from the paper's closed-form kernel: define `S(Pi,k)` / `sKernel[p,k]` with the same algebraic form; compute `S_shell` as the `k->0` limit; specialize to `k = pi/2` to get `S_q`; define `Pi_eq = M_s + M_q*S_q`; substitute the hardcoded `Pi_* = 1.50882951349316`; print the gain line. Function-by-function:

SymPy lines 21-25:
```
def S(Pi, kappa):
    return sp.simplify(
        Pi * (kappa * sp.tanh(kappa) + Pi * (sp.exp(-Pi) / sp.cosh(kappa) - 1))
        / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2))
    )
```

Mathematica lines 31-34:
```
sKernel[p_, k_] := FullSimplify[
  p*(k*Tanh[k] + p*(Exp[-p]/Cosh[k] - 1))/((1 - Exp[-p])*(k^2 - p^2)),
  Assumptions -> $Assumptions
];
```

Same kernel, character-for-character (modulo syntax). The same is true for the `S_shell = limit(..., k->0)` block (sympy line 29 ↔ math lines 36-39) and the `S_q(Pi_*)` evaluation block (sympy lines 39-43 ↔ math lines 55-61). Neither engine re-derives the kernel from the boundary-value problem.

**Why this matters:**
The two-engine policy requires both engines to derive the result independently. If the kernel formula in both scripts was copied from the same source, an error in the closed form (e.g., a sign on `Tanh` or a factor on `Cosh`) would be invisible to the cross-check because both engines would carry the same error.

**Required change:**
At minimum, have one of the two engines (recommend Mathematica, since it handles symbolic ODEs more crisply) derive `S_q` from the 1-D boundary-value problem `[-d^2/dx^2 + (pi/2)^2] u(x) = Sigma_Pi(x); u(0)=0; u'(1)=0` and then `S_q = u'(0) / Sigma_Pi`-norm, rather than typing the closed form. If a full independent derivation is out of scope, then at least cross-check the two engines' `S_q(Pi)` at a handful of numeric `Pi` values: have sympy compute `S_q(0.5), S_q(1.0), S_q(2.0)` and write them to a small text file; have mathematica read that file and assert agreement to 1e-15. That at least verifies the two transliterations of the formula agree across engines rather than across re-typings.

**Verification:**
After fix, the two scripts' derivations of `S_q` no longer share a verbatim algebraic core, or there is an explicit cross-engine numeric comparison in one of them.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration: same kernel formula, same `k->0` limit, same `k = pi/2` specialization, same `Pi_*` hardcoded literal, same gain-line computation. It adds two `expectZero` checks (one substantive — shell limit; one tautological — closed-form re-typing) on top of the same algebraic skeleton. See F4.

## Engine cross-check

Both engines report the same `S_shell = 1`, the same closed-form `S_q` (modulo presentation), the same `S_q(Pi_*) = 0.658075937605428...` to 28 digits, and the same gain line `1.50882951349316 - 0.658075937605428 * M_q`. The numbers agree, but since the closed form was typed identically into both engines (F4), agreement is not informative about correctness of the closed form itself.

## Verdict justification

Verdict: `findings`. The shell limit `S(Pi,0)=1` is genuinely verified in Mathematica (and not in SymPy). Everything else — the closed form of `S_q`, the fixed-point law, the numerical value `S_q(Pi_*)`, and the canonical gain line — is printed but not asserted on either side, except for one tautological mathematica check. The sympy script has no assertions at all. The paper's `\stagefield{Checks}` checklist contains items not exercised by either script. Outputs are fresh (script mtime May 11 11:56, sympy output May 11 12:45, mathematica output May 11 13:12), so this is not a `stale_output` issue — it's a substantive verification gap.

Attacks tried that the scripts withstood:
- Numerical check of `S_q(Pi_*)`: both engines independently evaluate to 0.658075937605428 to 15+ digits → mathematically consistent but doesn't verify the closed form because both share the same typed formula.
- Stage banner says "STAGE 117" (typo) — cosmetic, not a math finding, so not filed.

`stop_cold` not set: none of the findings render the math unfixable, and none mathematically propagate downstream in a way that requires halting. The paper_misalignment is a scope question, not a math error.

## Self-test notes

Trap checks walked through: (1) the proposed sympy assertion `simplify(S_shell - 1) == 0` does fail-test (the `S_shell = limit(...)` is a non-trivial computation, and substituting a perturbed kernel would break it); (2) the proposed numeric checks `S_q(Pi_*) ≈ 0.658…` carry independent target literals, not derived from `sKernel` itself, so they non-tautologically exercise the closed form; (3) the F3 paper_misalignment is correctly routed as a user-resolution question rather than a Codex auto-fix, since the appendix eqs. (`Sigma_m^* ≈ 0.45148…`, `M_s = 4*Sigma_m`, `M_q = -Sigma_m`) could plausibly belong to stage 134 or to a downstream stage, and only the user knows the intended scope split.
