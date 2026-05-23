---
unit_id: 055
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 055 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage055_explicit_reachability_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.txt`

## What the script claims to verify

The scripts check the explicit "exp + Robin" family for the lowest reachability lane in the moving-throat PDE. They build `Omega_exp(alpha) = pi*alpha*(2*alpha*e^alpha + pi)/((4*alpha^2 + pi^2)*(e^alpha - 1))` and a shape factor `A_K(y,x) = 1/(1 - x/4 + x*y^2/pi^2)`, form `zeta_family = Omega_exp^2 * A_K`, then check (i) the symmetric "twin" limit (`alpha -> 0`, `y = pi/2`) gives `zeta = 1`, (ii) the closure-maximum limit (`alpha -> oo`, `y = 0`) gives `zeta = pi^2/(4 - x)`, (iii) the compliance floor solving `zeta_max = zeta_req` for `x` yields `x_floor = 4 - pi^2/zeta_req`, and (iv) an equivalent stiffness-ratio identity `KX/KW = 1 - x/4 = pi^2/(4*zeta_req)` at the floor. The narrative "regime split" (overlap-only vs combined overlap+softening vs infeasible) is print-only — no assertion attached.

## Assertion inventory

| #  | Script      | Line | Form                                                                  | Anchored to claim? |
|----|-------------|------|-----------------------------------------------------------------------|--------------------|
| A1 | sympy       | 45   | `simplify(zeta_twin - 1) == 0`                                        | yes                |
| A2 | sympy       | 46   | `simplify(zeta_max - pi**2/(4 - x)) == 0`                             | yes                |
| A3 | sympy       | 52   | `simplify(x_floor - (4 - pi**2/zeta_req)) == 0`                       | yes                |
| A4 | sympy       | 56   | `simplify((1 - x/4).subs(x, x_floor) - pi**2/(4*zeta_req)) == 0`      | no (tautological)  |
| M1 | mathematica | 50   | `FullSimplify[zetaTwin - 1] == 0`                                     | yes                |
| M2 | mathematica | 51   | `FullSimplify[zetaMax - Pi^2/(4 - x)] == 0`                           | yes                |
| M3 | mathematica | 52   | `FullSimplify[xFloor - (4 - Pi^2/zetaReq)] == 0`                      | yes (definitional, but harmless given xFloor is hardcoded) |
| M4 | mathematica | 53   | `FullSimplify[(zetaMax /. x -> xFloor) - zetaReq] == 0`               | yes                |
| M5 | mathematica | 54   | `FullSimplify[((1 - x/4) /. x -> xFloor) - Pi^2/(4 zetaReq)] == 0`    | no (tautological)  |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py:56`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl:54`

**What's wrong:**
Both engines define `x_floor = 4 - pi^2/zeta_req` (SymPy via `solve`, then asserted equal to that closed form on line 52; Mathematica hardcoded on line 45 and re-asserted on line 52). They then assert the "KX/KW equivalence":

SymPy line 56:
```
expect_zero("KX/KW equivalence", (1 - x / 4).subs(x, x_floor) - pi**2 / (4 * zeta_req))
```

Mathematica line 54:
```
expectZero["KX/KW equivalence", ((1 - x/4) /. x -> xFloor) - Pi^2/(4 zetaReq)];
```

Once `x_floor = 4 - pi^2/zeta_req` is established (assertions A3 / M3), substituting it into `1 - x/4` yields `1 - (4 - pi^2/zeta_req)/4 = pi^2/(4*zeta_req)` by pure arithmetic on the prior definition. The residual is identically zero by construction — the assertion cannot fail unless A3 / M3 already failed first. It does not exercise a separate physical claim about `KX/KW`; in particular, the form `1 - x/4` is typed in by hand and never derived from `A_K` (e.g. as `1/A_K(y=0, x)`), so this check provides no independent verification of either the closure ceiling or the stiffness identity.

**Why this matters:**
The check looks like an independent cross-check of "stiffness-ratio form of the floor", but it's algebraically guaranteed by the substitution used to define `x_floor`. A future regression on the closure ceiling (e.g. an off-by-one in `A_K`'s denominator) would not be caught here — only A2/M2 would catch it. The script's print output advertises this as a separate equivalence (`# Equivalent stiffness-ratio form.`), so it is misleading as a verification line.

**Required change:**
Make the check non-tautological in both engines by deriving the LHS from `A_K` rather than from the hand-typed expression `1 - x/4`. Specifically, evaluate `1/A_K(y=0, x)` (which equals `1 - x/4` if and only if the AK definition is consistent with the closure-ceiling assertion), substitute `x -> x_floor`, and check the residual against `pi^2/(4*zeta_req)`. This anchors the check to the AK definition and the x_floor value simultaneously — a regression in either will break it.

SymPy: replace line 56 with
```
expect_zero("KX/KW equivalence", (1 / AK).subs(y, 0).subs(x, x_floor) - pi**2 / (4 * zeta_req))
```

Mathematica: replace line 54 with
```
expectZero["KX/KW equivalence", ((1/aK) /. y -> 0 /. x -> xFloor) - Pi^2/(4 zetaReq)];
```

**Verification:**
After the edit:
- SymPy: line 56 now invokes `AK` (not the hardcoded `1 - x/4`). Output line "KX/KW equivalence = 0" still appears; if `A_K` were changed inconsistently, this line would now flag it where before it would not.
- Mathematica: line 54 references `aK`. Same behavior. The "KX/KW equivalence = 0" / "PASS" lines still appear in the captured output.
- `redteam exec-sympy 055` and `redteam exec-mathematica 055` still exit 0.

## Independent-derivation check (Mathematica)

The Mathematica script is parallel to the SymPy one in structure (same `Omega_exp`, same `A_K`, same limits, same target identities), but the algebra is the natural expression of the underlying closed-form objects in each engine. It is not a line-by-line transliteration: Mathematica hardcodes `xFloor` rather than using `Solve` (SymPy uses `sp.solve`), adds an extra round-trip check (`(zetaMax /. x -> xFloor) - zetaReq`, line 53) that SymPy lacks, and uses `Limit[..., Direction -> "FromAbove"]` for the closure maximum. The two scripts agree on result but differ in path. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both engines report:
- `Omega_exp` = the same rational function (`pi*alpha*(2*alpha*exp(alpha) + pi)/((4*alpha^2 + pi^2)*(exp(alpha) - 1))`).
- `A_K` = `1/(1 - x/4 + x*y^2/pi^2)` (SymPy normalizes to `4*pi^2/(4*x*y^2 + pi^2*(4 - x))`; Mathematica keeps the unnormalized form). Algebraically identical.
- `zeta_twin = 1`, `zeta_max = pi^2/(4 - x)` (SymPy prints `-pi^2/(x - 4)`, same expression).
- `x_floor = 4 - pi^2/zeta_req`.
- All `expect_zero` / `expectZero` residuals are `0`.

No disagreement. The Mathematica `Limit::alimv` warnings on the alpha-limit are about ignoring assumptions involving the limit variable and do not affect the result (the assertions still pass).

## Verdict justification

The substantive limit checks (twin = 1, closure maximum = `pi^2/(4-x)`) hold up under attack: I verified by hand that `Omega_exp -> 1` as `alpha -> 0` (numerator and denominator both `~ pi^2 * alpha`) and `Omega_exp -> pi/2` as `alpha -> oo` (numerator `~ 2*pi*alpha^2*exp(alpha)`, denominator `~ 4*alpha^2*exp(alpha)`), giving `zeta_twin = 1*1 = 1` and `zeta_max = (pi^2/4)*(4/(4-x)) = pi^2/(4-x)`. The `x_floor` closed form follows by direct algebra from `zeta_max = zeta_req`. The only weakness is the trailing "KX/KW equivalence" check in both engines: it merely restates `x_floor`'s definition through arithmetic on `1 - x/4`, so it cannot fail independently. This is a real but cosmetic finding (severity low). Engines agree; outputs are fresh (sympy script mtime 2026-04-01, output 2026-05-11; mathematica script mtime 2026-05-11, output 2026-05-11 — output newer than script in both cases).

## Self-test notes

Checked: (i) variable independence — the proposed fix substitutes `y -> 0` and `x -> x_floor` into `1/AK`; AK depends on both, so neither substitution is trivial. (ii) Parity / symmetry — no integrals introduced. (iii) Trivial-case pre-check — substituting `1/AK = 1 - x/4 + x*y^2/pi^2` at `y = 0` gives `1 - x/4`, then `x -> 4 - pi^2/zeta_req` gives `pi^2/(4*zeta_req)`, so the new residual is zero exactly when `AK` is consistent with the closure-ceiling assertion (no spurious failure). (iv) Path specifications — no new files; edits are to existing `.py` (line 56) and `.wl` (line 54).
