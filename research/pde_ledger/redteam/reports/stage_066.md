---
unit_id: 066
batch: III.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage066_wall_figure_of_merit.md
  paper_appendix: present
---

# Audit unit 066 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_066.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage066_wall_figure_of_merit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 110)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage066_wall_figure_of_merit_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage066_wall_figure_of_merit_mathematica_audit.txt`

## What the paper claims

Stage 066 defines a single dimensionless wall figure of merit and shows that the matched thin-wall success/failure problem collapses to one comparison. From the paper card:

> `\stagefield{Output}`: "The wall control \eqref{eq:app-stage066-Wwall} and exact window \eqref{eq:app-stage066-Wwindow}."

The two boxed equations are:

1. Definition: `W_wall = 4 pi a^2 L^2 J_1 V_0^2 / (T_X ell)`.
2. Exact window: `W_wall <= Pe_req / Delta_inf  =>  fail`, `W_wall >= Pe_req / Delta_0  =>  succeed`; the intermediate band requires the full operator branch.

The source notes additionally state (a) the `kappa`-scaling identity `W_wall = kappa * G_eq^{tw}` (notes §1), (b) monotonicity directions in the six geometric parameters `V0, a, L, ell, T_X, J_1` (notes §3), and (c) a constant-compressibility form `W_H = 4 pi a^2 L^2 I_f V_0^2 / (H_w T_X ell)` with `J_1 -> I_f / H_w` (notes §4). The appendix row at part03 line 110 summarizes: "Single wall control W_wall and exact window."

## What the script claims to verify

The SymPy script defines `W_wall`, `W_fail = Pe_req/Delta_inf`, and `W_suff = Pe_req/Delta_0`, asserts that substituting the Stage-48/065 thresholds `V0_fail^2 = T_X ell Pe_req / (4 pi a^2 L^2 J_1 Delta_inf)` and the corresponding `V0_suff^2` into `W_wall` recovers `W_fail` and `W_suff`, *prints* four partial derivatives (`dW/dV0^2, dW/da, dW/dL, dW/dell`) without sign assertions, and confirms the `W_H` form arises by `J_1 -> I_f/H_w` plus the same threshold inversion. The Mathematica script mirrors all four threshold-inversion checks, additionally **asserts** the `kappa G_eq^{tw}` identity, and **asserts** all six monotonicity signs (`> 0` for V0^2, a, L, J_1; `< 0` for ell, T_X) via `expectTrue`.

## Paper ↔ script cross-check

| Deliverable | SymPy | Mathematica |
|---|---|---|
| `W_wall` definition (paper eq 1) | match (defined; W_H reduction also verified) | match (defined; explicit kappa relation also verified) |
| Exact window `W_wall<=Pe_req/Delta_inf -> fail`, `W_wall>=Pe_req/Delta_0 -> succeed` | match (substitution checks recover thresholds) | match (substitution checks recover thresholds) |
| `kappa`-scaling identity `W_wall = kappa G_eq^{tw}` (notes §1) | missing (no check) | match (`expectZero` at line 48) |
| Monotonicity in 6 directions (notes §3) | partial (derivatives printed only; no sign assertions; `dW/dJ_1` and `dW/dT_X` not even printed) | match (all 6 signs `expectTrue`-asserted) |
| Constant-compressibility form `W_H` (notes §4) | match (substitution chain) | match (substitution chain) |

`paper_alignment: partial` — the notes call out monotonicity as a stage deliverable, but the SymPy engine only prints derivatives without asserting signs and omits two of the six directions entirely. The Mathematica engine satisfies the requirement.

## Assertion inventory

| #  | Script      | Line  | Form                                              | Exercises which paper claim?                            | Anchored to claim?   |
|----|-------------|-------|---------------------------------------------------|---------------------------------------------------------|----------------------|
| A1 | sympy       | 63    | `expect_zero(W_wall(V0_fail^2) - W_fail)`         | Exact window — `V0_fail` ↦ `Pe_req/Delta_inf` map       | yes                  |
| A2 | sympy       | 64    | `expect_zero(W_wall(V0_suff^2) - W_suff)`         | Exact window — `V0_suff` ↦ `Pe_req/Delta_0` map         | yes                  |
| A3 | sympy       | 71-74 | `print(sp.diff(...))` only                        | Monotonicity (notes §3)                                 | no (no assertion)    |
| A4 | sympy       | 79    | `expect_zero(W_wall.subs(J1,I_f/H_w) - W_H)`      | Const-compressibility form definition                   | partial (tautological by construction) |
| A5 | sympy       | 80    | `expect_zero(W_H(V0_fail^2) - W_fail)`            | Exact window in `W_H` form                              | yes                  |
| A6 | sympy       | 81    | `expect_zero(W_H(V0_suff^2) - W_suff)`            | Exact window in `W_H` form                              | yes                  |
| B1 | mathematica | 48    | `expectZero[W_wall - kappa*G_eq^{tw}]`            | `kappa`-scaling (notes §1)                              | yes                  |
| B2 | mathematica | 53    | `expectZero[W_wall(V0_fail^2) - W_fail]`          | Exact window                                            | yes                  |
| B3 | mathematica | 54    | `expectZero[W_wall(V0_suff^2) - W_suff]`          | Exact window                                            | yes                  |
| B4 | mathematica | 73-78 | `expectTrue` on signs of dW/d(V0^2), da, dL, dell, dJ1, dTx | Monotonicity (notes §3, all 6 directions)      | yes                  |
| B5 | mathematica | 85    | `expectZero[W_wall(J1→I_f/H_w) - W_H]`            | Const-compressibility form definition                   | partial (tautological by construction) |
| B6 | mathematica | 86-93 | `expectZero` on `W_H` threshold inversions        | Exact window in `W_H` form                              | yes                  |

Comment on A4/B5: both engines define `W_H` directly as `4 pi a^2 L^2 I_f V_0^2 / (H_w T_X ell)` and then check that substituting `J_1 -> I_f/H_w` in `W_wall` reproduces `W_H`. Algebraically this is `(4 pi a^2 L^2 V_0^2/(T_X ell)) * (I_f/H_w) == W_H`, true by construction. It is mild bookkeeping rather than a substantive identity test. Not raised as a separate finding because both engines do it identically and the same substitution then feeds the non-tautological `W_H` threshold checks (A5/A6, B6).

Comment on A1/A2/B2/B3: `V0_fail^2` is defined inside the script as the Stage-48/065 threshold expression. The check that `W_wall(V0_fail^2) == Pe_req/Delta_inf` therefore verifies that the paper's compressed wall FOM picks up exactly the Stage-48 amplitude threshold — i.e., it is the concrete statement that the multi-parameter wall-amplitude threshold is equivalent to a single-number `W_wall` comparison. This is non-tautological because the constants `4 pi a^2 L^2 J_1 / (T_X ell)` could be mis-set (extra `2`, missing `pi`, swapped numerator/denominator) and the residual would not collapse to zero.

## Findings

### F1 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py:66-74`

**What's wrong:**
Notes §3 enumerates monotonicity as a stage deliverable:

> "So it is strictly monotone in the physically relevant directions: larger wall amplitude V0 increases W_wall, larger throat radius a increases W_wall, larger throat length L increases W_wall, thinner wall width ell increases W_wall, larger support tension T_X suppresses W_wall, and larger weighted wall-profile moment J_1 increases W_wall."

The SymPy script under the `MONOTONICITY` banner only **prints** four partial derivatives (`dW/dV0^2, dW/da, dW/dL, dW/dell`) without any assertion on their signs, and entirely omits `dW/dJ_1` and `dW/dT_X`:

```python
print("dW/d(V0^2) =", sp.diff(W_Vp, Vp))
print("dW/da =", sp.diff(W_wall, a))
print("dW/dL =", sp.diff(W_wall, L))
print("dW/dell =", sp.diff(W_wall, ell))
```

A sign flip or zero derivative would not be caught — the script would still print "STAGE 49 AUDIT PASSED" and exit 0. The Mathematica engine asserts all six signs via `expectTrue` (lines 73-78 of the `.wl`), so this is an asymmetric coverage gap where one engine certifies the claim and the other does not.

**Why this matters:**
The paper status for Stage 066 is `\StatusExactClosure{}`, meaning both engines should independently certify the closure. As written, the SymPy engine does not certify the monotonicity portion of the closure, and is also incomplete (missing `J_1`, `T_X`). If a future edit introduced a sign error (e.g., dropping a sign on `T_X` in the denominator), Mathematica would fail but SymPy would silently pass.

**Required change:**
After line 74 of the SymPy script, add explicit positivity / negativity assertions matching the Mathematica engine. Use the existing `positive=True` declarations so SymPy can resolve the inequalities. Specifically, augment the `MONOTONICITY` block as follows (replace lines 66-74; lines 66-69 of the section header and `Vp` substitution stay):

```python
banner("MONOTONICITY")

# W_wall is manifestly monotone in V0^2, a^2, L^2, J1 and inverse ell, T_X.
Vp = sp.symbols("Vp", positive=True, real=True)
W_Vp = sp.simplify(W_wall.subs(V0**2, Vp))

dW_dV0sq = sp.simplify(sp.diff(W_Vp, Vp))
dW_da    = sp.simplify(sp.diff(W_wall, a))
dW_dL    = sp.simplify(sp.diff(W_wall, L))
dW_dell  = sp.simplify(sp.diff(W_wall, ell))
dW_dJ1   = sp.simplify(sp.diff(W_wall, J1))
dW_dTX   = sp.simplify(sp.diff(W_wall, TX))

print("dW/d(V0^2) =", dW_dV0sq)
print("dW/da =", dW_da)
print("dW/dL =", dW_dL)
print("dW/dell =", dW_dell)
print("dW/dJ1 =", dW_dJ1)
print("dW/dT_X =", dW_dTX)

assert sp.simplify(dW_dV0sq > 0) is sp.true, "dW/d(V0^2) should be positive"
assert sp.simplify(dW_da    > 0) is sp.true, "dW/da should be positive"
assert sp.simplify(dW_dL    > 0) is sp.true, "dW/dL should be positive"
assert sp.simplify(dW_dell  < 0) is sp.true, "dW/dell should be negative"
assert sp.simplify(dW_dJ1   > 0) is sp.true, "dW/dJ1 should be positive"
assert sp.simplify(dW_dTX   < 0) is sp.true, "dW/dT_X should be negative"
```

Because all symbols `V0, ell, a, L, J1, T_X, Vp` are declared with `positive=True`, the simplified derivatives are manifestly signed monomials (e.g. `4*pi*J1*L**2*a**2/(T_X*ell)` for `dW/dVp`, and `-4*pi*J1*L**2*V0**2*a**2/(T_X**2*ell)` for `dW/dT_X`). SymPy `simplify(expr > 0)` returns `sp.true` for such manifestly-positive rationals over positive reals; if any sign actually fails it will return `sp.false` or remain a Relational, both of which trip the `assert ... is sp.true` guard.

**Verification:**
Re-run `redteam exec-sympy 066`. New lines `dW/dJ1 = ...` and `dW/dT_X = ...` must appear in the output, and the script must still exit 0. No existing pre-line-74 output line should change. The new assertions must hold; if any fail, the monotonicity claim in notes §3 is wrong and the directive must be revisited.

## Independent-derivation check (Mathematica)

Both scripts perform the same algebraic substitutions for the threshold inversions (define W_wall; substitute V0_fail^2 / V0_suff^2; reduce J_1 -> I_f/H_w; substitute again into W_H). The Mathematica engine **adds** an independent check tying W_wall to `kappa G_eq^{tw}` (line 48), using `K_X = T_X kappa / L^2` as the `kappa` definition from Stage 44 — a route the SymPy engine does not exercise. The Mathematica engine also derives sign decisions on six monotonicity directions via `FullSimplify[..., Assumptions]`, again a route not exercised in SymPy. Symbol naming differs (`v0/ell/a/len/j1/tx` vs. `V0/ell/a/L/J1/T_X`), and the assumption blocks are written in different idioms. The pattern of intermediate steps for the threshold inversions is similar (which is expected for trivial algebra), but the kappa-bridge and the monotonicity assertion block are genuinely independent in Mathematica. No `mathematica_transliteration` finding.

## Engine cross-check

The two engines agree on every overlapping check:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `W_wall` | `4*pi*J1*L**2*V0**2*a**2/(T_X*ell)` | `(4*a^2*j1*len^2*Pi*v0^2)/(ell*tx)` |
| `W_fail` | `Pe_req/Delta_inf` | `peReq/deltaInf` |
| `W_suff` | `Pe_req/Delta_0` | `peReq/delta0` |
| `W_wall(V0_fail^2) - W_fail` | `0` | `0` |
| `W_wall(V0_suff^2) - W_suff` | `0` | `0` |
| `W_H` | `4*pi*I_f*L**2*V0**2*a**2/(H_w*T_X*ell)` | `(4*a^2*ifMom*len^2*Pi*v0^2)/(ell*hw*tx)` |
| `J_1 = I_f/H_w` reduction | `0` | `0` |
| `W_H(V0_fail^2) - W_fail`, `W_H(V0_suff^2) - W_suff` | `0`, `0` | `0`, `0` |
| `dW/dV0^2, dW/da, dW/dL, dW/dell` (numerator/denominator structure) | matches | matches |

Mathematica covers additionally `W_wall - kappa G_eq^{tw} = 0` (`PASS`) and the six monotonicity signs (`PASS`). No disagreement; the asymmetry is a coverage gap in SymPy (see F1), not a numerical inconsistency.

## Verdict justification

The paper claims (W_wall definition and the exact two-sided wall-FOM window) are both verified in both engines, non-tautologically, via the substitution chain `V0_fail^2 -> W_fail` and `V0_suff^2 -> W_suff`. The constant-compressibility form is also covered. The notes-side deliverables (kappa-scaling and 6-direction monotonicity) are covered only by Mathematica. The SymPy engine prints derivatives without asserting their signs and omits two of the six directions entirely; this is a real but low-severity coverage gap, hence verdict `findings` with `paper_alignment: partial`. No engine disagreement, no tautological core claim (the W_H reduction is tautological by construction but is a bookkeeping step into the non-tautological W_H threshold checks), no hardcoded numeric constant, no missing script, outputs fresh (sympy output mtime > sympy script mtime; mathematica output mtime > mathematica script mtime, per `stat` check). Attacks tried that failed: (1) probe the V0_fail/V0_suff threshold inversion for hidden prefactor errors — robust; (2) attempt to find tautological structure in the central threshold-inversion chain — the substitution check is non-trivial; (3) probe for symbol-domain mistakes (positive-real assumptions on lengths/areas/times/dimensionless thresholds) — all consistent with physical setup; (4) probe whether the kappa-bridge identity in Mathematica hides a transliteration of SymPy algebra — it uses a route SymPy does not exercise. The one attack that succeeded: monotonicity coverage check on the SymPy side (F1).

## Self-test notes

(1) Variable independence: `W_wall = 4 pi a^2 L^2 J_1 V0^2 / (T_X ell)` depends on all six of {V0, a, L, ell, J_1, T_X}, so each `sp.diff(W_wall, var)` is generically nonzero; with `positive=True` declarations the derivatives are manifestly-signed rational monomials. (2) Sign check: under the existing `positive=True` declarations, `simplify(dW/dT_X < 0)` resolves to `sp.true` because `T_X` appears only as `1/T_X` with positive numerator, and similarly for the others; `assert ... is sp.true` is the correct comparison form for SymPy `Relational` evaluation. (3) Paper round-trip: the proposed F1 fix adds only sign assertions on derivatives of the same `W_wall` the paper defines, with no new constants and no new dependencies — no new `paper_misalignment` is introduced; the fix brings SymPy into line with Mathematica's already-asserted signs and with notes §3.
