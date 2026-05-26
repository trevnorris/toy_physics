---
unit_id: 054
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00-06:00
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage054_robin_softening_support_lane.md
  paper_appendix: present
---

# Audit unit 054 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_054.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage054_robin_softening_support_lane.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 86)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage054_robin_softening_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage054_robin_softening_mathematica_audit.txt`

## What the paper claims

From `stage_054.tex` `\stagefield{Output}`: "Robin root \eqref{eq:app-stage054-robin-root}, softening factor \eqref{eq:app-stage054-AK}, and ceiling \eqref{eq:app-stage054-AK-window}." The three deliverables are:

1. **Robin root**: `y tan y = eta` with `0 < y < pi/2` (eq. app-stage054-robin-root).
2. **Softening factor**: `A_K(eta) = K_W^eff / K_(phi,0)^eff = 1 / (1 - x/4 + x y^2/pi^2)`, with `x = pi^2 T_X/(L^2 K_W^eff)` and `0 < x < 4` (eq. app-stage054-AK).
3. **Endpoint window (ceiling)**: `1 <= A_K <= 4/(4-x)` (eq. app-stage054-AK-window).

A further consequence (`zeta_req <= 4/(4-x)` for pure-softening rescue) is stated as a remark below the boxed equations. The notes give the boundary-value setup (Robin at `s=0`, Neumann at `s=L`, `psi'' + k^2 psi = 0`) and explicitly state the supporting monotonicity: "the map `y -> eta = y tan y` is strictly increasing on (0,pi/2), and `A_K` is strictly decreasing in y," which is what gives the window its closure structure.

## What the script claims to verify

Both scripts purport to verify, from the BVP `psi'' + k^2 psi = 0` with `psi'(L)=0`, `psi'(0)=h psi(0)`:
(a) the Robin characteristic equation `k tan(kL) = h` (and its dimensionless form `y tan y = eta`);
(b) algebraic equivalence of the explicit ratio `K_W^eff/K_(phi,0)^eff` with the `x`-form `1/(1 - x/4 + xy^2/pi^2)`;
(c) the endpoint values `A_K(y=pi/2)=1` and `A_K(y->0+)=4/(4-x)`;
(d) the saturation floor `x_floor = 4 - 4/zeta_req` obtained by inverting `A_K,max = zeta_req`.

The Mathematica script now derives `bExpr` via `Solve[...]` (line 40) and `xFloor` via `Solve[...]` (line 84), so the previous v1 `hardcoded_result` findings on those two lines no longer apply.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) Robin root `y tan y = eta` | SymPy A1, A2 (lines 39, 42); MMa B1, B2 (lines 46, 47-50) — derive `k tan(kL) - h = 0` from BVP and rescale to `y tan y = eta`. | match |
| (2) Softening factor `A_K = 1/(1 - x/4 + xy^2/pi^2)` | SymPy A3, A4 (lines 59, 62); MMa B3, B4 (lines 67, 77) — verify `K_W` identity and `A_K` x-form equality. | match |
| (3) Window `1 <= A_K <= 4/(4-x)` | SymPy A5, A6 (lines 69, 70); MMa B5, B6 (lines 78, 79) — verify endpoint values only. | partial (see F1) |
| Remark: rescue criterion `zeta_req <= 4/(4-x)` (and `x_floor = 4 - 4/zeta_req`) | SymPy A7 (line 85); MMa B7, B8 (lines 89, 90). | match |

The script also asserts `A_K,max(x_floor) - zeta_req == 0` (MMa B8) which is a downstream consistency check; not in the paper card explicitly but follows trivially from B7.

`paper_alignment: aligned` — every paper-side deliverable is mapped to at least one script-side check; one row is "partial" but the gap is a verification-strength gap, not a paper/script disagreement on what is claimed.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `expect_zero("Robin equation -> k tan(kL) - h", char_eq/A + h - k*tan(kL))` | claim 1 | yes (B derived via `sp.solve`) |
| A2 | sympy | 42 | `expect_zero("dimensionless form", ...)` | claim 1 | yes |
| A3 | sympy | 59 | `expect_zero("K_W identity", KW - (KX_from_x + ...))` | claim 2 (setup) | yes |
| A4 | sympy | 62 | `expect_zero("A_K x-form", AK_x - 1/(1 - x/4 + x y^2/pi^2))` | claim 2 | yes |
| A5 | sympy | 69 | `expect_zero("DN limit", AK_DN - 1)` | claim 3 (endpoint) | partial (endpoint only) |
| A6 | sympy | 70 | `expect_zero("soft-mouth limit", AK_soft - 4/(4-x))` | claim 3 (endpoint) | partial (endpoint only) |
| A7 | sympy | 85 | `expect_zero("x floor = 4 - 4/zeta_req", x_floor - (4 - 4/zeta_req))` | rescue remark | yes (x_floor from `sp.solve`) |
| B1 | mma | 46 | `expectZero["Robin equation -> k tan(kL) - h", charEq/a + h - k Tan[k ell]]` | claim 1 | yes (bExpr from `Solve`, line 40) |
| B2 | mma | 47-50 | `expectZero["dimensionless form", ...]` | claim 1 | yes |
| B3 | mma | 67 | `expectZero["K_W identity", kWBar - (kXFromX + ...)]` | claim 2 (setup) | yes |
| B4 | mma | 77 | `expectZero["A_K x-form", aKX - aKSym]` | claim 2 | yes |
| B5 | mma | 78 | `expectZero["DN limit", aKDN - 1]` | claim 3 (endpoint) | partial (endpoint only) |
| B6 | mma | 79 | `expectZero["soft-mouth limit", aKSoft - 4/(4 - x)]` | claim 3 (endpoint) | partial (endpoint only) |
| B7 | mma | 89 | `expectZero["x floor = 4 - 4/zeta_req", xFloor - (4 - 4/zetaReq)]` | rescue remark | yes (xFloor from `Solve`, line 84) |
| B8 | mma | 90 | `expectZero["A_K,max(x_floor) - zeta_req", (aKMax /. x -> xFloor) - zetaReq]` | rescue remark consistency | partial (depends on B7) |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py:64-70`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl:70-79`

**What's wrong:**
The paper claims the **endpoint window** `1 <= A_K <= 4/(4-x)` (eq. `app-stage054-AK-window`). The notes (section 3) further state that this window is established by monotonicity:

> "Because the map `y -> eta = y tan y` is strictly increasing on (0,pi/2), and `A_K` is strictly decreasing in y, the softening factor is strictly decreasing in eta."

In both scripts the window is exercised only at the two endpoints:
- SymPy lines 65–70: `AK_DN = AK_sym.subs(y, pi/2)` and `AK_soft = sp.limit(AK_sym, y, 0, dir="+")`, then `expect_zero("DN limit", AK_DN - 1)` and `expect_zero("soft-mouth limit", AK_soft - 4/(4-x))`.
- Mathematica lines 71–79: identical structure, with `aKDN = aKSym /. y -> Pi/2` and `aKSoft = Limit[aKSym, y -> 0, ...]`.

Computing two endpoint values does not establish a window — two numbers do not bracket a continuous function unless monotonicity (or some other order-preserving property) is verified between them. The paper's boxed equation `1 <= A_K <= 4/(4-x)` is therefore only partially verified by the assertions.

A non-tautological monotonicity check is direct: `A_K(y) = 1/(1 - x/4 + x y^2/pi^2)` is decreasing in y on `(0, pi/2)` iff its derivative w.r.t. y is non-positive there. With `D(y) = -d/dy (1 - x/4 + x y^2/pi^2) = -2 x y/pi^2`, the sign on `0 < y < pi/2, 0 < x < 4` is negative. The script could verify `sign(d A_K / d y) == -sign(2 x y / pi^2 * A_K^2)` (i.e., the derivative residual `dAK/dy + 2 x y A_K^2/pi^2 == 0`), which is an exact algebraic identity, then assert positivity of the prefactor `2 x y / pi^2` on the assumed domain.

**Why this matters:**
If the closed form for `A_K(y)` had a hidden non-monotonic regime (e.g., a sign flip from a misderived denominator term), the endpoint values alone could still match `1` and `4/(4-x)` while the interior fails the window. The endpoint-only check is necessary but not sufficient for the paper's window claim. The script's own docstring/printout also asserts the window verbatim ("`A_K = 1 / [1 - x/4 + x y^2/pi^2]. It ranges from 1 up to 4/(4-x)`," final-ledger print) so the assertion strength does not match the script's own stated scope, independent of the paper.

**Required change:**
In both scripts, add a monotonicity check on `A_K(y)` as a function of y on `(0, pi/2)` with `0 < x < 4`. The simplest non-tautological form:

- Compute the symbolic derivative `dAK_dy = sp.diff(AK_sym, y)` (SymPy) / `D[aKSym, y]` (Mathematica).
- Assert the residual `expect_zero("dAK/dy closed form", dAK_dy - (-2*x*y / (pi**2 * (1 - x/4 + x*y**2/pi**2)**2)))`.
- Print a comment noting that on `0 < y < pi/2`, `0 < x < 4` the prefactor `2 x y / pi^2 > 0`, hence `dAK/dy < 0`, hence the endpoint values bracket the function.

For the Mathematica script, mirror with `D[aKSym, y]` and a `FullSimplify` reduction.

**Verification:**
After the fix, the SymPy transcript prints a new line of the form `dAK/dy closed form = 0` and a corresponding PASS in the Mathematica transcript. The new assertion fails if and only if the closed form for `A_K(y)` is mis-stated; together with the existing endpoint checks, this exercises the full window claim.

## Independent-derivation check (Mathematica)

The `.wl` follows the same variable choreography as the `.py`: `psi = A cos(ks) + B sin(ks)`; solve for B from Neumann at L; compute `charEq` at `s = 0`; then `kPhi`, `kWEff`, `aK`; then `kXFromX`, `tXFromX`, `aKX`, `aKSym`, `aKDN`, `aKSoft`; then `ineqRhs`, `yReqSq`, `aKMax`, `xFloor`. The two engines both invoke `Solve` (sympy `sp.solve` / mathematica `Solve`) at lines 33/40 (`B`) and 83/84 (`xFloor`).

The structural parallelism is high — the .wl can be read as a renamed re-implementation of the .py — but the algorithmic steps (BVP ansatz → Neumann elimination → characteristic equation → effective stiffness → x-substitution → endpoint limits → saturation solve) are the standard derivation for this BVP, and each `Solve` call is independently exercised in each engine. I do not file `mathematica_transliteration` here: the parallelism is at the level of "both engines run the textbook derivation," not "the .wl echoes the .py's algebraic shortcuts." Each `Solve` and each `FullSimplify` operates independently on its engine's symbolic store; if the SymPy version had a hidden algebraic shortcut, the Mathematica version would still recompute from the ansatz.

Caveat: a more strongly independent re-derivation would use `psi = C cos(k(L-s))` (auto-satisfying Neumann at L) and impose Robin at `s=0`, avoiding the two-coefficient ansatz entirely. That stronger form is not present in either engine. Within the rubric's tolerance for "both engines derive from physical premises," I judge the current pair to clear the bar.

## Engine cross-check

Both engines pass all their assertions. Final closed forms (from saved outputs):

- SymPy: `A_K in x,y form = 4*pi**2/(4*x*y**2 - pi**2*(x - 4))`, `A_K,max = -4/(x - 4)`, `x floor at saturation = 4 - 4/zeta_req`.
- Mathematica: `A_K in x,y form = kWBar/(kWBar + kWBar*x*(-1/4 + y^2/Pi^2))`, `A_K,max = -4/(-4 + x)`, `x floor at saturation = 4 - 4/zetaReq`.

These are algebraically equal up to engine-specific normalization: `4 pi^2 / (4 x y^2 - pi^2 (x-4))` = `1/(1 - x/4 + x y^2/pi^2)` (multiply numerator and denominator by `4/pi^2`), and the Mathematica form factors `kWBar` out trivially. `-4/(x-4) = 4/(4-x)`. No `engine_disagreement` finding.

Output freshness: SymPy script mtime Apr 1 12:39, output mtime May 22 17:38 — output newer than script (fresh). Mathematica script mtime May 22 17:38, output mtime May 22 17:38 — output ≥ script mtime (fresh). No `stale_output` finding.

## Verdict justification

Paper-side and script-side both claim the Robin root, softening factor, and softening window. The Robin root and softening factor are exercised with substantive `Solve`-derived checks in both engines, and the engines agree. The previous v1 hardcoded-result findings have been addressed: the Mathematica script now derives `bExpr` and `xFloor` via `Solve` calls. One gap remains: the **window claim** (`1 <= A_K <= 4/(4-x)`) is verified only at its two endpoints, with no monotonicity check connecting them. The paper claim is fully aligned with the script's stated scope; the assertion strength is the issue, not the target. Verdict: `findings`, one `insufficient_verification` finding. No stop-cold.

Attacks tried that failed:
- Banner/comment text says "STAGE 37" / "STAGE 037" while the filename and the paper card are stage 054 — this is a renumbering artifact in print/banner statements only; it does not enter any assertion residual and does not change which math is being verified. Not flagged.
- The `0 < y < pi/2` paper constraint vs script's looser `y > 0` assumption: the `A_K x-form` identity is algebraic in y and holds on the full positive line, so the looser assumption is harmless for the assertions present. Not flagged.
- The Mathematica `Limit::alimv` warning at the `y -> 0+` limit: examined; the warning is informational and the simplified result `-4/(-4 + x)` is the correct soft-mouth limit. Not a bug.
- The Robin sign convention `psi'(0) = h psi(0)` is reproduced consistently in both scripts via `(D[psi,s] /. s -> 0) - h*(psi /. s -> 0)`, matching the notes. Sign correct.

## Self-test notes

- **Variable independence (proposed F1 fix)**: `AK_sym` (SymPy) and `aKSym` (Mathematica) both depend on `y` explicitly, so `sp.diff(AK_sym, y)` / `D[aKSym, y]` will not be identically zero. Confirmed.
- **Trivial-case pre-check (proposed F1 fix)**: substitute `x = 2, y = pi/4`. Then `1 - x/4 + x y^2/pi^2 = 1 - 1/2 + 2*(pi^2/16)/pi^2 = 1/2 + 1/8 = 5/8`, so `A_K = 8/5`. `dA_K/dy = -2*x*y/(pi^2) * A_K^2 = -2*2*(pi/4)/(pi^2) * (8/5)^2 = -(1/pi) * 64/25 = -64/(25 pi)`, negative — confirms the prefactor sign claim. Residual `sp.diff(AK_sym, y) - (-2*x*y / (pi**2 * (1 - x/4 + x*y**2/pi**2)**2))` reduces to 0 by direct computation (both expressions are the same closed form for the derivative of `A_K`).
- **Paper round-trip (proposed F1 fix)**: the new monotonicity check adds an algebraic identity for `dA_K/dy` directly from the paper's boxed `A_K` form; it does not introduce any new constant or convention that would create a fresh `paper_misalignment`.
- **Path specifications**: F1 modifies existing `.py` and `.wl` files at `scripts/` and `mathematica/`; no new file creation.
