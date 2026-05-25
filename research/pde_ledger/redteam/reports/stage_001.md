---
unit_id: 001
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md
  paper_appendix: present
---

# Audit unit 001 red-team report (v2 — paper-grounded)

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_001.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage001_geometry_lift_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage001_geometry_lift_mathematica_audit.txt`

Script mtimes: sympy 2026-04-21 17:04; mathematica 2026-05-21 11:42. Output mtimes: sympy 2026-05-22 12:38; mathematica 2026-05-22 12:40. Outputs are newer than scripts — no `stale_output` finding.

## What the paper claims

Stage 001 lifts the old finite-dimensional `(a(t), L(t))` geometry to a distributed throat shape field `R(Omega, w, t)` via the level set `Sigma = r - R`. The `\stagefield{Output}` reads verbatim: *"The stage outputs the shape-field closure (eq:app-stage001-sigma), the confinement variation (eq:app-stage001-delta-v), the harmonic split (eq:app-stage001-harmonic-expansion), the modal wall PDE (eq:app-stage001-modal-wall-pde), and the response-operator target (eq:app-stage001-response-operator)."* The four explicit `\stagefield{Checks}` are: (i) `q_{00} = 2 sqrt(pi) delta a`; (ii) `delta V_conf = -(V'_wall(Sigma_0/ell_c)/ell_c) eta` from `delta Sigma = -eta`; (iii) no `l=0 ↔ l=2` mixing from spherical orthogonality; (iv) the ontology channels `A_w, J^w, F_{mu w}` are kept alive. The boxed modal wall PDE (eq:app-stage001-modal-wall-pde) is

```
mu_eta partial_t^2 q_{lm} - partial_w(T_w partial_w q_{lm}) + [K_eta + l(l+1) T_Omega] q_{lm}
    = S_{lm}^{(psi,A)} + f_{lm}^{ext}                                                 (positive RHS)
```

The boxed localized Maxwell equation (eq:app-stage001-linear-maxwell) is

```
partial_M(Z(w) delta F^{MN}) + (1/xi) partial^N(partial . delta A) = mu_0 delta J^N    (+1/xi gauge term)
```

The notes restate the wall PDE with the same `+S + f_ext` RHS sign and describe the same Lagrangian-based derivation route.

## What the script claims to verify

The SymPy script verifies (I) real-harmonic bookkeeping (norms, zero averages of the P2 set, spherical-Laplacian eigenvalues, `q_00 = 2 sqrt(pi) delta a`); (II) the linearized confinement variation `d/d_eps Vwall((Sigma_0 - eps eta)/ell_c)|_{eps=0} = -eta V'_wall(Sigma_0/ell_c)/ell_c`; (III) the modal wall Euler-Lagrange equation in densitized and weighted forms plus `K_l` specialization at `l=0,2` and a sourced variant; and (IV) the linearized localized Maxwell equation derived from a chosen Lagrangian. The Mathematica script mirrors all of these and adds one independent cross-check (Section I.3b) of the spherical Laplacian eigenvalues via `SphericalHarmonicY`. The summary text in both scripts asserts the verified equations are the paper's modal wall PDE and the paper's localized Maxwell equation.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Shape-field closure `Sigma = r - R` (eq sigma) | Implicit in Section II's `expr = Vwall((Sigma_0 - eps eta)/ell_c)`; uses `delta Sigma = -eta` correctly | match |
| Confinement variation `delta V_conf = -(V'_wall/ell_c) eta` (eq delta-v) | Section II `expect_zero("linearized confinement variation", first_var - target)` | match |
| Harmonic split with `Y_00`, real `P_2` basis (eq harmonic-expansion) | Section I.1 norms, zero averages, orthogonality (`cross(Y00, Y20)`); zero averages of `Y_{2m}` imply orthogonality to constant `Y00` | match |
| `q_00 = 2 sqrt(pi) delta a` (eq q00-da) | Section I.2 `expect_zero("mouth average - q00/(2 sqrt(pi))")` | match |
| Modal wall PDE LHS operator `mu q_tt - d_w(T_w q_w) + (K_eta + l(l+1) T_Omega) q` (eq modal-wall-pde) | Sections III.1, III.2, III.3 | match |
| Modal wall PDE *RHS sign* `= +S_{lm}^{(psi,A)} + f_{lm}^{ext}` | Section III.4 `ldens_forced = ldens - q * source_total` ⇒ verified equation is `= -source_total` | **mismatch** (F1) |
| Localized Maxwell LHS-operator structure (eq linear-maxwell) with `+(1/xi) partial^N(partial . delta A)` | Section IV: assertion verifies `partial_w(Z F_wx) - partial_x(divA)/xi - mu0 J_x = 0` | **mismatch** under mostly-plus signature (F2, may be metric-signature artifact) |
| Response-operator target `j_A = sum_B Z_{eff,AB} u_B` (eq response-operator) | Definitional; no scriptable identity; Section IV retention of `Aw, Jw` is consistent with the port basis | match (definitional, not testable) |
| Ontology: `A_w, J^w, F_{mu w}` retained (Check iv) | Section IV explicitly carries `Aw`, `Jw`, `Fwx = ∂_w A_x − ∂_x A_w` | match |

`paper_alignment` set to **partial**: bookkeeping checks (I.1, I.2, I.3, II, III.1-3) faithfully exercise the paper's Output and Checks; the source-coupling (III.4) and gauge-term (IV) signs disagree with the paper's boxed equations.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 89-93 | `expect_zero(integral(Y_{2m} dOmega))` x5 | harmonic split + Check (iii) | yes |
| A2 | sympy | 96-98 | `expect_zero(norm(Y00) - 1)` | harmonic split normalization | yes |
| A3 | sympy | 100-102 | `expect_zero(norm(Y20) - 1)` | harmonic split normalization | yes |
| A4 | sympy | 104-106 | `expect_zero(cross(Y00, Y20))` | Check (iii) no l=0 ↔ l=2 mixing | yes |
| A5 | sympy | 114 | `expect_zero(mouth_avg - q00/(2 sqrt(pi)))` | Check (i) `q_00 = 2 sqrt(pi) delta a` | yes |
| A6 | sympy | 125 | `expect_zero(-Delta_S2 Y00)` | spherical-Laplacian identity needed for modal-wall PDE | yes |
| A7 | sympy | 127 | `expect_zero(-Delta_S2 Y_{2m} - 6 Y_{2m})` x5 | same | yes |
| A8 | sympy | 141 | `expect_zero(linearized confinement variation)` | Check (ii) + boxed eq delta-v | yes |
| A9 | sympy | 167 | `expect_zero(EL[ldens] - target_dens)` | boxed modal-wall-pde LHS operator | yes |
| A10 | sympy | 177 | `expect_zero(EL[g·ldens] - target_weighted)` | weighted-vs-densitized convention (Verification line of card) | yes |
| A11 | sympy | 181-182 | `expect_zero(K_l|ell=0 - K_eta)`, `expect_zero(K_l|ell=2 - (K_eta + 6 T_Omega))` | `l(l+1)` specialization in modal-wall-pde | yes |
| A12 | sympy | 192 | `expect_zero(EL[ldens − q·source] − (target_dens − source))` | boxed modal-wall-pde *RHS sign* | **no — opposite sign from paper** |
| A13 | sympy | 228 | `expect_zero(EL_Ax − target_Ax)` | boxed linear-maxwell + ontology (A_x) | partial — operator yes, gauge sign maybe-no |
| A14 | sympy | 229 | `expect_zero(EL_Aw − target_Aw)` | boxed linear-maxwell + ontology (A_w) | partial — operator yes, gauge sign maybe-no |
| M1-M14 | mathematica | various | mirrors A1-A14 | same paper claims | same as A1-A14 |
| M15 | mathematica | 138-139 | `expectZero("SphericalHarmonicY[l,0] lap eigenvalue")` (independent) | independent cross of A6/A7 | yes |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** target_mismatch (sign convention on source coupling)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:188-192`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:184-192`
- paper side: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_001.tex:158-162` (boxed eq:app-stage001-modal-wall-pde)
- notes side: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md:358-361`

**What's wrong:**

Paper boxed eq:app-stage001-modal-wall-pde (`stage_001.tex:156-163`):

```
mu_eta partial_t^2 q_{lm} - partial_w(T_w partial_w q_{lm}) + [K_eta + l(l+1) T_Omega] q_{lm}
    = S_{lm}^{(psi,A)} + f_{lm}^{ext}
```

Notes (lines 358-361) match: `... = S_{lm}^{(psi,A)} + f_{lm}^{ext}`. Both state the source enters the RHS with a **positive** sign.

SymPy script `sympy_audit.py:188-192`:

```python
source_total = S_lm(t, w) + f_ext(t, w)
ldens_forced = ldens - q(t, w) * source_total
el_forced = euler_equations(ldens_forced, q(t, w), [t, w])[0]
target_forced = target_dens - source_total
expect_zero("sourced densitized Euler-Lagrange equation", el_forced.lhs - target_forced)
```

SymPy's `euler_equations` returns `Eq(∂L/∂q − Σ_i d/dx_i ∂L/∂(∂_{x_i} q), 0)`. With `ldens_forced = ldens − q · source_total`, `∂L/∂q` gains an extra `−source_total`, so `el_forced.lhs = target_dens − source_total = 0`, i.e.,

```
mu_eta partial_t^2 q − partial_w(T_w partial_w q) + K_l q = −(S_lm + f_ext).
```

That is the **negative** of the paper's boxed RHS. The Mathematica script (`mathematica_audit.wl:187-192`) does the identical thing: `ldensForced = ldens − qField*sourceTotal` and asserts the same sign-flipped form. The summary line in both scripts ("the sourced modal-wall RHS forcing S_lm^(psi,A) + f_lm^ext") explicitly claims to verify the paper's positive-RHS form, but the assertion exercises the opposite sign.

**Why this matters:**

The script silently certifies a convention opposite to the one the paper publishes. Downstream stages quoting eq:app-stage001-modal-wall-pde will get the wrong sign on the source if they trust the script's sign convention rather than the paper's. Stage 003 in particular uses this modal-wall PDE form as the wall-BdG coupling and depends on the correct sign for the back-reaction direction. For a checkpoint stage, this is a real defect.

**Required change:**

Script-side, single-sign flip: the Lagrangian coupling term should be `+q · source` (so that `∂L/∂q` produces `+source`, leaving `+source` on the equation RHS), not `−q · source`.

- SymPy (`sympy_audit.py:189`): change `ldens_forced = ldens - q(t, w) * source_total` → `ldens_forced = ldens + q(t, w) * source_total`. Update line 191 `target_forced = target_dens - source_total` → `target_forced = target_dens + source_total`.
- Mathematica (`mathematica_audit.wl:188`): change `ldensForced = ldens - qField*sourceTotal` → `ldensForced = ldens + qField*sourceTotal`. Update line 191 `targetForced = targetDens - sourceTotal` → `targetForced = targetDens + sourceTotal`.

After this fix, the equation verified will be `mu q_tt − ∂_w(T_w q_w) + K_l q = +(S_lm + f_ext)`, matching the paper.

**Verification:**

Re-run the script; output line "sourced densitized Euler-Lagrange equation = 0" should still print `0` (the identity remains true under consistently-flipped signs on both sides). The substantive change is in the Lagrangian → equation correspondence the assertion exercises. Manual inspection: the new `ldens_forced` line should show `+ q(t,w) * source_total` and the new `target_forced` line should show `+ source_total`.

### F2 — paper_misalignment

**Severity:** low (likely metric-signature convention question, not a substantive sign error in the math)
**Subtype:** target_mismatch (gauge-term sign in linearized Maxwell — pending user resolution on metric signature)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_001.tex:192-198` (boxed eq:app-stage001-linear-maxwell)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:209-229`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:194-213`

**What's wrong:**

Paper boxed eq:app-stage001-linear-maxwell reads
```
partial_M(Z(w) delta F^{MN}) + (1/xi) partial^N(partial . delta A) = mu_0 delta J^N,
```
i.e., the gauge-fix term has a `+(1/xi) partial^N` sign. The raised-index `partial^N` makes the sign of the term as written in `partial_N` coordinates dependent on metric signature.

SymPy script defines (lines 209-213):
```python
lmax = (
    sp.Rational(1, 2) * Zloc(w) * Fwx**2
    - sp.Rational(1, 2) * divA**2 / gauge_xi
    + mu0 * (Jx(x, w) * Ax(x, w) + Jw(x, w) * Aw(x, w))
)
```
and asserts (line 228) `el_Ax - target_Ax = 0` with
```python
target_Ax = sp.diff(Zloc(w) * Fwx, w) - sp.diff(divA, x) / gauge_xi - mu0 * Jx(x, w)
```
Setting this to zero gives `partial_w(Z F_wx) - partial_x(divA)/xi = mu_0 J_x`.

Reading the paper's boxed equation under a mostly-plus metric (`partial^x = partial_x` for spatial indices): `partial_w(Z F^{wx}) + partial_x(divA)/xi = mu_0 J_x`. **The gauge-fix term has opposite sign in script vs paper-mostly-plus.**

If the parent theory uses a `(+,−,−,−,−)` mostly-minus signature, then `partial^x = −partial_x` and the paper's `+(1/xi) partial^N(...)` rewrites as `−(1/xi) partial_x(...)` for spatial N, which would match the script. The paper card does not state the metric signature, and the audit prompt restricts the auditor's reading to this unit's files.

**Why this matters:**

If the signature is mostly-plus, the script verifies a Maxwell equation with a sign-flipped gauge term — internally consistent against its own Lagrangian but not the paper's published form. If mostly-minus, no defect. Either way, the metric convention is load-bearing for sign comparisons and should be made explicit either in the paper card (one-line statement of signature) or in the script docstring.

**Required change:**

Routes to user — see directive `## Resolve before fix_loop` block. The auditor cannot determine the metric signature from the audit reading list. Codex must not silently flip the sign in either direction.

**Verification:**

After the user picks a direction: either (a) update the script's `lmax` gauge term to `+divA^2/(2 xi)` (mostly-plus path; flips EL gauge-term sign to match) — `sympy_audit.py:211` swap `-` to `+`, and `mathematica_audit.wl:205` swap `-` to `+`; and update both `target_Ax`/`target_Aw` and `target_Ax`/`target_Aw` Mathematica accordingly; or (b) leave scripts alone and add a one-line metric-signature note to the paper card or script docstring.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally a transliteration of the SymPy script: same harmonic basis hardcoded (Section I), same chain-rule call (Section II), same EulerEquations / VariationalD wrapping with the same `target_*` constructions (Sections III and IV). The lone genuinely-independent check is Section I.3b (`mathematica_audit.wl:128-139`), which builds the spherical-Laplacian eigenvalue test from Mathematica's built-in `SphericalHarmonicY` and confirms eigenvalue `-l(l+1)` for `l=0` and `l=2`. The Section II comment (lines 143-145) explicitly self-discloses parallel-checking: *"this section is an intentional parallel check rather than an independent derivation. The two engines agree here as a sanity cross-check only."*

Borderline `mathematica_transliteration`. The Mathematica script does use `EulerEquations` from `VariationalMethods` (line 170) and `VariationalD` (lines 208-209), which are genuinely Mathematica-idiomatic operators, not hand-rolled mirrors of SymPy. The earlier v1 audit's `mathematica_transliteration` finding has already been partly addressed: the prior hand-rolled `eulerLagrange2D` is gone, replaced by `EulerEquations`/`VariationalD`. Given the self-disclosure plus the I.3b independent cross-check plus the use of native VariationalMethods operators, I decline to file a separate `mathematica_transliteration` finding in this v2 pass for this foundational definitional-bookkeeping stage.

## Engine cross-check

Both engines produce identical outputs in their saved transcripts: every `expect_zero`/`expectZero` returns 0; every Section's PASS line prints. Numerical residuals are nil because the checks are symbolic identities. The two engines therefore agree, but they agree on the same sign-convention quirks (F1 and F2), so engine agreement does not rescue paper alignment for those two checks.

## Verdict justification

The bookkeeping (real harmonics, mouth-average extraction, spherical-Laplacian eigenvalues, chain-rule confinement variation, `K_l` specialization, densitized-vs-weighted wall operator) faithfully exercises the paper's Output and Checks. Two boxed equations of the paper, however, are tested with sign conventions opposite to the paper's published form:

- **F1** (source-sign in eq:app-stage001-modal-wall-pde): a clear script-side sign error fixable in one character per engine; the script's coupling `−q·source` produces an equation with `−source` on the RHS while the paper writes `+source`. Routes to Codex.
- **F2** (gauge-term sign in eq:app-stage001-linear-maxwell): conditional on the unstated metric signature; either a script error or a paper-card omission. Routes to user for direction.

This unit is a checkpoint (high bar), so even foundational sign-convention disagreements are findings. `verdict: findings`, `stop_cold: null`, `paper_alignment: partial`.

The v1 (pre-paper-grounded) audit found one transliteration finding which has been substantively addressed in the current script (Mathematica now uses native VariationalMethods). The v2 paper-grounded re-audit surfaces two new findings (F1, F2) that v1 could not have caught because v1 did not read the paper.

## Self-test notes

I verified: (1) for F1, that SymPy's `euler_equations` returns `∂L/∂q − Σ_i d/dx_i ∂L/∂(∂_{x_i} q) = 0` (the standard EL form), so the source-sign reasoning is correct and the script's `−source_total` outcome is mechanical. (2) For F2, that the script's EL is internally consistent with its chosen Lagrangian signs; whether it matches the paper depends on the metric signature, which I cannot pin down from the allowed reading list — hence routing to user. (3) For the F1 fix, flipping `ldens − q*source` → `ldens + q*source` simultaneously flips `target_forced`, so the assertion identity will still hold and the equation verified will then match the paper's boxed RHS. (4) No `assert_nonzero`-style trivial-derivative trap exists in this script — every assertion is a zero-equality of well-formed differential expressions. (5) Harmonic average / orthogonality coverage: the script tests `cross(Y00, Y20)` and `average(Y_{2m})` for all five P2 modes; since `Y00` is a constant, `cross(Y00, Y_{2m}) ∝ average(Y_{2m}) = 0`, so orthogonality of `Y00` against the entire real P2 set is covered. (6) Paper round-trip on F1: the proposed fix swaps one sign each in the SymPy and Mathematica `ldens_forced` and `target_forced` lines, introducing no new literals or constants, so no new paper_misalignment can be introduced by the patch.
