---
unit_id: 158
batch: IV.6
created_at: 2026-05-27T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 158

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

---

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_158.tex:21-25` quote:
  ```
  \stagefield{Checks}{\begin{verificationchecklist}
  \item Check deviations are taken about the renormalized co-evolving canonical point.
  \item Check even-preservation constraints are imposed before reading the remaining odd defect.
  \item Check tangent motion on the parent compensation family gives \(\delta_\perp=0\).
  \end{verificationchecklist}}
  ```
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_158.tex:13` quote: "the transport of deviations into \(\Delta_Q\), D/N similarity slippage, and the final normal coordinate \(\delta_\perp\)."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py` — no assertion touches "even-preservation," "odd defect," "tangent motion," or "\(\delta_\perp\)".
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl` — same; transliteration of the sympy script's coverage.

**Notes side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage158_linear_defect_transport.md` — derives only the four boxed linearizations (`delta R`, `delta M_q`, `delta Pi`, `Delta_Q`). Contains a discussion of the compensated hybrid branch (§6.C) where the even branch stays canonical, but no `\delta_\perp` computation and no tangent-motion test.

## Resolve before fix_loop

The paper card promises three Checks; the scripts implement only the first (and only partially: linearization is taken about `g_* = r - \sqrt{1+r^2}/2`, the Family-1 canonical g-coordinate). The notes file itself does not contain a `\delta_\perp` computation or an even-preservation check. Which side is authoritative?

Possible directions (the user picks one):

- (a) **Paper card is authoritative; scripts under-cover.** Add new assertions to both scripts that test (a-2) even-preservation under the compensated-hybrid branch (`\rho_R = 4\sigma_W, \kappa_W = 1/3` per notes §6.C, with the resulting `\chi_Q^{hyb} = (1 - 9\sigma_W\gamma_W)/(1 - \sigma_W)`); and (a-3) a definition of `\delta_\perp` (presumably the perpendicular component of the deviation vector relative to the parent compensation family's tangent direction) and a check that tangent-aligned deformations give `\delta_\perp = 0`. The notes file would also need expansion to derive these so the script has anchors to cite. This is a real workload — three to four new boxed identities — and requires the user to specify the operational definitions of `\delta_\perp` and "parent compensation family" first.

- (b) **Notes/scripts are authoritative; paper card over-promises.** Trim items 2 and 3 from the `Checks` list in `stage_158.tex:23-24` so the card matches what notes and scripts actually verify. Also remove or soften the "the final normal coordinate \(\delta_\perp\)" phrase in `stage_158.tex:13` since this stage's notes do not compute that coordinate. No script changes needed.

- (c) **Both sides agree but the operational definitions are in an upstream stage.** Identify the upstream stage(s) where "even-preservation constraints" and "tangent motion / \(\delta_\perp\)" are operationally defined and verified, and update the paper card's Checks list to reference those upstream verifications by anchor (e.g., "carries even-preservation from stage_NNN") rather than claiming to verify them here. Script may need a single carry-forward comment but no new assertions.

The orchestrator will not invoke Codex on F1 until the user has chosen a direction. F2 and F3 below are independent script-side fixes that Codex may apply now.

---

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py:48-53`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl:41-45`

**Issue:** The `delta Ms law` check evaluates `(Sigma0 + dSigma0) - Sigma0 - dSigma0 = 0` — a pure rearrangement of the script's own definition `Ms = Sigma0 + dSigma0`. The check cannot fail and exercises no physical claim. The upstream identity it ostensibly tests (`M_s = \Sigma_0` from Stages 188-189) is hardcoded into the definition and never asserted.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py`, delete lines 48-53 inclusive (the `Sigma0, dSigma0, Rstar, dR` symbols block, the `Ms` definition, the `Mq`/`Mq_lin`/`Mq0` definitions, and the two `expect_zero` calls for `delta Ms law` and `delta Mq law`).

Wait — `Mq` and its check depend on these symbols. Do NOT delete the whole block. The precise edit is:

Before (sympy:48-54):
```python
Sigma0, dSigma0, Rstar, dR = sp.symbols("Sigma0 dSigma0 Rstar dR", real=True)
Ms = Sigma0 + dSigma0
Mq = -(Sigma0 + dSigma0) * (Rstar + dR)
Mq_lin = sp.expand(Mq).subs({dSigma0 * dR: 0})
Mq0 = -Sigma0 * Rstar
expect_zero("delta Ms law", Ms - Sigma0 - dSigma0)
expect_zero("delta Mq law", (Mq_lin - Mq0) - (-Rstar * dSigma0 - Sigma0 * dR))
```

After:
```python
Sigma0, dSigma0, Rstar, dR = sp.symbols("Sigma0 dSigma0 Rstar dR", real=True)
Mq = -(Sigma0 + dSigma0) * (Rstar + dR)
Mq_lin = sp.expand(Mq).subs({dSigma0 * dR: 0})
Mq0 = -Sigma0 * Rstar
expect_zero("delta Mq law", (Mq_lin - Mq0) - (-Rstar * dSigma0 - Sigma0 * dR))
```

(Remove lines 49 and 53.)

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl`, mirror the deletion at lines 41 and 45:

Before (wl:38-46):
```
Clear[sigma0, dSigma0, rStar, dR];
$Assumptions = Element[{sigma0, dSigma0, rStar, dR}, Reals];

mS = sigma0 + dSigma0;
mQ = -(sigma0 + dSigma0)*(rStar + dR);
mQLin = Expand[mQ] /. dSigma0*dR -> 0;
mQ0 = -sigma0*rStar;
expectZero["delta Ms law", mS - sigma0 - dSigma0];
expectZero["delta Mq law", (mQLin - mQ0) - (-rStar*dSigma0 - sigma0*dR)];
```

After:
```
Clear[sigma0, dSigma0, rStar, dR];
$Assumptions = Element[{sigma0, dSigma0, rStar, dR}, Reals];

mQ = -(sigma0 + dSigma0)*(rStar + dR);
mQLin = Expand[mQ] /. dSigma0*dR -> 0;
mQ0 = -sigma0*rStar;
expectZero["delta Mq law", (mQLin - mQ0) - (-rStar*dSigma0 - sigma0*dR)];
```

(Remove lines 41 and 45.)

**Verification:**

After Codex applies, the verifier will run `redteam exec-sympy 158` and `redteam exec-mathematica 158` and confirm:
1. Both transcripts no longer contain a `delta Ms law` line.
2. The remaining checks (`linear delta R law`, `delta Mq law`, `delta Pi law`, `linear Delta_Q law`, plus F3's two new checks) all show `= 0` and PASS.
3. Both scripts exit 0.

---

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py` (add new assertions after the existing `delta Pi law` check, before the numerical coefficients banner)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl` (mirror)

**Issue:** The script verifies four standalone linearizations but never asserts the two composed boxed identities from `notes/stages/moving_throat_pde_stage158_linear_defect_transport.md` §3-§4 that are this stage's actual carry-forward outputs:

```
\delta M_q = -(1/4)\delta\Sigma_0 + (\Sigma_0^{can}/\sqrt{1+r_{F1}^2})\delta g + O(2)
\delta \Pi = (1 - S_{can}/4)\delta\Sigma_0 - (\Sigma_0^{can}/4)\delta S + (\Sigma_0^{can} S_{can}/\sqrt{1+r_{F1}^2})\delta g + O(2)
```

The numerical-coefficients block prints values that match the notes, but those are independent literal arithmetic, not derived from the linearizations via assertion.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py`, after the `delta Pi law` check (currently line 64) and before the section-3 chi block (currently line 66), insert:

```python
# ---------------------------------------------------------------------------
# 3b. Composed boxed identities (notes/stages/moving_throat_pde_stage158_*.md §3-§4)
# ---------------------------------------------------------------------------
dg_sym, r_sym = sp.symbols("dg_sym r_sym", real=True, positive=True)
dR_from_dg = -dg_sym / sp.sqrt(1 + r_sym**2)

# Composed delta M_q in (dSigma0, dg) with R_* = 1/4
dMq_composed = -sp.Rational(1, 4) * dSigma0 - Sigma0 * dR_from_dg
dMq_boxed = -sp.Rational(1, 4) * dSigma0 + Sigma0 / sp.sqrt(1 + r_sym**2) * dg_sym
expect_zero("composed delta Mq law", sp.expand(dMq_composed - dMq_boxed))

# Composed delta Pi in (dSigma0, dS, dg) with R_* = 1/4, S_* = Sstar symbolic
dPi_composed = (1 - sp.Rational(1, 4) * Sstar) * dSigma0 \
    - Sigma0 * (sp.Rational(1, 4) * dS + Sstar * dR_from_dg)
dPi_boxed = (1 - Sstar / 4) * dSigma0 \
    - (Sigma0 / 4) * dS \
    + (Sigma0 * Sstar) / sp.sqrt(1 + r_sym**2) * dg_sym
expect_zero("composed delta Pi law", sp.expand(dPi_composed - dPi_boxed))
```

Note: `dg_sym` and `r_sym` are fresh symbols to avoid colliding with the section-1 `dg` and `r` symbols (which were used to derive `delta R`). The boxed identities use only the carried relation `dR = -dg/sqrt(1+r^2)` and the canonical value `R_* = 1/4`; `S_*` remains symbolic as `Sstar`.

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl`, after line 55 (the `delta Pi law` expectZero) and before line 57 (the `Clear[eps, ...]`), insert:

```mathematica
Clear[dgSym, rSym];
$Assumptions = Element[{sigma0, dSigma0, sStar, dS, dgSym, rSym}, Reals] && rSym > 0;
dRFromDg = -dgSym/Sqrt[1 + rSym^2];

dMqComposed = -(1/4)*dSigma0 - sigma0*dRFromDg;
dMqBoxed = -(1/4)*dSigma0 + (sigma0/Sqrt[1 + rSym^2])*dgSym;
expectZero["composed delta Mq law", Expand[dMqComposed - dMqBoxed]];

dPiComposed = (1 - sStar/4)*dSigma0 - sigma0*((1/4)*dS + sStar*dRFromDg);
dPiBoxed = (1 - sStar/4)*dSigma0 - (sigma0/4)*dS + (sigma0*sStar/Sqrt[1 + rSym^2])*dgSym;
expectZero["composed delta Pi law", Expand[dPiComposed - dPiBoxed]];
```

**Claim manifest:**

- M1: `\delta M_q = -(1/4)\delta\Sigma_0 + (\Sigma_0/\sqrt{1+r^2})\delta g` — symbolic identity with `R_* = 1/4` substituted, `\Sigma_0` and `r` left symbolic for canonical-point evaluation. (Sympy: `composed delta Mq law`. Mathematica: same.)
- M2: `\delta \Pi = (1 - S_*/4)\delta\Sigma_0 - (\Sigma_0/4)\delta S + (\Sigma_0 S_*/\sqrt{1+r^2})\delta g` — symbolic identity with `R_* = 1/4`, `S_*` symbolic. (Sympy: `composed delta Pi law`. Mathematica: same.)

**Self-test confirmation (per auditor instructions step 320):**

- Variable independence: `dMq_composed - dMq_boxed = -(1/4)*dSigma0 - Sigma0 * (-dg/sqrt(1+r^2)) - (-(1/4)*dSigma0 + Sigma0/sqrt(1+r^2)*dg) = -(1/4)*dSigma0 + Sigma0*dg/sqrt(1+r^2) + (1/4)*dSigma0 - Sigma0*dg/sqrt(1+r^2) = 0`. Holds symbolically. ✓
- `dPi_composed - dPi_boxed`: expand both — composed gives `(1 - Sstar/4)*dSigma0 - Sigma0/4*dS - Sigma0*Sstar*(-dg/sqrt(1+r^2)) = (1 - Sstar/4)*dSigma0 - Sigma0/4*dS + Sigma0*Sstar*dg/sqrt(1+r^2)`; boxed gives the same. ✓
- Symmetry/parity: no integrals; n/a.
- Trivial-case pre-check: at `dg = dSigma0 = dS = 0`, both composed and boxed are zero — assertion holds trivially. With `dg = 1, dSigma0 = 0, dS = 0, sigma0 = 4, Sstar = 0.5, r = 1`: composed `dMq = -4*(-1/sqrt(2)) = 4/sqrt(2)`; boxed `dMq = 4/sqrt(2)`. Match. With deliberately-wrong sign on the dg coefficient in `dMq_boxed`, mismatch is `8/sqrt(2)` — non-zero — finding fails. Non-tautological. ✓
- Path specifications: edits are to existing files; n/a.
- Paper round-trip: the boxed forms are quoted verbatim from `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage158_linear_defect_transport.md` §3 and §4. No new constants introduced; `R_* = 1/4` is the canonical value stated in notes line 41. ✓

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 158` and `redteam exec-mathematica 158` and confirm:
1. Both transcripts contain new lines `composed delta Mq law = 0` and `composed delta Pi law = 0` (and `PASS:` for the Mathematica versions).
2. Both scripts exit 0.
