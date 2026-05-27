---
unit_id: 096
batch: IV.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage096_geometry_lane_check_verdict.md]
  paper_appendix: present
---

# Audit unit 096 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_096.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage096_geometry_lane_check_verdict.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_096}` row at line 1226 is in scope)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.txt`

## What the paper claims

The stage card (`stage_096.tex` lines 7-25) calls this stage a "geometry-lane firewall ledger step" and quotes the conservative carrier as `\widehat Y_Q^{\rm cons}=3/4+(1/4)/(1-\omega^2/\Omega_Q^2)`. The narrative verdict is "the natural isotropic branch is conservatively clean; only genuine symmetry breaking or a second dynamic l=2 geometry pole can reopen the issue." The card's explicit Checks (lines 21-25) are (i) static limit `eps_2=eps_4=0` returns `c_pole=1/4`, (ii) `l=0`/`l=2` orthogonality before applying the geometry firewall, (iii) "any support/source success statement still carries the minimal-module hypothesis". The notes file lines 30-53 spell out the bottom-line deliverables: `eps_2=eps_4=0`, `c_pole=1/4`, `c_geom=3/4`, `Yhat_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)`, `rho_alpha=4/3`, `zeta_req=1/3`. Notes line 26-28 attribute the contamination-zero step to "Stage 77 proves exactly that the scalar/geometry l=0 lane and the grouped real l=2 bundle are block diagonal in the quadratic wall theory." This stage is flagged as a checkpoint.

## What the script claims to verify

The SymPy script (lines 2-18) states it (1) verifies the l=0/l=2 orthogonality on the isotropic branch, (2) collapses the obstruction formula to `3/4 + 1/4`, and (3) reproduces `rho_alpha=4/3`, `zeta_req=1/3`. Its non-trivial assertions are: (a) `<Y00|Y2m> = 0` for the five grouped real `P_2` harmonics; (b) `(-Delta_{S^2})Y2m - 6 Y2m = 0`; (c) `<Y00|(-Delta)Y2m> = 0`; (d) several arithmetic identities downstream of the hardcoded `eps_2 = eps_4 = 0`. The Mathematica script mirrors this set of checks using `FullSimplify` and Mathematica's own `Integrate`/Laplace stencil.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side coverage |
|---|---|
| `l=0`/`l=2` orthogonality before firewall (Check ii) | **match** — SymPy lines 73-84 and `.wl` lines 47-55 integrate `<Y00\|Y2m>` and `<Y00\|(-Delta)Y2m>` over `S^2` and assert vanishing for all five m-modes; Laplace eigenvalue 6 verified independently. |
| Static limit `eps_2=eps_4=0` returns `c_pole=1/4` (Check i) | **mismatch** — script does NOT verify the implication "block-diagonality ⇒ eps_2=eps_4=0"; it sets `eps_2 = sp.Integer(0)`, `eps_4 = sp.Integer(0)` by assignment (sympy line 88-89; `.wl` line 57-58) and then derives `c_pole = (1+0)/(4(1+0)^2)`. The `c_pole - 1/4` check then reduces to `1/4 - 1/4 == 0` arithmetic. The orthogonality block is not algebraically linked to the assignment. |
| `Yhat_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)` | **partial** — given the hardcoded eps's, the algebraic identity (c_geom + c_pole/(1-omega^2/Omega_Q^2)) − (3/4 + (1/4)/(1-omega^2/Omega_Q^2)) = 0 is non-trivial as a rational-function identity, but it inherits the hardcoded inputs. |
| `rho_alpha = 4/3`, `zeta_req = 1/3` | **partial** — same: arithmetic substitution after the hardcoded eps's; `expect_zero("zeta_req - c_pole/c_geom", ...)` on sympy line 117 is purely tautological since `zeta_req = c_pole/c_geom` was just defined on line 101. |
| Minimal-module hypothesis carried through support/source statement (Check iii) | **missing** — no script-side check or even a printed annotation that the rho_alpha=4/3, zeta_req=1/3 conclusions are conditioned on the minimal-module hypothesis. |

`paper_alignment` set to `partial` — bottom-line numbers match the paper; verification structure does not fully exercise the paper's "Checks" list.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 82 (5x) | `expect_zero("<Y00|Y2m>", overlap)` | Check (ii) orthogonality | yes |
| A2 | sympy | 83 (5x) | `expect_zero("(-Delta)Y2m - 6 Y2m", ...)` | Laplace eigenvalue for ell=2 (orthogonality prerequisite) | yes |
| A3 | sympy | 84 (5x) | `expect_zero("<Y00|(-Delta)Y2m>", ...)` | Check (ii) orthogonality | yes |
| A4 | sympy | 90 | `expect_zero("eps_2", eps_2)` with `eps_2 = sp.Integer(0)` | Check (i) precondition | **no — tautological** |
| A5 | sympy | 91 | `expect_zero("eps_4", eps_4)` with `eps_4 = sp.Integer(0)` | Check (i) precondition | **no — tautological** |
| A6 | sympy | 109 | `expect_zero("c_pole - 1/4", ...)` | Check (i) target | partial — follows by hardcoded substitution |
| A7 | sympy | 110 | `expect_zero("c_geom - 3/4", ...)` | derived from A6 | partial — substitution |
| A8 | sympy | 111-114 | `Yhat_Q^cons - [3/4 + (1/4)/(1 - omega^2/Omega_Q^2)] == 0` | conservative carrier identity | yes (genuine rational-function identity, modulo hardcoded eps's) |
| A9 | sympy | 115 | `expect_zero("rho_alpha - 4/3", ...)` | rho_alpha deliverable | partial — substitution |
| A10 | sympy | 116 | `expect_zero("zeta_req - 1/3", ...)` | zeta_req deliverable | partial — substitution |
| A11 | sympy | 117 | `expect_zero("zeta_req - c_pole/c_geom", ...)` | none — restates definition | **no — tautological** (`zeta_req` literally defined as `c_pole/c_geom` on line 101) |
| M1 | mathematica | 50 (5x) | `expectZero["<Y00|Y2m>", ...]` | Check (ii) | yes |
| M2 | mathematica | 51 (5x) | `expectZero["(-Delta)Y2m - 6 Y2m", ...]` | Laplace eigenvalue | yes |
| M3 | mathematica | 52 (5x) | `expectZero["<Y00|(-Delta)Y2m>", ...]` | Check (ii) | yes |
| M4 | mathematica | 74 | `expectZero["eps2", eps2]` | precondition | **no — tautological** |
| M5 | mathematica | 75 | `expectZero["eps4", eps4]` | precondition | **no — tautological** |
| M6 | mathematica | 76 | `expectZero["c_pole - 1/4", ...]` | Check (i) | partial — substitution |
| M7 | mathematica | 77 | `expectZero["c_geom - 3/4", ...]` | derived | partial — substitution |
| M8 | mathematica | 78 | `Yhat_Q^cons - [3/4 + (1/4)/(1 - omega^2/Omega_Q^2)] == 0` | conservative carrier identity | yes |
| M9 | mathematica | 79 | `expectZero["rho_alpha - 4/3", ...]` | rho_alpha | partial — substitution |
| M10 | mathematica | 80 | `expectZero["zeta_req - 1/3", ...]` | zeta_req | partial — substitution |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py:88-91` and `:117`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.wl:57-58` and `:74-75`

**What's wrong:**
The script defines

```
eps_2 = sp.Integer(0)
eps_4 = sp.Integer(0)
expect_zero("eps_2", eps_2)
expect_zero("eps_4", eps_4)
```

at sympy lines 88-91 (mirror at `.wl` 57-58, 74-75). The assertions `eps_2 == 0` and `eps_4 == 0` are trivially true by construction — they cannot fail no matter what the orthogonality block computes. The output transcript ("eps_2 = 0", "eps_4 = 0") cannot reveal anything.

Additionally, sympy line 117 has

```
expect_zero("zeta_req - c_pole/c_geom", zeta_req - c_pole / c_geom)
```

but `zeta_req` was defined on line 101 as `sp.simplify(c_pole / c_geom)`. The check `(c_pole/c_geom) - (c_pole/c_geom) == 0` is algebraically guaranteed.

**Why this matters:**
This is a checkpoint stage (higher bar). The notes lines 26-33 attribute `eps_2 = eps_4 = 0` to "Stage 77 proves exactly that the scalar/geometry l=0 lane and the grouped real l=2 bundle are block diagonal". The orthogonality block A1/A3 (sympy 73-84; .wl 47-55) is the substance the verdict rests on, but the script never builds an algebraic object whose vanishing IS `eps_2`/`eps_4`. As written, the script's verdict on Check (i) ("static limit eps_2=eps_4=0 returns c_pole=1/4") is "we type 0 and arithmetic substitution gives 1/4", not "block diagonality forces the obstruction numbers to vanish".

**Required change:**
The `paper_misalignment`-direction question (whether to call this `paper_misalignment`) is closed: paper and script agree on the bottom line. The script-side fix is to (a) delete the trivially-true assertions on sympy lines 90-91 and 117, and `.wl` lines 74-75, since they convey no information; OR (b) construct `eps_2`/`eps_4` as the actual L2 inner products that the orthogonality block already evaluates (e.g., the `<Y00|(-Delta)Y2m>`-style integrals that contamination numbers would correspond to in the obstruction formula), then assert those integrals vanish. Option (b) is preferred for a checkpoint. Either way, also remove the tautological `expect_zero("zeta_req - c_pole/c_geom", ...)` on sympy line 117.

**Verification:**
After the fix: sympy output should no longer contain bare "eps_2 = 0" / "eps_4 = 0" PASS lines that are guaranteed by construction. Either those lines are gone (cleanup option a) or they correspond to integrals whose integrand is non-trivially constructed from the Y_2m harmonics (option b — the integral expression appears in the script body, not a literal `Integer(0)`).

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py:86-116`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.wl:57-80`

**What's wrong:**
The paper card's Checks list (stage_096.tex lines 21-25) has three items:
1. Static limit `eps_2=eps_4=0` returns `c_pole=1/4`.
2. `l=0` and `l=2` orthogonality before applying the geometry firewall.
3. Any support/source success statement still carries the minimal-module hypothesis.

Check (ii) is exercised by the orthogonality block (sympy 73-84). Check (i) is exercised only via hardcoded substitution (see F1). Check (iii) — that the carried `rho_alpha=4/3`/`zeta_req=1/3` statement preserves the minimal-module hypothesis — is **not exercised at all**, even as a printed annotation in the FINAL LEDGER. The script's `FINAL LEDGER` (sympy lines 119-126) just lists the numerical results without flagging that they are conditional on the minimal-module hypothesis. Likewise the Mathematica "Stage 096 Mathematica audit passed." line at `.wl:83`.

**Why this matters:**
For a checkpoint, the downstream-use card (stage_096.tex line 27) is explicit: "When cited downstream, the status tag above must be carried with the result; the card is a derivation ledger entry, not an unconditional actual-branch theorem." A passing script transcript should make the conditioning visible so that downstream stages that quote `rho_alpha=4/3`/`zeta_req=1/3` know what assumption is propagating with them.

**Required change:**
Append an explicit annotation block (printed, not asserted) at the end of both scripts that reproduces the minimal-module hypothesis condition tag. For sympy, after line 126 add a print block that states: "These results assume the Part III minimal isotropic module and the grouped real P_2 carrier; they are not unconditional actual-branch theorems." Mirror the same block at the end of the `.wl` script before `Exit[0]`. This is a printed annotation, not a new assertion; the script keeps exit code 0.

**Verification:**
Run-output of both engines should include a "Hypothesis carried:" or equivalent line in their FINAL LEDGER section. Diff the new sympy `.txt` vs old: a single block of new print lines flagging the minimal-module hypothesis appears.

### F3 — stale_output (informational)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py` (mtime May 11 11:56)
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.txt` (mtime May 11 12:45)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.wl` (mtime May 11 11:56)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.txt` (mtime May 11 13:04)

**What's wrong:**
Outputs are newer than scripts, so this is not technically stale; but two cosmetic issues survived into the captured outputs:
- Both transcripts contain the banner "STAGE 079 — GEOMETRY-LANE CHECK VERDICT" (sympy script line 53; `.wl` line 26). The unit is `096`, not `079`. The numeric stale label is also present in the saved outputs at lines 11 of each `.txt`.
- The sympy docstring (sympy lines 2-18) refers internally to "Stage 75 obstruction formula" and "Stage 77" — those are upstream references, fine to retain — but the section banner labels this audit as `STAGE 079`, which is incorrect.

This is **informational only** — the math is unchanged — but should be fixed for traceability since this is a checkpoint.

**Why this matters:**
Auditors and downstream readers reading the captured PASS transcript could be confused about which stage card the output anchors to. Notes file line 12 already shows past renumbering history (it references `stage147` even though the current path is `stage096`); the banner should not amplify that confusion.

**Required change:**
- sympy line 53: change `banner("STAGE 079 — GEOMETRY-LANE CHECK VERDICT")` to `banner("STAGE 096 — GEOMETRY-LANE CHECK VERDICT")`.
- `.wl` line 26: change `banner["STAGE 079 — GEOMETRY-LANE CHECK VERDICT"];` to `banner["STAGE 096 — GEOMETRY-LANE CHECK VERDICT"];`.

After Codex applies, the verifier re-runs both engines, producing fresh `.txt` files with the corrected banner.

**Verification:**
Re-run captured outputs show "STAGE 096 — GEOMETRY-LANE CHECK VERDICT" at the banner line.

## Independent-derivation check (Mathematica)

The `.wl` file is structurally close to the SymPy script — same variable names (`y00`, `y20`, `y21c`, ...), same per-mode loop structure, same arithmetic chain `cPole -> cGeom -> rhoAlpha -> zetaReq`, same assertion list in the same order. The Laplacian implementation `lapS2` and integral `dOmegaIntegral` use Mathematica's own `D[…]` and `Integrate[…]`, so the spherical-harmonic and Laplace identities are independently derived by each engine's primitive (sympy `sp.integrate`/`sp.diff` vs Mathematica `Integrate`/`D`). I would not call this `mathematica_transliteration` because the core integrals are evaluated by independent symbolic engines and Mathematica reports the result of its own integration; the per-mode structural mirror is fine for a verdict-checkpoint where both engines should yield the same outcome via independent quadrature.

What is shared between the two — and is a finding regardless of engine — is the hardcoded `eps2 = 0; eps4 = 0` assignment. Both engines accept the same handed-down constant. That is F1 above, not `mathematica_transliteration`.

## Engine cross-check

Both engines produce the same numerical outputs:
- sympy `c_pole = 1/4`, `c_geom = 3/4`, `Yhat_Q^cons(omega) = (Omega_Q**2 - 3*omega**2/4)/(Omega_Q**2 - omega**2)`, `rho_alpha = 4/3`, `zeta_req = 1/3` (output lines 45-49).
- mathematica `c_pole = 1/4`, `c_geom = 3/4`, `Yhat_Q^cons(omega) = (3 + (1 - omega^2/omegaQ^2)^(-1))/4`, `rho_alpha = 4/3`, `zeta_req = 1/3` (output lines 51-55).

The two simplified forms of `Yhat_Q^cons` are algebraically equal (multiply out to confirm), and both engines pass the explicit identity `Yhat_Q^cons - [3/4 + (1/4)/(1 - omega^2/Omega_Q^2)] == 0`. No `engine_disagreement`.

## Verdict justification

The orthogonality block (A1–A3 / M1–M3) is solid: 15 non-trivial spherical integrals across both engines, plus the Laplace eigenvalue identity at ell=2, all genuine algebra. The `Yhat_Q^cons` rational-function identity (A8 / M8) is also genuine. But the **link** between the orthogonality verdict and the obstruction-numbers-vanish statement (which is what makes the stage a "verdict") is asserted by typing zero into a variable, not derived. For a non-checkpoint stage this would be a low-severity insufficiency note; for a checkpoint with the explicit "Check static limit eps_2=eps_4=0 returns c_pole=1/4" item in the card, it is a real gap. Add to that the entirely-missing minimal-module-hypothesis annotation (Check iii) and the stale "STAGE 079" banner, and the verdict is `findings` with three medium/low items, no stop-cold.

## Self-test notes

Trap 1 (variable independence): I did not propose any new `sp.diff` / `D[]` in F2 (the only literal addition is a `print` block). For the cleanup option (a) in F1, no new derivative or integral is added. For the substantive option (b) in F1, I left the implementation open for the user to decide (because choosing the right integral representation for `eps_2`/`eps_4` would require reading Stage 77's scripts, which is out of scope for this unit's audit) — that is why F1's directive presents both options. Trap 2 (parity): n/a — no new integrals proposed. Trap 3 (trivial-case): substituting `eps_2=0, eps_4=0` into `c_pole = (1+eps_4)/(4(1+eps_2)^2)` yields `1/4`, consistent with the current behavior; F2's print block has no assertion so no trivial-case risk. Trap 4 (paths): the directive's targets name files in `scripts/` and `mathematica/` explicitly. Trap 5 (paper round-trip): F2's annotation language is taken from the paper card's own line 27 ("the card is a derivation ledger entry, not an unconditional actual-branch theorem"), so the annotation cannot introduce a new paper_misalignment.
