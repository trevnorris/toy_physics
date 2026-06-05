---
unit_id: 024
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage024_overlap_isotropy.md]
  paper_appendix: present
---

# Audit unit 024 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_024.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage024_overlap_isotropy.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.txt`

## What the paper claims

Stage 024 takes the grouped Stage-023 coefficients and gives them their first explicit overlap-integral interpretation. `\stagefield{Output}` reads: "Stage~024 outputs the real-STF normalization \eqref{eq:app-stage024-Y-orthonormal}, the angular source-map identity \eqref{eq:app-stage024-source-map}, the isotropic overlap formulas \eqref{eq:app-stage024-B-moments}--\eqref{eq:app-stage024-D-isotropic}, the grouped isotropy collapse \eqref{eq:app-stage024-isotropic-collapse}, the weak-axisymmetric signature \eqref{eq:app-stage024-lambda-signature}, the line \(b=3a\), and the normalization transport coefficient \eqref{eq:app-stage024-PA-transport}." Distinct deliverables: (1) orthonormality `\int Y_A Y_B d\Omega = \delta_{AB}`; (2) angular source map `S_A^{port}=S_A`, hence `\hat m_{ang}=1`; (3) isotropic overlap moments `B_n`, `Z_n`, `N_n` and the conservative operator `D(\omega)`; (4) `O(3)` collapse `D_{20,n}=D_{21,n}=D_{22,n}=D_n`, same for `N`; (5) the axisymmetric triple-overlap matrix `M^{(20)}=\frac{\sqrt5}{7\sqrt\pi}\operatorname{diag}(1,1/2,1/2,-1,-1)` and grouped weights `\lambda=(1,1/2,-1)`; (6) the defect line `b=3a`; (7) the transport coefficient `P_1=(N_1 D_0 - N_0 D_1)/D_0^2`. The notes add detail, in particular the sixth-moment identity used to produce `M^{(20)}` and `\kappa_*=\sqrt5/(7\sqrt\pi)`. This is a checkpoint stage (`is_checkpoint: True`), so both engines are required and alignment must be exact.

## What the script claims to verify

Both scripts (SymPy + Mathematica) verify, section by section: Section I — the 5×5 Gram matrix of the normalized real-STF harmonics equals the identity (orthonormality), and `gram·svec = svec` (source-map identity, `\hat m_{ang}=1`). Section II — the Stage-022 grouped `(xbar,a,b)` decomposition under `(1,2,2)/5` weighting (witness substitutions plus a falsifiable round-trip reassembly). Section III — `B_n` from a two-mode BdG response series, then `Z_n`/`N_n`/`D_n` derived from a per-pair conservative 2×2 matrix inverse (a genuine physics anchor) and matched to the paper's closed forms. Section III.5 — per-lane `D`/`N` collapse under `O(3)` substitution plus a lane-breaking witness proving the collapse check is non-tautological. Section IV — the exact Y20 triple-overlap matrix equals `\frac{\sqrt5}{7\sqrt\pi}\operatorname{diag}(1,1/2,1/2,-1,-1)` and the grouped `b=3a` law. Section V — the first-order transport law for `P_A=N_A/D_A` and its grouped defects `b_P=3a_P`. The verdict applies to exactly these checks.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) `\int Y_A Y_B = \delta_{AB}` (eq:Y-orthonormal) | I.1 `expect_zero("Gram - I5", ...)` (py:164, wl:76) | match |
| (2) source map `S_A^{port}=S_A`, `\hat m_{ang}=1` (eq:source-map, eq:mhat-factor) | I.2 `projected - svec` (py:171, wl:81) | match |
| (3) `B_n`,`Z_n`,`N_n`,`D(\omega)` isotropic formulas (eq:B-moments..eq:D-isotropic) | III.1/III.2 series + matrix-inverse anchor (py:243-315, wl:119-183) | match |
| (4) collapse `D_{20,n}=D_{21,n}=D_{22,n}`, same `N` (eq:isotropic-collapse) | III.5.1/III.5.2 (py:363-370, wl:213-218) | match |
| (5) `M^{(20)}=\frac{\sqrt5}{7\sqrt\pi}\operatorname{diag}(1,1/2,1/2,-1,-1)`, `\lambda=(1,1/2,-1)` (eq:axisym-matrix, eq:lambda-signature) | IV.1 `M - M_target`; IV.2 weights (py:404,409-411, wl:236,240-243) | match |
| (6) `b=3a` (eq:b-3a) | IV.2 `b_x - 3 a_x` (py:425, wl:250) | match |
| (7) `P_1=(N_1 D_0 - N_0 D_1)/D_0^2` (eq:PA-transport), `b_P=3a_P` | V `P20/P21/P22` + grouped defects (py:458-470, wl:268-290) | match |
| Sixth-moment normalization `4\pi/105` (engine of #5, stated in notes as `4 pi / 122`) | py:128 `I6 = 4*pi*s/105`; wl:43 direct `Integrate` | mismatch (notes vs script) |

`paper_alignment: partial` — every boxed `\stagefield{Output}` deliverable is faithfully exercised by both engines, but the notes state a sixth-moment normalization factor that contradicts the script (and the notes' own downstream result).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 164 | `expect_zero(Gram - I5)` | (1) orthonormality | yes |
| A2 | sympy | 171 | `expect_zero(projected - svec)` | (2) source map | yes |
| A3 | sympy | 196-213 | grouped witnesses + round-trip reassembly | (4)/grouped coords | yes (round-trip is falsifiable) |
| A4 | sympy | 243-245 | `B0/B2/B4 - sum formula` (from series) | (3) B-moments | yes |
| A5 | sympy | 272-279 | matrix-inverse `Z/N` vs target (physics anchor) | (3) Z/N derivation | yes (non-tautological anchor) |
| A6 | sympy | 304-315 | `Z0..D4 - closed form` | (3) Z/N/D moments | yes (anchored by A5) |
| A7 | sympy | 363-370 | per-lane `D/N` collapse under iso subs | (4) collapse | yes |
| A8 | sympy | 383-384 | lane-breaking witness nonzero | (4) non-triviality guard | yes |
| A9 | sympy | 404 | `M - M_target` | (5) triple-overlap matrix | yes |
| A10 | sympy | 422-425 | `xbar/a/b` + `b - 3a` | (5)/(6) signature, `b=3a` | yes |
| A11 | sympy | 458-470 | `P20/P21/P22`, grouped defects | (7) transport | yes |
| B1 | mathematica | 76 | `Gram - I5` | (1) | yes |
| B2 | mathematica | 81 | `projected - svec` | (2) | yes |
| B3 | mathematica | 89-103 | grouped witnesses + reassembly | (4)/grouped | yes |
| B4 | mathematica | 119-121 | `B0/B2/B4 - sum` | (3) | yes |
| B5 | mathematica | 157-158 | matrix-inverse `Z/N` vs paper rational | (3) anchor | yes |
| B6 | mathematica | 175-183 | `Z0..D4 - closed form` | (3) | yes (anchored by B5) |
| B7 | mathematica | 213-218 | per-lane collapse | (4) | yes |
| B8 | mathematica | 222-223 | lane-breaking witness | (4) guard | yes |
| B9 | mathematica | 236 | `M - M_target` | (5) | yes |
| B10 | mathematica | 247-250 | `xbar/a/b`, `b - 3a` | (5)/(6) | yes |
| B11 | mathematica | 268-290 | quotient-rule transport + defects | (7) | yes |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** value_mismatch
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage024_overlap_isotropy.md:213`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:128`

**What's wrong:**
The notes state the exact sixth moment of the unit sphere with prefactor `4 pi / 122`:

> notes:213 — `` `int dOmega n_i n_j n_k n_l n_m n_n = (4 pi / 122) sum_pairings delta delta delta` ``

The SymPy script uses prefactor `4*pi/105`:

> py:128 — `return sp.simplify(4 * pi * s / 105)`

(The Mathematica engine sidesteps the explicit prefactor by integrating `n[i]...n[nn] Sin[theta]` directly, wl:42-46, and reproduces the same `M^{(20)}` and `\kappa_*`.) The correct mathematical value is `4\pi/105`: the script's `4\pi/105` reproduces the boxed downstream result `M^{(20)}=\frac{\sqrt5}{7\sqrt\pi}\operatorname{diag}(1,1/2,1/2,-1,-1)` (sympy output lines 126-145, mathematica output lines 85-86) and `\kappa_* = \sqrt5/(7\sqrt\pi)`. The notes' `4\pi/122` is internally inconsistent with the notes' OWN downstream line:

> notes:217 — `` `M^(20) = kappa_* diag(1, 1/2, 1/2, -1, -1)` `` with notes:221 `` `kappa_* = sqrt(5) / (7 sqrt(pi))` ``

If the prefactor were `4\pi/122`, `\kappa_*` would be rescaled by `105/122` and would not equal `\sqrt5/(7\sqrt\pi)`. So the `122` cannot be reconciled with either the script (105) or the notes' own boxed result. The `.tex` card (eq:app-stage024-fourth-moment is `4\pi/15`, and the sixth-moment paragraph, lines 160-166, refers to "the sixth moment" without stating the prefactor) is silent on this number, so the discrepancy lives between the notes and the script.

**Why this matters:**
This is the load-bearing intermediate that produces deliverable (5), the axisymmetric splitting matrix. A reader following the notes would compute the wrong `\kappa_*` and would not reproduce the boxed signature. It also undermines the notes' stated derivation chain (the quoted `4\pi/122 -> \kappa_*=\sqrt5/(7\sqrt\pi)` step does not actually hold).

**Required change:**
Per the paper_misalignment protocol, Codex applies nothing. Route to the user (see `## Resolve before fix_loop` in the directive). The discrepancy is almost certainly a notes-side typo (`122` -> `105`), but direction is the user's call; the script and the engines already agree on `105`/the direct integral.

**Verification:**
After the user picks a direction: if the notes are corrected to `4\pi/105`, the notes line 213 prefactor matches py:128 and the notes' own `\kappa_*` line 221 becomes consistent. No script change is needed (the engines already use the correct value and exit 0).

### F2 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.txt` (mtime 2026-05-26 00:38:42)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.txt` (mtime 2026-05-26 00:39:09)

**What's wrong:**
Both saved-output `.txt` files are OLDER than their scripts:
- sympy script mtime 2026-06-03 15:59:11 > sympy output mtime 2026-05-26 00:38:42
- mathematica script mtime 2026-06-03 15:59:11 > mathematica output mtime 2026-05-26 00:39:09

So the committed transcripts do not necessarily reflect the current script state. Reading both transcripts, every check still reads `= 0` / `PASS`, and the content matches what the current scripts would emit (all assertions are still present and the algebra is unchanged), so this is informational rather than a substantive disagreement.

**Why this matters:**
Freshness only; the captured output is consistent with the current source. The orchestrator's independent re-run will regenerate fresh transcripts.

**Required change:**
Re-run both scripts to refresh the saved outputs (the orchestrator does this on every batch). No source edit required.

**Verification:**
After re-run, output `.txt` mtimes exceed the script mtimes and content is unchanged (all PASS / `= 0`).

## Independent-derivation check (Mathematica)

The Mathematica script is NOT a transliteration of the SymPy script; the two engines use genuinely different machinery for the load-bearing steps:

1. **Sphere moments.** SymPy builds `I4`/`I6` from an explicit combinatorial pairing sum times `4\pi/15` and `4\pi/105` (py:109-128, with a hand-written `pairings` generator py:71-81). Mathematica instead does the literal surface integral `Integrate[n[i]...n[nn] Sin[theta], {theta,0,Pi},{phi,0,2Pi}]` with `n[1]=Sin[theta]Cos[phi]` etc. (wl:31-46). These are independent derivations of the same moments — the `.wl` never references a pairing sum or the `105` constant.

2. **Z/N response.** Both reach the closed forms, but SymPy's series anchor (py:265-279) and Mathematica's (wl:132-158) are independently coded; the `.wl` is not a line-by-line echo (e.g. wl uses `First[First[{coupling}.Inverse[mPair].Transpose[{coupling}]]]`).

3. **First-order transport.** SymPy uses `sp.series(expr, eps, 0, 2)` (py:448), Mathematica deliberately uses the quotient rule `D[pAFull[lam], eps] /. eps -> 0` (wl:259-261) "independent of SymPy's series-helper" (wl:256-258). Genuinely independent route.

No `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce identical final results at the level they claim:
- Gram - I5 = zero 5×5 (sympy out 19-28; mathematica out 9-10).
- M^(20) = `\sqrt5/(7\sqrt\pi) diag(1,1/2,1/2,-1,-1)` and `M - M_target = 0` (sympy out 126-155; mathematica out 85-86).
- All `B/Z/N/D` coefficient residuals = 0 (sympy out 73-93; mathematica out 35-62).
- All transport residuals = 0 (sympy out 176-186; mathematica out 99-114).
The lane-breaking witness coefficients are written in different but algebraically-equal forms (sympy out 116; mathematica out 79), as expected from the two different parametrizations — not a disagreement. `engines_agree: true`.

## Verdict justification

`verdict: findings` with two findings, neither stop-cold. Every boxed `\stagefield{Output}` deliverable is faithfully and non-tautologically exercised by BOTH engines, anchored where it could have been tautological (the matrix-inverse anchor in Section III converts the Z/N/D closed-form identities into real physics checks; the round-trip reassembly in Section II and the lane-breaking witness in Section III.5 guard against trivial passes). Attacks tried that failed: (a) checked whether the Gram/source-map checks are tautological — they are not, the quadrupole overlaps are computed from independent sphere moments; (b) checked whether the Z/N closed forms are self-referential — the matrix-inverse anchor breaks that; (c) checked the sixth-moment parity/normalization that produces `M^{(20)}` — the script's `4\pi/105` is correct and reproduces `\kappa_*`; (d) checked engine independence — confirmed (direct integral vs pairing sum). The two findings are: F1, a notes-side `value_mismatch` (sixth-moment prefactor written as `122` in notes:213, contradicting the script's `105` and the notes' own `\kappa_*`), routed to the user; and F2, `stale_output` (both transcripts predate their scripts but their content still matches), informational. I confirm I read the stage card, the notes file, and the part-02 appendix row (line 38) / summary (line 122), and the script's verified claims match the card's boxed Output exactly.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit against the `.tex` card and `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `\int Y_A Y_B d\Omega = \delta_{AB}` (Gram = I5) | py:160-164, wl:75-76; sympy out 9-28 | .tex:37-39 (eq:Y-orthonormal); notes:51-52 | MATCH |
| source map `S_A^{port}=S_A` ⇒ `\hat m_{ang}=1` | py:169-172, wl:81; sympy out 33-43 | .tex:52,59 (eq:source-map, eq:mhat-factor); notes:68,78 | MATCH |
| `B_0=\sum C_\alpha^2/\varpi^2`, `B_2`(/⁴), `B_4`(/⁶) | py:243-245, wl:119-121; sympy out 73-75 | .tex:93-98 (eq:B-moments); notes:132-136 | MATCH |
| `Z_0,Z_2,Z_4` closed forms | py:304-306, wl:175-177; sympy out 85-87 | .tex:118-120 (eq:ZN-moments); notes:155-159 | MATCH |
| `N_0,N_2,N_4` closed forms | py:307-312, wl:178-180; sympy out 88-90 | .tex:121-123; notes:163-167 | MATCH |
| `D_0=K-B_0-Z_0`, `D_2=-(M+B_2+Z_2)`, `D_4=-(B_4+Z_4)` | py:313-315, wl:181-183; sympy out 91-93 | .tex:131-134 (eq:D-isotropic); notes:175-179 | MATCH |
| collapse `D_{20,n}=D_{21,n}=D_{22,n}`, `N` same | py:363-370, wl:213-218; sympy out 102-111 | .tex:142-144 (eq:isotropic-collapse); notes:183-185 | MATCH |
| `M^{(20)}=\frac{\sqrt5}{7\sqrt\pi}\operatorname{diag}(1,1/2,1/2,-1,-1)` | py:402-404, wl:234-236; sympy out 126-155 | .tex:164 (eq:axisym-matrix); notes:217 | MATCH |
| `\kappa_* = \sqrt5/(7\sqrt\pi)` | py:402, wl:234; (embedded in M) | .tex:164 (prefactor); notes:221 | MATCH |
| splitting weights `\lambda_{20}=1,\lambda_{21}=1/2,\lambda_{22}=-1` | py:409-411,427-430, wl:240-243; sympy out 165-167 | .tex:171-175 (eq:lambda-signature); notes:231-235 | MATCH |
| grouped `xbar=x^{(0)}, a_x=\varepsilon x^{(1)}/4, b_x=3\varepsilon x^{(1)}/4, b=3a` | py:422-425, wl:247-250; sympy out 160-163 | .tex:190-196 (eq:b-3a); notes:247-255 | MATCH |
| `P_1=(N_1 D_0 - N_0 D_1)/D_0^2` | py:455, wl:264; sympy out 176-178 | .tex:212 (eq:PA-transport); notes:295 | MATCH |
| `a_P=\varepsilon P_1/4, b_P=3\varepsilon P_1/4, b_P=3a_P` | py:467-470, wl:288-290; sympy out 183-186 | .tex:218-222 (eq:P-defects); notes:299-305 | MATCH |
| fourth-moment prefactor `4\pi/15` | py:117, wl (direct integral); | .tex:32 (eq:fourth-moment); notes:48 | MATCH |
| **sixth-moment prefactor `4\pi/105`** | py:128 (`4*pi*s/105`); wl:42-46 (direct integral) | notes:213 states `4 pi / 122`; .tex silent | **MISMATCH** → F1 (value_mismatch) |

INTERNAL scaffolding (accounted for, no finding): the raw STF basis matrices `E20/E21c/E21s/E22c/E22s` and `NORM=\sqrt{15/(8\pi)}`; the per-pair 2×2 matrix `M_pair`/`mPair` and its inverse (anchor machinery); the two-mode/one-pair witness substitutions and the `(p,q,rr)` round-trip; the lane-breaking witness `delta`-coefficient; pass/fail flags and `= 0` residuals.

reconciliation: 16 deliverable values checked, 1 misaligned (sixth-moment prefactor: notes `122` vs script `105`). All other deliverables MATCH.

## Self-test notes

Trap checks: (1) Variable independence — the only derivatives are `sp.series`/`D[...,eps]` over `eps`/`omega`, and the expressions genuinely depend on those variables (e.g. `pAFull` depends on `eps`; the `B/Z/N` responses depend on `omega`), so no identically-zero derivative trap. (2) Symmetry/parity — the sphere integrals are even-order (4th, 6th) moments that are correctly nonvanishing; the script's pairing-sum/direct-integral both give nonzero, consistent results, and the odd-lane structure of `M^{(20)}` (signature `(1,1/2,1/2,-1,-1)`) is reproduced. (3) Trivial-case — the matrix-inverse anchor (Section III) and the lane-breaking witness (Section III.5) both confirm the collapse/closed-form checks are non-tautological. (4)/(5) No new script is prescribed (F1 routes to user, F2 is a re-run), so no path-spec or paper-round-trip trap applies; I confirmed the only numeric mismatch is notes-side and the script already uses the mathematically correct `4\pi/105`.
