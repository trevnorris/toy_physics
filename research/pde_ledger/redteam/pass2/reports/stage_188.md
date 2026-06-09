---
unit_id: 188
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage188_branch_observables_completion.md]
  paper_appendix: present
---

# Audit unit 188 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_188.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage188_branch_observables_completion.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (row at line 107; overview at line 265; subsection at 881; anchor at 1466)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage188_branch_observables_completion_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage188_branch_observables_completion_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` (stage_188.tex:15): "Compiles direct branch observables to \((\Theta_1,\Xi_1,\mathcal R_1)\) and records equivalent observable packets." The notes (the authoritative detail layer) enumerate four deliverables: (1) the exact coefficient identity `A_{tr,*}=B_* C_{tr,*}`; (2) an exact invertible first-order compiler `C_obs->quot = diag(-1/C_{tr,*},1,1)` from the branch-observable packet `(δln R_tr, δln N_*, δln ε_η)` to the Stage-187 tangent quotient packet `(Σ_tr,Σ_nt,Σ_η)`; (3) an exact invertible compiler `C_obs->def` (lower-triangular, with rows `Θ_1=δln R_tr`, `Ξ_1=δln N_* − B_* δln R_tr`, `R_1=−ε_η*/(1−ε_η*)·δln ε_η − Ξ_1`) and its inverse `C_def->obs`; (4) the branch-observable zero-defect theorem `Θ_1=Ξ_1=R_1=0 ⟺ δln R_tr=δln N_*=δln ε_η=0`, which follows from invertibility of both compilers. The notes also give the complementary identity `δln(1−ε_η)=R_1+Ξ_1=−ε_η*/(1−ε_η*)·δln ε_η` and the factorization `C_obs->def = C_quot->def · C_obs->quot`. The card is `\StatusExactClosure`, checkpoint: False; the card itself states "Mathematica audit: none yet" (stage_188.tex:11) — but a `.wl` is in fact present (this is the V.3 dual-engine retrofit batch); the card's prose lag is a numbering/prose-sync item, not a script defect.

## What the script claims to verify

Both scripts verify, symbolically over `χ0_*,δU_*>0`, `0<ε_η*<1`: (I) the two coherent-branch coefficient identities `C_* = 1/C_{tr,*}` and `A_{tr,*} = B_* C_{tr,*}`; (II) that the observable→quotient compiler equals `diag(−C_*,1,1)`, has determinant `−C_*≠0`, and is invertible; (III) that the quotient→defect compiler reproduces the Stage-183/185 triangular normal form and yields the three componentwise defect relations; (IV) that the direct observable→defect compiler equals the lower-triangular `C_obs->def`, factorizes as `C_quot->def·C_obs->quot`, has determinant `−ε_η*/(1−ε_η*)≠0`; (V) that its inverse equals the stated `C_def->obs` and round-trips the generic observable packet; (VI) the complementary-observable drift identity; (VII) that the zero-defect equivalence is exercised on the GENERIC packet (not the zero vector) via bijectivity, with both determinants confirmed nonzero. The SymPy script builds every compiler matrix as a hand-typed literal and checks products/inverses/dets. The Mathematica script independently re-derives each compiler from finite branch-ratio first-log-drifts (`Series`/`Log`) and `Solve` of the branch-drift equation systems, then differentiates with `Outer[D,...]`, and cross-checks against hard-typed SymPy reference matrices.

## Paper ↔ script cross-check

| Paper/notes deliverable | Script-side check | Status |
|---|---|---|
| `A_{tr,*}=B_* C_{tr,*}` (notes §2) | py:63 / wl:95 `expect_zero(A_tr,* − B_* C_tr,*)` | match |
| `C_obs->quot = diag(−1/C_{tr,*},1,1)`, det≠0, invertible (notes §3) | py:74,82-88 / wl:114-129 | match |
| `C_quot->def` triangular normal form + 3 defect relations (notes §4-5) | py:93-115 / wl:144-167 | match |
| `C_obs->def` lower-triangular + factorization + det (notes §4) | py:118-144 / wl:189-209 | match |
| `C_def->obs` inverse + round-trip (notes §4) | py:147-164 / wl:224-241 | match |
| complementary `δln(1−ε_η)=R_1+Ξ_1` (notes §1.3) | py:167-170 / wl:245-249 | match |
| zero-defect iff via invertibility on GENERIC packet (notes §6) | py:172-196 / wl:251-271 | match |

Every notes deliverable has a faithful, non-tautological script-side check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `expect_zero(C_* − 1/C_tr,*)` | coeff identity (precursor) | yes |
| A2 | sympy | 63 | `expect_zero(A_tr,* − B_* C_tr,*)` | claim 1 (hinge identity) | yes |
| A3 | sympy | 87-88 | `C_quot->obs·C_obs->quot − I` (both orders) | claim 2 (invertibility) | yes |
| A4 | sympy | 110-115 | 3× defect-relation residuals | claim 3 | yes |
| A5 | sympy | 130 | `factorized − expected compiler` | claim 3 (factorization) | yes |
| A6 | sympy | 139-144 | 3× `Δ_def` componentwise | claim 3 | yes |
| A7 | sympy | 159-164 | inverse + round-trip | claim 3 (inverse) | yes |
| A8 | sympy | 170 | `(R+Ξ) − δln(1−ε_η)` | complementary identity | yes |
| A9 | sympy | 177-196 | bijection on generic packet + det-nonzero | claim 4 (iff) | yes |
| B1 | mathematica | 94-95 | `expectZero` coeff identities | claims 1 (precursor+hinge) | yes |
| B2 | mathematica | 126-129 | obs→quot derived by Solve+D, == SymPy, inverse | claim 2 | yes |
| B3 | mathematica | 153-167 | quot→def derived by Solve+D, == SymPy, 3 relations | claim 3 | yes |
| B4 | mathematica | 201-209 | obs→def derived by Solve(branch drifts)+D, factorization, == SymPy, det | claim 3 | yes |
| B5 | mathematica | 234-241 | inverse by re-Solve, == SymPy, == matrix-inverse, round-trip | claim 3 | yes |
| B6 | mathematica | 246-249 | complementary drift from `Log[(1−η e^{s·etaObs})/(1−η)]` | complementary identity | yes |
| B7 | mathematica | 253-271 | bijection on generic packet + det-nonzero | claim 4 (iff) | yes |

Every row traces to a paper-side deliverable. No orphaned scaffolding.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.txt:3` (and :272)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py:35` (and :198)

**What's wrong:**
The committed SymPy output `.txt` is stale relative to the current `.py`. The `.py` mtime is `Jun 3 15:59`; the `.txt` mtime is `Jun 1 11:23` (`.txt` older → stale). The content also disagrees: the saved output banner reads `STAGE 171 — BRANCH-OBSERVABLE COMPLETION...` (`.txt:3`) and `STAGE 171 LEDGER` (`.txt:272`), whereas the current script prints `STAGE 188 ...` (`.py:35`) and `STAGE 188 LEDGER` (`.py:198`). The captured transcript predates the stage-banner renumber fix. The body math in the stale `.txt` otherwise matches the current script's identities (all the same residuals print `= 0`), so this is a label-only staleness, not a math discrepancy. The Mathematica output (`.txt:3` `STAGE 188 ...`) is already current.

**Why this matters:**
A committed transcript with a wrong stage banner is a numbering/audit-trail hazard (the same `STAGE 171`-vs-`188` self-label class the project is mid-remediation on) and means the saved SymPy artifact does not reflect the current script. It is not a math defect.

**Required change:**
Re-run the SymPy script to refresh the saved output so the banner reads `STAGE 188` and the mtime post-dates the `.py`. No source edit is required for this finding (the `.py` already prints `STAGE 188`); the verifier's standard re-run will regenerate the `.txt`.

**Verification:**
After the verifier re-runs `python3 scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py`, the saved `.txt` line 3 reads `STAGE 188 — BRANCH-OBSERVABLE COMPLETION...`, the ledger banner reads `STAGE 188 LEDGER`, and the `.txt` mtime is newer than the `.py` mtime.

## Independent-derivation check (Mathematica)

GENUINELY INDEPENDENT. The `.wl` does not type the compiler matrices as literals (the way the `.py` does); it derives each from the physical first-order-drift premises and only afterward cross-checks against hard-typed SymPy reference matrices. Three corresponding sections:

1. **Observable→quotient compiler.** `.py:74` constructs it as a literal: `C_obs_to_quot = sp.diag(-Cstar, 1, 1)`. `.wl:103-115` instead derives it: `Solve[{sigTr == firstLogDrift[Exp[-cStar*small*rObs]], sigNt == firstLogDrift[Exp[small*nObs]], sigEta == firstLogDrift[Exp[small*etaObs]]}, quotVars]` then `obsToQuot = linearCompiler[quotFromObs, obsVars]` (= `Outer[D, quotFromObs, obsVars]`). The diagonal `−C_*` entry emerges from the first-log-drift of the finite ratio `Exp[−C_* small rObs]`, not from being typed. Different route (Series/Log + autodiff vs. literal matrix).

2. **Direct observable→defect compiler.** `.py:118-124` types `C_obs_to_def_expected = sp.Matrix([[1,0,0],[-Bstar,1,0],[Bstar,-1,-epsetas/(1-epsetas)]])`. `.wl:171-189` derives the `nObs` row from `nCompositeDrift = firstLogDrift[Exp[small*xi]*Exp[small*theta]^bStar]` (i.e. it actually expands `δln(T²·R_tr^{B_*}) = Ξ + B_* Θ` from the composite definition `N_*=T² R_tr^{B_*}`) and the third row from a genuine finite ratio (see #3), then `Solve[...]`/`Outer[D,...]`. The `−B_*` entry is produced by the series expansion of the composite, not copied.

3. **Complementary drift / third defect row.** `.py:167` types `dln_E = -epsetas/(1-epsetas)*de` directly. `.wl:173-176` derives it from the genuine finite observable ratio: `etaComplementDrift = SeriesCoefficient[Log[(1 - eta*Exp[small*etaObs])/(1 - eta)], {small, 0, 1}]`, i.e. it expands `δln(1−ε_η)` from the actual `(1−ε_η)` ratio. This is the clearest evidence the `.wl` works from the physical premise rather than echoing the `.py`'s pre-collapsed coefficient.

The presence of `sympyObsToQuot`, `sympyQuotToDef`, `sympyObsToDef`, `sympyDefToObs` (wl:116-120,146-150,192-196,227-231) is a *cross-engine consistency layer* applied AFTER the independent Solve/Series/D derivation, not the derivation itself — the correct dual-engine pattern. Threshold applied consistently: the `.wl` feeds DIFFERENT expressions (finite ratios → Series/Solve) to DIFFERENT operations (autodiff/linear-solve) than the `.py` (hand-typed matrices → `.inv()`/`.det()`). Not a port. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree at the level they claim. Shared final forms from the two `.txt` transcripts:
- `C_tr,*`, `C_*`, `B_*` identical; `A_tr,*` printed as `2χ/(1+χ+δ+χδ)` (wl `.txt:12`) = `2χ/((1+χ)(1+δ))` (py `.txt:22-24`) — same expression, different factoring.
- `det(C_obs->quot) = −(1+χ)(1+δ)(1+χ+δ)/(χδ) = −C_*` on both.
- `det(C_obs->def) = ε/(ε−1) = −ε/(1−ε)` on both (wl `.txt:51`, py `.txt:188-191`).
- All residual checks print `0` / `PASS` in both transcripts (every `expect_zero`/`expectZero` line). No sign, factor-of-2, or convention disagreement.

## Verdict justification

`findings`. The single finding is a low-severity `stale_output`: the committed SymPy `.txt` still carries the pre-renumber `STAGE 171` banner while the current `.py` prints `STAGE 188`; the body math is unchanged and all residuals are zero. Attacks tried and failed: (a) I checked whether the `.wl` is a transliteration of the `.py` — it is not; it derives every compiler via finite-ratio first-log-drifts + `Solve`/`Outer[D]` and only cross-checks the SymPy literals afterward. (b) I checked the zero-set iff for the `M·0=0` triviality trap — both scripts deliberately exercise bijection on the GENERIC packet (py:177-184, wl:253-260) plus nonzero-determinant confirmation, so the iff is non-tautological. (c) I checked the coefficient identities `C_*=1/C_{tr,*}` and `A_{tr,*}=B_* C_{tr,*}` — these are real algebraic facts about distinct hand-built expressions, not definitional. (d) Symbol domains (`χ0_*,δU_*>0`, `0<ε_η*<1`) match the notes' physical setup and justify the `1/(1−ε)` and determinant-nonzero claims. (e) Paper alignment: every notes deliverable maps to a faithful script check. I read the card, notes, and appendix row before the scripts; the script's verified claim matches the paper's stated claim.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation (symbolic stage; no numeric figure-of-merit constants are emitted). All values are MATCH against the notes `.md` (the natural carrier; the `.tex` card is terse by design and legitimately omits intermediate forms).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `C_{tr,*}=χδ/((1+χ)(1+δ)(1+χ+δ))` | py:42-45 / wl:72-75 / out:9 | md:172-176 | MATCH |
| `C_* = 1/C_{tr,*}` | py:46-49 / wl:76-79 / out:10 | md:234 (as `1/C_{tr,*}` factor) | MATCH |
| `B_* = 2(1+χ+δ)/δ` | py:50 / wl:80-83 / out:11 | md:121, 187-190 | MATCH |
| `A_{tr,*} = 2χ/((1+χ)(1+δ))` | py:51 / wl:84-87 / out:12 | md:179-184 | MATCH |
| identity `A_{tr,*}=B_* C_{tr,*}` | py:63 / wl:95 / out:15 | md:191-196, 394-396 | MATCH |
| `C_obs->quot = diag(−C_*,1,1)` | py:74 / wl:114-124 / out:22 | md:242-253 | MATCH |
| `det C_obs->quot = −C_* ≠0` | py:83 / wl:125 / out:23 | md:254-257 | MATCH |
| `C_quot->def` (triangular) | py:93-99 / wl:144 / out:36 | md:377-384 | MATCH |
| `C_obs->def` (lower-tri, B_*,−ε/(1−ε)) | py:118-124 / wl:189 / out:50 | md:301-312 | MATCH |
| `det C_obs->def = −ε_η*/(1−ε_η*)` | py:136 / wl:200 / out:51 | md:349-353 | MATCH |
| `C_def->obs` (inverse) | py:147-153 / wl:224 / out:68 | md:336-346 | MATCH |
| `Θ_1=δln R_tr` | py:139 / wl:204 / out:58 | md:287, 315-317 | MATCH |
| `Ξ_1=δln N_* − B_* δln R_tr` | py:140 / wl:205 / out:60 | md:320-323 | MATCH |
| `R_1=−ε/(1−ε)δln ε_η − Ξ_1` | py:141-144 / wl:206-209 / out:62 | md:325-331 | MATCH |
| `δln(1−ε_η)=−ε/(1−ε)δln ε_η = R_1+Ξ_1` | py:167-170 / wl:245-249 / out:83 | md:155-163 | MATCH |
| zero-defect iff theorem | py:177-196 / wl:253-271 | md:405-414 | MATCH |

INTERNAL (scaffolding, no prose expectation): `obs`/`obsVars`, `quot`/`quotVars`, `def`/`defVars` symbol triples; `small` series variable; `firstLogDrift`/`linearCompiler`/`clean`/`expectZero` helpers; the `sympy*` reference matrices (cross-engine check arrays); all `… − I` / `… − 0` residuals; `det·(1/det)−1` well-definedness probes.

Note: the notes define `Λ_0 := 27π²Gc_s⁵/(20a⁵c⁵)` (md:144) and the `T²` definition (md:113), but the scripts do not compute or emit these — they are upstream/contextual quantities, not Stage-188 script deliverables, so they are out of reconciliation scope (no script value to reconcile).

reconciliation: complete; 16 deliverable values checked, 0 misaligned.

## Self-test notes

Checked the three required traps. (1) Variable independence: every `Outer[D, outputs, inputs]` in the `.wl` differentiates expressions (`quotFromObs`, `defFromQuot`, `defFromObsDirect`, `obsFromDefSolved`) that genuinely depend on the differentiation variables (`obsVars`/`quotVars`/`defVars`) — confirmed via the `Solve` antecedents, so no identically-zero-derivative trap. (2) No unbounded integrals here (linear-algebra/series stage), so the parity trap is N/A. (3) Trivial-case/zero-vector trap: explicitly avoided — both scripts exercise the zero-defect iff on the GENERIC packet plus nonzero-determinant, so the iff cannot pass trivially via `M·0=0`. The only finding (stale_output) is a transcript-refresh item with no source edit, so the paper round-trip trap does not apply.
