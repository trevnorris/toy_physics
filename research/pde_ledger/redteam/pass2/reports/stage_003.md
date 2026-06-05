---
unit_id: 003
batch: I.1
auditor_model: claude-opus-4-8
audit_date: 2026-06-04T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: false
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage003_bdg_coupling.md]
  paper_appendix: present
---

# Audit unit 003 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_003.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage003_bdg_coupling.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row at line 9)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage003_bdg_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage003_bdg_mathematica_audit.txt`

## What the paper claims

Stage 003 turns on the first conservative matter-support dynamics. It reduces the linearized BdG sector to stable normal coordinates, couples them to the wall via the Stage-001 confinement variation, and integrates them out exactly. `\stagefield{Output}`: "Stage~003 outputs the scalar Schur complement \eqref{eq:app-stage003-scalar-kernel}, the grouped kernels \eqref{eq:app-stage003-p2-kernel}, the low-frequency coefficients \eqref{eq:app-stage003-d0d2d4}, the exact pole shift \eqref{eq:app-stage003-pole-shift}, and the first microscopic grouped-isotropy criterion." Distinct deliverables: (1) the exact scalar effective kernel `D0_eff = K0 - omega^2 M0 - C(Omega0^2 - omega^2 I)^{-1} C^T`; (2) its low-frequency `K_eff/M_eff/N_eff` renormalization; (3) the grouped-P2 channel kernels `D_A^eff = K_A - M_A omega^2 - sum g_{Aa}^2/(varpi_{Aa}^2 - omega^2)`; (4) the channelwise low-frequency `d0/d2/d4` coefficients; (5) the exact two-pole spectrum `omega_pm^2` and the weak-coupling wall/matter pole shifts `delta Omega_eta^2 = -+ g^2/[M(varpi^2-Omega_eta^2)]`; (6) the grouped trace/anomaly diagnostic `bar x, a_x, b_x` with the isotropy theorem `a_x=b_x=0`. The notes (Sec. 6) additionally show the harmonic selection rule (l=0 / l=2 block diagonality). This is a checkpoint stage, so the bar is higher: both engines required, assertions substantive, alignment exact.

## What the script claims to verify

The SymPy script (4 sections) verifies: (I) Euler–Lagrange equations of the documented two-wall + two-scalar-BdG Lagrangian, the frequency-space elimination yielding `D0_eff`, its match to the explicit 2x2 manual form, and the low-frequency `K_eff/M_eff/N_eff` series; (II) the single-mode dispersion `(K-Mw2)(varpi2-w2)-g^2`, the two closed-form roots, and the perturbative wall/matter shifts; (III) the grouped-P2 channel kernels, the `d0/d2/d4` coefficients, the `d2bar/a2/b2` diagnostics, and `a2=b2=0` on the isotropic branch; (IV) the spherical-harmonic orthogonality selection rule. The Mathematica script claims to mirror all four sections with an independent EL route (explicit conjugate-momentum time derivative), `LinearSolve` for the kernel, an independent root-sum/product check, and a 4x4 overlap-matrix identity check.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) Scalar Schur kernel `D0_eff` | sympy I.2 `Deff` + manual; mathematica I `dEff` + manual | match (sympy); mathematica section I corrupted (see F1) |
| (2) `K_eff/M_eff/N_eff` low-freq | sympy I.3 `series match`; mathematica `series match` | match (both use independent fresh matrices, unaffected by F1) |
| (3) Grouped-P2 kernels `D_A^eff` | sympy III D20/D21/D22; mathematica III | match |
| (4) `d0/d2/d4` coefficients | sympy III.1; mathematica III | match |
| (5) Two-pole spectrum + pole shifts | sympy II; mathematica II | match |
| (6) `bar x, a_x, b_x`, `a2=b2=0` | sympy III.2/III.3; mathematica III | match (formulas agree with `.tex` eq:app-stage003-trace-anomaly) |
| selection rule (notes Sec. 6) | sympy IV; mathematica IV | match |
| notes Sec. 6 indices `d_{237/238/239}` | scripts use `d_{2,20/21/22}` | mismatch (F3 — notes garble) |

`paper_alignment: partial` — the SymPy script aligns exactly with all six paper deliverables; the Mathematica Section-I EL block is internally corrupt (F1) so the second-engine coverage of deliverable (1)'s EL derivation is not actually exercised by the current script, and the notes carry a garbled index (F3).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 113-128 | `expect_zero` EL residuals (4) | claim 1 (Lagrangian/EL) | yes |
| A2 | sympy | 156 | `derived D0_eff vs Deff == 0` | claim 1 | yes |
| A3 | sympy | 174 | `Deff - manual == 0` | claim 1 | yes |
| A4 | sympy | 200 | `series match == 0` | claim 2 | yes |
| A5 | sympy | 237 | derived dispersion vs `(K-Mw2)(varpi2-w2)-g^2` | claim 5 | yes |
| A6 | sympy | 249-250 | root closed forms == 0 | claim 5 | yes |
| A7 | sympy | 264-265 | wall/matter perturbative shift | claim 5 | yes |
| A8 | sympy | 339-342 | isotropic a2/b2/D-diffs == 0 | claim 6 | yes |
| A9 | sympy | 374-383 | harmonic overlaps / norms | selection rule | yes |
| B1 | mathematica | 85-94 | EL residuals (4) | claim 1 | NO — corrupted Lagrangian (F1) |
| B2 | mathematica | 117 | `dDerived - dEff == 0` | claim 1 | NO — `dDerived` from corrupted EL (F1) |
| B3 | mathematica | 130 | `dEff - manual == 0` | claim 1 | yes (dEff is independent) |
| B4 | mathematica | 145 | `series match == 0` | claim 2 | yes |
| B5 | mathematica | 168 | derived dispersion == 0 | claim 5 | yes |
| B6 | mathematica | 175-178 | roots solve + sum/product | claim 5 | yes |
| B7 | mathematica | 187-188 | wall/matter shift | claim 5 | yes |
| B8 | mathematica | 227-230 | isotropic a2/b2/D-diffs | claim 6 | yes |
| B9 | mathematica | 259 | overlap - identity == 0 | selection rule | yes |

## Findings

(Ordered: paper_misalignment first, then script-side by severity.)

### F1 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:61-67`
- (consequent assertions) `:85-94`, `:113-117`
- saved output `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage003_bdg_mathematica_audit.txt:9-18`

**What's wrong:**
Line 54-59 already assigns the complete reduced Lagrangian (qa/qL inertia, K-potential, two BdG kinetic+potential terms, four couplings). Lines 61-67 then execute

```
lRed = lRed + (
  -1/2 kaa qa^2 - kaL qa qL - 1/2 kLL qL^2
  + 1/2 D[xa, t]^2 - 1/2 wa^2 xa^2
  + 1/2 D[xb, t]^2 - 1/2 wb^2 xb^2
  + c1a qa xa + c1b qa xb + c2a qL xa + c2b qL xb
);
```

which ADDS a *second copy* of the potential, both BdG kinetic terms, both BdG potential terms, and all four wall–matter couplings, while leaving the qa/qL inertia terms (`maa, maL, mLL`) single. The effective Mathematica Lagrangian therefore has DOUBLED `kaa,kaL,kLL`, DOUBLED BdG kinetic/potential terms, and DOUBLED `c1a,c1b,c2a,c2b`, but single inertia. The git commit that introduced this (`10d4027`, batch I.1) describes it as a patch to a prior "`lRed = ...` continuation defect that had captured only kinetic terms" — but the new block over-corrects and doubles instead of restoring.

Consequence, derived by hand (no execution): the qa Euler–Lagrange equation of the doubled Lagrangian is
`elQa = maa qa'' + maL qL'' + 2 kaa qa + 2 kaL qL - 2 c1a xa - 2 c1b xb`,
so the assertion at lines 85-87
`expectZero["qa equation", elQa - (maa qa'' + maL qL'' + kaa qa + kaL qL - c1a xa - c1b xb)]`
has residual `kaa qaFun[t] + kaL qLFun[t] - c1a xaFun[t] - c1b xbFun[t] != 0`. Likewise `xa equation` gets residual `xaFun''[t] + wa^2 xaFun[t] - c1a qaFun[t] - c2a qLFun[t] != 0` (X-kinetic and X-potential doubled). And `dDerived` (lines 113-116) is built from this corrupted `elQa/elQL`; for the [1,1] entry the factor-2 on the X-kinetic cancels in the X-solution but the doubled `kaa` and doubled coupling do not, giving `dDerived[1,1] - dEff[1,1] = kaa - c1a^2/(wa^2-omega^2) - c1b^2/(wb^2-omega^2) != 0`. So a faithful run of the CURRENT `.wl` must FAIL (Exit[1]) at "qa equation" and at "derived D0_eff vs Deff".

The committed saved output (`.txt:9-18`) nonetheless reports `qa equation = 0 / PASS` and `derived D0_eff vs Deff = {{0,0},{0,0}} / PASS`. This is mathematically impossible for the current script (confirmed against the committed output via `git show 10d4027:...txt`), which proves the committed transcript was captured from a transient correct (pre-doubling) version during the I.1 fix-loop, not from the file as committed. The SymPy script (`.py:90-105`) defines each Lagrangian term exactly once and is correct; its EL/`dEff` checks genuinely hold. So the second engine does NOT independently verify deliverable (1)'s EL derivation, and the two engines disagree on the constructed Lagrangian.

**Why this matters:**
Stage 003 is a checkpoint; deliverable (1) (the scalar Schur kernel and its EL provenance) must be verified by both engines. As written, the Mathematica EL block verifies a doubled, physically wrong Lagrangian — and would not even run clean. The orchestrator's independent re-run will exit nonzero. The captured "PASS" transcript masks the defect.

**Required change:**
Delete the spurious additive block at lines 61-67 (the `lRed = lRed + ( ... );` statement and its `Clear[...]` of unused locals is fine to keep or remove). `lRed` is already complete after line 59. The `Clear[qa0, ...]` at line 61 may remain. After removal, `elQa/elQL/elXa/elXb` will match the single-coefficient EL forms in the assertions at lines 85-94, and `dDerived` will match `dEff`. Re-run `math -script` and confirm all Section-I checks pass (Exit[0]). Do NOT touch the SymPy script — it is correct.

**Verification:**
After the edit, a fresh `math -script` run prints `qa equation = 0`, `xa equation = 0`, `derived D0_eff vs Deff = {{0,0},{0,0}}`, and exits 0 with the same downstream Section II–IV PASSes. The regenerated `.txt` should be byte-comparable to the (correct) committed transcript except possibly for FullSimplify form of intermediate prints.

### F2 — stale_output

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage003_bdg_mathematica_audit.txt` (full file)
- script `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:61-67`

**What's wrong:**
mtime ordering alone looks fine (output 2026-05-21 11:50 newer than script 2026-05-21 01:03), so this is not a pure timestamp-staleness flag. The substantive problem is content disagreement: the saved transcript's Section-I results (`qa equation = 0`, `derived D0_eff vs Deff = {{0,0},{0,0}}`) cannot be produced by the current `.wl` (see F1). The captured output reflects a different, correct version of the script than the one committed. This is the standard "output's content disagrees with what the current script would produce" sub-case and is therefore blocking, not merely informational.

**Why this matters:**
A reader trusting the committed transcript would conclude the Mathematica Section I passes, when in fact the committed script fails. The transcript must be regenerated from the corrected script.

**Required change:**
This finding is resolved automatically once F1 is applied and the script is re-run; the verifier should regenerate `mathematica/output/moving_throat_pde_stage003_bdg_mathematica_audit.txt` from the corrected `.wl`.

**Verification:**
Regenerated transcript Section I shows the same PASS lines AND is produced by the corrected (non-doubled) script; the orchestrator's independent re-run exits 0.

### F3 — paper_misalignment

**Subtype:** notes_contradicts_script (numbering garble)
**Severity:** low
**Files:**
- notes `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage003_bdg_coupling.md:367-372`
- paper card `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_003.tex:94-99` and `:131-140`
- script `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py:300-323`

**What's wrong:**
The notes (Sec. 6) write the grouped trace/anomaly diagnostics with garbled channel indices:

```
\bar d_2=\frac{d_{237}^{\rm eff}+2d_{238}^{\rm eff}+2d_{239}^{\rm eff}}{5},
a_2=\frac{2d_{237}^{\rm eff}-d_{238}^{\rm eff}-d_{239}^{\rm eff}}{10},
b_2=\frac{d_{238}^{\rm eff}-d_{239}^{\rm eff}}{2}.
```

The intended quantities are `d_{2,20}^{eff}, d_{2,21}^{eff}, d_{2,22}^{eff}` (the omega^2 coefficient `d_{2A}` for channels `A in {20,21,22}`). The literals `237/238/239` are evidently the `2` subscript run together with channel labels `20/21/22` minus a digit — a renumbering/typesetting garble. The `.tex` card uses the correct `d_{2A}` form (eq:app-stage003-d0d2d4, line 96-98; trace-anomaly block lines 131-140), and both scripts use `d220/d221/d222` (SymPy) / `d2coeffRaw[[1..3]]` (Mathematica) with the formulas `(d220 + 2 d221 + 2 d222)/5`, `(2 d220 - d221 - d222)/10`, `(d221 - d222)/2` — matching the `.tex`. So the script and `.tex` agree; only the notes string is garbled.

**Why this matters:**
Purely a documentation defect — the algebra is unaffected (the coefficients/combinations are correct everywhere else). But the garbled indices are confusing and are an instance of the known numbering-drift family. Direction of resolution is the user's call (this is a notes prose edit; Codex does not edit notes/).

**Required change (route to user, no Codex edit):**
See `## Resolve before fix_loop` in the directive. Likely fix: notes lines 367-372 `d_{237}/d_{238}/d_{239}` → `d_{2,20}/d_{2,21}/d_{2,22}` (or `d_2^{(20)}` etc.), matching the `.tex` and scripts.

**Verification:**
Notes Sec. 6 indices read `d_{2,20}/d_{2,21}/d_{2,22}` (or equivalent), consistent with `stage_003.tex` and the script `d220/d221/d222`.

## Independent-derivation check (Mathematica)

Apart from the Section-I corruption (F1), the Mathematica script is a genuinely independent re-derivation, not a transliteration of the SymPy choreography:
- Section I EL: SymPy uses `sympy.calculus.euler.euler_equations`; Mathematica builds a hand-rolled `timeD[]` total-time-derivative operator over algebraized coordinates (lines 71-83). Different machinery.
- Kernel: SymPy uses `Matrix.inv()` (line 155); Mathematica uses `LinearSolve` (line 100). Different method.
- Roots/dispersion (Section II): SymPy `sp.solve` the dispersion and compares to closed forms (lines 238-250); Mathematica instead checks the closed forms *solve* the dispersion and verifies sum/product of roots (lines 175-178). Different, complementary checks.
- Selection rule (Section IV): SymPy enumerates pairwise integrals; Mathematica builds a 4x4 `overlap` matrix and checks `overlap - IdentityMatrix[4] == 0` (lines 255-260). Different aggregation.
No `mathematica_transliteration` finding. (The structural problem is the doubling bug F1, not echoing SymPy.)

## Engine cross-check

SymPy (correct) and the CURRENT Mathematica (corrupted Section I) disagree on the constructed Lagrangian / EL derivation: SymPy's EL residuals are genuinely zero; the current Mathematica's would be nonzero (F1). For all other sections (II–IV) the two engines agree, and the committed transcripts agree where the Mathematica script is uncorrupted:
- `D0_eff` (sympy `.txt:21-32` ↔ mathematica `.txt:19`): identical kernel up to sign-of-denominator presentation (`c1a^2/(omega^2-wa^2)` = `-c1a^2/(wa^2-omega^2)`).
- Pole roots, perturbative shifts, P2 d2bar/a2/b2, selection rule: identical between transcripts.
The disagreement is confined to Section I and is the F1 defect; once F1 is fixed the engines will agree end-to-end.

## Verdict justification

`verdict: findings`. The SymPy script is correct and exercises all six paper deliverables plus the selection rule non-tautologically; I attacked the dispersion-derivation substitution (`subs(omega**4, w2**2).subs(omega**2, w2)`), the `solve`→closed-form root check, the perturbative series, and the isotropy substitution, and all hold. The Mathematica script, however, contains a high-severity Section-I Lagrangian-doubling defect (F1): lines 61-67 add a second copy of the potential/BdG/coupling terms, so the committed `.wl` would fail its own "qa equation" and "derived D0_eff vs Deff" assertions on a faithful run — and the committed transcript (which shows PASS) cannot have come from the committed script, making the saved output stale/inconsistent (F2). The notes carry a low-severity index garble `d_{237/238/239}` vs the correct `d_{2,20/21/22}` (F3, routes to user). No `stop_cold`: F1 is a clean script fix (delete the spurious additive block), F2 follows from re-running, and F3 is a one-line notes correction that does not propagate (the algebra is correct everywhere it matters). Downstream stages (022/023 collect the B_{A0/A2/A4} moments) consume the kernel-moment FORMS, which SymPy already verifies correctly and the `.tex` states correctly, so F1's fix does not change any downstream-quoted value — not CRITICAL_DOWNSTREAM.

## Value Reconciliation (pass-2 augmentation)

Stage 003 is an entirely SYMBOLIC stage: the scripts emit no pinned numeric constants, benchmarks, or figures-of-merit. The deliverable "values" are closed-form symbolic results. Reconciliation is therefore form-level. The committed `.txt` transcripts are the authoritative record of emitted results; I base reconciliation on script source + transcripts and reason only (no execution). Note: the Mathematica Section-I transcript is inconsistent with its current script (F2), so for the Section-I deliverables I rely on the SymPy transcript and the script source.

| value (symbolic deliverable) | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `D0_eff = K0 - w^2 M0 - C(Omega0^2 - w^2 I)^{-1} C^T` | py:155 / py.txt:21-32; wl:100 | tex:41-44 (eq scalar-kernel); md:182-188 | MATCH |
| `D0_eff` explicit 2x2 manual form | py:162-173 / py.txt:21-32; wl:118-127 | (intermediate, implied by tex eq) — in md:182-193 | MATCH |
| `K_eff = K0 - C Omega0^{-2} C^T` | py:181-186 / py.txt:82; wl:133-136 | tex:59-60; md:208-213 | MATCH |
| `M_eff = M0 + C Omega0^{-4} C^T` | py:187-192 / py.txt:83; wl:137-140 | tex:61-62; md:215-220 | MATCH |
| `N_eff = C Omega0^{-6} C^T` | py:193-198 / py.txt:84; wl:141-144 | tex:63-64; md:222-227 | MATCH |
| `D_A^eff = K_A - M_A w^2 - sum g_{Aa}^2/(varpi^2 - w^2)` | py:285-287 / py.txt:135-137; wl:197-199 | tex:82-87; md:253-262 | MATCH |
| `d0A = K_A - sum g^2/varpi^2` | py:297-299 / py.txt:139-141 | tex:96; md:275-278 | MATCH |
| `d2A = -(M_A + sum g^2/varpi^4)` | py:300-302 / py.txt:143-145 | tex:97; md:279-285 | MATCH |
| `d4A = -sum g^2/varpi^6` | py:303-305 / py.txt:147-149 | tex:98; md:286-290 | MATCH |
| exact roots `omega_pm^2 = [Om^2+varpi^2 ± sqrt((Om^2-varpi^2)^2 + 4g^2/M)]/2` | py:247-248 / py.txt:94-105; wl:170-171 | tex:117-119; md:322-328 | MATCH |
| wall shift `delta Om^2 = -g^2/[M(varpi^2-Om^2)]` | py:264 / py.txt:112-115; wl:187 | tex:124-125; md:333-338 | MATCH |
| matter shift `delta varpi^2 = +g^2/[M(varpi^2-Om^2)]` | py:265 / py.txt:117-121; wl:188 | tex:126-127; md:340-345 | MATCH |
| `bar d2 = (d2,20 + 2 d2,21 + 2 d2,22)/5` | py:321 / py.txt:154-159; wl:213 | tex:133-134 (general `bar x`); md:367 | MATCH (.tex); MISMATCH index in md (F3) |
| `a2 = (2 d2,20 - d2,21 - d2,22)/10` | py:322 / py.txt:160-165; wl:214 | tex:135-138; md:370 | MATCH (.tex); MISMATCH index in md (F3) |
| `b2 = (d2,21 - d2,22)/2` | py:323 / py.txt:166-171; wl:216 | tex:139; md:372 | MATCH (.tex); MISMATCH index in md (F3) |
| isotropy theorem `a2 = b2 = 0` | py:339-340 / py.txt:176-177; wl:227-228 | tex:153 (boxed); md:386 | MATCH |
| selection rule (overlaps = delta_{l l'}) | py:374-383 / py.txt:186-195; wl:259 | md:68-74, 385-386 | MATCH (notes carrier; tex states block structure narratively) |

INTERNAL scaffolding (accounted for, no finding): pass/fail flags (`PASS:`/`FAIL:`), `expect_zero`/`expectZero` residual prints (all `= 0`), the Euler–Lagrange residual checks (qa/qL/xa/xb equations) which verify provenance not a reported value, `dDerived - dEff`, `series match`, `derived dispersion`, `sum/product of roots`, the `T0/Ta/Tb` projector matrices defined but unused in assertions (wl:204-206), phase/ansatz symbols, harmonic norm checks (`Y00 norm - 1` etc.).

reconciliation: 17 deliverable forms checked; 0 numeric mismatches; 1 documentation index-garble in the notes (`d_{237/238/239}` → `d_{2,20/21/22}`, folded as F3, a low-severity notes-side paper_misalignment). All symbolic deliverable forms reconcile between scripts and the `.tex` card; the only doc discrepancy is the F3 notes index garble.

## Self-test notes

Variable independence: F1's prescribed fix (deleting lines 61-67) leaves `lRed` depending on qa/qL/xa/xb exactly once each; the EL `timeD[D[lAlg, v..]]` derivatives are over genuine dependencies, so no identically-zero-derivative trap. Parity/symmetry: Section IV integrals are over the symmetric sphere; the cross-integrals (e.g. Y00·Y20, Y20·Y21c) vanish by genuine orthogonality (odd-m phi-integral or l-orthogonality), not by a domain artifact — confirmed the integrand parities. Trivial-case: for the corrected (single) Lagrangian, the qa EL coefficient of `kaa` is exactly 1, matching the assertion's subtracted `kaa qaFun[t]`, so the residual reduces to 0 — confirming the fix is correct and not merely accidentally-true. Path: F1/F2 edits target the existing `.wl` in `mathematica/`; no new-script path needed. F3 is routed to the user (notes prose), not Codex.

---

## Orchestrator review (exec-confirmed, 2026-06-04)

The audit agent's F1/F2 were resolved by direct execution + instrumentation (orchestrator, single Mathematica seat):

- **F1 `insufficient_verification` (HIGH) — FALSE POSITIVE.** The current `.wl` runs (exit 0) and passes ALL 19 checks; `expectZero` is sound (`Exit[1]` on any nonzero residual). Instrumenting the real file shows `lRed`'s `qa^2` coefficient is `0` BEFORE the `lRed = lRed + (...)` block and `-1/2 kaa` AFTER it — i.e. the block is NOT a doubling. Root cause: in the original multi-line assignment (lines 54-59), line 55 (`...+ 1/2 mLL D[qL,t]^2`) is a syntactically complete expression and lines 56-59 each begin with `-`, so the WL script parser evaluates lines 56-59 as SEPARATE, DISCARDED statements. The original `lRed` therefore captures only the kinetic terms; the parenthesized re-add (lines 62-67) supplies the stiffness/BdG/coupling exactly once. Net `lRed` is single-coefficient and correct; the derived EL equations match the hand-typed single-coefficient targets and the independent `dEff` matrix elimination.
- **F2 `stale_output` (MED) — FALSE POSITIVE.** A fresh run is byte-identical to the committed `.txt` (19 PASS / 0 FAIL, `diff` empty). The committed transcript IS reproducible by the current script.
- **F3 `paper_misalignment` (LOW) — VALID.** Notes index garble `d_{237}/d_{238}/d_{239}` vs the correct `d_{2,20/21/22}` in `.tex` + scripts. Routes to user (notes-side fix, canonical lives in card+scripts).

**New LOW-severity robustness observation (not a correctness finding):** the original multi-line `lRed` (lines 55-59) silently drops lines 56-59 due to newline-before-`-` statement splitting; correctness is rescued only by the parenthesized re-add. This is fragile — a future "remove the redundant duplicate" cleanup (precisely the agent's misread) would break the stage. Optional clarity fix: consolidate into a single parenthesized assignment. Defer to user.

**Net:** stage 003 has NO confirmed script-side defect. Only F3 (notes garble) is genuine, pending user resolution.
