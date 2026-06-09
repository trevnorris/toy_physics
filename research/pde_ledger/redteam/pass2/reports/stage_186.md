---
unit_id: 186
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage186_similarity_orbit_closure.md]
  paper_appendix: present
---

# Audit unit 186 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_186.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage186_similarity_orbit_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows/sections referencing stage 186: line 103 status row; §"Similarity-orbit theorem" lines 1049–1085, 1210 claim-status)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` (stage_186.tex:15): "Shows the zero-defect equations are the tangent equations of a five-parameter monomial-preserving similarity orbit \(\mathcal G_*\)." The notes (§§2–6) spell out the deliverables: (D1) the exact rank-3 monomial-drift matrix `M_*` (3×8), whose rows are the log-exponent vectors of `(C_tr, C_nt, eps_eta)`; (D2) the convenient `(τ_1,κ_η,μ_1)` minor determinant `= 1+χ_{0,*} > 0`, hence rank 3; (D3) `dim ker M_* = 8−3 = 5` (boxed, notes:196–199); (D4) the finite five-parameter orbit `G_*` with the three dependent scaling laws `K_η→e^{2C−U}`, `T_U→e^{U−((1+δ)/(1+χ))(Γ+C−U)}`, `μ_W→e^{M(Λ,C,Γ,U,W)}` and the explicit `M(·)` exponent (notes:229–260); (D5) exact preservation of the three direct monomials under `G_*`; (D6) the tangent-space relations `κ_η=2c_1−κ_U`, the `τ_1` formula, and the `μ_1` formula matching the Stage 185 compatibility ledger; (D7) `ker M_* = T_id G_*`. The appendix (lines 1052–1085, 1210) confirms the same finite-family content and the exact-closure claim status.

## What the script claims to verify

The SymPy/Mathematica scripts (a) write down `M_*` as a literal 3×8 matrix and the three `dlog` monomial drifts as literal coefficient combinations; (b) check the `(τ_1,κ_η,μ_1)` minor determinant equals `1+χ`; (c) `solve` the three drift equations `=0` for `(τ_1,κ_η,μ_1)` and assert the solutions equal the paper's hardcoded target formulas; (d) build the finite-orbit dependent scalings (sympy posits them by hand; the .wl additionally re-derives them by `Solve[mMat.orbitLogVector==0]`), assert they preserve all three monomials, and assert they match the paper's stated scaling laws; (e) linearize the orbit at the identity and assert the result reproduces the solved compatibility formulas. The docstring (py:8–15) states these five checks as the script's purpose.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 `M_*` rows = log-exponent vectors | Matrix hardcoded (py:52–56); .wl compares hardcoded `mMat` against `mMatDerived` built from hardcoded exponent lists | partial (rows asserted, not derived from monomial defs — see F1) |
| D2 minor det `= 1+χ`, rank 3 | py:61–64 / wl:73–77 `det == 1+chi` | match |
| D3 `dim ker M_* = 5` | not asserted (only rank-3 via nonzero minor + 8 cols implied) | missing-from-script (see Value Reconciliation; informational) |
| D4 finite-orbit scaling laws | py:101–103 (posited) + py:122–124; wl:108–130 (solved from kernel) → matched to paper forms (py:128–130 via expect_zero, wl:114–126) | match |
| D5 monomial preservation | py:110–116 / wl:132–138 `preserves C_tr/C_nt/eps_eta == 0` | match |
| D6 tangent relations = Stage 185 ledger | py:81–96 / wl:89–101 (each formula vs hardcoded target) | match |
| D7 `ker M_* = T_id G_*` | py:127–130 / wl:150–154 linearization == solved compatibility formulas | match (this is the operational content of D7) |

`paper_alignment: partial` — all load-bearing identities are present and the values agree, but D1's row-content is asserted rather than independently derived (F1, internal-strength concern), and the `dim ker = 5` boxed deliverable is never explicitly emitted (Value-Reconciliation note, low-severity).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 63 | `det_minor != 1+chi → raise` | D2 | yes |
| A2 | sympy | 81–84 | `expect_zero(tau_formula − paper τ_1)` | D6 | yes |
| A3 | sympy | 85 | `expect_zero(kEta − (2c1−kU))` | D6 | yes |
| A4 | sympy | 86–96 | `expect_zero(mu_formula − Stage185 form)` | D6 | yes |
| A5 | sympy | 114–116 | `expect_zero(preserves C_tr/C_nt/eta)` | D5 | yes |
| A6 | sympy | 121–124 | solve K_η scaling, match 2C−U | D4 (eps_eta leg) | partial (round-trip, see F2) |
| A7 | sympy | 128–130 | `expect_zero(linearized − solved formula)` | D7 | yes |
| A8 | mathematica | 69–71 | `expectZero(mMatDerived − mMat == 0)` | D1 | no (both hardcoded — F1) |
| A9 | mathematica | 77 | `detMinor === 1+chi` | D2 | yes |
| A10 | mathematica | 89–101 | `expectZero(solved − paper τ/κ/μ)` | D6 | yes |
| A11 | mathematica | 114–126 | `expectZero(kernel-solved scaling − paper law)` | D4 | yes (kernel-derived) |
| A12 | mathematica | 136–138 | `expectZero(preserves monomials)` | D5 | yes |
| A13 | mathematica | 146–148 | round-trip eta scaling | D4 (eps_eta leg) | partial (F2) |
| A14 | mathematica | 152–154 | `expectZero(linearized − solved)` | D7 | yes |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl:35-62`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py:44-56`

**What's wrong:**
The two engines share the same load-bearing object — the literal `M_*` matrix and the literal monomial log-drifts — rather than each deriving them independently from the physical monomial definitions in the notes (`C_tr=(γc/K_U)^{1+δ}(π²T_U/L²K_U)^{1+χ}`, etc.). In SymPy the drifts are written directly:
```
dlog_Ctr = (1+delta)*(gam1+c1-kU) + (1+chi)*(tau1-kU)   # py:44
M = sp.Matrix([[0, 1+delta, 1+delta, -(2+chi+delta), 0,0,0, 1+chi], ...])  # py:52
```
The `.wl` re-packages the *same numbers* as hardcoded exponent lists and dots them with the drift vector:
```
ctrCoreExponents = {0,1,1,-1,0,0,0,0}; thermalExponents = {0,0,0,-1,0,0,0,1};  # wl:35-36
monomialLogDrifts = {(1+delta)*logDriftFromExponents[ctrCoreExponents] + (1+chi)*logDriftFromExponents[thermalExponents], ...}  # wl:41-48
mMat = {{0,1+delta,1+delta,-(2+chi+delta),0,0,0,1+chi}, ...}  # wl:58-62
```
Then the `.wl`'s "M_* row N matches paper" checks (wl:69–71) compare `mMatDerived` (Coefficient-extracted from `monomialLogDrifts`) against `mMat` — but BOTH are built from the same hardcoded exponent constants, so the check is also tautological (it confirms `Coefficient` inverts the `.` product, not that the row content is physically correct). Every subsequent step (minor det, compatibility solve, monomial-preservation, linearization) is a line-for-line port of the SymPy script operating on this shared hardcoded matrix. The one genuinely independent move is the `.wl`'s `Solve[mMat.orbitLogVector==0]` for the scalings (wl:108–109) vs SymPy positing them (py:101–103) — useful, but it still consumes the same hardcoded `mMat`.

**Why this matters:**
The second-engine policy requires Mathematica to provide an *independent* derivation route. Here neither engine derives `M_*` from the monomial definitions; both transcribe the same answer matrix, so a transcription error in the exponent vectors would pass both engines identically. The "row matches paper" assertions give false confidence — they cannot fail.

**Required change:**
In the `.wl`, derive at least one `M_*` row from the actual symbolic monomial (e.g. construct `Log[C_tr]` from `C_tr=(γ c/K_U)^{1+δ}(π²T_U/L²K_U)^{1+χ}` in terms of the eight log-couplings and read off the exponent vector via `Coefficient`/`D[...,logvar]`), then compare that derived row to `mMat[[1]]`. This makes A8 non-tautological and gives the second engine a real independent check of D1. (If a full independent re-derivation is judged out of scope for this stage, downgrade — but at minimum the row-match check must compare a *derived* row to the hardcoded one, not two hardcoded representations.)

**Verification:**
After the fix, the `.wl` should construct row 1 (and ideally rows 2,3) from `Log` of the monomial expression and `expectZero[derivedRow − mMat[[N]]]` should still be 0, but now non-trivially. The output line "M_* row 1 matches paper = 0" must come from a `Log`-derived row, not from `Coefficient` of the same hardcoded `monomialLogDrifts`.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py:117-124`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl:139-148`

**What's wrong:**
The block labeled "Non-tautological ground check" is in fact a round-trip. It defines `eps_eta_logdrift = 2*C - U - eta_scaling` (py:120), solves `eps_eta_logdrift == 0` for `eta_scaling`, then asserts `solved_eta - (2*C - U) == 0` (py:122). The solve trivially returns `eta_scaling = 2C − U` by construction, so the assertion `solved − (2C−U) == 0` cannot fail regardless of the physics — it confirms the solver inverts a linear equation. The companion check `eps_eta_logdrift.subs(eta_scaling, Eta_exp) == 0` (py:124) is non-trivial only because `Eta_exp = 2C−U` was *defined* at py:101 to be exactly that, so it too is guaranteed. The block's comment claims it is "non-tautological," which is incorrect.

**Why this matters:**
A comment asserting non-tautology on a tautological check is a precise inversion of the auditor's signal and could mask a future regression. The genuine eps_eta-preservation content is already covered non-trivially by A5/A12 (`finite orbit preserves eps_eta`), so this block adds no coverage.

**Required change:**
Either (a) delete the round-trip block (py:117–124 / wl:139–148) since `finite orbit preserves eps_eta` already exercises the same claim, or (b) make it genuine by starting from the *physical* definition `eps_eta = c_etaU^2/(K_U K_eta)` and its actual log-drift `2*c1 - kU - kEta` (note: it is `2c_1`, not `2C`), solving for the `K_η` scaling that zeroes it, and confirming it matches `2C − U`. As written, it is scaffolding mislabeled as a ground check.

**Verification:**
After the fix, either the block is removed (output no longer prints "K_eta preserving scaling matches paper 2C-U") or the solved scaling is derived from the `2c1 − kU − kEta` drift rather than from the already-target-shaped `2C − U − eta_scaling`.

## Independent-derivation check (Mathematica)

The `.wl` is largely a port. The shared, load-bearing `mMat` matrix and monomial-drift definitions are hardcoded identically in both engines (wl:58–62 ↔ py:52–56; wl:41–48 ↔ py:44–46), and the compatibility solve, minor determinant, monomial-preservation, and linearization checks are step-for-step the same operations in Mathematica syntax. The single divergence is constructive: the `.wl` derives the three dependent orbit scalings by `Solve[mMat . orbitLogVector == 0, {kEtaScale, muScale, tauScale}]` (wl:108–112), whereas the `.py` posits them by hand (py:101–103). That is a genuine extra check on D4, but it operates on the same hardcoded `mMat`, so it does not lift the script out of "port" status for the core object `M_*`. Strongest single evidence: `wl:58-62` `mMat` is character-for-character the same literal matrix as `py:52-56`, and the `.wl`'s "row matches paper" check (wl:69–71) compares that hardcoded matrix to another hardcoded representation rather than to a `Log`-derived row. Call: **PORT** (with one constructive divergence on the orbit-scaling leg).

## Engine cross-check

Both engines pass every assertion (sympy `.txt` lines 24–26, 34–38, 43–45 all `= 0`; mathematica `.txt` all `PASS`). The final symbolic forms agree up to representation: e.g. SymPy `tau1 = (-c1*delta - c1 + chi*kU - delta*gam1 + delta*kU - gam1 + 2*kU)/(chi+1)` (sympy.txt:21) and Mathematica `tau1 = (-((1+delta)*(c1+gam1)) + (2+chi+delta)*kU)/(1+chi)` (math.txt:22) are the same rational function in factored vs expanded form. `det = 1+chi` in both. No disagreement. `engines_agree: true`.

## Verdict justification

The math is sound and the load-bearing identities all reconcile with the paper: the minor determinant `1+χ`, the three compatibility formulas, the monomial-preservation of the orbit, and the linearization→Stage-185 equivalence are all verified non-tautologically and agree across engines. However the verdict is `findings`, not `clean`, for two reasons: (F1) the second engine is essentially a port of the first for the central object `M_*` — both hardcode the matrix and its `.wl` "row matches paper" check compares two hardcoded representations, so it cannot detect a transcription error in the exponent vectors; and (F2) the block self-labeled "Non-tautological ground check" is in fact a round-trip that cannot fail. Both are script-internal-strength issues, not paper↔script value disagreements, so no user-resolution gate is triggered. I attacked the minor-column selection (correctly maps to (τ,κ_η,μ)), the sign of row 2 `2(1+E)=2+2E` (correct), the assumption set (`chi,delta>0` matches the paper's constructive coherent branch `χ_{0,*}>0`, and is needed for `det=1+χ>0` and the orbit positivity), and the linearization substitution basis (correct: `Λ,C,Γ,U,W → λ1,c1,γ1,κ_U,κ_W`); these held. I confirmed I read the stage card, the notes, and the appendix section; the script's verified claims match the paper's stated Output for D2/D4/D5/D6/D7.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 9 deliverable values checked, 0 misaligned (see D3 note below — a stated deliverable not emitted by the script, flagged informationally, not as a value MISMATCH).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `M_*` row 1 `[0,1+δ,1+δ,−(2+χ+δ),0,0,0,1+χ]` | py:53 / wl:59 / sympy.txt:10, math.txt:9 | notes:174 (`M_*` boxed) | MATCH |
| `M_*` row 2 `[2(1+E),0,2E,F−E,−1,−(2+E),1,−F]` | py:54 / wl:60 / sympy.txt:11 | notes:175 | MATCH |
| `M_*` row 3 `[0,2,0,−1,−1,0,0,0]` | py:55 / wl:61 / sympy.txt:14 | notes:176 | MATCH |
| minor det `(τ,κ_η,μ) = 1+χ` | py:62 / wl:76 / sympy.txt:16, math.txt:17 | notes:191 (`det = 1+χ_{0,*}>0`); appendix:1165 | MATCH |
| `τ_1 = κ_U − ((1+δ)/(1+χ))(γ_1+c_1−κ_U)` | py:82-84 / wl:91 / sympy.txt:24 PASS | notes:318-325 (boxed); appendix tangent | MATCH |
| `κ_η = 2c_1 − κ_U` | py:85 / wl:93 / sympy.txt:25 PASS | notes:313-316 (boxed) | MATCH |
| `μ_1 = 2c_1−κ_U+2κ_W−2λ_1 − E(2γ_1+2λ_1−κ_U−κ_W) − F((1+δ)/(1+χ))(γ_1+c_1−κ_U)` | py:86-96 / wl:94-101 / sympy.txt:26 PASS | notes:327-335 (boxed) | MATCH |
| orbit scaling laws: `K_η→2C−U`, `T_U→U−((1+δ)/(1+χ))(Γ+C−U)`, `μ_W→M(Λ,C,Γ,U,W)` with the explicit `M` exponent | py:101-107 / wl:114-130 / sympy.txt:31-33, math.txt:41-43 | notes:229-260 (boxed); appendix:1062-1084 | MATCH |
| monomial preservation `δ_{G*} ln C_tr=C_nt=eps_eta=0` | py:114-116 / wl:136-138 / sympy.txt:34-36 | notes §4 (271-299) | MATCH |
| `dim ker M_* = 5` | NOT emitted by either script (only rank-3 via nonzero minor + 8 columns) | notes:196-199 (boxed `dim ker M_* = 8−3 = 5`); §7.1 codim-3 | MISSING-from-script (informational, see note) |

**D3 / `dim ker M_* = 5` note:** This is a boxed deliverable in the notes (line 198) and the conceptual headline of §7.1 / the appendix codimension-3 statement. The scripts establish rank 3 (nonzero `(τ,κ_η,μ)` minor) on an 8-column matrix, which *implies* `dim ker = 5`, but neither script emits or asserts the kernel dimension explicitly (e.g. no `M.nullspace()` / `Length[NullSpace[mMat]] == 5`). Under the augmentation's "Guards," this is not a MISMATCH (no conflicting value exists) and the value lives correctly in the notes; rank-3 is the load-bearing fact and it *is* verified. I therefore record it as an informational MISSING-from-script gap rather than a `paper_misalignment` finding — the deliverable is supported by the verified rank-3 minor. (If the orchestrator wants strict deliverable-for-deliverable script coverage, an `assert M.nullspace()` length-5 check would close it; not raised as a blocking finding here.)

INTERNAL scaffolding (no finding): `dlog_Ctr/Cnt/eta` expanded prints, `subs_basis`/`basisSubs`, `eta_scaling` round-trip intermediate (F2 covers its mislabeling), all `expect_zero/expectZero = 0` residuals, PASS flags, banner strings, conclusion prose.

## Self-test notes

Checked: (1) Variable independence — no `diff`/`D` derivatives in either script, so the "zero-derivative" trap does not apply; the checks are linear-algebra solves and substitutions. (2) Symmetry/parity — no integrals; n/a. (3) Trivial-case — minor-column map `[7,4,6]`/`{8,5,7}` correctly selects `(τ_1,κ_η,μ_1)` from the drift order `(λ1,c1,γ1,κU,κη,κW,μ1,τ1)`; det `1+χ` is nonzero on `χ>0`, confirming rank 3. (4) Path specs — no missing-script finding; both engines present. (5) Paper round-trip — F1/F2 are script-internal-strength fixes (derive a row / drop a round-trip); neither introduces a new constant, so no new `paper_misalignment` is created. The `chi,delta>0` assumptions match the paper's constructive coherent branch and are required for `det=1+χ>0`.
