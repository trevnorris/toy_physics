---
unit_id: 229
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 229 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_229.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows 70, 809-865)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The card's `\stagefield{Output}` states verbatim: "Selected-branch classifier theorem: the pure-transfer defect is denominator-like near softening, with crossover controlled by the exact threshold $\delta=8/9$." The derivation ledger adds: factorize the selected branch into numerator-like and denominator-like pieces, derive the classifier threshold $\delta=8/9$, and prove the near-softening denominator-like signature. The notes (the authoritative derivation source, §2-§5) enumerate five deliverables: (1) exact reduction of the selected-mode product $N_-(x)$ to the universal dimensionless $F(\xi,\delta)=(9\delta+11\xi)^4/[81(1-\xi)(9\delta^2+18\delta\xi+11\xi^2)^2]$ and its factorization $F=F_{\rm num}F_{\rm den}$; (2) the exact log-slope classifier $\mathcal R_{ND}=72\delta^2(1-\xi)/[(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)]$; (3) onset law $\mathcal R_{ND}(0,\delta)=8/(9\delta)$ and near-softening limit $\to0$; (4) the exact crossover cubic $\mathcal P=0$ with $\partial_\xi\mathcal P>0$ monotonicity; (5) the always-denominator threshold $\delta=8/9$ and sample crossover depths for $\delta=1/4,1/2,3/4$. The appendix (lines 809-865) restates the $N_-$/$s_-$ forms, the threshold $\delta=8/9$, and the qualitative co-loading / near-softening-denominator-like theorem, but does NOT write out the cubic $\mathcal P$ explicitly.

## What the script claims to verify

The SymPy script verifies, with substantive (non-tautological) `assert simplify(...)==0` checks: (a) the exact reduction $N_-(A\xi,A\delta) = (8\beta_0/\pi^2 A) F(\xi,\delta)$ built from the genuine constants $\kappa_0^2=8/\pi^2,\ \kappa_1^2=16/(9\pi^2)$ (line 49); (b) the factorization $F=F_{\rm num}F_{\rm den}$ (line 53); (c) the log-slopes $L_{\rm num},L_{\rm den}$ computed via `sp.diff(sp.log(...))` and the classifier $R_{ND}$ matched against the closed forms (lines 71-73); (d) onset $R_{ND}(0,\delta)=8/(9\delta)$ (line 76) and $\lim_{\xi\to1^-}R_{ND}=0$ (line 79); (e) the crossover cubic obtained by `numerator = expand(-together(R_ND-1).as_numer_denom()[0])` matched against `expected_P` (line 97), its derivative (line 101), and $\mathcal P(0,\delta)=9\delta^2(9\delta-8)$ (line 104); (f) numerical sample crossover roots for $\delta=1/4,1/2,3/4$ with left/right sign probes (lines 122-135); (g) an always-denominator slice for $\delta=1$ (lines 141-144). These are well-anchored: the cubic and classifier are re-derived in-script from the constants and `diff`, then compared to hand-written targets, so the targets are genuinely exercised.

## Paper ↔ script cross-check

| Paper/notes deliverable | Script check | Status |
|---|---|---|
| (1) $N_-\to F$ reduction + $F=F_{\rm num}F_{\rm den}$ | lines 49, 53 | match |
| (2) classifier $\mathcal R_{ND}$ + $L_{\rm num},L_{\rm den}$ | lines 71-73 | match |
| (3) onset $8/(9\delta)$ + near-softening $\to0$ | lines 76, 79, 83 | match |
| (4) crossover cubic $\mathcal P$ + $\partial_\xi\mathcal P>0$ | lines 97, 101, 104 | match (script's 121, see F1) |
| (5) threshold $\delta=8/9$ + sample depths | lines 104, 122-135, 141-144 | match |
| Mathematica independent re-derivation | (none) | missing (F2) |

The script's verified mathematics matches the **paper card** and **appendix** in full. The misalignment is a **notes-side internal typo** in the cubic leading coefficient (F1, `notes_contradicts_script`): the notes write `189\xi^3` but their own derivative `363\xi^2` corresponds to `121\xi^3`. The script's `121` is the arithmetically correct value, so the *physics* aligns; only the notes' written coefficient is wrong. `paper_alignment: partial` because the notes prose disagrees with the (correct) script.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 49 | `reduced == 0` ($N_-\to F$) | claim 1 | yes |
| A2 | sympy | 53 | `simplify(F - F_num*F_den)==0` | claim 1 | yes |
| A3 | sympy | 71 | `simplify(L_num - expected)==0` | claim 2 | yes |
| A4 | sympy | 72 | `simplify(L_den - expected)==0` | claim 2 | yes |
| A5 | sympy | 73 | `simplify(R_ND - expected)==0` | claim 2 | yes |
| A6 | sympy | 76 | onset $=8/(9\delta)$ | claim 3 | yes |
| A7 | sympy | 79 | `soft_limit == 0` | claim 3 | yes |
| A8 | sympy | 83 | $L_{\rm num}$ softening limit | claim 3 | yes |
| A9 | sympy | 97 | `simplify(numerator - expected_P)==0` | claim 4 | yes (re-derived from R_ND) |
| A10 | sympy | 101 | `simplify(dP - expected_dP)==0` | claim 4 | yes |
| A11 | sympy | 104 | `P(0,delta)=9delta^2(9delta-8)` | claim 4/5 | yes |
| A12 | sympy | 125 | `assert_close(root, expected)` | claim 5 | yes (root of in-script poly) |
| A13 | sympy | 131-134 | left>1, right<1 sign probes | claim 5 | yes |
| A14 | sympy | 143 | always-denominator $\delta=1$ | claim 5 | yes |

All sympy assertions are non-tautological and anchored. A9 is the key adversarial test: `numerator` is recomputed from `R_ND` (itself from `diff(log F_num/F_den)`) and compared to the literal `expected_P`, so a wrong literal would fail. That is exactly what catches F1 — the script's `121` is forced by the recomputed `numerator`, proving `121` correct.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** notes_contradicts_script
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.md:292`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.py:96`

**What's wrong:**
The notes write the crossover cubic with leading coefficient `189`:

> notes:292 — `\mathcal P(\xi,\delta) = 189\xi^3+297\delta\xi^2+333\delta^2\xi+81\delta^3-72\delta^2 = 0`

The script asserts leading coefficient `121`:

> script:96 — `expected_P = 121 * xi**3 + 297 * delta * xi**2 + 333 * delta**2 * xi + 81 * delta**3 - 72 * delta**2`

The script is arithmetically correct. Expanding the classifier-equals-one condition $72\delta^2(1-\xi) - (9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)=0$ gives the $\xi^3$ term entirely from $11\xi\cdot11\xi^2 = 121\xi^3$. The script's A9 check (line 97) re-derives this numerator from `R_ND` and would fail if `121` were wrong; the saved output (line 23) confirms `P = ... + 121*xi**3`. Decisively, the notes are **internally inconsistent**: notes:302 give $\partial_\xi\mathcal P = 363\xi^2+594\delta\xi+333\delta^2$, and $\partial_\xi(121\xi^3)=363\xi^2$ — but $\partial_\xi(189\xi^3)=567\xi^2$. So the notes' own derivative (which matches the script and is correct) contradicts the notes' own cubic. This is the identical pattern flagged for sibling stages 222/223 (notes-side coefficient typo against a correct script).

**Why this matters:**
The notes are the authoritative derivation source for the .tex card. A reader transcribing the cubic from the notes would propagate `189`, which is wrong and would break the monotonicity/onset arguments downstream (the $\delta=8/9$ threshold derivation in notes §5.1 relies on $\mathcal P(0,\delta)=9\delta^2(9\delta-8)$, which is unaffected, but the cubic's $\xi$-dependence is). The .tex card and appendix do not restate the cubic, so the published paper is not directly affected — but the notes should match the verified script. Per protocol, the auditor flags and does not resolve.

**Resolution:** routed to user via `## Resolve before fix_loop` in the directive. The math strongly indicates the notes' `189` is the typo and `121` is correct, but notes are prose owned by Claude/user, not Codex; the direction (fix notes, not script) requires user sign-off.

**Verification:**
After user chooses direction (a), notes:292 leading term reads `121\xi^3`. No script change. The script already verifies `121` (output line 23).

### F2 — missing_verification_script

**Severity:** medium
**Subtype:** missing_mathematica
**Files:**
- `mathematica/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_mathematica_audit.wl` (does not exist)
- card admission: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_229.tex:11` "Mathematica audit: none yet."

**What's wrong:**
No Mathematica script exists. This unit is `is_status_only_candidate: False` and `is_checkpoint: False`, so the dual-engine rule applies: a `.wl` is REQUIRED wherever Mathematica CAN independently verify, and every claim here (closed-form rational reduction, symbolic log-derivatives, polynomial limits, a cubic obtained by clearing denominators, derivative-positivity, and real-root location) is squarely within Mathematica's native capability. The test is "is it possible," not "is it necessary." It is clearly possible.

**Why this matters:**
A single-engine stage has no cross-engine guard against transliteration-independent algebra errors. Notably, the F1 typo lived in the notes precisely because there was no second symbolic engine forcing agreement; an independent `.wl` deriving the cubic from native `Together`/`Numerator`/`D` would have surfaced `121` vs `189` independently.

**Required change:**
Codex creates the `.wl` per the claim manifest in the directive. It must NOT transliterate the `.py`: use native Mathematica primitives (`Together`, `Numerator`, `Cancel`, `D[Log[...],xi]`, `Limit[...,Direction->-1]`, `CoefficientList`, `Reduce`/`Resolve` for root location and positivity) with a different decomposition order, and derive the cubic by `Numerator[Together[RND-1]]` rather than copying `expected_P`.

**Verification:**
`redteam exec-mathematica 229` runs the new `.wl` to exit 0 with all `expectZero[...]`/`If[...,Exit[1]]` checks passing; the engine cross-check then confirms the `.wl`'s independently-derived cubic equals $121\xi^3+297\delta\xi^2+333\delta^2\xi+81\delta^3-72\delta^2$, agreeing with the SymPy residuals.

## Independent-derivation check (Mathematica)

No `.wl` exists; transliteration cannot yet be assessed. The directive's claim manifest specifies an independent derivation route and an explicit anti-transliteration guard so the new `.wl` is not a line-by-line port.

## Engine cross-check

Only SymPy is present; no cross-engine comparison possible. This is itself finding F2.

## Verdict justification

`verdict: findings`. The SymPy script's verified mathematics matches the paper card and appendix in full, and all 14 assertions are non-tautological and well-anchored (A9 is a genuine re-derivation of the cubic from the classifier, not a self-comparison). Attacks tried that FAILED to break it: (i) checked whether `expected_P` is hardcoded/tautological — no, it is forced by the recomputed `numerator` at line 97; (ii) verified the cubic's leading coefficient by hand-expanding the classifier-equals-one condition — got `121`, matching the script; (iii) checked the onset/threshold arithmetic $\mathcal P(0,\delta)=9\delta^2(9\delta-8)\Rightarrow\delta=8/9$ — correct, matches paper; (iv) checked symbol domains (`positive=True` on $\xi,\delta,A,\beta_0$) — consistent with the notes' $0\le\xi<1,\ \delta>0$ stable-branch setup, and the `dir='-'` one-sided limit at $\xi=1$ is the correct softening edge. Two findings remain: F1 (notes-side typo `189` vs correct script `121`, routed to user) and F2 (missing required second engine). Neither is `stop_cold`: F1 is a prose typo against a correct script (no downstream math propagation since the published cubic appears only in notes, and even there the threshold derivation it feeds is unaffected), and F2 is a coverage gap, not an inconsistency. The math itself is sound.

## Self-test notes

Checked: (1) Variable independence — `sp.diff(sp.log(F_num),xi)` and `sp.diff(sp.log(F_den),xi)`: `F_num` depends on both $\xi,\delta$ and `F_den` on $\xi$, so neither derivative is identically zero; the F2 manifest's `D[...,xi]` targets likewise depend on $\xi$. (2) Trivial-case pre-check — hand-substituted $\xi=0$ into the cubic giving $81\delta^3-72\delta^2=9\delta^2(9\delta-8)$, confirming the script's line-104 literal and the $\delta=8/9$ threshold. (3) Coefficient cross-check — independently expanded $(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)=81\delta^3+261\delta^2\xi+297\delta\xi^2+121\xi^3$ and the full crossover numerator, confirming `121` (script correct) and exposing the notes' `189` as the typo via its own inconsistent `363`-derivative. (4) The F2 directive Target path includes `_mathematica_audit` and lives under `mathematica/`.
