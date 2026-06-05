---
unit_id: 074
batch: III.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage074_family1_healing_lock.md]
  paper_appendix: present
---

# Audit unit 074 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_074.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage074_family1_healing_lock.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 122/126/339 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.txt`

## What the paper claims

Stage 074 fixes the Family-1 support scale by a healing-length lock. The card states the healing closure `\ell = \hbar/(2 m c_{s,w})` (eq:app-stage074-healing), which forces the dimensionless support scale `\chi_s = m c_{s,w} L/\hbar = L/(2\ell) = \Lambda_\ell/2`. With the Stage-073 input `L/\ell = \Lambda_\ell = 37`, this gives `\chi_s = 37/2` (eq:app-stage074-chi). Carrying forward the Stage-071 branch coefficient `\kappa = 4\chi_s^2 + (4/5)\Lambda_\ell^2` and inserting `\chi_s = \Lambda_\ell/2` reduces it to `\kappa = (9/5)\Lambda_\ell^2 = 12321/5`, with derived scale `\alpha = \sqrt\kappa = 111/\sqrt5` (eq:app-stage074-kappa). The verbatim `\stagefield{Output}{Family--1 support values \(\chi_s=37/2\), \(\kappa=12321/5\), \(\eta=37\).}` lists three deliverables; `\eta=37` is the Stage-073 carry-forward `\eta=\Lambda_\ell`, not a fresh derivation here. The notes add the numeric `\alpha \approx 49.6407091`. Appendix row 126 summarizes: `\(\ell=\hbar/(2mc_{s,w})\), \(\chi_s=37/2\), \(\kappa=12321/5\)`.

## What the script claims to verify

Both scripts perform the same short symbolic chain from the physical definition `\chi = m c_s L/\hbar`: substitute the healing width `c_s -> \hbar/(2 m \ell)` to obtain `\chi = L/(2\ell)`, substitute `L -> \Lambda_\ell \ell` to obtain `\chi_s = \Lambda_\ell/2`, then form `\kappa = 4\chi_s^2 + (4/5)\Lambda_\ell^2` and assert it equals `(9/5)\Lambda_\ell^2`. They then substitute `\Lambda_\ell = 37` and assert the reference-branch values `\chi_s = 37/2`, `\kappa = 12321/5`, `\alpha = 111/\sqrt5` (and print `\alpha` numerically `49.6407091...`). The SymPy script comment (lines 44–47) is explicit that the coefficients `4` and `4/5` are carried forward from the earlier Family-1 Euler–Lagrange branch and that this stage only verifies the reduction under the lock.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `\ell = \hbar/(2 m c_{s,w})` (healing closure) | applied as substitution `c_s -> hbar/(2 m ell)` (py:37, wl:33–36) | match |
| `\chi_s = L/(2\ell) = \Lambda_\ell/2` | `expect_zero("chi_s - Lambda_ell/2", ...)` (py:53, wl:42) | match |
| `\kappa = (9/5)\Lambda_\ell^2` | `expect_zero("kappa - (9/5) Lambda_ell^2", ...)` (py:54, wl:43) | match |
| `\chi_s = 37/2` (reference) | `expect_zero("chi_ref - 37/2", ...)` (py:68, wl:57) | match |
| `\kappa = 12321/5` (reference) | `expect_zero("kappa_ref - 12321/5", ...)` (py:69, wl:58) | match |
| `\alpha = 111/\sqrt5 \approx 49.6407091` | `expect_zero("alpha_ref - 111/sqrt(5)", ...)` + numeric print (py:70/66, wl:59/55) | match |
| `\eta = 37` | input identity `\eta=\Lambda_\ell`; carried as `Lambda_ref=37` (py:56, wl:45), not separately asserted | match (carry-forward; see note) |

Note on `\eta=37`: the notes (md:9–13, 100–103) define `\eta = \Lambda_\ell` as a Stage-073 carry-forward, not a new Stage-074 derivation. The scripts use `\Lambda_\ell = 37` as the substitution point, which is the same quantity. A separate `assert eta == 37` would be a tautological restatement of the input, so its absence is not a `script_missing_paper_claim`. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `expect_zero(chi_lock - Lambda_ell/2)` | `\chi_s=\Lambda_\ell/2` | yes |
| A2 | sympy | 54 | `expect_zero(kappa_lock - (9/5)Lambda_ell^2)` | `\kappa=(9/5)\Lambda_\ell^2` | yes |
| A3 | sympy | 68 | `expect_zero(chi_ref - 37/2)` | `\chi_s=37/2` | yes |
| A4 | sympy | 69 | `expect_zero(kappa_ref - 12321/5)` | `\kappa=12321/5` | yes |
| A5 | sympy | 70 | `expect_zero(alpha_ref - 111/sqrt(5))` | `\alpha=111/\sqrt5` | yes |
| A6 | mathematica | 42 | `expectZero[chiLock - lambdaEll/2]` | `\chi_s=\Lambda_\ell/2` | yes |
| A7 | mathematica | 43 | `expectZero[kappaLock - (9/5)lambdaEll^2]` | `\kappa=(9/5)\Lambda_\ell^2` | yes |
| A8 | mathematica | 57 | `expectZero[chiRef - 37/2]` | `\chi_s=37/2` | yes |
| A9 | mathematica | 58 | `expectZero[kappaRef - 12321/5]` | `\kappa=12321/5` | yes |
| A10 | mathematica | 59 | `expectZero[alphaRef - 111/Sqrt[5]]` | `\alpha=111/\sqrt5` | yes |

All ten assertions are non-tautological. The load-bearing pair (A1/A2, A6/A7) derive `chi_lock`/`kappa_lock` by independent substitutions (`c_s -> hbar/(2 m ell)` then `L -> Lambda_ell ell`) and then check them against the claimed closed forms — the LHS is not defined as the RHS, so the check can genuinely fail if the healing substitution were wrong. The reference-branch checks (A3–A5, A8–A10) confirm the `Lambda_ell=37` numerical specialization of those derived forms; `chi_ref`/`kappa_ref`/`alpha_ref` come from substitution into the derived symbolic results, not from hardcoded literals.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.txt` (mtime 2026-05-27 02:15:34) vs script `..._sympy_audit.py` (mtime 2026-05-27 02:38:20)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.txt` (mtime 2026-05-27 02:18:07) vs script `..._mathematica_audit.wl` (mtime 2026-05-27 02:38:22)

**What's wrong:**
Both saved outputs predate their scripts by ~20–23 minutes. The content disagreement is confined to the banner label: the SymPy output line 3 reads `STAGE 57 — HEALING-LENGTH LOCK AND SUPPORT SCALE`, whereas the current script line 27 banner is `STAGE 074 — ...`; the Mathematica output line 3 reads `STAGE 057 — ...`, whereas the current `.wl` line 26 banner is `STAGE 074 — ...`. All numeric/symbolic result lines in the outputs (chi_s, kappa, alpha, the `= 0` residuals, `49.640709100495331260`) match what the current scripts emit, so only the stale banner differs.

**Why this matters:**
The committed transcript a reader inspects mislabels the stage as 57/057. The mismatch is cosmetic for the math but is exactly the script/output-band self-label drift this pipeline tracks; a fresh run will correct it.

**Required change:**
Re-run both scripts to refresh the committed `.txt` transcripts so the banner reads `STAGE 074`. (Done together with F2 so the docstring fix is captured in the same refresh.)

**Verification:**
After refresh, SymPy output line 3 and Mathematica output line 3 both read `STAGE 074 — HEALING-LENGTH LOCK AND SUPPORT SCALE`; all result lines unchanged (`chi_s = 37/2`, `kappa = 12321/5`, `alpha (numeric) = 49.640709100495331260`).

### F2 — paper_misalignment (subtype: notes_contradicts_script) — self-label only

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:3`

**What's wrong:**
The SymPy module docstring's first line is a stale self-label: `moving_throat_pde_stage57_family1_healing_lock_sympy_audit.py`. The file's actual name and the docstring's own next line (`SymPy audit for Stage 074:`, line 5) are `074`. This is the known EM-extension renumber drift (a `stage57`/`STAGE 57` remnant). It is an unambiguous self-label — the file is unmistakably stage 074 — so it is in-scope for the Reading-2 in-loop policy (verdict:findings ⇒ fix unambiguous self-labels), not a content/value disagreement.

**Why this matters:**
A reader scanning the docstring sees the wrong stage number. No math is affected. Note: this is *self-label* drift, not a genuine paper↔script value disagreement; it does not require user resolution and does not set `needs_user_resolution`. It is filed here for completeness of the self-label inventory and is corrected mechanically.

**Required change:**
At `..._sympy_audit.py:3` replace `moving_throat_pde_stage57_family1_healing_lock_sympy_audit.py` with `moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`.

**Verification:**
`..._sympy_audit.py:3` reads `moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`; `grep -n "stage57\|STAGE 57" ..._sympy_audit.py` returns nothing. (The `.wl` has no stale self-label: its banner and final print already read `STAGE 074` / `Stage 074`.)

## Independent-derivation check (Mathematica)

The `.wl` runs the same substitution chain as the `.py`: `chiFromHealing = (mpsi cSw len/hbar) /. cSw -> hbar/(2 mpsi ell)` (wl:33–36) mirrors `chi_in_ell = chi_def.subs(c_s, hbar/(2 m_psi ell))` (py:37); `chiLock = ... /. (len/ell) -> lambdaEll` (wl:37) mirrors `chi_lock = chi_in_ell.subs(L, Lambda_ell*ell)` (py:41); `kappaLock = 4 chiLock^2 + (4/5) lambdaEll^2` (wl:38) mirrors py:48. This is the *same* algebra in both engines. I considered `mathematica_transliteration` but did not file it: the stage's entire content is a single closed-form reduction `\chi=L/(2\ell)=\Lambda_\ell/2 \Rightarrow \kappa=(9/5)\Lambda_\ell^2`, for which there is no genuinely distinct derivation route — both engines must apply the same two substitutions to the same physical definition. The Mathematica side is a legitimate independent confirmation of an elementary identity, not an echo of nontrivial SymPy-specific scaffolding. (One implementation note, not a finding: wl:37 uses the rewrite rule `(len/ell) -> lambdaEll`; the committed output line 5 `chi_s (locked) = lambdaEll/2` confirms the rule fired and the result is correct.)

## Engine cross-check

Both outputs agree exactly on every result:

| quantity | SymPy output | Mathematica output |
|---|---|---|
| chi_s (locked) | `Lambda_ell/2` (l.6) | `lambdaEll/2` (l.5) |
| kappa(Lambda_ell) | `9*Lambda_ell**2/5` (l.7) | `(9*lambdaEll^2)/5` (l.6) |
| chi_s residual | `0` PASS | `0` PASS |
| kappa residual | `0` PASS | `0` PASS |
| chi_s (ref) | `37/2` (l.13) | `37/2` (l.13) |
| kappa (ref) | `12321/5` (l.14) | `12321/5` (l.14) |
| alpha | `111*sqrt(5)/5` (l.15) | `111/Sqrt[5]` (l.15) |
| alpha numeric | `49.640709100495331260` (l.16) | `49.64070910049533126...` (l.16) |
| ref residuals | all `= 0` PASS | all `= 0` PASS |

`111*sqrt(5)/5` and `111/Sqrt[5]` are the same number; the residual check `alpha_ref - 111/sqrt(5) = 0` passes in both. `engines_agree: true`.

## Verdict justification

The math holds up under attack. The healing-lock substitution chain is genuine (not tautological: the LHS quantities are derived by substitution, then checked against independently-stated closed forms), the carried-forward `4` and `4/5` coefficients are correctly anchored to Stage 071 in both the script comment and the notes and match the paper card, symbol positivity is physically justified in both engines, and every paper-side deliverable (`\chi_s=37/2`, `\kappa=12321/5`, `\alpha=111/\sqrt5`, plus the carried `\eta=37`) is exercised and reconciles to both `.tex` and `.md`. The two engines agree exactly. The only defects are low-severity script/output-band self-label drift: a stale `stage57` docstring line in the SymPy script (F2) and stale `STAGE 57`/`STAGE 057` banners in both committed transcripts (F1) — both are remnants of the known EM-extension renumber and self-correct on a docstring fix plus a fresh run. No `paper_misalignment` value/target disagreement exists, so `needs_user_resolution: false`. Verdict is `findings` solely on these two mechanical label/refresh items.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `chi_s = Lambda_ell/2` (symbolic lock) | py:50/53, wl:40/42; sympy.txt l.6, math.txt l.5 | tex:23 (locked form), md:71 | MATCH |
| `kappa = (9/5) Lambda_ell^2` (symbolic) | py:51/54, wl:41/43; sympy.txt l.7, math.txt l.6 | tex:29, md:91 | MATCH |
| `chi_s = 37/2` (reference) | py:63/68, wl:52/57; sympy.txt l.13, math.txt l.13 | tex:23, tex:34, md:75/99/129 | MATCH |
| `kappa = 12321/5` (reference) | py:64/69, wl:53/58; sympy.txt l.14, math.txt l.14 | tex:29, tex:34, md:95/103/130 | MATCH |
| `alpha = 111/sqrt(5) ≈ 49.6407091` | py:65/66/70, wl:54/55/59; sympy.txt l.15–16, math.txt l.15–16 | tex:31, md:117 | MATCH |
| `eta = 37` (carry-forward from Stage 073) | input `Lambda_ref=37` (py:56, wl:45) | tex:34, md:11/101/128 | MATCH |

INTERNAL scaffolding (no finding): banner strings; `chi (after healing-length substitution) = L/(2*ell)` intermediate; per-check residual `= 0` lines and PASS flags; `Lambda_ell = 37` echo print; `Final ledger:` recap; the `Stage 074 Mathematica audit passed.` flag.

All six emitted deliverable values reconcile to both the `.tex` card and the `.md` notes with no value disagreement. The only doc/script mismatches are the stale stage-number self-labels (docstring line + transcript banners), captured as F1/F2 above; those are label-only, not value misalignments, and do not trigger a `value_mismatch` paper_misalignment.

## Self-test notes

Checked: (1) variable independence — no `diff`/`D` in this stage, only substitutions, so the zero-derivative trap does not apply. (2) Symmetry/parity — no integrals. (3) Trivial-case pre-check — hand-evaluated the chain: `m c_s L/hbar` with `c_s=hbar/(2 m ell)` gives `L/(2ell)`, with `L=Lambda_ell ell` gives `Lambda_ell/2`; `kappa=4(Lambda_ell/2)^2+(4/5)Lambda_ell^2=Lambda_ell^2+(4/5)Lambda_ell^2=(9/5)Lambda_ell^2`; at 37: `kappa=9*1369/5=12321/5`, `sqrt=111/sqrt5≈49.6407091` — all match the asserted RHS, confirming the `expect_zero` residuals are genuinely zero (not zero-by-construction). (4) Path specifications — no missing-script finding. (5) Paper round-trip — F2 replaces a stale filename literal in a docstring with the correct one and introduces no new constant; F1 is a refresh-only directive; neither alters any asserted value, so no new paper_misalignment is created.
