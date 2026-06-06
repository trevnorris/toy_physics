---
unit_id: 088
batch: III.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage088_loading_ratio_from_minimal_module.md]
  paper_appendix: present
---

# Audit unit 088 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_088.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage088_loading_ratio_from_minimal_module.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 154; section header row 4)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.txt`

## What the paper claims

Stage 088 reads the minimal isotropic conservative quadrupole precursor `Y_Q^cons(omega) = c0 + c1/(1 - omega^2/Omega_Q^2)` with the fixed coefficients `c0 = 3/4`, `c1 = 1/4` as a contact-plus-pole decomposition, and extracts the loading ratio. The `\stagefield{Output}` is the minimal-module ratio box, eq:app-stage088-rho: "\(\rho_\alpha=\frac{1}{c_0}=\frac43, \zeta_{\rm req}=\frac{c_1}{c_0}=\frac13, \Pi_{\rm tr}=\frac43C_{\rm mix}\)." The distinct deliverables are: (1) the inverse map `rho_alpha = 1/c0`, `zeta_req = c1/c0`; (2) the numeric results `rho_alpha = 4/3`, `zeta_req = 1/3` from `c0=3/4, c1=1/4`; (3) the product-language statement `Pi_tr = (4/3) C_mix`; and (4) the regime classification `C_mix < Pi_tr < 2 C_mix`, equivalently `0 < zeta_req = 1/3 < 1` (symmetric-lowest-twin). The notes (sections 1-5) add the contact-plus-pole construction `c0 = 1/rho_alpha`, `c1 = (rho_alpha-1)/rho_alpha` with `c0 + c1 = 1`, and the equivalence to the `alpha_req/alpha_mix` loading variable. `Omega_Q = 3 c_s/(2a)` is an INPUT carried forward (not a stage-088 deliverable); the loading-ratio result is independent of its value.

## What the script claims to verify

The SymPy script (docstring + assertions) verifies: (1) the loading-variable form equals the rho-parameterized form under `alpha_req -> rho_alpha*alpha_mix` (line 41-44); (2) `(c0,c1)` *extracted* by pole-residue limit at u=1 plus subtraction at u=0 equal `1/rho_alpha` and `(rho_alpha-1)/rho_alpha` (lines 63-67) — explicitly engineered to avoid the "define then assert" tautology (comment lines 46-58); (3) contact-plus-pole reconstruction and `c0+c1=1` (lines 71-75); (4) the inverse maps round-trip (lines 85-90); (5) the paper-form extraction `c0_paper=3/4, c1_paper=1/4` → `rho_min=4/3, zeta_min=1/3` (lines 95-109); (6) `Pi_tr = (4/3) C_mix` and the regime bounds `1 < rho_min < 2` (lines 116-119). The Mathematica script independently extracts `(c0,c1)` directly from the paper `yQpaper` form (different substitution path), checks `rho_min=4/3, zeta_min=1/3`, reconstructs both the coefficient form and the rho-parameterized form against the paper precursor, and checks the product identity and regime bound.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `rho_alpha = 1/c0` (inverse map) | py:77,85 `rho_from_c0`, round-trip; wl:64 `rhoMin=1/c0Paper` | match |
| `zeta_req = c1/c0` (inverse map) | py:79,87-90; wl:65 `zetaMin=c1Paper/c0Paper` | match |
| `rho_alpha = 4/3` (numeric) | py:108 `rho_min - 4/3`; wl:70 `rho_min - 4/3` | match |
| `zeta_req = 1/3` (numeric) | py:109 `zeta_min - 1/3`; wl:71 `zeta_min - 1/3` | match |
| `Pi_tr = (4/3) C_mix` (product) | py:116; wl:90 | match |
| `C_mix < Pi_tr < 2 C_mix` / `0<zeta<1` (regime) | py:118-119 `1<rho_min<2`; wl:91 `1<rhoMin<2` | match |
| `c0=3/4, c1=1/4` (input coeffs) | py:99-100; wl:59-60 | match |
| contact-plus-pole `c0=1/rho`, `c1=(rho-1)/rho`, `c0+c1=1` (notes) | py:63-75; wl:61,77,83 | match |

`paper_alignment: aligned`. Every paper-side deliverable maps to a non-tautological script-side check that is faithfully covered by both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41-44 | `expect_zero(Y_loading.subs - Y_rho)` | loading-variable ≡ rho-form (notes §1) | yes |
| A2 | sympy | 63 | `c0_extracted - 1/rho_alpha == 0` | contact fraction = 1/rho (notes §1) | yes (extraction ≠ definition) |
| A3 | sympy | 64-67 | `c1_extracted - (rho_alpha-1)/rho_alpha == 0` | pole residue = (rho-1)/rho (notes §1) | yes |
| A4 | sympy | 71-74 | `Y_rho - (c0_ext + c1_ext/(1-u)) == 0` | contact-plus-pole reconstruction | yes |
| A5 | sympy | 75 | `c0_ext + c1_ext - 1 == 0` | normalized static limit `c0+c1=1` | yes |
| A6 | sympy | 85-90 | inverse-map round-trips | `rho=1/c0`, `zeta=c1/c0` inverse maps | partial (round-trip on independently-derived c) |
| A7 | sympy | 99-100 | `c0_paper-3/4==0`, `c1_paper-1/4==0` | input coeffs from paper form | yes |
| A8 | sympy | 108-109 | `rho_min-4/3==0`, `zeta_min-1/3==0` | numeric loading ratio (Output) | yes |
| A9 | sympy | 116 | `Pi_tr_from_rho - (4/3)C_mix == 0` | product-language `Pi_tr=(4/3)C_mix` | yes (anchored to derived rho_min) |
| A10 | sympy | 118-119 | `rho_min>1`, `rho_min<2` | regime classification | yes |
| B1 | mathematica | 59-61 | `c0_paper-3/4`, `c1_paper-1/4`, `c0+c1-1` | input coeffs + static-limit | yes |
| B2 | mathematica | 70-71 | `rho_min-4/3`, `zeta_min-1/3` | numeric loading ratio (Output) | yes |
| B3 | mathematica | 77 | paper form − reconstruction from (c0,c1) | contact-plus-pole reconstruction | yes |
| B4 | mathematica | 83 | rho-parameterized form (rho→rho_min) − paper form | rho-form ≡ paper precursor | yes |
| B5 | mathematica | 90 | `Pi_tr_from_rho - (4/3)C_mix` | product identity | yes (anchored to derived rho_min) |
| B6 | mathematica | 91 | `1 < rho_min < 2` | regime classification | yes |

A6 is the only "partial" row: lines 85-90 compose the construction map with its inverse (`1/c0` with `c0=c0_extracted`). Because `c0_extracted` was *independently derived* (residue limit + subtraction, not defined as `1/rho_alpha`), the round-trip does exercise the inverse-formula claim, but it is weaker than A2/A3, which already establish the load-bearing identity directly. Not a finding.

## Findings

### F1 — stale_output (stale self-label, numbering)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py:3`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py:5`

**What's wrong:**
The SymPy docstring carries the pre-renumber self-label. Line 3: `moving_throat_pde_stage71_loading_ratio_from_minimal_module_sympy_audit.py`; line 5: `SymPy audit for Stage 71.` The file's canonical stage number is 088 (filename, paper card `stage:088`, manifest, banner on line 29 all say 088). "Stage 71" / "stage71" is a stale OWN-number self-label from the known +17 EM-extension renumber drift (088 − 17 = 071). This is the file labeling itself, not a cross-reference, so it is in-scope.

Note (deferred, NOT a finding): the `stage085` references at py:112 and wl:87 are CROSS-references to another stage (the upstream Stage-085 product identity) and are correctly numbered; cross-reference numbering is out of scope per the deferred numbering plan.

**Why this matters:**
A reader cross-checking the script against the paper card or manifest sees "Stage 71" and may distrust the file or mis-route it. No math impact.

**Required change:**
At py:3 replace `moving_throat_pde_stage71_loading_ratio_from_minimal_module_sympy_audit.py` with `moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py`. At py:5 replace `SymPy audit for Stage 71.` with `SymPy audit for Stage 088.`

**Verification:**
`grep -n "stage71\|Stage 71" scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py` returns nothing; the banner (line 29) and filename already read 088. Re-running the script must still exit 0 with the identical output transcript (docstring change is non-functional).

## Independent-derivation check (Mathematica)

The `.wl` is a genuine independent derivation, not a transliteration. The SymPy script first builds the `alpha_req/alpha_mix` loading form and the `rho_alpha`-parameterized `Y_rho`, extracts `(c0,c1)` via `Y_rho.subs(omega**2, u*Omega_Q**2)` (substituting `omega**2`, py:56), then maps to the paper form. The Mathematica script skips the loading-variable / rho-parameterized scaffolding entirely and works straight from the paper `yQpaper = 3/4 + (1/4)/(1 - omega^2/omegaQ^2)` (wl:46), using the *combined-ratio* substitution `yQpaper /. omega^2/omegaQ^2 -> u` (wl:47) — the opposite substitution choice. Both then extract `c1` via `Limit[(1-u)*…, u->1]` and `c0` by subtraction, but on different starting expressions. The `.wl` additionally checks the rho-parameterized rebuild (wl:82-84) as a separate confirmation. Quoted correspondence:
- SymPy: `Y_rho_u = sp.simplify(Y_rho.subs(omega**2, u * Omega_Q**2))` (py:56) vs Mathematica: `yQpaperU = yQpaper /. omega^2/omegaQ^2 -> u` (wl:47) — different operand and different starting form.
- SymPy derives from `rho_alpha` symbolically then specializes to the paper numbers; Mathematica derives the paper numbers directly. Not line-by-line parallel.
No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on every shared deliverable:
- `c0 = 3/4`, `c1 = 1/4` (sympy out lines 20-21; math out lines 11,13).
- `rho_alpha = 4/3`, `zeta_req = 1/3` (sympy out lines 22-25; math out lines 16-21).
- `Pi_tr = (4/3) C_mix` (sympy out line 26; math out lines 26-28).
- regime `1 < rho < 2` (sympy assert lines 118-119; math out lines 29-30).
`engines_agree: true`.

**Special-care verification (first-pass fragilities, both confirmed FIXED):**
1. SymPy vacuous-substitution risk: the previous failing combined-ratio sub `omega**2/Omega_Q**2 -> u` was replaced by `Y_rho.subs(omega**2, u*Omega_Q**2)` (py:56), with the prior failure documented in the comment (py:53-54). After `sp.simplify`, `Y_rho` is `(Omega_Q**2*rho_alpha - omega**2)/(rho_alpha*(Omega_Q**2 - omega**2))` (out line 6) — `omega**2` appears atomically, so the substitution lands. The output confirms a genuine extraction: `c0_extracted = 1/rho_alpha`, `c1_extracted = (rho_alpha-1)/rho_alpha` (out lines 8-9), not `0` (which is what a no-op sub would have produced for the residue). NOT vacuous.
2. Mathematica comment-close risk: the old bug was a comment containing the literal `stage085_*)` which closed `(* … *)` early. The current wl:87 reads `in the stage 085 Mathematica audit files). Substitute rho_min. *)` — the upstream reference is written "stage 085" (space, no `_*`), so no embedded `*)`. `grep` for `*)` shows only the seven legitimate terminal comment-closers (wl:38,45,51,63,75,81,87). No early-close. Confirmed by output completeness below.
3. PASS-count vs assertion-count: the `.wl` contains exactly 9 assertion calls (expectZero at wl:59,60,61,70,71,77,83,90; expectTrue at wl:91). The committed output contains exactly 9 `PASS:` lines (out:11,13,15,19,21,23,25,28,30) plus the final "Stage 088 Mathematica audit passed." NO assertion is silently skipped; NO partial run. (The `Limit::alimv` warning on out:6 is the standard benign "assumptions involving the limit variable ignored" notice for the `u->1` limit; the residue still evaluates correctly to 1/4 on out:9.)

## Verdict justification

`findings` — one low-severity stale self-label only; no math defect. I attacked the obvious failure modes and all held: (a) the SymPy `(c0,c1)` extraction is a real residue/subtraction probe, not a "define then assert" tautology — confirmed by the non-trivial extracted forms in the output; (b) the round-trip checks (A6) compose construction with inverse on *independently-derived* coefficients, and are backed by the stronger A2/A3/B-row direct checks, so not self-cancelling; (c) the two flagged first-pass fragilities are both genuinely fixed (vacuous substitution avoided via `omega**2` operand; comment no longer self-closes), and the PASS-line count matches the assertion count, ruling out a silent partial run; (d) the Mathematica derivation is independent (different starting form and substitution), not a transliteration; (e) `Omega_Q = 3c_s/(2a)` is an input, not a stage-088 deliverable, and is correctly carried as an opaque positive symbol. Paper alignment is exact on all four deliverables across both engines. I read the paper card, notes, and appendix row; the script's verified claim matches the paper's stated claim.

## Value Reconciliation (pass-2 augmentation)

Every result value the scripts emit, located in the `.tex`/`.md`:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `c0 = 3/4` | py:99, wl:59; sympy out:20, math out:9,11 | `.tex:20` (`c_0=\frac34`), `.md:83` (`c0 = 3/4`) | MATCH |
| `c1 = 1/4` | py:100, wl:60; sympy out:21, math out:9,13 | `.tex:22` (`c_1=\frac14`), `.md:84` (`c1 = 1/4`) | MATCH |
| `c0 + c1 = 1` | py:75, wl:61; sympy out:13, math out:15 | `.md:65` (`c0 + c1 = 1`) | MATCH |
| `rho_alpha = 4/3` | py:108, wl:70; sympy out:22,24, math out:16,19 | `.tex:37` (`\rho_\alpha=\frac{1}{c_0}=\frac43`), `.md:88,167` | MATCH |
| `zeta_req = 1/3` | py:109, wl:71; sympy out:23,25, math out:17,21 | `.tex:39,48` (`\zeta_{\rm req}=\frac{c_1}{c_0}=\frac13`), `.md:90,171` | MATCH |
| `Pi_tr = (4/3) C_mix` | py:116, wl:90; sympy out:26, math out:26,28 | `.tex:41` (`\Pi_{\rm tr}=\frac43C_{\rm mix}`), `.md:103,172` | MATCH |
| regime `C_mix < Pi_tr < 2 C_mix` (`1<rho<2`) | py:118-119, wl:91; math out:29-30 | `.tex:46` (`C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}`), `.md:121` | MATCH |
| inverse maps `rho_alpha = 1/c0`, `zeta_req = c1/c0` | py:14-19, wl:64-68; sympy out:14,16, math out:16-17 | `.tex:37,39` (`=\frac{1}{c_0}`, `=\frac{c_1}{c_0}`), `.md:71,73` | MATCH |

INTERNAL scaffolding (accounted for, no finding expected in prose): `Y_loading(omega)` and `Y_rho(omega)` simplified rational forms (sympy out:5-6) — intermediate constructions; `c0_extracted`/`c1_extracted` symbolic forms (sympy out:8-9) — intermediate extraction outputs feeding asserts; `rho_from_c0 = 1/c0`, `rho_from_c1 = -1/(c1-1)`, `zeta_from_c = c1/c0` printouts (sympy out:14-16) — symbolic map definitions; all residual `= 0` / `PASS` lines and the `Limit::alimv` warning — verification scaffolding. `Omega_Q = 3c_s/(2a)` is a carried-forward INPUT, not emitted as a stage-088 result, and is held as an opaque symbol; not a reconciliation item.

reconciliation: complete; 8 deliverable values checked, 0 misaligned

## Self-test notes

I checked: (1) variable-dependence/substitution traps — the SymPy `omega**2 -> u*Omega_Q**2` and Mathematica `omega^2/omegaQ^2 -> u` substitutions both genuinely land (operands present in their respective starting forms; outputs show non-trivial extracted residues, not the `0` a no-op would yield); (2) the Mathematica comment block never self-closes (no embedded `*)`; 9 assertions ↔ 9 PASS lines, no partial run); (3) the `(c0,c1)` extraction is non-tautological (residue-limit + subtraction, not "define then assert") and the round-trip A6 is backstopped by direct A2/A3/B-row checks; (4) the proposed fix is docstring-only (py:3,5 stale `Stage 71` → `Stage 088`) and is non-functional, so it cannot introduce a new paper_misalignment and the output transcript is unchanged. Conclusion: one low-severity stale self-label; math is sound and paper-aligned across both engines.
