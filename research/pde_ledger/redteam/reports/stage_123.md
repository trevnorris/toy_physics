---
unit_id: 123
batch: IV.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: misaligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage123_parent_normalized_branch_values.md
  paper_appendix: present
---

# Audit unit 123 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_123.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage123_parent_normalized_branch_values.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_123}` row near line 1280; the appendix carries no narrative summary of its own for this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py`
- mathematica: (missing — stage card line 11 explicitly states "Mathematica audit: none yet"; not filed as a finding per audit instructions)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.txt`
- mathematica output: (missing — no .wl)

## What the paper claims

The stage card (stage_123.tex) is terse: it advertises the "Output" only as "Converts the Family-1 compensated target to normalized parent flow and traction numbers (Xi_v, Xi_T)." The detailed boxed claims live in the notes file. Specifically the notes state (in their own LaTeX):

1. The mixed-flow definition `Xi_v := q_* sqrt(mu_0 m_psi rho_w) a L^{3/2} ell^{3/2} v_{w0} / (hbar sqrt(Z_q) c_s)` plus the exact branch law `Xi_v = -(3*sqrt(30)*pi^{3/2}/228) * r` (notes lines 25-31).
2. The numerical Family-1 value `Xi_v^{F1} ≈ -1.01675633282526` (notes lines 37-41).
3. The mouth-traction definition `Xi_T := T_m a sqrt(mu_0 rho_w L ell)/sqrt(Z_q m_psi)` plus the branch law `Xi_T = 3*sqrt(30)/(10*sqrt(pi)) * (1/g)` (notes lines 56-72).
4. The numerical values `Xi_T^nat ≈ 0.927058084855655`, `Xi_T^(-) ≈ 1.22297517701464`, `Xi_T^(+) ≈ 0.331334521644609` (notes lines 76-101).

There is an internal arithmetic inconsistency inside the notes themselves: the symbolic boxed formula on line 30 has denominator 228, but the numerical Xi_v^F1 quoted on line 41 is only consistent with denominator 160 (see Findings).

## What the script claims to verify

The script symbolically constructs `r_from_parent = lam / sqrt(Ks*Kq)` using stated upstream formulas `Ks = 3*pi*a^2 hbar^2/(5 m rho ell)`, `Kq = Z_q/mu_0 * pi^2 c_s^2/(4 L^2)`, `lam = -(8*sqrt(2)/3)*q*v_w0*a^2*ell*sqrt(L)`, solves `r = r_from_parent` for `v_w0`, substitutes into Xi_v_def, and asserts (line 46) `Xi_v_expr + 3*sqrt(30)*pi^(3/2)*r/160 == 0`. It then constructs `g_from_parent = sqrt(2 Z_q Ks)/(T_m J_s c_s sqrt(mu_0 L))` with `c_s -> hbar/(2 m ell)` applied, solves for T_m, substitutes into Xi_T_def, and asserts (line 47) `Xi_T_expr - 3*sqrt(30)/(10*sqrt(pi)*g) == 0`. It then prints the four Family-1 numeric values.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Xi_v = -(3*sqrt(30)*pi^{3/2}/228) * r` (notes:30) | `Xi_v_expr + 3*sqrt(30)*pi^(3/2)*r/160 == 0` (script:46) | **mismatch (228 vs 160)** |
| `Xi_v^F1 ≈ -1.01675633282526` (notes:41) | printout `Xi_v(F1) ≈ -1.0167563328252594644` (script:59; output:17) | match (numerically agrees with script's 160-denominator law, NOT with notes' 228 symbolic) |
| `Xi_T = 3*sqrt(30)/(10*sqrt(pi)*g)` (notes:70) | `Xi_T_expr - 3*sqrt(30)/(10*sqrt(pi)*g) == 0` (script:47) | match |
| `Xi_T^nat ≈ 0.92706` (notes:82) | printout `Xi_T(nat) ≈ 0.92705808485565499282` (script:60; output:18) | match |
| `Xi_T^(-) ≈ 1.22298` (notes:97) | printout `Xi_T(-)` ≈ 1.22297517701463916 (script:61; output:19) | match |
| `Xi_T^(+) ≈ 0.33133` (notes:99) | printout `Xi_T(+)` ≈ 0.33133452164460908 (script:62; output:20) | match |

Set `paper_alignment: misaligned` because the load-bearing symbolic identity for Xi_v in the notes disagrees with the script's symbolic identity; the rest matches.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 46 | `expect_zero("Xi_v law", Xi_v_expr + 3*sqrt(30)*pi^(3/2)*r/160)` | claim 1 (Xi_v law) | partial — script's denominator 160 disagrees with notes' 228 |
| A2 | sympy | 47 | `expect_zero("Xi_T law", Xi_T_expr - 3*sqrt(30)/(10*sqrt(pi)*g))` | claim 3 (Xi_T law) | yes |
| A3 | sympy | 59-62 | `print(...)` for `Xi_v_F1`, `Xi_T_nat`, `Xi_T_minus`, `Xi_T_plus` | claims 2, 4 | print-only (no `assert`); informational |

A3 is print-only, so the Family-1 / branch numerical values are not under an executable check. The numerics match the notes' boxed numbers though, so this is a low-severity gap (informational), not a separate failure.

## Findings

### F1 — paper_misalignment

**Subtype:** value_mismatch (paper-side notes have an internal typo that also disagrees with the script)
**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage123_parent_normalized_branch_values.md:25-31` (symbolic claim)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage123_parent_normalized_branch_values.md:37-41` (numerical claim)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py:46` (script assertion)

**What's wrong:**

The notes file claims, in a boxed equation (line 30):
> `Xi_v = -(3*sqrt(30)*pi^{3/2}/228) * r`

but later boxes (lines 37-41) state numerically:
> `Xi_v^{F1} = -(3*sqrt(30)*pi^{3/2}/228) * r_{F1} ≈ -1.01675633282526`

The script's symbolic derivation (lines 25-36) re-derives the law from the same upstream Ks, Kq, lam definitions and obtains
- `Xi_v(r) = -3*sqrt(30)*pi**(3/2)*r/160`

(see output line 13). Its assertion (line 46) therefore uses denominator 160. The notes' own numerical value -1.01675633282526 is consistent with the 160-denominator form only:

- Denominator 160: `-3*sqrt(30)*pi^(3/2) * r_{F1} / 160` with `r_{F1} = sqrt(12*(37/20)^2/pi^2 - 1)` evaluates to `-1.01675...` (matches the boxed numeric).
- Denominator 228: same product divided by 228 evaluates to roughly `-0.713...` (does NOT match the boxed numeric).

So the notes are internally inconsistent: the symbolic constant `228` is incompatible with the numeric `-1.01675633282526`. The script's value `160` is self-consistent with both its own derivation and with the notes' numeric.

**Why this matters:**

This is the load-bearing constant for the parent-flow branch law. Any downstream stage that quotes the symbolic form `Xi_v = -(3*sqrt(30)*pi^{3/2}/228)*r` will compute the wrong numeric. If `160` is the correct constant (as the derivation and the boxed numeric both indicate), the notes' `228` is a typo that risks propagating into the paper's main text or into a downstream stage's hardcoded numbers.

**Required change:**

This is a paper-side disagreement (notes ↔ script). Codex must NOT auto-edit. The directive contains a `## Resolve before fix_loop` block. The expected resolution is paper-side: the notes' boxed denominator should be corrected to `160` (matching the derivation and the boxed numeric). After user confirmation, no script change is needed.

**Verification:**

Once the user picks a direction, the verifier confirms either (a) `notes ... /stage123_*.md:30` now reads `/160` and the script assertion at line 46 is unchanged, OR (b) the script assertion at line 46 is changed to `/228` and the new sympy output disagrees with the notes' numeric `-1.01675...` (which would be the wrong direction — re-flag).

### F2 — paper_misalignment

**Subtype:** notes_contradicts_script (minor cosmetic banner mismatch)
**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py:16` (banner)
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.txt:11` (mirrors banner)
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_123.tex:1` (paper says "Stage~140")

**What's wrong:**

The script's banner reads
> `banner("STAGE 106 — PARENT-NORMALIZED BRANCH VALUES")`

The paper card section header (stage_123.tex:1) calls this "Stage~140". The label is `stage:123` and the file is `stage_123.tex`. So three different numbers are floating around: 106 (script banner), 123 (file/label), 140 (paper card section heading). This is cosmetic — no math depends on it — but the banner "STAGE 106" is dead wrong by any of the live conventions.

**Why this matters:**

Low impact, but the banner is the first thing in the script's output and is mirrored in the saved output file. If a future reader greps for "stage 123" or "stage 140" in the captured outputs they will not match this file. This is a label drift, not a math error.

**Required change:**

In `scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py:16`, replace
- `banner("STAGE 106 — PARENT-NORMALIZED BRANCH VALUES")`

with the matching stage id from the paper card and file name (e.g. `"STAGE 123 — PARENT-NORMALIZED BRANCH VALUES"`).

**Verification:**

After Codex applies, the sympy output's banner line will read `STAGE 123 — ...` and match the file name / paper card label. The assertions and numeric values must not change.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` script exists for this stage. The paper card explicitly states `"Mathematica audit: none yet"` (stage_123.tex:11). Per the audit instructions, no `missing_mathematica` finding is filed because no upstream-referenced result requires the Mathematica mirror for this unit.

## Engine cross-check

Not applicable — only one engine present.

## Verdict justification

The script's symbolic derivation of Xi_v from upstream Ks, Kq, lam is correct: starting from the stated Ks = 3πa^2ℏ^2/(5mρℓ), Kq = Z_q π^2 c_s^2 / (4 μ_0 L^2), lam = -(8√2/3) q v_w0 a^2 ℓ √L, the substitution Xi_v_def = q√(μ_0 m ρ) a L^{3/2} ℓ^{3/2} v_w0 / (ℏ √Z_q c_s) and solving v_w0 from r = lam/√(Ks Kq) yields exactly `-3√30 π^{3/2} r / 160`. I verified this longhand and the script's printed expression agrees. The Xi_T side, after the c_s → ℏ/(2mℓ) substitution applied to g_from_parent, likewise yields the stated `3√30/(10√π g)` cleanly. The four Family-1 numbers all match the notes to displayed precision.

So the script is mathematically clean. The ONLY substantive finding is that the notes' symbolic denominator `228` in the boxed Xi_v law contradicts BOTH the derivation AND the notes' own numeric. The cleanest reading is that `228` is a typo for `160` in the notes. This is a paper_misalignment requiring user resolution before any code change.

Verdict: `findings` with `stop_cold: null`. The pipeline is gated on user resolution of F1 before Codex applies F2.

## Self-test notes

- Variable independence: no `sp.diff` in the proposed changes; F2 is a banner string only. N/A.
- Symmetry/parity: no integrals affected by either finding. N/A.
- Trivial-case pre-check: re-evaluated `Xi_v(r_{F1})` for both denominators 160 and 228 with `r_{F1} = sqrt(12*(37/20)^2/pi^2 - 1) ≈ 1.778` — 160 gives ≈ -1.0168 (matches notes' numeric), 228 gives ≈ -0.713 (does NOT match). Confidence in F1 is high.
- Path specifications: no missing-script findings; no path manifest needed.
- Paper round-trip: F1 is a paper-side fix (user-routed). F2 is a script banner only; it does not introduce a new paper_misalignment because it merely aligns the script banner with the paper's existing labeling.
