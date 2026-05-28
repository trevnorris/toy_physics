---
unit_id: 138
batch: IV.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage138_normalized_mouth_gain_family.md]
  paper_appendix: present
---

# Audit unit 138 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_138.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage138_normalized_mouth_gain_family.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only line 1310 references this stage via `\input{stages/stage_138}`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage138_normalized_mouth_gain_family_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage138_normalized_mouth_gain_family_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage138_normalized_mouth_gain_family_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage138_normalized_mouth_gain_family_mathematica_audit.txt`

## What the paper claims

The stage card's body-equation block states verbatim: "Gain ratio \(R_q=(\mathfrak g_c-\mathfrak r)^2/(1+\mathfrak r^2)\); exact compensation gives \(R_q=1/4\)." The notes back this with the explicit derivation in normalized parent variables \(\mathfrak r=\lambda/\sqrt{K_sK_q}\), \(\mathfrak g_c=g_q\sqrt{K_s}/(g_s\sqrt{K_q})\), \(\Sigma_0=Lg_s^2/(K_s\Theta_\sigma)\), giving boxed deliverables \(M_s=\Sigma_0\), \(M_q=-\Sigma_0(\mathfrak g_c-\mathfrak r)^2/(1+\mathfrak r^2)\), and the exact ratio \(R_q=(\mathfrak g_c-\mathfrak r)^2/(1+\mathfrak r^2)\). Inserting the Stage 98 compensation surface \(\mathfrak g_c=\mathfrak r\pm\frac12\sqrt{1+\mathfrak r^2}\) yields \(R_q=1/4\). A corollary identification \(M_s=4\Sigma_m,\ M_q=-\Sigma_m\) with \(\Sigma_m:=\Sigma_0/4\) is presented as bookkeeping equivalence with Stage 237, not an independent claim. The Inputs and Checks fields list procedural imports/checklist items (outlet consistency, susceptibility closure, numerical fixed points), not new identities the script must verify in-stage.

## What the script claims to verify

Both engines start from the unnormalized \(M_s=Lg_s^2/(K_s\Theta)\), \(M_q=-L(K_sg_q-\lambda g_s)^2/(K_s(K_sK_q+\lambda^2)\Theta)\); apply substitutions \(\lambda\to\hat r\sqrt{K_sK_q}\) and \(g_q\to\hat g_c g_s\sqrt{K_q/K_s}\); compute \(R_q=-M_q^{\rm norm}/M_s^{\rm norm}\); then substitute the two compensation branches \(\hat g_c=\hat r\pm\frac12\sqrt{1+\hat r^2}\) and assert each yields \(R_q=1/4\). The Mathematica script additionally asserts \(R_q-(\hat g_c-\hat r)^2/(1+\hat r^2)=0\) explicitly via `expectZero`. The SymPy `Rq_raw.subs(Ms_norm, Sigma0)` step is a no-op on the simplified ratio (the \(\Sigma_0\)-factor has already cancelled out of \(R_q\)).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| \(R_q=(\mathfrak g_c-\mathfrak r)^2/(1+\mathfrak r^2)\) | Mathematica L45 `expectZero["R_q exact formula", rQ - (gCore - rHat)^2/(1 + rHat^2)]`; SymPy implicit via printed `Rq` (no `assert`) | match (Mathematica) / partial (SymPy: prints only) |
| \(R_q=1/4\) on \(+\) branch | SymPy L31 `assert sp.simplify(sol_plus - sp.Rational(1,4)) == 0`; Mathematica L51 `expectZero["R_q on + branch - 1/4", solPlus - 1/4]` | match |
| \(R_q=1/4\) on \(-\) branch | SymPy L32 `assert sp.simplify(sol_minus - sp.Rational(1,4)) == 0`; Mathematica L52 `expectZero["R_q on - branch - 1/4", solMinus - 1/4]` | match |
| \(M_s=\Sigma_0\) (definitional) | not asserted; both scripts print \(\Sigma_0\) after textual `subs(Ms_norm, Sigma0)` | partial — definitional, not load-bearing |
| \(M_q=-\Sigma_0(\hat g_c-\hat r)^2/(1+\hat r^2)\) | printed by both engines, derived from substitution; Mathematica's `R_q exact formula` check implicitly verifies this | match (via \(R_q\)) |
| \(M_s=4\Sigma_m,\ M_q=-\Sigma_m\) Stage 237 identification | not exercised | partial — trivial corollary of \(R_q=1/4\) plus \(\Sigma_m:=\Sigma_0/4\) definition; no independent claim |

Overall: paper-side load-bearing deliverables map to script-side checks. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 31 | `assert sp.simplify(sol_plus - sp.Rational(1,4)) == 0` | \(R_q=1/4\) on \(+\) branch | yes |
| A2 | sympy | 32 | `assert sp.simplify(sol_minus - sp.Rational(1,4)) == 0` | \(R_q=1/4\) on \(-\) branch | yes |
| A3 | mathematica | 45 | `expectZero["R_q exact formula", rQ - (gCore - rHat)^2/(1 + rHat^2)]` | closed-form \(R_q=(\hat g_c-\hat r)^2/(1+\hat r^2)\) | yes |
| A4 | mathematica | 51 | `expectZero["R_q on + branch - 1/4", solPlus - 1/4]` | \(R_q=1/4\) on \(+\) branch | yes |
| A5 | mathematica | 52 | `expectZero["R_q on - branch - 1/4", solMinus - 1/4]` | \(R_q=1/4\) on \(-\) branch | yes |

All assertions are non-tautological: each substitutes a concrete branch into a previously derived ratio and verifies a specific rational closed-form. None of them shortcut via assumption tricks.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage138_normalized_mouth_gain_family_mathematica_audit.wl:26`

**What's wrong:**
The Mathematica banner on L26 reads `banner["STAGE 121 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO"];`, but the unit ID is 138 (see filename, stage card `\label{stage:138}` at L2 of `paper/stages/stage_138.tex`, and Stage Title `Stage 138: Normalized Mouth-Gain Family and Compensation Ratio` at L1). The Mathematica banner therefore mis-identifies the audit's owning stage; the output transcript echoes the wrong stage number (`STAGE 121 — ...` at L11 of the `.txt`). Subtype: `notes_contradicts_script` is closest fit — the stage card and notes both say 138, the script's own banner says 121.

**Why this matters:**
A future reader scanning audit transcripts to triage failures could be misdirected to Stage 121 when the failing assertion lives in Stage 138's audit. It also weakens the bookkeeping trust that script identifiers line up with paper identifiers across the ledger.

**Required change:**
Edit `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage138_normalized_mouth_gain_family_mathematica_audit.wl:26`: change the banner string literal from `"STAGE 121 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO"` to `"STAGE 138 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO"`. This is a single-token in-place fix; no math changes.

**Verification:**
After Codex applies, re-running `redteam exec-mathematica 138` should produce a transcript whose banner line (currently `STAGE 121 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO`) reads `STAGE 138 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO`, all three `PASS:` lines remain present, and `EXIT_CODE: 0` is preserved.

## Independent-derivation check (Mathematica)

Both scripts traverse the same algebraic path: declare base symbols, define `M_s` and `M_q` in unnormalized form, apply the substitutions `lam -> rHat*Sqrt[kS*kQ]` and `gQ -> gCore*gS*Sqrt[kQ/kS]` (SymPy: L13-16; Mathematica: L36), then compute `-M_q^norm / M_s^norm`. The variable choreography is parallel (`Ms`/`mS`, `Mq`/`mQ`, `subs`, `Ms_norm`/`mSNorm`, `Mq_norm`/`mQNorm`, `Rq_raw`/`rQRaw`, `Rq`/`rQ`). Compare SymPy L10-19 against Mathematica L33-40 — the steps are 1:1 with renamed identifiers. However, the underlying derivation is essentially a one-line algebraic simplification (the cancellation \(K_s g_q-\lambda g_s = g_s\sqrt{K_sK_q}(\hat g_c-\hat r)\) is the only non-trivial step), so an "independent re-derivation" cannot meaningfully diverge from this path without introducing artificial complexity. Additionally, Mathematica adds a check the SymPy lacks (the explicit closed-form `expectZero` for `R_q - (gCore - rHat)^2/(1 + rHat^2)` at L45). Verdict: not flagged as `mathematica_transliteration` — the algebra is short enough that structural parallelism is forced, and the engines exercise the simplification independently via different simplifier internals (`sp.simplify` vs `FullSimplify[Together[Expand[...]]]`).

## Engine cross-check

SymPy output reports `R_q = (g_c - r)**2/(r**2 + 1)`, branches both `1/4`. Mathematica output reports `R_q = (gCore - rHat)^2/(1 + rHat^2)`, `R_q exact formula = 0`, both branches `1/4` with the residual `R_q on ± branch - 1/4 = 0` printed before `PASS:`. The two engines agree on every quoted residual and final ratio. `engines_agree: true`.

## Verdict justification

The script faithfully verifies the paper card's load-bearing deliverable (`R_q = (g_c - r)^2/(1 + r^2)`, collapsing to `1/4` on both compensation branches). Both engines agree, both outputs are fresh relative to script mtimes (output mtime > script mtime on both pairs), assertions are non-tautological, and no symbol-assumption issues hide branch errors (`r`/`g_c` declared `real`, all positives consistent). Attacks tried: (i) checked whether `Rq_raw.subs(Ms_norm, Sigma0)` could mask a `Σ_0`-dependent residual — it does not because the simplified `Rq_raw` has already cancelled `M_s_norm`; (ii) checked whether `g_c = r ± (1/2)sqrt(1+r^2)` and its alternate-sign cousin yield identical `R_q` for symmetry reasons (they do — `(g_c - r)^2 = (1+r^2)/4` for both signs, which is the substantive content); (iii) checked notes for additional boxed claims (`M_s = 4Σ_m, M_q = -Σ_m`) — these are trivial corollaries of `R_q = 1/4` and the definition `Σ_m := Σ_0/4`, not independent identities. The only finding is the Mathematica banner mislabel `STAGE 121` (should be `STAGE 138`); the math itself is correct. Verdict: `findings`, but no `stop_cold` flag — the issue is cosmetic, with no downstream propagation risk.

## Self-test notes

Traps checked: (1) `Rq_raw.subs(Ms_norm, Sigma0)` as potentially-tautological textual substitution — confirmed it is a no-op because the `Σ_0`-factor has already cancelled; (2) parity/symmetry of the two compensation branches — confirmed `(g_c - r)^2 = (1+r^2)/4` is identical for both `±` signs, so the two assertions A1/A2 (and A4/A5) are symmetric but each independently exercises a substitution; (3) banner-fix round-trip — confirmed editing a single string literal at L26 leaves all `expectZero[...]` calls and `Exit[0]` path untouched, no new paper-misalignment introduced (paper card and notes both refer to Stage 138).
