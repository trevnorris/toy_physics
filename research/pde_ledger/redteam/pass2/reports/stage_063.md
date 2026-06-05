---
unit_id: 063
batch: III.3
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
  notes_stage_files: [moving_throat_pde_stage063_parent_thresholds.md]
  paper_appendix: present
---

# Audit unit 063 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_063.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage063_parent_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 104; `\input` line 244)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage063_parent_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.txt`

## What the paper claims

Stage 063 (`\StatusExactClosure{}`) converts the Stage-062 parent-projected microscopic gain into exact fail/succeed surfaces. `\stagefield{Output}` reads: "Amplitude thresholds \eqref{eq:app-stage063-gphi-thresholds}, coherence thresholds \eqref{eq:app-stage063-coherence-thresholds}, and Cauchy no-go \eqref{eq:app-stage063-cauchy-nogo}." The three boxed deliverables are: (1) amplitude thresholds `g_{φ,fail}² = m c_{s,*}² K_X N_σσ G_fail/(ρ_* O_σφ²)` and the `G_suff` twin; (2) coherence thresholds `C_fail² = m c_{s,*}² K_X G_fail/(ρ_* g_φ² N_φφ)` and the `G_suff` twin; (3) the Cauchy no-go `C_fail² > 1 ⟹` no source profile in the channel reaches even the fail threshold. The notes add the load-bearing supporting identities: the parent gain `G_micro = ρ_* g_φ² O_σφ²/(m c_{s,*}² K_X N_σσ)`, the overlap factorization `O_σφ² = N_σσ N_φφ C_σφ²`, the best-case `G_max(g_φ) = ρ_* g_φ² N_φφ/(m c_{s,*}² K_X)` at `C_σφ²=1`, the ordering `C_fail² ≤ C_suff²`, and the `κ = K_X L²/T_X` insertion of `G_{fail,suff} = Pe_req/(κ Δ_{∞,0})` showing the explicit `K_X` cancels (surviving only inside `Δ` through `κ`). The notes attribute `G_fail/G_suff = Pe_req/(κ Δ)` and the phase diagram to Stage 061.

## What the script claims to verify

The SymPy script (banner "STAGE 63") pins the two amplitude thresholds and two coherence thresholds as symbolic expressions, then asserts a battery of structural identities: the overlap→coherence substitution `O_σφ² = C² N_σσ N_φφ` turns each `g²` threshold into its coherence-form twin (lines 51–60); the coherence-threshold ratio `C_suff²/C_fail² = G_suff/G_fail` (line 67); the factorization `G_micro = G_max·C²` (lines 71–74); the four `κ`-insertion identities for both fail/suff in both overlap and coherence form (lines 80–100); that an independent `sp.solve(G_micro = G_{fail,suff})` for `g_φ²` reproduces the hand-rearranged thresholds (lines 105–115); and the Cauchy saturation `G_max = G_micro|_{O²=N_σσ N_φφ}` (lines 119–122). The Mathematica script mirrors the same identity set but derives the `g_φ²` roots via `Reduce[... && gphiSq>0, …, Reals]` with positive-root `Cases` extraction instead of `Solve`, and constructs `c2Def = O²/(N_σσ N_φφ)` directly.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Amplitude thresholds `g_{φ,fail/suff}²` (eq gphi-thresholds) | py lines 45–48 print exactly `m c² K_X N_ss G/(ρ O²)`; reproduced by `solve` (105–115); wl 36–37 | match |
| Coherence thresholds `C_{fail/suff}²` (eq coherence-thresholds) | py lines 63–66 print exactly `m c² K_X G/(ρ g² N_pp)`; wl 51–52 | match |
| Overlap↔coherence equivalence (`O²=N_σσ N_φφ C²`) | py 51–60 / wl 42–49 substitution checks | match |
| `C_fail² ≤ C_suff²` ordering | py 67 ratio `= G_suff/G_fail` (with `G_fail≤G_suff`) / wl 55 | match (ordering inferred from positive ratio = `G_suff/G_fail`) |
| Cauchy no-go `C_fail²>1 ⟹ unreachable` | py 71–74 (`G_micro=G_max·C²`) + 119–122 (saturation at C²=1, i.e. `G_max` is max gain); wl 58, 102–105 | match (mathematical content: max reachable gain = `G_max`; the `>1` restatement is a logical corollary of `C²≤1`, not a separate identity) |
| `κ`-insertion + `K_X` cancellation | py 80–100 (4 checks) / wl 64–83 | match |
| `G_max = ρ_* g_φ² N_φφ/(m c² K_X)` | py 70, 119–122 / wl 57, 102–105 | match |

Every paper deliverable maps to a substantive, non-tautological script check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 51–55 | `expect_zero(g_fail² with O²→C²N_ssN_pp − coherence form)` | coherence thresholds / equivalence | yes |
| A2 | sympy | 56–60 | same for g_suff² | coherence thresholds | yes |
| A3 | sympy | 67 | `C_suff²/C_fail² − G_suff/G_fail == 0` | ordering `C_fail²≤C_suff²` | yes |
| A4 | sympy | 71–74 | `G_micro(O²→C²N_ssN_pp) − G_max·C² == 0` | Cauchy no-go (factorization) | yes |
| A5 | sympy | 80–84 | κ-insert fail (overlap form) `== 0` | κ-insertion / K_X cancel | yes |
| A6 | sympy | 85–89 | κ-insert suff (overlap form) `== 0` | κ-insertion / K_X cancel | yes |
| A7 | sympy | 91–95 | κ-insert fail (coherence form) `== 0` | κ-insertion / K_X cancel | yes |
| A8 | sympy | 96–100 | κ-insert suff (coherence form) `== 0` | κ-insertion / K_X cancel | yes |
| A9 | sympy | 105–111 | `solve(G_micro=G_fail) − g_fail² == 0` | amplitude thresholds (indep. solve) | yes |
| A10 | sympy | 113–115 | `solve(G_micro=G_suff) − g_suff² == 0` | amplitude thresholds | yes |
| A11 | sympy | 119–122 | `G_micro(O²→N_ssN_pp) − G_max == 0` | Cauchy saturation / no-go | yes |
| B1–B11 | mathematica | 42–105 | same identity set; `Reduce`+`Cases` for the two root checks | all of the above | yes |

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage063_parent_thresholds_sympy_audit.txt` (mtime 2026-05-22 19:43) vs `.py` (mtime 2026-06-03 15:59)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.txt` (mtime 2026-05-22 19:44) vs `.wl` (mtime 2026-06-03 15:59)

**What's wrong:**
Both committed transcripts predate their scripts and their content disagrees with what the current scripts emit. The current SymPy banner is `STAGE 63 — PARENT-OVERLAP THRESHOLD THEOREM` (`.py:31`) but the committed transcript opens `STAGE 46 — PARENT-OVERLAP THRESHOLD THEOREM` (`...sympy_audit.txt:3`). The current Mathematica banner is `STAGE 063 — PARENT THRESHOLDS` (`.wl:26`) but the committed transcript opens `STAGE 046 — PARENT THRESHOLDS` (`...mathematica_audit.txt:3`). The transcripts were captured before the banner-renumber edit, so they no longer reflect the current scripts.

**Why this matters:**
The committed transcript is the auditable record; a stale banner means the saved output cannot be trusted to reflect the current script state. The numeric/symbolic residuals themselves are unchanged (all `= 0`), so the math is not affected, but the freshness gate fails.

**Required change:**
Re-run both scripts and overwrite the two `.txt` files so the banner and footer match the current scripts.

**Verification:**
After re-run, `...sympy_audit.txt:3` reads `STAGE 063 …` (per F2 fix) and footer line reads `All Stage 063 symbolic checks passed.`; `...mathematica_audit.txt:3` reads `STAGE 063 — PARENT THRESHOLDS` and footer reads `Stage 063 Mathematica audit passed.`

### F2 — stale numbering self-labels (SymPy)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:3`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:31`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:124`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:76` (cross-reference, see note)

**What's wrong:**
The SymPy script carries the known `+17` pre-renumber drift in its own self-labels. The file is `..._stage063_...` and the banner intends stage 63, yet:
- `.py:3` docstring: `Stage 46 SymPy audit.` (self-label, should be 063; 46+17=63)
- `.py:124` footer: `print("\nAll Stage 46 symbolic checks passed.")` (self-label, should be 063)
- `.py:31` banner: `banner("STAGE 63 — …")` uses un-padded `63`; the canonical/`.wl` form is `063`.

These are unambiguous self-labels of the *current* stage and fall under the in-loop Reading-2 policy (verdict:findings ⇒ fix unambiguous self-labels + refresh outputs). The Mathematica `.wl` self-labels are already correct (`063` at `.wl:26` and `.wl:108`).

Additionally `.py:76` is a *cross-reference*: `# Insert Stage-44 threshold formulas`. The notes attribute `G_fail/G_suff = Pe_req/(κ Δ)` and the phase diagram to Stage 061 (notes lines 63, 92, 182–186), and 44+17=61, so the correct cross-reference is Stage 061. (The docstring at `.py:9`/`.py:8` also says "Stage-44 threshold formulas" implicitly via the comment chain — there is no number in the docstring item itself.)

**Why this matters:**
Stale stage numbers in script labels and the captured footer make the transcript self-contradictory (file says 063, body says 46) and propagate the +17 numbering drift the program is actively reconciling.

**Required change:**
- `.py:3`: `Stage 46 SymPy audit.` → `Stage 063 SymPy audit.`
- `.py:124`: `All Stage 46 symbolic checks passed.` → `All Stage 063 symbolic checks passed.`
- `.py:31`: `banner("STAGE 63 — PARENT-OVERLAP THRESHOLD THEOREM")` → `banner("STAGE 063 — PARENT-OVERLAP THRESHOLD THEOREM")` (3-digit padding to match `.wl` and filename).
- `.py:76`: `# Insert Stage-44 threshold formulas` → `# Insert Stage-061 threshold formulas` (cross-reference corrected per notes attribution).

**Verification:**
`grep -nE "Stage 46|Stage-44|STAGE 63[^0-9]" ...sympy_audit.py` returns no hits; the refreshed `.txt` (F1) shows `STAGE 063` banner and `All Stage 063 …` footer.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. It uses `FullSimplify[Together[Expand[...]]]` (`.wl:20–24`) where SymPy uses `simplify(expand(...))`, sets `$Assumptions` as a global positivity block (`.wl:29–33`) rather than per-symbol `positive=True`, and — most decisively — derives the `g_φ²` roots via `Reduce[(gMicro /. gPhi^2 -> gphiSq) == gFail && gphiSq > 0, gphiSq, Reals]` followed by `Cases[…, HoldPattern[gphiSq == rhs_] :> rhs]` (`.wl:88–91`), with an explicit comment that this differs from the SymPy `Solve` path (`.wl:85–86`). It also builds `c2Def = O²/(N_σσ N_φφ)` (`.wl:41`) and substitutes that, rather than SymPy's reverse substitution `O² → C² N_σσ N_φφ`. Same target identities, genuinely different machinery. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both engines assert the identical identity set and both committed transcripts show every residual `= 0` with `PASS`/passed footers (sympy `.txt:7–19`, mathematica `.txt:7–30`). The printed closed forms agree symbol-for-symbol: `g_(phi,fail)^2 = G_fail*K_X*N_ss*cs_star_sq*m/(O_sp**2*rho_star)` (sympy `.txt:5`) ≡ `(csStarSq*gFail*kX*m*nSS)/(oSP^2*rhoStar)` (mathematica `.txt:5`); `C_fail^2 = G_fail*K_X*cs_star_sq*m/(N_pp*g_phi**2*rho_star)` (sympy `.txt:9`) ≡ `(csStarSq*gFail*kX*m)/(gPhi^2*nPP*rhoStar)` (mathematica `.txt:11`). Engines agree. The only discrepancy between the two transcripts is the stale banner stage number (F1), which is a transcript-freshness issue, not a math disagreement.

## Verdict justification

The math is sound and exactly aligned with the paper. I attacked each assertion for tautology and found none: A1/A2 exercise the genuine algebraic content of `O_σφ² = N_σσ N_φφ C_σφ²` (the two sides are built from independent expressions, not defined-then-checked); A3 is a non-trivial ratio of two independently pinned thresholds; A4/A11 exercise the factorization and Cauchy saturation that carry the no-go's mathematical content; A9/A10 re-derive the thresholds by an independent `solve`/`Reduce` rather than restating them. Every paper deliverable has a matching, anchored check, and the printed closed forms match the boxed paper equations symbol-for-symbol. The only defects are housekeeping: the two committed transcripts are stale (banner says 46/046 while scripts say 63/063) and the SymPy script still carries +17 pre-renumber self-labels (lines 3, 31, 124) plus one stale cross-reference (line 76 → Stage 061). Both are low-severity, script-side, and mechanically fixable. Verdict: `findings`; no `material_change` to the verified physics.

## Value Reconciliation (pass-2 augmentation)

The scripts emit **no numeric constants** — every result is a closed-form symbolic threshold in free symbols. The labeled symbolic deliverables and their reconciliation against the paper card and notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `g_{φ,fail}² = m c² K_X N_σσ G_fail/(ρ_* O_σφ²)` | py:45–47, sympy.txt:5; wl:36, math.txt:5 | tex eq:app-stage063-gphi-thresholds (`.tex:17–18`); notes lines 30–31, 104–105 | MATCH |
| `g_{φ,suff}² = m c² K_X N_σσ G_suff/(ρ_* O_σφ²)` | py:46,48, sympy.txt:6; wl:37, math.txt:6 | tex eq (`.tex:20–21`); notes lines 33–34, 107–108 | MATCH |
| `C_fail² = m c² K_X G_fail/(ρ_* g_φ² N_φφ)` | py:63,65, sympy.txt:9; wl:51, math.txt:11 | tex eq:app-stage063-coherence-thresholds (`.tex:27`); notes line 130 | MATCH |
| `C_suff² = m c² K_X G_suff/(ρ_* g_φ² N_φφ)` | py:64,66, sympy.txt:10; wl:52, math.txt:12 | tex eq (`.tex:28`); notes line 132 | MATCH |
| `G_max = ρ_* g_φ² N_φφ/(m c² K_X)` (best-case gain) | py:70 (asserted, not printed); wl:57 | notes lines 60, 156 | MATCH (notes carrier; intermediate, card omits legitimately) |
| `G_micro = ρ_* g_φ² O_σφ²/(m c² K_X N_σσ)` | py:41; wl:35 | notes lines 7, 100 | MATCH (notes carrier) |
| κ-form fail `g²= m c² T_X N_σσ Pe_req/(ρ_* L² O_σφ² Δ_∞)` | py:83; wl:67 | tex prose absent; notes lines 71–72, 196–200 | MATCH (notes carrier) |
| κ-form suff `g²= m c² T_X N_σσ Pe_req/(ρ_* L² O_σφ² Δ_0)` | py:88; wl:71 | notes lines 74–75, 204–206 | MATCH (notes carrier) |

INTERNAL scaffolding (no prose expectation, no finding): `gphi_sq`/`gphiSq` solve dummy, `G_fail_sub`/`G_suff_sub` substitution placeholders, `kappa = K_X L²/T_X` substitution, all `expect_zero`/`expectZero` residuals (all `= 0`), pass/fail flags, `c2Def`.

reconciliation: complete; 8 deliverable values checked, 0 misaligned.

## Self-test notes

Checked the four self-test traps: (1) Variable independence — no `sp.diff`/`D[]` derivatives in either script, so the zero-derivative trap does not apply; the prescribed fixes (F1 refresh, F2 label edits) introduce no derivatives. (2) Symmetry/parity — no integrals; not applicable. (3) Trivial-case — F2 edits are pure string/label changes to docstring/banner/comment/footer and do not touch any assertion expression, so all residuals remain `= 0`; I confirmed the changed lines (3, 31, 124, 76) are non-executable strings/comments. (4) Paper round-trip — the label fixes carry no constant and the cross-reference correction (Stage-44 → Stage-061) matches the notes' own attribution (notes lines 63, 182), introducing no new `paper_misalignment`. No directive trap surfaced.
