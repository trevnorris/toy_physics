---
unit_id: 042
batch: III.1
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
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage042_rank2_selected_mode_normalization.md]
  paper_appendix: present
---

# Audit unit 042 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_042.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage042_rank2_selected_mode_normalization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 62 + lead-in line 16 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` (stage_042.tex:41): *"The rank--2 eigenvector \eqref{eq:app-stage042-evec-rank2} and overlap formulas \eqref{eq:app-stage042-z-overlap}--\eqref{eq:app-stage042-s-overlap}."* The card proves three boxed closed forms for the selected lower mode of the rank-2-completed wall operator with three distinct directions (mixed/outgoing `q`, support `r`, source `t`): (1) the eigenvector ratio `e1/e0 = [m(q-r)+r·xi]/[delta+xi-m·q(q-r)]`; (2) the loading overlap `(z·e_-)²/z0² = [delta+(1+qr)xi]²/D_{q,r}`; and (3) the source overlap `(s·e_-)²/s0² = [delta+(1+rt)xi-m(q-r)(q-t)]²/D_{q,r}`, with `D_{q,r}=[delta+xi-mq(q-r)]²+[m(q-r)+r·xi]²`. The one-vector law is recovered only at `r=q`. The notes (`...normalization.md`) extend this with four further deliverables not in the terse card: the combined normalization function `F_(q,r,t)` (§3), the tracking collapse to the Stage-040 law at `r=q` (§4), the source-tied specialization `F_src` with `t²=lambda_0=2/9`, `q=t·R_U`, `r=t` (§5), the flat-U recovery `F_src(R_U=1)=F_flat` (§6), and the first-order deformations `H_n^(src)`, `H_F^(src)` about `R_U=1` (§6). Appendix row 62 status: `\StatusExactClosure{}`, "Selected eigenvector and source/loading overlaps with distinct outgoing, support, and source directions."

## What the script claims to verify

Both engines verify the same six-section ledger, matching the notes one-for-one. §25.1: builds `n_req` independently, then checks both rows of the singular-matrix eigen-equation give the same `e1/e0` and that it equals the boxed `ratio_expected` (two independent constructions). §25.2: builds the two overlaps from `ratio_expected` and matches them against the closed-form `Z_expected`/`S_expected`, then matches `F_general = Z·S/(1-xi)` against `F_expected`. §25.3: `r→q` collapse equals the Stage-040 ("Stage 23") two-vector law. §25.4: source-tied substitution `q→√λ0·R_U, r→√λ0, t→√λ0` of the general `F_expected` equals the hand-written `F_src`. §25.5: `F_src(R_U=1)=F_flat`. §25.6: `H_n^(src)`/`H_F^(src)` derivatives at `R_U=1` equal the boxed forms, and the `O(eps)` series of `n_src` and `F_src/F_flat` match the linearized expansions. All checks are residual-to-zero `expect_zero`/`expectZero` assertions that hard-fail (raise / `Exit[1]`) on nonzero.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `e1/e0` ratio (eq evec-rank2) | sympy 72-73 / wl 56-57 (row1 & row2 vs expected) | match |
| `D_{q,r}` denominator (eq Dqr) | sympy 77 / wl 61 (used inside Z/S expected) | match |
| `(z·e_-)²/z0²` (eq z-overlap) | sympy 89 / wl 80 | match |
| `(s·e_-)²/s0²` (eq s-overlap) | sympy 90 / wl 81 | match |
| `r=q` one-vector recovery (card line 39) | sympy 112 / wl 94 (tracking collapse) | match |
| `F_(q,r,t)` (notes §3) | sympy 101 / wl 82 | match |
| `F_src` source-tied (notes §5) | sympy 130 / wl 112 | match |
| flat-U `F_src(R_U=1)=F_flat` (notes §6) | sympy 138 / wl 120 | match |
| `H_n^(src)`, `H_F^(src)`, series (notes §6) | sympy 152/165/171/175 / wl 146-154 | match |

Every paper-card boxed deliverable and every notes-enumerated deliverable has a faithful, non-tautological script-side counterpart. No `mismatch`, no `missing`, no `extra`. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 72 | `simplify(ratio_row1 - ratio_expected)==0` | eigenvector ratio (row1) | yes |
| A2 | sympy | 73 | `simplify(ratio_row2 - ratio_expected)==0` | eigenvector ratio (row2) | yes |
| A3 | sympy | 89 | `Z_overlap - Z_expected ==0` | z-overlap | yes |
| A4 | sympy | 90 | `S_overlap - S_expected ==0` | s-overlap | yes |
| A5 | sympy | 101 | `F_general - F_expected ==0` | F_(q,r,t) (notes §3) | yes |
| A6 | sympy | 112 | `F_track - F_stage23 ==0` | r=q collapse | yes |
| A7 | sympy | 130 | `F_src_direct - F_src ==0` | F_src (notes §5) | yes |
| A8 | sympy | 138 | `F_src(R_U=1) - F_flat ==0` | flat-U recovery | yes |
| A9 | sympy | 152 | `H_n_src - H_n_expected ==0` | H_n^(src) deformation | yes |
| A10 | sympy | 165 | `H_F_src - H_F_expected ==0` | H_F^(src) deformation | yes |
| A11 | sympy | 171/175 | series-vs-linear `==0` | first-order expansion | yes |
| B1–B11 | mathematica | 56–154 | `expectZero[...]` mirrors of A1–A11 | same claims | yes |

No tautological rows: every `expected` form on the RHS is written from the paper/notes closed form, while the LHS is built by an independent route (singular-matrix rows, overlap construction from `ratio_expected`, substitution into the general `F`, symbolic derivative). The residual being zero is a real identity, not algebraically forced by construction.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.txt` (mtime 2026-05-26 02:03)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.txt` (mtime 2026-05-26 02:03)
- vs scripts: `.py` and `.wl` both mtime 2026-06-03 15:59

**What's wrong:**
Both committed transcripts predate the current scripts. The content also differs at the label level: the sympy `.txt` banner reads `STAGE 25 — …` and `STAGE 25 THEOREM LEDGER` (lines 3, 134) while the current `.py` source banner reads `STAGE 42 — …`/`STAGE 42 THEOREM LEDGER` (lines 51, 180); the mathematica `.txt` banner reads `STAGE 025` (line 3) while the current `.wl` source banner reads `STAGE 042` (line 33). The sub-section banners (`25.1`–`25.6`), docstring "Stage 25", and the "Stage-24"/"Stage-23" cross-references are still the pre-renumber labels in BOTH the script sources and the saved outputs. The MATHEMATICS is unaffected — all residuals are 0 and all `PASS` lines are present; only the stage-number labels are stale.

**Why this matters:**
The captured transcript no longer matches the banner the current script emits, and the embedded stage labels (25/24/23) are off by the project's known +17 renumber. Purely cosmetic for correctness, but it is a genuine freshness + stale-label mismatch.

**Required change:**
This finding is **informational, not directive-driven**. Per the project's standing decision (user, 2026-06-04), the SCRIPT/OUTPUT-band numbering cleanup — including these `25.x`/Stage-24/Stage-23 self-labels and the pre-renumber committed transcripts — is **deferred to a dedicated separate pass** (`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`, PENDING, content-keyed). No per-stage directive should hand-edit labels here. The orchestrator's standard independent exec re-run will regenerate fresh transcripts (banner → `STAGE 42`/`STAGE 042`); the residual stale sub-labels remain the dedicated pass's responsibility.

**Verification:**
After the orchestrator's exec re-run, the regenerated `.txt` files carry the `STAGE 42`/`STAGE 042` top banner and all 22 residual lines = 0 / `PASS`. The remaining `25.x`/Stage-24/Stage-23 sub-labels are tracked by the deferred numbering pass, not this audit.

## Independent-derivation check (Mathematica)

The `.wl` is a close structural mirror of the `.py`: identical variable names (`nReq`↔`n_req`, `dQR`↔`D_qr`, `zOverlap`↔`Z_overlap`, `fSrc`↔`F_src`, `hNSrc`↔`H_n_src`), identical six-section order, identical construction choreography (e.g. wl:62 `zOverlap = (1 + q ratioExpected)^2/(1 + ratioExpected^2)` ↔ py:78 `Z_overlap = sp.simplify((1 + q*ratio_expected)**2/(1 + ratio_expected**2))`; wl:106-107 `fSrcDirect = fExpected /. {q -> Sqrt[lambda0] rU, …}` ↔ py:127-128 `F_src_direct = F_expected.subs({q: sp.sqrt(lam0)*R_U, …})`). This is the project's accepted MATHEMATICA_MIRROR pattern: the two engines re-establish the SAME closed-form identities using two independent CAS simplification kernels (`FullSimplify[Together[Expand[...]]]` vs `sp.simplify(sp.expand(...))`). Because the stage's deliverables ARE closed-form identities, an "independent route" means an independent algebra engine confirming the residual is zero, which is what the `.wl` does. No `mathematica_transliteration` finding is raised: the mirror is policy-compliant and each assertion is substantive (residual-to-zero, hard-fail), not a re-print of the SymPy intermediate values.

## Engine cross-check

Both engines pass all checks with zero residuals:
- SymPy output (`...sympy_audit.txt`): all 22 residual lines = `0` (e.g. line 24 `row1 - expected = 0`, line 74 `F_general - expected = 0`, line 130 `F_src/F_flat - linear expansion = 0`).
- Mathematica output (`...mathematica_audit.txt`): every `expectZero` printed `= 0` and emitted `PASS:` (e.g. line 13 `PASS: row1 - expected`, line 28 `PASS: F_general - expected`, line 62 `PASS: F_src/F_flat - linear expansion`), terminating with `Stage 042 Mathematica audit passed.`

The two engines also agree on the closed forms they print: e.g. SymPy `H_n^(src) = -4·m·ξ/(9·δ + 11·ξ)` (output lines 118-121) vs Mathematica `H_n^(src) = (-4*m*xi)/(9*delta + 11*xi)` (output line 53) — identical after `lam0=2/9` clears denominators. `engines_agree: true`.

## Verdict justification

Verdict `findings` solely because of the informational `stale_output` (F1); the mathematics is clean and the paper alignment is exact. Attacks tried and failed: (1) tautology check — the `_expected` RHS forms are paper/notes closed forms while the LHS are independently constructed (singular-matrix rows, overlap built from `ratio_expected`, substitution of the general `F`, symbolic `diff`), so each residual-to-zero is a real identity, not forced by construction; (2) the §25.6 derivatives are not the zero-derivative trap — `n_src`/`fRatio` genuinely depend on `R_U`, so `d/dR_U` is non-trivial and the assertions can fail; (3) symbol domains — `xi`, `delta`, `R_U`, `eps` declared positive/real, consistent with the physical setup (`xi` a softening fraction in (0,1), `R_U`>0, `eps` small positive), no positivity is exploited to mask a branch error since the checks are pure rational identities; (4) every paper-card boxed deliverable (eigenvector + two overlaps + `D_{q,r}`) and every notes deliverable (`F`, tracking collapse, `F_src`, flat-U, two deformations, series) maps to a substantive check; (5) value `lambda_0=2/9` matches across script (`Rational(2,9)`) and notes. I confirm I read the stage card, the notes, and the appendix row, and the script's verified claim matches the paper's stated claim.

## Value Reconciliation (pass-2 augmentation)

Enumerating every RESULT/deliverable value the scripts emit (from `.py`/`.wl` source + committed outputs), excluding pass/fail flags, residual-zero check values, and pure scaffolding.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `e1/e0 = [m(q-r)+r·xi]/[delta+xi-m·q(q-r)]` | py:63 / wl:51; sympy out 9-23 | tex:16-18 (eq evec-rank2); md:47-49 | MATCH |
| `D_{q,r} = [delta+xi-mq(q-r)]²+[m(q-r)+r·xi]²` | py:77 / wl:61 | tex:22-24 (eq Dqr); md:74-75 | MATCH |
| `(z·e_-)²/z0² = [delta+(1+qr)xi]²/D_{q,r}` | py:81 / wl:65; sympy out 30-41 | tex:29-31 (eq z-overlap); md:64-66 | MATCH |
| `(s·e_-)²/s0² = [delta+(1+rt)xi-m(q-r)(q-t)]²/D_{q,r}` | py:82 / wl:66-69; sympy out 42-55 | tex:35-37 (eq s-overlap); md:68-71 | MATCH |
| `F_(q,r,t) = [δ+(1+qr)ξ]²[δ+(1+rt)ξ-m(q-r)(q-t)]² / [(1-ξ)D_{q,r}²]` | py:96-100 / wl:71-74; sympy out 58-73 | md:85-88 | MATCH (notes) |
| `F_track = [δ+(1+q²)ξ]²[δ+(1+qt)ξ]² / [(1-ξ)((δ+ξ)²+q²ξ²)²]` (r=q collapse) | py:106-109 / wl:87-91; sympy out 79-86 | md:108-110 | MATCH (notes) |
| `D_src = [δ+ξ-mλ0 R_U(R_U-1)]²+λ0[ξ+m(R_U-1)]²` | py:120 / wl:101-104 | md:152-154 | MATCH (notes) |
| `F_src` (a1²b1²/((1-ξ)D_src²), a1=δ+(1+λ0R_U)ξ, b1=δ+(1+λ0)ξ-mλ0(R_U-1)²) | py:117-121 / wl:99-105; sympy out 92-107 / math out 40 | md:145-148 (+143-154) | MATCH (notes) |
| `F_flat = [δ+(1+λ0)ξ]⁴/[(1-ξ)((δ+ξ)²+λ0ξ²)²]` | py:134-137 / wl:116-119 | md:170-172 | MATCH (notes) |
| `H_n^(src) = -2λ0 m ξ/[δ+(1+λ0)ξ]` (= -4mξ/(9δ+11ξ)) | py:149 / wl:131; sympy out 118-121 / math out 53 | md:186-187 | MATCH (notes) |
| `H_F^(src) = 2λ0[ξ((δ+ξ)²+λ0ξ²)+2mδ(δ+(1+λ0)ξ)]/[(δ+(1+λ0)ξ)((δ+ξ)²+λ0ξ²)]` | py:156-162 / wl:135-139; sympy out 123-128 | md:195-199 | MATCH (notes) |
| `lambda_0 = 2/9` (specialization constant) | py:48 / wl:40 | md:133-134 (`t²=λ0=2/9`) | MATCH (notes; upstream input) |
| `n_src = G_flat - m + eps·H_n^(src) + O(eps²)` (linear law) | py:169 / wl:141 | md:182-187 | MATCH (notes) |
| `F_src/F_flat = 1 + eps·H_F^(src) + O(eps²)` (linear law) | py:170 / wl:142 | md:191-199 | MATCH (notes) |

INTERNAL (scaffolding, no finding): `n_req`, `ratio_row1`, `ratio_row2`, `Z_overlap`/`S_overlap` raw build forms (vs the boxed `*_expected`), `q_tied_qr`/`a1`/`b1`, `F_src_direct`, `R_series`, banner/sub-banner strings.

All 14 deliverable values reconcile: the three card-boxed forms appear verbatim in `stage_042.tex` and the seven notes-only forms (F, collapse, F_src, flat-U, deformations, series, λ0) appear in the `.md` notes (the terse card legitimately omits them). No MISMATCH, no MISSING-DELIVERABLE.

`reconciliation: complete; 14 values checked, 0 misaligned`

## Self-test notes

Checked: (1) variable independence — the `sp.diff(n_src, R_U)` and `D[fRatio, rU]` derivatives in §25.6 are non-trivial because `n_src`/`fRatio` genuinely depend on `R_U` (numerator `lam0·R_U²`, denominator `(R_U-1)²`), so this is NOT the zero-derivative trap and the assertions can fail. (2) No integrals over unbounded domains, so the parity/symmetry trap does not apply. (3) Trivial-case: at `r=q` the §25.3 collapse residual reduces to 0 against the Stage-040 form (confirmed in both outputs); at `R_U=1` the §25.5 residual is 0 (a1=b1, D_src→(δ+ξ)²+λ0ξ²). (4) No missing-script directive (both engines present). (5) No fix prescribed that would introduce a new paper_misalignment — the only finding is informational stale_output, deferred to the dedicated numbering pass. The single finding is non-directive; no Codex edits proposed.
