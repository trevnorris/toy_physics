---
unit_id: 065
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
  notes_stage_files: [moving_throat_pde_stage065_thin_wall_confinement.md]
  paper_appendix: present
---

# Audit unit 065 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_065.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage065_thin_wall_confinement.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 108)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.txt` (present; STALE)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.txt` (present; STALE)

## What the paper claims

Stage 065 translates the parent gain thresholds (from Stage 064) into wall-amplitude thresholds for a concrete thin-wall confinement profile `V_conf(r;a) = V0 f((r-a)/ell)`. The card's bottom-line `\stagefield{Output}` is terse: "A physical wall-amplitude interpretation of `g_phi`," with the single boxed derivation equation `g_phi = V0/ell` (eq:app-stage065-gphi). The notes (the authoritative carrier here) enumerate the full deliverable chain: (1) `g_phi = V0/ell`; (2) exact shell integral `I1 = 4 pi ell (a^2 J1 + 2 a ell J2 + ell^2 J3)` with `Jn = int xi^n f'^2/H`; (3) `J2 = 0` for a centered symmetric layer; (4) exact gain `G_eq = 4 pi V0^2 (a^2 J1/ell + 2 a J2 + ell J3)/K_X`; (5) thin-wall leading gain `G_eq^tw = 4 pi a^2 V0^2 J1/(K_X ell)`; (6) wall-amplitude thresholds `V0_fail^2 = K_X ell G_fail/(4 pi a^2 J1)`, similarly suff; (7) after inserting `kappa = K_X L^2/T_X`, K_X cancels, leaving `V0_fail^2 = T_X ell Pe_req/(4 pi a^2 L^2 J1 Delta_inf)` (suff with Delta_0); (8) under constant compressibility `J1 = I_f/H_w`, giving `V0_fail^2 = H_w T_X ell Pe_req/(4 pi a^2 L^2 I_f Delta_inf)`. The appendix row (108) summarizes: "Wall amplitude thresholds for `V_conf = V0 f((r-a)/ell)`." Status: ExactClosure/Reduced.

## What the script claims to verify

Both scripts assert: (a) `g_phi = V0/ell` (assigned; the SymPy 1/ell scaling is independently derived from `V_conf` via `sp.diff` at one e-fold, lines 108-123); (b) the `(1,2,1)` polynomial coefficient structure of `I1` via a direct shell integral against the concrete Gaussian profile `f(u)=exp(-u^2)` with constant `h'=1`; (c) `J2 = 0` by parity for that symmetric profile; (d) the thin-wall remainder `G_eq^sym - G_eq^tw = 4 pi V0^2 ell J3/K_X` (i.e. the dropped `O(ell/a)` term); (e) the threshold solves for `V0^2` from `G_eq^tw = G_fail/G_suff`, the `K_X` cancellation under `kappa = K_X L^2/T_X`, and the constant-compressibility reduction `J1 -> I_f/H_w`. The Gaussian profile is an audit-only anchor (the moment integrals `J1_num`, `J3_num` are scaffolding, not stage deliverables). All eight notes deliverables have a corresponding script-side symbolic check.

## Paper ↔ script cross-check

| Paper/notes deliverable | Script check | Status |
|---|---|---|
| (1) `g_phi = V0/ell` | py:61 assign + py:114-123 derived 1/ell scaling; wl:66 | match |
| (2) `I1 = 4 pi ell (a^2 J1 + 2 a ell J2 + ell^2 J3)` | py:66 + py:128-133 direct shell integral; wl:46-50 | match |
| (3) `J2 = 0` (symmetric) | py:70,100; wl:43 | match |
| (4) `G_eq = 4 pi V0^2 (a^2 J1/ell + 2 a J2 + ell J3)/K_X` | py:74,78; wl:74,78 | match |
| (5) `G_eq^tw = 4 pi a^2 V0^2 J1/(K_X ell)` | py:76,81-84; wl:76,81-84 | match |
| (6) `V0_fail^2 = K_X ell G_fail/(4 pi a^2 J1)` | py:140-144; wl:92-96 | match |
| (7) K_X cancellation -> `T_X ell Pe_req/(4 pi a^2 L^2 J1 Delta_*)` | py:147-160; wl:98-112 | match |
| (8) const-H `V0_fail^2 = H_w T_X ell Pe_req/(4 pi a^2 L^2 I_f Delta_inf)` | py:165-176; wl:116-128 | match |

`paper_alignment: aligned` — every notes deliverable maps to a non-tautological script check, and the symbolic forms agree verbatim with the notes. The `.tex` card is terse but consistent (its only equation `g_phi=V0/ell` matches).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 81-84 | `expect_zero(G_eq_sym - G_eq_tw - 4πV0²ell·J3/KX)` | claim 5 (thin-wall remainder) | yes |
| A2 | sympy | 100 | `expect_zero(J2_num)` (Gaussian) | claim 3 | yes |
| A3 | sympy | 105-106 | `expect_zero(rel correction - (ell/a)²J3/J1)` | claim 5 scaling | yes |
| A4 | sympy | 122-123 | `expect_zero(ampl + 2V0·e^-1/ell)` | claim 1 (1/ell scaling) | yes |
| A5 | sympy | 132-133 | `expect_zero(I1_full - I1_poly)` (1,2,1) | claim 2 | yes |
| A6 | sympy | 153-160 | `expect_zero(K_X cancellation fail/suff)` | claim 7 | yes |
| A7 | sympy | 169-176 | `expect_zero(const-H fail/suff)` | claim 8 | yes |
| A8 | sympy | 182-183 | `expect_zero(J1_num - If_num/Hw_num)` | claim 8 reduction | yes |
| A9 | math | 43 | `expectZero[j2Num]` | claim 3 | yes |
| A10 | math | 49-50 | `expectZero[i1Direct - i1Poly]` | claim 2 | yes |
| A11 | math | 55-56 | `expectZero[j1Num - ifMomDirect/hwSym]` | claim 8 | yes |
| A12 | math | 81-84 | `expectZero[gEqSym - gEqTw - 4πv0²ell·j3/kx]` | claim 5 | yes |
| A13 | math | 105-112 | `expectZero[K_X cancellation fail/suff]` | claim 7 | yes |
| A14 | math | 121-128 | `expectZero[const-H fail/suff]` | claim 8 | yes |

No tautological rows: each `expect_zero`/`expectZero` compares an independently-built quantity (a symbolic `solve`/division result, a direct `integrate`, or a `diff`) against the target closed form. A1/A12 could in principle fail had the `O(ell/a)` term been mis-derived; A4 genuinely exercises the `1/ell` power via `sp.diff` on `V_conf` (a wrong power would surface). A5/A10 exercise the `(1,2,1)` coefficient structure via direct shell integration.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.txt` (mtime 2026-05-22 19:52:15) vs script `.py` (mtime 2026-06-03 15:59:11)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.txt` (mtime 2026-05-22 19:52:21) vs script `.wl` (mtime 2026-06-03 15:59:11)

**What's wrong:**
Both committed transcripts predate their scripts by ~12 days. The content disagreement is the stage self-label: the SymPy transcript header reads `STAGE 48 — THIN-WALL CONFINEMENT BRANCH` (line 3) and `STAGE 48 AUDIT PASSED` (line 48), whereas the current `.py` banner emits `STAGE 65` (py:51) and `STAGE 65 AUDIT PASSED` (py:185). The Mathematica transcript header reads `STAGE 048 — THIN-WALL CONFINEMENT BRANCH` (line 16) whereas the current `.wl` emits `STAGE 065` (wl:58). The committed SymPy docstring still carries the stale label too (`Stage 48 SymPy audit`, py:2; `Stage-44 thresholds`, py:22), but that is docstring prose, not transcript content.

The numeric/symbolic math content of both transcripts otherwise matches what the current scripts produce (g_phi, I1, G_eq, all `= 0` residuals, J1_num = sqrt(2)*sqrt(pi)/2 = Sqrt[Pi/2], J3_num = 3*sqrt(2)*sqrt(pi)/8). Only the stage-number self-labels are stale.

**Why this matters:**
The committed transcript is the auditable record. A reader cross-referencing the Stage-065 card against a transcript that announces "STAGE 48 / STAGE 048" gets a numbering contradiction. This is the known SCRIPT/OUTPUT-band stale-label pattern.

**Required change:**
Refresh both transcripts by re-running the current scripts (`python3 .../stage065_..._sympy_audit.py` and `math -script .../stage065_..._mathematica_audit.wl`) and committing the regenerated `.txt`. No script-logic edit is required for this finding; the banner code already emits the correct `STAGE 65 / STAGE 065`.

**Verification:**
After refresh, SymPy transcript line 3 reads `STAGE 65 — THIN-WALL CONFINEMENT BRANCH`, line for the footer reads `STAGE 65 AUDIT PASSED`; Mathematica transcript header reads `STAGE 065 — THIN-WALL CONFINEMENT BRANCH`. Output mtimes newer than script mtimes.

### F2 — numbering self-label (stale docstring/comment labels)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py:2` — `Moving-Throat PDE — Stage 48 SymPy audit`
- `.../stage065_..._sympy_audit.py:22` — `5. Inserting the Stage-44 thresholds and kappa = K_X L^2/T_X gives`

**What's wrong:**
The SymPy docstring carries stale stage self-labels: the title says `Stage 48` (should be `Stage 65`) and the threshold-insertion step cites `Stage-44 thresholds`. The notes/appendix make clear the fail/succeed surfaces being inserted come from **Stage 061** (notes line 47: "comparing this with the Stage-061 fail/succeed surfaces"; notes section 5 header "the Stage-061 operator formulas"). `Stage 48` is the pre-renumber (+17) self-label for this unit; `Stage-44` is the pre-renumber self-label for the Stage-061 source. These are unambiguous self/cross labels for THIS unit and its named upstream source, so they fall inside the in-loop Reading-2 scope (verdict:findings stage ⇒ fix its unambiguous self-labels).

The `.wl` script header self-label is already correct (`STAGE 065`, wl:58), so no `.wl` label edit is needed.

**Why this matters:**
Same numbering-contradiction concern as F1, but at the source level rather than the transcript. The docstring is the human-readable statement of which stage the script audits; "Stage 48" contradicts the filename `stage065` and the card.

**Required change:**
- py:2 — change `Stage 48 SymPy audit` to `Stage 65 SymPy audit`.
- py:22 — change `Inserting the Stage-44 thresholds` to `Inserting the Stage-061 thresholds` (the source of `G_fail`/`G_suff` per the notes).
No assertion or math changes. Re-run after the edit so the refreshed transcript (F1) also picks up the corrected docstring is not echoed (docstring is not printed, so transcript is unaffected, but keep edits in one pass).

**Verification:**
py:2 reads `Stage 65 SymPy audit`; py:22 reads `Inserting the Stage-061 thresholds`. Script still exits 0 with all checks `= 0`.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent** re-derivation, not a transliteration. Evidence:
- Block ordering differs: the `.wl` front-loads an "INDEPENDENT SHELL-INTEGRAL DERIVATION" section (wl:26-56) computing `j1Num/j2Num/j3Num` and `i1Direct` via `Integrate` BEFORE the symbolic stage block; the `.py` interleaves the concrete moments in the middle (py:89-133).
- Threshold method differs: `.py` uses `sp.solve(sp.Eq(G_eq_tw, G_fail), V0**2)` (py:140-141); the `.wl` instead forms `gEqCoeff = gEqTw/v0^2` and divides `v0FailSq = gFail/gEqCoeff` (wl:90-93) — an algebraic coefficient division, a structurally different route to the same `V0^2`.
- The `.wl` does NOT carry the `.py`'s `g_phi` 1/ell `sp.diff` derivation (py:108-123) nor the relative-correction `(ell/a)^2 J3/J1` check (py:105-106); each engine chose its own anchor set.
- The Gaussian moment integrals are evaluated by each engine's own integrator (SymPy `sp.integrate` vs Mathematica `Integrate`), and the two return the same value in different surface forms (`sqrt(2)*sqrt(pi)/2` vs `Sqrt[Pi/2]`).

## Engine cross-check

Both engines agree. Final symbolic forms (modulo each CAS's print formatting):
- `g_phi = V0/ell` (py txt:5) ≡ `g_phi = v0/ell` (wl txt:18).
- `I1 = 4*pi*ell*(J1*a^2 + 2*J2*a*ell + J3*ell^2)` (py txt:6) ≡ `I_1 = 4*ell*(a^2*j1 + 2*a*ell*j2 + ell^2*j3)*Pi` (wl txt:19).
- `G_eq^tw = 4*pi*J1*V0^2*a^2/(K_X*ell)` (py txt:18) ≡ `(4*a^2*j1*Pi*v0^2)/(ell*kx)` (wl txt:23).
- `V0_fail^2 (kappa) = Pe_req*T_X*ell/(4*pi*Delta_inf*J1*L^2*a^2)` (py txt:33) ≡ `(ell*peReq*tx)/(4*a^2*deltaInf*j1*len^2*Pi)` (wl txt:32).
- All `expect_zero`/`expectZero` residuals are `0` / `PASS` in both transcripts.
`J1_num`: `sqrt(2)*sqrt(pi)/2` (py) = `Sqrt[Pi/2]` (wl) — equal. `J3_num`: `3*sqrt(2)*sqrt(pi)/8` (py) = `(3*Sqrt[Pi/2])/4` (wl) — equal. No engine disagreement.

## Verdict justification

`findings: 2`, both low-severity numbering/freshness items; no math defect. I attacked: (a) the `g_phi = V0/ell` derivation — the comment at py:116-119 wrongly states `f'(0)=1` for the Gaussian (whose `f'(0)=0`), but the script explicitly acknowledges this and pivots to a genuine `1/ell`-scaling check at one e-fold (`-2 V0 e^-1/ell`), which is non-tautological and correct, so not a finding; (b) the moment integrals (J1=sqrt(pi/2), J3=(3/4)sqrt(pi/2)) — recomputed by hand, both correct; (c) the threshold `solve` and `K_X` cancellation — verified algebraically; (d) the `(1,2,1)` coefficient check — genuinely exercises the shell-weight expansion via direct `integrate`; (e) symbol assumptions (`positive=True` on V0/ell/a/K_X/J1/J3; J2 set to 0 in the symmetric branch) — all physically justified, none invalidates a `simplify`. Every one of the eight notes deliverables maps to a non-tautological check in both engines, and the symbolic outputs match the notes verbatim. The only defects are the stale transcripts (F1) and the stale `Stage 48`/`Stage-44` docstring self-labels (F2), both inside the in-loop Reading-2 scope. `material_change: false`.

## Value Reconciliation (pass-2 augmentation)

All deliverable values are symbolic closed forms; they are carried by the **notes** (`.md`), which the terse `.tex` card legitimately summarizes (card states only the boxed `g_phi=V0/ell`). Per the augmentation guards, a value living correctly in the `.md` is a MATCH.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `g_phi = V0/ell` | py:61 / wl:66; py-txt:5, wl-txt:18 | `.tex:16` (eq:app-stage065-gphi); `.md:30,106` | MATCH |
| `I1 = 4πell(a²J1 + 2a·ell·J2 + ell²J3)` | py:66 / wl:67; py-txt:6, wl-txt:19 | `.md:33,133` | MATCH |
| `I1\|J2=0 = 4πell(a²J1 + ell²J3)` | py:70 / wl:68; py-txt:11, wl-txt:20 | `.md:40,149` | MATCH |
| `G_eq = 4πV0²(a²J1/ell + 2aJ2 + ell·J3)/K_X` | py:74 / wl:74; py-txt:16, wl-txt:21 | `.md:43,158` | MATCH |
| `G_eq^tw = 4πa²V0²J1/(K_X·ell)` | py:76 / wl:76; py-txt:18, wl-txt:23 | `.md:46,166` | MATCH |
| `V0_fail² = K_X·ell·G_fail/(4πa²J1)` (pre-kappa) | py:140 / wl:92; py-txt:31, wl-txt:30 | `.md:49,188` | MATCH |
| `V0_suff² = K_X·ell·G_suff/(4πa²J1)` | py:141 / wl:93; py-txt:32, wl-txt:31 | `.md:51,190` | MATCH |
| `V0_fail²(kappa) = T_X·ell·Pe_req/(4πa²L²J1·Delta_inf)` | py:148 / wl:99; py-txt:33, wl-txt:32 | `.md:62,210` | MATCH |
| `V0_suff²(kappa) = T_X·ell·Pe_req/(4πa²L²J1·Delta_0)` | py:149 / wl:100; py-txt:34, wl-txt:33 | `.md:64,212` | MATCH |
| `V0_fail²\|H~const = H_w·T_X·ell·Pe_req/(4πa²L²I_f·Delta_inf)` | py:165 / wl:116; py-txt:41, wl-txt:42 | `.md:77,242` | MATCH |
| `V0_suff²\|H~const = H_w·T_X·ell·Pe_req/(4πa²L²I_f·Delta_0)` | py:166 / wl:117; py-txt:42, wl-txt:43 | `.md:79,244` | MATCH |

INTERNAL (audit-only Gaussian-anchor scaffolding, not stage deliverables — the profile `f=exp(-u²)` is chosen by the script to give concrete numbers, not a paper claim): `J1_num = sqrt(pi/2)`, `J2_num = 0`, `J3_num = (3/4)sqrt(pi/2)`, `If_num`, `Hw_num = 1`, `ampl = -2V0·e^-1/ell` (one-e-fold scaling anchor), and all residual/`PASS` flags.

reconciliation: complete; 11 deliverable values checked, 0 misaligned.

## Self-test notes

Variable-independence: the only `sp.diff` is `dV_dr = diff(V0·exp(-((r-a)/ell)²), r)` (py:115) — `EXPR` depends on `r`, so the derivative is non-trivial; the `D[fProf[u],u]` (wl:30) depends on `u` — both genuine, no identically-zero trap. Parity: J2 = ∫ xi·f'²/H over symmetric domain with `f'²` even ⇒ odd integrand ⇒ 0; the `expect_zero(J2_num)` correctly vanishes, and J1/J3 (even integrands) are correctly nonzero — confirmed by hand (J1=sqrt(pi/2), J3=(3/4)sqrt(pi/2)). Trivial-case: the threshold residuals reduce to 0 by the `solve`/cancellation algebra I re-derived. No paper round-trip issue: F1/F2 are label/freshness only and prescribe no math change.
