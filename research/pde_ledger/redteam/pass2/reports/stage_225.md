---
unit_id: 225
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 225 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_225.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (row 62 + Part VII narrative paras 684–745 reviewed)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` (stage_225.tex:15): "Microscopic compiler $\Xi_1=N_{01}/N_0-D_{01}/D_0$ together with first-order even-preserving compensation conditions and the first surviving mixed-sector family." The card's purpose is the first microscopic formula for the weak-axisymmetric prefactor slope $\Xi_1$; the derivation ledger expands $u_2=-D_2/D_0$, $u_4=(D_2^2-D_0D_4)/D_0^2$, $P_0=N_0/D_0$ to first weak-axisymmetric order and solves the first conservative compensation surface, sorting wall-only / BdG-only / mixed-sector mechanisms. The notes (`.md` §§2–6, §8) enumerate the full deliverable set: (1) arbitrary-base first-order formulas for $u_2^{(1)},u_4^{(1)},\Xi_1$; (2) the exact conservative compensation surface $D_{21}=-u_2 D_{01}$, $D_{41}=(D_4/D_0)D_{01}=(u_2^2-u_4)D_{01}$, reducing on a one-pole branch ($u_4=4u_2^2$) to $D_{41}=-3u_2^2 D_{01}$; (3) the primitive logarithmic-slope compiler (B/Δ/S2/H/Q/P-raw drifts → $D_{01},D_{21},D_{41},N_{01}$); (4) the Stage-223 compatibility-point numbers ($D_0,D_2,D_4,u_2,u_4,D_4/D_0,P_{0,\text{target}}$, one-pole identity $u_4-4u_2^2=0$); (5) the concrete $\Xi_1$ linear form (9 coeffs); (6) wall-only and pure-BdG no-go results ($\Delta_{\rm BdG}\approx-5.119\times10^{-5}$); (7) the mixed/U rank-2/nullity-3 compensation matrix, its null basis $v_1,v_2,v_3$, and $\Xi_1(v_i)$; (8) the four transported Stage-224 amplitude windows on $v_1$. The appendix row 62 summarizes the same. **Note:** the card's `\stagefield{Verification}` line says "Mathematica audit: none yet," and notes §8 lists only the SymPy file — both predate the now-present `.wl`. This is a documentation lag, not a value defect (see Verdict justification; not raised as a blocking finding).

## What the script claims to verify

Both scripts verify the full deliverable list above. SymPy (docstring/banner line 14) derives the first-order slopes by `sp.diff(...).subs(eps,0)`, solves the compensation surface via `sp.solve(u2_1==0,D21)` then `sp.solve(u4_1==0,D41)`, reduces the one-pole case by substituting $D_4=-3D_0u_2^2$, builds the primitive compiler by differentiating dressed (exponentially parameterized) primitives, pins the seven Stage-223 numbers via `assert_close`, asserts the 9 $\Xi_1$ coefficients, proves wall-only triviality via `sp.solve`, proves the BdG no-go via the symbolic determinant + sample value, constructs the mixed nullspace by block elimination (`A11.inv()`), and divides the four budgets by $\sigma_1$. Mathematica mirrors each deliverable with `expectZero/expectClose/expectTrue/expectVectorClose`, plus two extra controls: an `expectNonZero` negative control on the one-pole reduction (M2) and `NullSpace[wallMatrix]===\{\}` plus an `expectNonZero` on the BdG sample determinant (M6).

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| Output $\Xi_1=N_{01}/N_0-D_{01}/D_0$ | py:51 `Xi1-(N01/N0-D01/D0)==0`; wl:151 M1 Xi1 | match |
| $u_2^{(1)}$ formula | py:41-42; wl:139-142 M1 | match |
| $u_4^{(1)}$ formula | py:43-50; wl:143-150 M1 | match |
| Compensation surface $D_{21}=-u_2D_{01}$, $D_{41}=(D_4/D_0)D_{01}$ | py:64-66; wl:158-160 M2 | match |
| One-pole reduction $D_{41}=-3u_2^2D_{01}$ | py:70-72; wl:162-168 M2 (+ neg. control) | match |
| Primitive log-slope compiler | py:161-193; wl:267-288 M3 (16 checks) | match |
| Stage-223 numbers ($D_0,D_2,D_4,u_2,u_4,D_4/D_0,P_0$) | py:233-241; wl:323-330 M4 | match |
| 9 $\Xi_1$ coefficients at sample | py:273-281; wl:340-348 M5 | match |
| Wall-only no-go | py:297-300; wl:359-361 M6 | match |
| Pure-BdG no-go ($\Delta_{\rm BdG}$ formula + sample) | py:313-320; wl:368-383 M6 | match |
| Mixed/U matrix, rank 2 / nullity 3, basis, $\Xi_1(v_i)$ | py:329-383; wl:387-442 M7 | match |
| 4 transported amplitude windows | py:404-409; wl:444-461 M8 | match |
| (card Verification line "none yet") | `.wl` now present | extra (doc lag, non-blocking) |

`paper_alignment: aligned` — every deliverable has a faithful, non-tautological script-side check on both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 41-43 | `simplify(u2_1 - …)==0` (and `+(D21+u2 D01)/D0`) | $u_2^{(1)}$ | yes |
| A2 | sympy | 43-50 | `simplify(u4_1 - …)==0` | $u_4^{(1)}$ | yes |
| A3 | sympy | 51 | `simplify(Xi1-(N01/N0-D01/D0))==0` | Output $\Xi_1$ | yes |
| A4 | sympy | 64-66 | compensation surface residuals | comp. surface | yes |
| A5 | sympy | 71-72 | `one_pole_D41-(-3u2^2)D01==0` | one-pole reduction | yes |
| A6 | sympy | 161-193 | 16 drift/compiler residuals | primitive compiler | yes |
| A7 | sympy | 233-241 | `assert_close` × 7 + one-pole identity | Stage-223 numbers | yes |
| A8 | sympy | 273-281 | `assert_close` × 9 coeffs | $\Xi_1$ linear form | yes |
| A9 | sympy | 297-300 | solve → `{xK:0,xM:0}` | wall no-go | yes |
| A10 | sympy | 315-320 | det formula + sample value | BdG no-go | yes |
| A11 | sympy | 351-383 | rank==2, basis, $\Xi_1(v_i)$ | mixed survivor | yes |
| A12 | sympy | 406-409 | 4 window `assert_close` | transported windows | yes |
| B1 | math | 139-151 | M1 expectZero × 3 | $u_2^{(1)},u_4^{(1)},\Xi_1$ | yes |
| B2 | math | 158-164 | M2 expectZero × 4 | comp. surface + one-pole | yes |
| B3 | math | 165-168 | M2 expectNonZero (−2 control) | one-pole load-bearing guard | yes |
| B4 | math | 267-288 | M3 expectZero × 16 | primitive compiler | yes |
| B5 | math | 323-330 | M4 expectClose × 7 + expectZero | Stage-223 numbers | yes |
| B6 | math | 340-348 | M5 × 9 coeffs | $\Xi_1$ linear form | yes |
| B7 | math | 360-383 | M6 rank/nullspace/det + nonZero | wall + BdG no-go | yes |
| B8 | math | 397-442 | M7 rank/nullity/in-nullspace/vec/$\Xi_1$ | mixed survivor | yes |
| B9 | math | 454-461 | M8 expectClose × 4 | transported windows | yes |

## Findings

None. (Zero findings.)

## Independent-derivation check (Mathematica)

The `.wl` is **genuinely independent**, not a transliteration. Three distinct routes confirm it:

1. **First-order slopes.** SymPy uses analytic differentiation: `u2_1 = sp.diff(u2A, eps).subs(eps, 0)/lam` (py:36-38, 143-149). Mathematica uses Taylor-series coefficient extraction: `epsSlope[expr_] := Coefficient[Normal[Series[expr,{eps,0,1}]], eps]` (wl:92, applied wl:134-136, 227-238). Different algorithms (symbolic `D` vs. truncated series) producing the same first-order forms.

2. **Surviving-family nullspace.** SymPy hand-rolls block elimination: partitions the 2×5 matrix into `A11=A[:,:2]`, `A12=A[:,2:]`, inverts `A11.inv()`, and builds each basis vector as `-A11_inv*A12*free` (py:354-363). Mathematica calls the built-in `NullSpace[mixedMatrix]` for a raw basis (wl:398), then re-gauges it into the same free-variable normalization with `LinearSolve[Transpose[freeBlock], UnitVector[3,j]]` (wl:401-408) — and additionally asserts `mixedMatrix . v == 0` for each normalized vector (wl:410-416), a residual check the SymPy script does not perform.

3. **Wall-only no-go.** SymPy solves `sp.solve([...],[xK,xM])` and asserts the lone solution is `{0,0}` (py:297-300). Mathematica instead proves it structurally: `MatrixRank[wallMatrix]==2` AND `NullSpace[wallMatrix]===\{\}` (wl:360-361). Different proof strategies for the same conclusion.

The `.wl` also carries genuine extra controls absent from the `.py`: the M2 `expectNonZero["…wrong -2 coefficient…"]` negative control and the M6 `expectNonZero` BdG-determinant nonvanishing check. The variable names (`Kwall`/`mass`/`omU` vs. `K`/`M`/`OmU`) and the assertion *order within M7* also differ from the `.py`. **Independence call: independent.**

## Engine cross-check

Both engines agree to ≥14 significant figures on every numeric deliverable and produce identical symbolic residuals (all `0`). Spot comparisons (sympy.txt ↔ mathematica.txt):
- $D_0$: 24.23730998862225 ↔ 24.23730998862225050 (residual 4.9e-14, tol 1e-12) — agree.
- $\Xi_1$ coeff $x_K$: −1.0097554097702972 ↔ −1.00975540977029715 — agree.
- $\Delta_{\rm BdG}$: −5.118869961200109e-05 ↔ −0.00005118869961200108 — agree.
- $\Xi_1(v_2)$: −14.431027813975458 ↔ −14.43102781397545563 — agree.
- both_30 window: 2.1678897807090363 ↔ 2.16788978070903489 — agree.
All M1–M8 `expectZero` residuals are exact `0`; all `expectClose` residuals are ≤ ~5e-14 (under their 1e-12 tolerances). No `engine_disagreement`.

## Verdict justification

**clean.** I read the paper card, the full notes file, and the relevant Part VII appendix rows/narrative before opening the scripts, and the scripts' claims match the paper's claims deliverable-for-deliverable. Attacks tried that failed: (a) I checked the one-pole reduction for tautology — it is genuinely load-bearing on both engines (see below); (b) I checked the `.wl` for transliteration across three subsystems and found independent routes each time; (c) I checked variable-independence in every `diff`/`Series` slope (each dressed quantity genuinely depends on the relevant `eps`-coupled log-slope, so no slope collapses to identically zero); (d) I checked the sign convention $D_2,D_4<0$ at the sample is consistent between notes (−1.1856, −0.17399) and both outputs; (e) I checked the wall-only and BdG no-go are non-trivial (the BdG determinant is a nonzero closed form, output wl:119); (f) I confirmed outputs are fresher than scripts (17:05 vs 16:46/16:51 — no `stale_output`); (g) I confirmed the renumber landed — no surviving 240/241/242 stage labels in the notes (the lone "0.241…/0.242…" hits are matrix-coefficient digits, not labels). The only paper-side imperfection is the documentation lag on the card's `\stagefield{Verification}` line ("Mathematica audit: none yet") and notes §8 (SymPy-only), both of which predate the now-present `.wl`; this understates verification coverage rather than misstating any result value, touches no `\stagefield{Output}` quantity, and is squarely the deferred prose-pointer-sync class (the orchestrator's post-batch tracker/doc-pointer sync), not a new red-team blocker. I am not raising it as a blocking finding.

### Re-audit confirmations (per re-audit prompt)
- **(a) `.wl` genuinely independent:** YES — independent (series-coefficient slopes vs. analytic-diff; built-in `NullSpace`+re-gauge vs. hand block-elimination; rank+empty-nullspace vs. `solve`). Not a port.
- **(b) one-pole 0==0 → $D_4=-3D_0u_2^2$ fix holds:** YES. The M2 check `onePoleD41 - (-3 u2Base^2)D01` is a substantive identity (passes, residual 0), and the negative control `expectNonZero["M2 wrong -2 coefficient…", onePoleD41-(-2u2^2)D01]` genuinely *can* fail and *does* return a nonzero residual `-((D01*D2^2)/D0^2)` (output line 27) — equal to $-u_2^2 D_{01}$, confirming the −3 coefficient is uniquely pinned and the check is not 0==0. The SymPy side mirrors this with the substantive `assert one_pole_D41-(-3*u2**2)*D01==0` (py:72). Fix holds on both engines.
- **(c) every emitted deliverable value reconciles to `.tex`+`.md`:** YES — see Value Reconciliation below.

## Value Reconciliation (pass-2 augmentation)

The `.tex` card is terse and carries only the symbolic Output formula plus the structural deliverable names; the per-stage `.md` notes are the natural numeric carrier and report every emitted number. Both `.txt` outputs are fresh (newer than their scripts), so reconciliation is grounded in source + committed outputs.

| value | source (py/wl + output) | .tex/.md location | status |
|---|---|---|---|
| $\Xi_1=N_{01}/N_0-D_{01}/D_0$ | sympy.txt:6; wl M1 | tex:15 Output; md:220 | MATCH |
| $u_2^{(1)}=(-D_0D_{21}+D_2D_{01})/D_0^2$ | sympy.txt:4; wl M1 | md:206-210 | MATCH |
| $u_4^{(1)}$ closed form | sympy.txt:5; wl M1 | md:213-216 | MATCH |
| $D_{21}=-u_2D_{01}$ (=$D_2D_{01}/D_0$) | sympy.txt:9; wl M2 | md:241 | MATCH |
| $D_{41}=(D_4/D_0)D_{01}$ | sympy.txt:10; wl M2 | md:244 | MATCH |
| one-pole $D_{41}=-3u_2^2D_{01}$ | sympy.txt:11; wl M2 | md:260 | MATCH |
| $K_{\rm compat}=24.4737548792910$ | sympy.txt:16; wl M4 ctx | md:123 | MATCH |
| $D_0=24.2373099886223$ | sympy.txt:17; wl M4 | md:130 | MATCH |
| $D_2=-1.18562046858190$ | sympy.txt:18; wl M4 | md:131 | MATCH |
| $D_4=-0.173991572849491$ | sympy.txt:19; wl M4 | md:132 | MATCH |
| $u_2=0.0489171640391802$ | sympy.txt:20; wl M4 | md:135 | MATCH |
| $u_4=0.00957155575054425$ | sympy.txt:21; wl M4 | md:136 | MATCH |
| $D_4/D_0=-0.00717866681290820$ | wl M4 D4/D0 | md:138 | MATCH |
| $P_{0,\rm target}=0.00206979231806289$ | sympy.txt:22; wl M4 | md:149 | MATCH |
| $\Xi_1$ coeff $x_K=-1.00975540977030$ | sympy.txt:25; wl M5 | md:362 | MATCH |
| $\Xi_1$ coeff $x_M=0$ | sympy.txt:26; wl M5 | md:373 | MATCH |
| coeff $x_{\lambda_B}=0.00418038073077834$ | sympy.txt:27; wl M5 | md:363 | MATCH |
| coeff $x_\varpi=-0.00418038073077834$ | sympy.txt:28; wl M5 | md:364 | MATCH |
| coeff $x_{\lambda_U}=0.324464020216766$ | sympy.txt:29; wl M5 | md:366 | MATCH |
| coeff $x_{\lambda_W}=1.69086641859305$ | sympy.txt:30; wl M5 | md:366 | MATCH |
| coeff $x_{\lambda_R}=0.423379354382463$ | sympy.txt:31; wl M5 | md:367 | MATCH |
| coeff $x_{\Omega_U}=-0.747843374599229$ | sympy.txt:32; wl M5 | md:368 | MATCH |
| coeff $x_{\Omega_W}=-4.11424577297551$ | sympy.txt:33; wl M5 | md:369 | MATCH |
| $\Delta_{\rm BdG}=-5.11886996120011\times10^{-5}$ | sympy.txt:40; wl M6 | md:446 | MATCH |
| mixed matrix row1 (5 entries) | sympy.txt:43; wl M7 | md:464 | MATCH |
| mixed matrix row2 (5 entries) | sympy.txt:43; wl M7 | md:465 | MATCH |
| null basis $v_1,v_2,v_3$ | sympy.txt:45-47; wl M7 | md:472-478 | MATCH |
| $\Xi_1(v_1)=1.36026097049402$ | sympy.txt:49; wl M7 | md:484 | MATCH |
| $\Xi_1(v_2)=-14.4310278139755$ | sympy.txt:50; wl M7 | md:485 | MATCH |
| $\Xi_1(v_3)=-5.01037421295998$ | sympy.txt:51; wl M7 | md:486 | MATCH |
| Stage-224 $\Xi_1$ budgets (4) | py:399-403; wl M8 in | md:155-168 | MATCH |
| window both_10 = 0.270485102839510 | sympy.txt:54; wl M8 | md:520 | MATCH |
| window one_10 = 0.542262903708006 | sympy.txt:55; wl M8 | md:525 | MATCH |
| window both_30 = 2.16788978070904 | sympy.txt:56; wl M8 | md:530 | MATCH |
| window one_30 = 3.40747461278373 | sympy.txt:57; wl M8 | md:535 | MATCH |

INTERNAL scaffolding (no prose expected, no finding): `assert_close`/`expectClose` tolerances; the −2-coefficient negative-control residual; `expectZero` residual values (all 0); rank/nullity boolean flags; `eps`/`lam` symbol declarations; the dressed-primitive intermediate slopes ($B_{i,1},Z_{i,1},\Delta_1,S_{2,1},H_1,Q_1,P_1^{\rm raw},N_{0,1}$ — these appear as compiler intermediates and ARE additionally documented in md §4.1–4.2, but are not bottom-line deliverables); pass/fail print labels.

reconciliation: complete; 37 deliverable values checked, 0 misaligned.

## Self-test notes

Checked all five traps. (1) Variable independence: every `sp.diff(...,eps)` / `Series[...,{eps,0,1}]` acts on a dressed quantity that genuinely contains `eps` via at least one `exp(eps·x)` factor, so no slope is identically zero (e.g. `B0d` depends on `xLB,xV`; `N0d` on the full mixed set) — confirmed against the nonzero coefficients in output M3/M5. (2) Parity/integrals: no unbounded-domain integrals in this stage (purely algebraic), N/A. (3) Trivial-case: the negative control `(-3)−(-2)=-1` times `u2^2 D01` = `-(D2^2/D0^2)D01` ≠ 0, matching output line 27 — the load-bearing one-pole check cannot be 0==0. (4) Path specs: N/A (no missing-script finding). (5) Paper round-trip: no fix prescribed; the only paper imperfection is a non-blocking documentation lag on the card Verification line, routed to ordinary doc-pointer sync, not a finding. Zero findings; no directive written.
