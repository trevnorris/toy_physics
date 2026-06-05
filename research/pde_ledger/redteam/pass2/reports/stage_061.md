---
unit_id: 061
batch: III.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage061_microscopic_gain_thresholds.md]
  paper_appendix: present
---

# Audit unit 061 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_061.tex`
- notes: `/var/projects/toy_projects/../notes/stages/moving_throat_pde_stage061_microscopic_gain_thresholds.md` → actual path `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage061_microscopic_gain_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 100)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage061_microscopic_gain_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage061_microscopic_gain_thresholds_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` is: *"The microscopic gain phase diagram \eqref{eq:app-stage061-G-thresholds}--\eqref{eq:app-stage061-phase-diagram}."* The stage rewrites the Stage-059 support/source verdict in terms of a microscopic gain. Its load-bearing deliverables are: (1) the factorization `\Xi_{\rm micro}=\kappa G_{\rm micro}` with `G_{\rm micro}:=\chi_\sigma\Lambda_\phi^2/K_X` and `\kappa=K_X L^2/T_X` (eq. app-stage061-Xi-micro); (2) the boxed gain thresholds `G_{\rm fail}=\Xi_{\rm fail}/\kappa`, `G_{\rm suff}=\Xi_{\rm suff}/\kappa` (i.e. `Pe_req/[\kappa\Delta_\infty]` and `Pe_req/[\kappa\Delta_0]`) and the boxed three-region phase diagram (eq. app-stage061-G-thresholds / app-stage061-phase-diagram). The notes add four further deliverables beyond the terse card: (3) the inverted threshold surfaces for `\chi_\sigma` and `\Lambda_\phi^2`; (4) the soft-support limit `\kappa\to0^+` giving `\Delta_0\to1/2`, `\Delta_\infty\to1`, hence `G_{\rm fail}\sim Pe_req/\kappa`, `G_{\rm suff}\sim 2Pe_req/\kappa`; (5) the compliant-mouth limit `\eta\to+\infty` giving `\Delta_0^{(\infty)}=(1-\mathrm{sech}\sqrt\kappa)/\kappa`, `\Delta_\infty^{(\infty)}=\tanh\sqrt\kappa/\sqrt\kappa`, hence `G_{\rm fail}^{(\infty)}=Pe_req/[\sqrt\kappa\tanh\sqrt\kappa]`, `G_{\rm suff}^{(\infty)}=Pe_req/[1-\mathrm{sech}\sqrt\kappa]`; (6) the combined stiff-support/compliant-mouth limit `\kappa\gg1` giving `\sqrt\kappa\,G_{\rm fail}^{(\infty)}\to Pe_req`, `G_{\rm suff}^{(\infty)}\to Pe_req`. The appendix row 100 summarizes it as the "Dimensionless gain phase diagram and geometry-controlled success/failure surfaces" with status `\StatusExactClosure{}`.

## What the script claims to verify

Both scripts pin the Stage-058 endpoint functions `Delta_0`, `Delta_inf` (parametrized by `kappa`, `eta`) as upstream inputs, then: (a) assert `Xi_micro - kappa*G_micro == 0` after substituting `T_X = K_X L^2/kappa` — the factorization (1); (b) build `G_fail = Pe_req/(kappa*Delta_inf)`, `G_suff = Pe_req/(kappa*Delta_0)` and the inverted surfaces `chi_fail/chi_suff/Lambda^2_fail/Lambda^2_suff`, asserting each inversion is consistent with `G_fail`/`G_suff` (deliverables 2,3); (c) take `kappa->0+` limits, asserting `Delta_0->1/2`, `Delta_inf->1`, and `kappa*G_fail->Pe_req`, `kappa*G_suff->2Pe_req` (4); (d) take `eta->oo` limits, asserting the `sech`/`tanh` endpoint forms and the closed-form `G_fail^(inf)`, `G_suff^(inf)` (5); (e) take the `sqrt(kappa)->oo` stiff-support limit, asserting `sqrt(kappa)*G_fail^(inf)->Pe_req`, `G_suff^(inf)->Pe_req` (6). The SymPy module docstring/theorem-ledger and the Mathematica banner enumerate exactly these five ledger claims.

## Paper ↔ script cross-check

| Paper/notes deliverable | Script-side check | Status |
|---|---|---|
| (1) `Xi_micro = kappa*G_micro`, `G_micro=chi*Lam^2/K_X` | py:40-41 / wl:48 `Xi_micro - kappa G_micro == 0` (subs `T_X=K_X L^2/kappa`) | match |
| (2) `G_fail=Pe_req/[kappa Delta_inf]`, `G_suff=Pe_req/[kappa Delta_0]` | py:44-49 / wl:39-47 build + print; exercised by surface-inversion asserts | match |
| (3) `chi`/`Lam^2` threshold surfaces | py:51-66 / wl:50-57, four `== 0` consistency asserts | match |
| (4) soft-support `kappa->0+`: `Delta_0->1/2`,`Delta_inf->1`, `G_fail~Pe_req/kappa`, `G_suff~2Pe_req/kappa` | py:69-79 / wl:59-66, four `== 0` asserts | match |
| (5) compliant-mouth `eta->oo`: sech/tanh forms + `G_fail^(inf)`,`G_suff^(inf)` | py:82-98 / wl:68-80, four `== 0` asserts | match |
| (6) stiff-support `kappa>>1`: `sqrt(kappa)G_fail->Pe_req`,`G_suff->Pe_req` | py:100-107 / wl:82-86, two `== 0` asserts | match |

Every paper-side deliverable has a non-tautological script-side check in both engines. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 40-41 | `expect_zero(Xi_micro.subs(TX,KX L^2/kappa) - kappa*G_def)` | claim 1 | yes |
| A2 | sympy | 59-62 | `expect_zero(chi_{fail,suff}.subs(...) - (KX/Lam^2)*G_{fail,suff})` | claim 3 | yes |
| A3 | sympy | 63-66 | `expect_zero(Lam2_{fail,suff}.subs(...) - (KX/chi)*G_{fail,suff})` | claim 3 | yes |
| A4 | sympy | 73-74 | `expect_zero(Delta0_k0 - 1/2)`, `(Delta_inf_k0 - 1)` | claim 4 | yes |
| A5 | sympy | 76-79 | `expect_zero(lim kappa*G_fail - Pe_req)`, `(lim kappa*G_suff - 2 Pe_req)` | claim 4 | yes |
| A6 | sympy | 86-89 | `expect_zero(Delta0_eta_inf - (1-sech)/kappa)`, `(Delta_inf_eta_inf - tanh/alpha)` | claim 5 | yes |
| A7 | sympy | 95-98 | `expect_zero(G_fail_inf - Pe_req/(alpha tanh))`, `(G_suff_inf - Pe_req/(1-sech))` | claim 5 | yes |
| A8 | sympy | 104-107 | `expect_zero(lim z*G_fail_inf_z - Pe_req)`, `(lim G_suff_inf_z - Pe_req)` | claim 6 | yes |
| B1 | math | 48 | `expectZero[Xi_micro/.tX->kX ell^2/kappa - kappa gMicro]` | claim 1 | yes |
| B2 | math | 54-57 | `expectZero[{chiFail,chiSuff,lambda2Fail,lambda2Suff}/.... - ...gFail/gSuff]` | claim 3 | yes |
| B3 | math | 63-66 | `expectZero[delta0K0-1/2]`,`[deltaInfK0-1]`,`[lim kappa gFail-peReq]`,`[lim kappa gSuff-2peReq]` | claim 4 | yes |
| B4 | math | 77-80 | `expectZero[delta0EtaInf-(1-Sech)/kappa]`,`[deltaInfEtaInf-Tanh/alpha]`,`[gFailInf-...]`,`[gSuffInf-...]` | claim 5 | yes |
| B5 | math | 85-86 | `expectZero[lim z gFailInfZ-peReq]`,`[lim gSuffInfZ-peReq]` | claim 6 | yes |

No tautological rows. Each `expect_zero`/`expectZero` subtracts an independently-built quantity from an independently-stated closed form; the printed residuals being `0` while the *printed* intermediate forms differ from the asserted forms (e.g. SymPy prints `G_suff^(inf)=Pe_req*cosh/(cosh-1)` then asserts equality with `Pe_req/(1-sech)`; Mathematica prints `gFailInf=peReq*Coth/Sqrt[kappa]` then asserts equality with `peReq/(alpha Tanh)`) confirms the simplifier actually performs work — the checks can fail.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.txt` (mtime 2026-05-11 12:44:14)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage061_microscopic_gain_thresholds_mathematica_audit.txt` (mtime 2026-05-11 12:54:35)

**What's wrong:**
Both saved `.txt` outputs predate their scripts. The `.py` (mtime 2026-06-03 15:59:11) and `.wl` (mtime 2026-06-03 15:59:11) were both edited after the captured transcripts. The content also disagrees on the stage banner: the current SymPy script emits `banner("STAGE 61 — MICROSCOPIC GAIN THRESHOLDS")` (py:22) and `banner("STAGE 61 THEOREM LEDGER")` (py:109), but the saved transcript still shows `STAGE 44 — MICROSCOPIC GAIN THRESHOLDS` (txt:11) and `STAGE 44 THEOREM LEDGER` (txt:46). The current Mathematica script emits `banner["STAGE 061 — MICROSCOPIC GAIN THRESHOLDS"]` (wl:26), but the saved transcript shows `STAGE 044 — MICROSCOPIC GAIN THRESHOLDS` (txt:11). All numeric/symbolic result lines and every `PASS:` line in the saved transcripts match what the current scripts would produce; only the stage-number banner labels are stale. This is the known pass-2 "outputs predate the banner fix" signal — the math is current, the label is not.

**Why this matters:**
Low impact: no result value is wrong, and exit code is 0 in both transcripts. But the banner label drift means the committed transcript misidentifies the stage, and a future reader diffing transcript-vs-script sees a spurious mismatch. The orchestrator's independent re-run refreshes these.

**Required change:**
No script edit required. The orchestrator/verifier re-runs `python3 scripts/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.py` and `math -script mathematica/moving_throat_pde_stage061_microscopic_gain_thresholds_mathematica_audit.wl` and overwrites the two `.txt` files. After refresh, the banners read `STAGE 61`/`STAGE 61 THEOREM LEDGER` (sympy) and `STAGE 061` (mathematica).

**Verification:**
After re-run, sympy output line 11 reads `STAGE 61 — MICROSCOPIC GAIN THRESHOLDS` and line 46 reads `STAGE 61 THEOREM LEDGER`; mathematica output line 11 reads `STAGE 061 — MICROSCOPIC GAIN THRESHOLDS`. Both still `EXIT_CODE: 0`.

## Independent-derivation check (Mathematica)

Independent (acceptable). The stage is intrinsically an algebraic restatement plus limits of upstream-quoted endpoint formulas (Stage 058 `Delta_0`/`Delta_inf`), so both engines must quote the same closed forms and apply the same substitutions — that shared structure is dictated by the math, not by transliteration. Evidence of genuine independence: (a) the two engines parametrize the stiff-support limit differently — SymPy substitutes `alpha -> z` (py:102) while Mathematica substitutes `kappa -> z^2` (wl:83) — yet both reach `-> Pe_req`; (b) the engines emit *different* intermediate closed forms for the same quantity and each verifies against the asserted target independently, e.g. for `G_fail^(inf)` SymPy prints `Pe_req/(sqrt(kappa)*tanh(sqrt(kappa)))` (txt:38) while Mathematica prints `peReq*Coth[Sqrt[kappa]]/Sqrt[kappa]` (txt:49); for `G_suff^(inf)` SymPy prints `Pe_req*cosh/(cosh-1)` (txt:39) while Mathematica prints `peReq*(1+(-1+Cosh[Sqrt[kappa]])^(-1))` (txt:50). A line-by-line port would not produce these divergent intermediate forms. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on every deliverable. Final residuals are `0` for all asserted identities in both transcripts (sympy txt lines 17,24-27,30-33,36-37,40-43; mathematica txt lines 19-28,33-46,51-62), and both exit 0. The printed `Delta_0`, `Delta_inf`, `G_fail`, `G_suff`, `Xi_micro`, `G_micro` forms are algebraically identical across engines (mathematica writes `kappa^(3/2)` where sympy writes `sqrt(kappa)*sinh(...)`/`alpha^2*...`; these are the same expression). No `engine_disagreement`.

## Verdict justification

`verdict: findings` with a single low-severity `stale_output`. Attacks that failed: (1) tautology — every `expect_zero` subtracts an independently-constructed quantity from an independently-stated closed form, and the printed-vs-asserted-form divergence proves the simplifier does real work, so none is guaranteed-zero by construction; (2) hardcoded result — no literal numeric answer is pinned; `1/2`, `1`, `2*Pe_req`, `Pe_req` are derived as limits, not asserted as inputs; (3) wrong-branch — the soft-support and compliant-mouth limits use the correct one-sided/`Infinity` directions and the printed limits (`1/2`, `1`, `tanh/sqrt(kappa)`, `(1-sech)/kappa`) match the notes exactly; (4) paper misalignment — all six deliverables map to non-tautological checks in both engines, with `G_fail`/`G_suff` and `G_micro` matching the boxed card equations and the notes; (5) symbol-domain — all symbols are `positive=True, real=True` (py:25-26) / `> 0` Reals (wl:30-32), consistent with the physical setup (gains, stiffnesses, Peclet number all positive). The only defect is the stale banner label in the committed transcripts; the math content is current and correct, so the finding is informational and non-blocking. I read the paper card, the notes file, and the appendix row before opening the scripts; the scripts' claims match the paper's claims.

## Self-test notes

Checked: variable-independence — no `sp.diff`/`D[]` derivatives appear in this unit; all checks are algebraic identities and symbolic limits, so the zero-derivative trap does not apply. Symmetry/parity — no unbounded integrals; not applicable. Trivial-case — substituting the soft-support concrete reductions (`Delta_0->1/2`, `Delta_inf->1`) into `kappa*G_fail`/`kappa*G_suff` gives `Pe_req`/`2 Pe_req`, matching the asserted residual-zero; the `eta->oo` and stiff-support limits likewise reduce to the asserted closed forms. Paper round-trip — the single finding (stale_output) prescribes no script edit, so it cannot introduce a new paper_misalignment.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 12 values checked, 0 misaligned.

All emitted deliverable values are symbolic closed forms (this stage pins no numeric benchmark constants beyond the limiting literals `1/2`, `1`, `2 Pe_req`, `Pe_req`, which are themselves derived limits). Each is carried in the paper card and/or notes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Xi_micro = kappa*G_micro`, `G_micro = chi_sigma Lambda_phi^2/K_X` | py:38-39 / wl:45-48; sympy txt:15-17, math txt:15-16,19 | tex:14-17 (eq Xi-micro), md:100-108 | MATCH |
| `G_fail = Pe_req/[kappa Delta_inf]` | py:46,48 / wl:39; sympy txt:18, math txt:17 | tex:19-25 (boxed), md:124 | MATCH |
| `G_suff = Pe_req/[kappa Delta_0]` | py:47,49 / wl:40; sympy txt:19, math txt:18 | tex:19-25 (boxed), md:126 | MATCH |
| phase diagram (fail/succeed/root-sensitive regions) | py:111 (ledger) / wl banner; sympy txt:49 | tex:27-35 (boxed), md:130-132 | MATCH |
| `chi_fail`/`chi_suff` threshold surfaces | py:51-52,55-56 / wl:50-51; sympy txt:20-21 | md:144-146 | MATCH (notes carrier; card terse) |
| `Lambda^2_fail`/`Lambda^2_suff` threshold surfaces | py:53-54,57-58 / wl:52-53; sympy txt:22-23 | md:150-152 | MATCH (notes carrier) |
| soft-support `Delta_0->1/2`, `Delta_inf->1` | py:71-74 / wl:61-64; sympy txt:28-29, math txt:33-34 | md:172-174 | MATCH |
| soft-support `G_fail~Pe_req/kappa`, `G_suff~2Pe_req/kappa` | py:76-79 / wl:65-66; sympy txt:32-33 | md:178-180, ledger md:60-62 | MATCH |
| `Delta_0^(inf)=(1-sech)/kappa`, `Delta_inf^(inf)=tanh/sqrt(kappa)` | py:84-89 / wl:73-78; sympy txt:34-35, math txt:47-48 | md:194-196 | MATCH |
| `G_fail^(inf)=Pe_req/[sqrt(kappa)tanh]`, `G_suff^(inf)=Pe_req/[1-sech]` | py:93-98 / wl:75-80; sympy txt:38-39, math txt:49-50 | md:200-202 | MATCH |
| stiff-support `sqrt(kappa)G_fail->Pe_req` | py:104,106 / wl:85; sympy txt:42, math txt:59 | md:232, 78 | MATCH |
| stiff-support `G_suff->Pe_req` | py:105,107 / wl:86; sympy txt:43, math txt:61 | md:234, 78 | MATCH |

INTERNAL (scaffolding, no finding): `alpha=sqrt(kappa)`, `Xi_fail=Pe_req/Delta_inf`, `Xi_suff=Pe_req/Delta_0` (intermediate), the substitution `T_X=K_X L^2/kappa`, the `z`/`alpha->z`/`kappa->z^2` change of variable, and all `... = 0` residual print lines and `PASS:` flags.

No MISMATCH and no MISSING-DELIVERABLE. Reconciliation is consistent with the standard verdict; the sole finding remains the informational `stale_output`.
