---
unit_id: 047
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage047_coherent_kernel_map.md]
  paper_appendix: present
---

# Audit unit 047 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_047.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage047_coherent_kernel_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 72)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.txt`

## What the paper claims

Stage 047 compresses the first coherent local D/N kernel (from Stages 045/046) into a dimensionless placement map and isolates how the support lane enters. The two boxed deliverables in the card are: (1) the **support-enhancement factor**
`\boxed{S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\epsilon\zeta}}` (eq. `app-stage047-S`), entering through `M_tr = M_mix S(\zeta;\epsilon)` (eq. `app-stage047-Mtr`); and (2) the **target ratio**
`\boxed{R_{\rm target}=\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}}` (eq. `app-stage047-Rtarget`), which contains no `\zeta` — so the support lane supplies loading but does not move the target locus. The `\stagefield{Output}` is verbatim: "Support enhancement factor \eqref{eq:app-stage047-S} and target independence \eqref{eq:app-stage047-Rtarget}." The notes add the full supporting map: the dimensionless ratio definitions (`eps_eta, eps_W, Z_W, delta_0, chi_0, zeta, Lambda`), the proportionality identities `eps_phi = zeta eps_W`, `Z_phi = zeta Z_W`, the split ratios (`eps = eps_W[1-(2/11)delta_U/(1+delta_U)]`, `R_tr`, `delta`), the baselines `M_mix`, `M_supp`, monotonicity `dS/dzeta=(1-eps)/(1-zeta eps)^2>0`, `S(0;eps)=1`, and the product law `R_target M_tr = 8 Lambda (1-eps)/pi^2 * S`.

## What the script claims to verify

Both engines build the dimensionless ratios from the microscopic definitions in notes §2 (`chi_0, eps_eta, eps_W, Z_W, zeta_def, delta_0, delta_U, Lambda, eps_phi, Z_phi`), then assert: (a) the proportionality identities `eps_phi - zeta_def*eps_W = 0` and `Z_phi - zeta_def*Z_W = 0`; (b) the factorization `M_tr - M_mix*S = 0` where `M_tr = M_mix + M_supp` with `M_mix`, `M_supp`, `S` independently constructed from their closed forms; (c) the product law `R_target*M_tr - 8 Lambda (1-eps)/pi^2 * S = 0`; (d) ζ-independence of `R_target` via `dR_target_loaded/dzeta = 0`, `R_target_loaded(zeta)-R_target_loaded(0)=0`, and `support-loaded R_target reconstruction - R_target = 0`; (e) the monotonicity-derivative identity `dS/dzeta - (1-eps)/(1-zeta eps)^2 = 0` and `S(zeta=0)-1=0`. The sympy script adds a numeric **negative control** that spoils the support-loading coefficient and confirms ζ-blindness then breaks (nonzero derivative `-44.866…`). The `chi_0 = sigma_0 = rho_0` rename is deliberately NOT asserted (commented as tautological, deferred to upstream Stage 28).

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| `S(\zeta;\epsilon)=1+\zeta(1-\eps)/(1-\eps\zeta)` (eq. -S) | `S` def line 91 (py)/89 (wl); `S(zeta=0)-1=0`; `dS/dzeta-(1-eps)/(1-zeta eps)^2=0`; monotone via positive numerator | match |
| `M_tr = M_mix S(\zeta;\epsilon)` (eq. -Mtr) | `expect_zero("M_tr - M_mix*S")` line 101 (py)/99 (wl), with `M_tr=M_mix+M_supp` independently built | match |
| `R_{\rm target}=\Lambda(1-\eps_\eta)(1-\eps)^2/[Z_W(1+\chi_0)^2]` (eq. -Rtarget) | `Rtarget` def line 93 (py)/92 (wl); reconciled to product law via `support-loaded R_target reconstruction = 0` | match |
| `R_target` independent of `\zeta` (Output) | `dR_target_loaded/dzeta=0`, `R_target_loaded(zeta)-R_target_loaded(0)=0` lines 113,114-117 (py)/108,109-112 (wl) | match |
| product law `R_target M_tr = 8\Lambda(1-\eps)/\pi^2 S` (notes §5) | `expect_zero("product law")` line 111 (py)/106 (wl) | match (notes-level) |
| proportionality `eps_phi=\zeta eps_W`, `Z_phi=\zeta Z_W` (notes §2) | lines 75,76 (py)/64,65 (wl) | match (notes-level) |
| split `eps`, `R_tr`, `delta` (notes §3) | computed + printed (lines 79-81 py / 67-69 wl); no standalone assert, but fed into the asserted M/R identities | partial (printed, exercised indirectly) |

`paper_alignment: aligned` — every boxed/Output deliverable has a faithful, non-tautological script-side check; both engines agree.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 75 | `eps_phi - zeta_def*eps_W == 0` | proportionality (notes §2) | yes |
| A2 | sympy | 76 | `Z_phi - zeta_def*Z_W == 0` | proportionality (notes §2) | yes |
| A3 | sympy | 101 | `M_tr - M_mix*S == 0` | eq. -Mtr / -S | yes |
| A4 | sympy | 111 | `R_target*M_tr - 8Λ(1-eps)/π²·S == 0` | product law (notes §5) | yes |
| A5 | sympy | 112 | `R_target_loaded - R_target == 0` | eq. -Rtarget | yes |
| A6 | sympy | 113 | `dR_target_loaded/dzeta == 0` | ζ-independence (Output) | yes |
| A7 | sympy | 114-117 | `R_target_loaded(zeta) - R_target_loaded(0) == 0` | ζ-independence (Output) | yes |
| A8 | sympy | 118 | `dS/dzeta - (1-eps)/(1-zeta eps)² == 0` | monotonicity (notes §4) | yes |
| A9 | sympy | 119 | `S(zeta=0) - 1 == 0` | `S(0;eps)=1` (notes §4) | yes |
| A10 | sympy | 155-156 | spoiled `dR_target/dzeta` numerically nonzero | negative control for ζ-blindness | yes |
| B1-B9 | mathematica | 64,65,99,106,107,108,109-112,113,114 | mirror of A1-A9 (independent FullSimplify) | same | yes |

All load-bearing assertions are non-tautological: each combines independently-built closed forms (`M_mix`, `M_supp`, `S`, `R_target`, `product_expected`) and asks whether two distinct constructions coincide. The factorization A3 is the strongest — `M_tr` is built as a *sum* of two baselines while `M_mix*S` is built from a *product*; their equality is the genuine algebraic content. A6/A7 are not tautological because `R_target_loaded = product_expected/M_tr` and both `product_expected` (via `S`) and `M_tr` carry `zeta`; the cancellation to a ζ-free expression is the real claim.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.txt` (mtime 2026-05-26 02:04:43)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.txt` (mtime 2026-05-26 02:04:49)
- vs scripts mtime 2026-06-03 15:59:11 (both `.py` and `.wl`).

**What's wrong:**
Both committed output transcripts predate the current scripts by ~8 days. The staleness is also visible in content: the sympy transcript banner reads `STAGE 30 — COHERENT KERNEL MAP AUDIT` (line 3) and closes `All Stage-30 symbolic checks passed.` (line 53); the Mathematica transcript banner reads `STAGE 030 — COHERENT KERNEL MAP` (line 4) and closes `Stage 047 Mathematica audit passed.` (line 40, inconsistent within the same file). The current scripts print `STAGE 47` (py line 32) / `STAGE 047 — COHERENT KERNEL MAP` (wl line 26). So the transcripts were captured against an earlier banner state. The *substantive* results in the transcripts are nonetheless consistent with what the current scripts would produce (all identities `= 0`, PASS lines present, spoiled control `-44.866…` nonzero), so this is informational, not a math defect.

**Why this matters:**
A reader auditing the committed transcript sees "STAGE 30" / "STAGE 030" labels that contradict the stage-047 filename, undermining traceability. The orchestrator's independent re-run will regenerate fresh transcripts; this finding records that the committed ones are stale.

**Required change:**
No script-math change. Re-run both scripts to refresh the two `output/*.txt` transcripts so banners read STAGE 047 and the close lines agree. (This is the orchestrator's standard refresh, not a Codex edit.)

**Verification:**
After refresh, `scripts/output/...sympy_audit.txt` line 3 reads `STAGE 47 ...` and the close line reads `All Stage-47 ...` (or whatever the current script prints), and the mtime is newer than the `.py`. Same for the `.wl` transcript.

NOTE on stale self-labels in the LIVE scripts (not a separate finding here — deferred per the script/output-band numbering plan): the live `.py` still carries `"""Stage 30 SymPy audit.` (line 2), `print("\nAll Stage-30 symbolic checks passed.")` (line 158), and comments "Stage 28" (line 45) / "Stage-30 support-loading" (line 121); the live `.wl` carries comment "Stage 28" (line 41). These are label-only numbering drift on a clean stage — exactly the class the user gated to the dedicated `NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` pass (content-keyed, never offset-swept). Flagged here for the band tracker; not folded into F1's required change to avoid pre-empting that deferred pass.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration. It re-declares all symbols under a single global `$Assumptions` block (Reals + positivity), uses `FullSimplify[Together[Expand[...]]]` as its zero-test, and constructs each closed form (`mMix`, `mSupp`, `sEnhance`, `rTarget`, `productExpected`) directly from the dimensionless definitions — independently from the SymPy `sp.simplify` route. Corresponding sections justify independence:
- SymPy `expect_zero` (py 25-29) raises a Python `AssertionError` on nonzero; Mathematica `expectZero` (wl 20-24) calls `Exit[1]` via `fail[]` — different control flow, different simplifier (`FullSimplify` vs `simplify(expand())`).
- SymPy builds `Mtr = sp.simplify(Mmix + Msupp)` (py 92) and tests `Mtr - Mmix*S`; Mathematica builds `mTr = FullSimplify[mMix + mSupp]` (wl 86) and tests `mTr - mMix*sEnhance` — same physical identity, independent algebra. The two transcripts even print `M_tr` in *different normal forms* (py line 36 fully cleared denominator; wl line 21 as a sum of two reciprocal terms), confirming independent simplification paths that nonetheless both yield `0` for the difference.
- The `.wl` omits the SymPy-only numeric negative control (py 121-156); that is a permissible asymmetry (the control is a sanity probe, not a paper deliverable), not a transliteration tell.

No `mathematica_transliteration` finding.

## Engine cross-check

Both engines emit the same nine `= 0` / PASS identities (sympy txt lines 20,21,38,45,46,47,48,49,50; wl txt lines 11-13,23-24,27-38). The printed closed forms agree up to normal form: e.g. sympy `R_target` (txt line 37) and wl `R_target` (txt line 22) are the same rational function with sign absorbed into the `(KU*Keta_eff - c_etaU^2)` vs `(cEtaU^2 - kEtaEff*kU)` denominator orientation. The sympy-only spoiled control (`-44.866…`, txt line 51) has no Mathematica counterpart by design. `engines_agree: true`.

## Verdict justification

`verdict: findings` solely on the low-severity `stale_output` (F1); the math is clean. Attacks tried and failed: (1) I checked whether `M_tr - M_mix*S = 0` is a construction tautology — it is not, because `M_tr` is built as a sum and `M_mix*S` as a product, so their equality is real algebra (and the two engines reach it via different normal forms). (2) I checked whether the ζ-independence assertions are vacuous — `R_target` (line 93) literally has no `zeta`, but the load-bearing check is `R_target_loaded = product_expected/M_tr` (both ζ-bearing) reducing to ζ-free; that is non-trivial and additionally tied back to `R_target` via A5. (3) I checked the negative control — its `.subs` dict mixes expression-keys (`eps`, `eps_eta`, `chi0`, `ZW`, `Lambda`) with the underlying microscopic-symbol assignments; the expression-keys are largely inert after `sp.simplify` expands to microscopic symbols, but the microscopic assignments alone drive the nonzero `-44.866…`, so the control still fires correctly (minor code-smell, not a defect — it does not weaken any paper claim). (4) I confirmed the `chi_0=sigma_0=rho_0` non-assertion is correctly justified as a tautological rename deferred upstream, not a silently-dropped deliverable. (5) Paper alignment is exact: card eqs -S, -Mtr, -Rtarget, the Output ζ-independence line, the appendix row (line 72), and the notes product law/monotonicity all have faithful checks. I read the stage card, the notes file, and the appendix row in full before the scripts.

## Self-test notes

Traps checked: (1) **Variable independence** — `sp.diff(Rtarget_loaded, zeta)` and `D[loadedRTarget, zeta]`: `Rtarget_loaded` genuinely depends on `zeta` (through `S` in the numerator and `M_tr` in the denominator), so the derivative-is-zero is a real cancellation, not a trivial "var not present" zero; `sp.diff(S, zeta)` likewise depends on `zeta`. (2) **Symmetry/parity** — no unbounded integrals in this stage. (3) **Trivial-case** — `S(zeta=0)=1` checked by both engines and consistent with the closed form; the spoiled control gives a concrete nonzero literal (`-44.866…`), confirming the `assert_nonzero` path fires. No directive with a Codex edit is warranted (the only finding is informational stale_output handled by the orchestrator's refresh + the deferred band pass).

## Value Reconciliation (pass-2 augmentation)

All emitted values are symbolic closed forms (no pinned numeric constants except the structural rational `sigma = 88/(9 pi^2)` and the spoiled-control probe number). The deliverable-level reconciliation:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `S(\zeta;\eps)=1+\zeta(1-\eps)/(1-\zeta\eps)` | py:91 / wl:89; sympy txt:35, wl txt:20 (expanded form, equals closed form) | `.tex:21-23` (boxed eq -S); `.md:166-167` | MATCH |
| `M_tr = M_mix S(\zeta;\eps)` | py:92,101 / wl:86,99; txt `M_tr - M_mix*S = 0` (sympy:38, wl:23-24) | `.tex:14-17` (eq -Mtr); `.md:161-163,220` | MATCH |
| `M_mix = 8 Z_W(1+\chi_0)^2/[\pi^2(1-\eps_\eta)(1-\eps)]` | py:89 / wl:75; sympy txt:33, wl txt:18 | `.md:153` | MATCH (notes carrier) |
| `M_supp = 8\zeta Z_W(1+\chi_0)^2/[\pi^2(1-\eps_\eta)(1-\zeta\eps)]` | py:90 / wl:82; sympy txt:34, wl txt:19 | `.md:157` | MATCH (notes carrier) |
| `R_{target}=\Lambda(1-\eps_\eta)(1-\eps)^2/[Z_W(1+\chi_0)^2]` | py:93 / wl:92; sympy txt:37, wl txt:22 | `.tex:27-29` (boxed eq -Rtarget); `.md:189` | MATCH |
| product law `R_target M_tr = 8\Lambda(1-\eps)/\pi^2 \cdot S` | py:104,111 / wl:101,106; sympy txt:43, wl txt:25 | `.md:194-202` | MATCH (notes carrier) |
| `R_tr = [1+\chi_0/(1+\delta_U)]/(1+\chi_0)` | py:80; sympy txt:26, wl txt:15 | `.md:126,218` | MATCH (notes carrier) |
| `eps = eps_W[1-(2/11)\delta_U/(1+\delta_U)]` | py:79; sympy txt:27, wl txt:16 | `.md:136-137` | MATCH (notes carrier) |
| `delta = [\delta_0+\eps_\eta\delta_U/(1+\delta_U)]/(1-\eps_\eta)` | py:81; sympy txt:28, wl txt:17 | `.md:143` | MATCH (notes carrier) |
| `eps_phi=\zeta eps_W`, `Z_phi=\zeta Z_W` | py:75,76 / wl:64,65; sympy txt:20-21, wl txt:11-13 | `.md:116-117,139` | MATCH (notes carrier) |
| `dS/d\zeta=(1-\eps)/(1-\zeta\eps)^2` | py:118 / wl:113; sympy txt:49, wl txt:35 | `.md:177` | MATCH (notes carrier) |
| `\chi_0 = \gamma c_{\eta U}/K_U` | py:50 / wl:45; sympy txt:9, wl txt:5 | `.md:104` | MATCH (notes carrier) |
| `\Lambda = 27\pi^2 G c_s^5 K_W^{eff}/(20 a^5 c^5 \mu_W)` | py:63 / wl:55 | `.md:108` | MATCH (notes carrier) |
| `\sigma = 88/(9\pi^2)` | py:39 / wl:35 | `.md:88` | MATCH (notes carrier) |

INTERNAL (scaffolding, no prose expected): `chi0/eps_eta/epsW/ZW/zeta_def/delta0/deltaU/eps_phi/Zphi` intermediate symbol bindings used only to assemble the asserted identities; `product_expected`, `product_actual`, `Rtarget_loaded` reconstruction intermediates; `eta_bad` spoil parameter and the `-44.866…` negative-control probe number (verification scaffolding, not a deliverable); pass/fail flags and `= 0` residuals.

reconciliation: complete; 14 deliverable values checked, 0 misaligned. Every emitted deliverable reconciles to the card (boxed eqs -S/-Mtr/-Rtarget) or the notes §1-§5; no MISMATCH or MISSING-DELIVERABLE, so no `value_mismatch`/`script_missing_paper_claim` folded into Findings.
