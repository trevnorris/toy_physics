---
unit_id: 067
batch: III.3
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage067_sech_gaussian_resonance.md]
  paper_appendix: present
---

# Audit unit 067 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_067.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage067_sech_gaussian_resonance.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows at lines 112, 114, 252 reference this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage067_sech_gaussian_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.txt`

## What the paper claims

Stage 067 benchmarks the parent source/support coherence against an explicit independent-profile family: a sech source `chi_sigma = sech(y/w_f)` and a Gaussian support `chi_phi = exp(-y^2/w_g^2)`. The card's boxed deliverables are: (i) an exact self-dual stationary point `w_g/w_f = sqrt(pi)` (`eq:app-stage067-self-dual`); (ii) the maximal coherence there, `\boxed{C_{\rm res}^2 \simeq 0.994418836451529}` (`eq:app-stage067-Cres`); and (iii) the penalty factor `\boxed{P_{\rm res} = 1/C_{\rm res}^2 \simeq 1.005612487760576}` (`eq:app-stage067-Pres`). `\stagefield{Output}` reads: "The near-perfect benchmark coherence \eqref{eq:app-stage067-Cres} and penalty \eqref{eq:app-stage067-Pres}." The notes add the underlying structure the card omits: exact norms `N_(sigma sigma)=2 w_f`, `N_(phi phi)=w_g sqrt(pi/2)`, the closed form `C^2(r)=I(r)^2/[r sqrt(2 pi)]`, the exact overlap duality `I(r)=(r/sqrt(pi)) I(pi/r)` ⇒ `C^2(r)=C^2(pi/r)`, the shortfall `1-C_res^2 ≈ 0.56%`, and the explicit caveat that this benchmark does NOT prove threshold survival (which still depends on Stages 061–066). The appendix row (line 112) summarizes: "Coherence resonance at `w_g/w_f=sqrt(pi)`."

## What the script claims to verify

Both scripts (1) derive the exact norms by direct symbolic integration and confirm them against the declared `2 w_f` / `w_g sqrt(pi/2)`; (2) verify the *algebraic* implication that the overlap duality `I(r)=(r/sqrt(pi))I(pi/r)` forces `C^2(r)=C^2(pi/r)` (honestly scoped in comments as holding for any `I`); (3) verify that a function symmetric under `r<->pi/r` has zero derivative at `r_*=sqrt(pi)` (honestly flagged as a calculus tautology, with a non-tautological falsification guard that a perturbed slope breaks stationarity); (4) numerically evaluate the real non-elementary overlap to produce `C_res^2`, `P_res`, `1-C_res^2`; (5) numerically verify the sech-Gaussian-specific duality at sample points; (6) numerically verify strict monotone increase up to `r_*` and decrease after, establishing the global maximum on the constructive branch. The `.wl` additionally cross-checks its `NIntegrate` `C_res^2`/`P_res` against the sympy mpmath values (honestly labeled as cross-engine numerical agreement, not closed-form).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Exact norms `N_ss=2 w_f`, `N_pp=w_g sqrt(pi/2)` (notes §1) | py L70-75 / wl L45-53: direct integration + equality check | match |
| Closed form `C^2(r)=I(r)^2/[r sqrt(2 pi)]` (notes §1) | py L79-80 / wl L56-57: prints simplified form | match |
| Exact duality `C^2(r)=C^2(pi/r)` (notes §2; appendix) | py L88-94 (algebraic) + L177-184 (numeric) / wl L61-67 + L139-147 | match |
| Self-dual stationary point `r_*=sqrt(pi)` (eq self-dual; notes §3) | py L109-144 (calculus, self-flagged tautological) + monotonicity scan L190-209 / wl L71-96 + L149-171 | match |
| `C_res^2 = 0.994418836451529` (eq Cres; Output) | py L163-167 / wl L119,134: numeric eval | match |
| `P_res = 1.005612487760576` (eq Pres; Output) | py L164-168 / wl L120,135 | match |
| `1-C_res^2 ≈ 0.56%` (notes §4) | py L169 / wl L132 | match |
| "does NOT prove threshold survival" caveat (notes §5-6) | py L223-225 / wl interpretation print | match (interpretive) |

All paper-side deliverables have a faithful script-side check; no `mismatch`, no `missing`, no `extra`. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 74 | `expect_zero(Nss_integral - 2 w_f)` | norm N_ss | yes (integral derived in-script) |
| A2 | sympy | 75 | `expect_zero(Npp_integral - w_g sqrt(pi/2))` | norm N_pp | yes |
| A3 | sympy | 94 | `expect_zero(C2_dual - C2_target)` | duality (algebraic) | partial (any-`I`; self-labeled; numeric backs it) |
| A4 | sympy | 118-121 | `expect_zero(self-dual overlap-slope)` | stationary point | no (self-flagged tautological calculus) |
| A5 | sympy | 131-134 | `expect_zero(dC2 at self-dual)` | stationary point | no (self-flagged tautological calculus) |
| A6 | sympy | 141-142 | `if dC2_broken==0: raise` (perturbation guard) | stationary point | yes (falsification: perturbed slope ≠ 0) |
| A7 | sympy | 183-184 | `if diff > 1e-40: raise` (numeric duality) | sech-Gaussian duality | yes (real overlap integral) |
| A8 | sympy | 198-200 | strict-increase monotonicity | r_* is max | yes (real C^2 evals) |
| A9 | sympy | 207-209 | strict-decrease monotonicity | r_* is max | yes |
| B1 | mathematica | 52-53 | `expectZero(nssDirect - 2 wf / nppDirect - ...)` | norms | yes (Integrate in-script) |
| B2 | mathematica | 67 | `expectZero(c2Dual - c2Target)` | duality (algebraic) | partial (any-`I`; self-labeled) |
| B3 | mathematica | 93-96 | `expectZero(c2PrimeLeft from Solve)` | stationary point | no (self-flagged tautological) |
| B4 | mathematica | 134-135 | `expectApprox(C_res^2 / P_res vs sympy target)` | C_res^2, P_res cross-engine | yes (NIntegrate vs mpmath) |
| B5 | mathematica | 145 | `expectTrue(duality diff <= 1e-35)` | sech-Gaussian duality | yes |
| B6 | mathematica | 157-160 | `expectTrue(increase up to r_*)` | r_* is max | yes |
| B7 | mathematica | 168-171 | `expectTrue(decrease after r_*)` | r_* is max | yes |

The non-anchored rows (A4/A5/B3) are explicitly self-labeled in the script comments as calculus tautologies whose substantive counterpart is the numerical monotonicity scan (A8/A9/B6/B7), which IS anchored. A3/B2 are honestly scoped as any-`I` algebraic implications backed by the numeric duality A7/B5. The load-bearing physics (the actual sech-Gaussian overlap, its duality, and the global max at `r_*`) is exercised by real numerical integration in both engines — non-tautological.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage067_sech_gaussian_sympy_audit.txt:3`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.txt:3`

**What's wrong:**
Both committed `.txt` outputs predate their scripts and carry a stale stage-number banner. mtimes:
- sympy `.py` mtime `2026-05-26 13:11:12`; sympy `.txt` mtime `2026-05-22 19:56:14` (output older).
- mathematica `.wl` mtime `2026-05-26 13:11:29`; mathematica `.txt` mtime `2026-05-22 19:56:29` (output older).

Content corroborates the staleness. The current sympy script (L53) prints `banner("STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK")` (hyphen), but the saved output L3 reads `STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE BENCHMARK` (en-dash, "50"). The current mathematica script (L38) prints `banner["STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"]`, but the saved output L3 reads `STAGE 050 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` ("050"). Note the mathematica output's *final* line (L78) already reads "Stage 067 Mathematica audit passed", so the capture caught the corrected final label but not the corrected top banner — consistent with the banner having been edited 50/050→067 after the `.txt` was committed. The numeric/symbolic results in the stale outputs still match what the current scripts would emit (only the banner label changed), so this is informational, not a math defect.

**Why this matters:**
The committed transcript misrepresents the stage number ("50"/"050") it belongs to, and does not reflect the current banner. It is the standard pass-2 script/output-band staleness signal.

**Required change:**
No script edit. The orchestrator's independent re-run (`exec-sympy 067` / `exec-mathematica 067`) refreshes both `.txt` files; the refreshed banner will read `STAGE 067` and the mtimes will post-date the scripts.

**Verification:**
After refresh, both `.txt` L3 banners read `STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK`, output mtimes exceed script mtimes, and all PASS lines / numeric values are unchanged (C_res^2 = 0.994418836451529348…, P_res = 1.00561248776057621695…).

## Independent-derivation check (Mathematica)

Not a transliteration. Section 3 (stationary point) uses a genuinely different construction in each engine: the `.py` (L109-126) hand-writes the tangent `Iprime_left - (Istar/sqrt(pi) - sqrt(pi)*Iprime_dual/r)` and a hand-written `dC2_selfdual` formula, then substitutes the slope; the `.wl` (L75-96) instead applies real symbolic differentiation `D[c2Fn[r] - c2Fn[Pi/r], r]` with `Derivative[1]` replacement and `Solve[... == 0, c2PrimeLeft, Reals]` to *solve* for the stationary slope. The numerical sections use each engine's native quadrature — mpmath `mp.quad` over `[-inf, inf]` (`.py` L154-156) vs. Mathematica `2*NIntegrate[..., {x,0,Infinity}, WorkingPrecision->80]` exploiting evenness (`.wl` L101-114). The `.wl`'s L126-135 explicit cross-check of its `NIntegrate` result against the sympy mpmath value is honestly labeled (L122-125) as cross-engine numerical agreement on the same definite integral, not a closed-form benchmark; that is legitimate corroboration, not an echo of the `.py` algebra. Independent.

## Engine cross-check

The engines agree. Symbolic side: both derive `N_ss=2 w_f`, `N_pp=w_g sqrt(pi/2)`, `C^2(r)=I(r)^2/[r sqrt(2 pi)]`, and the algebraic duality `C^2(r)=C^2(pi/r)` (both `=0` residuals). Numeric side: sympy `C_res^2 = 0.994418836451529348706428351608877628170873348983716948813464`; mathematica `C_res^2 = 0.99441883645152934870642835160887762817087334898371694998969514…` — agree to ~32 significant digits, well inside the `.wl`'s `10^-35`/`10^-34` cross-engine tolerances (mpmath dps60 vs. NIntegrate WP80 with AccuracyGoal/PrecisionGoal 32). `P_res`, `1-C_res^2`, and all monotonicity-grid `C^2` values agree to displayed precision. No sign/factor disagreement.

## Verdict justification

The scripts faithfully and non-tautologically exercise every paper-side deliverable. Attacks tried that failed: (1) I checked whether the stationary-point assertions are vacuous — they ARE self-flagged calculus tautologies (A4/A5/B3), but the substantive maximum claim is independently carried by the real-integral monotonicity scans (A8/A9/B6/B7) and the perturbation falsification guard (A6), so the deliverable is genuinely exercised, not faked. (2) I checked the duality assertion for being an empty any-`I` algebraic identity — it is, but it is honestly scoped and backed by the numeric sech-Gaussian duality (A7/B5). (3) I checked the norms for hardcoding — they are derived by direct in-script integration and checked against the declared forms, so a wrong norm would fail. (4) I checked symbol domains — all widths/`r` declared `positive=True`/`>0`, matching the physical setup; the Parseval/Fourier duality the notes invoke is exercised numerically rather than assumed symbolically. (5) I verified engine independence (different section-3 constructions, native quadratures). The only defect is the stale committed transcripts carrying a "STAGE 50"/"STAGE 050" banner versus the corrected "STAGE 067" scripts — a low-severity `stale_output` resolved by the orchestrator's independent re-run, with no math change. Paper alignment is exact; all deliverable values reconcile (see below). Verdict: `findings` (one informational stale_output).

## Value Reconciliation (pass-2 augmentation)

`reconciliation: complete; 8 values checked, 0 misaligned`

All emitted deliverable values appear in the `.tex` card and/or `.md` notes and agree to displayed precision. The paper `.tex` quotes truncated/rounded forms of the higher-precision script values, which is expected for a card; the notes carry the longer forms, and the scripts carry the full precision.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `r_* = w_g/w_f = sqrt(pi)` | py L104,166 / wl L71,118; sympy out L35 | `.tex:17` `w_g/w_f=\sqrt\pi`; `.md:104` `w_g/w_f = sqrt(pi)` | MATCH |
| `N_(sigma sigma) = 2 w_f` | py L62,213 / wl L47; sympy out L9 | `.md:45` `N_(sigma sigma) = … = 2 w_f` | MATCH |
| `N_(phi phi) = w_g sqrt(pi/2)` | py L63,214 / wl L48; sympy out L10 | `.md:47` `N_(phi phi) = w_g sqrt(pi/2)` | MATCH |
| `C^2(r) = I(r)^2/[r sqrt(2 pi)]` | py L79,215 / wl L56; sympy out L15 | `.md:60-61`; `.tex` (implicit via boxed values) | MATCH |
| duality `I(r)=(r/sqrt(pi))I(pi/r)` ⇒ `C^2(r)=C^2(pi/r)` | py L94,216-217 / wl L67; sympy out L20-21 | `.md:82,86` | MATCH |
| `C_res^2 = 0.994418836451529348…` | py L163,167,221 / wl L119,134; sympy out L36 | `.tex:23` `\simeq0.994418836451529`; `.md:117` `0.9944188364515293487…` | MATCH (card truncated) |
| `P_res = 1.005612487760576217…` | py L164,168,222 / wl L120,135; sympy out L37 | `.tex:29` `\simeq1.005612487760576`; `.md:128` `1.0056124877605762169…`; appendix L114 `\simeq1.005612` | MATCH (card/appendix truncated) |
| `1 - C_res^2 = 0.005581163548470651… (≈0.56%)` | py L169 / wl L132; sympy out L38 | `.md:121,123` `0.0055811635484706513…`, "about 0.56%" | MATCH |

INTERNAL (scaffolding, no finding expected in prose): `dC2_broken` perturbation residual, the `Iprime_left/Istar/delta_bad` symbolic slope tokens, monotonicity-grid sample `C^2(x)` values (0.55/0.7/…/4.0), per-sample numeric-duality `|I(r)-…|` diffs, cross-engine `expectApprox` diffs, all PASS/FAIL flags, tolerances (`1e-40`, `10^-35`, `10^-34`), and precision/grid parameters.

## Self-test notes

Variable-independence: no zero-by-construction derivatives — the only symbolic `D[]` is `.wl` L76 `D[c2Fn[r]-c2Fn[Pi/r], r]`, where `c2Fn` genuinely depends on `r`, so the derivative is non-trivial; the `.py` uses hand-written tangents, not `sp.diff` on a constant. Symmetry/parity: the overlap integrand `sech(x) exp(-x^2/r^2)` is even, so `mp.quad` over `[-inf,inf]` and `2*NIntegrate[..,{x,0,inf}]` correctly agree (the `.wl`'s factor-2 evenness shortcut is valid and matches the `.py` full-line integral to ~32 digits). Trivial-case/anchoring: the norm checks integrate the actual profiles in-script (not hardcoded), the monotonicity/duality checks use the real non-elementary overlap, and the perturbation guard (A6) confirms a nonzero slope breaks stationarity — so the substantive claims are exercised, while the section-3 calculus identities are honestly self-labeled tautologies backed by the numerics.
