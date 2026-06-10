---
unit_id: 222
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
  notes_stage_files: [moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 222 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_222.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows/eqs 55-56, 518-572, 679, 1433)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_mathematica_audit.txt`

## What the paper claims

Stage 222 inserts the lowest concrete finite-throat primitive branch into the resonance/linewidth theorem. `\stagefield{Output}`: "Concrete branch pole census, residue/linewidth cancellation, and first branch-level dynamic survival test. The pole polynomial and survival gate are exact in the primitive closure; the displayed slice is a numerical placement." The notes enumerate six deliverables: (1) the primitive branch from the lowest N/N wall mode `u_0=1/sqrt(L)` and lowest D/N half-wave `f_0=sqrt(2/L) sin(pi s/2L)` with exact overlap `kappa=2sqrt(2)/pi`; (2) the exact quartic pole polynomial `F(y)` with `D(omega)=F(omega^2)/((varpi^2-omega^2)Delta(omega))`; (3) the exact residue/linewidth cancellation `R_{Q,*}=27 c_s^5/(a^5 omega_*^5 N(omega_*))`; (4) the exact low-loss survival inequality `R_{Q,*} >= 2 DeltaV_req (1+eta^2)/eta x^6`; (5) one explicit numerical pole census on the admissible sample slice; (6) the first static/dynamic tension: stronger outgoing coupling lambda_W raises P0 but lowers the upper-wall R_{Q,*}. The notes also report a concrete numerical census table (sec. 6.2) and the lambda_W scan table (sec. 8). The published card and Part-VII appendix carry only the abstract/symbolic forms (eqs. app-part07-primitive-kappa, -quartic-pole, -residue-linewidth-primitive); the numeric slice lives in the notes (this is the natural carrier and is consistent with the card's `\StatusNumerical{}` flag).

## What the script claims to verify

The SymPy script verifies, in order: the exact overlap constant `kappa=2sqrt(2)/pi` (symbolic integral, line 56); the exact quartic identity `D - F(omega^2)/((varpi^2-omega^2)Delta)==0` and `deg_y F = 4` (lines 97-101); the residue/linewidth cancellation `R_{Q,*}=27 c_s^5/(a^5 omega_*^5 N_star)` (line 117); the low-loss peak and threshold closed forms (lines 119-124); the static slice data (Delta0/D0/N0/P0, lines 167-170); the uncoupled and coupled pole census against pinned reference roots (lines 189-194); the four pure-Q R_{Q,*} figures (lines 204-208); the eta=0.1/0.3 thresholds (lines 219-220); and the full lambda_W scan rows incl. monotone-increasing P0 / monotone-decreasing upper-wall R_{Q,*} (lines 259-264). The Mathematica `.wl` mirrors the same six deliverables via an independent route (see Independent-derivation check).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) overlap kappa=2sqrt(2)/pi | py L55-56 / wl L113-114 (symbolic integral) | match |
| (2) quartic F(y), D=F/((varpi^2-omega^2)Delta), deg 4 | py L92-101 / wl L138-163 | match |
| (3) residue/linewidth R_{Q,*}=27 c_s^5/(a^5 omega_*^5 N) | py L113-117 / wl L167-172 | match |
| (4) low-loss survival inequality / threshold | py L119-124 / wl L174-183 | match |
| (5) numerical pole census (slice 6.2) | py L167-208 / wl L210-259 | match |
| (6) static/dynamic tension under lambda_W scan (sec 8) | py L229-264 / wl L272-306 | match |

`paper_alignment: aligned`. Each notes-side numeric (overlap, static data, pole census, R_Q figures, thresholds, scan table incl. the corrected `145.483858657863`) maps to a script-side check and a saved-output line. The published card/appendix carry only symbolic forms and are unaffected.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 56 | `simplify(kappa - 2sqrt2/pi)==0` | claim 1 | yes |
| A2 | sympy | 100 | `cancel(D - F/((varpi^2-omega^2)Delta))==0` | claim 2 | yes |
| A3 | sympy | 101 | `Poly(F,y).degree()==4` | claim 2 | yes |
| A4 | sympy | 117 | `simplify(RQ_star - 27 c_s^5/(a^5 omega_*^5 N_star))==0` | claim 3 | yes |
| A5 | sympy | 167-170 | `assert_close` Delta0/D0/N0/P0 | claim 5 | yes |
| A6 | sympy | 189-194 | `assert_close` uncoupled+coupled roots | claim 5 | yes |
| A7 | sympy | 204-208 | `assert_close` four R_{Q,*} | claim 5 | yes |
| A8 | sympy | 219-220 | `assert_close` thresholds eta=0.1/0.3 | claim 4 | yes |
| A9 | sympy | 259-261 | `assert_close` full lambda_W scan rows | claim 6 | yes |
| A10 | sympy | 263-264 | monotone inc P0 / monotone dec R_{Q,*} | claim 6 | yes |
| M1 | mathematica | 114 | `expectZero[kapIntegral - 2sqrt2/Pi]` | claim 1 | yes |
| M2 | mathematica | 153-163 | native-numerator reconstructs D; equals F; deg 4 | claim 2 | yes |
| M3 | mathematica | 172,182,183 | residue cancellation; peak; Solve-derived threshold | claims 3,4 | yes |
| M4 | mathematica | 215-218 | `expectClose` Delta0/D0/N0/P0 | claim 5 | yes |
| M5 | mathematica | 236-259 | `expectVectorClose` roots + R_{Q,*} | claim 5 | yes |
| M6 | mathematica | 269-306 | thresholds + full scan + monotone | claims 4,6 | yes |

Every assertion traces to a paper deliverable; none are orphaned. None are tautological: A2/M2 cross-check two *independently-constructed* objects (D built from the K_B - Q/Delta bundle vs. F built as the product/native-numerator); A4/M3 derive R_{Q,*} from `(1/D0')/(Gamma5 omega^5 N/D0')` so the D0' cancellation is exercised, not assumed; A10/M6 are genuine sign tests on numerically-distinct values.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is a **genuinely independent route**, not a transliteration. Justification by three corresponding sections:

1. **Quartic (M2).** The `.py` writes the product form by hand (`F_y = (((K-M y)(varpi^2-y)-C^2)((Omega_U^2-y)(Omega_W^2-y)-R^2) - (varpi^2-y)(...))`, py L92-95) and asserts `D - F/denominator == 0`. The `.wl` instead extracts the quartic from D's *own* cleared denominator — `nativeNumerator = Numerator[Together[dBranch]]`, then `quarticFromD = Collect[Expand[PowerExpand[nativeNumerator /. omega -> Sqrt[y]]], y, FullSimplify]` (wl L138-143) — and only *then* checks it against the product form `fProductY` (wl L157) AND reconstructs `Together[D]` from it (wl L153-156) AND extracts `CoefficientList`/`Exponent` (wl L147-148, L162-163). The `.py` never forms the native numerator nor the coefficient list; the `.wl` derives the polynomial from the bundle rather than asserting a pre-written form. Different route.

2. **Survival threshold (M3).** The `.py` writes the threshold as a direct closed form: `survival_threshold = simplify(2*DeltaV_req*(1+eta^2)*x^6/eta)` (py L122-124). The `.wl` *solves* for it: `thresholdFromPeak = stripConditional[rqRatio /. First[Solve[stagePeak == deltaVreq, rqRatio]]]` (wl L179) and then checks equality with the expected form (wl L183). Inversion vs. direct assertion — different route.

3. **Root finding (M5/M6).** The `.py` uses `sp.nroots(expand(poly), n=30)` (py L9); the `.wl` uses `NSolve[poly==0, y, WorkingPrecision->50]` (wl L67). Different solvers; the `.wl` residuals in its output (e.g. M5 ~1e-17, M6 ~1e-14) are independently computed, not echoed.

Variable naming also diverges throughout (`lam_B`/`K`/`a`/`omega_star` in py vs `lamB`/`kWall`/`aScale`/`omegaStar` in wl), and the `.wl` is organized into M1-M6 banner blocks vs the `.py`'s flat flow. No shared-op signature indicating a port. INDEPENDENCE CALL: **independent**.

## Engine cross-check

Both engines agree at the level claimed:
- overlap: both `2 Sqrt[2]/Pi` (py out L4 / wl out L11).
- static slice: Delta0=1.90933940817883, D0=2.76355510933127, N0=0.0501661980249591, P0=0.0181527764203328 (py out L28-31 / wl out L45-52, residuals ~1e-17).
- coupled roots: {0.9382727..., 1.3914108..., 1.7204537..., 2.0453994...} (py out L35 / wl out L59, residuals ~1e-17).
- four R_{Q,*}: 18.7069..., 0.380740..., 16.0250..., 32.0025... (py out L38-41 / wl out L60, residuals ~1e-13).
- lambda_W scan upper-wall R_{Q,*}: 145.4838586578641 / 32.0025... / 13.6885... / 7.5809... / 4.8273... — identical to 13+ figures in both engines (py out L48-52 / wl out L78-82). Both report monotone P0 up / R_{Q,*} down (py out L54-55 / wl out L84-87).

No engine disagreement.

## Verdict justification

`clean`. I read the card, the notes, and the Part-VII appendix rows, then attacked the scripts. Attacks tried and failed: (a) checked A2/M2 for tautology — they cross two independently-built objects, not a define-then-assert; (b) checked A4/M3 for hidden cancellation assumption — the D0' cancellation is genuinely exercised by constructing both `|A_*|=1/D0'` and `gamma_*=Gamma5 omega^5 N/D0'` and dividing; (c) checked the `assume positive=True` domains (py L64-68, wl L89-105) against the physical setup — all couplings/frequencies/scales are positive by construction, consistent with the card's admissibility conditions; (d) checked symmetry/parity of the overlap integral — `u_0 f_0` on [0,L] integrates to the stated `2sqrt(2)/pi` with no sign trap; (e) checked monotone tests A10/M6 for triviality — values are numerically distinct and the direction is the load-bearing claim. The script's claim matches the paper's claim across all six deliverables.

## Value Reconciliation (pass-2 augmentation)

Authoritative record: script source + committed `.txt` outputs (both fresh; see mtime note). The natural carrier for numeric slice values is the `.md` notes (the published card/appendix carry only symbolic forms, consistent with `\StatusNumerical{}`).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| kappa = 2sqrt(2)/pi | py L56 / wl L114; py out L4, wl out L11 | notes L83-84; appendix eq L531 | MATCH |
| F(y) quartic + D=F/((varpi^2-omega^2)Delta) | py L92-101; py out L13-14 | notes L143-154; appendix eq L545-553 | MATCH |
| R_{Q,*}=27 c_s^5/(a^5 omega_*^5 N) | py L117; py out L17 | notes L210-212,227-229; appendix eq L566-570 | MATCH |
| threshold = 2 DeltaV_req (1+eta^2) x^6/eta | py L124; py out L21 | notes L263-269 | MATCH |
| C=0.450158158078553 | py out L24 | notes L298 | MATCH |
| G_U=0.3 | py out L25 | notes L300 | MATCH |
| G_W=0.360126526462842 | py out L26 | notes L301 | MATCH |
| R=0.225079079039277 | py out L27 | notes L302 | MATCH |
| Delta0=1.90933940817883 | py out L28 / wl out L45 | notes L309 | MATCH |
| D0=2.76355510933127 | py out L29 / wl out L47 | notes L310 | MATCH |
| N0=0.0501661980249591 | py out L30 / wl out L49 | notes L313 | MATCH |
| P0=0.0181527764203328 | py out L31 / wl out L51 | notes L314 | MATCH |
| uncoupled wall roots 1.68143182591478, 2.04274007519334 | py out L33 / wl out L57 | notes L323-324 | MATCH |
| uncoupled internal roots 0.974601723746314, 1.41779810977117 | py out L34 / wl out L58 | notes L326-327 | MATCH |
| coupled poles 0.938272741746754 / 1.39141087653804 / 1.72045371048003 / 2.04539948783659 | py out L35 / wl out L59 | notes L339-342 | MATCH |
| R_{Q,*} census 18.7069287828307 / 0.380740659074003 / 16.0250330226177 / 32.0025481088465 | py out L38-41 / wl out L60 | notes L339-342 | MATCH |
| thresholds eta=0.1: 21.8545662963584; eta=0.3: 7.86187368416853 | py out L44-45 / wl out L73,75 | notes L362,365 | MATCH |
| benchmark V_known(1)=1.181909222592, DeltaV_req(1)=1.081909222592 | py L213-215 | notes L351,355-356 | MATCH |
| lambda_W scan P0 column 0.00594740531769 / 0.01815277642033 / 0.03800016314041 / 0.06717078268091 / 0.10847330811048 | py out L48-52 / wl out L78-82 | notes L395-399 | MATCH |
| lambda_W scan D0 column 2.82723442158450 / 2.76355510933127 / 2.66591349720966 / 2.53430958521967 / 2.36874337336129 | py out L48-52 / wl out L78-82 | notes L395-399 | MATCH |
| lambda_W scan upper omega_* 2.04402272302752 / 2.04539948783659 / 2.04793277506821 / 2.05190668889211 / 2.05778339035510 | py out L48-52 / wl out L78-82 | notes L395-399 | MATCH |
| **lambda_W=0.2 upper-wall R_{Q,*} = 145.483858657863** | py out L48 / wl out L78,83 | notes L395 | **MATCH** (corrected value confirmed) |
| lambda_W scan upper R_{Q,*} 32.0025481088465 / 13.6885356356808 / 7.58097126746582 / 4.82738925564702 | py out L49-52 / wl out L79-82 | notes L396-399 | MATCH |

**reconciliation: complete; 24 deliverable rows checked, 0 misaligned.**

INTERNAL scaffolding (no finding): pass/fail flags and `PASS:` prints; residual/diff magnitudes (~1e-17 to 1e-13); tolerances (1e-12, 1e-10); the intermediate `Q(omega)`, `K_B(omega)`, `Delta(omega)` symbolic invariant prints (py out L7-10); CoefficientList of the quartic (wl out L17); the native-numerator reconstruction residuals.

### R_Q reconfirmation (explicit)

- SymPy emits `145.4838586578641` (assertion py L253; output py-out L48). NOT 213.x.
- Mathematica emits `145.48385865786412765610...` (assertion wl L291; output wl-out L78 and L83 `lambda_W=0.2 upper-wall R_Q,* = 145.483...`). NOT 213.x.
- Notes line 395 reads `145.483858657863` — the corrected value. No surviving `213.483858657863` in the stage-222 `.md` or in the card/appendix.
- The only files containing the string `213.483858657863` are the process trackers (`STAGE_VERIFICATION_COVERAGE.md`, `CHECKPOINT_CONSTANT_PROVENANCE.md`, `PAPER_CLEANUP_TRACKER.md`, `MATHEMATICA_MIRROR_POLICY.md`), each in historical-record prose describing the typo correction itself — these are out-of-scope process docs, not stage content, and the correct `→145...` resolution is recorded there. No action needed.

### Stale-output note

Not stale. SymPy output (2026-06-02 16:15:50) is newer than `.py` (2026-05-11 11:58). Mathematica output (2026-06-02 16:15:50) is newer than `.wl` (2026-06-02 14:30). Reconciliation is based on script source + committed outputs (no execution).

## Self-test notes

Variable-independence: no `sp.diff`/`D[...]` derivatives in either script, so the zero-derivative trap does not apply. Symmetry/parity: the only integral is the bounded overlap `int_0^L u_0 f_0` over a finite interval — value `2sqrt(2)/pi` is positive and well-posed, no symmetric-domain vanishing trap. Trivial-case: A2/M2 are non-trivial (two independently-built objects equated); A4/M3 exercise the D0' cancellation by explicit construction; A10/M6 monotone tests run on numerically-distinct values. No directive written (zero findings).
