---
unit_id: 086
batch: III.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T23:20:35Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage086_family1_loading_ratio_window.md]
  paper_appendix: present
---

# Audit unit 086 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_086.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage086_family1_loading_ratio_window.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 150)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage086_family1_loading_ratio_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage086_family1_loading_ratio_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage086_family1_loading_ratio_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage086_family1_loading_ratio_window_mathematica_audit.txt`

## What the paper claims

Stage 086 converts the Family-1 support/source verdict into a pure window in the loading ratio `rho_alpha = alpha_req/alpha_mix`. The mechanism is the Stage-081 product inversion `Pi_tr = C_mix Q(zeta;eps_blk)` with `Q(zeta;eps_blk) = [1 + (1 - 2 eps_blk) zeta]/[1 - eps_blk zeta]`, combined with Stage-085's reduction `Pi_tr/C_mix = rho_alpha`. The card's `\stagefield{Output}` is "The pure Family--1 loading-ratio window \eqref{eq:app-stage086-success-window}--\eqref{eq:app-stage086-fail-window}", with the two boxed equations stating, at `lambda_mu=1`, `eps_blk=0`: `rho_alpha <= 3.46622291347846` => guaranteed success; `rho_alpha >= 3.46752913273870` => guaranteed failure; `rho_alpha < 3.46752922945601` => absolute constructive ceiling. The notes add detail beyond the .tex: the unblocked reduction `Q(zeta;0) = 1 + zeta`; the conservative lower envelope `rho_suff^(J) = 3.44257571477179`; the input zeta values (`zeta_suff^(chi)=2.46622291347846`, `zeta_fail^(chi)=2.46752913273870`, `zeta_suff^(J)=2.44257571477179`, `zeta_max^(F1)=2.46752922945601`, carried from Stages 080/081); the blocking cap `eps_blk < 1/zeta_max ≈ 0.405263689711371`; and the monotonicity `d rho_max/d eps_blk = zeta_max(zeta_max-1)/(1 - eps_blk zeta_max)^2 > 0`. The appendix row (150) is a one-line status: "Pure loading-ratio success/failure window", \StatusExactClosure{}.

## What the script claims to verify

Both scripts define the Q-map verbatim, then: (1) assert the exact unblocked reduction `Q(zeta;0) - (1+zeta) == 0` (symbolic, non-tautological — the actual eps->0 collapse); (2) take the four Family-1 zeta values as carry-forward inputs and evaluate the rho thresholds via `rho = Q(zeta, eps=0) = 1 + zeta`, comparing each against the paper's literal rho values to tolerance; (3) print/assert the derivative form `dQ/dzeta = (1 - eps)/(1 - eps zeta)^2`; (4) compute the blocking cap `1/zeta_max` and verify it against the notes value; (5) print/assert the closed form `d rho_max/d eps = zeta_max(zeta_max-1)/(1 - eps_blk zeta_max)^2`. The Mathematica script additionally asserts the `dQ` and `d rho_max/d eps` closed forms symbolically (against independently supplied `dQExpected`/`dQMaxExpected`, using the symbol `zetaMax`, not the numeric), and asserts that each input zeta literal matches the paper value — checks the SymPy script only prints or only consumes as input.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Q(zeta;eps_blk)` defining map | sympy L30 / wl L39 define it verbatim | match |
| Unblocked reduction `Q(zeta;0)=1+zeta` | sympy L35 `expect_zero`; wl L46 `expectZero` | match |
| `rho_suff^(chi) = 3.46622291347846` (success) | sympy L53; wl L68 numeric check | match |
| `rho_fail^(chi) = 3.46752913273870` (failure) | sympy L54; wl L69 | match |
| `rho_max^(F1) = 3.46752922945601` (ceiling) | sympy L56; wl L71 | match |
| `rho_suff^(J) = 3.44257571477179` (notes lower envelope) | sympy L55; wl L70 | match |
| blocking cap `1/zeta_max ≈ 0.405263689711371` (notes) | sympy L59 print; wl L77 `expectApprox` | match |
| `d rho_max/d eps_blk = zeta_max(zeta_max-1)/(1-eps_blk zeta_max)^2 > 0` (notes monotonicity) | sympy L64 print only; wl L74-76 symbolic `expectZero` | match |
| `dQ/dzeta = (1-eps)/(1-eps zeta)^2` (notes) | sympy L32 print only; wl L40-45 symbolic `expectZero` | match |

`paper_alignment: aligned`. Every boxed window value and every notes-stated auxiliary value has a corresponding script-side check; no script-side assertion is orphaned from the paper/notes; no value differs.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 35 | `expect_zero(Q(zeta;0) - (1+zeta))` | unblocked reduction | yes |
| A2 | sympy | 53 | `expect_close(rho_suff_chi, 3.46622291347846)` | success threshold | partial (weak corroboration; see note) |
| A3 | sympy | 54 | `expect_close(rho_fail_chi, 3.46752913273870)` | failure threshold | partial |
| A4 | sympy | 55 | `expect_close(rho_suff_J, 3.44257571477179)` | conservative envelope | partial |
| A5 | sympy | 56 | `expect_close(rho_max, 3.46752922945601)` | ceiling | partial |
| A6 | mathematica | 45 | `expectZero(dQ - dQExpected)` | dQ/dzeta closed form | yes |
| A7 | mathematica | 46 | `expectZero(Q(.,0) - (1+zeta))` | unblocked reduction | yes |
| A8 | mathematica | 53-56 | `expectApprox(zeta_* vs paper)` | input zeta carry-forward | yes (validates inputs) |
| A9 | mathematica | 68-71 | `expectApprox(rho_* numeric)` | window thresholds | partial (round-trip of 1+zeta; see note) |
| A10 | mathematica | 76 | `expectZero(dQMax - dQMaxExpected)` | d rho_max/d eps monotonicity form | yes |
| A11 | mathematica | 77 | `expectApprox(1/zeta_max vs 0.40526...)` | blocking cap | yes |

Note on A2-A5/A9 ("partial"): each rho check evaluates `rho = Q(zeta, eps=0)`, which by the (independently asserted) identity A1/A7 equals `1 + zeta`, then compares against the paper rho literal. Since the input zeta literal and the target rho literal both originate in the same notes and differ by exactly 1, these numeric checks are weak corroboration rather than independent re-derivations of the thresholds — they confirm the typed zeta inputs round-trip through the +1 reduction to the documented rho. They are NOT tautological (they do exercise the `eps=0` collapse on concrete numbers and would catch a mistyped literal), and the load-bearing physics content lives in A1/A7 (the symbolic reduction) plus A8 (the zeta-vs-paper validation). No finding: the substantive identities are genuinely exercised, and the threshold values themselves are derived upstream (Stage 080 solves them from transcendental equations); Stage 086 legitimately consumes them as carry-forward inputs.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is not a transliteration of the `.py`. Both must define the same `Q` map and the same `Q(zeta;0)=1+zeta` reduction (forced — these are the paper's defining formulas), but the Mathematica script diverges substantively from there: it supplies its OWN closed forms `dQExpected = (1 - epsBlk)/(1 - epsBlk*zeta)^2` (wl L41) and `dQMaxExpected = zetaMax*(zetaMax - 1)/(1 - epsBlk*zetaMax)^2` (wl L75) and ASSERTS the SymPy-derived derivatives equal them symbolically (wl L45, L76), whereas SymPy merely PRINTS `dQ/dzeta` and `d rho_max/d eps` (py L32, L64) with no assertion. The Mathematica script also asserts each input zeta literal against the paper value (wl L53-56), which SymPy does not. The `d rho_max/d eps` check uses the symbol `zetaMax` (general, assumed `>1`), not the numeric — a fully symbolic, independent verification of the monotonicity form. This is an independent second engine, not an echo.

## Engine cross-check

Both engines agree at full precision. SymPy output line 8: `rho_suff^(chi) = 3.46622291347846012143918414949`; Mathematica numeric target (wl L68) = `3.46622291347846012143918414949` — identical to all 30 digits, and the same for fail/J/max (py output 9-11 vs wl L69-71). The unblocked reduction prints `0` in both (py output L7, wl output L9). The `dQ/dzeta` printed forms agree: py `(1 - eps)/(eps**2*zeta**2 - 2*eps*zeta + 1)` = wl `(1 - epsBlk)/(-1 + epsBlk*zeta)^2`. The eps cap agrees: py `0.4052636897113714997686884` (output L16) = wl target `0.4052636897113714997686884` (wl L77). No residual, sign, or factor disagreement.

## Verdict justification

`clean`. I read the card, the notes, and the appendix row before opening the scripts. The paper's `\stagefield{Output}` — the three-value loading-ratio window plus the notes' auxiliary deliverables (the `J` envelope, the blocking cap, the two derivative closed forms) — is fully and faithfully covered by the script-side checks, with every value matching to full precision across both engines. Attacks tried that failed: (1) tautology — the rho checks looked like candidates, but the load-bearing identities (the symbolic `Q(zeta;0)=1+zeta` reduction and the two symbolic derivative closed forms) are genuine and non-tautological, and the rho threshold values are derived upstream (Stage 080, from transcendental solves) so consuming them as carry-forward inputs is correct rather than self-referential; (2) symbol-assumption — `eps` real-unrestricted in sympy and `0<=epsBlk<1/zetaMax` in Mathematica both leave the asserted identities valid for all `eps != 1/zeta`, hiding nothing; (3) derivative-by-hand recomputation — both `dQ/dzeta = (1-eps)/(1-eps zeta)^2` and `d rho_max/d eps = zeta_max(zeta_max-1)/(1-eps zeta_max)^2` match my own quotient-rule computation, and the signs (`>0` for both, given `eps<1` and `zeta_max>1`) match the notes' monotonicity claims; (4) transliteration — the `.wl` adds independent symbolic assertions the `.py` lacks. Output transcripts are fresh (both `.txt` mtimes 2026-05-27 10:24-10:26, newer than their scripts' 10:17). The self-labels read "STAGE 086" correctly (py L26, wl L32/L80). One stale CROSS-reference is present (py L37 comment "carried from Stages 63-64"; notes attribute these to Stages 080/081 — consistent with the known +17 pre-renumber drift, since 63+17=80, 64+17=81), but per the in-loop numbering policy CROSS-references to other stages are DEFERRED, not findings; I note it here only.

## Value Reconciliation (pass-2 augmentation)

Authoritative record: script source + committed `.txt` outputs (both fresh; not run). Deliverable-level table:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Q(zeta;eps) = (1 + (1-2 eps) zeta)/(1 - eps zeta)` | py L30 / wl L39; py out L5, wl out L5 | .tex eq lines 14-18 (Q via boxed window) ; .md L46-47 | MATCH |
| `Q(zeta;0) = 1 + zeta` | py L35 / wl L46; py out L7 | .md L85 | MATCH |
| `dQ/dzeta = (1-eps)/(1-eps zeta)^2` | py L32 / wl L40-41; py out L6, wl out L6 | .md L133-134 (monotonicity context); implied | MATCH (notes) |
| `rho_suff^(chi) = 3.46622291347846` | py L43,L48,L53 / wl L58,L68; py out L8, wl out L19 | .tex L17 (boxed) ; .md L107,L154 | MATCH |
| `rho_fail^(chi) = 3.46752913273870` | py L44,L49,L54 / wl L59,L69; py out L9, wl out L20 | .tex L23 (boxed) ; .md L109,L156 | MATCH |
| `rho_suff^(J) = 3.44257571477179` | py L45,L50,L55 / wl L60,L70; py out L10, wl out L21 | .md L111 (not in .tex card) | MATCH (notes) |
| `rho_max^(F1) = 3.46752922945601` | py L46,L51,L56 / wl L61,L71; py out L11, wl out L22 | .tex L27 (boxed ceiling) ; .md L113,L160 | MATCH |
| `eps_blk cap = 1/zeta_max ≈ 0.405263689711371` | py L59 / wl L77; py out L16, wl out L33 | .md L127 (not in .tex card) | MATCH (notes) |
| `d rho_max/d eps = zeta_max(zeta_max-1)/(1-eps zeta_max)^2` | py L64 / wl L74-76; py out L18, wl out L31 | .md L133-134 | MATCH |
| input `zeta_suff^(chi) = 2.46622291347846` | py L38 / wl L48,L53 | .md L97 | MATCH |
| input `zeta_fail^(chi) = 2.46752913273870` | py L39 / wl L49,L54 | .md L99 | MATCH |
| input `zeta_suff^(J) = 2.44257571477179` | py L40 / wl L50,L55 | .md L101 | MATCH |
| input `zeta_max^(F1) = 2.46752922945601` | py L41 / wl L51,L56 | .md L103 | MATCH |

INTERNAL (scaffolding, no finding expected): banner strings; `expect_zero`/`expect_close`/`expectZero`/`expectApprox`/`pass`/`fail`/`fmt` helpers; the `expect_close` residual-diff prints (~1e-16/1e-17, near-zero check values); tolerances (`1e-13`, `10^-14`); the unsimplified printed `rho_max^(F1)(eps)` rational (py out L17) and the unsimplified `d rho_max/d eps` cubic-denominator form (py out L18), which are intermediate display forms of the asserted closed form; the FINAL LEDGER prose lines (py out L20-23).

reconciliation: complete; 13 values checked, 0 misaligned

## Self-test notes

Checked: (1) variable independence — `D[qMax, epsBlk]` (wl L74) operates on `qMap /. zeta -> zetaMax`, which genuinely depends on `epsBlk`, so the derivative is non-trivial (not identically zero); `diff(Q, zeta)` likewise depends on `zeta`. (2) Parity/symmetry — no unbounded integrals in this stage; n/a. (3) Trivial-case — substituting `eps=0` into `Q` gives `1+zeta` (the asserted reduction) and the rho values equal `1 + zeta_input` numerically, confirmed against outputs; the two derivative closed forms reproduce my hand quotient-rule computations and have the correct sign (`>0`). (4) Path specs — n/a (no missing-script directive; both engines present). (5) Paper round-trip — no fix prescribed; the clean verdict introduces no new misalignment. No directive written (zero findings).
