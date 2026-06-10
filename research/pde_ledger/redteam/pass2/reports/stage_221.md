---
unit_id: 221
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
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 221 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_221.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows: line 54 status table; lines 461-513 theorem block + claim status; line 1431 stage input)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.txt`

## What the paper claims

Stage 221 is the Part VII resonance/linewidth checkpoint (anchor MTDC-T11.1). It expands the Stage-220 dynamic one-port mixed bundle near a simple conservative pole `F_0(omega_*)=0`, `F_0'(omega_*)!=0` (and the wall specialization `D_0(omega_*)=0`, `Delta_0(omega_*)!=0`), reduces the reduced susceptibility to a universal Breit–Wigner form, and derives an exact dispersive/absorptive no-free-lunch tradeoff. The card's `\stagefield{Output}` is verbatim: "Dispersive/no-free-lunch theorem: $|\Re\chi_*|/|\Im\chi_*|=|\delta|/\gamma_*$, so resonant conservative gain is inseparable from absorptive loading, with maximal balanced leverage at equal conservative and absorptive magnitudes." The notes (Purpose, sections 2-9) enumerate the full deliverable set: (1) the simple-pole Breit–Wigner normal form `chi ~ A_*/(delta - i gamma_*)` with `A_* = N_s/F_0'`, `gamma_* = Gamma_* Z_*/|F_0'|`; (2) the Stage-220 derivative identity `dD_Pi/dPi = -N`, `N=(A G_W + R G_U)^2/Delta_Pi^2`; (3) the wall-like specialization `chi_qq ~ (1/D_0') / (delta - i gamma_wall)`; (4) the exact Re/Im line shapes and the ratio `r`; (5) the conservative-max identity (`r=1`, `max|Re chi| = |A_*|/(2 gamma_*)`); (6) the low-loss factorization identity and bound `|A_*|/gamma_* * eta/(1+eta^2)`; (7) the barrier/absorbed-power ratio `gamma_*/|delta|`; (8) the quality-factor detuning bound `>= 1/(2 Q_* eta)`; (9) the linear survival-window inequality in residue-to-linewidth form. The appendix theorem block (lines 461-513) carries the same equations and marks Stages 220-221 ExactClosure for the algebra with the actual-branch realization Open.

## What the script claims to verify

Both scripts verify the seven theorem-block identities plus the survival window, all symbolically (no numeric proof path). SymPy: section I builds the local denominator `F_lin = Fprime*delta - Pi*Z_star` and checks the Breit–Wigner collapse (`chi_lin == chi_bw`, passive form); section II posits the perfect-square `N` and checks `diff(D_Pi,Pi)+N==0`; section III checks the wall-like specialization on the `D0prime>0` slice; section IV does the real/imag split via `expand_complex(...).as_real_imag()` and the two factorization identities; section V the barrier/power ratio; section VI the quality-factor bound and the two survival-window connecting identities. Mathematica mirrors this theorem path but with native independent operations (see independence check). A probe-only numeric slice (`Aabs=7, gamma=2, r=3, eta=1/4, ...`) runs in both engines and is explicitly flagged non-proof (notes lines 79, 156-158). All checks are substantive (`assert simplify(...)==0` / `expectZero[..., ...]` with `Exit[1]` on nonzero) — no print-only deliverables.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) Breit–Wigner normal form, A_*, gamma_* | py L18-33 / wl L62-85 (+ Residue) | match |
| (2) dD_Pi/dPi = -N, N perfect square | py L43-54 / wl L87-102 | match |
| (3) wall-like specialization, gamma_wall | py L56-66 / wl L104-117 | match |
| (4) Re/Im line shapes, |Re|/|Im| = r | py L71-84 / wl L119-134 | match |
| (5) max identity, r=1, |A_*|/(2 gamma_*) | py L89-91,99 / wl L136-140 | match |
| (6) low-loss factorization + bound | py L93-96,100 / wl L141-144,167-168 | match |
| (7) barrier/absorbed-power ratio 1/r | py L104-114 / wl L146-156 | match |
| (8) quality-factor detuning bound | py L119-131 / wl L158-165 | match |
| (9) linear survival window (residue/linewidth) | py L136-151 / wl L167-183 | match |

`paper_alignment: aligned`. Every Output/deliverable item in the card, notes, and appendix theorem block has a substantive script-side check in BOTH engines, with the same symbolic forms.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 24 | `simplify(chi_lin - chi_bw)==0` | (1) BW collapse | yes |
| A2 | sympy | 28 | `simplify(chi_passive - chi_passive_expected)==0` | (1) passive form | yes |
| A3 | sympy | 52 | `simplify(diff(D_Pi,Pi)+Nfun)==0` | (2) derivative id | yes |
| A4 | sympy | 62 | `simplify(chiqq_lin.subs(...) - chiqq_expected)==0` | (3) wall | yes |
| A5 | sympy | 77-79 | `re_r==re_exp`, `im_r==im_exp`, `re_r/im_r==r` | (4) line shape + ratio | yes |
| A6 | sympy | 91 | `max_identity + (r-1)^2/(2(1+r^2))==0` | (5) max id | yes |
| A7 | sympy | 96 | low-loss factorization == 0 | (6) low-loss id | yes |
| A8 | sympy | 109 | `ratio_barrier - 1/r==0` | (7) power ratio | yes |
| A9 | sympy | 124,127 | detuning q-form + low-loss boundary | (8) Q bound | yes |
| A10 | sympy | 143,144-147 | survival-window connecting identities | (9) survival window | yes |
| B1 | mathematica | 72,77 | `chiLin-chiBW`, passive == 0 | (1) | yes |
| B2 | mathematica | 78-81 | `Residue[chiPassive,{delta,I gammaStar}] - Astar == 0` | (1) residue (extra, independent) | yes |
| B3 | mathematica | 98-102 | `D[QPi/DeltaPi,portPi]` perfect square + `D[DPi,portPi]+N==0` | (2) derivative id (native D, derived N) | yes |
| B4 | mathematica | 111-114 | wall specialization == 0 | (3) | yes |
| B5 | mathematica | 130-134 | generic + substituted Re/Im + ratio | (4) | yes |
| B6 | mathematica | 137-144 | max + low-loss identities | (5),(6) | yes |
| B7 | mathematica | 152 | barrier ratio - 1/r == 0 | (7) | yes |
| B8 | mathematica | 164-165 | detuning + low-loss boundary | (8) | yes |
| B9 | mathematica | 173-180 | survival-window connecting identities | (9) | yes |

All rows are non-tautological: each LHS is built from one expression chain (the line-shape / denominator derivation) and compared against an independently written closed-form deliverable, so a wrong derivation would produce a nonzero residual.

## Findings

None. (See verdict justification.)

## Independent-derivation check (Mathematica)

**INDEPENDENCE CALL: independent.** The re-authored `.wl` derives every load-bearing quantity by a genuinely different operation than the `.py`, and adds independent checks the `.py` lacks. Three load-bearing quantities, each examined:

1. **The residue `A_*` (line-shape residue).** `.py` only rearranges algebraically: `chi_lin = together(Num_star/F_lin)` and asserts it equals the rewritten `chi_bw` (L19-24) — pure `together`. `.wl` instead extracts the residue with the native pole operation `Residue[chiPassive, {delta, I gammaStar}] - Astar` (L78-81), which the `.py` never does. Different operation: SymPy `together`-equality vs Mathematica `Residue` at the explicit pole `delta = I gamma_*`. Independent.

2. **The linewidth / derivative identity `dD_Pi/dPi = -N`.** This is the sanctioned native-`D` route. `.py` (L50-52) POSITS the closed form `Nfun = (A*G_W + R*G_U)^2/Delta_Pi^2` and only checks `diff(D_Pi,Pi)+Nfun==0`. `.wl` (L96-102) DERIVES `NfunDerived = Together[D[QPi/DeltaPi, portPi]]` from `Q_Pi/Delta_Pi` natively, THEN checks (a) that the derived `N` equals the perfect square, AND (b) `D[DPi,portPi]+NfunDerived==0`. So the perfect-square form is an *output* on the Mathematica side and an *input* on the SymPy side — strictly more independent, not a shared posited form. The identity itself holds: `D_Pi = K_B - Q_Pi/Delta_Pi` with `K_B` independent of `Pi`, `dW/dPi=-1`, giving numerator `A*Q_Pi - G_U^2*Delta_Pi = (A G_W + R G_U)^2`, a perfect square; hence `dD_Pi/dPi = -(A G_W + R G_U)^2/Delta_Pi^2 = -N`. Confirmed by hand. Independent, identity sound.

3. **The real/imaginary line-shape split (the line-shape "derivative" of the Breit–Wigner form).** `.py` (L72-73) substitutes `delta = r*gamma` FIRST, then splits with `sp.expand_complex(chi_r).as_real_imag()`. `.wl` (L121-125) splits the GENERIC `delta`-dependent form FIRST with `ComplexExpand[Re[...]]`/`ComplexExpand[Im[...]]`, checks the uncollapsed Breit–Wigner forms `Aabs*delta/(delta^2+gamma^2)` and `Aabs*gamma/(delta^2+gamma^2)` (L130-131, an extra check with no `.py` counterpart), THEN substitutes `delta -> r gamma`. Different mechanism (expand_complex/as_real_imag vs ComplexExpand) AND different order (substitute-then-split vs split-then-substitute). Independent.

Conclusion: this is NOT the precedent failure mode where both engines computed the load-bearing quantity by the same operation despite a re-author. The residue is derived (Mathematica `Residue` vs SymPy algebra), `N` is derived rather than posited on the Mathematica side, and the line shape is split by a different operation in a different order. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce identical symbolic deliverables and identical probe numerics. Spot comparison of saved outputs:
- `A_* = Num_star/Fprime`, `gamma_* = Gamma_out*Z_star/Fprime` — py out L5-6 vs wl out L16-17. Agree.
- `Re chi = Aabs*r/(gamma*(1+r^2))`, `Im chi = Aabs/(gamma*(1+r^2))`, `|Re|/|Im| = r` — py out L16-18 vs wl checks IV (all PASS). Agree.
- `U_disp = -Aabs*Sfam^2*r/(2 gamma(1+r^2))`, `P_abs = Aabs*Sfam^2*omega_drive/(2 gamma(1+r^2))`, ratio `1/r` — py out L25-27 vs wl out L58-60. Agree.
- Probe slice: `Re chi=1.05`, `Im chi=0.35`, `|Re|/|Im|=3`, `|U_disp|=13.125`, `P_abs=48.125`, `P_abs/(omega|U|)=0.3333`, `low-loss |U|_max=10.294...`, `required |A|/gamma=0.204` — identical in py out L38-45 and wl out L79-86.

Note: the committed Mathematica `.txt` is STALE. Its banner prints "STAGE 204" (wl out L3), whereas the current `.wl` source prints "STAGE 221" (`.wl` L35). The `.wl` mtime (Jun 3 15:59) is newer than its `.txt` (Jun 2 16:15), so the saved transcript predates the re-author. The CONTENT of the checks is unchanged (all the same PASS lines and identical symbolic forms appear), so the staleness is the deferred-cosmetic banner drift the prompt already sanctioned — not a content disagreement. Flagged informationally; the orchestrator's independent re-run will refresh it.

## Verdict justification

`clean`. Attacks tried that failed: (a) checked every assertion for tautology — each compares a derivation chain against an independently written closed form, so none is `x==x`; (b) checked the line-shape `.wl` split for shared-operation transliteration — it uses a different mechanism and order, plus a `Residue` extraction and a native-`D` derived `N` that the `.py` lacks, so the re-author is genuinely independent on all three load-bearing quantities; (c) re-derived the `dD_Pi/dPi=-N` perfect-square identity by hand and confirmed it; (d) checked symbol domains — both engines fix the positive sign-fixed slice (`Fprime>0`, `D0prime>0`, etc.) matching notes lines 18-24, so `|F_0'|=F_0'` is legitimate and the absolute-value forms restore at the end as the notes state; (e) confirmed F1 (survival round-trips) and F2 (deliverable #9) fixes hold — both round-trips are non-tautological connecting identities present as real `assert`/`expectZero` checks in BOTH engines, not print-only. The only blemish is the stale Mathematica `.txt` banner ("STAGE 204"), which is the sanctioned deferred-cosmetic drift, not a content defect. All nine deliverables reconcile to the card, notes, and appendix theorem block. I read the card, the full notes, and the appendix theorem rows; the script's verified claim matches the paper's claim exactly.

## Value Reconciliation (pass-2 augmentation)

Enumerated every symbolic deliverable and labeled numeric the scripts emit (from `.py`/`.wl` source + committed outputs), excluding verification scaffolding (residual==0 PASS flags, the probe-only numeric slice).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `chi ~ A_*/(delta - i gamma_*)`, `A_*=Num_star/Fprime`, `gamma_*=Gamma_out*Z_star/Fprime` | py L31-33 / out L4-6; wl L83-85 / out L15-17 | notes L179-187 (boxed); appendix `eq:app-part07-breit-wigner` L471-475 | MATCH |
| `dD_Pi/dPi = -N`, `N=(A G_W+R G_U)^2/Delta_Pi^2` | py L53-54 / out L8-9; wl L98-102 / out L22-25 | notes L207-211 (boxed) | MATCH |
| `chi_qq ~ (1/D0prime)/(delta - i gamma_wall)`, `gamma_wall=Gamma_out*N_star/D0prime` | py L65-66 / out L12-13; wl L116-117 / out L32-33 | notes L233-242 | MATCH |
| `Re chi=Aabs*r/(gamma(1+r^2))`, `Im chi=Aabs/(gamma(1+r^2))`, `|Re|/|Im|=r` | py L82-84 / out L16-18; wl L130-134 | notes L282-293 (boxed); appendix `eq:app-part07-re-im-ratio` L480-484 | MATCH |
| `f(r)=r/(1+r^2)`; `f-1/2 = -(r-1)^2/(2(1+r^2))`; `f-eta/(1+eta^2) = -(r-eta)(eta r-1)/((1+r^2)(1+eta^2))` | py L98-100 / out L20-22; wl L136-144 | notes L303-309 (max), L352-355 (low-loss) | MATCH |
| `U_disp=-Aabs*Sfam^2*r/(2 gamma(1+r^2))`, `P_abs=Aabs*Sfam^2*omega_drive/(2 gamma(1+r^2))`, ratio `1/r` | py L112-114 / out L25-27; wl L154-156 / out L58-60 | notes L384-409 (boxed ratio) | MATCH |
| `|omega-omega_*|/omega_* = r/(2 Qfac)`; low-loss boundary `1/(2 Qfac eta)` | py L130-131 / out L30-31; wl L164-165 | notes L427-429, L436-440 (boxed) | MATCH |
| `|U_disp|_max(low-loss) = Aabs*Sfam^2*eta/(2 gamma(1+eta^2))` | py L150 / out L34; wl L182 | notes L450-453 | MATCH |
| required `|A|/gamma >= 2*DeltaV_req*(1+eta^2)/(Sfam^2*eta)` | py L151 / out L35; wl L183 | notes L467-473 (boxed); appendix `eq:app-part07-low-loss-bound` L494-498 + `linear-survival-window` L501-506 | MATCH |

INTERNAL (scaffolding, no doc carrier expected, no finding): probe-only numeric slice values `Aabs=7, gamma=2, r=3, eta=1/4, Sfam=5, omega_drive=11, Qfac=40, omega_star=80, DeltaV_req=3/5` and the derived probe numerics (`Re chi=1.05`, `Im chi=0.35`, `|U_disp|=13.125`, `P_abs=48.125`, `P_abs/(omega|U|)=0.3333`, `low-loss |U|_max=10.294`, `required |A|/gamma=0.204`) — explicitly flagged probe-only/non-proof in notes lines 79, 156-158; all residual==0 PASS flags.

reconciliation: complete; 9 deliverable values checked, 0 misaligned.

## Self-test notes

Trap 1 (variable independence of derivatives): the load-bearing `D[QPi/DeltaPi, portPi]` and `sp.diff(D_Pi, Pi)` are taken w.r.t. `Pi`/`portPi`, on which both `Q_Pi` and `Delta_Pi` genuinely depend (`W = OmegaW^2 - omega^2 - portPi`); the derivative is non-trivially nonzero and yields the perfect square — not an identically-zero derivative. Trap 2 (symmetry/parity): no unbounded-domain integrals in this stage. Trap 3 (trivial-case): the survival round-trips (py L143-147 / wl L173-180) connect two independently written expression chains (line-shape `re_expected` vs the boxed deliverable; `residue_requirement` vs the saturated survival LHS) and are non-tautological; substituting the boundary `r=1/eta` and `Aabs->residue_requirement*gamma` reduces the residual to literal 0 only because the derivations are consistent, exactly the intended check. Trap 5 (paper round-trip): no script-side change prescribed (zero findings), so no new paper_misalignment introduced. The only residual signal is the stale Mathematica `.txt` banner ("STAGE 204" vs source "STAGE 221"), the sanctioned deferred cosmetic — informational, not blocking.
