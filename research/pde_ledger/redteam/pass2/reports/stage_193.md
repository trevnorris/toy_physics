---
unit_id: 193
batch: V.3
auditor_model: Opus 4.8 (1M context)
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage193_isotropic_grouped_p2_target_surface.md]
  paper_appendix: present
---

# Audit unit 193 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_193.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 117, 1216-1261)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Freezes the isotropic one-pole conservative target surface and proves the scalar/geometry lane enters the grouped carrier only at quadratic anisotropy order." The notes (§2-§6) and the appendix (eqs. `app-part05-conservative-target-surface`, `-one-pole-surface-D`, `-conservative-one-pole-module`, `-geometry-firewall-schur`, `-geometry-firewall-result`) enumerate the deliverables precisely:
1. exact grouped trace/anomaly map with its exact inverse, and the isotropic common-lane collapse `a2=b2=a4=b4=0`;
2. the one-pole identity `Delta_pole = bar_nu4 - 4 bar_nu2^2 = -(3 D2^2 + D0 D4)/D0^2`, equivalent to `D0 D4 + 3 D2^2 = 0`;
3. the one-parameter conservative carrier `Y_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` with `Omega_Q^2 = -D0/(4 D2)`, matching `1 + bar_nu2 omega^2 + bar_nu4 omega^4` through `O(omega^4)`;
4. the scalar/geometry firewall: for a block operator with LINEAR `chi` `l=0 <-> l=2` coupling, the exact Schur complement of the `l=2` block is `D2 I - chi^2 C^T C / D0scalar`, hence `partial_chi D_{2,eff}|_{chi=0} = 0` — no `O(chi)` contamination, contamination begins only at `O(chi^2)`.

This is a `\claimstatus{\StatusExactClosure{}}` stage; checkpoint: no.

## What the script claims to verify

Both scripts verify exactly the four deliverables above. (I) the grouped trace/anomaly map (`xbar, ax, bx`), its inverse reconstructs the lanes, and the common-lane substitution `nu2->nu2->nu2` (resp. `nu4`) drives the anomalies `a,b` to zero. (II) the one-pole identity in `D`-space and its equivalence form. (III) the `3/4 + 1/4(1-omega^2/Omega_Q^2)^{-1}` series matches the one-pole target series through `O(omega^4)`. (IV) the Schur complement of a linearly-`chi`-coupled block operator equals `D2 I - chi^2 C^T C / D0scalar`, the `chi`-linear part vanishes, and `D_eff|_{chi=0} = D2 I`. The script's docstring/ledger prose matches the paper's claim verbatim.

## Paper <-> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Grouped trace/anomaly map + exact inverse | py L57-68 (`grouped_trace_anomaly`/`grouped_inverse`, `inverse x20/21/22 == 0`); wl L80-98 (`projector`, `Inverse[projector] == inverseTarget`, reconstructs lanes) | match |
| Isotropic common-lane collapse `a2=b2=a4=b4=0` | py L71-76; wl L100-103 | match |
| One-pole identity `Delta_pole = -(3D2^2+D0D4)/D0^2` and `D0D4+3D2^2=0` | py L84-103; wl L107-123 | match |
| One-parameter carrier series through `O(omega^4)` | py L111-120; wl L127-138 | match |
| Firewall: linear-`chi` block -> `chi^2` Schur correction, `partial_chi|_0=0` | py L137-162; wl L142-162 | match |

`paper_alignment: aligned`. Every paper deliverable has a faithful, non-tautological script-side check in BOTH engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 66-68 | `inverse x2k - x2k == 0` | claim 1 (inverse) | yes |
| A2 | sympy | 73-76 | `a2/b2/a4/b4 (iso) == 0` | claim 1 (collapse) | yes |
| A3 | sympy | 94-97 | `Delta_pole + (3D2^2+D0D4)/D0^2 == 0` | claim 2 | yes |
| A4 | sympy | 100-103 | `Delta_pole(D4=onepole) == 0` | claim 2 | yes |
| A5 | sympy | 120 | `Y_pole_series - Y_expected == 0` | claim 3 | yes |
| A6 | sympy | 153-156 | `Deff - (D2 I - chi^2 C^T C/D0s) == 0` | claim 4 | yes |
| A7 | sympy | 158-161 | `d/dchi Deff|_0 == 0`; `Deff|_0 - D2 I == 0` | claim 4 (firewall) | yes |
| A8 | math | 97-98 | `Inverse[projector]-inverseTarget==0`; reconstructs lanes | claim 1 | yes |
| A9 | math | 102-103 | `commonNu2/4[[2;;3]] == 0` | claim 1 (collapse) | yes |
| A10 | math | 117-119 | `nu2/nu4/Delta_pole` from `D[Y,{w,n}]` match closed forms | claim 2 | yes |
| A11 | math | 121-123 | `Solve[Delta_pole==0]` gives `d0d4+3d2^2=0`; `Delta_pole|_surface==0` | claim 2 | yes |
| A12 | math | 136-138 | carrier `w^2`/`w^4` coeffs match; series - target == 0 | claim 3 | yes |
| A13 | math | 157-158 | `D[offdiag,chi]==C`; `D[offdiag,{chi,2}]==0` | claim 4 (input is linear) | yes |
| A14 | math | 159-162 | `schurEff - schurTarget==0`; `det` identity; `d/dchi|_0==0`; `\|_0 - D2 I==0` | claim 4 (firewall) | yes |

No row is tautological: each assertion compares an independently-constructed quantity against a target it was not defined to equal. In §I both engines verify the inverse against an independently-stated `inverseTarget`/`grouped_inverse`. In §II `Delta_pole` is built from `nu2,nu4` (py hand-coded, wl from Taylor derivatives) then compared to the closed form `-(3D2^2+D0D4)/D0^2` — non-trivial. In §IV the Schur complement of a LINEARLY-coupled block (off-diagonals `= chi*C`) is computed via matrix inverse and shown to be quadratic in `chi`; the linear→quadratic promotion is genuine, not assumed.

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.txt` (mtime 2026-06-01 11:43:38)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py` (mtime 2026-06-03 15:59:11)

**What's wrong:**
The committed SymPy output `.txt` is older than the `.py` script (output 06-01, script 06-03). The content confirms it: the saved output's banner reads `STAGE 176 — EXACT ISOTROPIC GROUPED-P2 TARGET SURFACE...` and `STAGE 176 LEDGER` (lines 3 and 127), whereas the current `.py` source emits `STAGE 193` (py L49) and `STAGE 193 LEDGER` (py L164). So the committed output predates the banner renumber fix and does not reflect the current script. All numeric/symbolic check results in the stale output are nonetheless correct and match what the current script would print (the F1 first-pass firewall fix — explicit block operator with linear `chi` and the Schur complement at py L137-162 — is fully reflected in the saved output's §IV, lines 69-124). The Mathematica output (06-01 11:43:38) is also older than... no: the `.wl` is 06-01 11:41:10 and its `.txt` is 06-01 11:43:38, so the Mathematica output is FRESH relative to its script.

**Why this matters:**
Only the SymPy stage label in the captured transcript is stale; the math is current. The verifier's independent re-run will regenerate a `STAGE 193`-bannered transcript. This is the known SCRIPT/OUTPUT-band lag, not a math defect.

**Required change:**
Re-run `python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py` and overwrite the saved `.txt`. No source edit.

**Verification:**
New SymPy `.txt` banner reads `STAGE 193 ...` / `STAGE 193 LEDGER`; all `expect_zero` lines still print `= 0`.

## Independent-derivation check (Mathematica)

GENUINELY INDEPENDENT — and the firewall (the first-pass DIRTY section) is exercised independently in both engines, with the `.wl` carrying a STRONGER linear→quadratic exercise. Evidence:

- **§I inverse.** SymPy: hand-written `grouped_inverse(xbar,ax,bx)` returns `xbar+4ax`, `xbar-ax+bx`, `xbar-ax-bx` (py L42-46) and asserts `inverse - x == 0`. Mathematica: builds a 3x3 `projector` matrix (wl L81-85) and computes `Inverse[projector]` (wl L88), then checks it against an independently stated `inverseTarget` and also `Inverse.projector.lane - lane == 0` (wl L97-98). Different method (native matrix inversion vs hand-coded algebraic inverse).
- **§II coefficients.** SymPy hand-codes `nu2_common = -D2/D0`, `nu4_common = (D2^2 - D0 D4)/D0^2` (py L84-85). Mathematica derives them as Taylor coefficients: `nu2FromD = (D[responseShape[w],{w,2}]/2) /. w->0`, `nu4FromD = (D[responseShape[w],{w,4}]/24)/.w->0` (wl L110-111), and `onePoleD4 = d4 /. First[Solve[deltaPole==0,d4]]` (wl L121) — a `Solve`, where SymPy hand-substitutes `D4_onepole = -3 D2^2/D0` (py L99). Genuinely different routes to the same closed forms.
- **§IV firewall (the load-bearing re-check).** SymPy: `Dblock` via `BlockMatrix`, then `Deff = D2blk*I3 - (chi*Cvec.T)*(scalar.inv())*(chi*Cvec)` (py L137-142). Mathematica: `blockOperator` via `ArrayFlatten`, then `schurEff = l2Block - lowerMix.Inverse[scalarBlock].upperMix` (wl L148-152). Both build the off-diagonals as genuinely `chi`-LINEAR (`upperMix = chi*mixSeed`, `chi*Cvec`) and take the Schur complement, which produces the `chi^2` form — so linear→quadratic is exercised, not assumed, in BOTH. The `.wl` ADDS two checks SymPy lacks: (a) `D[blockOperator offdiag, chi] == C` and `D[..,{chi,2}] == 0` (wl L157-158), explicitly certifying the INPUT coupling is exactly linear before showing the output is quadratic; (b) the Schur determinant identity `Det[block] - Det[scalar]*Det[schurEff] == 0` (wl L160), an independent cross-validation of the Schur complement that SymPy does not perform. This is the opposite of a transliteration: the firewall fix is non-tautological and the `.wl` is the stronger of the two.

No CONSISTENT-threshold transliteration trigger: there is no section where the same expression is computed by the same operation. The choreography differs (algebraic inverse vs `Inverse[]`; hand-coded coefficients vs `D[]`/`Solve[]`; `BlockMatrix` vs `ArrayFlatten`+det-identity).

## Engine cross-check

Both engines PASS every check and agree on every emitted form (saved outputs confirm): `Delta_pole = -((3 d2^2 + d0 d4)/d0^2)` (wl txt L25) == `(-D0 D4 - 3 D2^2)/D0^2` (py txt L42-47); `Omega_Q^2 = -1/4 d0/d2` (wl txt L40) == `-D0/(4 D2)` (py txt L54-57); Schur complement `D2 I - chi^2 C^T C/D0scalar` identical in both (wl txt L53, py txt L77-91); all firewall residual matrices `{{0,0,0},...}`. No `engine_disagreement`.

## Verdict justification

`findings` — single low/informational `stale_output` on the SymPy `.txt` (label-band lag: committed output still bannered `STAGE 176`; the math content is current and correct, including the first-pass firewall fix). Attacks tried and failed: (1) firewall tautology — refuted: the off-diagonals are genuinely `chi`-linear and the Schur complement via matrix inverse is what promotes them to `chi^2`; the `.wl` additionally certifies linearity via `D[..,chi]` and validates via the det identity. (2) Transliteration — refuted: independent methods in all four sections (algebraic vs `Inverse[]`; hand-coded coeffs vs Taylor `D[]`/`Solve[]`; `BlockMatrix` vs `ArrayFlatten`+det-identity). (3) Paper misalignment — refuted: card Output, notes §2-§6, and appendix eqs all match the script's load-bearing assertions exactly. (4) `symbol_assumption_error` — `D0,D2,D4,D0scalar,D2blk` are declared `nonzero` (py L83,127-128; wl L75), which is required for the divisions/inverses and is physically the non-degenerate-block assumption; no over-assumption hides a branch. Paper card `\stagefield{Verification}` says "Mathematica audit: none yet" but a passing `.wl` exists — noted below as a paper-side card-text lag (paper-cleanup class), NOT a script finding.

## Self-test notes

Variable-independence: every `D[...,chi]`/`.diff(chi)` acts on `Deff`/`schurEff`, which genuinely depend on `chi` (off-diagonals `= chi*C`), so the `d/dchi|_0 == 0` checks are non-trivially satisfied (the `chi^2` term's derivative is `2 chi C^T C/D0`, which vanishes at `chi=0` — a real result, not a zero-derivative-of-a-constant trap). Parity/integration: no unbounded integrals in this stage. Trivial-case: substituting `C=0` collapses `Deff -> D2 I` correctly; substituting the one-pole `D4 = -3D2^2/D0` gives `Delta_pole = 0` as asserted. No directive trap triggered.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation. All emitted RESULT values are symbolic (no benchmark/figure-of-merit numbers in this stage). Sources cited as `py`/`wl` source line + saved-output line; doc locations are the stage card (`.tex`), notes (`.md`), and part appendix.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Isotropy surface `a2=b2=a4=b4=0` | py L73-76 / py-txt L24-27; wl L102-103 / wl-txt L15-18 | tex L20; md §2 L113, §6 L335; appdx eq `conservative-target-surface` L1226 | MATCH |
| One-pole defect `Delta_pole = -(3 D2^2 + D0 D4)/D0^2` | py L86,94-97 / py-txt L42-48; wl L112,119 / wl-txt L25,30 | md §3 L177; appdx (via eq L1228 + L1233) | MATCH |
| One-pole surface `D0 D4 + 3 D2^2 = 0` | py L99-103 / py-txt L49; wl L121-122 / wl-txt L32 | md §3 L183-185; appdx eq `one-pole-surface-D` L1233 | MATCH |
| `Omega_Q^2 = -D0/(4 D2)` | py L111 / py-txt L54-57; wl L127 / wl-txt L40 | md §4 L221-224; appdx (implicit in carrier eq) | MATCH |
| Carrier `Y_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` | py L112,120 / py-txt L58-64; wl L128,138 / wl-txt L41,46 | tex L15 (prose); md §4 L235-239, §6 L344-349; appdx eq `conservative-one-pole-module` L1243-1245 | MATCH |
| Firewall Schur `D_eff,l=2 = D2 I - chi^2 C^T C / D0scalar` | py L142,153-156 / py-txt L77-112; wl L152,159 / wl-txt L53,58 | md §5.2 L294-301; appdx eq `geometry-firewall-schur` L1250-1254 | MATCH |
| Firewall `partial_chi D_eff\|_{chi=0} = 0` (no O(chi)) | py L158-161 / py-txt L113-118; wl L161 / wl-txt L62 | tex L15; md §5.2 L308-310, §6 L357-361; appdx eq `geometry-firewall-result` L1259 | MATCH |

INTERNAL scaffolding (no finding): `xbar, ax, bx` projector expressions; `inverseTarget`/`grouped_inverse` reconstruction; `nu2_common`, `nu4_common` / `nu2FromD`, `nu4FromD`; `D4_onepole`/`onePoleD4`; `Y_expected`/`targetSeries`; carrier `w^2`/`w^4` coefficient cross-checks; `Det[block]-Det[scalar]*Det[schur]` identity; per-check residual matrices.

reconciliation: complete; 7 deliverable values checked, 0 misaligned.

Note (paper-cleanup class, NOT a script finding): the stage card `\stagefield{Verification}` (tex L11) states "Mathematica audit: none yet", but a passing `.wl` exists (`mathematica/moving_throat_pde_stage193_..._mathematica_audit.wl`). The notes §8 "Script-backed status" likewise lists only the SymPy audit. This is card-text lag from the dual-engine retrofit; resolution is paper-side (user/Codex paper-edit lane), not a script change.
