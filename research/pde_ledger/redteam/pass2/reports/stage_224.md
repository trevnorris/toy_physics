---
unit_id: 224
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 224 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_224.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows 60, 599-679 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Actual-branch kill test: the branch survives only if the transported quantity $(\Delta_{\rm norm}+T_{\rm quad})(1+|\varepsilon\Xi_1|)/\hat m_0^2$ stays below the critical prefactor budget." The card + notes + appendix (Theorem `thm:app-part07-actual-branch-ceiling`) enumerate the deliverables: (1) the grouped inverse map $(\bar P_0,a_{P_0},b_{P_0})\!\mapsto\!(P_{20},P_{21},P_{22})$ with $P_{20}=\bar P_0+4a_{P_0}$, $P_{21}=\bar P_0-a_{P_0}+b_{P_0}$, $P_{22}=\bar P_0-a_{P_0}-b_{P_0}$; (2) the normalization compiler $\bar P_0=(\Delta_{\rm norm}+T_{\rm quad})/\hat m_0^2$, $T_{\rm quad}=54Gc_s^5/(5a^5c^5)$; (3) the isotropic ceiling $\Delta_{\rm norm}\le\hat m_0^2 P_{\rm crit}-T_{\rm quad}$; (4) the weak-axisymmetric lane law $P_A=\bar P_0(1+\epsilon\lambda_A\Xi_1)$ with $\lambda=(1,\tfrac12,-1)$; (5) the trace/anomaly compiler $a_{P_0}=\epsilon\bar P_0\Xi_1/4$, $b_{P_0}=3\epsilon\bar P_0\Xi_1/4$, hence $b_{P_0}=3a_{P_0}$; (6) the robust all-lane collapse $\max\{P_{20},P_{21},P_{22}\}=\bar P_0(1+|\epsilon\Xi_1|)=\bar P_0+4|a_{P_0}|$; (7) the calibrated ($\Delta_{\rm norm}=0$) lower bounds $\hat m_0^2\ge T_{\rm quad}/P_{\rm crit}$ and $\hat m_0^2\ge T_{\rm quad}(1+|\epsilon\Xi_1|)/P_{\rm crit}$; (8) four explicit headroom budgets at the Stage-223 compatibility point $\bar P_0\approx0.002069792318062885$ for the four carried ceilings. The numerical budgets and ceilings are `\StatusNumerical{}`; the symbolic compiler is exact.

## What the script claims to verify

The SymPy script asserts, symbolically: the grouped inverse projection recovers $(\bar P_0,a_{P_0},b_{P_0})$ from the three lanes (M1, lines 34-36); the isotropic-ceiling rearrangement $(\bar P_0-P_{\rm crit})\hat m_0^2=\Delta_{\rm norm}-(\hat m_0^2 P_{\rm crit}-T_{\rm quad})$ (line 52); the WA lane law and its trace/anomaly inversion $a_{P_0}=\epsilon\bar P_0\Xi_1/4$, $b_{P_0}=3\epsilon\bar P_0\Xi_1/4$, $b_{P_0}=3a_{P_0}$ (lines 76-78), plus a re-expansion round-trip (lines 84-86); the robust max-lane collapse in both sign cases via an explicit $z_{\rm abs}$ substitution for $|\epsilon\Xi_1|$ (lines 106-120); and four numerical headroom budgets each checked against its **ceiling** via a defining relation (lines 153-156). The Mathematica script verifies the same eight deliverables by an independent route (Solve / LinearSolve / Refine, see below) and additionally proves which lane dominates in each sign case.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) grouped inverse map | sympy M1 L34-36; wl M1 Solve L95-110 | match |
| (2) normalization compiler $\bar P_0=(\Delta+T)/\hat m_0^2$ | sympy L25 (used in all); wl L91 | match |
| (3) isotropic ceiling $\Delta\le\hat m_0^2 P_{\rm crit}-T$ | sympy L52; wl M2 L114-117 | match |
| (4) WA lane law, $\lambda=(1,1/2,-1)$ | sympy L65-71; wl L121-122 | match |
| (5) $a_{P_0},b_{P_0}$ compiler, $b=3a$ | sympy L76-78; wl M3 L129-131 | match |
| (6) robust collapse $=\bar P_0(1+|\epsilon\Xi_1|)=\bar P_0+4|a_{P_0}|$ | sympy L106-120; wl M5 L152-181 | match |
| (7) calibrated lower bounds on $\hat m_0^2$ | sympy L54,122 (printed); wl printed via forms | match (printed deliverable, also implied by ceiling checks) |
| (8) four headroom budgets at compat point | sympy L153-156; wl M6 L204-211 | match |

`paper_alignment: aligned`. Every paper deliverable has a load-bearing, dual-engine check. The single finding below is a card verification-status staleness (the card says no `.wl` exists, but one does), routed to the user.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 34-36 | `simplify(inv_* - *) == 0` | (1) inverse map | yes |
| A2 | sympy | 52 | `expand((Pbar-Pcrit)*mhat0^2) == expand(Delta - rhs)` | (3) isotropic ceiling | yes |
| A3 | sympy | 76-78 | `simplify(aP_wa - ...)==0`, `bP=3aP` | (4)(5) WA compiler | yes |
| A4 | sympy | 84-86 | `simplify(P_from_ab - P_wa)==0` | (1)+(4) round-trip | yes |
| A5 | sympy | 106-120 | sign-case diffs + `subs(Abs→zabs)` collapse | (6) robust collapse | yes |
| A6 | sympy | 155-156 | `assert_close((budget+1)*barP0, Pcrit)` etc | (8) headroom vs ceiling | yes |
| B1 | wl | 108-110 | `expectZero[(Solve.../packetLaneRules) - input]` | (1) inverse map (via Solve) | yes |
| B2 | wl | 114-117 | `expectZero[ceiling residual]` | (3) isotropic ceiling | yes |
| B3 | wl | 129-131 | `expectZero[aPwa-...]`, `bP=3aP` (LinearSolve) | (4)(5) WA compiler | yes |
| B4 | wl | 136-138 | `expectZero[reexpand - waLanes]` | (1)+(4) round-trip | yes |
| B5 | wl | 152-181 | `expectTrueUnder[lane dominates]` + `expectZeroUnder[max = Pbar(1+|.|)]` | (6) robust collapse (+ domination proof) | yes |
| B6 | wl | 204-211 | `expectZero[(budget+1)*barP0 - Pcrit]` etc | (8) headroom vs ceiling | yes |

No tautological or unanchored rows. Every row traces to a paper deliverable.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** notes_contradicts_script (card verification-status staleness)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_224.tex:11`

**What's wrong:**
The card states the verification status as:
> `\stagefield{Verification}{SymPy audit: \StageFile{scripts/...sympy_audit.py}.  Mathematica audit: none yet.}`

But a complete, passing Mathematica audit now exists and is committed:
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_mathematica_audit.wl` (6 sections M1-M6, all PASS in the committed `.txt`). The card understates the verification coverage that exists. This is paper-side text only; Codex cannot edit it.

**Why this matters:**
The card is the citable record of what is verified. Saying "Mathematica audit: none yet" when one exists and passes misrepresents the stage's dual-engine status to any downstream reader/citer.

**Required change:**
Routed to user (paper-side). See `## Resolve before fix_loop` in the directive. The natural fix is to update line 11 to reference the existing `.wl`. No script change.

**Verification:**
Card line 11 names the `.wl` path; no longer says "none yet".

## Independent-derivation check (Mathematica)

The `.wl` is **genuinely independent**, not a transliteration. Three corresponding sections, contrasted:

1. **Inverse map.** SymPy M1 *hardcodes* the inverse projection coefficients: `inv_bar = (P20 + 2*P21 + 2*P22)/5`, `inv_a = (2*P20 - P21 - P22)/10`, `inv_b = (P21 - P22)/2` (L30-32). Mathematica M1 instead *derives* the inverse by symbolically solving the forward linear system: `inverseRule = First[Solve[Thread[{p20Lane,p21Lane,p22Lane} == laneMatrix . {pbarUnknown,aUnknown,bUnknown}], {pbarUnknown,aUnknown,bUnknown}]]` (L95-100), then substitutes the lane packet. The `Solve` output (committed `.txt` L9: `pbarUnknown -> p20/5 + 2p21/5 + 2p22/5`, etc.) independently *reproduces* the SymPy hardcoded coefficients — that agreement is the cross-engine confirmation, not a shared derivation.

2. **WA trace/anomaly compiler.** SymPy M3 recovers $a_{P_0},b_{P_0}$ with the same hardcoded projection (`(2*P20_wa - P21_wa - P22_wa)/10`, `(P21_wa - P22_wa)/2`, L73-74). Mathematica M3 uses `LinearSolve[laneMatrix, waLanes]` (L123) — a different inversion mechanism on the forward lane matrix.

3. **Robust collapse.** SymPy M5 only asserts algebraic differences between lanes after substituting `zabs` for `Abs` (L106-120); it never proves *which* lane is the max. Mathematica M5 additionally *proves the domination* with `expectTrueUnder["...P20 dominates", waLanes[[1]]>=waLanes[[2]] && waLanes[[1]]>=waLanes[[3]], positiveCase]` (L152-156) and the symmetric negative case, using `Refine[Abs[eps xi1], ...]` to resolve the absolute value. This is strictly more than the SymPy route does.

The only shared items are the physical *inputs* (the lane matrix `{{1,4,0},{1,-1,1},{1,-1,-1}}` and the signature `{1,1/2,-1}`), which are paper premises, not algebra to be re-derived. **INDEPENDENCE CALL: independent.**

## Engine cross-check

Both engines emit identical symbolic forms and identical numerical budgets:
- Inverse map: sympy `.txt` L7-9 vs wl `.txt` L9-15 (Solve coefficients identical) — agree.
- WA compiler: sympy L20-21 `aP=Xi1*eps*(Delta+T)/(4*mhat0^2)`, `bP=3*...` vs wl L27 `(eps*(deltaNorm+tQuad)*xi1)/(4*mhat0^2)`, `3*...` — agree.
- Headroom budgets: sympy L31-34 vs wl L66/71/76/81 — agree to all printed digits (e.g. both `0.36793032849264599...`, `0.000190384841874108125`). Engines agree.

## Verdict justification

Both engines verify all eight paper deliverables; the assertions are non-tautological, well-anchored, and exercised over the right symbolic/numeric space. Attacks tried that failed: (a) checked the M1/M3 inversions for being value-vs-itself — they are not; the `.wl` derives the inverse independently via Solve/LinearSolve and the agreement with SymPy's hardcoded coefficients is genuine cross-engine confirmation; (b) checked the M5 `Abs` handling for hidden branch errors — the `.wl` resolves `Abs` via `Refine` under sign assumptions AND proves lane domination, closing the gap the SymPy `zabs` substitution leaves open; (c) checked the M6 headroom checks for the prior hardcoded_result defect — confirmed the fix holds (see Self-test). The single finding is a paper-card verification-status staleness ("Mathematica audit: none yet" while a passing `.wl` exists), which is paper-side and routes to the user; it does not affect the math. Verdict: `findings` (one low-severity paper_misalignment), `needs_user_resolution: true`.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 12 deliverable values checked, 0 misaligned.

The `.tex` card is intentionally terse and carries the symbolic constant $T_{\rm quad}=54Gc_s^5/(5a^5c^5)$ and the structural law $b_{P_0}=3a_{P_0}$; all numeric deliverables live correctly in the notes `.md` (a legitimate carrier per the augmentation guards). Scripts keep $T_{\rm quad}$ as a free symbol, so no numeric scale is emitted for it — correct.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| inverse map $P_{20}=\bar P_0+4a_{P_0}$, $P_{21}=\bar P_0-a_{P_0}+b_{P_0}$, $P_{22}=\bar P_0-a_{P_0}-b_{P_0}$ | py L26-28, sympy.txt L4-6; wl Solve, wl.txt L9 | notes L86-93; appendix L618-625 | MATCH |
| $\bar P_0=(\Delta_{\rm norm}+T_{\rm quad})/\hat m_0^2$ | py L25, sympy.txt L7,12 | notes L74-80; appendix L613-616 | MATCH |
| isotropic ceiling $\Delta_{\rm norm}\le\hat m_0^2 P_{\rm crit}-T_{\rm quad}$ | py L58, sympy.txt L13 | notes L167-171 | MATCH |
| calibrated bound $\hat m_0^2\ge T_{\rm quad}/P_{\rm crit}$ | py L59, sympy.txt L14 | notes L183-187 | MATCH |
| WA law $P_A=\bar P_0(1+\epsilon\lambda_A\Xi_1)$, $\lambda=(1,1/2,-1)$ | py L65-71, sympy.txt L17-19; wl.txt L26 | notes L196-217; appendix L632-643 | MATCH |
| $a_{P_0}=\epsilon\bar P_0\Xi_1/4$ | py L76, sympy.txt L20; wl.txt L27 | notes L222; appendix L650 | MATCH |
| $b_{P_0}=3\epsilon\bar P_0\Xi_1/4$; $b_{P_0}=3a_{P_0}$ | py L77-78, sympy.txt L21-22; wl.txt L27 | notes L224,228; appendix L652 | MATCH |
| robust collapse $=\bar P_0(1+|\epsilon\Xi_1|)=\bar P_0+4|a_{P_0}|$ | py L115-120,127, sympy.txt L27 | notes L247-251,255-259 | MATCH |
| calibrated weak bound $\hat m_0^2\ge T_{\rm quad}(1+|\epsilon\Xi_1|)/P_{\rm crit}$ | py L122, sympy.txt L28 | notes L286-290; appendix L671-675 | MATCH |
| carried ceilings 0.0028313316855593175 / 0.0035965105896846573 / 0.00817339430971383 / 0.0116633929790174 | py L137-141; wl L188-191 | notes L103,105,110,112 | MATCH |
| compat point $\bar P_0=0.002069792318062885$ | py L134; wl L186 | notes L13,301 | MATCH |
| four headroom budgets ($|\epsilon\Xi_1|$ and $|a_{P_0}|$, all 4 ceilings) | sympy.txt L31-34; wl.txt L66,71,76,81 | notes L316,319,324,327,335,338,343,346 | MATCH |

Internal scaffolding (no prose expected, no finding): pass/fail flags, `residual = 0` lines, `laneMatrix`, `inverseRule` intermediate, `zabs`/`z` substitution symbols, `Refine` outputs, `assert_close` tol `1e-12`, the round-trip re-expansion intermediates.

## Self-test notes

(1) Variable-independence / derivative traps: N/A — this stage has no `sp.diff`/`D[...]`; all checks are linear-algebra inversions and rearrangements. (2) The M6 hardcoded_result→defining-relation fix HOLDS: budgets are computed (`epsXiBudget = pcrit/barP0 - 1`, `aBudget = (pcrit - barP0)/4`) and each is checked back against the **ceiling** `pcrit` (`(budget+1)*barP0 - pcrit == 0`, `barP0 + 4*aBudget - pcrit == 0`), not against itself; a wrong rearrangement (e.g. `/2` instead of `/4`) would fail, so the check is load-bearing. (3) Trivial-case sanity: with `eps*xi1 > 0`, `waLanes[[1]] - waLanes[[2]] = pbarIso*eps*xi1/2 > 0` so M5's domination `expectTrueUnder` is genuinely substantive (output confirms `= True`, not an unevaluated conditional). Independence confirmed via distinct inversion routes (hardcoded projection in py vs Solve/LinearSolve in wl) plus the wl-only domination proofs.
