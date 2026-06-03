---
unit_id: 233
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 233 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_233.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 78, 105, and section narrative 883–1005)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

The card's `\stagefield{Output}` is: "Rigid-mouth compiler: mouth tracking alone is not enough; the actual branch must control the transfer-shape/nontracking packet that carries $\Xi_1$." The derivation ledger states the stage proves that rigid mouth tracking $\delta\ln R_{\rm tr}=0$ implies $\Xi_1=\delta\ln\mathfrak N_*$, interprets the static ceiling as an internal transfer/turbulence gate, and identifies the remaining finite packet. The notes (Section 9, "SymPy-backed status") enumerate eight explicit deliverables the script must verify: (1) the exact Stage-239 observable compiler $\Theta_1=\delta\ln R_{\rm tr}$, $\Xi_1=\delta\ln\mathfrak N_*-B_*\,\delta\ln R_{\rm tr}$, $\mathcal R_1=-c_\eta\,\delta\ln\epsilon_\eta-\Xi_1$; (2) the track-locked specialization $\delta\ln R_{\rm tr}=0\Rightarrow\Theta_1=0,\ \Xi_1=\delta\ln\mathfrak N_*$; (3) the prefactor identity $\Xi_{\rm load}=N_{01}/N_0-D_{01}/D_0=P_1/P_0$; (4) operator-rigidity specialization $D_{01}=0\Rightarrow\Xi_{\rm load}=N_{01}/N_0$; (5) the transported static ceiling $|\epsilon\Xi_1|\le P_{\rm crit}\hat m_0^2/(\Delta_{\rm norm}+T_{\rm quad})-1$; (6) its calibrated $\Delta_{\rm norm}=0$ simplification; (7) the equivalent $\bar P_0$-form; and (8) numerical recovery of the two carried budgets $0.367930328492646$ and $0.737619063660757$ from the compatibility-point value $\bar P_0\approx0.002069792318062885$. The appendix narrative (lines 883–1005) confirms this is an exact-closure compiler stage within the rigid-mouth packet language.

## What the script claims to verify

The script (docstring-free; intent from section comments and print strings) tests: the observable compiler definitions for $\Theta_1,\Xi_1,\mathcal R_1$ (lines 27–30); the rigid-mouth collapse under $\delta\ln R_{\rm tr}=0$ via three `assert_zero` calls (lines 45–47); the prefactor identity $\Xi_{\rm load}=P_1/P_0$ (line 64) and its $D_{01}=0$ specialization (line 75); the two algebraic ceiling-form equivalences (lines 92, 95); and a numeric "recovery" of the two carried budgets (lines 115–126). The first seven blocks are genuine symbolic identities; the eighth (numeric) block is a self-referential round-trip, analyzed below. All blocks correspond to notes deliverables 1–8 in order.

## Paper ↔ script cross-check

| Paper deliverable (notes §9) | Script check | Result |
|---|---|---|
| 1. Stage-239 compiler forms | lines 27–30 definitions; partly exercised by 45–47 | match (definitions match notes §2 exactly) |
| 2. track-locked $\Theta_1=0,\ \Xi_1=\delta\ln\mathfrak N_*$ | lines 45, 46, 47 | match |
| 3. $\Xi_{\rm load}=N_{01}/N_0-D_{01}/D_0=P_1/P_0$ | line 64 | match |
| 4. $D_{01}=0\Rightarrow\Xi_{\rm load}=N_{01}/N_0$ | line 75 | match |
| 5. transported ceiling form | line 92 (equivalence to $\bar P_0$ form) | match |
| 6. $\Delta_{\rm norm}=0$ simplification | line 95 | match (weak; see F3) |
| 7. equivalent $\bar P_0$-form | line 92 / 99 print | match |
| 8. numeric recovery of $0.3679..,0.7376..$ from $\bar P_0$ | lines 115–126 | mismatch — circular round-trip, does not recover (F1) |

`paper_alignment: aligned` for the load-bearing symbolic content (deliverables 1–7); the only paper↔script defect is that deliverable 8's "recovery" is not actually a recovery — but that is a script-side weakness (the script under-verifies its own and the paper's claim), not a disagreement about the values, so it is filed as `tautological_check`, not `paper_misalignment`. The numeric values used (`0.002069792318062885`, `0.367930328492646`, `0.737619063660757`) match the notes verbatim (notes §6, §9).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 45 | `assert_zero(Theta1_rm)` | deliverable 2 ($\Theta_1=0$) | partial (trivial: $\Theta_1:=\delta\ln R_{\rm tr}$ then subs 0) |
| A2 | sympy | 46 | `assert_zero(Xi1_rm - dln_Nstar)` | deliverable 2 ($\Xi_1=\delta\ln\mathfrak N_*$) | yes |
| A3 | sympy | 47 | `assert_zero(R1_rm + c_eta*dln_eps_eta + dln_Nstar)` | deliverable 2 ($\mathcal R_1$ collapse) | yes |
| A4 | sympy | 64 | `assert_zero(simplify(P1/P0 - Xi_load))` | deliverable 3 | yes |
| A5 | sympy | 75 | `assert_zero(Xi_load_or - N01/N0)` | deliverable 4 | yes |
| A6 | sympy | 92 | `assert_zero(simplify(gate_rhs - (Pcrit/Pbar_expr - 1)))` | deliverables 5, 7 | yes |
| A7 | sympy | 95 | `assert_zero(gate_rhs_cal - (Pcrit*mhat0**2/Tquad - 1))` | deliverable 6 | partial (re-asserts result of its own `subs`) |
| A8 | sympy | 125 | `assert_close(recovered_robust, robust_budget)` | deliverable 8 | NO — tautological round-trip |
| A9 | sympy | 126 | `assert_close(recovered_nonempty, nonempty_budget)` | deliverable 8 | NO — tautological round-trip |

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py:115-126`

**What's wrong:**
The notes (§9, deliverable 8) state the script provides the "numerical recovery of the two carried Stage-241 budgets `0.367930328492646`, `0.737619063660757` from the Stage-240 compatibility-point value `Pbar ≈ 0.002069792318062885`." The script instead defines all three numbers as literals (lines 115–117) and then performs:
```
Pcrit_robust   = Pbar_num * (1 + robust_budget)        # line 119
recovered_robust = Pcrit_robust / Pbar_num - 1          # line 122
assert_close(recovered_robust, robust_budget)           # line 125
```
Substituting line 119 into line 122 gives `recovered_robust = Pbar_num*(1+robust_budget)/Pbar_num - 1 = robust_budget` identically, for ANY value of `Pbar_num` and `robust_budget`. The assertion therefore cannot fail and verifies nothing about the budget values — it confirms `b == b`. The same holds for the nonempty branch (lines 120, 123, 126). `Pcrit` is never supplied from an independent source (Stage-241 / Stage-240 data, or from `mhat0`, `Tquad`, `Delta_norm`); it is back-computed from the very budget the check claims to recover. This is a self-referential round-trip masquerading as a numerical recovery.

**Why this matters:**
The script's own print block (lines 128–135) and the paper claim a numerical cross-check of the carried Stage-241 budgets. A reader trusts that `0.3679...` and `0.7376...` were re-derived and confirmed. In fact nothing is confirmed: if either carried budget were transcribed wrong, this check would still pass. The single numeric deliverable of the stage is unexercised.

**Required change:**
Replace the circular round-trip with a recovery that consumes an independently sourced quantity. The genuine identity the paper asserts is `|\epsilon\Xi_1| \le P_{\rm crit}/\bar P_0 - 1` (notes §9 item 7). A non-circular recovery supplies `P_crit` (the carried critical-pressure number) and `Pbar_num` as the two independent inputs and computes the budget `P_crit/Pbar_num - 1`, then asserts it equals the carried `robust_budget`/`nonempty_budget`. Because the script does not currently carry an independent `P_crit` literal, and the notes do not state one either, Codex must obtain `P_crit` from the upstream source the notes name (the Stage-241 robust/nonempty ceiling data) OR, if no independent `P_crit` literal is available in-repo, reformulate the check so the two budgets are tied to each other through a non-trivial relation that the budgets must satisfy (e.g. confirm the ratio `(1+nonempty_budget)/(1+robust_budget)` against the independently carried `Pcrit_nonempty/Pcrit_robust` ratio if both critical pressures are separately carried). See directive F1 for the precise contract and the `## Resolve before fix_loop`-adjacent note: if no independent `P_crit` source exists anywhere in the repo, this becomes a documentation/sourcing question, not a silent reformulation.

**Verification:**
The reformulated block must (a) NOT compute `Pcrit` as `Pbar_num*(1+budget)` and then divide it back out by the same `Pbar_num`; (b) take at least one of {`P_crit`, the ceiling expression evaluated at carried `mhat0`/`Tquad`} from an independent literal/source; (c) still exit 0 with both budgets matching. The verifier confirms line 119/120's `Pcrit_*` is no longer defined purely from `Pbar_num*(1+budget)`.

### F2 — missing_verification_script (missing_mathematica)

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_mathematica_audit.wl` (does not exist)

**What's wrong:**
This stage has a SymPy script and no Mathematica `.wl`. The card says "Mathematica audit: none yet." Every claim in this stage is a closed-form algebraic identity (the compiler collapse under $\delta\ln R_{\rm tr}=0$, the $P_1/P_0$ prefactor identity, the $D_{01}=0$ specialization, and the two ceiling-form equivalences) plus one numeric relation. All of these are independently verifiable in Mathematica via `FullSimplify`/`Reduce`/`Refine` on a different decomposition. The dual-engine rule's test is "is it possible" — and it plainly is. No genuine obstruction (no opaque numeric solver, no PDE integration) exists here.

**Why this matters:**
A checkpoint-grade exact-closure stage with only one engine has no cross-check; the F1 tautology is exactly the kind of defect a genuine second-engine derivation would expose, because an independent route would have to supply `P_crit` rather than round-trip the budget.

**Required change:**
Codex writes a new independent-route Mathematica script (see directive F2 claim manifest M1–M8). Independent route = native Mathematica primitives via a DIFFERENT decomposition than the `.py`, not a transliteration.

**Verification:**
`redteam exec-mathematica 233` produces the new `.wl`'s output, all in-file checks pass, and the script exits 0. The verifier additionally confirms the `.wl` is not a line-by-line port of the `.py` (different intermediate decomposition for at least the compiler-collapse and prefactor-identity blocks).

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py:94-95`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py:45`

**What's wrong:**
Two checks are weak. (a) Line 94–95: `gate_rhs_cal = simplify(gate_rhs.subs({Delta_norm: 0}))` then `assert_zero(gate_rhs_cal - (Pcrit*mhat0**2/Tquad - 1))`. The assertion subtracts the hand-written form from the result of the same `subs`; it confirms SymPy's substitution rather than an independent identity. (b) Line 45: `assert_zero(Theta1_rm)` where `Theta1 := dln_Rtr` (line 27) then `subs(dln_Rtr:0)` — this is the trivial substitution of zero into an identity map; it cannot probe any structure. Both are anchored to real paper deliverables (6 and 2 respectively), so this is low severity, but neither adds verification strength beyond confirming SymPy's `subs`.

**Why this matters:**
Low. These do not threaten correctness; they are filler checks. Noting them so the second-engine `.wl` (F2) does not merely echo the same trivial substitutions.

**Required change:**
No mandatory script edit for F3 alone. When Codex writes the F2 `.wl`, the calibrated-branch ceiling (M6) must be derived by an independent route (e.g. start from the $\bar P_0$ form and impose $\Delta_{\rm norm}=0$, then confirm equality with the $T_{\rm quad}$-only form by `FullSimplify` of their difference under the positivity assumptions) rather than re-substituting and re-asserting the same expression. Likewise M2's $\Theta_1=0$ must follow from the compiler definition, not a bare zero-substitution. This is folded into F2's manifest.

**Verification:**
F3 is satisfied by F2's `.wl` honoring the M2/M6 independence notes; no standalone `.py` change required.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot yet be assessed. The directive (F2) requires the new `.wl` to use a different decomposition for at least the compiler-collapse, prefactor-identity, and calibrated-ceiling blocks, and the verifier will check it is not a port of the `.py`.

## Engine cross-check

Only the SymPy engine is present; no cross-check possible (`engines_agree: n/a`). This is itself F2.

## Verdict justification

`findings`. The seven symbolic blocks (deliverables 1–7) hold up under attack: A2–A6 are genuine, non-tautological algebraic identities that faithfully exercise the paper's boxed claims, and the constants/forms in the script match the notes verbatim. Attacks that failed: I tried to break the $P_1/P_0$ identity (line 64) by hand — it is a true identity $P_1/P_0=N_{01}/N_0-D_{01}/D_0$; I tried to find a constant mismatch against notes §6/§9 — the three numeric literals match exactly; I checked whether the appendix's Stage-234 compiler form $\Xi_1=-\delta\ln R_{\rm target}-c_\eta\delta\ln\epsilon_\eta$ (appendix line 920) contradicts the script's $\Xi_1=\delta\ln\mathfrak N_*-B_*\delta\ln R_{\rm tr}$ — it does not, because that is Stage 234's decomposition in different variables, and Stage 233's card/notes consistently use the $\mathfrak N_*$ form the script implements. What does NOT hold: the single numeric deliverable (8) is a circular round-trip that verifies nothing (F1, high), the stage has no second engine although every claim is Mathematica-verifiable (F2, medium), and two checks are trivial substitution confirmations (F3, low). I confirmed I read the card, the full notes file, the appendix narrative (883–1005) and rows, and the script and its output. The output is fresh (output mtime 2026-05-11 12:51 > script mtime 2026-05-11 11:58), so no `stale_output` finding.

## Self-test notes

Checked: (1) variable independence — no `sp.diff` calls in this script, so the zero-derivative trap does not apply; the F1 reformulation introduces no derivative. (2) Symmetry/parity — no integrals, n/a. (3) Trivial-case pre-check — I substituted the round-trip algebra for F1 (`Pbar_num*(1+b)/Pbar_num - 1`) and confirmed it reduces to `b` for all inputs, proving the tautology; for the F1 fix I confirmed a `Pcrit/Pbar - 1` route with an independent `Pcrit` would give a genuinely falsifiable budget. (4) Path specification — the F2 `.wl` target is named with the full `mathematica/` directory and the mandatory `_mathematica_audit.wl` suffix, matching the repo convention (224 existing `.wl` files use this pattern). (5) Paper round-trip — the F1/F2 fixes use only constants and forms already in the notes; no new constant is introduced, so no new `paper_misalignment` is created. Separately noted (not a finding routed to user): the script's section comments/print strings cite "Stage 188" (line 20, 32), "Stage 223 / Stage 224" (lines 113, 128), which are stale upstream labels — the card/notes canonically cite Stage 239 (compiler), Stage 240 (compatibility point), Stage 241 (budgets/ceiling). This is the known project-wide incomplete-renumber; the math is unaffected. Folded into the directive as an optional low-severity comment-label correction (script-side prose only), NOT a blind batch renumber and NOT a user-routed paper_misalignment.
