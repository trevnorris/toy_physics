---
unit_id: 135
batch: IV.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage135_outlet_consistent_mouth_closure.md"]
  paper_appendix: present
---

# Audit unit 135 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_135.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage135_outlet_consistent_mouth_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_135}` line at 1304; this part appendix is a manifest of `\input` rows, no per-stage prose row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.txt`

## What the paper claims

Stage 135 imposes the Stage-115 compensated-outlet weighting `4 : -1` onto the dimensionless mouth-layer gains, giving the **gain pair** `M_s = 4 Σ_m`, `M_q = -Σ_m` (one residual amplitude `Σ_m > 0`). Under that substitution the two-gain Family-1 mouth-layer equation **collapses to the one-parameter fixed-point law** `Π = Σ_m[4 - S_q(Π)]`. Requiring the imported canonical compensation point `Π_* ≈ 1.50882951349316` (carried from Stage 131) then selects the **exact outlet-consistent gain** `Σ_m^* = Π_*/(4 - S_q(Π_*)) ≈ 0.451485277739090`, with `M_s^* ≈ 1.80594111095636`, `M_q^* ≈ -0.451485277739090`, kernel value `S_q(Π_*) ≈ 0.6581 < 1`, and a mixed-lane correction `M_q^* S_q(Π_*) ≈ -0.297111597463199`. The card quotes verbatim: "Inheriting the compensated outlet ratio gives \(M_s=4\Sigma_m\), \(M_q=-\Sigma_m\), and \(\Sigma_m^*\approx0.451485\)." The card's `\stagefield{Checks}` adds a third deliverable: the fixed points are recorded as numerically located, not as closed-form constants — consistent with the `\StatusExactClosure` tag (a ledger entry, not an unconditional branch theorem).

## What the script claims to verify

The SymPy script verifies, in four steps: (1) the algebraic reduction `M_s + M_q S_q == Σ_m(4 - S_q)` under `M_s→4Σ_m, M_q→-Σ_m` (asserted residual 0); (2) the kernel-value inequality `0 < S_q(Π_*) < 1` (asserted, exercises the actual `S_q` evaluation); (3) the numeric anchors `Σ_m^*, M_s^*, M_q^*` against the notes values within tolerances 1e-12/1e-11/1e-12; (4) the mixed-lane correction `M_q^* S_q(Π_*)` against `-0.297111597463199` within 1e-11. It also prints (no assert) the closure residual `Π_* - Σ_*(4 - S_*)`, explicitly labelled "was the original tautological closure residual." The Mathematica script mirrors all of these and additionally asserts the closure residual ≈ 0 (line 78) and `3 Σ_m^* < Π_* < 4 Σ_m^*` (line 79).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Gain pair `M_s = 4Σ_m, M_q = -Σ_m` (card:16, notes:33-37) | py:36 / wl:52 substitution; reduction asserted py:40, wl:56 | match |
| One-parameter law `Π = Σ_m(4 - S_q(Π))` (card body 13-17 intent; notes:56-60,144-148) | py:38-40 reduction residual==0; wl:56 expectZero | match |
| `S_q(Π_*) ≈ 0.6581 < 1` (notes:68) | py:56-58 `0<S_star<1` asserted; wl:67 expectTrue | match |
| `Σ_m^* ≈ 0.451485` (card:16; notes:89,153) | py:69 anchor < 1e-12; wl:74 expectApprox 1e-15 | match |
| `M_s^* ≈ 1.80594111095636` (notes:97,121) | py:71 anchor < 1e-11; wl:75 expectApprox 1e-14 | match |
| `M_q^* ≈ -0.451485277739090` (notes:100) | py:73 anchor < 1e-12; wl:76 expectApprox 1e-15 | match |
| Mixed-lane correction `≈ -0.297111597463199` (notes:126) | py:81 anchor < 1e-11; wl:77 expectApprox 1e-14 | match |
| `Π_*` imported from Stage 131 (notes:81,109-114) | py:46 / wl:58 literal `1.50882951349316` | match (carry-forward) |
| Fixed points numerically located, not closed-form (card Checks 24) | scripts use `sp.Float`/numeric `Solve`, no symbolic closed form claimed | match |
| (none) `3 Σ_m^* < Π_* < 4 Σ_m^*` | wl:79 expectTrue | extra (consistent with notes:70-72 prose; not asserted in py) |

`paper_alignment: aligned`. Every paper-side deliverable has a faithful, non-tautological script-side anchor. The single `extra` row (wl:79) is a true and harmless restatement of the notes §2 prose bound; it is not a paper_misalignment.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 40 | `simplify(reduced_law - expected_law) == 0` | one-parameter reduction | partial — algebraic-identity-only (distributive law), independent of `S_q` value |
| A2 | sympy | 58 | `assert 0 < S_star < 1` | `S_q(Π_*)<1` deliverable | yes — exercises actual kernel value |
| A3 | sympy | 69 | `abs(Sigma_star - target) < 1e-12` | `Σ_m^*` | yes |
| A4 | sympy | 71 | `abs(M_s_star - target) < 1e-11` | `M_s^*` | yes |
| A5 | sympy | 73 | `abs(M_q_star - target) < 1e-12` | `M_q^*` | yes |
| A6 | sympy | 81 | `abs(mixed_correction - target) < 1e-11` | mixed-lane correction | yes |
| A7 | mathematica | 56 | `expectZero[reducedLaw - expectedLaw]` | one-parameter reduction | partial — same algebraic-identity-only as A1 |
| A8 | mathematica | 67 | `expectTrue[0 < sStar < 1]` | `S_q(Π_*)<1` | yes |
| A9 | mathematica | 74-77 | `expectApprox[..., 1e-15/1e-14]` | `Σ_m^*, M_s^*, M_q^*`, mixed | yes |
| A10 | mathematica | 78 | `expectApprox["closure residual", residual, 0, 1e-14]` | none — residual is 0 by construction | **no — tautological (see F1)** |
| A11 | mathematica | 79 | `expectTrue[3 σ* < Π_* < 4 σ*]` | notes §2 prose bound (extra) | yes |

A1/A7 are the algebraic-reduction checks. They confirm `4Σ - Σ·S_q ≡ Σ(4-S_q)` for ANY `S_q`, i.e. they exercise distributivity, not the kernel. They are weak but they faithfully verify exactly the algebraic step the paper claims (the two-gain → one-parameter collapse is itself a distributive identity), so they are not flagged as a defect. A10 is the only assertion that traces to no real claim and cannot fail — flagged below.

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.wl:78`

**What's wrong:**
The Mathematica script asserts the closure residual is zero:
```
residual = N[piStar - sigmaStar*(4 - sStar), 30];   (* line 64 *)
...
expectApprox["closure residual", residual, 0, 10^-14];   (* line 78 *)
```
But `sigmaStar` was DEFINED two lines earlier as the solution of exactly that equation:
```
sigmaStar = N[sigmaVar /. First[Solve[piStar == sigmaVar*(4 - sStar), sigmaVar, Reals]], 30];   (* line 60 *)
```
i.e. `sigmaStar = piStar/(4 - sStar)`, so `piStar - sigmaStar*(4 - sStar) ≡ 0` algebraically, independent of any physics or of the kernel value. The assertion cannot fail no matter what `S_q(Π_*)` is. The SymPy author already recognized this exact construction as circular and **removed the assert**, leaving only a print with the explicit comment (sympy:84): "Numerical sanity probe (was the original tautological closure residual)." The two engines are inconsistent: SymPy demoted this to a print, Mathematica still asserts it as a PASS line (`PASS: closure residual` in the saved output).

**Why this matters:**
A tautological PASS line masquerades as verification and inflates the apparent check count of the Mathematica engine. It contributes nothing to confirming the stage's claim; the substantive Mathematica anchors that DO exercise the kernel are lines 74-77 (`Σ_m^*, M_s^*, M_q^*`, mixed) and the range checks 67/79. Removing the dead assert keeps the two engines in the intended parity that SymPy already established.

**Required change:**
Bring the Mathematica side into parity with the SymPy side: at `.wl:78`, remove the `expectApprox["closure residual", residual, 0, 10^-14];` assertion. Retain the print of the residual already emitted at line 72 (`Print["Pi_star - Sigma_star*(4 - S_star) = ", residual];`) as a sanity probe, matching the SymPy script's sympy:85-86. Do not add any new assertion in its place — the kernel value is already independently exercised by the numeric anchors at lines 74-77.

**Verification:**
After the fix, `.wl:78` no longer contains an `expectApprox`/`expect*` assertion for "closure residual"; the saved Mathematica transcript no longer contains the `PASS: closure residual` line; all remaining `expect*` checks still PASS and the script still `Exit[0]`. The print at line 72 (`Pi_star - Sigma_star*(4 - S_star) = 0...`) is unchanged.

## Independent-derivation check (Mathematica)

The `.wl` is structurally a close port of the `.py`: identical kernel formula (`sKernel[p,k]` vs `S(Pi,kappa)`, character-for-character the same `p*(k*Tanh[k] + p*(Exp[-p]/Cosh[k]-1))/((1-Exp[-p])*(k^2-p^2))`), identical `genericLaw = ms + mq*sQ`, identical substitution `ms→4σ, mq→-σ`, identical hardcoded `Pi_star = 1.50882951349316`, identical `Solve` for the gain, and identical numeric anchors. I considered raising `mathematica_transliteration`. I am NOT raising it as a hard finding, for two reasons: (1) the kernel `S_q` is an **imported result from upstream stages (114-116/131)**, not something this stage derives — both engines must legitimately use the same imported formula, so the kernel parallelism is forced, not an echo; (2) the operative computation here is a one-line algebraic substitution plus a 1-D numeric root, which leaves essentially no room for an "alternate route." The Mathematica engine does provide genuine independent value: it re-evaluates the kernel and re-solves with native `Solve[..., Reals]`/`FullSimplify` and confirms the same numeric anchors to 1e-14/1e-15, and adds the `3σ*<Π_*<4σ*` cross-check. Net: the second engine is doing real cross-validation of the numeric anchors, which is the meaningful content of this stage. I note the structural parallelism for the record but it does not meet the bar for an actionable transliteration finding given the imported kernel.

## Engine cross-check

Both engines agree to full displayed precision:

| quantity | SymPy output | Mathematica output |
|---|---|---|
| reduction residual | `0` (py out:5) | `0`, PASS (wl out:6-7) |
| `S_q(Π_*)` | `0.658075937605428494...` (py out:16) | `0.6580759376054294867...` (wl out:8) |
| `Σ_m^*` | `0.451485277739089696...` (py out:18) | `0.4514852777390898089...` (wl out:11) |
| `M_s^*` | `1.80594111095635878...` (py out:19) | `1.8059411109563592359...` (wl out:12) |
| `M_q^*` | `-0.451485277739089696...` (py out:20) | `-0.4514852777390898089...` (wl out:13) |
| mixed corr | `-0.297111597463198745...` (py out:22) | `-0.2971115974631992675...` (wl out:14) |

Agreement is to ~14-15 significant figures (the engines differ only in the last ~2 digits of their 28-30 digit displays, consistent with independent finite-precision evaluation of the same kernel). I independently recomputed the kernel at 40 dps with mpmath and obtained `S_q(Π_*)=0.65807593760542948827...`, `Σ_m^*=0.45148527773908981865...`, `M_s^*=1.80594111095635927460...`, mixed `=-0.29711159746319927460...` — all consistent with both engines. No engine disagreement.

## Verdict justification

Attacks tried: (1) checked A1/A7 for tautology — they are algebraic-identity-only but they map to the paper's actual one-parameter-collapse claim, so defensible, not flagged; (2) checked whether `Pi_star` is an unanchored hardcoded result — it is a legitimate Stage-131 carry-forward stated in notes:81, so not `hardcoded_result`; (3) independently recomputed the kernel and every gain value at 40 dps and confirmed all numeric anchors and the `0<S_q<1`, `3σ*<Π_*<4σ*` bounds; (4) checked symbol domains (`positive=True, real=True` / `$Assumptions pi>0, sigmaM>0`) — consistent with `Σ_m>0` and `Π_*>0`, and the kernel denominator `(κ²-Π²)` with `κ=π/2≈1.571, Π_*≈1.509` is nonzero, no spurious branch; (5) confirmed engine agreement to display precision and verified the `Solve` is real-branch-restricted. The one real defect is the tautological `closure residual` assertion left in the Mathematica script (`.wl:78`) that the SymPy author already demoted to a print — a low-severity inconsistency, not a math error. Paper ↔ script alignment is exact across all eight deliverables. Verdict: `findings` (one low-severity, script-side, Codex-fixable), no stop-cold, no paper_misalignment.

## Self-test notes

I checked: (1) variable-independence — no `diff`/`D` derivatives appear in either script, so the zero-derivative trap does not apply; the only `Solve` is over the genuinely-present `sigmaVar`. (2) Trivial-case pre-check on F1 — substituting `sigmaStar = piStar/(4-sStar)` into `piStar - sigmaStar*(4-sStar)` collapses to literal 0, confirming A10 is tautological and the proposed removal (not replacement) is correct. (3) Paper round-trip — the F1 fix only removes a dead assert and adds nothing, so it introduces no new constant and cannot create a paper_misalignment; the retained print is unchanged. No symmetry/parity or integration traps apply (no integrals in this unit).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 8 deliverable values checked, 0 misaligned

Deliverable-level table (authoritative record = script source + committed saved outputs; both outputs are fresh, so values below are quoted from the committed `.txt`):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| gain pair `M_s = 4Σ_m`, `M_q = -Σ_m` (symbolic) | py:36, wl:52 | tex:16 `\(M_s=4\Sigma_m\), \(M_q=-\Sigma_m\)`; md:33-37 (boxed) | MATCH |
| one-parameter law `Π = Σ_m(4 - S_q(Π))` (symbolic) | py:42-44, wl:53-55; py out:6-15, wl out:5 | md:56-60 and md:144-148 (boxed Result) | MATCH |
| `S_q(Π_*) = 0.658075937605428...` | py out:16, wl out:8 | md:68 `\mathcal S_q(\Pi_*)\approx 0.6581<1` | MATCH (agrees to quoted 4 sig figs) |
| `Σ_m^* = 0.451485277739090` | py out:18, wl out:11 | tex:16 `\Sigma_m^*\approx0.451485`; md:89 and md:153 `\approx 0.451485277739090` | MATCH |
| `M_s^* = 1.80594111095636` | py out:19, wl out:12 | md:97 and md:121 `\approx 1.80594111095636` | MATCH |
| `M_q^* = -0.451485277739090` | py out:20, wl out:13 | md:100 `\approx -0.451485277739090` | MATCH |
| `M_q^* S_q(Π_*) = -0.297111597463199` | py out:22, wl out:14 | md:126 `\approx -0.297111597463199151` | MATCH (notes trailing digits differ at ~4e-16, far below the 1e-11 tolerance; computed value `-0.2971115974631987...` is the correct product of the quoted `M_q^*` and `S_q`) |
| `Π_* = 1.50882951349316` (imported anchor) | py:46, wl:58 | md:81 `\Pi=\Pi_*\approx 1.50882951349316` (carry-forward from Stage 131, md:109-114) | MATCH |

INTERNAL (scaffolding, no prose expected; raise no finding): `residual_sub`/`outlet-consistent reduction`==0 (algebraic check flag); `s_in_range`/`0<S_q<1` boolean; `residual`/`closure residual` near-zero probe (the tautological one, F1); `3 Σ_m^* < Π_* < 4 Σ_m^*` boolean (a true restatement of the notes §2 prose bound, not a numeric deliverable); all PASS/FAIL flags, tolerances (1e-11..1e-15), and `expectApprox` diffs.

All eight emitted deliverable values reconcile against the `.tex` card and/or the `.md` notes. No MISMATCH, no MISSING-DELIVERABLE. The only sub-1e-15 numeric jitter is in the notes' last-quoted digits of the mixed-lane correction (md:126), which is round-off in the prose, not a real value disagreement, and is well inside the script tolerance — not raised as a finding.
