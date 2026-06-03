---
unit_id: 234
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 234 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_234.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 80, 905-938)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

The card's `\stagefield{Output}` reads: "Direct static gate: on the rigid-mouth branch, the static same-charge test is decided by the realized pair $(R_{\rm target},\epsilon_\eta)$, with tracking data factored out." The `\stagefield{Derivation ledger}` says the derivation "shows two cancellations, reduces the static gate to $R_{\rm target}$ and $\epsilon_\eta$, and states the two-observable kill test." The notes (section 9) enumerate six SymPy-backed deliverables: (1) the exact finite quotient chart $(q_{\rm tr},q_{\rm nt},q_\eta)$ in $(R_{\rm tr},R_{\rm target},\epsilon_\eta)$ and its inverse; (2) the first-order linearization of the chart; (3) the exact cancellation of $R_{\rm tr}$ out of $\Xi_1$, i.e. $\partial_{\delta\ln R_{\rm tr}}\Xi_1=0$; (4) the rigid-mouth reduction $q_{\rm nt}=\Xi_1$ at $\delta\ln R_{\rm tr}=0$; (5) the strip form of the robust ($0.367930328492646$) and nonempty ($0.737619063660757$) static gates on the $(R_1,E_1)$ plane; and (6) the canonical direct-branch families (pure-target, pure-dressing, balanced minimal-norm). The triangular compiler is $\Theta_1=\delta\ln R_{\rm tr}$, $\Xi_1=-\delta\ln R_{\rm target}-c_\eta\,\delta\ln\epsilon_\eta$, $\mathcal R_1=\delta\ln R_{\rm target}$ with $c_\eta:=\epsilon_{\eta,*}/(1-\epsilon_{\eta,*})$. The two ceiling constants are not derived here; they are imported branch constants carried forward from upstream stages (224/225/226/230/231/233).

## What the script claims to verify

The script (docstring/print at line 17) audits the "direct branch-observable static gate and the two-observable kill test." Its assertions test: the inverse-map roundtrip of the finite log chart (lines 55-57); the first-order drifts $q_{\rm tr}^{(1)}=-C_*r_1$, $q_{\rm nt}^{(1)}=B_*r_1-c_\eta E_1-R_1$, $q_\eta^{(1)}=E_1$ (lines 88-90); the compiler rows $\Theta_1=r_1$, $\Xi_1=-R_1-c_\eta E_1$, $\mathcal R_1=R_1$ and the inverse $\delta\ln\epsilon_\eta$ row (lines 105-110); the $R_{\rm tr}$ cancellation $\partial_{r_1}\Xi_1=0$ and the rigid-mouth reduction $\Xi_1|_{r_1=0}=q_{\rm nt}^{(1)}|_{r_1=0}$ (lines 122, 134-135); the strip-endpoint algebra (lines 160-163); and the three canonical families incl. the Lagrange-multiplier balanced minimal-norm solve (lines 177-200). All assertions use `assert_zero` = `simplify(expr) != 0 -> raise`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (1) finite chart + inverse roundtrip | lines 49-57 (`qtr/qnt/qeta_roundtrip - q_* == 0`) | match |
| (2) first-order linearization | lines 82-90 | match |
| compiler triangular map ($\Theta_1,\Xi_1,\mathcal R_1$) + $\epsilon_\eta$ inverse | lines 101-110 | match |
| (3) $R_{\rm tr}$ cancels: $\partial_{\delta\ln R_{\rm tr}}\Xi_1=0$ | line 122 | match (load-bearing, genuine cancellation) |
| (4) rigid-mouth $q_{\rm nt}=\Xi_1$ | lines 124-135 | match |
| (5) robust + nonempty strip form | lines 146-163 | partial (form verified, ceiling values never exercised — tautological endpoints) |
| (6) canonical families incl. balanced minimal-norm | lines 173-200 | match |
| second engine (Mathematica) | none | missing |

`paper_alignment: aligned` — every paper deliverable has a corresponding script check; the only weakness is the strip-endpoint tautology, which still tests the correct *form*.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 55 | `simplify(qtr_roundtrip - q_tr) == 0` | (1) chart inverse | yes |
| A2 | sympy | 56 | `simplify(qnt_roundtrip - q_nt) == 0` | (1) chart inverse | yes |
| A3 | sympy | 57 | `simplify(qeta_roundtrip - q_eta) == 0` | (1) chart inverse | yes |
| A4 | sympy | 88 | `qtr1 + C_star*r1 == 0` | (2) linearization | yes |
| A5 | sympy | 89 | `qnt1 - (B_star*r1 - c_eta*E1 - R1) == 0` | (2) linearization | yes |
| A6 | sympy | 90 | `qeta1 - E1 == 0` | (2) linearization | yes |
| A7 | sympy | 105 | `Theta1 - r1 == 0` | compiler row 1 | yes |
| A8 | sympy | 106 | `Xi1 + R1 + c_eta*E1 == 0` | compiler row 2 ($\Xi_1$) | yes |
| A9 | sympy | 107 | `Rcal1 - R1 == 0` | compiler row 3 | yes |
| A10 | sympy | 110 | `E1_inv - E1 == 0` | inverse drift map | yes |
| A11 | sympy | 122 | `diff(Xi1, r1) == 0` | (3) $R_{\rm tr}$ cancellation | yes (genuine; r1 truly present pre-simplify) |
| A12 | sympy | 126 | `qnt_rigid_finite != desired -> raise` | (4) rigid finite form | partial (structural-equality guard, not `simplify`) |
| A13 | sympy | 134 | `Theta1_rigid == 0` | (4) rigid mouth | yes |
| A14 | sympy | 135 | `Xi1_rigid - qnt1_rigid == 0` | (4) rigid $q_{\rm nt}=\Xi_1$ | yes |
| A15 | sympy | 160-163 | `Xi1.subs(R1, strip_edge) ∓ const == 0` (x4) | (5) strip form | no — tautological (const cancels by construction) |
| A16 | sympy | 177 | `Xi_target - xi == 0` | (6) pure-target family | yes |
| A17 | sympy | 181 | `Xi_dressing - xi == 0` | (6) pure-dressing family | yes |
| A18 | sympy | 190-200 | Lagrange solve + `R1_bal/E1_bal/Xi` checks | (6) balanced minimal-norm | yes (genuine `solve`) |

## Findings

### F1 — missing_verification_script (missing_mathematica)

**Severity:** medium
**Files:**
- `(missing)` — no `.wl` exists; `scripts/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.py` is the only engine.

**What's wrong:**
The unit is `is_checkpoint: false` but also `is_status_only_candidate: false`, so the dual-engine rule applies: both engines are required wherever Mathematica CAN independently verify. Every deliverable here is squarely in Mathematica's wheelhouse via native primitives and a different decomposition than the SymPy script uses:
- The finite chart + inverse roundtrip (claim 1) is a `Simplify`/`PowerExpand` log identity.
- The first-order linearization (claim 2) — the SymPy script differentiates an `exp(eps*x)` substitution at `eps=0`; Mathematica can instead use `Series[..., {eps,0,1}]` / `Normal` / `Coefficient`, an independent route.
- The $R_{\rm tr}$ cancellation (claim 3) is `D[Xi1, r1]` or, independently, forming $\Xi_1=q_{\rm nt}^{(1)}+(B_*/C_*)q_{\rm tr}^{(1)}$ and `FullSimplify`-ing to confirm no `r1`.
- The strip form (claim 5) is direct substitution.
- The balanced minimal-norm family (claim 6) — the SymPy script uses an explicit Lagrange-multiplier `solve`; Mathematica can use `Minimize[{R1^2+E1^2, -R1-c_eta E1==xi}, {R1,E1}]` (a genuinely different route) and confirm the same closed form.

So Mathematica can independently verify the stage, which makes the missing `.wl` a finding (not a legitimate single-engine carve-out).

**Why this matters:**
Single-engine coverage means a SymPy-specific simplification quirk (e.g. a `simplify`/`expand_log(..., force=True)` branch choice on the log chart) could mask a real defect with no cross-check. The cancellation in claim 3 is the marquee physics result; it deserves a second, independently-routed confirmation.

**Required change:**
Codex writes a NEW independent-route Mathematica script (see directive F1 for path, claim manifest, and acceptance criteria). Independent route = native Mathematica primitives via a DIFFERENT decomposition than the `.py` (Series/Normal for the linearization, Minimize for the balanced family), NOT a transliteration of the SymPy `diff(..., eps).subs(eps,0)` choreography.

**Verification:**
`redteam exec-mathematica 234` produces a `.wl` output transcript that exits 0 with all in-file checks passing; the new `.wl` independently reproduces $\Xi_1=-R_1-c_\eta E_1$, $\partial_{r_1}\Xi_1=0$, the rigid-mouth $q_{\rm nt}=\Xi_1$, and the balanced minimal-norm $(R_1,E_1)=(-\xi/(1+c_\eta^2),-c_\eta\xi/(1+c_\eta^2))$.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage234..._sympy_audit.py:146-163`

**What's wrong:**
The Section 5 strip checks are algebraically guaranteed by construction and never exercise the ceiling values. The script defines (line 150) `robust_lower = -c_eta*E1 - robust/Abs(vareps)`, then (line 155) substitutes it into `Xi1_twoobs = -R1 - c_eta*E1` and (line 160) asserts `Xi1_twoobs.subs(R1, robust_lower) - robust/Abs(vareps) == 0`. Expanding: `-(-c_eta*E1 - robust/Abs(vareps)) - c_eta*E1 = robust/Abs(vareps)`, so the residual is identically zero for ANY value of `robust` — the literal `0.367930328492646` (line 146) cancels and is never tested. The same holds for the upper edge (line 161) and both `nonempty` edges (lines 162-163). Replacing the constants with arbitrary placeholders would leave all four assertions passing. The checks confirm only the trivial rearrangement "if $R_1=-c_\eta E_1\mp B/|\varepsilon|$ then $\Xi_1=\pm B/|\varepsilon|$," which restates the strip-equivalence the notes already present as "Equivalently."

Note: the constants themselves are legitimate *carried imports* (they appear identically in notes for stages 224/225/226/230/231/233), so this is NOT a `hardcoded_result` against the paper and NOT a `paper_misalignment` — the value provenance is upstream. The defect is purely that the check is tautological in the value.

**Why this matters:**
The strip-form deliverable (notes section 9 item 5, card "two-observable kill test") is reported as SymPy-backed, but the only thing actually backed is the `Xi1 = -R1 - c_eta*E1` rearrangement (already covered by A8). If the ceiling values were ever mis-transcribed (e.g. a copy error swapping robust/nonempty, or a digit drop), these assertions would still pass green, giving false confidence.

**Required change:**
Replace the tautological endpoint identities with a check that actually pins the two ceiling constants and their ordering. Concretely, add assertions that (a) the strip half-width for the robust gate equals `robust/Abs(vareps)` AND that `robust < nonempty` (the nonempty corridor must be the looser bound), and (b) confirm the robust and nonempty constants are distinct nonzero literals carried as stated (e.g. assert their numeric difference equals the documented `nonempty - robust` and that neither is the other). Do not invent a derivation for the constants — they are carried imports; just make the check fail if the values or their ordering are corrupted. (Codex: keep the existing `Xi1 = -R1 - c_eta*E1` rearrangement print; only the four `assert_zero` endpoint identities at lines 160-163 are the tautology to strengthen.)

**Verification:**
After the fix, swapping `robust` and `nonempty` literals (or zeroing one) must make the script exit nonzero; the new assertion appears near line 160 and the saved output reflects the ordering check.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration is not yet a risk — but the directive (F1) must pre-empt it: the new script must use `Series`/`Normal` for the linearization and `Minimize` for the balanced family, NOT a port of the SymPy `sp.diff(expr, eps).subs(eps, 0)` and explicit Lagrange `sp.solve` choreography.

## Engine cross-check

N/A — only the SymPy engine is present. Cannot compare residuals. This is exactly the gap F1 closes.

## Verdict justification

`findings`. The SymPy script is paper-aligned: all six notes-enumerated deliverables map to non-tautological, well-anchored assertions, with the marquee $R_{\rm tr}$-cancellation (A11) being a genuine algebraic cancellation (I verified `r1` truly appears in both `qnt1` and `qtr1` before `simplify`, so `diff(Xi1, r1) == 0` is substantive, not a zero-derivative trap), and the balanced minimal-norm family verified via a real `solve` (A18). Attacks that failed: (i) I checked the chart roundtrip is not circular — the inverse map is built from independent closed forms and substituted back, not by re-using the forward expression; (ii) I checked `Rcal1 = -Xi1 - c_eta*qeta1` (A9) is a legitimate compiler row recovering $R_1=\delta\ln R_{\rm target}$, not a self-referential definition; (iii) I checked symbol domains — `eps_eta`, `eps_eta_star`, `Rtr/Rtarget` etc. are `positive=True`, which is justified by the physical setup (radii and a dressing fraction $\in(0,1)$) and makes the `expand_log(..., force=True)` log-splits valid. Two real defects remain: the missing second engine (F1, dual-engine rule applies — Mathematica clearly CAN verify) and the tautological strip-endpoint checks (F2, the carried ceiling values are never actually exercised). Neither is stop-cold; both are script-side and Codex-fixable.

## Self-test notes

Checked variable-independence trap: all three `eps`-derivatives (lines 82-84) act on expressions that genuinely contain `eps` via the `exp(eps*...)` substitutions, and `diff(Xi1, r1)` (line 122) operates on a pre-simplify expression where `r1` is genuinely present in both summands — so no identically-zero derivative masquerading as a pass. Checked trivial-case for A11: substituting concrete $B_*, C_*$ confirms the `r1` term cancels only because the $B_*/C_*$ weight is exact, so the `assert_zero` is load-bearing. Checked the F2 fix round-trip against the paper: the strengthened ordering/value check introduces no new constant and no new paper claim (the two ceilings are already stated verbatim in the notes), so no `paper_misalignment` is created. No symmetry/parity integrals appear in this unit.
