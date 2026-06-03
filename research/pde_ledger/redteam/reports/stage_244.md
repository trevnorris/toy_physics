---
unit_id: 244
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T00:00:00-06:00
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 244 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_244.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (row at line 86)
- sympy: `/var/projects/toy_projects/.../scripts/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.py` → `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 244 compiles the first two open-system observables that switch on once the `J^w ≠ 0` lane (declared in Stage 243) is opened: the projected leakage source `S_leak` and the scalar-photon work channel `J^w E_w`. Using a normalized Gaussian projector `W_λ = e^{-w²/λ²}/(λ√π)` and an odd profile `φ_λ = 2w e^{-w²/λ²}/(√π λ³)` with `E_w = -E₀ φ_λ`, `j^w = μ_w ρ₀ E_w`, `J^w = q j^w`, the card states the exact closed forms: `S_leak = √2 μ_w ρ₀ /(2√π λ³) · E₀` (eq:sleak-e0), `W_w^bulk = √2 μ_w q ρ₀/(2√π λ³) · E₀²` (eq:bulk-work), `W_w^bulk = q E₀ S_leak` (eq:work-leak-relation), `W_w^sess = 2 μ_w q ρ₀/λ² · E₀² = 2√(2π) λ W_w^bulk` (eq:session-work), and the quadratic law `W_w^sess = 4π q λ⁴/(μ_w ρ₀) · S_leak²` (eq:quadratic-law). It then pulls the amplitude back through the Stage-242 selected-support demand `Π_tr = (4/3)C_mix = 32Λ(1-ε)/(3π²) = (16/π²)Λϱ` with `ϱ=(2/3)(1-ε)` (eq:pitr) via `E₀ = η_leak Π_tr` (eq:e0-pullback), giving e.g. `W_w^sess = 2048 η²_leak μ_w q ρ₀/(9π⁴λ²) · Λ²(1-ε)²` (eq:work-epsilon). The two structural deliverables are the **support-versus-orbit split** `∂_{R_tr} S_leak = ∂_{R_target} S_leak = ∂_{ε_η} S_leak = 0` (eq:support-orbit-split, "the formulas depend only on selected-support variables") and the **orientation parity** `S_leak(-η)=-S_leak(η)`, `W_w(-η)=W_w(η)` (eq:parity). The appendix row (line 86) summarizes exactly this list.

## What the script claims to verify

The docstring enumerates seven items: (1) the projected leakage source by Gaussian integration of `W'·j^w`, (2) the one-mode bulk work scalar by integrating `J^w E_w`, (3) the Session-I scalar as the thickness-rescaled bulk scalar, (4) the Π_tr pullback, (5) the support-versus-orbit split, (6) orientation parity, (7) recovery at `η_leak=0`. The assertions: Section 1 verifies the boundary term vanishes and the `S_leak` integral equals the closed form; Section 2 verifies `W_bulk` integral = closed form, `W_bulk = q E₀ S_leak`, `W_sess = 2√(2π)λ W_bulk` = closed form, and `W_sess = 4π q λ⁴ S_leak²/(μ_w ρ₀)`; Section 3 verifies `Π_tr`, its ϱ form, and the four pulled-back observables against independently written closed forms; Section 4 checks the support/orbit split via `diff(expr, R_tr/R_target/eps_eta)==0` and recovery at `η_leak=0`; Section 5 checks parity under `η_leak → -η_leak`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `S_leak = √2 μ_w ρ₀/(2√π λ³) E₀` (eq:sleak-e0) | L44-56 integrate `W'·j^w`, compare `S_expected` | match |
| boundary term vanishes (notes §2) | L43,55 `boundary==0` | match |
| `W_bulk = √2 μ_w q ρ₀/(2√π λ³) E₀²` (eq:bulk-work) | L62-63,76 | match |
| `W_bulk = q E₀ S_leak` (eq:work-leak-relation) | L66,77 | match |
| `W_sess = 2 μ_w q ρ₀/λ² E₀² = 2√(2π)λ W_bulk` (eq:session-work) | L64-65,78 | match |
| `W_sess = 4π q λ⁴ S_leak²/(μ_w ρ₀)` (eq:quadratic-law) | L67,79 | match |
| `Π_tr = 32Λ(1-ε)/(3π²) = 16Λϱ/π²` (eq:pitr) | L87-91,118-119 | match |
| `E₀ = η_leak Π_tr` (eq:e0-pullback) | L93,113 | match |
| compiled `S_leak`, `W_bulk`, `W_sess` (incl. eq:work-epsilon `2048/9π⁴`) | L96-108,120-125 | match |
| support/orbit split `∂_{R_tr}=∂_{R_target}=∂_{ε_η}=0` (eq:support-orbit-split) | L133-142 `diff(expr,var)==0` | **mismatch (vacuous)** — see F1 |
| orientation parity (eq:parity) | L160-171 | match |
| recovery at `η_leak=0` (notes §5.2) | L144-154 | match |

`paper_alignment: aligned` — every constant the script carries matches the paper card verbatim, and every assertion traces to a paper deliverable. The single weak point is the support/orbit-split *method*, which is vacuous (F1), not misaligned in claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 55 | `boundary == 0` | boundary vanishes (notes §2) | yes |
| A2 | sympy | 56 | `simplify(S_leak - S_expected) == 0` | eq:sleak-e0 | yes |
| A3 | sympy | 76 | `simplify(W_bulk - W_bulk_expected) == 0` | eq:bulk-work | yes |
| A4 | sympy | 77 | `simplify(W_bulk - W_bulk_from_S) == 0` | eq:work-leak-relation | yes |
| A5 | sympy | 78 | `simplify(W_sess - W_sess_expected) == 0` | eq:session-work | yes |
| A6 | sympy | 79 | `simplify(W_sess - W_sess_from_S) == 0` | eq:quadratic-law | yes |
| A7 | sympy | 118 | `simplify(Pi_tr - Pi_expected) == 0` | eq:pitr | yes |
| A8 | sympy | 119 | `simplify(Pi_tr_varrho - Pi_expected_varrho) == 0` | eq:pitr (ϱ form) | yes |
| A9-A14 | sympy | 120-125 | pulled-back `S/W_bulk/W_sess` vs closed forms | compiled formulas incl. eq:work-epsilon | yes |
| A15-A23 | sympy | 140-142 | `diff(expr, R_tr/R_target/eps_eta) == 0` (×3 exprs) | eq:support-orbit-split | **no (vacuous)** |
| A24-A26 | sympy | 152-154 | `*_rec == 0` (η_leak=0) | recovery slice (notes §5.2) | yes |
| A27 | sympy | 169 | `simplify(S_flip + S_pull_varrho) == 0` | eq:parity (odd) | yes |
| A28-A29 | sympy | 170-171 | `simplify(W_*_flip - W_*_pull_varrho) == 0` | eq:parity (even) | yes |

## Findings

### F1 — insufficient_verification (VARIABLE-INDEPENDENCE SELF-TEST TRAP)

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.py:131-142`

**What's wrong:**
This is the textbook support-blindness self-test trap. Section 4 declares fresh orbit symbols on line 131:
```
R_tr, R_target, eps_eta = sp.symbols("R_tr R_target eps_eta", real=True)
```
then for each compiled observable asserts:
```
d_Rtr     = sp.simplify(sp.diff(expr, R_tr));     assert d_Rtr == 0
d_Rtarget = sp.simplify(sp.diff(expr, R_target)); assert d_Rtarget == 0
d_epseta  = sp.simplify(sp.diff(expr, eps_eta));  assert d_epseta == 0
```
But `expr` is one of `S_pull_varrho`, `W_bulk_pull_varrho`, `W_sess_pull_varrho`. Tracing their construction: `S_pull = S_expected.subs(E0, E0_pull)` → free symbols `{lam, mu_w, rho0, eta_leak, Lam, eps}`; after the `1-eps → (3/2)varrho` subs, `{lam, mu_w, rho0, eta_leak, Lam, varrho}` (plus `q` for the works). At no point do `R_tr`, `R_target`, or `eps_eta` enter these expressions — they are brand-new symbols created two lines earlier and never substituted into anything. Therefore `sp.diff(expr, R_tr)` is **identically zero by construction**, exactly as the saved output confirms (lines 45-53: all nine print `0`). The assertion `== 0` can never fail no matter what the physics is — even if `S_leak` secretly *did* depend on an orbit variable, that variable would not be the symbol `R_tr` declared here, so the derivative would still be zero. The paper's central structural deliverable (eq:app-part08-stage244-support-orbit-split, "The formulas depend only on selected-support variables") is the headline result of Section 5 of the notes, yet the script's "proof" of it is vacuous.

**Why this matters:**
The support-versus-orbit split is one of only two structural theorems the stage proves (the other being parity); the appendix row (line 86) names it explicitly. A vacuous check provides zero assurance that a future edit (e.g. an amplitude pullback that accidentally pulls in an orbit coordinate) would be caught. The correct statement of support-blindness is a *free-symbol containment* claim, not a derivative.

**Required change:**
Replace the derivative-based checks with a non-vacuous structural test that the compiled observables' free symbols exclude the orbit variables. Concretely, after computing the three pulled-back observables, assert for each:
```
orbit_syms = {R_tr, R_target, eps_eta}
for expr in (S_pull_varrho, W_bulk_pull_varrho, W_sess_pull_varrho):
    assert orbit_syms.isdisjoint(expr.free_symbols)
    # and confirm they DO depend on the support variables, so the test is not vacuous:
    assert {Lam, varrho, eta_leak}.issubset(expr.free_symbols)  # (drop eta_leak for the works? no — both depend on eta_leak)
```
The `isdisjoint` test genuinely fails if any orbit symbol leaks into a compiled formula; the `issubset` positive control confirms the observables are non-trivial functions of the support variables (so the disjointness is meaningful, not "expr is a constant"). Keep the derivative prints if desired for human readability, but the load-bearing assertion must be the symbol-set test. See directive F1 for exact code.

**Verification:**
After the fix, `expr.free_symbols` for each of the three observables must contain `{Lam, varrho, eta_leak}` (and `mu_w, rho0, lam`, plus `q` for works) and must be disjoint from `{R_tr, R_target, eps_eta}`. The script must still exit 0.

### F2 — missing_verification_script (missing_mathematica)

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_mathematica_audit.wl` (does not exist)

**What's wrong:**
Stage 244 is SymPy-only. The dual-engine rule (test = "is it possible," not "is it necessary") clearly applies: every claim is a Gaussian integral over `(-∞,∞)`, an algebraic substitution, or a parity/free-symbol statement — all of which Mathematica verifies natively with `Integrate`, `Limit`, `FullSimplify`, `D`, and `FreeQ`. There is no transcendental obstruction, no special-function gap, no reason a second engine cannot independently re-derive `S_leak`, `W_bulk`, `W_sess`, the pullback, the parity, and the support-blindness. This stage is a checkpoint: `False` per the prompt, but `is_status_only_candidate: False` and not a carry-forward — both engines are required.

**Why this matters:**
A second-engine confirmation of the two Gaussian integrals (`∫ W'·j^w` and `∫ J^w E_w`) is exactly the kind of cross-check the policy exists to provide, since a single CAS could in principle mis-handle the odd-Gaussian moment. No Mathematica corroboration currently exists.

**Required change:**
Create the independent-route `.wl` named in directive F2. It must NOT transliterate the `.py` (different decomposition: e.g. compute the leakage moment via `Integrate[Wprime jw,{w,-Infinity,Infinity}]` directly using native Gaussian-moment evaluation, and independently confirm `W_bulk` via the second-moment integral, rather than reusing the same intermediate-variable choreography). Must verify all six claim-manifest items below.

**Verification:**
`redteam exec-mathematica 244` produces a `.wl` output that exits 0 with all `expectZero`/`expectTrue` checks passing, and the leakage/work closed forms match the SymPy script's (engine cross-check then becomes available).

### F3 — paper_misalignment? (NOT raised — notes typo only; informational)

**Severity:** low (informational; NOT a script defect, NOT routed to user gate as a blocking item)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.md:366`

**What's wrong:**
The notes §4.2 line 366 writes the ϱ-form bulk work scalar as `196√2 η²_leak μ_w q ρ₀ Λ²ϱ²/(π^{9/2}λ³)`. The correct coefficient, obtained from the `(1-ε)` form `512√2/9` (notes line 355, which the script confirms at L105) via `(1-ε)²=(9/4)ϱ²`, is `512√2/9 · 9/4 = 128√2`. The script line 106 correctly uses `128√2` and the saved output (line 63) prints `128√2`. So **the script and the paper card are both correct**; the `196√2` in the notes is a stray typo. The published paper card (`stage_244.tex`) does not carry this ϱ-form coefficient at all (it only states the `(1-ε)` forms), so there is **no paper_misalignment in the published card** and no script change is warranted.

**Why this matters:**
Only flagged so the orchestrator can, separately from the script red-team, correct the prose notes typo (`196 → 128`). It does not block, does not require user resolution for a value choice (the value is unambiguously `128√2`), and the red-team does not edit notes. No directive entry.

**Required change:** none for scripts. (Prose-side note for the orchestrator: notes line 366 coefficient should read `128√2`, matching the script and the derivable value.)

**Verification:** n/a (informational).

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration is moot for now. The directive F2 prescribes an independent decomposition (native Gaussian-moment integration and a distinct ordering of the pullback substitutions) explicitly to avoid producing a transliteration.

## Engine cross-check

Not applicable — only one engine present. After F2 lands, the leakage/work closed forms (`S_leak`, `W_bulk`, `W_sess`, `Π_tr`, and the `2048/9π⁴` work-epsilon coefficient) become cross-checkable between engines.

## Verdict justification

The arithmetic core of the stage holds up under attack: the two Gaussian integrals are evaluated genuinely (not resubstituted), every closed form is written independently and compared via `simplify(... - ...) == 0`, every load-bearing constant in the script matches the paper card verbatim (including the `2048/(9π⁴)` work-epsilon coefficient and the `Π_tr = 32Λ(1-ε)/(3π²)` pullback), the parity checks substitute `η→-η` into genuinely η-dependent expressions, and the recovery slice is a real `η_leak=0` substitution into nonzero expressions. The verdict is `findings` for two reasons: (F1) the support-versus-orbit split — one of the stage's two structural theorems — is "proved" by differentiating w.r.t. orbit symbols that the compiled expressions never contain, making all nine of those derivatives identically and vacuously zero (the support-blindness self-test trap); and (F2) the stage is single-engine despite being fully Mathematica-verifiable, so a second-engine `.wl` is required. No `paper_misalignment` against the published card (the only mismatch is a prose typo in the notes, `196√2` vs the correct `128√2`, which the script gets right). No stop-cold condition: F1 is mechanically fixable and F2 is additive; neither propagates a changed value downstream.

## Self-test notes

Checked the variable-independence trap explicitly and it IS the live defect (F1): `R_tr/R_target/eps_eta` are absent from all three differentiated expressions, so `diff(...)==0` is vacuous — my prescribed `isdisjoint(free_symbols)` plus an `issubset` positive control is non-vacuous (it fails if an orbit symbol leaks in and the positive control rules out the "expr is constant" degenerate pass). Parity: `S_pull_varrho` is degree-1 in `eta_leak` (odd) and the works are degree-2 (even); `η→-η` substitution checks are genuine, not vacuous. Symmetry/parity of the integrands: `W_λ` is even, `φ_λ` is odd, so `W'·j^w` ∝ (odd)·(odd)=even ⇒ nonzero leakage integral (matches `√2/(2√π)` moment), and `J^w E_w` ∝ (odd)² = even ⇒ nonzero work integral; both assertion targets are correctly nonzero. Round-trip: the `2√(2π)λ·W_bulk` and `q E₀ S_leak` relations are checked against independently-written closed forms, not against themselves. Paper round-trip: confirmed the `128√2`/`512√2`/`2048` constants are mutually consistent under `(1-ε)=(3/2)ϱ` and match the card.
