---
unit_id: 195
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch.md]
  paper_appendix: present
---

# Audit unit 195 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_195.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 121, 135; narrative 1306–1318; anchor MTDC-T9.5 line 1467)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The card's `\stagefield{Output}` states verbatim: "Factorizes the observable odd closure as \(\widehat m_0^{2}\chi_QN_Q=1\) and collapses \(\Delta_{\rm norm}\)." The appendix narrative (line 1306–1310) repeats this as the load-bearing headline `\widehat m_0^{2}\chi_QN_Q=1`, plus the natural source map `\widehat m_0\to1,\;N_Q=\chi_Q^{-1}` (line 1311–1317). The notes enumerate five deliverables: (1) the exact isotropic retarded odd ratio `Γ̄₅/Γ̄₅^target = χ_Q N_Q`; (2) the boxed observable odd-closure factorization `m̂₀²χ_Q N_Q = 1`, derived from the observable odd condition `m̂₀²Γ̄₅ = Γ̄₅^target` together with deliverable (1); (3) the exact collapse `Δ_norm = P0_target(1/χ_Q − 1) = −P0_target Δ_Q/(1+Δ_Q)`; (4) the natural source-map reduction `m̂₀→1 ⟹ N_Q = 1/χ_Q`; (5) the explicit DtN-deformation-algebra forms for χ_Q, N_Q, Δ_norm^pt, their linearization in (β,Σ₀,Σ₅), and the canonical-branch specialization (β=1,Σ₀=Σ₅=0 ⟹ χ_Q=N_Q=1, Δ_norm=0). The card is the "original-stage audit card" and explicitly "restores original stage order"; the notes body uses the alias numbering 246 (with upstream 242/244/245). This stage imports the DtN fingerprint and source map; it introduces no new constitutive law.

## What the script claims to verify

The SymPy script (banner "STAGE 178", a stale internal alias) walks the same five blocks. It checks: (I) the two boxed forms of Γ̄₅ agree, the N_Q definition, and the odd ratio `Γ̄₅/Γ̄₅^target = χ_Q N_Q`; (II) the Packet-A residual `Δ_norm = P0_target(m̂₀²N_Q − 1)`, an "odd closure factorization" line, the collapse `Δ_norm = P0_target(1/χ_Q − 1)`, and the Δ_Q form; (III) the point-particle reductions; (IV) the DtN-deformation forms for χ_Q, N_Q, Δ_norm^pt and the first-order linearization in (eps_beta, dSigma0, dSigma5); (V) the canonical-branch specialization. All checks route through `expect_zero`, which `simplify(expand(...))`s and raises on nonzero. Output is fresh (output mtime 12:48 > script mtime 11:58) and exits 0.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) `Γ̄₅/Γ̄₅^target = χ_Q N_Q` | lines 64–67 | match (genuine substitution P0→N_Q·P0_target) |
| (2) `m̂₀²χ_Q N_Q = 1` (boxed Output) | line 89 `odd_closure − (m̂₀²χ_Q N_Q − 1)` | **mismatch/tautological** — `odd_closure` is *defined* as that expression (line 76) then asserted equal to itself; the derivation `m̂₀²Γ̄₅ = Γ̄₅^target` + (1) ⟹ factorization is never composed |
| (3) `Δ_norm = P0_target(1/χ_Q − 1)`, Δ_Q form | lines 90, 94–97, 112 | match |
| (4) `m̂₀→1 ⟹ N_Q = 1/χ_Q` | lines 111, 113–116 | match (weakly definitional, but substitution m̂₀→1 is the content) |
| (5) χ_Q, N_Q, Δ_norm^pt deformation forms + linearization + canonical | lines 139–183 | match (Δ_norm^pt and linearization are the strongest, fully non-tautological checks) |

Dominant pattern is `match`, but the single headline deliverable (2) — the card's `\stagefield{Output}` — is verified by a self-comparison, so `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `simplify(Gamma5 − Gamma5_alt) == 0` | del 1 (Γ̄₅ two forms agree) | yes |
| A2 | sympy | 63 | `simplify(N_Q_def − P0/P0_target) == 0` | definition only | no (tautological: N_Q_def := P0/P0_target) |
| A3 | sympy | 64–67 | `simplify(Gamma5/Gamma5_target).subs(P0,N_Q·P0_target) − χ_Q·N_Q == 0` | del 1 | yes |
| A4 | sympy | 88 | `Delta_norm_NQ − P0_target(m̂₀²N_Q−1) == 0` | del 3 form | no (algebraically guaranteed by the subs that built Delta_norm_NQ) |
| A5 | sympy | 89 | `odd_closure − (m̂₀²χ_Q N_Q − 1) == 0` | del 2 (boxed Output) | **no — tautological self-echo** |
| A6 | sympy | 90 | `Delta_norm_from_odd − P0_target(1/χ_Q−1) == 0` | del 3 | yes |
| A7 | sympy | 94–97 | `Delta_norm in Delta_Q == 0` | del 3 | yes |
| A8 | sympy | 111 | `NQ_pt − 1/χ_Q == 0` | del 4 | partial (m̂₀→1 subst is content; NQ_from_odd hand-defined) |
| A9 | sympy | 112 | `Delta_norm_pt − P0_target(1/χ_Q−1) == 0` | del 3/4 | yes |
| A10 | sympy | 113–116 | `(NQ_pt−1) in Delta_Q == 0` | del 4 | yes |
| A11 | sympy | 139–142 | `NQ_from_def − (3S−Σ₀)/(3(Sβ⁵+9Σ₅)) == 0` | del 5 | partial (reciprocal of defined χ_from_def; can still catch a typo in the RHS) |
| A12 | sympy | 143 | `Delta_norm_from_def − expected == 0` | del 5 | yes (strong) |
| A13 | sympy | 160–163 | `linearized N_Q − 1 == 0` | del 5 | yes (strong; genuine series) |
| A14 | sympy | 164–167 | `linearized Delta_norm == 0` | del 5 | yes |
| A15 | sympy | 174–183 | canonical-branch specializations == 0 | del 5 | yes |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py:76` and `:89`

**What's wrong:**
The card's headline `\stagefield{Output}` is the boxed factorization `m̂₀²χ_Q N_Q = 1`. The notes (section 2, lines 130–141) derive it: from the observable point-particle odd normalization condition `m̂₀²Γ̄₅ = Γ̄₅^target`, substituting the exact ratio `Γ̄₅/Γ̄₅^target = χ_Q N_Q` (deliverable 1) gives `m̂₀²χ_Q N_Q = 1`. The script never composes these. It defines

```
odd_closure = sp.simplify(mhat0**2 * chi_Q * N_Q - 1)   # line 76
```

and then "verifies" it with

```
expect_zero("odd closure factorization", odd_closure - (mhat0**2 * chi_Q * N_Q - 1))   # line 89
```

i.e. `simplify(m̂₀²χ_Q N_Q − 1) − (m̂₀²χ_Q N_Q − 1) == 0`. This is `X − X == 0`; it cannot fail regardless of the physics, and it does not use the observable odd condition `m̂₀²Γ̄₅ = Γ̄₅^target` or the Γ̄₅ definition at all. The one deliverable the whole stage is named for is therefore asserted, not derived. (The ratio deliverable 1 is genuinely checked at lines 64–67, but the *composition into the factorization* — which is the Output — is missing.)

**Why this matters:**
The stage's load-bearing claim, the box that the appendix cites at `\resultanchor{MTDC-T9.5}` and that Stages 196/197 build on (`m̂₀²χ_Q N_Q=1` feeding `Δ_branch=0 ⟺ χ_Q=1`), is unverified. If the upstream Γ̄₅ or Γ̄₅^target normalization had a hidden factor error, the genuine derivation would catch it; the current self-echo would still pass. This is exactly a "passes for the wrong reason" hole on the most important assertion in the unit.

**Required change:**
Replace the self-echo with the actual derivation of the factorization from the observable odd condition. Insert a check that the observable odd condition residual factors through the verified ratio:

```python
odd_condition_residual = sp.simplify(
    (mhat0**2 * Gamma5 - Gamma5_target).subs(P0, N_Q * P0_target)
)
expect_zero(
    "observable odd condition factorizes as Gamma5_target*(mhat0^2 chi_Q N_Q - 1)",
    odd_condition_residual - Gamma5_target * (mhat0**2 * chi_Q * N_Q - 1),
)
```

Keep the existing `odd_closure` symbol if downstream lines reference it, but the load-bearing assertion at line 89 must be the non-tautological one above (it may replace line 89 outright). This derives `m̂₀²Γ̄₅ − Γ̄₅^target = Γ̄₅^target(m̂₀²χ_Q N_Q − 1)` from the carried Γ̄₅ = χ_Q a⁵P0/(27c_s⁵) and Γ̄₅^target = 2G/(5c⁵), so the box now follows from the observable condition rather than being declared.

**Verification:**
A new `expect_zero("observable odd condition factorizes ...", ...) = 0` line appears in the output between the section-II banner and the `Delta_norm ...` lines; the script still exits 0. Confirm by hand: `m̂₀²·χ_Q a⁵P0/(27c_s⁵)` with `P0=N_Q·54Gc_s⁵/(5a⁵c⁵)` equals `m̂₀²χ_Q N_Q·2G/(5c⁵) = m̂₀²χ_Q N_Q·Γ̄₅^target`, so the residual is `Γ̄₅^target(m̂₀²χ_Q N_Q − 1)` and the check is zero only because of the correct normalization, not by construction.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py:63`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py:88`

**What's wrong:**
Two secondary checks are definitional echoes that cannot fail:
- Line 63: `expect_zero("N_Q definition", N_Q_def - P0 / P0_target)` where `N_Q_def = sp.simplify(P0 / P0_target)` (line 52). This is `simplify(P0/P0_target) − P0/P0_target == 0`.
- Line 88: `expect_zero("Delta_norm - P0_target*(mhat0^2 N_Q - 1)", Delta_norm_NQ - P0_target * (mhat0**2 * N_Q - 1))` where `Delta_norm_NQ` (line 75) is exactly `mhat0**2*P0 - P0_target` with `P0 → N_Q*P0_target`, i.e. `P0_target*(mhat0**2*N_Q - 1)` by construction. So the subtraction is identically zero.

Neither exercises any physics; each restates a definition already made a few lines above.

**Why this matters:**
Low impact on its own — these are bookkeeping restatements, not the headline. But they pad the PASS transcript with checks that can never fail, which obscures that the real headline (F1) is the one that wasn't tested. They should either be removed or upgraded to actually anchor a paper form.

**Required change:**
Either delete lines 63 and 88, or make them substantive. The minimal, safe edit is deletion of the two `expect_zero` calls at lines 63 and 88 (the surrounding `print`/`pprint` lines may stay for display). Do not introduce any new physics here; if F1's fix already establishes the Δ_norm form, line 88 is redundant. If you prefer to keep line 88 as a display anchor, leave the `print` but drop the `expect_zero`.

**Verification:**
The two `... = 0` echo lines ("N_Q definition = 0" at output line 39, "Delta_norm - P0_target*(mhat0^2 N_Q - 1) = 0" at output line 65) no longer appear (or are clearly marked display-only); the script still exits 0.

## Independent-derivation check (Mathematica)

No `.wl` exists for this unit, so the `mathematica_transliteration` category does not apply.

## Engine cross-check

Only the SymPy engine is present; no cross-check possible. `engines_agree: n/a`.

On `missing_mathematica` (line-114 judgment): NOT flagged. Every deliverable of this stage is a closed-form symbolic algebra identity — a rational-function ratio (del 1), a factorization of a normalization condition (del 2), a reciprocal/substitution collapse (del 3,4), and a first-order Taylor series plus specializations of a rational function (del 5). SymPy fully and genuinely settles all of these; there is no transcendental evaluation, no nontrivial integral, no branch-cut ambiguity, and no numerical tolerance where a second CAS would add independent confidence. The only verification gap (F1) is a SymPy-side tautology to be fixed in SymPy, not a "needs a second engine" gap. Single-engine is acceptable here, consistent with the established SymPy-only non-status-only precedent (stages 121/122/123). I could not point to any specific claimed result that a corrected SymPy fails to settle, which is the bar the prompt sets for a valid `missing_mathematica` finding.

## Verdict justification

The mathematics holds up: I hand-verified the odd ratio (lines 64–67), the Δ_norm collapse and its Δ_Q form (lines 90, 94–97), the point-particle reductions (lines 111–116), the DtN-deformation χ_Q/N_Q/Δ_norm^pt forms (lines 139–143), the first-order linearization (lines 160–163, numerator first-order `15S·eps_beta + dSigma0 + 27dSigma5` over `3S` gives `−5eps_beta − dSigma0/(3S) − 9dSigma5/S`), and all canonical-branch specializations (lines 174–183) — these are genuine, non-tautological, and match the notes. Symbol domains are physically sound (no `symbol_assumption_error`); denominators `3S−Σ₀` and `Sβ⁵+9Σ₅` stay generically nonzero. The defect is that the stage's named headline — the boxed `\stagefield{Output}` `m̂₀²χ_Q N_Q=1` — is verified only by a self-comparison (`odd_closure − odd_closure`, line 89) instead of being derived from the observable odd condition `m̂₀²Γ̄₅=Γ̄₅^target` plus the verified ratio (F1, medium), with two further low-severity definitional echoes (F2). The fix is a small, mechanical SymPy edit that composes existing in-script quantities, introducing no new constant and no paper drift. Hence `verdict: findings`, `paper_alignment: partial`, no stop_cold.

## Self-test notes

Variable-independence: checked the only `sp.diff`-equivalent — the `series(t,0,2)` linearization (lines 147–153); `t` is wired into β, Σ₀, Σ₅ via `.subs`, and `NQ_from_def` genuinely depends on all three, so the first-order coefficient is non-vacuous (confirmed against the hand value). Trivial-case/round-trip: F1's proposed residual reduces by hand to `Γ̄₅^target(m̂₀²χ_Q N_Q − 1)` using the in-script Γ̄₅ and P0_target, so it is zero only for the correct normalization, not by construction, and uses no constant absent from the paper (P0_target and Γ̄₅^target are the boxed notes values). Stale labels ("STAGE 178" banner, "Stage 194 deformation algebra") are cosmetic print strings, not assertions, so they do not gate the math; noted but not raised as a separate finding.
