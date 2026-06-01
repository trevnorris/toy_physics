---
unit_id: 192
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md]
  paper_appendix: present
---

# Audit unit 192 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_192.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 115, 265, 1139-1210, 1561)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.py`
- mathematica: (missing — SymPy-only stage)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage builds exact orbit/quotient projectors on the 8-dimensional log-ratio drift space ordered `(λ, c, γ, U, K_η, W, μ, T)`. `\stagefield{Output}`: *"Builds exact projectors \(Q_{\rm quot}\), \(O_{\rm orb}\), showing quotient failure lives only on \((T_U,K_\eta^{\rm eff},\mu_W)\)."* The notes (the authoritative source, more detailed than the terse card) enumerate five deliverables: (1) the exact pivot block `P_(T,K_η,μ)` selected from the carried Stage 238 monomial-drift matrix `M_*` with `det = 1+χ0_* > 0`; (2) the canonical quotient section `S = E P^{-1}` satisfying `M_* S = I_3`; (3) the complementary projectors `Q_quot = S M_*`, `O_orb = I_8 - Q_quot` with `Q^2=Q`, `O^2=O`, `QO=OQ=0`, `M_*O=0`, `M_*Q=M_*`; (4) the quotient-failure support result with closed forms `(Δ_T)_fail = q_tr/(1+χ0)`, `(Δ_Keta)_fail = -q_η`, `(Δ_μ)_fail = F_*/(1+χ0) q_tr + q_nt - q_η`, all five free coordinates vanishing; and (5) the orbit-lock equivalence `Q_quot Δx = 0 ⟺ M_* Δx = 0 ⟺ Δx ∈ ker M_*`, plus the orbit-law closed forms recovered by `O_orb`. The appendix (lines 1139-1208) restates `M_*`, the pivot block, projectors, failure components and equivalence verbatim, with the claim-status that Stages 186-192 are exact closures on the positive coherent reduced state space.

## What the script claims to verify

The SymPy script hardcodes the carried `M_*` (a Stage-238 input, not derived here), selects columns 7/4/6 (T/K_η/μ) to form the pivot block, inverts it, builds the embedding `E`, the section `S = E P^{-1}`, and the projectors `Q = S M_*`, `O = I_8 - Q`. It then asserts the full ledger of identities: pivot block and its inverse match the paper's closed forms AND satisfy `P^{-1}P = P P^{-1} = I`; `M_* S = I_3`; the six projector identities; that rows 0,1,2,3,5 of `Q` vanish (failure support); the closed-form `Δx_fail` and `Δx_orbit` components; the reconstruction `Δx = orbit + fail`; `M_* Δx_orbit = 0`; and (Section VI) for an independent `q`-packet, `M_*(Sq) = q`, `Q(Sq) = Sq`, `O(Sq) = 0`. All over fully symbolic parameters (chi0_*, deltaU_*, E_*, F_* and eight free real drift components).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `M_*` carried matrix (notes §1, app. eq:q-packet) | lines 43-49 hardcoded; matches notes/app entry-by-entry | match |
| Pivot block `P_(T,K_η,μ)` and `det=1+χ0` (notes §2) | lines 59-72: column-selected, compared to expected, det factored | match |
| `P^{-1}` (notes §2) | lines 74-86: computed `.inv()` + self-checks `P^{-1}P=I`, `PP^{-1}=I` | match |
| Section `S = E P^{-1}`, `M_*S=I_3` (notes §3) | lines 88-110 | match |
| Projectors `Q,O` and 6 identities (notes §4, app eq:projector-identities) | lines 112-124 | match |
| Failure support only on dependent triple (Output, notes §5) | lines 126-129 (rows 0,1,2,3,5 zero) + 141-157 closed forms | match |
| Orbit-law closed forms + free coords preserved (notes §6) | lines 159-178 | match |
| Orbit-lock equivalence `QΔx=0 ⟺ M_*Δx=0` (notes §7, app eq:orbit-lock-projector) | pinned by `Q=S M_*` construction + `M_*Q-M_*=0` (line 124) + Section VI q-packet checks | match |

Dominant pattern: aligned.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 69 | `Pdep - Pdep_expected == 0` | pivot block (column-select of M_*) | yes |
| A2 | sympy | 72 | `factor(det) - (1+chi0s) == 0` | `det = 1+χ0` | yes |
| A3 | sympy | 84 | `Pdep_inv - expected == 0` | `P^{-1}` closed form | yes |
| A4 | sympy | 85-86 | `P^{-1}P - I == 0`, `P P^{-1} - I == 0` | `P^{-1}` is true inverse (self-validating) | yes |
| A5 | sympy | 109 | `Sdep - expected == 0` | section closed form | yes |
| A6 | sympy | 110 | `M_* S - I_3 == 0` | `S` is right inverse | yes |
| A7 | sympy | 119-124 | `Q^2-Q`, `O^2-O`, `QO`, `OQ`, `M_*O`, `M_*Q-M_*` all `==0` | 6 projector identities | yes |
| A8 | sympy | 128-129 | `Q rows {0,1,2,3,5} == 0` | failure support on dependent triple | yes |
| A9 | sympy | 156 | `QΔx - Δx_fail_expected == 0` | failure closed forms | yes |
| A10 | sympy | 157 | `M_* Δx_fail - q == 0` | failure carries full q | yes |
| A11 | sympy | 177 | `OΔx - Δx_orbit_expected == 0` | orbit-law closed forms + free coords | yes |
| A12 | sympy | 178 | `M_* Δx_orbit == 0` | orbit piece ∈ ker M_* | yes |
| A13 | sympy | 179 | `Δx - orbit - fail == 0` | unique decomposition | yes |
| A14 | sympy | 188-190 | `M_*(Sq)-q`, `Q(Sq)-Sq`, `O(Sq)` all `==0` | equivalence via independent q-packet | yes |
| A15 | sympy | 194 | `M_* Δx_orbit_expected == 0` | paper's orbit-law lies in ker M_* (cross-check) | yes |
| A16 | sympy | 193 | `S·0 == 0` | (none — vacuous) | no (harmless) |

## Findings

None. See observations in Verdict justification for the one non-finding (cosmetic stage-number banner) and the one vacuous-but-harmless check (A16, `S·0=0`).

## Independent-derivation check (Mathematica)

No `.wl` exists; transliteration question is moot.

## Engine cross-check

Single engine; not applicable.

## Verdict justification

Clean. I read the paper card, notes, and appendix first, built the model, then attacked the script. Attacks tried and failed:

- **Tautology hunt on the `*_expected` literals.** Every hardcoded `expected` block (pivot, inverse, section, failure, orbit) is compared against a quantity *independently computed bottom-up from `M_*`* (column selection, `.inv()`, matrix products), not substituted back into its own defining expression. The comparisons would fail if column indices, the inverse, or the projector products were wrong, so they are substantive, not X−X==0. The truly self-validating checks (`P^{-1}P=I`, `M_*S=I_3`, `Q^2=Q`, `M_*Q=M_*`) confirm the closed forms independently of the literals.
- **Column-index sabotage.** Verified by hand that the basis order `(λ,c,γ,U,K_η,W,μ,T)` makes T=7, K_η=4, μ=6, and that `hstack(M[:,7],M[:,4],M[:,6])` and the `Edep` embedding both reproduce the paper's `(T,K_η,μ)` ordering exactly. Recomputed `det P = 1+χ0` and `P^{-1}` by hand; both match.
- **Headline equivalence under-verification.** The `⟺` is not merely printed: `Q=S M_*` (construction) gives one direction, and the asserted `M_*Q-M_*=0` (line 124) gives `M_*QΔx = M_*Δx = q`, so `QΔx=0 ⟹ q=0`; Section VI's independent q-packet (`M_*(Sq)=q`, `Q(Sq)=Sq`, `O(Sq)=0`) pins the other direction. The equivalence is genuinely settled.
- **Symbol-domain / branch attack.** `chi0_*,deltaU_*,E_*,F_*` are `positive,real`; the only branch quantity is `det=1+χ0>0`, exactly the paper's "constructive coherent branch" condition that justifies `Pdep.inv()`. Drift symbols are `real` only (correct for arbitrary log-ratios). No over-strong assumption hides a branch error; no special-case substitution — verification is fully symbolic over the entire parameter space.
- **Staleness.** Output (mtime May 11 12:48) is newer than the script (May 11 11:58), exit 0, every printed residual matrix/vector is identically zero.

On the missing-Mathematica question: I did NOT flag `missing_mathematica`. This stage is exact linear algebra over a rational-function field in four symbols. SymPy's exact matrix arithmetic drives every residual to an identically-zero matrix; there is no numerical approximation, integral, or transcendental simplification where a second CAS would plausibly catch a SymPy error. There is no claimed result that SymPy fails to genuinely verify — all five notes deliverables and both Output projectors are fully and exactly settled. Per the prompt's line-114 judgment and SymPy-only precedent (stages 121/122/123; non-status-only allowed), single-engine is acceptable here.

## Self-test notes

- **Variable independence:** No `sp.diff` appears; the script is matrix algebra only, so the diff-against-uninvolved-symbol trap does not arise.
- **Symmetry/parity:** No integrals; not applicable.
- **Trivial-case pre-check:** Mentally substituted the load-bearing identities (`P^{-1}P=I`, `det=1+χ0`, `M_*S=I_3`, the failure/orbit closed-form rows against the column-selected `M_*`) by hand and confirmed each reduces to the asserted zero. The one vacuous check is A16 (`S·0=0`, trivially true) — harmless redundant scaffolding, not load-bearing on any paper claim, so not filed as a finding. Noted cosmetic-only mismatch: the print banner says "STAGE 175"/"Stage 187"/"Stage 238" (pre-renumber IDs) while the unit/card is Stage 192; these strings feed no assertion and carry no math, so they are below the bar for a script-side finding and do not constitute paper_misalignment of any verified claim.
