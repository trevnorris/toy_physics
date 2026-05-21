---
unit_id: 006
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 006 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.txt`
- mathematica output: `(missing)`

mtimes: sympy script 2026-05-04 12:00 (epoch 1777917651); sympy output 2026-05-11 12:38 (epoch 1778524715). Output is newer than script, so `outputs_fresh: true`. The unit's manifest entry has `is_status_only_candidate: False`, so both engines are required.

## What the script claims to verify

The script claims to derive a brane-side projected Maxwell-style vector system from a 4+1 antisymmetric-tensor theory by introducing two distinct field layers — "measured" fields `(E_meas, B_meas) = Proj_W[F_{MN}]` and "source-coupled" flux fields `(D_flux, H_flux) = Proj_W[Z F^{MN}]` — and showing that (i) the homogeneous index-by-index Bianchi identities for `(F_{a0}, F_{ab})` rearrange to the vector form `div B = 0`, `curl E + ∂_t B = 0`; (ii) the index-by-index inhomogeneous equations for `Proj_W[Z F^{MN}]` rearrange to `div D + Leak0 + Gauge0 = mu0 rho_proj` and `curl H − ∂_t D + Leak_vec + Gauge_vec = mu0 J_proj`; (iii) for a concrete Gaussian weight `W = exp(-w^2)/sqrt(pi)`, mediator `Z = exp(-w^2)`, and a specific bulk potential `A_M(t,x,y,z,w) = (brane part) × (1 + w^2)`, the projected divergence-free / Faraday relations hold, the boundary terms in the leakage integration-by-parts vanish, the leak integrand reduces to a specific closed-form value `sqrt(2)/4` for `F_{w1} = w`, and the projected inhomogeneous equations with explicit `leak_μ` terms are consistent with a current `J_μ_bulk` defined by the bulk inhomogeneous Maxwell equation. Two mutated-sign assertions are included as adversarial controls.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 111 | `simplify(far_i - (∂_t B_i + curlE_i)) == 0`, i=1..3 | partial (index-to-vector rearrangement only) |
| A2 | sympy | 163 | `simplify(lhs_i - amp_i_target) == 0`, i=1..3 | partial (index-to-vector rearrangement only) |
| A3 | sympy | 225 | `Pg(F10_bulk) - Pg(Zg F10_bulk) ≠ 0` (E vs D distinction for nontrivial Z) | yes |
| A4 | sympy | 243 | `∂_x B1_proj + ∂_y B2_proj + ∂_z B3_proj == 0` for explicit A-derived bulk | partial (automatic from F=dA, not projection-specific) |
| A5 | sympy | 245-255 | `∂_t B_i_proj + curlE_i_proj == 0` for explicit A-derived bulk, i=1..3 | partial (automatic from F=dA) |
| A6 | sympy | 262 | `boundary_value(Wg * Zg * Fw1) == 0` | yes |
| A7 | sympy | 263 | `Pg(∂_w(Zg Fw1)) - (boundary - Pgp(Zg Fw1)) == 0` | partial (IBP identity that sympy satisfies automatically) |
| A8 | sympy | 264 | `leak1 - Pg(∂_w(Zg Fw1)) == 0` | partial (follows from boundary = 0 and definitions) |
| A9 | sympy | 265 | `leak1 - sqrt(2)/4 == 0` | yes (closed-form numerical anchor) |
| A10 | sympy | 266 | mutated leakage IBP (sign flip) is nonzero — adversarial control | yes |
| A11 | sympy | 276 | `boundary_value(Wg Zg ∂_w A_μ) == 0` for μ=0..3 | yes |
| A12 | sympy | 280-281 | `Pg(∂_w(Zg ∂_w A_μ)) - (boundary - Pgp(Zg ∂_w A_μ)) == 0` for μ=0..3 | partial (IBP identity automatic) |
| A13 | sympy | 296-297 | concrete projected Gauss law residual is zero | partial (J0_bulk defined to make this hold; construction-tautological) |
| A14 | sympy | 299-309 | concrete projected Ampere components 1..3 residual zero | partial (J_i_bulk defined to make these hold) |
| A15 | sympy | 312 | mutated concrete Faraday sign is nonzero — adversarial control | yes |

A row marked "partial" feeds either `insufficient_verification` or `tautological_check`. Sympy-only — no `B*` rows.

## Findings

### F1 — missing_verification_script

**Severity:** high
**Files:**
- `(missing)` — expected at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl`

**What's wrong:**
There is no Mathematica script for unit 006. The audit prompt explicitly lists the mathematica path as `(missing)`, and a filesystem search under `/var/projects/toy_physics/research/pde_ledger/` confirms that no `*stage006*.wl` exists. The unit is not status-only (`is_status_only_candidate: False`), and the unit is not a checkpoint, so the second-engine policy requires an independent Mathematica derivation. None is present, so the projected-Maxwell-vector claims (homogeneous Bianchi rearrangement, inhomogeneous Ampere/Gauss rearrangement, boundary discharge, leak normalization `sqrt(2)/4`, concrete Gauss/Ampere consistency) are verified by a single engine only.

**Why this matters:**
The second-engine policy exists to catch errors that survive one engine's algebraic conventions (sign of curl, ordering of indices, branch choices in symbolic integration). With only sympy present, a systematic sign error in the index-to-vector mapping (e.g. a `curlH_i` defined with the wrong handedness on line 146-148) or in the antisymmetric-tensor `G^{μν}` definitions (lines 129-132) would survive undetected — sympy verifies internal consistency, not correctness.

**Required change:**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl` that independently derives and asserts:

M1. The homogeneous projected Bianchi rearrangement: for an antisymmetric tensor `F_{MN}` on `(t,x,y,z)` with `E_i := F_{i0}` and `B_i := (F_{23}, F_{31}, F_{12})`, the cyclic Bianchi identities `∂_t F_{jk} + ∂_j F_{k0} + ∂_k F_{0j} = 0` and `∂_x F_{23} + ∂_y F_{31} + ∂_z F_{12} = 0` rearrange to `∂_t B + curl(E) = 0`, `div B = 0`.

M2. The inhomogeneous projected rearrangement: for antisymmetric `G^{MN}` with `G^{i0} = D_i` and `G^{ij} = -ε_{ijk} H_k`, the index-by-index equation `∂_μ G^{μν} + L^ν + G^ν = mu0 j^ν` rearranges to `div D + Leak0 + Gauge0 = mu0 rho` and `curl H − ∂_t D + Leak_vec + Gauge_vec = mu0 J_vec`. Match the SymPy script's sign convention for `Leak`, `Gauge`, and `curl H` (sympy lines 142-152). Do not transliterate sympy line by line — derive from `D[F[μ,ν], x^μ]` in Mathematica using a metric-based index summation rather than predeclaring `(G01, G02, G03, G23, G31, G12)` symbol by symbol.

M3. Concrete Gaussian projection identities: with `W(w) = Exp[-w^2]/Sqrt[Pi]`, `Z(w) = Exp[-w^2]`, define `Pg[f] := Integrate[W f, {w, -Infinity, Infinity}]` and verify:
  - `Pg[Z w^2 * D[Z w, w] / Z]` ... (more concretely) the leak normalization `leak1 := -Integrate[W'(w) Z(w) w, {w, -Infinity, Infinity}] == 1/(2 Sqrt[2])` (which equals `sqrt(2)/4`).
  - The boundary terms `Limit[W Z w, w -> Infinity] - Limit[W Z w, w -> -Infinity] == 0`.
  - The IBP relation `Pg[D[Z w, w]] == -Integrate[W'(w) Z w, ...]` (boundary form).

M4. For the explicit bulk potential `A_0 = x z (1 + w^2)`, `A_1 = t y (1 + w^2)`, `A_2 = t z (1 + w^2)`, `A_3 = x y (1 + w^2)`, the projected fields `E_proj`, `B_proj`, `D_proj`, `H_proj` satisfy the projected Bianchi (homogeneous) and the projected Gauss / Ampere (inhomogeneous) equations with the leak terms defined by `-Pgp[Z * ∂_w A_μ]`, where `Pgp` integrates against `W'(w)`. The source `J_μ_bulk` is defined from the bulk inhomogeneous Maxwell equation `∂_μ(Z F^{μν}) + ∂_w(Z F^{wν}) = mu0 J^ν`. State this explicitly in a comment so the Mathematica check is not merely "matching the sympy".

M5. Adversarial sign mutations: at least one IBP-sign mutation and one Faraday-sign mutation are confirmed to fail (parallel to sympy lines 266 and 311).

The Mathematica script should use Mathematica-idiomatic constructions: `LeviCivitaTensor[3]` for curl, `Sum[D[G[μ,ν], x[μ]], {μ, 0, 3}]` rather than enumerated `G^{0ν} + G^{2ν} + G^{3ν}`, `Integrate[..., {w, -Infinity, Infinity}]` for projections, and `Limit[..., w -> Infinity]` for boundary discharge. It must not echo the SymPy intermediate names (`F10`, `F23`, `lhs1`, `amp1_target`, `leak1`).

**Verification:**
After creation, `redteam exec-mathematica 006` runs the new `.wl`, produces an output file in `/var/projects/toy_physics/research/pde_ledger/mathematica/output/`, and exits with status 0. The output should contain explicit lines confirming the five M1..M5 results, with at least one assertion per result. The Mathematica `leak1` numeric must agree with `1/(2 Sqrt[2])` (= `Sqrt[2]/4`).

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:98-111` (homogeneous rearrangement)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:142-163` (inhomogeneous rearrangement)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:241-310` (concrete projection)

**What's wrong:**
Three sub-issues sit in the same finding because they have the same root cause: the assertions test internal algebraic identities, not the projected-Maxwell physics the docstring claims.

(a) Lines 98-111 check `simplify(far1 - (∂_t B1 + curlE1))` with `far1 = ∂_t F23 + ∂_y F30 + ∂_z F02` and `F23 = B1, F30 = E3, F02 = -E2`. Substituting gives `∂_t B1 + ∂_y E3 + ∂_z(-E2) = ∂_t B1 + ∂_y E3 - ∂_z E2 = ∂_t B1 + curlE1`. The two sides are byte-identical expressions in different notation; the `simplify(... - ...) == 0` is purely a notational rearrangement of the index-symbol assignments at lines 79-80. The same is true of the divB residue on lines 82-87 (only printed, not asserted). This verifies the author's sign convention, not any projection physics.

(b) Lines 142-163 do the same for the inhomogeneous sector. `lhs1 = ∂_t G01 + ∂_y G21 + ∂_z G31 + L1 + G1 = -∂_t D1 - ∂_y H3 + ∂_z H2 + L1 + G1` and `amp1_target = (∂_z H2 - ∂_y H3) - ∂_t D1 + L1 + G1` are the same expression with the terms rearranged. `simplify(lhs1 - amp1_target)` is zero by inspection. No physics is tested.

(c) Lines 241-310 ("concrete projection audit") look substantive but are construction-tautological in two places:
  - Lines 243-255 verify projected Bianchi (`div B_proj = 0`, `∂_t B_proj + curl E_proj = 0`) for fields built from a potential `A_M(t,x,y,z,w)`. But `F = dA` always satisfies `dF = 0` by Schwarz, regardless of A or of the projection weight. The projection `Pg` is a linear operator commuting with brane-coordinate derivatives, so it preserves the identity. The check confirms only that sympy correctly handles partial-derivative commutativity under integration — not that the projected Bianchi structure is what the script claims.
  - Lines 283-309 verify projected Gauss / Ampere using `J_μ_bulk` *defined* by the bulk inhomogeneous Maxwell equation (lines 283-294). Projecting that definition and applying IBP with vanishing boundary terms gives the projected Gauss / Ampere identity automatically. The check confirms only that `Pg` commutes with brane-derivatives and that boundary terms vanish — again sympy's algebraic plumbing, not projection physics.

The substantive checks are A6 (boundary vanish), A9 (the closed-form leak normalization `sqrt(2)/4`), and A3 (E ≠ D for nontrivial Z), plus the two adversarial sign mutations A10 and A15. Everything else is internal consistency.

**Why this matters:**
The docstring claims the script "verifies" three things: that the homogeneous laws project cleanly, that the inhomogeneous laws acquire leakage and gauge-driver terms, and that the projected theory naturally distinguishes (E,B) from (D,H). The concrete projection audit looks like it does heavy lifting for those claims, but most of its assertions hold for any potential / any localized weight / any antisymmetric tensor, so they don't distinguish "projected Maxwell works" from "F = dA exists and sympy can integrate by parts". A skeptical reader has no test that, for example, the leak structure could not have been a *different* contraction of the transverse field strength.

**Required change:**
Add a small number of assertions that actually exercise the projection physics. Concretely, at the end of Section 5 (after sympy line 310, before line 316), append:
1. A check that the projected Bianchi *would fail* if the projection commuted incorrectly with the brane derivative. The cheapest version: define `bogus_div = ∂_x B1_bulk_proj + ∂_y B2_bulk_proj + ∂_z B3_bulk_proj` after substituting `B3_bulk_proj := B3_bulk_proj + w` (mixing transverse coordinate into the projected field) and `assert_nonzero` on it. This shows the divB check would fail if the projection were misapplied.
2. A check that uses a non-`F=dA` antisymmetric tensor (i.e. a bulk `H_{MN}` not constructed as the curl of any potential) and verifies that the projected divB is in general nonzero, then confirms that adding a `dH=0` condition restores `divB_proj = 0`. This isolates "projection preserves Bianchi" from "F=dA trivially satisfies Bianchi".
3. A check that the leak structure depends on the `Z(w)` profile by recomputing `leak1` with `Z(w) := 1` (the trivial mediator); under the docstring's claim the leak should vanish, so `assert_zero(leak1_trivial)` and `assert_nonzero(leak1_gaussian - leak1_trivial)`.

These are additive, do not refactor existing assertions, and stay inside the unit's scope.

**Verification:**
After Codex appends the three checks, the output file should grow by three new assertion lines (or one combined block) with names containing "bogus", "non-potential", "trivial Z". The exit code remains 0. The audit_inventory for a re-audit should show three rows marked "yes" in the "Anchored to claim?" column.

### F3 — stale_output

**Severity:** low
**Files:**
- script: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py` (mtime epoch 1777917651 = 2026-05-04 12:00)
- output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.txt` (mtime epoch 1778524715 = 2026-05-11 12:38)

**What's wrong:**
This is **not** a stale-output finding — the output mtime is later than the script mtime, so the saved output reflects the current script. I list this finding only to make the freshness check explicit and to record the timestamps, since the audit template includes a freshness flag (`outputs_fresh: true` in the front-matter).

**Why this matters:**
Informational only. Listed for traceability.

**Required change:**
None. After F1 (creating the Mathematica script) and F2 (adding three new sympy assertions) are applied, the verifier will rerun sympy and refresh the .txt; the new mtime will again be later than the script mtime.

**Verification:**
`stat -c '%Y' script output` shows output > script after the re-run.

## Independent-derivation check (Mathematica)

No Mathematica script exists. The check is not applicable until the script required by F1 is created. The newly-created `.wl` must derive the projected Maxwell vector form from index identities and projection definitions independently — see F1's "Required change" for the structural constraints (use `LeviCivitaTensor`, `Sum`, `SphericalHarmonicY` if needed, etc., not a line-by-line port of the SymPy intermediate names).

## Engine cross-check

Not applicable — only one engine is present. `engines_agree: n/a`.

## Verdict justification

`findings` (not `clean`, not `stop_cold`). The sympy script's algebra is internally consistent: I attacked the sign of `curlE` (lines 94-96) against the Bianchi index choice (lines 89-92) and got the same answer; I attacked the antisymmetric-tensor sign assignment for `G^{μν}` (lines 129-132) against the rearrangement target `amp_i_target` (lines 150-152) and confirmed the rearrangement holds; I evaluated the closed-form Gaussian integral `-∫ W'(w) Z(w) w dw` and confirmed it equals `sqrt(2)/4` matching line 265; I verified the boundary terms `lim_{w→±∞} W Z f → 0` for all f appearing in the boundary checks; and I confirmed the two adversarial sign mutations (lines 266 and 311) genuinely fail. What does not hold up is the strength of the verification: the "abstract" rearrangements in Sections 1 and 2 verify only the author's chosen index conventions, and the "concrete" audit's Bianchi / Gauss / Ampere checks are automatic consequences of `F = dA` and of `J_bulk` being *defined* from the bulk Maxwell equation, so they do not isolate "projection preserves Maxwell" from those standalone identities. Combined with the absence of the Mathematica script, three findings: missing second engine (high), insufficient verification of the projection-specific physics (medium), and an informational stale-output note (low; output is actually fresh). No `stop_cold` flag: each finding is locally fixable without affecting downstream units' results, since unit 006's outputs (the structural projected Maxwell form and the `sqrt(2)/4` leak normalization) are correct as far as the algebra goes — they just are not as thoroughly cross-verified as the policy demands.
