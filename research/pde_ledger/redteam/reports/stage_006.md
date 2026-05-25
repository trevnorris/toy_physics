---
unit_id: 006
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 006 red-team report (v2 paper-grounded re-audit)

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_006.tex`
- notes: (none) — em_projected per-stage notes (notes/em_projected/step_NN_*.md) were not committed to the repo per orchestrator note
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 34 references stage 006)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.txt`

Also consulted (for sign-convention inheritance only): `paper/stages/stage_005.tex`, the parent law `eq:stage005-projected-maxwell-expanded` referenced from the stage 006 card.

Freshness check: sympy script mtime epoch 1779383264; sympy output 1779384369 (newer, fresh). Mathematica script 1779348128; Mathematica output 1779385859 (newer, fresh). `outputs_fresh: true`.

## What the paper claims

The stage 006 card defines two field layers built from the projection kernel `W`: the measured fields `E_meas = ∫ W E dw, B_meas = ∫ W B dw` (eq:stage006-measured-fields) and the source-coupled flux fields `D_flux = ∫ W Z E dw, H_flux = ∫ W Z B dw` (eq:stage006-flux-fields). It then states the projected vector laws under the convention `F^{i0} = E_i`: the homogeneous pair `∂_t B_meas + ∇×E_meas = 0` and the inhomogeneous Ampere law `∇×H_flux − ∂_t D_flux = μ₀ J_proj + L_mix` (eq:stage006-ampere), where `L_mix` is the vector form of the transverse leakage term carried forward from stage 005 (eq:stage005-projected-maxwell-expanded). The card's Output paragraph (verbatim): "Stage 006 exports the field split (eq:stage006-measured-fields)–(eq:stage006-flux-fields) and the vector-law sign convention used by the rest of the projected Maxwell block." The appendix row (line 34 of stage_appendix_part01.tex) summarizes the deliverables as "Measured-field homogeneous equations, source-flux inhomogeneous equations, and the open-system charge balance." No "gauge-driver" symbolic term appears in the stage 006 card or its inherited stage 005 equations; the only inhomogeneity beyond `μ₀ J_proj` is `L_mix`.

## What the script claims to verify

Both scripts verify: (a) the homogeneous Bianchi rearrangement `div B = 0` and `∂_t B + ∇×E = 0` from the antisymmetric tensor identities (symbolic); (b) the inhomogeneous projected rearrangement in vector form, which the script writes as `div D + Leak0 + Gauge0 = μ₀ ρ_proj` (Gauss) and `∇×H − ∂_t D + Leak_vec + Gauge_vec = μ₀ J_proj` (Ampere), again symbolic; (c) a concrete projection audit on a fixed bulk potential `A^μ = ...·(1+w²)` with `W(w) = exp(−w²)/√π`, `Z(w) = exp(−w²)`, exercising the projected Bianchi, Gauss, and Ampere identities numerically and the leakage IBP relation; (d) the leakage normalization `leak1 = √2/4 = 1/(2√2)` from `−∫ W'(w) Z(w) w dw`, which both engines independently compute and agree on; (e) sign-mutation guards (Faraday, IBP) that flip a sign and confirm the assertion would fail; (f) projection-discrimination tests (bogus projection mixing in `w`, asymmetric weight on a non-potential 2-form, antisymmetric-Z mediator killing the leak). The SymPy script docstring (lines 7–12) and comment (line 57) frame the extra symbols `Leak_μ` (transverse leakage) and `Gauge_μ` (gauge-driver) as additive open-system terms.

## Paper ↔ script cross-check

| Paper deliverable | Script-side coverage | Status |
|---|---|---|
| Measured fields `E_meas, B_meas` (eq:stage006-measured-fields) | SymPy lines 42–47 declare measured-field symbols; concrete projection L241–242 builds `E_bulk_proj, B_bulk_proj` from the bulk potential via `Pg`. Mathematica L153–154 (`Eprojected`, `Bprojected`). | match |
| Flux fields `D_flux, H_flux` (eq:stage006-flux-fields) | SymPy lines 50–55 declare flux symbols; L272–273 build `D_bulk_proj, H_bulk_proj` via `Pg(Z·...)`. L225 explicitly checks measured ≠ flux under nontrivial Z. Mathematica L155–156. | match |
| Homogeneous Faraday `∂_t B + ∇×E = 0` | SymPy L98–111 (symbolic), L244–255 (concrete); Mathematica M1 (L82–97), M4 (L165) | match |
| Homogeneous `div B = 0` | SymPy L82–87 (definitional print only), L243 (concrete assertion); Mathematica M1 (L93–98), M4 (L164) | match |
| Inhomogeneous Ampere `∇×H − ∂_t D = μ₀ J + L_mix` (eq:stage006-ampere) | SymPy L150–163 (symbolic), L299–310 (concrete); Mathematica M2 (L120–125), M4 (L170–174). Sign-convention bookkeeping: the script puts the leakage on the LHS with opposite sign (`+ Leak`); this is algebraically equivalent to the paper's `= μ₀ J + L_mix` with `L_mix = −Leak`. SymPy's `leak = −Pgp(Z F^{w,μ})` while the paper's `L_mix` traces back to `+∫(∂_w W) Z F^{wμ} dw = +Pgp(Z F^{wμ})` via stage 005 with boundary=0. Equivalent. | match |
| "Open-system charge balance" (appendix row) — Gauss law `div D + leak0 = μ₀ ρ_proj` | SymPy L135–139 (symbolic), L295–298 (concrete); Mathematica M2 (L119, L124), M4 (L166–169). Not redisplayed in the body of the stage 006 .tex but inherited from stage 005's general-μ projected law. | match (anchored via stage 005 inheritance) |
| Leakage `L_mix` (vector form of stage 005 transverse leakage) | SymPy L259, L264–265 (`leak1 = √2/4`); Mathematica M3 (L134, `vectorLeakMoment = 1/(2√2)`) — both engines agree numerically | match |
| Sign convention `F^{i0} = E_i` | SymPy L77–79 hardcode `F10=E1, F20=E2, F30=E3, F01=−E1, ...`; Mathematica L77–80 define `fieldF[i,0] = Evec[[i]]`, `fieldF[0,i] = −Evec[[i]]` | match |
| "Gauge-driver" terms `G_μ` (script-side only) | SymPy L66–69, L135, L142–144; Mathematica L104, L117. Not exercised by any concrete check (cancel between LHS and RHS in every assertion they appear in). No paper-side counterpart. | extra (paper does not mention a separate gauge-driver term) |

Dominant pattern: every paper-side deliverable has a faithful script-side check on both engines, plus the script carries an extra symbolic placeholder (`Gauge_μ`) that the paper does not reference. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 110–111 | `sp.simplify(far_i − (∂_t B_i + curlE_i)) == 0` | Faraday rearrangement (homogeneous) | yes (algebraic identity from antisymmetric tensor; the Mathematica M1 independently re-derives it via Levi-Civita contraction so the sign convention itself is cross-verified) |
| A2 | sympy | 162–163 | `sp.simplify(lhs_i − amp_target_i) == 0` | Ampere rearrangement (inhomogeneous) | yes; Gauge_μ enters trivially and is not anchored (see F1) |
| A3 | sympy | 225 | `simplify(E_meas − D_flux) != 0` for concrete `Z` | (E,B) ≠ (D,H) — paper's "two pairs" claim | yes |
| A4 | sympy | 243 | concrete `divB == 0` | homogeneous div B | yes (non-tautological: built from explicit potential, projection actually integrated) |
| A5 | sympy | 244–255 | concrete `∂_t B + ∇×E == 0` (3 components) | homogeneous Faraday | yes |
| A6 | sympy | 262 | `[Wg·ZFw1]_∂ == 0` | localized boundary discharge (stage 005 carry-forward) | yes |
| A7 | sympy | 263 | IBP identity with explicit boundary | leakage IBP from stage 005 | yes |
| A8 | sympy | 264 | `leak1 == Pg(∂_w(Z·Fw1))` | leakage definition vs IBP | yes |
| A9 | sympy | 265 | `leak1 == √2/4` | numerical leakage normalization | yes (non-tautological literal anchored to specific Z, W choice) |
| A10 | sympy | 266 | mutated IBP sign != 0 | sign-mutation guard | yes |
| A11 | sympy | 295–298 | concrete projected Gauss `div D + leak0 = μ₀ Pg(J⁰)` | inhomogeneous Gauss (open-system charge balance) | yes |
| A12 | sympy | 299–310 | concrete projected Ampere (3 components) | inhomogeneous Ampere | yes |
| A13 | sympy | 320–324 | bogus-projection divB fails | projection-must-respect-3D-coords guard | yes |
| A14 | sympy | 338–341 | non-potential 2-form, symmetric weight: divB = 0 | projection parity tracking | yes (informational) |
| A15 | sympy | 344–348 | non-potential 2-form, asymmetric weight: B₁ ≠ 0 | shows divB=0 was due to weight symmetry | yes |
| A16 | sympy | 361–365 | antisymmetric-Z mediator kills leak | parity-of-Z controls leak | yes |
| A17 | sympy | 366–370 | Gaussian-Z leak ≠ antisymm-Z leak | discriminates mediator parities | yes |
| A18 | sympy | 371–374 | mutated concrete Faraday sign != 0 | sign-mutation guard | yes |
| M1a | math | 97 | M1 Faraday rearrangement = 0 | homogeneous Faraday (independent Levi-Civita derivation) | yes |
| M1b | math | 98 | M1 divB rearrangement = 0 | homogeneous div B | yes |
| M2a | math | 124 | M2 Gauss rearrangement = 0 | inhomogeneous projected Gauss form | yes; Gauge term cancels trivially |
| M2b | math | 125 | M2 Ampere rearrangement = 0 | inhomogeneous projected Ampere form | yes; Gauge term cancels trivially |
| M3a | math | 132 | M3 boundary discharge = 0 | boundary term for Gaussian weight | yes |
| M3b | math | 133 | M3 projected IBP relation = 0 | IBP identity | yes |
| M3c | math | 134 | `vectorLeakMoment − 1/(2√2) = 0` | numerical leakage normalization (engine cross-check vs A9) | yes |
| M4a | math | 164 | M4 projected Bianchi divB = 0 | concrete homogeneous div B | yes |
| M4b | math | 165 | M4 projected Bianchi Faraday = 0 | concrete homogeneous Faraday | yes |
| M4c | math | 166–169 | M4 projected Gauss law residual = 0 | concrete inhomogeneous Gauss | yes |
| M4d | math | 170–174 | M4 projected Ampere law residual = 0 | concrete inhomogeneous Ampere | yes |
| M5a | math | 179 | M5 IBP sign mutation ≠ 0 | sign-mutation guard | yes |
| M5b | math | 180 | M5 concrete Faraday sign mutation ≠ 0 | sign-mutation guard | yes |

All non-trivial assertions; both engines independently confirm `leak1 = √2/4 = 1/(2√2)`.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_006.tex:31-43` (the only inhomogeneity named in the stage 006 card is `L_mix`; no "gauge-driver" symbol appears)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:66-69, 135, 142-144, 166-167, 193-195` (introduces and carries `Gauge0..Gauge3`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl:104, 117` (parallel `gauge` symbols in `projectedInhom[nu]`)

**What's wrong:**
Both scripts carry a set of symbolic placeholder terms `Gauge_μ` (SymPy) / `gauge` (Mathematica) as additive contributions to the projected inhomogeneous equations. The SymPy docstring (lines 7–12) describes them as "gauge-driver terms" and the inline comment at line 57 says "Projected source, leakage, and gauge-driver terms." However, the paper card for stage 006 (and the upstream stage 005 card it inherits from) name only the transverse leakage `L_mix` / `−[WZF^{wμ}]_∂ + ∫(∂_w W) Z F^{wμ} dw` as the open-system inhomogeneity. There is no separate `G_μ` / gauge-driver term on the paper side.

Paper quote (`paper/stages/stage_006.tex:38-43`):
```
\nabla\times\mathbf H_{\rm flux}-\partial_t\mathbf D_{\rm flux}
=\mu_0\mathbf J_{\rm proj}+\mathbf L_{\rm mix}.
```
followed by: "Here `\mathbf L_{\rm mix}` is the vector form of the transverse leakage term in eq:stage005-projected-maxwell-expanded."

Script quote (`scripts/...sympy_audit.py:165-167`):
```
print("Compact inhomogeneous vector form:")
print("  div D + Leak0 + Gauge0 = mu0 rho_proj")
print("  curl H - partial_t D + Leak_vec + Gauge_vec = mu0 J_proj")
```

Mathematica analogue (`mathematica/...mathematica_audit.wl:114-117`):
```
projectedInhom[nu_Integer] := Sum[
  D[fluxG[mu, nu], braneCoords[[mu + 1]]],
  {mu, {0, 1, 2, 3}}
] + leak[[nu + 1]] + gauge[[nu + 1]];
```

Because `Gauge_μ` is added equally on both sides of every assertion that references it (e.g. SymPy `lhs1 = ... + L1 + G1` vs `amp1_target = ... + L1 + G1` at lines 142, 150; Mathematica's M2 rearrangement at lines 119, 121 likewise has the `gauge` symbol on both sides of the residual), the symbol cancels trivially. No concrete-projection check (Section 5 of SymPy, M4 of Mathematica) ever uses `Gauge_μ` — those checks use only the leak term and the bulk current.

**Why this matters:**
Low-severity documentation drift. The scripts' Section 4 / Compact-summary form `div D + Leak + Gauge = μ₀ ρ` and `curl H − ∂_t D + Leak + Gauge = μ₀ J` advertise an extra inhomogeneity term that the paper card does not derive or even mention. A reader following the paper would wonder whether the scripts are exercising additional physics that has not been documented. Conversely, if the gauge-driver terms are intended (e.g., as placeholders for projected gauge-fixing contributions that could appear in extended derivations), the paper card needs to introduce them, since the script's Section 4 summary makes them look load-bearing. Either way, the paper and the scripts should agree on what symbols appear in the projected inhomogeneous equation.

**Required change:**
See "Resolve before fix_loop" block in the directive. Codex must not edit this autonomously.

**Verification:**
After user-chosen direction is applied, the script's Section 4 compact summary and the paper card's body equations should reference exactly the same set of inhomogeneity symbols (`Leak` alone, or `Leak` + `Gauge`).

## Independent-derivation check (Mathematica)

The Mathematica script is not a transliteration of the SymPy script. Concrete evidence:

- SymPy enumerates components by hand: `F10, F20, F30 = E1, E2, E3`, `F23, F31, F12 = B1, B2, B3` (lines 77–79), and writes the Faraday cyclic-Bianchi components explicitly as `far1 = sp.diff(F23, t) + sp.diff(F30, y) + sp.diff(F02, z)` etc. (L90–92). Curl components are unrolled: `curlE1 = sp.diff(E3, y) - sp.diff(E2, z)` (L94–96).
- Mathematica defines a generic antisymmetric tensor function `fieldF[i_Integer, j_Integer] /; 1 <= i <= 3 && 1 <= j <= 3 := Sum[eps3[[k, i, j]] Bvec[[k]], {k, {1,2,3}}]` (L77–80) and computes the cycle residual using `Sum[eps3[[i,j,k]] (D[fieldF[j,k], t] + ...), ...]` (L82–92). The curl is generic: `curl3[v_] := Table[Sum[eps3[[i,j,k]] D[v[[k]], spaceCoords[[j]]], ...], {i, 1, 3}]` (L61–64).
- For the inhomogeneous sector, Mathematica defines a parameterized `projectedInhom[nu_]` that sums `D[fluxG[mu, nu], braneCoords[[mu + 1]]]` over all four indices (L114–117), then forms the rearrangement residual. SymPy enumerates `lhs0, lhs1, lhs2, lhs3` by hand (L135, L142–144).
- For the concrete bulk-potential audit, Mathematica defines `potential = {...}` as a list and uses `twoForm[a, b]` as a generic function (L137–146); the projected fields are then `Eprojected = Table[Pg[twoForm[i, 0]], {i, 1, 3}]` (L153). SymPy builds the corresponding components by hand, one at a time (L229–242).

These are different algebraic styles for the same physical claim; they re-derive the projected Maxwell rearrangements independently. Not transliteration.

## Engine cross-check

The two engines independently compute the leakage normalization:

- SymPy: `assert_zero("nonzero vector leakage value", leak1 - sp.sqrt(2) / 4)` (`leak1 = √2/4`, with `leak1 = -Pgp(Zg * Fw1)` and `Fw1 = w`).
- Mathematica: `expectZero["M3 leakage normalization", vectorLeakMoment - 1/(2 Sqrt[2])]` (`vectorLeakMoment = 1/(2√2)`, with `vectorLeakMoment = -Pgp[ZZ[w] w]`).

Since `√2/4 = 1/(2√2)`, the engines agree exactly. All other engine-comparable checks (divB, Faraday, Gauss, Ampere on the concrete bulk potential) yield residual 0 in both engines (confirmed in the saved outputs).

The IBP sign-mutation guards: SymPy expects `projected_transverse_derivative - (boundary + Pgp)` to be nonzero (line 266); Mathematica expects `Pg[D[gaussianSample, w]] - (gaussianBoundary + Pgp[gaussianSample])` to be nonzero (line 177) and reports residual `1/Sqrt[2]`. Both engines correctly detect the sign mutation.

## Verdict justification

The scripts pass the v2 paper-grounded audit on every load-bearing claim of the stage 006 card: measured/flux field definitions, homogeneous Faraday and `div B`, inhomogeneous Ampere and Gauss with leakage, the `F^{i0}=E_i` sign convention, and the numerical leakage normalization. The Mathematica script independently re-derives the same identities (using generic Levi-Civita / Sum-over-indices machinery rather than enumerated components) and agrees numerically. The mutation tests (sign flips on Faraday and on IBP) demonstrate that the assertions would fail under the wrong sign convention, confirming the checks are non-tautological.

The only discrepancy is the scripts' inclusion of a `Gauge_μ` placeholder symbol that the paper card does not introduce. Because the placeholder cancels in every assertion it appears in, it does not affect any verified physics — it is documentation drift, not a math error. Verdict: `findings` with a single low-severity `paper_misalignment` (subtype `paper_missing_script_claim`), routed to the user for direction (either drop the symbol from the scripts, or introduce/justify it in the paper card).

Attacks tried that all failed (scripts held up):

- Re-derived the leakage `−∫ W'(w) exp(−w²) w dw` by hand → `√2/4`, agrees with both engines.
- Checked the sign convention chain: paper Ampere `= μ₀ J + L_mix` vs script `+ Leak + Gauge = μ₀ J`. With `L_mix = +Pgp(...)` (from stage 005 eq:stage005-projected-maxwell-expanded with boundary=0) and `Leak = −Pgp(...)` (from script line 259), the two forms are algebraically equivalent. No sign error.
- Checked whether the symmetric Gaussian weight `Wg = exp(−w²)/√π` projects out the bulk potential's `(1+w²)` factor non-trivially → `Pg(1+w²) = 3/2`, nonzero. The concrete checks really do exercise the integration.
- Checked whether the Mathematica `bulkCurrent[nu]` sums over all five bulk indices including `w` (otherwise the Gauss/Ampere checks would be trivially zero). Line 160: `Sum[..., {mu, {0, 1, 2, 3, 4}}]` — yes, all five. ✓
- Checked whether the SymPy `J0_bulk` includes `∂_w(Zg Fw0_bulk)` (line 284) — yes. ✓
- Tried to identify a tautological assertion: A1, A2 are tensor-rearrangement identities (true by construction of antisymmetric F), but they are still useful as type-checks of the vector form and are independently re-derived by Mathematica via Levi-Civita contraction; the concrete checks A11–A12 and M4c–M4d provide substantive verification by integrating against an explicit bulk potential.
- Confirmed the v1 audit's earlier insufficient_verification concern was addressed: the script now includes assertions A13–A18 (bogus projection, non-potential 2-form with asymmetric weight, antisymmetric-Z mediator) that isolate projection physics from F=dA-trivial identities.

## Self-test notes

(1) Variable independence — the bulk potential `A^μ = ...·(1+w²)` actually depends on `(t, x, y, z, w)`, so every `sp.diff(..., x)` and `D[..., t]` is meaningful. The bogus-projection test `B3_bogus = B3_bulk_proj + z * w**2` produces a `z`-dependent term, so `∂_z B3_bogus` is nonzero — correctly catches the leak. (2) Parity — `Wg(w) = exp(−w²)/√π` is even, `Wgp(w)` is odd, `Z(w) = exp(−w²)` is even, `Fw1 = w` is odd; the product `Wgp · Z · Fw1 = odd·even·odd = even`, integrable on `(−∞,∞)` and nonzero — consistent with `leak1 = √2/4`. The antisymmetric-Z test uses `Z=w` (odd), so `Wgp · w · w = odd·even = odd`, integral = 0 — consistent with the assertion. (3) Trivial-case substitution — substituting the concrete bulk potential `A^0 = x z (1+w²)` etc. into the projected Gauss law: `D1 = Pg(Zg · F10_bulk) = Pg(exp(−w²) · z(1+w²)) = z · Pg(exp(−w²)(1+w²))` (using `Pg(c·f(w)) = c·Pg(f)` for `c` independent of `w`), and the residual checks then reduce to algebraic identities that simplify cleanly. (4) The single finding's resolution does not require any new derivation, just a paper-vs-script symbol agreement; it routes to the user.
