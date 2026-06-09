---
unit_id: 167
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage167_bundle_transport_tangent_compensation.md]
  paper_appendix: present
---

# Audit unit 167 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_167.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage167_bundle_transport_tangent_compensation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 65, 334–361)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` is verbatim: "Proves arbitrary first-order isotropic bundle drift is tangent to the exact compensated parent family." The derivation ledger says the stage "removes tangent-compensated isotropic motion, proves no linear grouped feed-down into scalar slippages, and reduces the direct grouped outlet map to the physical even-response and prefactor slopes." The card is terse; the notes are authoritative on the algebra. The notes enumerate the full set of deliverables: (1) the bundle transport laws for the remaining mouth/background variables `c_{s,w}, ell, L_W, a, c_s, Z_q, v_{w0}, T_m, g_s, g_q, lambda` as explicit linear images of `(dΘ_w, dK_s, dK_q, dP_0)`; (2) the three parent invariants `δln r_c = δln 𝔯 = δln 𝔤 = 0`; (3) the exact vanishing of the two Stage-163 off-family channels and of the normal coordinate `δ⊥ = 0`; (4) the cancellation of `dP_0` from `δln v_{w0}`. The appendix (eq. `app-part05-tangent-parent-family`) restates only the invariant triple `δln r_c = δln 𝔯 = δln 𝔤 = 0`.

## What the script claims to verify

Both scripts define the four bundle drifts `dTheta, dKs, dKq, dP`, build the Stage-166 inversion primitives (`drho, da, dcs, dZ`) and the frozen-EOS / healing-lock identities (`dcsw=2drho, dell=-dcsw, dLW=da`), then assemble the Stage-165 branch laws `dv, dT`, the parent couplings `dgq, dgs, dlam` (via the intermediate `dIsq = 2da+dell+dLW/2`), and finally the invariants. The SymPy script `expect_zero`s only the three invariants (`drc, dr, dg`), the two channels, and `delta_perp`; it merely PRINTS the carry-forward laws. The Mathematica script additionally `expectZero`s the carry-forward laws as primitive-route-vs-boxed-result closures: `dv` equals `-3dKs/4+dKq/2+13dTheta/8`, `dT` equals `-5dKs/4+dKq/2+15dTheta/8-2dP/5`, `dgs`, `dgq`, `dlam` against their boxed forms, plus `(dv-dT)=(2dcs-da)` and `(dv+dT)=(dZ+5drho-4da)`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Output: invariance / tangency `δln r_c = δln 𝔯 = δln 𝔤 = 0` | `expect_zero(drc/dr/dg)` (py 77–79, wl 76–78) | match |
| δ⊥ = 0 and channel vanishing (notes §5) | `expect_zero(chan1/chan2/delta_perp)` (py 84–89, wl 83–87) | match |
| Bundle transport laws v_{w0}, T_m (notes §2 boxes) | wl `expectZero` carry-forward (wl 53–54); py prints only (py 56–57, 92–93) | match (wl) / partial (py) |
| Transport laws g_s, g_q, lambda (notes §3 boxes) | wl `expectZero` (wl 68–70); py prints only (py 67–70, 94–96) | match (wl) / partial (py) |
| dP_0 cancels from δln v_{w0} (notes box line 144) | implied by `dv` having no `dP` term; not separately asserted, but visible in output (`dKq/2 - 3dKs/4 + 13dTheta/8`) | partial |
| v/T, v*T sum/difference forms (notes §2 boxes) | wl asserts `(dv-dT)`, `(dv+dT)` against alt forms (wl 55–56); py prints (py 58–59) | match (wl) |
| primitives ρ_w, c_{s,w}, ell, a, L_W, c_s, Z_q (notes §1) | printed in both engines (py 44–50, wl 39–45) | match (input-encoded) |

Dominant pattern: aligned. Every paper/notes deliverable is exercised by at least the Mathematica engine; the SymPy engine under-asserts (prints the carry-forward laws rather than asserting them) but still asserts the load-bearing invariants/channels/δ⊥, which is the card's stated Output. No identity mismatch anywhere.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 77 | `expect_zero(2*dlam-dKs-dKq)` | Output: δln r_c=0 | yes |
| A2 | sympy | 78 | `expect_zero(dlam-(dKs+dKq)/2)` | Output: δln 𝔯=0 | yes |
| A3 | sympy | 79 | `expect_zero(dgq+dKs/2-dgs-dKq/2)` | Output: δln 𝔤=0 | yes |
| A4 | sympy | 84 | `expect_zero(dgq+dKs-dgs-dlam)` | notes §5 channel 1 | yes |
| A5 | sympy | 85 | `expect_zero(dKs+dKq-2*dlam)` | notes §5 channel 2 | yes |
| A6 | sympy | 89 | `expect_zero(gstar*chan1+chan2/(4√(1+rstar²)))` | notes §5 δ⊥=0 | yes |
| A7 | math | 53 | `expectZero(dv-(-3dKs/4+dKq/2+13dTheta/8))` | notes §2 v_{w0} box | yes |
| A8 | math | 54 | `expectZero(dT-(-5dKs/4+dKq/2+15dTheta/8-2dP/5))` | notes §2 T_m box | yes |
| A9 | math | 55 | `expectZero((dv-dT)-(2dcs-da))` | notes §2 v/T form | yes |
| A10 | math | 56 | `expectZero((dv+dT)-(dZ+5drho-4da))` | notes §2 v*T form | yes |
| A11 | math | 68 | `expectZero(dgs-(-dKs/4+dKq/2+3dTheta/8-2dP/5))` | notes §3 g_s box | yes |
| A12 | math | 69 | `expectZero(dgq-(-3dKs/4+dKq+3dTheta/8-2dP/5))` | notes §3 g_q box | yes |
| A13 | math | 70 | `expectZero(dlam-(dKs+dKq)/2)` | notes §3 lambda box | yes |
| A14 | math | 76–78 | `expectZero(drc/dr/dg)` | Output: invariants | yes |
| A15 | math | 83–87 | `expectZero(chan1/chan2/deltaPerp)` | notes §5 | yes |

None of these are tautological. Each invariant is built from the primitive transport laws (`dlam = dv + dIsq`, with `dv`, `dIsq` assembled from `drho/da/dcs/dZ/dcsw/dell/dLW`), so a wrong upstream coefficient would make the residual nonzero. The carry-forward `expectZero`s in the .wl are genuine primitive-route-vs-boxed-result closures (not `x==x`).

## Findings

### F1 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl:1-87`
- (compare) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.py:17-89`

**What's wrong:**
The `.wl` mirrors the `.py` line-for-line in scaffolding, variable choreography, and assembly order rather than re-deriving the result through an independent route. Corresponding excerpts:

Primitive block — py 31–37:
```
dTheta, dKs, dKq, dP = sp.symbols(...)
drho = sp.Rational(1, 2) * dTheta
da = sp.Rational(1, 2) * dKs - sp.Rational(1, 4) * dTheta
dcs = sp.Rational(1, 2) * dKs - sp.Rational(1, 4) * dTheta + sp.Rational(1, 5) * dP
dZ = dKq - sp.Rational(2, 5) * dP
dcsw = 2 * drho ; dell = -dcsw ; dLW = da
```
wl 31–37 (identical order, identical RHS):
```
drho = dTheta/2; da = dKs/2 - dTheta/4;
dcs = dKs/2 - dTheta/4 + dP/5; dZ = dKq - 2*dP/5;
dcsw = 2*drho; dell = -dcsw; dLW = da;
```
Branch-law assembly — py 53–54 vs wl 47–48 use the identical decomposition
`(dZ-drho)/2 + 3*dcsw/2 ± dcs - (5|3)*da/2`, and the parent-coupling block (py 62–65 `dgq/dgs/dIsq/dlam`) is reproduced verbatim in wl 58–61 including the non-obvious internal name `I_sq`/`dIsq`. The print-string blocks (py 91–98, wl 90–97) are character-for-character identical.

The shared primitive definitions are partly dictated by the physics (they ARE the Stage-165/166 carry-forward identities), so the overlap is not entirely avoidable. The mitigating factor is that the `.wl` carries strictly MORE independent verification than the `.py`: it `expectZero`s the carry-forward laws (A7–A13) and the v/T, v*T alternate forms (A9–A10) as primitive-route-vs-boxed-result closures, which the `.py` only prints. So the second engine does add real cross-checking value and is not a pure echo — hence low severity rather than a blocking transliteration.

**Why this matters:**
The second-engine policy wants the two engines to reach the result by independent algebra so that a shared transcription error cannot pass both. Here the assembly path is shared, so a sign/coefficient error in, e.g., the `dv` decomposition would be reproduced identically in both engines and would not be caught by engine cross-check. The extra `.wl` carry-forward assertions partly compensate by checking the assembled `dv`/`dT`/`dgs`/`dgq`/`dlam` against the independently-stated boxed laws from the notes.

**Required change:**
Have the `.wl` derive `dv`, `dT`, `dgq`, `dgs`, `dlam` from a route that does not replay the SymPy decomposition — e.g., start from the boxed carry-forward laws (notes §2–§3) as the primary definitions and then `expectZero` them against the primitive-assembled forms (the reverse of the current direction), or substitute concrete numeric drifts `(dTheta, dKs, dKq, dP)` and confirm each invariant/channel vanishes numerically as a second route. Either keeps a genuinely independent second derivation while leaving the SymPy script as-is. Cite `mathematica/...stage167...wl:47-61`.

**Verification:**
After Codex revises the `.wl`, the script still exits 0 with all PASS lines, but the assembly of `dv/dT/dgq/dgs/dlam` no longer mirrors `scripts/...stage167...py:53-65`; the new closures should be visible as additional `expectZero` lines or a numeric-substitution block in the refreshed `mathematica/output/...stage167...txt`.

## Independent-derivation check (Mathematica)

See F1. The `.wl` is a structural port of the `.py` (identical variable set, identical RHS in identical order, identical internal `dIsq` name, identical print strings), but it is not a pure echo: it adds A7–A13 and A9–A10 as independent closure assertions that the `.py` omits. Confidence the `.wl` is substantially a transliteration: high for the scaffolding/assembly, but moderate overall because the extra assertions give it genuine standalone verification content. Net call: low-severity `mathematica_transliteration`.

## Engine cross-check

Both engines agree exactly. Compare the primitive block (sympy output lines 5–11 vs mathematica output lines 5–11): identical. The parent couplings (sympy lines 18–21: `dgq = dKq - 3dKs/4 - 2dP/5 + 3dTheta/8`, `dgs = dKq/2 - dKs/4 - 2dP/5 + 3dTheta/8`, `dIsq = 5dKs/4 - 13dTheta/8`, `dlam = (dKq+dKs)/2`) match the mathematica output (lines 30–33) verbatim. Both report `δln r_c = δln 𝔯 = δln 𝔤 = 0`, both channels = 0, `delta_perp = 0`. The branch laws agree: sympy `dv = dKq/2 - 3dKs/4 + 13dTheta/8` (line 13) equals mathematica `(4dKq - 6dKs + 13dTheta)/8` (line 16). No residual, sign, or factor disagreement. `engines_agree: true`.

## Verdict justification

The stage's stated Output — first-order tangency of bundle drift to the compensated parent family — is faithfully and non-tautologically verified by both engines: the three invariants, the two off-family channels, and `δ⊥` all reduce to zero from primitive-assembled quantities, and every detailed transport law in the notes matches the printed/asserted output exactly. I attacked the invariants for tautology (they are built from upstream primitives, not pinned to their answers — not tautological), checked the `delta_perp` combination (depends on `chan1`/`chan2` which are independently derived — a nonzero channel would surface), verified the `dv` coefficient by hand (`13/8 dTheta`, `dP` cancels — matches), and confirmed the v/T and v*T forms against the notes boxes. No `paper_misalignment`: the card is terse but every deliverable lives in the notes and reconciles. The single finding is the borderline `mathematica_transliteration` (low) — the `.wl` shares the `.py`'s assembly path but adds genuine extra closure assertions, so it is not a blocking echo. Verdict: findings (one low-severity script-side finding); not stop_cold.

## Self-test notes

Checked: (1) no `sp.diff`/`D[]` derivatives present — all checks are linear-combination residuals, so the "derivative of independent variable is zero" trap does not apply; (2) no integrals — parity trap N/A; (3) trivial-case pre-check — substituting `dTheta=dKs=dKq=dP=0` gives all-zero (trivially passes) and substituting a single nonzero drift (e.g. `dKs=1`, rest 0) gives `dlam=1/2`, `drc=1-1=0`, `dr=1/2-1/2=0`, confirming the invariants vanish for a genuine nonzero drift, not just at the origin; (4) the F1 required change names the correct `.wl` path under `mathematica/`; (5) paper round-trip — the suggested `.wl` revision changes only the derivation route, not any constant, so it introduces no new paper_misalignment.

## Value Reconciliation (pass-2 augmentation)

Every emitted deliverable value was located in the notes `.md` (the natural carrier; the `.tex` card is intentionally terse and the appendix carries only the invariant triple). No stale `168π²`/`100π²` constant and no Family-1 radius literal `√(4107−100π²)/(10π)` appears anywhere in this stage — it is purely symbolic first-order transport algebra with no `π`, no `4107`. The "Family-1" mentions in card/notes are textual branch names, not the radius constant.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| δln ρ_w = dΘ/2 | py44 / wl39, sympy.txt:5 | notes.md:419 (§1, §7) | MATCH |
| δln c_{s,w} = dΘ | py48 / wl43, sympy.txt:9 | notes.md:80 (boxed) | MATCH |
| δln ell = −dΘ | py49 / wl44, sympy.txt:10 | notes.md:91 (boxed) | MATCH |
| δln a = dKs/2 − dΘ/4 | py45 / wl40, sympy.txt:6 | notes.md:426 | MATCH |
| δln L_W = dKs/2 − dΘ/4 | py50 / wl45, sympy.txt:11 | notes.md:102–107 (boxed); appendix:337 | MATCH |
| δln c_s = dKs/2 − dΘ/4 + dP/5 | py46 / wl41, sympy.txt:7 | notes.md:429 | MATCH |
| δln Z_q = dKq − 2dP/5 | py47 / wl42, sympy.txt:8 | notes.md:432 | MATCH |
| δln v_{w0} = dKq/2 − 3dKs/4 + 13dΘ/8 | py53 / wl53, sympy.txt:13 | notes.md:137–139 (boxed) | MATCH |
| δln T_m = −5dKs/4 + dKq/2 + 15dΘ/8 − 2dP/5 | py54 / wl54, sympy.txt:14 | notes.md:152–157 (boxed) | MATCH |
| δln(v/T) = dKs/2 − dΘ/4 + 2dP/5 | py58 / wl55, sympy.txt:15 | notes.md:165–172 (boxed) | MATCH |
| δln(v·T) = −2dKs + dKq + 7dΘ/2 − 2dP/5 | py59 / wl56, sympy.txt:16 | notes.md:175–181 (boxed) | MATCH |
| δln g_s = −dKs/4 + dKq/2 + 3dΘ/8 − 2dP/5 | py63 / wl68, sympy.txt:19 | notes.md:222–227 (boxed) | MATCH |
| δln g_q = −3dKs/4 + dKq + 3dΘ/8 − 2dP/5 | py62 / wl69, sympy.txt:18 | notes.md:232–237 (boxed) | MATCH |
| δln λ = (dKs+dKq)/2 | py65 / wl70, sympy.txt:21 | notes.md:243–247 (boxed); appendix | MATCH |
| δln r_c = 0 | py77 / wl76, sympy.txt:26 | notes.md:276; card:15; appendix:353 | MATCH |
| δln 𝔯 = 0 | py78 / wl77, sympy.txt:27 | notes.md:299; appendix:355 | MATCH |
| δln 𝔤 = 0 | py79 / wl78, sympy.txt:28 | notes.md:322; appendix:357 | MATCH |
| δln(g_q K_s/(g_s λ)) = 0 | py84 / wl83, sympy.txt:33 | notes.md:345–349 (boxed) | MATCH |
| δln(K_s K_q/λ²) = 0 | py85 / wl84, sympy.txt:34 | notes.md:345–349 (boxed) | MATCH |
| δ⊥ = 0 | py89 / wl87, sympy.txt:35 | notes.md:364 (boxed) | MATCH |

INTERNAL (scaffolding / intermediate, not stated deliverables, no finding): `δln I_sq = 5dKs/4 − 13dΘ/8` (parent-action `J_s ∝ a²ℓL_W^{1/2}` intermediate used to build λ; printed but not a boxed result), and the symbolic placeholders `gstar`, `rstar` (Stage-163 evaluation-point constants carried symbolically in δ⊥).

reconciliation: complete; 20 deliverable values checked, 0 misaligned
