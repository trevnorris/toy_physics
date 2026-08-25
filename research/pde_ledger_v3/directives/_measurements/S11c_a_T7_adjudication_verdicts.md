# Measurements twin — S11c_a_T7_adjudication_verdicts.md

Inputs (hash-locked): PY `~/.s11_build/S11c_a_sympy_engine.out` sha
`6386471555b1e99d0aeb0f716eea30f839d59be50c0cedd4677ea7b376b79129`; WL committed `.out` sha
`82062bd36cfb07b1f18631077f0c63ac1cbce7834967686f680fa9f30019e4ec`. Spec `2926c71c`.
⚠ rule-5: no literal output-expression fragment appears below (both legs flagged the earlier draft's
fragments); divergences are described structurally, ⛔ never by their computed value.

CONVENTIONS (all commands below): run from repo root `/var/projects/toy_physics` with
`PY=/home/trevnorris/.s11_build/S11c_a_sympy_engine.out`,
`WL=research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`,
`SPEC=research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`.

## VERDICT B — density
Script + literal stdout: `_measurements/s11ca_t7_adjudication/density_adjudication.{py,stdout}` (the script
reads both streams and pairs cases by (branch,face,dof), diffing the two rep payloads; prints DIFFER/
IDENTICAL counts, no expression values).
Result: WL KINEMATIC_BALANCE `DIFFER=0 IDENTICAL=8`; WL RELATIVE_FLUX `DIFFER=0 IDENTICAL=8`; PY
PROJECTION_SHAPE_DERIV / _DYNAMIC_OPERAND / _RESIDUAL `DIFFER=4` each; PROJECTION_STATIC_OPERAND `IDENTICAL=4`.
PROJECTION_TERM_ORIGINS (fold) CMD: same pairing over `^PY_S11CA_PROJECTION_TERM_ORIGINS:` ⇒ `DIFFER=4`.
The projection DIFFER is a STRUCTURAL difference in the density-carrying term (the RHO4 and RHOBR integrands
are built from different density-profile factorizations, §2b), ⛔ not a Dummy-index artifact — confirmed by
inspecting that the differing region is inside the integrand, not a leading label (value not reproduced here).
WL KINEMATIC value is velocity-based (its variation is carried by the bulk-velocity trace and the
geometric half-width, ⛔ not a rep-dependent density factor — which is why the two rep payloads coincide);
its `RESIDUAL_A_MINUS_B` field is the definitional identity.
Spec: `sed -n '351p;352p' "$SPEC"` ⇒ `J_s^α ≡ ρ_m (v_bulk − v_face)·n̂`; `n̂·v = V + J/ρ_m`. `sed -n '76p' "$SPEC"`
(§1b projection uses ρ_4D). `sed -n '315p' "$SPEC"` (§3a ρ_4D anchoring-dep). `sed -n '398p' "$SPEC"`
("both reps wherever density enters"). ⚠ `sed -n '229p' "$SPEC"` names only T-g/T-h/T-i (NOT T-f). `ρ_m` is a
bound S11b KNOB (`LEDGER['rho_m']`), rep-independent.

## VERDICT C — virtual work
Script + stdout: `_measurements/s11ca_t7_adjudication/virtualwork_offdiagonal.{py,stdout}` (WL-only test —
it reads the WL stream and reports, per (physical,virtual) DOF, whether the `SHAPE_DERIVATIVE.EXPRESSION`
field is the literal `0`; it prints booleans, ⛔ no expression values). Result: all 8 off-diagonal
(phys≠virt) are `zero=False`; all 8 diagonal `zero=False`. WL total keys = 16; PY = 8 (diagonal).
Physical-DOF redundancy (fold, both legs): for fixed (branch,rep,virtual DOF) the two WL physical-DOF rows
are byte-identical across the FULL emitted record (Codex `s11ca_t7_independent_check`: "WL complete emitted
records distinguish physical rows: False"). Producer: PY `:919-924`; WL `:1051`. Spec: `sed -n '417,419p' "$SPEC"`
(T-d "…which virtual-displacement pairings occur is part of the computation"). ⚠ PY off-diagonals are
UNMEASURED (PY emits diagonal only) — the PY patch supplies the cross-engine check.

## VERDICT H.1 — coverage (executable command in this twin, per Codex rule-2 note)
CMD (complete quantity vocab, both streams):
```
QUANTS='FACE_NORMAL|CONORMAL_DERIV|FACE_MEASURE_SHAPE_DERIV|FACE_VELOCITY|RELATIVE_FLUX|KINEMATIC_BALANCE|TRACTION|VIRTUAL_WORK_SHAPE_DERIV|CLOSURE_SHAPE_DERIV|VIRTUAL_CONSTRAINT|EVOLUTION_MASS_BALANCE|EVOLUTION_TERM_ORIGINS|PROJECTION_SHAPE_DERIV|PROJECTION_STATIC_OPERAND|PROJECTION_DYNAMIC_OPERAND|PROJECTION_RESIDUAL|PROJECTION_TERM_ORIGINS|FACE_SHIFT'
comm -23 <(grep -m1 '^PY_S11CA_CONTROL_FORM_BASE_OPERAND:' "$PY" | grep -oE "Str\('($QUANTS)'\)" | sed "s/Str('//;s/')//" | sort -u) \
         <(grep -m1 '^WL_S11CA_CONTROL_FORM_BASE_OPERAND:' "$WL" | grep -oE "\"($QUANTS)\|" | sed 's/[|"]//g' | sort -u)
```
⇒ PY-only = {EVOLUTION_TERM_ORIGINS, PROJECTION_STATIC_OPERAND, PROJECTION_DYNAMIC_OPERAND,
PROJECTION_RESIDUAL, PROJECTION_TERM_ORIGINS}; the reverse `comm -13` ⇒ ∅. Same 5 absent from WL
UNIFORM_LIMIT (grep the WL uniform histogram). Spec: `sed -n '424,429p;435,438p' "$SPEC"` (§4 T-f/T-h enumerate the
five separately); `sed -n '490p' "$SPEC"` (§5b every T-object); `sed -n '497p' "$SPEC"` (§5c each S11c-a object).

## VERDICT BG — background state
CMD (WL loads): `grep -m1 '^WL_S11CA_BACKGROUND_STATE:' "$WL" | grep -oE '"[A-Z_]+"' | sort -u` ⇒ incl.
BOUNDARY_LOADS. CMD (PY zeros present — Codex rule-2 note):
`grep -m1 '^PY_S11CA_BACKGROUND_STATE:' "$PY" | grep -oE "Equality\(Symbol\('(theta_0|V_0_[A-Z_]+|J_0_[A-Z_]+|A_0_[A-Z_]+)'\), Integer\(0\)\)"`
⇒ theta_0, V_0_*, J_0_*, A_0_* all `=0` present. CMD (PY loads absent):
`grep -m1 '^PY_S11CA_BACKGROUND_STATE:' "$PY" | grep -ocE "load|BOUNDARY|t_hold|f_hold"` ⇒ 0. PY holds live
only in `ADMISSIBILITY_PREMISE` (`f_hold_0, t_hold_±_0`). Spec: `sed -n '246,272p' "$SPEC"` — §2d:251 `𝔅⁰ ≡ {…,
boundary loads}`, :254 `θ⁰≡0`, :261 zeros, :267-272 "Emit the state … S11CA_BACKGROUND_STATE".

## BACKGROUND_DENSITY_MAP adjudication (legs split → §2b governs)
CMD: pair PY density-map cases by rep, compare LAB_HELD vs MATERIAL_ADVECTED payload ⇒ identical for both
reps (branch-independent). Spec: `sed -n '228,231p' "$SPEC"` — "construct Σ_E⁰(y) … Emit the two computed maps"
(two = per REPRESENTATIVE, on pre-anchoring y); `sed -n '242,244p' "$SPEC"` (§2c both branches — applies to
anchoring-dependent profiles). ⇒ 2-per-rep spec-faithful (WL); PY branch axis redundant (shallow).

## Discipline
Each verdict = computation + spec citation (rule 4/rule 2). No expected OUTPUT value is stated (rule 5);
the patch directives specify WHAT to compute only. Two legs verified every engine assignment before any
builder (rule 7/rule 13); corrections folded above.
