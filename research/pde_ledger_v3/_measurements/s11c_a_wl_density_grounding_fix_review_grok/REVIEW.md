# Independent build review — WL density-grounding fix (Grok)

Artifact: `research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl`
(+ regenerated `.out`). Diff vs `HEAD`.

## Verdict

**Mostly §3c-correct and complete across the three shifted-trace sites**, with **one finding**.

The free-premise `rhoBulkBackground(x,w,t)` is gone from all three FACE_SHIFT-bearing packages;
background face values come from the §2b representatives; `∂_w ρ⁰ = 0` *emerges* from
differentiation (not asserted); the density axis is present at all three sites; the form control
bites; pressure exact/shape match the pre-fix baseline; collateral is confined to the expected
shifted-trace families.

**Finding (incomplete §2c branching on the density background):**
`traceSource` always evaluates `fieldZero` at `spatialCoordinates`
(`S11c_a_interface_geometry_mathematica_audit.wl:428`), while
`rhoBulkRepresentativeZero` ignores the branch
(`:446–447`). For `MATERIAL_ADVECTED|RHOBR_CONSTANT` the EXACT_TRACE background is
`rhoBr/widthProfile[x1,x2,x3]`, not `ρ_4D,bg⁰(χ)` as required by §2c/§3a and as already wired
elsewhere via `branchedRho4Source` (`:666–667`). Script
`probe_material_branch_density.wl` shows the missing shape contribution is the material transport
of `W_bg(χ)` — nonzero residual vs the branched path. LAB_HELD and RHO4_CONSTANT are unaffected
(constant / χ→x at background order for the O(1) face value).

## 1. Independent §3c derivation

Script:
`research/pde_ledger_v3/_measurements/s11c_a_wl_density_grounding_fix_review_grok/derive_s3c_density_trace.py`

Literal stdout (`derive_s3c_density_trace.stdout`):

```
=== background normal derivative by differentiation (must emerge) ===
RHO4_CONSTANT: D[rho4_bg, w] = 0
RHO4_CONSTANT: D[rho4_bg, w] == 0 ? True
RHOBR_CONSTANT: D[rho4_bg, w] = 0
RHOBR_CONSTANT: D[rho4_bg, w] == 0 ? True

=== §3c shift term δh · ∂_w ρ⁰ |_{h0}  (DOF=DELTA_W, ζ_c=0) ===
RHO4_CONSTANT: delta_h * D_w(rho0)|_h0 = 0
RHOBR_CONSTANT: delta_h * D_w(rho0)|_h0 = 0

=== composed first-order δ[ρ] at LAB_HELD/MINUS/DELTA_W ===
RHO4_CONSTANT: first-order delta[rho] = -W_0*delta_rho_4D_face_minus_dw*eta_bg*w1_profile(x1, x2, x3)/2 + delta_rho_4D_face_minus
RHOBR_CONSTANT: first-order delta[rho] = -W_0*delta_rho_4D_face_minus_dw*eta_bg*w1_profile(x1, x2, x3)/2 + delta_rho_4D_face_minus
```

Engine match (LAB_HELD|RHO4|FACE_MINUS|DOF_DELTA_W|FIELD_BULK_DENSITY):
EXACT `rhoBr/W0 + waveOrder*rhoBulkPerturbation[..., height]`;
SHAPE `rhoBulkPerturbation[..., -W0/2] - (etaBg*W0*w1Jet[0,0,0]*Derivative[...][rhoBulkPerturbation])/2`
— free-premise `Derivative[rhoBulkBackground]` absent; perturbation retains first-shape-order η.

## 2. FORM ablation (mandatory)

Copy under `/tmp/.../ablate_density_bg_form.wl` (also saved in `_measurements/...`).
Replaced grounded `rhoBr/W0` with `rhoBr/W0 + alphaTest*w^2` (form change, not rescale).

Literal stdout:

```
=== BEFORE (grounded RHO4_CONSTANT; w-independent) ===
D[bg,w] = 0
SHAPE_DERIVATIVE = rhoBulkPerturbation[x1, x2, x3, {-1/2*(W0*(1 + etaBg*w1Jet[0, 0, 0])), time}]
=== AFTER (FORM ablation: bg -> rhoBr/W0 + alphaTest*w^2) ===
D[bg,w] = 2*alphaTest*normalCoordinate
AFTER - BEFORE = (alphaTest*W0*deltaWidth[x1, x2, x3, time] + alphaTest*etaBg*W0*deltaWidth[x1, x2, x3, time]*w1Jet[0, 0, 0])/2
shift_reappeared? True
```

The shift is wired through a real `D[·,w]`.

## 3. Independence / not a hardcode

`ablate_representative_independence.wl` stdout:

```
RHO4_CONSTANT face = rhoBr/W0
RHOBR_CONSTANT face = rhoBr/(W0*(1 + etaBg*w1Jet[0, 0, 0]))
faces_differ? True
```

Transcript: RHO4 EXACT carries `rhoBr/W0`; RHOBR EXACT carries `rhoBr/widthProfile[x1,x2,x3]`.

## 4. Collateral

Changed top-level tags vs `HEAD` (only these seven):

| Tag | Why expected |
|---|---|
| `WL_S11CA_FACE_SHIFT` | primary T-e + DENSITY axis |
| `WL_S11CA_CONTROL_FORM_{BASE,ABLATED,RESIDUAL}` | form-control FACE_SHIFT slice + DENSITY axis |
| `WL_S11CA_UNIFORM_LIMIT_{S11CA,S11B}_OPERAND`, `..._RESIDUAL` | uniform-limit FACE_SHIFT slice + DENSITY axis |

Within form-control and uniform-limit sections, **only** `FACE_SHIFT|*` keys were added/replaced;
no non-FACE_SHIFT stem moved. No other of the 40 tags changed.

## 5. All three sites

| Site | Loop / keys | Lines |
|---|---|---|
| Primary `FACE_SHIFT` | `{branch, densityNames, sign, traceFieldNames, dof}` via `densityFaceCaseKey` + `\|FIELD_` | `:1120–1136` |
| Form-control slice | same iterators; keys `FACE_SHIFT\|` <> `densityFaceCaseKey` <> `\|FIELD_` | `:1883–1898` |
| Uniform-limit slice | same; keys `FACE_SHIFT\|` <> `densityFaceCaseKey` <> `\|FIELD_` | `:2097–2109` |

Inventory: `traceFieldInventory[density]` with
`BULK_DENSITY -> {rhoBulkRepresentativeZero[density], ...}` (`:625–628`).

## 6. Blindness

No `Get`/`Needs`/`Import` of the SymPy sibling. Only `FileNameJoin` for form-control temp files
(`:1698`). Engine re-derives from local ansatz + §2b `rho4Profile` (`:305–310`).

## Finding detail — MATERIAL_ADVECTED RHOBR

`probe_material_branch_density.stdout`:

```
SHAPE FACE_SHIFT path: 0
SHAPE branched path:   (etaBg*rhoBr*uThree*... + ...)/(W0*(1 + etaBg*w1)^2)
material_transport_missing_from_FACE_SHIFT_path? True
```

Fix scope for a follow-up: evaluate the density zero through `branchedProfileCoordinates[branch]`
(same object other density-bearing families already use), not bare `spatialCoordinates`.
