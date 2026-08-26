# Independent physics review — S11c-a WL background-current fix (Grok)

Artifact: `research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` vs `a7459cb8` (~lines 435–448).
Spec authority: `directives/S11c_a_SHARED_PHYSICS.md`. Directive not reviewed.

## 1 · Independent derivation (from the spec)

**§1b** (mass current):

```text
j = ρ_4D v_bulk ,           ∂_tρ_4D + ∇₄·j = 0 .
```

**§1** (drain is scope limit only): `v_bulk_normal_0` “remains only the inherited rest-frame scope limit.”

**§2d** (background state): `𝔅⁰` carries `ρ_4D,bg⁰` (nonzero density map) and supplies `V_s⁰ = 0`, `J_s⁰ = 0`.

**§3c** (shifted traces; free-premise ban):

> Every background face value or normal derivative appearing in this law is obtained by differentiating a
> member of the supplied background state `𝔅⁰` (§2d); **none may be introduced as a free premise.** In this scope
> the traced bulk velocity, the perturbation pressure, and the bulk current have zero background — `V_s⁰=J_s⁰=0`
> (§2d), the drain `v_bulk_normal_0` is the inert rest-frame scope limit of §1, `δp` has background value zero
> (§3b), and **the rest-frame background current `ρ_4D⁰v_bulk⁰` vanishes** — and the supplied density background
> depends on the in-plane anchor, not on `w`.

Given the engine’s rest-frame `bulkVelocityZero ≡ {0,0,0,0}` and §1b’s product law,

```text
j⁰ = ρ_4D⁰ · v_bulk⁰ = ρ_4D⁰ · 0 = 0 ,
```

while `ρ_4D⁰` remains the nonzero symbolic background density. Vanishing of `j⁰` is a **consequence** of `v_bulk⁰=0` under §1b, not a free `current*Background` premise (banned by §3c).

CAS script + stdout:

- `/tmp/s11ca_wl_bgcurrent_review_grok/derive_j0_from_spec.py`
- `/tmp/s11ca_wl_bgcurrent_review_grok/derive_j0_from_spec.stdout`
- copies under this `_measurements/` directory

Literal stdout:

```text
OPERAND_rho0 = rho_4D^0
OPERAND_v_bulk0 = [0, 0, 0, 0]
PRODUCT_j0 = rho0 * v_bulk0 = [0, 0, 0, 0]
RESIDUAL_j0_minus_zero = [0, 0, 0, 0]
j0_is_identically_zero = True
rho0_stays_symbolic_nonzero = True
PROBE_FORM_j = [rho_4D^0*vp1, rho_4D^0*vp2, rho_4D^0*vp3, rho_4D^0*vpW]
PROBE_currentW = rho_4D^0*vpW
PROBE_currentX = [rho_4D^0*vp1, rho_4D^0*vp2, rho_4D^0*vp3]
```

## 2 · Mandatory form ablation — background-velocity probe

Isolated copy of the engine’s `bulkCurrentZero` / `currentWZero` / `currentXZero` forms; `bulkVelocityZero` replaced by `{vp1,vp2,vp3,vpW}`. Full audit engine **not** run (`timeout 600`; seconds-scale `wolframscript -file`).

Scripts/stdout:

- `ablate_velocity_probe.wl` / `.stdout`
- `ablate_real_zero_velocity.wl` / `.stdout`
- `ablate_pkspec1_check.wl` / `.stdout`

### Nonzero probe (faithfulness)

Literal probe output (abridged):

```text
PROBE_bulkCurrentZero = {vp1*rhoBulkBackground[...], ..., vpW*rhoBulkBackground[...]}
PROBE_currentWZero = vpW*rhoBulkBackground[...]
PROBE_currentXZero[1] = vp1*rhoBulkBackground[...]
PROBE_currentXZero[2] = vp2*rhoBulkBackground[...]
PROBE_currentXZero[3] = vp3*rhoBulkBackground[...]
TRACKS_probe_W = True
TRACKS_probe_X1 = True
TRACKS_probe_X2 = True
TRACKS_probe_X3 = True
STAYS_IDENTICALLY_ZERO_UNDER_PROBE = False
```

**Result:** background current **tracks** the probe (`j⁰ → ρ⁰·v_probe`). Not a disguised `:= 0`.

### Real zero velocity (reverse)

```text
REAL_bulkVelocityZero = {0, 0, 0, 0}
REAL_bulkCurrentZero = {0, 0, 0, 0}
REAL_currentWZero = 0
REAL_currentXZero[1..3] = 0
SYMBOLIC_INDEX_currentXZero[k] = currentXZero[k][{x1, x2, x3}, w]   (inert)
Part_pkspec1_message_count = 0
rhoBulkZero_still_nonzero_head = True
```

**Result:** with real `bulkVelocityZero`, background current is `0`; density head `rhoBulkBackground` preserved; no `Part::pkspec1`; symbolic index stays inert under `index_Integer`.

## 3 · No-regression / scope

Source diff vs `a7459cb8`: **one hunk**, +5/−4, only the background-current block.

| Object | vs baseline |
|---|---|
| `rhoBulkZero` / `rhoBulkWave` | identical |
| `currentWWave` / `currentXWave` | identical |
| `bulkVelocityZero` / `pressureZero` / `heightSource` | identical |
| `currentWBackground` / `currentXBackground*` in `.wl` | **0** matches |
| `current*Background` in `/tmp/s11ca_wl_bgcurrent_full.out` | **0** (baseline out had 307+547) |
| `rhoBulkBackground` counts in outs | 63 = 63 |
| Unique `S11CA_*` tag set | equal (40 `WL_S11CA_*` payloads; 41 name strings) |

Pure face-geometry payloads identical BASE↔NEW: `FACE_MAP_*`, `FACE_NORMAL`, `FACE_VELOCITY`, `FACE_MEASURE_SHAPE_DERIV`, `CONORMAL_DERIV`, `BACKGROUND_STATE`, `BACKGROUND_DENSITY_MAP`.

Tags that **consumed** free `current*Background` in baseline (`FACE_SHIFT`, projection dynamic/shape/residual/term-origins, control-form, uniform-limit operands) shrink and lose those symbols — expected under the §3c product law, not a geometry regression. Perturbation-current **definitions** unchanged; relative-flux / kinematic-balance / traction / evolution payloads byte-identical.

Commands supporting the greps: `regression_greps.txt`, `tag_payload_compare.txt`, `current_block.diff` in this directory.

## 4 · Physics findings

**None.** The prior free `current*Background` symbols were the §3c free-premise defect; replacing them by `ρ⁰·v⁰` with `v⁰=0` is the required construction. Density background stays nonzero. Perturbation currents untouched. No disguised zeroing (probe tracks).

## 5 · Verdict

**Faithful to §1b / §2d / §3c.** The change is the correct rest-frame background-current construction: `j⁰ = ρ_4D⁰ v_bulk⁰` vanishes because `v_bulk⁰ = 0`, while `ρ_4D⁰` remains.
