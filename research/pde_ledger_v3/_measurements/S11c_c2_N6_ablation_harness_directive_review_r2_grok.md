I'll start by reading the review-legs skill and the directive, then check each knife against the engines without touching the working tree.The cited directive path is missing; I'll locate the actual file and the engines it wraps.The artifact lives under `research/pde_ledger_v3`. Next I’ll read the directive, engines, template, and clearance records.I'll verify each knife against the live engines, starting with the WL and covariance sources.Next I’ll pull the construction sites the knives patch and confirm the APIs they call actually exist.I’ll read `face_factory`, `flatten`, `build_increment`, and the S11c-a APIs the patches call.I’ll run uniqueness and parse checks on the quoted fragments and the patched Python/WL replacements.Reviewed `research/pde_ledger_v3/directives/S11c_c2_N6_ablation_harness_directive.md` against the four live engines. Sites, patches, FORM/COEFFICIENT class, one-sided cones, K_normal scoping, PRINT-not-PASS, and certified-object coverage all hold. No blocking defect.

Check script: `/tmp/n6_r2_review/knife_checks.py`. Literal stdout (all checks completed; a trailing leftover `PY` then raised `NameError`):

```
FaceSource frozen dataclass: True
unique fragment counts:
  1  WL materialNormalKnife=0
  1  WL actualJunkCoefficient=0
  1  WL sourceChannel False
  1  COV ACTUAL_JUNK=0
  1  COV image assignment
  1  COV mu_pred block
  1  DIAG flatten E_W
  1  DIAG slots
  1  REC face_factory MATERIAL
  1  REC ms.get
  1  REC source=closed_response
APIs: True True True True True True
PARSE_OK K_EW_pop
PARSE_OK K_rank_form
PARSE_OK K_normal_wrapper
PARSE_OK K_source_unpack
PARSE_OK K_slotdrop_coeff
flatten_pop_keys ['U0', 'U1', 'U2', 'THETA']
velocity_before_carrier True
ms_uses_m_v True
WL_junk_case actualJunkCase = {"MATERIAL_ADVECTED", "RHO4_CONSTANT"};
WL_activeJunk activeJunk = If[case === actualJunkCase, actualJunkCoefficient, 0];
leak PASS 0
leak FAIL 0
leak must vanish 0
leak non-trivial 0
leak verdict 0
build_increment_slots_in_body 0
```

---

### H1 WL — `K_carrier`
- **Site unique.** `materialNormalKnife = 0;` occurs once (wl:17). FORM `= 1` / ×2 `= 2` / IDENTITY `= 0` are all parseable assignments.
- **Construction.** `invT[[1, 4]] += knife jet["WBg", {2}]` (wl:289) mixes a different WBg slope into the mapped covector. Carrier path `materialGeometry[..., materialNormalKnife]` (wl:860) is live; source path is a separate builder with knife `0` (wl:861).
- **FORM vs ×2.** Enabling that off-block addend leaves the uncontaminated Jacobian family; `= 2` only scales it.
- **Cone.** `cm` (wl:868) and therefore `R_N6` / `SPLIT_CHECK` move. DEAD list is `CARRIER_EULERIAN`, `SOURCE_*`, `R_COV`, `FROZEN_PHI`. `R_COV_INCREMENT` (wl:896–900, uses `cm`) and `ACTUAL_CONTROL_PARAMETERS` (embeds `materialNormalKnife` at wl:980) are correctly excluded. Bites at `LAB_HELD × RHOBR_CONSTANT` (`jet["WBg", {2}]` is a live symbol; mapped x-component picks up `s * knife * jet`).

### H1 WL — `K_junk`
- **Site unique.** `actualJunkCoefficient = 0;` once (wl:20). FORM/×2 enable `jk junkMu eW` at wl:247.
- **Case gate is engine-owned.** `activeJunk = If[case === actualJunkCase, actualJunkCoefficient, 0]` (wl:845) with `actualJunkCase = {"MATERIAL_ADVECTED", "RHO4_CONSTANT"}` (wl:21). The directive keeps that gate. H1 runs all four cases (wl:1011–1014), so the junk case is in the transcript. Not a pinned-case no-op of the dropped-advection kind.
- **DEAD is right.** `faceLaws` puts `mu` in `affinity` (wl:317) but `carrier[rows_] := D[rows, p] /. slotZero` (wl:336) kills the pressure-independent junk addend, so `CARRIER_*` stay independent. `SOURCE_PREDICTED` / `FROZEN_PHI` / `PHI_DOMAIN_CENSUS` never see `activeJunk`.

### H1 WL — `K_split_route`
- **Site unique.** Full `sourceChannel = joinFaceMaps[buildContraction[cm[#], es[#] - ms[#], response[#], #, False] &];` occurs once (wl:892). Bare `es[#] - ms[#]` also appears on `crossChannel` (wl:893); the directive scopes the edit to this assignment.
- **FORM.** `affine=True` selects `kernelFamilies` whose first entry is `"LOCAL_BARE"` (wl:485, 576) and writes `entry = gScale[cg, -pressureSlots[[slot]]]` (wl:571–572) — i.e. `-C_M·p` into the source channel. That extra bare family is exactly what `allFamilies` currently zero-fills (wl:718–721, 894), so `SPLIT_CHECK` (wl:895) leaves the identity. ×2 keeps `False` and only rescales `es-ms`.
- **DEAD.** `CARRIER_*`, `SOURCE_EULERIAN/MATERIAL`, `R_COV` are all built at wl:874–886, before the patched assignment.

### H2 covariance — `K_junk`
- **Site unique.** `ACTUAL_JUNK = sp.Integer(0)` once (cov:35). FORM/×2 feed `actual.append(kappa_j * junk * b.e_W)` (cov:145), unconditional at the pinned `LAB_HELD × RHOBR_CONSTANT`.
- **DEAD.** `SOURCE_PREDICTED` is `source_terms(..., mu_pred, ev)` (cov:191); `FROZEN_PHI` / `PHI_DOMAIN_CENSUS` are emitted inside `prolonged_phi` (cov:105–124) before `actual_amplitudes`.

### H2 covariance — `K_circular`
- **Site unique.** The three-line block at cov:186–188 matches the quoted old fragment exactly (count 1). FORM reorders to `mu_pred = mu_actual` and leaves `predicted_amplitude` unpatched; line 191 still builds `SOURCE_PREDICTED`. That is an operand identification, not a coefficient. ×2 is `[2 * term for term in predicted_amplitude(...)]`; `predicted_amplitude` returns a one-element list (cov:134).
- **DEAD.** `SOURCE_ACTUAL` still comes from `mu_actual` + `mv` (cov:192). `FROZEN_PHI` already emitted.

### H2 covariance — `K_rank`
- **Site unique.** The `images[atom] = b.total_derivative(...)` assignment inside nested `image_of` is unique (cov:92–96). `paths[atom][1]` is the multi-index (cov:78, 85). `b.total_derivative(..., background_depth=)` exists and is keyword-only (S11c-b:782–786).
- **FORM.** Rank ≥ 2 is identified with the rank-1 parent image — a jet-order collapse, not a rescale. Coverage is unchanged: every `atom in paths` still gets an `images[atom]` entry, so `if uncovered: raise` (cov:125–126) does not newly fire. `background_depth` / BFS truncation are not used.
- **DEAD.** `SOURCE_ACTUAL` does not consume `phi`. Certified set includes `FROZEN_PHI` and `PHI_DOMAIN_CENSUS`. Rank-2 jets are present in imported `μ_E` (clearance structure: `theta_didi` / `e_W_didi`, max rank 2), so the knife is live at the pinned case.

### H3 diagnostic — `K_EW_rowdrop`
- **Site unique.** Flatten line at diag:311, bound to `m_rows` at diag:806.
- **`pop` is a real row-family removal.** `flatten` zips `ROWS = ('U0','U1','U2','THETA','E_W')` (diag:44, 95–96). After flatten the key `'E_W'` exists; `rows.pop('E_W')` leaves `['U0','U1','U2','THETA']` (script stdout). Not a scalar.
- **No crash against `build_increment`.** Output domain is `BLOCKS × GRADES × signatures × faces` filled by `setdefault` (diag:510–512), not row keys. `slots` is unused in the body (count 0). `ROW_DIM[row]` is only hit for remaining coeff keys (`U0`…`THETA` are in `ROW_DIM`, diag:47–48). `unflatten` (diag:99–101) is only used inside `template` on a fresh `dict.fromkeys(ROWS)` (diag:434–435), not on `m_rows`. Eulerian imported rows (diag:796) are untouched.
- **Path to the certified residual.** `m_rows` → `m_coeff` (809) → `M` (851) → `REP_INVARIANCE_RESIDUAL` (853). Imported `E` is DEAD. Sign-flip is correctly labeled COEFFICIENT and does not replace the ×2 companion.

### H3 diagnostic — `K_slotdrop`
- **Site unique.** `slots=tuple(...)` once in `run` (diag:788). FORM remaining names are original order minus `delta_p_plus`. Same tuple feeds coeff extraction (797, 809) and both increment calls (850–851).
- **×2 companion.** `globals()['pressure_coefficients']` is the module function (diag:410–412); the nested wrapper parses; `key[1]` is the slot symbol. All four slots kept; only the extracted `delta_p_plus` column is scaled.
- **DEAD.** `MU_RECONSTRUCTION_*` come from `constitutive` (diag:819–824); `N6_MU_AMPLITUDE` / `N6_FACE_VELOCITY` from `mu` / factory velocity (825–828), not from `slots`. No advection mutation. At `RHOBR_CONSTANT` the engine’s own `N4_ADVECTION` triple is emitted (diag:860–870).

### H4 reconcile — `K_normal`
- **Site unique.** `n.face_factory(..., 'MATERIAL', mu_slot)` once, inside `build_material_carrier` (rec:91).
- **Patch APIs exist.** `a.build_material_face_source` (S11c-a:703), `a._FACE_CACHE` (a:626), `a.grad_W` (a:129), `a.dot` (a:390–391), `n.replace` = `dataclasses.replace` (diag:24) on `@dataclass(frozen=True) class FaceSource` (a:604–605). Wrapper parses. `build_face_source` looks up `build_material_face_source` in the S11c-a global namespace at call time (a:832–835), so the monkeypatch is seen.
- **Carrier-only.** `m_v = build_material_velocity(...)` (rec:204) runs *before* `build_material_carrier` (rec:205). Velocity uses `face_velocity_raw(source.normal_exact)` (a:890–891) on an unpatched builder (rec:84–85). `ms` is `source_terms(..., mu_m, m_v[s])` (rec:208). Wrapper is installed only around the carrier `face_factory` call and restored in `finally`, including cache. Does not patch Eulerian `normal_exact` (a:850) or `material_inverse_transpose` (a:694). Does not globally patch S11c-a.
- **FORM vs ×2.** `components[0] += a.grad_W[1]` mixes a different slope channel into the already-normalized material 4-normal; ×2 only scales existing component 0. Live at `LAB_HELD` (`grad_W[1]` is a live jet; component 3 of the unit normal stays nonzero, so renormalize does not divide by 0).
- **DEAD.** `SOURCE_*` and `EULERIAN_OPERAND` do not read the wrapped `normal_exact`.

### H4 reconcile — `K_source_route`
- **Site unique.** `source = closed_response(comp, inputs, m_coeff, ds, kernels)` once (rec:260). FORM unpacks `n.build_increment(...)` which returns `(output, dimensions)` (diag:513) and takes `slots` (in scope at rec:192–193). That injects signature-0 `-C·p` (diag:491–493) into a channel that `closed_response` restricts to 6/9/12 (rec:115, 124). Line 259 (`carrier`) is not edited. DEAD `CARRIER_*` are already built at rec:250–254.

### H4 reconcile — `K_operand_swap`
- **Site unique.** `ms[s].get(w, sp.S.Zero)` once (rec:211), inside `ds = es - ms`. FORM replaces only that read with `es[s].get(...)` → `es − es`. Not the whole-operand `ms↔es` exchange (that is `ds → −ds`, coefficient ×−1, correctly excluded). ×2 is `2 * ms[s].get(...)`, a rescale of the existing `ms` addend. `ds` feeds `source`/`cross` (rec:260–261), not `CARRIER_*`. At pinned `RHOBR_CONSTANT`, `es ≠ ms` (engine’s N4 advection is live), so `SOURCE_CHANNEL` moves.

### PRINT-not-PASS, bounds, completeness
- Zero hits for `PASS` / `FAIL` / `must vanish` / `non-trivial` / `verdict` in the builder-facing directive. Contract is tagged `{baseline, corrupted, diff}` plus named DEAD prints, no value claim.
- Builder sequence is build → run once → report paths → stop; no other model; no commit.
- Certified sets cover: WL `R_N6` / `SPLIT_CHECK` / `R_COV` under the engine’s own primes; covariance `R_COV` plus non-circularity (`K_circular`) and rank-2 `Φ` (`FROZEN_PHI`, `PHI_DOMAIN_CENSUS`); diagnostic `REP_INVARIANCE_*` (the engine’s `residual(E,M)` at diag:850–853) plus both `TILT` and `N4_ADVECTION` control triples; reconcile `R_N6` / `SPLIT_CHECK` / both channels / both bridges.

No finding survives the physics filter.

**SOUND**
