I'll review this bounded native pressure-unit build from the fixed packet only, starting with the build spec, worker, library, and evidence guide.The review is still open. I am continuing through the fixed packet only, starting with the native C2 bridge and the unread face operands.The remaining chemical-unit records are the part that can still change the verdict. I am reading those stored units and the unread face operands from the packet only.The chemical expression uses several gamma families, and their stored units are not interchangeable. I am reading each used gamma’s four registry records before judging the expected chemical unit.The used gamma records agree. I am checking the face source, inverse, velocity, and density operands that the certificates actually select.CLEAR FOR THIS BOUNDED NATIVE PRESSURE-UNIT BUILD

This verdict covers only the unbound native pressure-source/closure unit bridge in the fixed packet. It does not accept a runtime result, and it does not clear source-grade/field transport, H/J/direct-profile or measure units, the all-544 address proof, or a numerical evaluator. No READY gate is present in this packet.

## Substantive blockers

None. The certificates attach to the original unbound constructors. The stored units of the symbols those constructors actually use are homogeneous at the units the worker requires. That conclusion is a structural reading of the saved constructors against `native-c2.py` line 59 and the gamma records below. The interpreter itself was not executed.

## Source locations

Selection and certificates are in `worker.py`. `join_native_case` (`worker.py` 116–133) parses the original restore text, requires one label match, and records the label census, index, full payload, and VALUE before the equality guard. `certify` (`worker.py` 102–113) walks that constructor, saves a refusal decision and the partial walk before the exception leaves, and returns a unit only after the expected-unit check. Controls (`worker.py` 249–257) emit the walk, then the decision, then require a message beginning `inhomogeneous native Add`.

| Join | Original constructor | Saved operand the certificate uses |
|---|---|---|
| Chemical `LAB_HELD` / `RHO4_CONSTANT`, index 0 | `native-b-exports.py` key `mu_theta_operator` at line 5628; receipt in `saved/native/source.json` 31–48 | Second leaf of `Tuple(Str('mu_theta'), Add(...))`, `worker.py` 172–178. The embedded dimension matrix stays inside the payload identity check. |
| Live density `RHO4_CONSTANT`, index 0 | `background_density_map` receipt, `source.json` 406–436 | VALUE leaf 1, `Mul(rho_br, Add(Mul(eta_bg, w1_profile), Integer(1)))`, equal to `binding-context.json` `densityMap.rho_br_bg_rho4_constant` |
| Velocity, both faces | `face_velocity` cases index 0 and 1, `source.json` 157–188 | `Mul(Rational(1, 2), W_0, e_W_t, epsilon_shape)`, equal to both `plus-source-input.json` and `minus-source-input.json` `nativeVelocity` |
| Raw c1 source | `DELTA_P` Add inside each face VALUE, `source.json` 87 and 112 | `epsilon_shape**(-1)` times that same Add, in both source-input files |
| Flat diagonal | `_LEDGER` begins at `native-c1-exports.py` line 76; key `dtn_kernel` is line 85. Display line 86 has four faces. `Inputs.self.kernel = cases(values['dtn_kernel'])` is `native-c2.py` line 199 | `saved/native/dtn-operands.json` literal: each `FLAT_DIAGONAL` is `omega*rho_m/q_out` times three top-level `DiracDelta` factors. Dropping those three factors is the flat-join Mul, residual `Integer(0)` |
| Dynamic Z | `kernel_bridge` lines 374–377 name `FLAT_DIAGONAL`, collect its deltas, replace them, and assign `DIMENSION_SCHEMA[z.name]` | `worker.py` 224–237 assigns the returned flat unit, after the `(-3,-1,1)` check, to `s11cc1_dtn_operator_lab_held_{plus,minus}`. Static schema line 59 keeps both symbols at `[0,0,0]` |
| Inverse and pressure | `RESOLVENT_DEFINITION` leaf 1 and `DELTA_P` on both faces, `source.json` 87 and 112 | Feedback coefficient without Z; full inverse Add on the face registry; `DELTA_P` product on that same registry |

The chemical gammas used by the unbound Add agree across `LEFTNativeSource`, `RIGHTNativeSource`, `LEFTPairing`, and `RIGHTPairing` in `saved/units/merged.json`:

- `gamma_s11cb_w_bg_04` is `(-2,-2,1)` (line 1552). It multiplies `u_i`, which the schema gives `(1,0,0)`.
- `gamma_s11cb_w_bg_06`, `_07`, `_08`, `_12`, `_13`, and `_14` are `(0,-2,1)` (lines 1724, 1810, 1896, 2240, 2326, 2412).
- `gamma_s11cb_mu_r_bg_04` is `(0,0,0)` (line 262). It multiplies `mu_R*u_i/W_0`.
- `gamma_s11cb_mu_r_bg_06`, `_07`, `_08`, `_12`, `_13`, and `_14` are `(2,0,0)` (lines 434, 520, 606, 950, 1036, 1122).

With those stored triples, every term family in the chemical Add — the `B_rho_3`, `C*W_0`, `G_theta_u`, `kappa_theta`, `kappa_theta_W`, and gamma families above — has exponent sum `(-1,-2,1)`. `worker.py` 169–170 loads these records and sets `newInference` false. `native-c2.py` `infer_dimensions` is not called.

Both raw sources are `epsilon_shape**(-1)` times an Add whose `Lambda_V_0/rho_m` pole sits beside `Integer(1)`. With schema `rho_m = (-4,0,1)` that inner Add is dimensionless and the raw source is `(1,-1,0)`. Reassigning `rho_m` to `(-3,0,1)` makes that inner Add inhomogeneous. The inverse Add is `Lambda_A_0*z/rho_m**2` plus the identity. The flat return `(-3,-1,1)` makes the feedback coefficient `(3,1,-1)` and the inverse dimensionless; pressure is then source times a dimensionless resolvent times that Z, which is `(-2,-2,1)`. Putting the static `(0,0,0)` on that same Z symbol makes the inverse Add inhomogeneous. Both controls therefore meet a real inhomogeneous Add on the actual constructors.

Inherited consumer rows are restored only. `consumer-unit-joins.json` has four `THETA_BALANCE` slots, `delta_p_plus`, `delta_p_minus`, `d_w_delta_p_plus`, and `d_w_delta_p_minus`, each with `total` equal to `expected`. The four `E_W_BALANCE` rows stay outside that filter (`worker.py` 261–264). `pressureSummandUnitsComplete`, `numericalEvaluatorReady`, and `sourceGradeFieldUnitTransportComplete` are returned false (`worker.py` 265–266).

Authority scope in `execution-authority.json` matches the manifest scope: real frequency 3, rest bulk, `LAB_HELD` / `RHO4_CONSTANT`, edges `1/5` and `1/10`, one scientific execution, no automatic retry, no deadline. The launcher command tail after `-u` is the worker argument tail (`launcher.py` 36–41, `worker.py` 272–273). `subprocess.run` of that command has no timeout. The 30-second `select` and the 50-step arming poll (`launcher.py` 47–49 and 94–104) are startup handshakes that refuse to launch when the hook is not armed.

## Coverage limits

The full node walk, label census on the untouched export strings, and hash identity of the chemical value line (161136 bytes) and face-response value line (206655 bytes) remain runtime checks. `worker.py` 142–147 performs those hashes. The displayed `dtn_kernel` faces match the saved literal’s four labels and three-delta `FLAT_DIAGONAL`; byte equality of the restore argument is the runtime `literalMatches` check (`worker.py` 197–200).

On this literal the three deltas are top-level factors, so dropping them yields the flat-join Mul. `kernel_bridge` instead replaces every `DiracDelta` by `One`. Those `One` factors are dimensionless, so the unit of the dropped Mul is the unit this bridge applies. The two constructors were not executed.

`merged.json` records `staticEntries` 823. The four-way agreement check in `worker.py` 161–168 covers every generated symbol, including symbols absent from the chemical Add. Only the gammas named above were read here. All four origin tails record `sourceReexecuted` false: native-source registries have 259 entries, pairing registries 271. Opaque pickle bytes were not opened.

`binding-context.json` stores parameter `omega` as the string `"1"` and frequency as `Integer(3)`. The worker requires that pair (`worker.py` 140). The unit walks use the unbound symbols.

The containment function body, guard, and supervisor are hash-pinned and are not in this packet. `worker.py` 274–275 executes the single `containment` function from the pinned helper before `run`. Their resource behavior is a runtime obligation. Synthetic tooling tests are not native unit evidence.

## Runtime obligations

A future gate must carry `READY_FOR_ONE_PACKET_NATIVE_UNITS`, independent build clearance, and both literal verdicts equal to the sentence above (`worker.py` 80–94). This review does not create that gate. The run must save every constructor, registry record, node walk, control, and posthash under containment, with the native memory cap, zero swap, and no automatic retry. Exit status `NATIVE_PRESSURE_UNIT_BRIDGE_COMPLETE_PENDING_INSPECTION` is inspection evidence only. Grade, jet, field, profile, measure, and all-544 address transport remain required before any evaluator can be called ready.