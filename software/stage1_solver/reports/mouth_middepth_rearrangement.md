# Restoring force of a rearranged, flowing patch at mid-depth

**The compressed material formed a zigzag that repeated every two arriving sections. After rearrangement it pushed back against a weak sideways ripple in both directions.** With the finite-displacement wall contact law, a ripple of amplitude +0.0001 bare bond lengths gave a mean restoring force of −0.125488 in model units. Reversing the step gave +0.125488. [Measurements and independent local review](../../../research/pde_ledger_v3/_measurements/mouth_middepth_rearrangement.json)

**Preparation.** For a radius-one mouth at intake 0.25, the earlier flowing strip supplied the mean section spacing 0.67952, the hoop stretch 1.90500 and the arrival time 0.0541596 at depth 1.175. Local shear order was 0.87530. We held those conditions fixed in a cold patch, fed centered material through its entrance and let sites move sideways with unit drag. Nearby sites interacted through the existing Morse bond force and could change neighbors. Wall attachment ends stayed fixed. We tested both the linearized wall law and the finite-displacement Morse law averaged over the orientations of nine measured offsets. [Feed and compression](mouth_joint_strip_hold.md) · [Material inputs](../tools/mouth_middepth_rearrangement.py) · [Patch evolution and readout](../tools/mouth_middepth_rearrangement.cpp)

**Recurrence.** After six full material replacements, the finite-displacement wall run had a sideways RMS offset of 0.35963 bond lengths. At a fixed location, successive arrivals changed the pattern by 0.71926 RMS. Arrivals two sections apart reproduced it to within 1.91 × 10⁻¹² RMS over the final eight comparisons. All 16 runs showed the two-arrival period. Different initial disturbances could change which neighboring rows tilted the same way.

**Force measurement.** From each released state we held a localized ripple at zero, positive and negative amplitude and both half amplitudes, keeping the ripple fixed in space as material passed. All other sideways motions remained free. The table gives the late force **exerted by the patch** relative to the flowing, undisplaced copy, averaged over eight arrivals. Each ripple had the same envelope. The stiffness is the external holding-force difference divided by the displacement difference.

| Wall contact | Ripple wavelength in section spacings | Patch force at +0.0001 | Patch force at −0.0001 | Stiffness |
|---|---:|---:|---:|---:|
| Finite-displacement Morse | 8 | −0.125488 | +0.125488 | 1254.879 |
| Linearized | 8 | −0.139485 | +0.139485 | 1394.853 |
| Linearized | 16 | −0.047724 | +0.047724 | 477.237 |

With the finite-displacement wall, a second initial disturbance gave 1274.141. Halving the time step in the reference preparation gave 1254.879 to the precision shown. Across these three runs, halving the ripple amplitude reproduced the full-step stiffness to within 0.000037. Material bonds contributed 1014.985 to the reference stiffness, wall contacts 25.802 and maintaining the pattern as material flowed through it 214.092.

**Implication for 0PN.** The rearranged mid-depth material can transmit a local sideways load under controlled axial flow. To use this result on the first bend, project its actual displacement onto the rearranged material and predict the force reaching the wall. The wavelength dependence and the separate transport contribution matter for that prediction. The freely moving mouth's acceleration remains open.

Working base `49dca23e`. Independent local verdict for the finite-displacement wall runs: “Scoped clear for two-arrival recurrence and a positive local restoring coefficient under the nine-offset Morse wall law. The whole-mouth force remains open.” Formal project review and commit are pending.
