# MT15-05 Spike-3 Clean R_norm Direct-Derived Packet

Generated from `software/stage1_solver/mathematica/mt15_05_spike3_clean_rnorm.wls`.

**LOUD LABEL:** CLEAN structure packet: no posited `mixed_ports` in the `N0` path. The computed residual is a MACHINERY number on the M1b smoke background, not a physical falsification result until M1c.

## Inputs

| term | value | provenance |
| --- | --- | --- |
| K | 1.4830544424409784 | M1b mt15_02 wall weak form |
| M | 0.46742489897033024 | M1b mt15_02 wall weak form |
| B0 | 0.0009701851412361938 | M1b mt15_02 BdG moments |
| B2 | 0.00002229517295898049 | M1b mt15_02 BdG moments |
| B4 | 5.229950077800426e-7 | M1b mt15_02 BdG moments |
| Z0 | 2.6733784164555218e-6 | Spike-2 mt15_04 Green transfer |
| Z2 | -1.3784738025851937e-6 | Spike-2 mt15_04 Green transfer |
| Z4 | 4.4921775805827964e-7 | Spike-2 mt15_04 Green transfer |
| N0 | 4.947131039276808e-6 | Spike-2 mt15_04 retarded transfer |
| N2 | -2.5567729980147607e-6 | Spike-2 mt15_04 retarded transfer |
| N4 | 8.774217407344629e-7 | Spike-2 mt15_04 retarded transfer |

## Packet

- Direct-derived packet: `/var/projects/toy_physics/software/stage1_solver/mathematica/runs/mt15_05_spike3_clean_rnorm/mt15_05_v2_22b_direct_derived_packet.json`
- Forbidden-field teeth packet: `/var/projects/toy_physics/software/stage1_solver/mathematica/runs/mt15_05_spike3_clean_rnorm/mt15_05_v2_22b_direct_derived_forbidden_R_norm_packet.json`
- Diagnostics: `/var/projects/toy_physics/software/stage1_solver/mathematica/runs/mt15_05_spike3_clean_rnorm/mt15_05_direct_derived_diagnostics.json`
- `mixed_ports` key present in direct packet: `False`
- `derived_maxwell_transfer.status`: `derived_green_function_transfer`
- `derived_bdg_wall_coefficients.status`: `derived_bdg_wall_m1b`
- `Gamma_port`: `0.14814814814814817`

## Direct Formula Preview

These values are not exported in the V2-22B packet; the Python V2 chain remains the external residual judge.

- `D0 = K - B0 - Z0 = 1.4820815839213257`
- `R_norm = -0.1105891059875127`
- `R_pole = -0.6555153308620991`
- `P2 = 3.804508041557425e-7`
- `P4 = 4.99964161716701e-7`

## Term/Provenance Map

- M1b (`mt15_02_bdg_wall_derivation.wls`): `K`, `M`, `B0`, `B2`, `B4`, wall/profile/BdG provenance.
- Spike-2 (`mt15_04_spike2_transfer_n0.wls`): `Z0`, `Z2`, `Z4`, `N0`, `N2`, `N4`, `Gamma_port`, Green/gauge validation checks.
- V2-22B extension: additive direct-derived branch. Legacy `mixed_ports` packets still follow the old path; this packet omits `mixed_ports` and generates V2-21 `direct_coefficients`.

## V2 Chain

Run the V2-22C chain on the direct-derived packet to produce the authoritative observable packet. This report records the packet-side machinery preview only until that run is completed.