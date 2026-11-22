# Does the independently prepared receiver lose its sideways restoring force?

For Newtonian free fall, the material leaving the brane must be able to receive momentum. Its shear response constrains the load that can return to the mouth. The interface traction, outgoing momentum and nearby storage also enter the force balance. A rigid crystal and a dense material whose constituents exchange neighbors can give different answers even at the same temperature. We tested whether preparing the same observable material state on an independent geometry reproduced the earlier pressure and shear response. [Earlier measurement](mouth_pressure_matched_receiver.md)

**Independent preparation.** We used the 128-site sample that retained crystalline order in the earlier test. Before measuring the shear response of the selected preparation, we prescribed three heating excursions, a return to thermal energy 1.28 at zero imposed external pressure, and a choice based on density, structural order, site exchange and connectedness. Only the excursion to 1.408 bond depths met the criteria. Its returned density was 1.06691, all 128 sites belonged to one connected cluster at the final four readings, and its strongest density-wave intensity went from 0.08502 to 0.07882. The excursions to 1.536 and 1.664 returned at densities 0.59059 and 0.17687, respectively, below the prespecified 0.95 minimum. The original sample and this sample supply two founding geometries. [Predeclared state test](../tools/mouth_receiver_precursor_match.py)

**Same material under load.** We cooled that sample to 1.088 and compared opposite sideways strains in the laboratory plane spanning the brane direction and the fourth direction. For each strength, subtracting the negative-strain stress from the positive one measures the part that changes sign with the strain. The table gives that late response divided by the initial response. The pressure tests measured the change in log volume per unit applied pressure. [Signed force readout](../tools/mouth_receiver_matched_response.py)

| Protocol on the same starting geometry | Full strain ±0.02 | Half strain ±0.01 | Compression compliance |
|---|---:|---:|---:|
| Standard thermal path A | +0.05302 | −0.33843 | 0.11243 |
| Standard thermal path B | +0.12984 | −0.04530 | 0.09373 |
| Half integration step | −0.02928 | −0.12797 | 0.15078 |
| Half pressure-controller mobility | +0.05978 | +0.10834 | 0.10799 |

All 28 arms satisfied the finite-window density, connectivity, site-motion and endpoint order tests. The earlier crystalline preparation of this same founding sample retained about four-fifths of the initial shear response. The measurable state of the sample therefore distinguishes two mechanically different preparations at the same bath temperature and imposed pressure. The half-strain readings still do not bound the residual within ±0.25. A negative individual estimate carries the uncertainty of subtracting fluctuating stresses. The four rows share a geometry, and changing the time step or pressure controller can change the microscopic path as well. [Measurement and independent audit](../../../research/pde_ledger_v3/_measurements/mouth_independent_receiver.json)

**Follow-up at one nominal transit.** The short force holds lasted 48.88171 native units. Four receiver preparations and two crystalline controls subsequently retained their respective structures at every saved reading through an undriven hold of 840.17332 native units. This represents one mouth transit if the earlier calibration of site mass, pitch and binding energy transfers. The longer signed measurements found small residual receiver shear at both strengths. The physical conversion mechanism and force at a mobile interface remain open. [Full-transit results](mouth_receiver_transit.md) · [Free-fall roadmap](../../../docs/zero_pn_free_fall_roadmap.md)

Working baseline `a91a633f`. Formal repository review pending.
