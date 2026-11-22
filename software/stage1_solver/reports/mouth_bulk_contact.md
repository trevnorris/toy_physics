# Does converted material release shear when it meets another patch?

**Two converted patches formed 48 contacts and continued to transmit a sideways force.** In four unstrained runs no constituent was observed crossing into the other half. Contact alone failed to produce a shear-free receiver on the measured timescale.

Each patch had 32 mobile sites, two periodic tangent directions and a free normal direction. We prepared them separately on the radius-one conversion profile at intake 0.25, then brought their faces together at temperature 0.008, zero applied outer normal traction and 2.32093 sites per projected area per patch. The same pair law acted across the contact. Disconnected controls used the same noise. We allowed one downstream transit for contact formation and one for the measurement, each lasting 48.8817124312 model time units. Intake set the preparation and clock. There was no imposed throughflow during contact. [Passing-patch preparation](mouth_brace_grip.md) · [Contact calculation](../tools/mouth_bulk_contact.py)

| Area-preserving stretch of the whole sample | Joined | Disconnected |
|---|---:|---:|
| Final-quarter incremental stiffness relative to initial | 87.03–95.20% | 90.96–93.02% |
| Cross-patch contacts throughout the sampled hold | 48 | 0 |

All 336 starting neighbor pairs in each joined sample were retained at every reported hold frame. Contacts formed below separation 1.15 and persisted up to 1.30. Crossing required reaching 0.2 past the plane through the combined center of mass, checked at intervals of at most 0.128. Shorter excursions remain unresolved. The trials include two pair laws, a change of alignment and seeds, and a smaller time step with a different noise history.

**Force through the contact.** At one common alignment and seed pair, we shifted the two original populations sideways relative to each other by ±0.0001 and ±0.00005. Only their sideways centroids were held. The signed derivative of the upper holding traction with relative displacement gives the stiffness in patch force per area per unit displacement. Disconnecting the patches at the same positions gave zero transmitted force. Both tables average the final 12 of 48 post-step readings. [Relative-displacement calculation](../tools/mouth_bulk_contact_slip.py)

| Depth-dependent mixing weight | Initial stiffness | Final-quarter stiffness | Retained |
|---|---:|---:|---:|
| First law, $f$ | 74.94879 | 49.90095 | 66.58% |
| Second law, $f^2$ | 62.00835 | 21.83864 | 35.22% |

**Next test for 0PN.** Both laws retain a preferred axis after conversion. A unit-length bond has energy −1 along the old slab and −0.4778574362 along its normal. Test an isotropic bulk interaction on independently prepared condensed and disordered samples at the same bulk conditions. Measure whether it resists compression while allowing shear to decay. A passing law can then be tested for transmission of lip force to the flare. The extended bulk’s density and pressure, its optical response, and free-mouth acceleration remain open.

Working checkpoint `a91a633f`. Separate local reviewers checked the whole-sample and contact-force measurements. Formal repository review is pending. [Measurements and independent checks](../../../research/pde_ledger_v3/_measurements/mouth_bulk_contact.json)
