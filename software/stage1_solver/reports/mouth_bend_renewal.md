# Does the restoring force survive material replacement?

**It can.** On the later-opening, deeply ordered flare, fresh material can sustain a sideways restoring force if its contacts carry their relaxed configuration downstream. If the contacts renew in place, their relaxed configuration gradually conforms to the held bend and the flow eventually amplifies it. Both tests used radius 1, inward intake 0.25, zero wave occupation, ordered shear modulus 1 and tangential attachment strength 1. [Static preparation](mouth_bend_static_shear.md)

An ideal marker moving with the archived carrier along the wall took **26.0591 model time units** to travel from the waist to the end of the deep bend. Set the mean time for replacing the relaxed contact configuration to one such transit for the comparison. Positive holding force in the table means that the mouth pushes back toward alignment. [Measured work and carrier transit](../../../research/pde_ledger_v3/_measurements/mouth_bend_renewal.json)

| Holding the deeper bend at +0.01, with the entrance bend fixed at zero | Deep holding force |
|---|---:|
| Immediately after elastic adjustment, either law | +0.00007070 |
| After one transit, contacts renew in place | −0.00003478 |
| After one transit, contacts move with the carrier, entering with the straight reference | +0.00006490 |

We let the shear field equilibrate in the four-dimensional guide and both reservoirs as the contact reference changed. In the local orthonormal wall frame, a rest offset \(p\) relaxes toward the actual tangential slip \(s\) at rate \((s-p)/\tau\). The two trials use the local time derivative or the derivative following a parcel. They keep the static carrier force and the same elastic response throughout. Fixed sites first lost a restoring direction at **0.53638** renewal times. For traveling contacts the steady force changed from amplifying to restoring in the weakest direction near **0.08403** transits of memory with the straight entry reference, or **0.08720** with zero traction at entry. These numerical crossings concern directional work on the two held bends. The latter entry comparison was evaluated at steady state. [Calculation](../tools/mouth_bend_renewal.py) · [Force histories](figures/mouth_bend_renewal.png)

**What to measure.** Mark contacts as they enter and measure how far their relaxed configuration travels. Halving the transport speed and doubling the renewal time preserved the steady response. A second discriminator is the deep holding force from a +0.01 entrance bend after settling. At one transit of memory it was +0.00009439 for straight-prepared entry and −0.00002735 for traction-free entry with the same meridional reference. [Readout](../tools/mouth_bend_renewal_readout.py)

For 0PN this supplies a possible mechanism for sustained internal alignment during drainage. Contact speed, renewal and birth references remain hypotheses. The transport law balances the loss of mean elastic energy through prescribed reference motion. Its power supply, reciprocal flow response and free acceleration remain open. Working base `49dca23e`. Two local reviewers checked the tracer clock and the held-force comparison within the recorded scopes. Formal review and commit are pending.
