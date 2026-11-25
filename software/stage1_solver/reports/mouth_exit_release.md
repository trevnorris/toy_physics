# What happens when the downstream exit is free?

**The bonded mouths at radii 1 and 2 narrowed and continued to drain when four coordinates were released. Allowing independent changes in crystalline order exposed lower-energy configurations at both radii. Whether either exit ultimately stays open or heals remains unresolved.**

The reversible order and carrier equations leave the identity of newly bonded lattice points unspecified. We therefore retained each core’s original bond to its collar as a supplied preparation. The entrance radius, reservoir phase difference and stored wave action of one were held fixed. [Bond audit](mouth_exit_bond_audit.md)

| Entrance | Fine-mesh exit before → after | Depth after | Inward current after | Cavity share of wave kinetic energy before → after |
|---:|---:|---:|---:|---:|
| 1 | 0.73437 → 0.44978 | 0.69226 | 0.28639 | 92.87% → 90.76% |
| 2 | 1.28925 → 1.03697 | 0.82834 | 1.32343 | 86.69% → 82.81% |

Both meshes found a local minimum in exit radius, depth, waist width and one order variation. The coarse exits were 0.46202 and 1.01986. At the fine minima the elastic energy slopes per unit log exit radius were −2.84427 and −21.66036, while the wave slopes were −0.15456 and −0.39925. The strained crystal supplied most of the outward force. [Signed work and wave losses](../../../research/pde_ledger_v3/_measurements/mouth_exit_release.json)

### What happened when interior order could vary independently?

Each run started from its four-coordinate minimum. All four stopped during order relaxation, before the first shape update. The work columns give the change in energy in separate tests that narrowed or widened the starting exit by 0.05%, holding the saved order fixed. A negative entry is an energy decrease. The order residual target was 0.00002.

| Entrance / mesh | Last saved order iteration | Order residual | Work on narrowing | Work on widening |
|---|---:|---:|---:|---:|
| 1 / 48 | 49 | 4.28842 | −0.00036094 | +0.00036271 |
| 1 / 72 | 100 | 19.29603 | −0.00041082 | +0.00041251 |
| 2 / 48 | 1,190 | 19.78155 | −0.00708096 | +0.00709627 |
| 2 / 72 | 31,475 | 0.63932246 | Unmeasured at stop | Unmeasured at stop |

The two radius-1 runs exhausted their order line searches. On the coarse radius-2 mesh the shear-bearing region split, leaving the continuation of the occupied mode ambiguous. On the fine radius-2 mesh the user stopped the run at exit 1.03697383. Its wave still occupied one connected shear-bearing region. Inward current was 1.3190232126 and terminal current was 1.3190232153. Across the last 1,000 order iterations the interior energy fell by only 1.53 × 10⁻⁸, with the sampled order residual between 0.53512 and 0.67656. [Final measurements and checkpoint hashes](../../../research/pde_ledger_v3/_measurements/mouth_exit_release.json)

Separate signed tests locate the physical evidence for narrowing. At fine radius 1 a joint shape and order change reached exit 0.44828 and lowered energy by 0.01635. At fine radius 2, **at iteration 2,325**, narrowing by about 0.1% lowered energy by 0.01672 with order held and by 0.01681 with a local order change. The final checkpoint has no signed exit measurement. The model fixes the terminal order at one half and uses a strictly positive exit coordinate. Establishing arrest requires restoring work about an open equilibrium under the remaining material and shape variations. Establishing closure requires a description that admits contact and tracks the occupied wave through loss of shear connectivity. The free-release and blind-radius gates remain closed for this branch. [Wave and closure audit](mouth_trapping_closure_gate.md) · [0PN roadmap](../../../docs/zero_pn_free_fall_roadmap.md)

**Independent review.** “Scoped clear for the archived bonded-exit readout at radii 1 and 2. Both four-coordinate minima drained, and the signed narrowing tests lowered energy at the sampled states. All four unrestricted runs stopped before a geometry update. The final fine radius-2 checkpoint establishes continued current and unresolved order stationarity. Neither self-support nor healing was established.”
