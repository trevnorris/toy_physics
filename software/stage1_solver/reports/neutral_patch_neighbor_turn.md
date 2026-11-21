# How neighboring brane material resists a local turn

The same stationary mouth at inward intake 0.24999999999821823 supplied four fixed-order preparations. Its material can turn and transmit torque while the opening stays still.

**Direct turn at finite temperature.** We applied equal and opposite moments to six inner sites per layer and their mobile neighbors. After one transit the mean relative angle per unit torque at orders 0.95, 0.50, 0.10 and 0.01 was 0.001542, 0.002131, 0.006334 and 0.012829 for central bonds. For layered bonds it was 0.000701, 0.000798, 0.001332 and 0.001643. All eight paired bath histories at each gate gave a positive increment. They share four parent labels, identical at entry. The two bond models have different six-site selections. The zero-torque control removes spontaneous motion. Half-strength drives and a halved timestep checked the signal. These readings measure finite-time angular compliance inside the patch. [Thermal torque drive](../tools/neutral_patch_warm_neighbor_torque.py)

**Internal torque.** Opposite small turns of four labeled inner sites per plane against their mobile neighbors gave these cold equilibrium stiffnesses, in patch energy per radian squared. Labels were identical across bond laws. [Microscopic torque test](../tools/neutral_patch_neighbor_torque.py)

| Material and turn prescribed | 95% order | 50% order | 10% order | 1% order |
|---|---:|---:|---:|---:|
| Central bonds, both layers measured together | 190.32 | 114.14 | 26.11 | 2.70 |
| Central bonds, relative turn measured within each layer | 683.56 | 679.20 | 675.28 | 674.40 |
| Layered bonds, both layers measured together | 281.66 | 156.21 | 32.81 | 3.32 |
| Layered bonds, relative turn measured within each layer | 661.99 | 651.92 | 638.98 | 635.71 |

The within-layer drive is insensitive to rigid motion of either layer. Only interlayer bonds weaken by assumption. Under sideways strain, the neighbor torque sometimes opposed and sometimes aided the induced turn across the eight thermal histories at each order (six inner sites per plane). [Direct torque on the saved states](../tools/neutral_patch_turning_warm.py)

**Force after one transit.** We bonded the whole patch’s mean angle to an overdamped sheet fixed at eight patch radii. A triangular-lattice shear test set the entry modulus. The sheet modulus scales with order, and site drag is one. The table gives the cold internal layer force per unit imposed strain after 5.29667 time units at half order. [Exterior torque calibration](../tools/neutral_patch_turning_neighbor.py) · [Coupled force calculation](../tools/neutral_patch_turning_hold.py)

| Surroundings | Central bonds | Layered bonds |
|---|---:|---:|
| Free edge | 4.961 | 5.691 |
| Sheet with the order-weighted curl energy | 9.235 | 5.801 |
| Sheet with symmetric-strain energy | 10.104 | 5.814 |

**Geometry.** In an infinite planar sheet, an edge turn can relax to a curl-free displacement proportional to $1/r$, even though the sheet resists over a finite transit. In three dimensions, a shell covering the poles of a sphere of radius $a$ supplies static torque $8\pi\mu a^3\theta/3$ under the same curl law. A co-rotating core of equal modulus adds $16\pi\mu a^3\theta/3$. Both controls suppress variation in the fourth coordinate. [Dimension and core checks](../tools/neutral_patch_turning_depth.py) · [Rotational reference requirement](../../../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md#r-s8-05--the-frame-the-branes-rotational-stiffness-is-measured-against)

**0PN requirement.** Determine which neighboring material stays shear-connected across the whole patch’s edge and faces as order falls. Measuring its relative turn, slip and transmitted torque would supply an independent input for the companion-driven load. The whole-mouth force and inertia remain open. [Comparison plots](figures/neutral_patch_neighbor_turn.png) · [Guarded results and review](../../../research/pde_ledger_v3/_measurements/neutral_patch_neighbor_turn.json)

Working base `49dca23e`. Independently reviewed locally within the scopes recorded in the measurements. Formal review and commit are pending.
