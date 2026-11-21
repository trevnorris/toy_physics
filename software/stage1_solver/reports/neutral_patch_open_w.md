# Does releasing the patch toward the bulk remove its shear stiffness?

**The patch rearranged when freed along $w$ and could still support a static sideways load.** The ordered bilayer originally required an outward holding force of 3.96969644 per projected area. On release the separation between the upper and lower halves, sorted by height, fell from 1 to 0.748–0.826 across the 13 trials. Individual sites could move farther apart. This cold patch has two periodic in-plane coordinates and one open coordinate $w$. Site count and projected area stay fixed. Volume density can change. There is no receiving bath, intake or imposed normal pressure. [Earlier confined test](neutral_patch_structural_endpoints.md)

**Material assumption.** We kept both existing Morse pair potentials. Each pair now takes a weighted mean of them based on its current separation in $w$, independent of its original layer. In units of the old gap let $t=\min(|\Delta w|,1)$ and $f=t^3(10-15t+6t^2)$. We compared weights $f$ and $f^2$ for the cross-layer potential. Both reproduce the old confined energy, in-plane forces and shear stiffness. They make different assumptions about material at intermediate separations. [Release calculation](../tools/neutral_patch_open_w.py)

| Initial 72-site patch | Confined | Released with $f$ | Released with $f^2$ |
|---|---:|---:|---:|
| Ordered | 87.091 | 85.006 | 59.622 |
| Defect preparation 191 | 75.426 | 19.464 | 9.384 |
| Defect preparation 817 | 17.666 | 34.445 | 24.289 |
| Defect preparation 1321 | 92.084 | 36.629 | −5.881 |

Entries are changes in shear stress per unit simple shear, in energy per projected area, after all internal forces settle. Each row uses the same small initial perturbation for the two release laws. Repeating with a different perturbation changed 817's response to 88.184 and 18.001, respectively. For the ordered start it gave 85.006 and 49.852. The negative slope in 1321 occurred with positive internal stiffness at fixed cell shape and a pre-existing shear stress of +0.60903649. It identifies an instability when the cell is allowed to shear. Each slope was checked against both stress differences and energy curvature at two strain amplitudes.

**Zero applied shear.** We let the cell itself shear at fixed projected area for preparation 1321. Under $f^2$ it settled at strain −0.06653563. The cell-averaged shear stress returned to zero within numerical precision and the new stiffness was +41.461. With $f$, the same preparation settled at −0.00443667 with stiffness +35.309. The negative slope identified a change of shape. The rearranged state again resisted shear. [Shape-release calculation](../tools/neutral_patch_shear_shape.py)

**Size and gravity.** An independent calculation prepared a flat crystal with each site above a gap in the opposite layer. Its small-cell stiffness was 85.006 for $f$ and 84.243 for $f^2$. In a 24-by-24 periodic cell a longer wave lowered that flat state's energy under both laws. Direct energy differences confirmed the instability. More than 98% of each mode's squared displacement lies along $w$. The next discriminating test is to let that wave settle and shear the buckled patch. The present tests supply neither a shear-free bulk reference nor the slow force transmitted to a draining mouth. Each phase received an independent local review. Formal review and commit are pending on working base `1dc45d09`. [Guarded measurements and independent reviews](../../../research/pde_ledger_v3/_measurements/neutral_patch_open_w.json)
