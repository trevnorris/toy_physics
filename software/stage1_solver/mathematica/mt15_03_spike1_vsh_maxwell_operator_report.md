# MT15-03 Spike-1 VSH Maxwell Operator Report

Status: Spike-1 only. No outgoing DtN/open-impedance solve, no wall-to-current transfer, no N0/Z0 extraction, and no V2 packet emission.

## Five-lane operator

Lane order: `(phi, a_r, a_E, a_B, a_w)` with `ell=2`, `L=ell(ell+1)=6`.

Target-blind placeholders: `R0(w)=a+(R_exit-a)(3x^2-2x^3)` with `a=1`, `R_exit=1.65`, and `Z(w)=Z_floor+Z_amp Exp[-((w-w0)/sigma)^2]` with `Z_floor=0.08`, `Z_amp=0.92`, `w0=L_w/2`, `sigma=0.38 L_w`.  These are geometry/localization placeholders only; the closed test grid is rectangular and avoids the origin (`r_min=0.35`) so Spike-1 can isolate the VSH operator without origin-parity or open-exit machinery.

The assembled action operator is `AAction=W_lane^-1 K`, where `K=sum sign M^dagger W M` over the D1 component matrices `G`, `E_r/E_E/E_B/E_w`, `C_r/C_E/C_B`, and `B_B/B_r/B_E`.  The primary MMS applies this weak Gram assembly directly to `u_exact`; the collocation `strongOp` generated from the Euler variation is retained as a secondary cross-reference.  The physical quotient norm uses the positive field-strength form `|E|^2+|B|^2+|C|^2`.

Cross-term carrying lanes: `a_r` row contains `-2 a_r/r^2 + 2 L a_E/r^2`; `a_E` row contains `+2 a_r/r^2`; `a_B` has the scalar vector-harmonic block only.  The smoke test deletes exactly this cross block.

## MMS

Primary gate: `AAction.u_exact-f_exact` is measured directly against the analytic Euler-Lagrange source on the interior rows where the centered weak Gram assembly represents a pointwise strong residual; two Dirichlet closure layers per side are excluded from this pointwise norm.  The full-grid `strongOp` cross-reference below is unchanged.  Manufactured solution has all five lanes nonzero and has coupled nonzero `a_r` and `a_E`; the source is derived by Euler variation of the same VSH weak density used by the matrix assembly.

| nr | nw | h | AAction phi | AAction a_r | AAction a_E | AAction a_B | AAction a_w | ASmoke phi | ASmoke a_r | ASmoke a_E | ASmoke a_B | ASmoke a_w |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 8 | 0.222222 | 0.059604801580567815 | 0.10115321848269036 | 0.06997364779399394 | 0.14060207593725887 | 0.189013068729257 | 0.059604801580567815 | 0.48041518135767336 | 0.1686727831458823 | 0.14060207593725887 | 0.189013068729257 |
| 12 | 12 | 0.153846 | 0.031019504588791033 | 0.05507346575224177 | 0.04219619097324425 | 0.06816604828460233 | 0.10983523700418085 | 0.031019504588791033 | 0.44988548026917463 | 0.13370875170350824 | 0.06816604828460233 | 0.10983523700418085 |
| 18 | 18 | 0.105263 | 0.014873798600273455 | 0.02781012143758397 | 0.021776180647156935 | 0.03298888331557206 | 0.059371219954522526 | 0.014873798600273455 | 0.4361822562237557 | 0.11916827305665606 | 0.03298888331557206 | 0.059371219954522526 |
| 26 | 26 | 0.074074 | 0.007496189591954471 | 0.014278190932017985 | 0.011164281038538256 | 0.01680091517022353 | 0.03164688219999131 | 0.007496189591954471 | 0.42878223866127746 | 0.1138516251323567 | 0.01680091517022353 | 0.03164688219999131 |

| lane | AAction orders | AAction final rel err | ASmoke orders | ASmoke final rel err |
| --- | --- | --- | --- | --- |
| phi | 1.776111, 1.936851, 1.949944 | 0.007496189591954471 | 1.776111, 1.936851, 1.949944 | 0.007496189591954471 |
| a_r | 1.653324, 1.800492, 1.897185 | 0.014278190932017985 | 0.178551, 0.081512, 0.048694 | 0.42878223866127746 |
| a_E | 1.375455, 1.743165, 1.901256 | 0.011164281038538256 | 0.631713, 0.303374, 0.129883 | 0.1138516251323567 |
| a_B | 1.968829, 1.912505, 1.920152 | 0.01680091517022353 | 1.968829, 1.912505, 1.920152 | 0.01680091517022353 |
| a_w | 1.476198, 1.621050, 1.790477 | 0.03164688219999131 | 1.476198, 1.621050, 1.790477 | 0.03164688219999131 |

Weak-operator teeth result: `aaction_teeth_pass=True`.  The genuine weak operator converges; deleting the D1 cross block from `AAction` leaves an O(1) residual on `a_r`/`a_E`, while `phi`/`a_B`/`a_w` remain second-order.

Weak-vs-strong consistency ratio `||AAction.u_exact-strongOp.u_exact||/||strongOp.u_exact||` decreases under refinement:

| nr | nw | h | ratio phi | ratio a_r | ratio a_E | ratio a_B | ratio a_w |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 8 | 0.222222 | 0.045754665888954694 | 0.08282229801936834 | 0.057329523614906494 | 0.10762211907264407 | 0.1530814917228946 |
| 12 | 12 | 0.153846 | 0.023790609470503057 | 0.04460872849246017 | 0.03368356124757608 | 0.05187627313386096 | 0.08604775751958267 |
| 18 | 18 | 0.105263 | 0.01143239896195668 | 0.02242560520893338 | 0.017228453171654303 | 0.02501161798530964 | 0.045885078429677845 |
| 26 | 26 | 0.074074 | 0.005770700625201258 | 0.011490377224321252 | 0.00879567180366958 | 0.012714430988643551 | 0.02439107619379279 |

| lane | weak-vs-strong ratios | ratio orders | final ratio |
| --- | --- | --- | --- |
| phi | 0.045755, 0.023791, 0.011432, 0.005771 | 1.778512, 1.931119, 1.945538 | 0.005770700625201258 |
| a_r | 0.082822, 0.044609, 0.022426, 0.011490 | 1.682693, 1.812240, 1.902952 | 0.011490377224321252 |
| a_E | 0.057330, 0.033684, 0.017228, 0.008796 | 1.446206, 1.766709, 1.913223 | 0.00879567180366958 |
| a_B | 0.107622, 0.051876, 0.025012, 0.012714 | 1.984540, 1.922374, 1.925461 | 0.012714430988643551 |
| a_w | 0.153081, 0.086048, 0.045885, 0.024391 | 1.566574, 1.656863, 1.798311 | 0.02439107619379279 |

Secondary `strongOp` cross-reference MMS (kept to preserve the original collocation Euler check):

| nr | nw | h | strongOp phi | strongOp a_r | strongOp a_E | strongOp a_B | strongOp a_w | strong smoke phi | strong smoke a_r | strong smoke a_E | strong smoke a_B | strong smoke a_w |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 8 | 0.222222 | 0.01597486780266993 | 0.02627949252476174 | 0.02276432508006703 | 0.037428304169532384 | 0.07384434008157693 | 0.01597486780266993 | 0.4292417817217282 | 0.11069368163188238 | 0.037428304169532384 | 0.07384434008157693 |
| 12 | 12 | 0.153846 | 0.007817647476418834 | 0.012925210705357652 | 0.011100250255929445 | 0.01834537699047001 | 0.03646487897912949 | 0.007817647476418834 | 0.42142233927411515 | 0.10934128975512929 | 0.01834537699047001 | 0.03646487897912949 |
| 18 | 18 | 0.105263 | 0.003708629893931507 | 0.006146274435323357 | 0.0052464957680725555 | 0.008699616511631135 | 0.01733993617125115 | 0.003708629893931507 | 0.41685939168262726 | 0.10908498955787577 | 0.008699616511631135 | 0.01733993617125115 |
| 26 | 26 | 0.074074 | 0.0018508846243665243 | 0.003070554764287982 | 0.002611088226880347 | 0.0043380937656308195 | 0.008652540490169685 | 0.0018508846243665243 | 0.41431637581708364 | 0.10903300907236665 | 0.0043380937656308195 | 0.008652540490169685 |

| lane | strongOp orders | strongOp final rel err | strong smoke orders | strong smoke final rel err |
| --- | --- | --- | --- | --- |
| phi | 1.943391, 1.965063, 1.977812 | 0.0018508846243665243 | 1.943391, 1.965063, 1.977812 | 0.0018508846243665243 |
| a_r | 1.929729, 1.958772, 1.974935 | 0.003070554764287982 | 0.049996, 0.028687, 0.017414 | 0.41431637581708364 |
| a_E | 1.953164, 1.974777, 1.985764 | 0.002611088226880347 | 0.033429, 0.006184, 0.001356 | 0.10903300907236665 |
| a_B | 1.939085, 1.966058, 1.980217 | 0.0043380937656308195 | 1.939085, 1.966058, 1.980217 | 0.0043380937656308195 |
| a_w | 1.918853, 1.958782, 1.978268 | 0.008652540490169685 | 1.918853, 1.958782, 1.978268 | 0.008652540490169685 |

## Gauge Deflation

- Gauge test grid: 10x9, raw dimension `450`.
- Pure-gauge rank: `90`; Lorenz/Gauss residual rank: `90`.
- Scalar gauge anchor: `chi(anchor_cell)=0` at cell `45`.
- Physical quotient dimension (`G=0` modulo pure gauge): `270`; expected `3N=270`.
- Pure-gradient max `|F|`: `1.887379141862766e-15`; positive physical field norm: `1.1556837592353772e-15`.
- Non-deflated gauge-fixed operator norm: `59.81320135249373`; after W-orthogonal gauge deflation, vector norm `9.612737051010784e-16`, operator norm `3.775450008221674e-14`.

The nonzero pre-constraint `G` norm (`11.011808393381854`) is expected: this is a pure gauge field-strength test, not a Lorenz-satisfying representative.  The physical quotient deflates it.

## Assembly Hermiticity Identity

`WLane.AAction = KAction = sum M^dagger W M` is Hermitian by construction, so this is an assembly/bookkeeping identity check for `WLaneInv`, not a physics gate.  Operator correctness is covered by the primary `AAction` MMS.  At `omega=0.37`, relative Frobenius norm of `W A - A^dagger W` is `6.217396463978332e-17`.

## Cross-Term Audit

- `a_r` cross-block norm on the MMS profile: `2.5303051904714082`.
- `a_E` cross-block norm on the MMS profile: `1.0265225657313926`.
- `a_B` cross-block norm on the MMS profile: `0.`.

## Term Map

| operator term | implementation | parent/design source |
| --- | --- | --- |
| parent localized Maxwell action | weighted Z(w) weak density with H=Z gauge fixing | notes/moving_throat_pde_program_compact.md:592-620 |
| exact localized Maxwell equation | operator is the Hessian of the VSH-reduced weak form, A=W^-1 K | notes/moving_throat_pde_program_compact.md:674-689 |
| mixed gauge invariants E_w,C_a | E_w=i omega a_w-d_w phi; C_r/C_E/C_B include d_w vector lanes | notes/moving_throat_pde_program_compact.md:769-786 |
| five-lane VSH decomposition | lane order (phi,a_r,a_E,a_B,a_w) | software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9829-9843 and decisions/05_mixed_maxwell_spike_design.md:12-19 |
| G,E,C,B reduced components | component matrices mG,mEr,mEE,mEB,mEw,mCr,mCE,mCB,mBB,mBr,mBE | software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9845-9862 |
| quadratic VSH weak form and angular weights | K=sum sign M^dagger W M with L weights on Psi/Phi components | software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9864-9872 and decisions/05_mixed_maxwell_spike_design.md:14-16 |
| pure-gauge quotient | GaugeMap chi=(i omega chi,d_r chi,chi/r,0,d_w chi); W-orthogonal deflation | software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9874-9879 |
| vector-Laplacian cross terms | CrossBlock rows a_r,a_E delete-test block; smoke A=A-crossBlock | software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9881-9890 and decisions/05_mixed_maxwell_spike_design.md:27-29 |
| closed Spike-1 BC and physical norm | closed Dirichlet MMS/assembly identity; positive F norm for quotient diagnostics | software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9912-9925 and software/stage1_solver/directives/mt15_03_spike1_vsh_maxwell_operator.md:29-48 |

Diagnostics JSON: `/var/projects/toy_physics/software/stage1_solver/mathematica/runs/mt15_03_spike1_vsh_maxwell_operator/mt15_03_spike1_vsh_maxwell_operator_diagnostics.json`.