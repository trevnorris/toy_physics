# Order stress at a curved mouth with an open, flared receiver

At 0PN the equal-depth order energy gives a curved conversion layer tension. A stationary passive layer takes power from material that crosses it. Its tension and conversion drag describe the same normal force in the thin, locally equilibrated limit. Using that force again as the pressure that sustains the flow would count it twice. The accepted order balance leaves the conversion drive, relative order flux and bulk traction unspecified. [Order law](../../../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md) · [Earlier energy balance](../../../docs/_measurements/rim_curvature_head.txt)

## What fails at a tangent lip

Take one stationary, regular surface of fixed shear-order fraction $\chi$, with ordered feed on one side and shear-free material on the other. Let $n>0$ be constituent density, $\mu_\chi=\delta F/\delta \chi$ the order-energy derivative, and $M>0$ the mobility. Test the passive closure $\mathbf u\cdot\nabla\chi=-M\mu_\chi$, with conversion acting back on the same phase-carrying material by force density $\mu_\chi\nabla\chi$. Its power density is $-M\mu_\chi^2$. The specified bare order energy has equal minima, so ordered material entering steadily brings no bare endpoint energy advantage.

A smooth stationary potential velocity, with barotropic and quantum pressure, requires a single-valued specific head $B$. Its momentum equation is
$$
dB=(\mu_\chi/n)\,d\chi.
$$
Taking a second exterior derivative shows that $\mu_\chi/n$ must be constant along each connected regular order level. At a tangent reflecting lip the material velocity along the level can be finite, but its normal component vanishes. A bounded order gradient and finite positive limiting $Mn$ therefore give $\mu_\chi/n\to0$. Passive inward crossing at regular interior points would make that same quantity nonzero. This rules out the proposed simultaneous steady solution for this attachment and force law. The general order balance also permits independent conversion or transport relative to the carrier. Those terms require their own momentum and energy partners.

On the saved radius-one mouth at intake 0.25, a diagnostic logistic order layer of width $0.02R$ would require specific force work of 0.01558055 around a closed loop, where a scalar head gives zero. At width $0.04R$ the work was 0.01558100. These are evaluations on the original carrier, with the same tension and curvature mobility across widths. Changing width rescales the trial gradient coefficient and mobility at fixed tension and curvature mobility. The calculation evaluates the order force that the passive law would require. [Loop calculation](../tools/mouth_order_phase_work.py)

For outward normal pointing into the bulk and curvature $H=\nabla\cdot\mathbf N$, the thin-layer order force on the carrier is $-\sigma H\mathbf N$. The previous $5/24$ came from the small-slope law relating curvature to number current, with its trial coefficient $K=0.4$. Passive relaxation relates curvature to normal velocity, which also depends on the local density. The exact axial tension resultant on that same trial graph has magnitude 0.20807614. An opposing pressure or other traction and the energy to cross the layer must be accounted for independently. [Stress derivation](../tools/mouth_order_jump_check.py) · [Previous curved-mouth preparation](mouth_curved_phase.md)

## How the downstream shape changes the answer

The connected carrier solve uses $r=\sqrt{x^2+y^2+z^2}$ and bulk depth $w$, with spatial volume $4\pi r^2\,dr\,dw$. The feeding slab, rounded lip, radius-one opening and intake 0.25 are identical in all three rows. In the two new receivers a common collar expands from the inherited waist to radius 1. The conduit then flares to exit radius 1.5 over axial length 1.85, with either the project's historical cubic profile or a delayed flare. Both open across the same outlet shoulder into the same meshed finite chamber. Changing from immediate discharge to the conduit also changes the collar and outlet shoulder. The two conduit rows alone isolate the axial flare profile. [Historical geometry](../decisions/07_gate_a_freeze_sheet.md) · [Solver](../tools/mouth_flared_receiver.py)

| Receiver | Bath phase difference (rad) | Current crossing at center | Current crossing at $0.99R$ |
|---|---:|---:|---:|
| Previous rounded bend, immediate discharge | 0.41696727 | 0.18679798 | 0.000918654 |
| Earlier flare | 0.45014697 | 0.18613935 | 0.000918721 |
| Delayed flare | 0.45882304 | 0.18614168 | 0.000918721 |

All 2,000 sampled interior radii and all 51 shells carried inward flow in both new receivers, with whole-mesh single-phase certificates. On the finest mesh six axial sections per receiver recovered the target intake to within 0.114%. Moving the chamber roof from $10R$ to $7R$ changed the phase difference by about 0.00017 radians on the middle mesh. The paired change in center current between the flares was smaller than the change with mesh refinement in either absolute center reading. The phase values are stationary reservoir offsets required at fixed intake. They measure neither a pressure difference nor energy delivered to the order layer.

The rounded entrance has its narrowest section half a mouth radius down the $w$ axis, where its radius is $0.524837R$. In four spatial dimensions a round axial cut has three-volume $4\pi a(w)^3/3$. At fixed intake this waist raises the mean axial number current by a factor of 6.917 relative to a radius-one section. The physical location of the minimum radius remains uncalibrated. If the opening is itself the minimum, the present downstream numbers do not apply directly. [Geometry and section measurements](figures/mouth_order_and_flare.png)

**Next 0PN discriminant.** On an isolated mouth, determine the radius as a function of $w$ and where the surface of half shear order meets it. A freely varying square-gradient order field with zero wall energy gives $\partial_{\nu}\chi=0$ on the wall, hence orthogonal contact at a regular intersection. A tangent attachment needs independent support from wall physics. A measured contact angle, normal current, and traction on the same section would distinguish a changed attachment from power or stress supplied by an additional conversion channel. A two-mouth acceleration remains open.

Working result on base `49dca23e`, uncommitted. Independent local checks covered stress, phase compatibility and the geometry comparison. Formal project review is pending. [Guarded measurements](../../../research/pde_ledger_v3/_measurements/mouth_order_and_flare.json)
