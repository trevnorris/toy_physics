# Does it cost work to deepen a wave-supported mouth?

**Result.** In a conservative one-mouth test with the rear support anchored to the bulk, deepening the surrounding brane compressed the standing wave and increased its energy. The slab, wall and rear support together released more energy than the wave gained. The reverse drive cost work. A rear reference moving exactly with the surrounding brane gave zero excess work. The physical rear support is unknown. [SymPy and Mathematica runs](../../../docs/_measurements/supported_mouth_bend.json)

**Test.** Take positive \(w\) into the throat. Impose a common displacement \(b\) on both faces at a spherical surface of radius \(R\) well outside the mouth radius \(a\). Keep the slab flat at infinity. The mouthless control has the same exterior displacement and a flat, translated interior. Subtract its work and the zero-drive baseline. For a midplane gradient energy with stiffness \(\kappa>0\), the annulus \(a<r<R\) contributes \(k(z-b)^2/2\), where \(z\) is the rim depth and \(k=4\pi\kappa/(1/a-1/R)\). [Two-face slab model](two_face_rim_bend.md)

Let the reflecting end at \(Y\) move and put \(L=Y-z\). Use the straight-wall energy \(TL\) and a weakly exchanging standing mode at fixed action \(J\), with energy \(J\sqrt{\omega_*^2+\nu^2/L^2}\). Their sum is \(V(L)\). Choose a rear spring of stiffness \(k_b>0\) with reference \(Y_*+\eta b\) as an explicit test assumption. Relax both \(z\) and \(Y\) at every drive. This tests stored work at fixed radius and wave preparation. The calculation includes no discharge, material conversion or pump work. [Wall load](common_face_work.md) · [Wave and endpoint forces](bend_front_attachment.md)

At equilibrium the net tensile load is \(P=V'(L)=k(z-b)\). The excess external work satisfies \(dW/db=(\eta-1)P\), counting work at both controls when the rear reference moves. An isolated positive rim depth \(X\) gives \(P_0=kX>0\). A bulk-fixed reference (\(\eta=0\)) therefore gives a negative slope for **every** such equilibrium in this model. In that case \(d^2W/db^2=(1/k+1/V''+1/k_b)^{-1}>0\). Positive stiffness and a negative initial work slope coexist because the system is already loaded.

**Four finite drives.** Dimensionless example with \(k=k_b=1,\ T=2,\ J\nu=1,\ \omega_*=0,\ Y_*=3,\ \eta=0\), initially \(z=L=1\). SymPy minimized the stored energy locally by solving both force balances on the stable branch. Mathematica integrated the independently calculated drive reaction. All four rows agree to ten printed decimal places.

| Broad bend \(b\) | Rim shift | Reflector shift | Wave energy change | Excess work |
|---:|---:|---:|---:|---:|
| +0.02 | +0.0119904408 | +0.0080095592 | +0.0039967924 | −0.0199199362 |
| −0.02 | −0.0120096411 | −0.0079903589 | −0.0040031922 | +0.0200799358 |
| +0.01 | +0.0059976051 | +0.0040023949 | +0.0019991990 | −0.0099799920 |
| −0.01 | −0.0060024051 | −0.0039975949 | −0.0020007990 | +0.0100199920 |

**Physical question.** Does the throat constrain a material deformation that requires positive *total* work when an incident bend deepens it? The present solver holds the screen geometry fixed and supplies neither the axial reflector motion nor the matched rim and bulk work. This one-mouth calculation determines the sign for the specified supports. The interaction between two physical mouths remains open. Working calculation on base e0a94594, uncommitted pending independent review. [Solver readouts](rim_bend_readouts.md)
