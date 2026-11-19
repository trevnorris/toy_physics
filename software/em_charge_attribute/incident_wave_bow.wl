(* Vary a displaced radial potential by translating the integration variable
   into the wall frame.  Integrating by parts in the radius then gives the
   kernel multiplying V'(r). Boundary terms vanish for a localized potential
   contrast (subtract its constant exterior value). *)
ClearAll["Global`*"];
emit[name_, v_] := Print[name, " = ", ToString[v, InputForm]];
assum = r > 0 && k > 0 && Element[{alp, dx}, Reals];
sphereSquare = Integrate[2 Pi Cos[k (r u + dx) + alp]^2, {u,-1,1},
    Assumptions->assum];
radialForceWeight = FullSimplify[-r^2 (D[sphereSquare,dx]/.dx->0)/2,
    Assumptions->assum];
(* F = integral V(r) weight dr = -integral V'(r) primitive dr. *)
kernel = FullSimplify[-Integrate[radialForceWeight /. r->v,{v,0,r},
    Assumptions->assum], Assumptions->assum];
angle = FullSimplify[(2 kernel/r^2 /. k->z/(2 r)),
    Assumptions->z>0 && r>0 && Element[alp,Reals]];
emit["angular kernel from a rigid shift of the incoming trace",angle];
emit["kernel at z=Pi and phases +Pi/4, -Pi/4, 0",
    (FullSimplify[angle /. {z->Pi,alp->#}]& /@ {Pi/4,-Pi/4,0})];
emit["long wavelength angular kernel",Limit[angle,z->0]];

(* Solve the local dipole equation at rest, then Fourier analyze the answer.
   Its two participating Neumann modes can be checked in the full equation. *)
profile[y_] := Cos[Pi y/2]^2;
local[y_] := peak profile[y]/inertia;
av[expr_] := Integrate[expr,{y,-1,1}]/2;
shift = av[local[y]];
bow = av[local[y] Cos[Pi y]]/av[Cos[Pi y]^2];
emit["uniform and first even Neumann accelerations",{shift,bow}];
emit["midplane minus mean end acceleration",local[0]-(local[1]+local[-1])/2];
emit["uniform depth control modes",{av[peak/inertia],
    av[peak/inertia Cos[Pi y]]/av[Cos[Pi y]^2]}];
residual=FullSimplify[local[y]-shift-bow Cos[Pi y]];
emit["remaining force density divided by inertia",residual];

(* Use homogeneous two-boson polynomials in the Bargmann representation.
   Occupation basis entries multiply p^2/Sqrt[2], p m, m^2/Sqrt[2]. *)
quadrature[f_] := Expand[(p D[f,m]-m D[f,p])/(2 I)];
basis={p^2/Sqrt[2],p m,m^2/Sqrt[2]};
coords[f_] := {Sqrt[2] Coefficient[f,p,2],Coefficient[Coefficient[f,p,1],m,1],
    Sqrt[2] Coefficient[f,m,2]};
qmat = Transpose[coords[quadrature[#]]& /@ basis];
q2 = qmat.qmat;
emit["one photon each: odd quadrature mean and second moment",{qmat[[2,2]],q2[[2,2]]}];
(* The three spectral projectors follow from interpolation of the measured
   characteristic polynomial, without choosing eigenvector phases. *)
roots=Sort[lambda/.Solve[CharacteristicPolynomial[qmat,lambda]==0,lambda]];
projector[t_] := Fold[Dot,IdentityMatrix[3],
   ((qmat-# IdentityMatrix[3])/(t-#)& /@ DeleteCases[roots,t])];
weights=Table[{t,projector[t][[2,2]]},{t,roots}];
emit["spectrum and one-each weights",weights];
lower[f_] := Exp[I th] D[f,p]+Exp[-I th] D[f,m];
amp=FullSimplify[lower[lower[p m]]];
emit["two-photon coherent absorption amplitude, common envelope removed",amp];
emit["its angular dipole",Integrate[2 Pi u (amp /. th->kr u),{u,-1,1}]];
If[!TrueQ[FullSimplify[D[kernel,r]+radialForceWeight==0,Assumptions->assum] &&
    residual==0 && Total[weights[[All,2]]]==1],Exit[1]];
