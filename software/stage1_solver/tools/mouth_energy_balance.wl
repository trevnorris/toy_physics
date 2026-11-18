(* Independent Cartesian construction from the translated potential and
   unsteady Bernoulli pressure. No normal-motion condition at the mouth. *)
ClearAll["Global`*"];
$HistoryLength = 0;
Print["Wolfram version = ", $Version];
emit[label_, value_] := Print[label, " = ", ToString[Factor[value],InputForm]];
check[label_, lhs_, rhs_] := Module[{residual = FullSimplify[lhs-rhs]},
 Print[label, " | ", ToString[lhs,InputForm], " | ", ToString[rhs,InputForm],
       " | residual ", ToString[residual,InputForm]];
 If[residual =!= 0, Exit[1]]];

exterior[n_] := Module[{coords,rad,omega,normal,phi,velocity,dtphi,pressure,
 cutSpeed,energyDensity,energyNormal,momentumNormal,onSphere,angularArea,
 integrate,flux0,flux2,momentum1,dipole,normalDerivative,green,stored,storedRate,
 wrongPressure,wrongFlux,defect,cross},
 coords = Table[Unique["coordinate"],n];
 rad = Sqrt[coords.coords];
 omega = 2 Pi^(n/2)/Gamma[n/2];
 normal = PadRight[{Cos[theta],Sin[theta]},n];
 phi = q rad^(2-n)/((2-n) density omega)+eps (u First[coords]+b First[coords]/rad^n);
 velocity = Grad[phi,coords];
 (* Uniform lab stream is stationary. Only the source and decaying dipole
    translate; b changes with rate bd. The far lab pressure remains a datum. *)
 dtphi = eps bd First[coords]/rad^n-eps v D[phi-eps u First[coords],First[coords]];
 pressure = pinf+density (eps^2 u^2/2-dtphi-velocity.velocity/2);
 onSphere[expr_] := FullSimplify[expr/.Thread[coords->r normal],
                                Assumptions->r>0 && 0<theta<Pi];
 velocity = onSphere[velocity];
 pressure = onSphere[pressure];
 cutSpeed = velocity.normal-eps v Cos[theta];
 energyDensity = density (internal+velocity.velocity/2);
 energyNormal = Expand[energyDensity cutSpeed+pressure velocity.normal];
 momentumNormal = Expand[density First[velocity] cutSpeed+pressure Cos[theta]];
 angularArea = 2 Pi^((n-1)/2)/Gamma[(n-1)/2] r^(n-1);
 integrate[expr_] := FullSimplify[angularArea Integrate[Expand[expr] Sin[theta]^(n-2),
                                                   {theta,0,Pi}],Assumptions->r>0];
 flux0 = integrate[Coefficient[energyNormal,eps,0]];
 flux2 = integrate[Coefficient[energyNormal,eps,2]];
 momentum1 = integrate[Coefficient[momentumNormal,eps,1]];
 Print["dimension = ",n];
 emit["energy base",flux0];
 emit["quadratic total energy flux",flux2];
 emit["linear momentum flux",momentum1];

 (* Green's identity for the harmonic decaying dipole computes its kinetic
    energy from the two shell boundaries. Cross term with the uniform stream
    vanishes on every sphere; the radial base cross term is angularly odd. *)
 dipole = b Cos[theta]/r^(n-1);
 normalDerivative = D[dipole,r];
 green = integrate[density dipole normalDerivative/2];
 stored = Factor[(green/.r->outer)-(green/.r->a)];
 storedRate = D[stored,b] bd;
 cross = integrate[(Coefficient[First[velocity],eps,1]-u) u density];
 check["uniform-dipole kinetic cross term",cross,0];
 emit["localized kinetic energy in shell",stored];
 emit["localized kinetic energy rate",storedRate];
 check["shell energy conservation",(flux2/.r->outer)-(flux2/.r->a),-storedRate];
 wrongPressure = pressure+eps density bd Cos[theta]/r^(n-1);
 wrongFlux = integrate[Coefficient[Expand[energyDensity cutSpeed+wrongPressure velocity.normal],eps,2]];
 defect = Factor[(wrongFlux/.r->outer)-(wrongFlux/.r->a)+storedRate];
 emit["omit unsteady pressure, shell defect",defect];
 If[defect===0,Exit[1]];
 {flux0+(flux2/.bd->0),momentum1/.bd->0}
];

feed = exterior[3];
bulk = exterior[4];
(* Replace specific internal energy plus pressure/density by a still-open
   specific enthalpy. No common pressure or conversion energy is assumed. *)
feed = Simplify[feed/.{q->-md,u->ub,internal->hf-pinf/density}];
bulk = Simplify[bulk/.{q->md,u->uw,internal->hb-pinf/density}];
arrival = Factor[-First[feed]-First[bulk]];
force = Factor[-Last[feed]-Last[bulk]];
power = Factor[v force];
remaining = Factor[arrival-power];
emit["net energy arrival",arrival];
emit["ideal steady mouth force",force];
emit["mechanical power delivered to center or holder",power];
emit["remaining energy rate",remaining];
check["mouth-frame head difference",remaining,md (hf-hb+((ub-v)^2-(uw-v)^2)/2)];
emit["held brane stream, remaining rate",remaining/.{uw->0,v->0}];
emit["tow through resting media, remaining rate",remaining/.{ub->0,uw->0}];
emit["equal head, W=0, V=U diagnostic",remaining/.{hb->hf,uw->0,v->ub}];
