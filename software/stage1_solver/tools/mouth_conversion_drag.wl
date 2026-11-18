(* A compatibility calculation with supplied ideal exterior momentum rates.
   This does not solve the interface PDE or determine a physical inertia. *)
ClearAll["Global`*"];
Print["Wolfram version = ", $Version];
check[label_, lhs_, rhs_] := Module[{residual = FullSimplify[lhs-rhs]},
  Print[label, " | ", ToString[lhs, InputForm], " | ", ToString[rhs, InputForm],
        " | residual ", ToString[residual, InputForm]];
  If[!TrueQ[And @@ (#===0& /@ Flatten[{residual}])], Exit[1]]];

(* An oriented incidence matrix for the two material subsystems. Each column
   is a directed momentum transfer; the conversion and inter-medium traction
   each enter one subsystem and leave the other. Outer entries are supplied. *)
incidence = {{1,0,-1,-1,-1,0,-1,0},{0,-1,1,1,0,-1,0,-1}};
rates = {md (ub-v),md (uw-v),carriage,transfer,ff,fb,kf,kb};
budget = incidence.rates;
Print["oriented subsystem budgets = ", ToString[budget,InputForm]];
localForces = First@Solve[Thread[budget=={0,0}],{ff,fb}];
Print["required subsystem reactions = ",ToString[localForces,InputForm]];
total = Factor[ff+fb/.localForces];
Print["mouth force = ",ToString[total,InputForm]];
check["internal exchange cancellation",Take[Total[incidence],{3,4}],{0,0}];
check["tow through resting media",total/.{ub->0,uw->0,kf->0,kb->0},0];
check["held brane stream",total/.{uw->0,v->0,kf->0,kb->0},md ub];

(* The local first variation includes the change in the area measure, too. *)
localTransfer = (j0+eps j1) (vx0+eps vx1) (area0+eps area1);
firstTransfer = D[localTransfer,eps]/.eps->0;
Print["live carried first variation including measure = ",ToString[firstTransfer,InputForm]];

(* Construct the upward oriented area of a graph in four coordinates by its
   3x3 minors, before choosing any pressure or geometry perturbation. *)
tangents = {{1,0,0,hx},{0,1,0,hy},{0,0,1,hz}};
areaVector = Table[(-1)^(i+4) Det[Delete[#,i]& /@ tangents],{i,4}];
Print["upward graph area vector = ",ToString[areaVector,InputForm]];
pressureX = Expand[(-pressure areaVector)[[1]]];
check["normal pressure projected along x",pressureX,pressure hx];

(* Differentiate the live graph in Cartesian coordinates, then restrict to
   a sphere. Differentiate pressure times graph slope at eps=0. *)
rad = Sqrt[x^2+y^2+z^2];
height = h0[rad]+eps h1[rad] x/rad;
pressureField = p0[rad]+eps p1[rad] x/rad;
cartesianForce = D[pressureField D[height,x],eps]/.eps->0;
angularForce = FullSimplify[cartesianForce/.{x->r mu,y->r Sqrt[1-mu^2],z->0},
                           Assumptions->r>0 && -1<mu<1];
shellForce = FullSimplify[2 Pi r^2 Integrate[angularForce,{mu,-1,1}],Assumptions->r>0];
Print["projected pressure per radial interval = ",ToString[shellForce,InputForm]];
endpoint = 4 Pi p0[r] r^2 h1[r]/3;
remaining = 4 Pi r^2 (p1[r] h0'[r]-p0'[r] h1[r])/3;
check["pressure moment with endpoint",shellForce,D[endpoint,r]+remaining];
flat = shellForce/.{h0'[r]->0,h1[r]->0,h1'[r]->0};
check["flat unchanged face",flat,0];

(* A possible force moment, without asserting that these functions solve
   either fluid problem: integrate Cartesian pressure*x-slope directly. *)
momentDensity = kp x D[kh (x^2+y^2+z^2)/2,x];
example = Integrate[
  (momentDensity/.{x->r Cos[th]}) r^2 Sin[th],
  {r,0,a},{th,0,Pi},{phi,0,2 Pi},Assumptions->a>0];
Print["curved face pressure moment = ",ToString[example,InputForm]];

(* Same-family scalar slow-response limit only; inertia is independent data. *)
sourceVelocity = ds/(4 Pi density distance^2);
acceleration = dt sourceVelocity/mt;
coupling = Factor[acceleration distance^2/ms];
Print["mass-normalized inverse-square coefficient = ",ToString[coupling,InputForm]];
Print["common intake per inertia condition = ",
      ToString[Factor[coupling/.{dt->lambda mt,ds->lambda ms}],InputForm]];
