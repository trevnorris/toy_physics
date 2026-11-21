(* Edge-supported order front in a three-ball. Pressures are local normal
   tractions, and the free surface spans a mouth in four spatial dimensions. *)
ClearAll["Global`*"];
$HistoryLength = 0;
require[statement_] := If[!TrueQ[statement], Print["failed: ", statement]; Exit[1]];

(* First extract the depth coefficients by integrating the equilibrium kink.
   z is measured in units of its width delta. *)
(* dchi/dz=-chi(1-chi) gives a finite polynomial integral on [0,1]. *)
shapeNorm = Integrate[chi (1-chi), {chi,0,1}];
Print["dimensionless depth integral = ", shapeNorm];
require[shapeNorm == 1/6];

(* Steady radial virtual work: energy under an allowed trial displacement
   h=b(1-r^2/R^2), verified below against the full EL equation. *)
assumptions = R>0 && sigma>0 && nu>0 && rho>0 && dp>0;
hTrial = b (1-r^2/R^2);
functional = 4 Pi Integrate[(sigma D[hTrial,r]^2/2-dp hTrial)r^2,{r,0,R}];
amplitude = b /. First[Solve[D[functional,b] == 0,b]];
h = Simplify[hTrial /. b->amplitude];
lapRad[profile_] := D[profile,{r,2}] + 2 D[profile,r]/r;
flow = Simplify[-nu lapRad[h]];
area = 4 Pi R^3/3;
edge = Simplify[-4 Pi R^2 sigma (D[h,r]/.r->R)];
mode = Sin[Pi r/R]/(Pi r/R);
eigenResidual = Simplify[nu lapRad[mode]+nu (Pi/R)^2 mode, assumptions];
tDecay = R^2/(Pi^2 nu);
sagOverTime = Simplify[flow tDecay / amplitude, assumptions];
Print["center sag = ",InputForm[amplitude], "  normal speed = ",InputForm[flow],
      "  edge normal force = ", InputForm[edge],"  speed times tau / sag = ",InputForm[sagOverTime]];
require[Simplify[dp+sigma lapRad[h],assumptions] == 0];
require[Simplify[edge - dp area] == 0 && eigenResidual == 0 && sagOverTime == 6/Pi^2];
(* Independent radial discretization in y=r h, y(0)=y(R)=0, linear
   forcing on uniform grid. Read sag from a one-sided cubic-exact derivative
   of y, and decay time by an actual symmetric tridiagonal eigensolve. *)
sigma0=1/6; nu0=2/5; rho0=6/5; pressure0=1/5000;
Do[
  dr = N[radius/cells];
  stencil = SparseArray[{Band[{1,1}]->2.,Band[{1,2}]->-1.,Band[{2,1}]->-1.},
                        {cells-1,cells-1}]/dr^2;
  y = Prepend[Append[LinearSolve[stencil,Table[pressure0 i dr/sigma0,{i,cells-1}]],0.],0.];
  sag = (-11 y[[1]]+18 y[[2]]-9 y[[3]]+2 y[[4]])/(6 dr);
  rate = Min[Eigenvalues[nu0 Normal[stencil]]];
  dRead = rho0 (4 Pi radius^3/3) (6/Pi^2) sag rate;
  dTheory = N[rho0 (4 Pi radius^3/3) nu0 pressure0/sigma0];
  Print["R=",radius," cells=",cells," sag=",InputForm[sag]," tau=",InputForm[1/rate],
        " inferred_D=",InputForm[dRead]," predicted_D=",InputForm[dTheory],
        " relative_error=",InputForm[dRead/dTheory-1]];
  require[Abs[dRead/dTheory-1]<2/cells^2],
 {radius,{24,48}}, {cells,{32,64,128}}];
(* Tangential control. If pressure stays spatially uniform, a uniform normal
   lab velocity is incompatible with a stationary non-flat front and nonzero
   tangential advection. If normal lab velocity stays uniform instead, the
   first-order change of shape and pressure is odd across the mouth. *)
cartH = dp (R^2-x^2-y1^2-y2^2)/(6 sigma);
shapeOdd = tVel dp x (R^2-x^2-y1^2-y2^2)/(30 sigma nu);
pressureOdd = -sigma Laplacian[shapeOdd,{x,y1,y2}];
mixedTest = Simplify[nu Laplacian[shapeOdd,{x,y1,y2}] - tVel D[cartH,x]];
Print["uniform pressure: required shift in normal velocity = ",InputForm[tVel D[cartH,x]],
      "  uniform normal velocity: pressure dipole = ", InputForm[Simplify[pressureOdd]]];
require[mixedTest == 0 && Simplify[pressureOdd+(pressureOdd/.x->-x)] == 0];
(* A compliant edge changes both the absolute sag and the decay eigenvalue.
   Bi=1 yields y'(R)=0 at y=r h during decay; sag is measured relative to
   the normal anchor position, not the displaced edge. *)
softH = h + dp R^2/(3 sigma);
softMode = Sin[Pi r/(2 R)]/(Pi r/(2 R));
softBoundary[f_] := Simplify[(D[f,r]+f/R)/.r->R];
gainRatio = Simplify[(6/Pi^2) (softH/.r->0)*nu (Pi/(2 R))^2/flow];
require[softBoundary[softH] == 0 && softBoundary[softMode] == 0 && gainRatio==3/4];
Print["Bi=1: uncorrected absolute-sag intake / true intake = ",gainRatio];
(* If h vanishes on the edge, integration by parts on B3 turns the
   projected pressure force -integral p h_x into +integral h p_x. *)
pressureForceRatio = Simplify[4 Pi Integrate[h r^2,{r,0,R}]/(area amplitude)];
Print["prescribed pressure-gradient force per gradient, area and sag = ",pressureForceRatio];
require[pressureForceRatio == 2/5];
(* A flux concentrated near the center probes the shape identifiability.
   Use s=r/R and dimensionless maximum through-speed one. *)
profile[s_] := (1-s^2)^3;
cu = Integrate[s^2 profile[s],{s,0,1}];
centerWeight = Integrate[(s-s^2) profile[s],{s,0,1}];
Print["weighted normal flow mean = ",3 cu,"  dimensionless center sag = ",centerWeight,
      "  center estimate / total flow = ",6 centerWeight/(3 cu)];
require[6 centerWeight/(3 cu) == 187/64];
Exit[0];
