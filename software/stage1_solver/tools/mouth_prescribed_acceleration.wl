(* Independently differentiate the incoming momentum flux at fixed carrier
   field on the actual 3D face. The grid and translated reference density
   match the time integrator's trigonometric interpolant. *)
raw = "_scratch/s11c/mouth_rim_push_sympy01/resting_face.f64";
If[FileHash[raw,"SHA256","HexString"] =!=
  "e76812f953ed1af4fd216a5d7fb9d70da6801536dfc786a48fcea165c8ca8053",Exit[1]];
{re,im,j} = ArrayReshape[BinaryReadList[raw,"Real64",ByteOrdering->-1],{3,32,32,32}];
dx = 1/8; length = 32 dx; axis = N[Range[-31,31,2]/16];
ball = Table[i^2+k^2+l^2,{i,axis},{k,axis},{l,axis}];
mask = Map[Boole[#<1]&,ball,{3}];
core = Map[Max[0,1-#]&,ball,{3}];
g = core^3;
dg = 6 MapThread[#1*#2&,{axis,core^2}];
den = re^2+im^2;
integ[arr_] := dx^3 Total[Flatten[arr]];
der[arr_] := (RotateLeft[arr]-RotateRight[arr])/(2 dx);
spectralD = N[Table[If[i==k,0,Pi/length*(-1)^(i-k)*Cot[Pi*(i-k)/32]],
    {i,32},{k,32}]];
densityGradient = ArrayReshape[spectralD.ArrayReshape[den,{32,1024}],{32,32,32}];
rate = integ[j];
dRaw = dg+g*densityGradient/den;
dj = rate*(dRaw-g*integ[dRaw]/integ[g])/integ[g];
vx = (re*der[im]-im*der[re])/den;
dNormalLog = -2*dg*(1-1/Sqrt[den]);
frozenSlope = integ[mask*dj*vx]+integ[dNormalLog*der[den*mask]]/4;
Print["fixed-field incoming load derivative per sideways displacement = ",InputForm[frozenSlope]];
If[!TrueQ[Abs[rate-.25]<10^-13 && Abs[integ[dj]]<10^-13 && frozenSlope!=0],Exit[1]];

(* Compare control-subtracted loads for four independently evolved fields.
   Only this cross-check consumes the run outputs. *)
readCase[tag_,acc_] := Module[{dir,out,data},
  dir = "_scratch/s11c/mouth_prescribed_accel_"<>tag<>"02/";
  out = Import[dir<>"outcome.json","RawJSON"];
  If[!TrueQ[out["limitsVerified"] && out["exitCode"]==0 &&
     out["childOutcome"]["exitCode"]==0 && out["stderrBytes"]==0],Exit[1]];
  data = ImportString[ReadString[dir<>"stdout"],"RawJSON"];
  If[!TrueQ[data["acceleration"]==acc && data["grid"]=={65,32,32,32}],Exit[1]];
  data];
amps = {0,.02,-.02,.01,-.01};
cases = MapThread[readCase,{{"zero","positive","negative","positive_half","negative_half"},amps}];
loads[item_,steps_] := Lookup[item[steps],"feed_load"];
times = Lookup[First[cases]["32"],"t"];
If[!TrueQ[Max[Abs[times-{.00025,.0005,.001}]]<10^-15],Exit[1]];
increments = Table[loads[item,"32"]-loads[First[cases],"32"],{item,Rest[cases]}];
odds = {(increments[[1]]-increments[[2]])/(2*amps[[2]]),
        (increments[[3]]-increments[[4]])/(2*amps[[4]])};
evenError = Max[Abs[Flatten[{increments[[1]]+increments[[2]],increments[[3]]+increments[[4]]}]]];
halfError = Max[Abs[odds[[1]]-odds[[2]]]];
stepError = Max[Abs[Flatten[(loads[#,"32"]-loads[#,"16"])& /@ cases]]];
rateError = Max[Abs[Flatten[Lookup[#["32"],"inward_rate"]& /@ cases]-.25]];
q = 2*First[odds]/times^2;
extrapolated = 2*q[[1]]-q[[2]];
Print["times = ",InputForm[times]];
Print["incremental incoming flux at final time for +A,-A,+A/2,-A/2 = ",
  ExportString[Last /@ increments,"RawJSON","Compact"->True]];
Print["full and half odd incoming slopes per acceleration at final time = ",
  ExportString[Last /@ odds,"RawJSON","Compact"->True]];
Print["max paired even residual, full/half slope discrepancy, timestep change, intake error = ",
  ExportString[{evenError,halfError,stepError,rateError},"RawJSON","Compact"->True]];
Print["2*odd incoming slope/t^2, first two samples extrapolated = ",
  ExportString[{q,extrapolated},"RawJSON","Compact"->True]];
If[!TrueQ[evenError<10^-14 && halfError<10^-12 && stepError<10^-13 && rateError<10^-13 &&
   Abs[extrapolated/frozenSlope-1]<.01],Exit[1]];

(* A known trajectory gives the exterior load L. Treat the remaining
   momentum rate B on that path as an independent unknown in the force
   equation, with no bare kinetic term or force input in the exterior. *)
row = Last[CoefficientArrays[{force+incoming-otherRate},{force,otherRate}]];
Print["inverse force balance nullspace (actuator force, remaining momentum rate) = ",
  InputForm[NullSpace[Normal[row]]]];
If[NullSpace[Normal[row]]=!={{1,1}},Exit[1]];
