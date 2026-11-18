(* The surface density of virtual work is varied at fixed solid angle.
   Test the centroid and total load for arbitrary finite angular variations. *)
rpert = a + e*(h0+h1*u+h2*(3*u^2-1)/2);
vol = 2 Pi Integrate[rpert^3/3,{u,-1,1}];
mom = 2 Pi Integrate[u*rpert^4/4,{u,-1,1}];
shift = Simplify[D[mom/vol,e] /. e->0, Assumptions->a>0];
traction = 3 f u/(4 Pi a^2);
work = Integrate[2 Pi a^2 traction D[rpert,e],{u,-1,1}];
Print["center variation and rim work = ",InputForm[{shift,work}]];
If[Simplify[work-f*shift]=!=0 || shift=!=h1,Exit[1]];

(* With a locked field, impose y=cX and calculate power directly. With a
   compliant field use the retarded Green function for a unit-kinetic mode.
   Both terms can be added to the same otherwise fixed feed equations. *)
path[t_] := drive*t^2/2;
lockedEnergy = D[c*path[t],t]^2/2;
lockedF = Simplify[D[lockedEnergy,t]/D[path[t],t]];
response = Integrate[Sin[om*(t-u)]/om * om^2*c*path[u],{u,0,t},
   Assumptions->om>0 && t>0];
springF = Simplify[c*om^2*(c*path[t]-response)];
energy = (D[response,t]^2+om^2*(response-c*path[t])^2)/2;
powerResidual = FullSimplify[D[energy,t]-springF*D[path[t],t]];
delta = FullSimplify[lockedF-springF];
initial = Simplify[{response,D[response,t],lockedF,springF} /. t->0];
Print["locked force from power, compliant force from field propagation = ",
  InputForm[{lockedF,springF}]];
Print["initial field, field speed and two extra actuator forces = ",InputForm[initial]];
Print["difference in required force and power residual = ",InputForm[{delta,powerResidual}]];
If[!TrueQ[initial=={0,0,c^2*drive,0} && powerResidual==0 &&
    Simplify[delta-c^2*drive*Cos[om*t]]==0],Exit[1]];

(* Use the accepted refinements of the four exterior experiments as data.
   Applying identical X in either interface example leaves their fluxes fixed.
   Only the known additive term in the momentum budget can be recovered. *)
manifest = "docs/_measurements/mouth_prescribed_acceleration.txt";
If[FileHash[manifest,"SHA256","HexString"] =!=
  "db7ae4e93ffb1ddf105d2afd02194f29009464d9a642721a18e3edd3aeeeae73",Exit[1]];
manifestLines = StringSplit[ReadString[manifest],"\n"];
readRun[tag_,expected_] := Module[{dir,file,record,out,last},
  dir = "_scratch/s11c/mouth_prescribed_accel_"<>tag<>"02/";
  Do[file=dir<>leaf;
     If[!MemberQ[manifestLines,FileHash[file,"SHA256","HexString"]<>"  "<>file],Exit[1]],
     {leaf,{"outcome.json","stdout"}}];
  out=Import[dir<>"outcome.json","RawJSON"];
  If[!TrueQ[out["limitsVerified"] && out["exitCode"]==0 &&
     out["childOutcome"]["exitCode"]==0],Exit[1]];
  record=ImportString[ReadString[dir<>"stdout"],"RawJSON"];
  last=Last[record["32"]];
  If[!TrueQ[record["acceleration"]==expected && last["t"]==.001 &&
     Abs[last["inward_rate"]-.25]<10^-13],Exit[1]];
  last["feed_load"]];
amps={1/50,-1/50,1/100,-1/100};
control=readRun["zero",0];
known=control-MapThread[readRun,{{"positive","negative","positive_half","negative_half"},amps}];
Print["known contribution minus change in incoming flux, in drive order = ",
  ExportString[known,"RawJSON","Compact"->True]];
Print["four force differences per c^2 at onset = ",InputForm[amps*(delta/(drive*c^2) /. t->0)]];
