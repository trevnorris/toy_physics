(* Use the azimuthal vorticity b(r) sin(theta) rather than a displacement
   potential. Laplace s>0 and the exterior domain select the decaying branch. *)
b[r_] := Exp[-k*(r-a)]*(bb/r^2+hh/r);
radial = Factor[Exp[k*(r-a)]*(D[r^2 b'[r],r]-(2+k^2*r^2)*b[r])];
ut[r_] := (b'[r]+b[r]/r)/k^2; (* coefficient of +sin(theta) *)
constants = First[Solve[{radial==0,ut[a]==-d},{bb,hh}]];
bwall = Factor[b[a] /. constants];
(* Work on the exterior, whose outward normal at its inner boundary is -e_r. *)
force = Integrate[(-mod*bwall*Sin[th])*(-Sin[th])*2 Pi a^2 Sin[th],{th,0,Pi}];
transfer = Factor[force/(s^2*d) /. {k->s/c,mod->rho*c^2}];
slow = Limit[transfer,s->0];
response = Factor[transfer/slow /. c->a/tau];
Print["outer radial-equation residual and tangential trace = ",
  InputForm[Simplify[{radial,ut[a]} /. constants]]];
Print["traction integral divided by field-trace acceleration = ",InputForm[transfer]];
Print["slow coefficient and normalized transfer = ",InputForm[{slow,response}]];

(* Check the step response by solving the time-domain initial-value equation
   from the transfer denominator. Initial data follow from the numerator. *)
ts = Cancel[response/s];
g0 = Limit[s*ts,s->Infinity];
g1 = Limit[s*(s*ts-g0),s->Infinity];
sol = DSolveValue[{tau^2 g''[t]+tau g'[t]+g[t]==1,g[0]==g0,g'[0]==g1},g[t],t];
gstep = FullSimplify[ComplexExpand[sol],Assumptions->tau>0 && t>0];
shape = 1-Exp[-t/(2 tau)]*(Cos[Sqrt[3]*t/(2 tau)]-Sin[Sqrt[3]*t/(2 tau)]/Sqrt[3]);
Print["normalized force on acceleration step = ",InputForm[gstep]];
Print["step onset, slope, decay limit = ",InputForm[{g0,g1,
  Limit[shape,t->Infinity,Assumptions->tau>0]}]];
If[!TrueQ[Simplify[(radial /. constants)]===0 && Simplify[(ut[a] /. constants)+d]===0 &&
  Simplify[slow-8 Pi rho*a^3/3]===0 && g0===0 && Simplify[tau*g1]===1 &&
  FullSimplify[gstep-shape,Assumptions->t>0 && tau>0]===0],Exit[1]];

(* Compare two different terminations of the same shear equation. A direct
   tangential tie sets d=X, while an unforced traction-free rim with zero
   initial field has b(a)=d=0. This does not impose a normal fluid condition. *)
free = Solve[{radial==0,b[a]==0},{bb,hh}];
freeTrace = Simplify[ut[a] /. First[free]];
Print["unforced traction-free tangential displacement = ",InputForm[freeTrace]];
If[freeTrace=!=0,Exit[1]];
file = "docs/_measurements/mouth_prescribed_acceleration.txt";
If[FileHash[file,"SHA256","HexString"]=!=
  "db7ae4e93ffb1ddf105d2afd02194f29009464d9a642721a18e3edd3aeeeae73",Exit[1]];
lines = Select[StringSplit[ReadString[file],"\n"],StringStartsQ[#,"{\"acceleration\":"]&];
cases = (ImportString[#,"RawJSON"]&) /@ lines;
amps = Lookup[cases,"acceleration"];
If[amps!={0,.02,-.02,.01,-.01} || !And@@(KeyExistsQ[#,"32"]& /@ cases),Exit[1]];
loads = (Last[#["32"]]["feed_load"]&) /@ cases;
Print["four drives and known changes of feeding-side flux = ",
 ExportString[{Rest[amps],Rest[loads]-First[loads]},"RawJSON","Compact"->True]];
Print["conditional tied-minus-free shear force divided by slow coefficient and step shape = ",
 InputForm[Rest[amps]]];
