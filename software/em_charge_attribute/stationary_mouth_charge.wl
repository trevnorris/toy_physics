(* Calculate the angular current directly on a moving, distorted collar. *)
ClearAll["Global`*"];
check[name_,value_,expected_]:=Module[{v=FullSimplify[value]},
  Print[name," = ",ToString[v,InputForm]];
  If[!TrueQ[FullSimplify[v==expected]],Exit[1]]
];
(* A small rectangle in turn/order state space measures path dependence
   of chi dB + f dchi, and fixes how f must depend on B. *)
loop = c db + f[b+db,c] dc - (c+dc) db - f[b,c] dc;
check["infinitesimal state-space loop per oriented area",
  D[loop,db,dc]/.{db->0,dc->0},Derivative[1,0][f][b,c]-1];
primitive=c b+h[c];
check["general endpoint charge", (primitive/.c->0)-(primitive/.{c->1,b->0}),h[0]-h[1]];
check["vertical-layer local contribution",Integrate[s-kk,{c,1,0}],kk-s];
(* Use a smooth spatial transition rather than order as the coordinate. *)
drop=3u^2-2u^3;
chi=1-drop;
turn=s(drop+alpha drop(1-drop));
Do[
 raw=Integrate[chi D[turn,u],{u,0,1}]/.{alpha->drive,s->1};
 total=Integrate[chi D[turn,u]+(turn-sector)D[chi,u],{u,0,1}]/.{alpha->drive,s->1};
 Print["shape, ordered, neutral total, unit total = ",
   ToString[{drive,raw,total/.sector->0,total/.sector->1},InputForm]];
 check["sum at general sector",total,sector],
 {drive,{-1/2,0,1/2}}];

(* Compute the pullback from a unit director, including a time-dependent
   twist. x=cos(theta); use angular orientation dphi wedge dx. *)
er={Sqrt[1-x^2]Cos[phi+tau u],Sqrt[1-x^2]Sin[phi+tau u],k x};
director=er/Sqrt[1+(k^2-1)x^2];
area[i_,j_]:=FullSimplify[director.Cross[D[director,i],D[director,j]],
  -1<x<1 && k>0 && Element[{phi,tau,u},Reals]];
angular=area[phi,x];
check["normal solid-angle density",angular,k/(1+(k^2-1)x^2)^(3/2)];
check["time twist contributes to pulled-back solid angle",area[tau,x],u angular];
check["closed angular two-form on the moving (u,time,x) chart",
  D[area[tau,x],u]+D[area[x,u],tau]+D[area[u,tau],x],0];
(* For -dchi wedge Omega with chi=chi(u), the spatial density is
   proportional to angular and the transverse time flux to area[tau,x].
   The temporal twist rotates the density without changing its integral. *)
Do[
 signed = sg Integrate[angular/.k->stretch,{x,-1,1}]/2;
 check["orientation and stretch "<>ToString[{sg,stretch},InputForm],signed,sg],
 {sg,{-1,1}},{stretch,{1/2,1,2}}];
check["surface density double-radius ratio",
  (1/(4Pi (2a)^2))/(1/(4Pi a^2)),1/4];
modeOverlap=Integrate[-D[chi,u],{u,0,1}];
check["trial layer response to constant potential",
  -D[-probe strength s modeOverlap,probe],strength s];

(* Same throughput at both cuts. Count occupants by elapsed time of flight
   through a stationary profile rather than a constituent's charge loss. *)
travel=Integrate[(1+eps Sin[Pi y/len]^2)/vref,{y,0,len},Assumptions->len>0&&vref>0];
check["trapped versus transit control, excess carrier charge",q flux (travel-len/vref),
  q flux eps len/(2vref)];
check["same charge in and out",q flux-q flux,0];
check["dropping the carrier charge at the front, no other outlet",q flux-0,q flux];
