(* Calculate the wall stress using normalized eigenfunctions on the moving side.
   The two modes use the first radial root in their own angular sector. *)
ClearAll["Global`*"];
say[key_,value_] := Print[key," = ",ToString[value,InputForm]];
archive = Import["software/stage1_solver/reports/pathA_26_derrick_results.yaml","Text"];
value[pattern_] := ToExpression[First[StringCases[archive,RegularExpression[pattern]->"$1"]]];
rad0 = value["phase_A:\\s+q_star:\\s+- ([0-9.]+)"];
len0 = value["phase_A:\\s+q_star:\\s+- [0-9.]+\\s+- ([0-9.]+)"];
param[name_] := value["(?m)^  "<>name<>": ([0-9.]+)$"];
speed = param["c_w"]; gap = param["mu_in"]; cutoff=param["mu_out"];
root[l_,n_] := BesselJZero[l+1/2,n];
freq[l_,n_,j_,r_,len_] := Sqrt[gap^2+speed^2 ((root[l,n]/r)^2+((j+1/2)Pi/len)^2)];
freqs=N[{freq[0,1,0,rad0,len0],freq[1,1,0,rad0,len0],freq[0,1,1,rad0,len0],
         freq[0,2,0,rad0,len0]}];
say["reference frequencies and trapping threshold",{freqs,cutoff}];
zeta=N[root[1,1],28];
aspect=ll /. FindRoot[zeta^2+(Pi/2/ll)^2==Pi^2+(3 Pi/2/ll)^2,{ll,1.4}];
rcrit = rr /. FindRoot[gap^2+speed^2 (zeta^2+(Pi/(2 aspect))^2)/rr^2==cutoff^2,{rr,2.}];
say["crossing aspect, radius and length at loss of binding",{aspect,rcrit,aspect*rcrit}];
axial[j_,v_] := Sqrt[2] Sin[(j+1/2) Pi v];
az0 = 1/Sqrt[4 Pi]; az1=Sqrt[3/(4 Pi)] t;
ang = Integrate[2 Pi az0 az1 t,{t,-1,1}];
azTerms = Integrate[2 Pi t {az0^2, az1^2},{t,-1,1}];
ax = Integrate[axial[0,v] axial[1,v] Cos[Pi v],{v,0,1}];
common = Integrate[axial[0,v] axial[1,v],{v,0,1}];
(* Normalize each radial eigenmode by volume integration, then read the
   two surface derivatives. Scaling lengths by r restores dimensions. *)
(* Use the derivative recurrence at integer order independently of root finding. *)
slopes=Table[Module[{k=N[root[l,1],30],nr},
 nr=Sqrt[NIntegrate[v^2 SphericalBesselJ[l,k v]^2,{v,0,1},WorkingPrecision->26]];
 (l SphericalBesselJ[l,k]-k SphericalBesselJ[l+1,k])/nr],{l,0,1}];
g[r_] := -speed^2 ang ax (Times @@ slopes)/r^3;
gmatrix={{0,-1},{-1,0}};
state={1,Exp[I phase]}/Sqrt[2];
signedProjection=ComplexExpand[Conjugate[state].gmatrix.state]//Simplify;
say["angular diagonals, bow and translation overlaps, phase projection",
    {azTerms,ax,common,signedProjection}];
(* At an exact crossing, slowly following the lower excited eigenfrequency
   gives equal weights on the two straight-tube occupations at b = 0+. *)
readings=Table[Module[{len=aspect r,om,wL},
 om=freq[1,1,0,r,len];
 wL= D[(freq[0,1,1,r,t]+freq[1,1,0,r,t])/2,t] /. t->len;
 {om,Abs[g[r]]/(2 om),-wL}],{r,{2,4}}];
say["tuned frequencies, bow load/J and anchored depth load/J (ell=1)",readings];
(* At fixed aspect and relative bow eta, use stationarity against a common
   dilation to eliminate the wave action from the generalized lateral force. *)
wave = Sqrt[mass^2+c^2 lam[eta]/r^2];
geom = surf r^3+vol r^4+curv r;
eliminated=Factor[-iW D[wave,eta]/r /. First[Solve[r D[geom,r]==-iW r D[wave,r],iW]]];
say["wave load along bow at size equilibrium",eliminated];
audit=N[Flatten[{freqs,aspect,aspect*rcrit,readings}]];
Print["comparison ",StringRiffle[ToString[#,InputForm]& /@ audit," "]];
If[!TrueQ[Max[Abs[slopes+Sqrt[2]{Pi,zeta}]]<10^-18 &&
  azTerms=={0,0} && common==0 && ax==1/2 && signedProjection==-Cos[phase] &&
  Min[Flatten[readings]]>0],Exit[1]];
