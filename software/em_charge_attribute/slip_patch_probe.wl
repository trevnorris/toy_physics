(* Exact surface kinematics and an independent standing-mode readout on saved graphs. *)
ClearAll["Global`*"];
write[tag_, val_] := Print[tag, " = ", ToString[val,InputForm]];
height = (r^2-1)/10 - r u/2;
rs = 5/2;
hr = D[height,r] /. r->rs;
ht = -Sqrt[1-u^2]/rs (D[height,u] /. r->rs);
w = FullSimplify[Sqrt[1+hr^2+ht^2],-1<=u<=1];
speed2 = 1+hr^2; (* unit coordinate radial speed *)
mean = (2/w+1/w^3)/15; (* Hessian of height is the Cartesian identity / 5 *)
pair[expr_] := FullSimplify[{expr /. u->1,expr /. u->0}];
write["slip/equator: W, surface n/n0, speed^2/v0^2, RH, height/R (common offset removed)",
      pair /@ {w,1/w,speed2,mean,height /. r->rs}];
(* Conditional carrier check with a common normal thickness and fixed velocity.
   Only the local EOS and convective kinetic energy are included. *)
chemical = D[mm n vv/2 + kap n^5/4,n] /. {n->n0/w,vv->v0^2 speed2};
contrast = Factor[Subtract@@pair[chemical]];
write["slip-equator EOS and kinetic insertion difference",contrast];
write["slope of conditional chemical potential in 1-cos(theta) at pole",
      FullSimplify[D[chemical /. u->1-y,y] /. y->0]];
(* In the actual collar the wave and geometric terms are independent of probe number. *)
held[dn_,site_] := area+bend+j omega+dn q pot[site];
write["direct first variation of held geometry/wave plus electric probe, slip minus equator",
      FullSimplify[D[held[dn,slip]-held[dn,equator],dn]]];

(* Chebyshev polynomials on each radial half, cosine series in polar angle.
   Compute mass and stiffness using the inverse induced metric in polar coordinates. *)
gauss[nn_] := gauss[nn] = Module[{mat=ConstantArray[0.,{nn,nn}],evec,eval},
  Do[mat[[i,i+1]]=mat[[i+1,i]]=N[i/Sqrt[4i^2-1]],{i,nn-1}];
  {eval,evec}=Eigensystem[mat]; SortBy[Transpose[{eval,2 evec[[All,1]]^2}],First]];
step[x_] := x^3 (10-15x+6x^2);
radial[h_,x_,p_] := Module[{b=Table[x(1-x) ChebyshevT[i,2x-1],{i,0,p-2}],z},
  z=ConstantArray[0,Length[b]]; If[h==1,Join[{x,0},b,z],Join[{1-x,x},z,b]]];
tensor[a_,b_] := Flatten[Outer[Times,a,b]];
gram[a_] := Transpose[a].a;
probe[path_,label_] := Module[{state,bb,coeff,nr,na,pow,p=10,lmax=9,rq=gauss[26],
    tq=gauss[32],hfun,gradient,rr,ang,dr,dt,wt,chi,basis,bp,ap,dap,rows,
    mblock,kr,kt,kcross,mass,stiff,eigen,vec,omega,values,part,energy,
    rv,rg,av,ag,pr,pt,volume,weight,loc,ur,ptval,sens,res,scale},
  state=SelectFirst[path["states"],TrueQ[#["beta"]==2.5]&];
  bb=state["beta"]; coeff=state["coeff"]; {nr,na}=path["shape"];
  pow=If[path["free_angle"],1,2]; ap=Table[Cos[k angle],{k,0,lmax}]; dap=D[ap,angle];
  {mblock,kr,kt,kcross}=Developer`ToPackedArray /@ Transpose[Reap[Do[
    hfun=(r^2-1)/10 - bb r u If[h==1,step[r-1],1]/5;
    If[h==1,hfun+=(r-1)^pow (2-r)^2 Sum[coeff[[1+i*nr+k]]
         LegendreP[i,u] LegendreP[k,2r-3],{i,0,na},{k,0,nr-1}]];
    gradient=Compile[{{r,_Real},{u,_Real}},Evaluate[HornerForm[#,{r,u}]& /@ {D[hfun,r],D[hfun,u]}]];
    basis=radial[h,x,p]; bp=D[basis,x];
    Do[rr=h+(1+s[[1]])/2; rv=basis /. x->rr-h; rg=bp /. x->rr-h;
      chi=If[h==1,1.,step[3-rr]];
      Do[ang=Pi(1+t[[1]])/2; av=ap /. angle->ang; ag=dap /. angle->ang;
        {pr,pt}=gradient[rr,Cos[ang]]; pt=-Sin[ang]pt/rr;
        volume=Sqrt[1+pr^2+pt^2]; weight=Pi^2 s[[2]] t[[2]] rr^2 Sin[ang] chi/2;
        dr=tensor[rg,av]; dt=tensor[rv,ag]/rr; scale=Sqrt[weight/volume];
        Sow[{Sqrt[weight volume] tensor[rv,av],scale dr,scale dt,scale (pt dr-pr dt)}],
        {t,tq}],{s,rq}],{h,2}]][[2,1]]];
  mass=gram[mblock]; stiff=gram[kr]+gram[kt]+gram[kcross];
  {eigen,vec}=Eigensystem[{stiff,mass},-1]; omega=Sqrt[First[eigen]];
  vec=First[vec]; vec=vec/Sqrt[vec.mass.vec];
  res=Norm[stiff.vec-omega^2 mass.vec]/Norm[stiff.vec];
  values=Table[av=ap /. angle->pangle; loc=tensor[radial[2,.5,p],av].vec;
    ur=tensor[D[radial[2,x,p],x] /. x->.5,av].vec;
    {loc^2,2loc ur,-path["action"]omega loc^2/2},{pangle,{0.,N[Pi/2]}}];
  write[label<>" omega, residual, each site {U^2, radial slope of U^2, mass sensitivity}",
        {omega,res,values}];
  If[!TrueQ[res<10^-7 && Abs[omega-state["parts"][[3]]]<0.002],Exit[1]]];
paths=Import["docs/_measurements/relaxed_neck_paths.json","RawJSON"];
Do[probe[paths[label],label],{label,{"free","clamped"}}];
