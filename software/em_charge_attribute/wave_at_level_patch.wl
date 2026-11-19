(* A strip with linear order loss, then an independent meridional Ritz solve. *)
ClearAll["Global`*"];
j=N[BesselJZero[0,1],30];
u[v_]:=BesselJ[0,j (1-v)];
norm=NIntegrate[(1-v) u[v]^2,{v,0,1},WorkingPrecision->25];
load=N[u'[0]^2/(2 j^2 norm)];
Print["strip dimensionless frequency, midpoint amplitude, load times length per energy = ",
 ToString[N[{j,u[1/2],load},12],InputForm]];
If[Abs[load-1]>10^-12 || Abs[u[0]]>10^-20 || Abs[u'[1]]>10^-20,Exit[1]];

gauss[n_]:=gauss[n]=Module[{a=ConstantArray[0.,{n,n}],ev,vec},
 Do[a[[i,i+1]]=a[[i+1,i]]=N[i/Sqrt[4 i^2-1]],{i,n-1}];
 {ev,vec}=Eigensystem[a]; SortBy[Transpose[{ev,2 vec[[All,1]]^2}],First]];
gate[x_]:=x^3(10-15 x+6 x^2);
(* Endpoint hats plus interval-local polynomial bubbles, joined at r=2.
   Their outer endpoint is a free degree of freedom. *)
radial[half_,s_,p_]:=Module[{b=Table[s (1-s) ChebyshevT[k,2s-1],{k,0,p-2}],z},
 z=ConstantArray[0,Length[b]];
 If[half==1,Join[{s,0},b,z],Join[{1-s,s},z,b]]];
collar[p_,a_,beta_]:=Module[{ruleR=gauss[2p+4],ruleT=gauss[3a+6],rr,rt,theta,
  angular,da,rv,dr,bv,db,rad,mu,chi,hR,hT,metric,weight,field,grad1,grad2,
  wr,wt,mat,ks,ms,sol,vec,lam,rin,der,load,mid,eval,
  x,xt,yl,ylp,angfun,angder,rows,vals,weightedGram},
 rr=(1+ruleR[[All,1]])/2; wr=ruleR[[All,2]]/2;
 theta=Pi (1+ruleT[[All,1]])/2; wt=Pi ruleT[[All,2]]/2;
 angfun=Table[Cos[k xt],{k,0,a}]; angder=D[angfun,xt];
 angular=Table[angfun/.xt->v,{v,theta}]; da=Table[angder/.xt->v,{v,theta}];
 rows=Reap[Do[
   yl=radial[h,x,p]; ylp=D[yl,x];
   Do[rad=h+rr[[i]]; rv=yl/.x->rr[[i]]; dr=ylp/.x->rr[[i]];
    bv=If[h==1,gate[rr[[i]]],1.]; db=If[h==1,30 rr[[i]]^2(1-rr[[i]])^2,0.];
    chi=If[h==1,1.,gate[1-rr[[i]]]];
    Do[mu=Cos[theta[[k]]]; hR=0.2(rad-beta mu (bv+rad db));
      hT=0.2 beta bv Sin[theta[[k]]];
      metric=Sqrt[1+hR^2+hT^2];
      weight=2 Pi wr[[i]] wt[[k]] Sin[theta[[k]]] rad^2 chi;
      field=Flatten[Outer[Times,rv,angular[[k]]]];
      grad1=Flatten[Outer[Times,dr,angular[[k]]]];
      grad2=Flatten[Outer[Times,rv,da[[k]]]]/rad;
      (* The inverse metric quadratic form is the sum of these three squares. *)
      Sow[{Sqrt[weight metric] field,Sqrt[weight/metric] grad1,
           Sqrt[weight/metric] grad2,Sqrt[weight/metric](hT grad1-hR grad2)}],
      {k,Length[theta]}],{i,Length[rr]}],{h,2}]][[2,1]];
 mat=Developer`ToPackedArray[N[#]]& /@ Transpose[rows];
 weightedGram[m_]:=Transpose[m].m;
 ms=weightedGram[mat[[1]]]; ks=Total[weightedGram /@ Rest[mat]];
 {vals,sol}=Eigensystem[{ks,ms},-1];
 lam=First[vals]; vec=First[sol]; vec=vec/Sqrt[vec.ms.vec];
 der=(D[radial[1,x,p],x]/.x->0.);
 rin=Table[Flatten[Outer[Times,der,v]].vec,{v,angular}];
 load=2 Pi Total[wt Sin[theta] rin^2]/(2 lam (1+0.2^2));
 mid=Abs[Flatten[Outer[Times,radial[2,0.5,p],ConstantArray[1.,a+1]]].vec];
 eval=Norm[ks.vec-lam ms.vec]/Norm[ks.vec];
 Print["p, angular, beta, omega, rim load per energy, midpoint, residual = ",
   ToString[{p,a,beta,Sqrt[lam],load,mid,eval},InputForm]];
 If[!TrueQ[lam>0 && load>0 && eval<10^-7],Exit[1]];
 Sqrt[lam]];
Do[Module[{freq=Table[collar[level[[1]],level[[2]],b],{b,{2.49,2.5,2.51}}]},
  Print["tilt force per action = ", ToString[-(Last[freq]-First[freq])/0.02,InputForm]]],
 {level,{{7,6},{11,10}}}];
