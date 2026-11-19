(* Independent energy readout on the supplied physical graphs. *)
ClearAll["Global`*"];
paths=Import["docs/_measurements/relaxed_neck_paths.json","RawJSON"];
gauss[n_]:=gauss[n]=Module[{a=ConstantArray[0.,{n,n}],val,vec},
 Do[a[[i,i+1]]=a[[i+1,i]]=N[i/Sqrt[4i^2-1]],{i,n-1}];
 {val,vec}=Eigensystem[a];SortBy[Transpose[{val,2vec[[All,1]]^2}],First]];
step[x_]:=x^3(10-15x+6x^2);
radial[h_,s_,p_]:=Module[{b=Table[s(1-s)ChebyshevT[k,2s-1],{k,0,p-2}],z},
 z=ConstantArray[0,Length[b]];If[h==1,Join[{s,0},b,z],Join[{1-s,s},z,b]]];
(* A Cartesian graph has three-mean-curvature div(grad h / W).
   Expand the divergence in r, cos(theta), independently of the Hessian frame. *)
kernel[h_] := Module[{fields},
 fields=HornerForm[Expand[#],{r,u}]& /@
    {D[h,r],D[h,u],D[h,{r,2}],D[h,r,u],D[h,{u,2}]};
 Compile[{{r,_Real},{u,_Real}},Evaluate[fields]]];
waveRows[h_,order_,p_,ell_]:=waveRows[h,order,p,ell]=Module[{rg=gauss[order],tg=gauss[order+4],
 b,db,y,dy,out,rad,theta,ar,fr,fdr,fa,fda,rval,mu,measure,v,vr,vt},
 b=radial[h,x,p];db=D[b,x]; y=Table[Cos[i angle],{i,0,ell}];dy=D[y,angle];
 out=Flatten[Table[rad=h+(1+qr[[1]])/2;theta=Pi(1+qt[[1]])/2;mu=Cos[theta];
   fr=b/.x->rad-h;fdr=db/.x->rad-h;fa=y/.angle->theta;fda=dy/.angle->theta;
   measure=Pi^2 qr[[2]]qt[[2]] rad^2 Sin[theta]/2;
   {rad,mu,measure,If[h==1,1.,step[3-rad]],Flatten[Outer[Times,fr,fa]],
    Flatten[Outer[Times,fdr,fa]],Flatten[Outer[Times,fr,fda]]/rad},
   {qr,rg},{qt,tg}],1];Developer`ToPackedArray /@ Transpose[out]];
energy[path_,state_] := Module[{beta=state["beta"],coef=state["coeff"],nr,na,pow,
  height,s=0.,bend=0.,km=0.,mm=0.,spec,shape,tab,rd,mu,weight,chi,v,vr,vt,
  jets,pr,hm,hrr,hrm,hmm,pt,w,mean,cv,scale,weighted,mode,lam},
 {nr,na}=path["shape"];pow=If[path["free_angle"],1,2];
 Do[
   height=(r^2-1)/10-beta*r*u*If[h==1,step[r-1],1]/5;
   If[h==1,height+= (r-1)^pow(2-r)^2 Sum[
       coef[[1+l*nr+m]] LegendreP[l,u] LegendreP[m,2r-3],{l,0,na},{m,0,nr-1}]];
   shape=kernel[height];{rd,mu,weight,chi,v,vr,vt}=waveRows[h,28,10,9];
   {pr,hm,hrr,hrm,hmm}=Transpose[MapThread[shape,{rd,mu}]];
   pt=-Sqrt[1-mu^2]hm/rd;w=Sqrt[1+pr^2+pt^2];
   mean=((2rd pr+rd^2 hrr+(1-mu^2)hmm-2mu hm)/w-
      (rd^2 pr*(pr*hrr+(1-mu^2)hm*hrm/rd^2-(1-mu^2)hm^2/rd^3)+
      (1-mu^2)hm*(pr*hrm+(1-mu^2)hm*hmm/rd^2-mu*hm^2/rd^2))/w^3)/(3rd^2);
   s+=Total[weight*w];bend+=Total[weight*w*mean^2];
   cv=pt*vr-pr*vt;scale=weight*chi/w;
   weighted[a_,q_]:=Transpose[a].(q*a);
   km+=weighted[vr,scale]+weighted[vt,scale]+weighted[cv,scale];
   mm+=weighted[v,weight*chi*w],
   {h,2}];
 lam=First[Eigenvalues[{km,mm},-1]];
 {s,bend,Sqrt[lam],s+path["stiffness"]bend+path["action"]Sqrt[lam]}];
Do[Module[{path=paths[label],points,base,z,res},
 points=Select[path["states"],MemberQ[{0.,0.025,1.1,2.,2.5,3.},N[#["beta"]]]&];
 If[Length[points]!=6,Exit[1]];
 Do[res=energy[path,z];
    If[N[z["beta"]]==0.,base=Last[res]];
    Print[label," beta=",z["beta"]," A,B,omega,E,dE = ",ToString[Append[res,Last[res]-base],InputForm]];
    If[!VectorQ[res,NumericQ] || !TrueQ[res[[3]]>0],Exit[1]],
    {z,points}]],{label,{"free","clamped"}}];
