(* Independent consistent-mass finite elements for the longer positive benchmark.
   Eliminate the driven wall coordinate. All other wall and packet nodes are free. *)
ClearAll["Global`*"]; $HistoryLength=0;
cut=1.85/4; depth=5 cut/4; nc=64; ne=5 nc/4; nn=ne+1; dz=depth/ne;
nodes=Range[0,ne] dz; mid=Most[nodes]+dz/2; ix=nc+1;
f2=Cos[Pi mid/(2*1.85)]^2;
assemble[weight_,local_]:=Module[{diag,off},
  diag=PadRight[weight local[[1,1]],nn]+PadLeft[weight local[[2,2]],nn];
  off=weight local[[1,2]];
  SparseArray[{Band[{1,1}]->diag,Band[{2,1}]->off,Band[{1,2}]->off},{nn,nn}]];
ones=ConstantArray[1.,ne]; zero=SparseArray[{},{nn,nn}];
massForm=dz {{2,1},{1,2}}/6; gradForm={{1,-1},{-1,1}}/dz;
mw=assemble[ones,massForm]; my=assemble[f2,massForm];
kw=assemble[ones,gradForm]; ky=assemble[0.5625 f2,gradForm];
lock=assemble[2 f2,massForm];
mm=SparseArray[ArrayFlatten[{{mw,zero},{zero,my}}]];
kk=SparseArray[ArrayFlatten[{{kw+lock,-lock},{-lock,ky+lock}}]];
rot=RotationMatrix[0.3];
zz=rot.DiagonalMatrix[{1.08/1.2,0.4/0.8}].Transpose[rot];
cc=SparseArray[Flatten[Table[{i nn,j nn}->zz[[i,j]],{i,2},{j,2}]],{2 nn,2 nn}];
free=Range[2,2 nn]; m=mm[[free,free]]; k=kk[[free,free]]; damp=cc[[free,free]];
amp=0.02; dt=1./1024; drive[t_]:=amp t^2/2;
force[t_]:=-amp mm[[free,1]]-amp t cc[[free,1]]-drive[t] kk[[free,1]];
solve=LinearSolve[m+dt damp/2+dt^2 k/4];
u=ConstantArray[0.,Length[free]]; v=u;
acc=LinearSolve[m,Normal[force[0]]]; impulse={0.,0.}; traction=impulse;
tCut={1.,0.5625 Cos[Pi cut/3.7]^2};
Do[
  time=j dt; up=u+dt v+dt^2 acc/4; vp=v+dt acc/2;
  acc=solve[Normal[force[time]-damp.vp-k.up]];
  u=up+dt^2 acc/4; v=vp+dt acc/2;
  pos=Partition[Prepend[u,drive[time]],nn];
  current=-tCut (pos[[All,ix+1]]-pos[[All,ix-1]])/(2 dz);
  impulse+=dt (traction+current)/2; traction=current,
  {j,1024}];
vel=Partition[Prepend[v,amp time],nn];
(* Element midpoint integration on the truncated interval, with independent
   spatial quadrature for momentum. Flux uses a centered gradient at the cut. *)
p=MapThread[dz Total[Take[#1 (Most[#2]+Rest[#2])/2,nc]]&,
            {{ones,f2},vel}];
computed={p,traction,impulse};
Print["finite element {stored, cut flux, cut impulse} = ",InputForm[computed]];
log=FileNameJoin[{DirectoryName[$InputFileName],"..","..","..","_scratch","s11c","signed_mouth_collar01","stdout"}];
rows=ImportString[#,"RawJSON"]& /@ Import[log,"Lines"];
sample=SelectFirst[rows,Lookup[#,"case",""]=="benchmark_continued_128_positive"&];
observed=Lookup[sample,{"stored","cut_flux","cut_impulse"}];
rel=Abs[(computed-observed)/observed];
Print["signed benchmark readings = ",InputForm[observed]," fractional differences = ",InputForm[rel]];
If[!TrueQ[Dimensions[rel]=={3,2} && Max[rel]<0.005],Exit[1]];
Exit[0];
