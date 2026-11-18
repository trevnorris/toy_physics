(* Independent signed virtual-work and onset checks on the archived feeding face. *)
ClearAll["Global`*"];
file="_scratch/s11c/mouth_rim_push_sympy01/resting_face.f64";
If[FileHash[file,"SHA256","HexString"]=!="e76812f953ed1af4fd216a5d7fb9d70da6801536dfc786a48fcea165c8ca8053",Exit[1]];
{re,im,j0}=ArrayReshape[BinaryReadList[file,"Real64",ByteOrdering->-1],{3,32,32,32}];
dx=1/8; dz=1/64; axis=N[Range[-31,31,2]/16]; den=re^2+im^2; u=Sqrt[den];
h[s_]:=Table[Max[0.,1-(x-s)^2-y^2-z^2],{x,axis},{y,axis},{z,axis}];
integral[arr_]:=dx^3 Total[Flatten[arr]];
energy[s_]:=integral[h[s]^3 (u-1)^2]/2;
weight=h[0]; g=weight^3;
g1=6 MapThread[Times,{axis,weight^2}];
g2=24 MapThread[Times,{axis^2,weight}]-6 weight^2;
slope=-integral[g2 (u-1)^2]/2;
fd=Table[-(energy[e]+energy[-e]-2 energy[0])/e^2,{e,{0.0001,0.00005}}];
Print["energy curvature, analytic force per displacement = ",InputForm[fd]," ",InputForm[slope]];
If[!TrueQ[Max[Abs[(fd-slope)/slope]]<10^-5],Exit[1]];

(* At onset d_X (normal slope) changes the density through the imposed phase
   flux. The finite-grid ghost plane contributes i d_X(normal)/dz. *)
fourierDx=N[Table[If[i==j,0,Pi/(32 dx) (-1)^(i-j) Cot[Pi (i-j)/32]],{i,32},{j,32}]];
gradn=ArrayReshape[fourierDx.ArrayReshape[den,{32,1024}],{32,32,32}];
dr=g1+g gradn/den;
dj=integral[j0] (dr-g integral[dr]/integral[g])/integral[g];
fieldCubic=integral[g1 (u-1) dj/u]/(6 dz);
Print["field-dependent virtual force per A*t^3 at onset = ",InputForm[fieldCubic]];
If[!TrueQ[Abs[integral[dj]]<10^-13 && fieldCubic!=0],Exit[1]];

(* Every equation for a prescribed-motion, incrementally source-free receiver
   is unchanged up to a common multiplier. Currents and momenta carry it. *)
local={mw aw-tw qw2-k (y-q), my ay-ty yy2+k (y-q)};
end={-tw qw1-z11 vw-z12 vy,-ty yy1-z12 vw-z22 vy};
observers={mw vw,my vy,-tw qw1,-ty yy1};
scale=Table[c->lambda c,{c,{mw,my,tw,ty,k,z11,z12,z22}}];
res=Simplify[(Join[local,end,observers]/.scale)-lambda Join[local,end,observers]];
Print["common load scale residuals = ",InputForm[res]];
If[res=!=ConstantArray[0,8],Exit[1]];
Exit[0];
