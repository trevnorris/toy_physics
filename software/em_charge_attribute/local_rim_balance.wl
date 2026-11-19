(* Minimize exterior volume energy in the decaying dipole radial family.
   Test the result in the full vector PDE, independently of minimization. *)
ClearAll["Global`*"];
emit[key_, v_] := Print[key," = ",ToString[v,InputForm]];
ass=a>0 && m>0 && b>0 && Element[{zn,zt,al,be},Reals];
(* s=r/a. Both trial functions have the imposed traces for all al,be. *)
rad[s_]:=zn/s^3+al(1/s-1/s^3);
tan[s_]:=zt/s^3+be(1/s-1/s^3);
dd[s_]=D[rad[s],s]+2(rad[s]-tan[s])/s;
cc[s_]=(rad[s]-tan[s])/s-D[tan[s],s];
integrand=Expand[a s^2 (b dd[s]^2+2 m cc[s]^2) 2 Pi/3];
ener=Integrate[integrand,{s,1,Infinity},Assumptions->ass];
best=First[Solve[{D[ener,al]==0,D[ener,be]==0},{al,be}]];
estar=Factor[ener/.best];
emit["on-shell volume energy",estar];
(* m curl curl u - b grad div u =0; radial and polar components multiplied
   by a^2, and with their angular factors Cos[theta],Sin[theta] removed. *)
res=Factor /@ ({2 m cc[s]/s-b D[dd[s],s],
   -m D[s cc[s],s]/s+b dd[s]/s}/.best);
emit["two bulk residuals",res];
q=Factor /@ (D[estar,#]& /@ {zn,zt});
net=Factor[Total[q]];
emit["boundary normal and tangential loads",q];
emit["exterior resultant",net];
emit["full-follow stiffness",Factor[net/.{zn->1,zt->1}]];
emit["response Hessian rank",MatrixRank[Table[D[estar,i,j],{i,{zn,zt}},{j,{zn,zt}}]]];
free=First[Solve[q[[2]]==0,zt]];
emit["free tangential trace, its load, its energy",{zt/.free,net/.free,estar/.free}//Simplify];
(* A Cartesian harmonic potential constructs the zero-load field directly. *)
xyz={xx,yy,zz}; rr=Sqrt[xyz.xyz]; pot=-zn a^3 xx/(2 rr^3);
vec=D[pot,#]& /@ xyz;
emit["harmonic-gradient divergence and curl",
    FullSimplify[{Div[vec,xyz],Curl[vec,xyz]},Assumptions->rr>0]];
(* Independent signs on both sides of a diffuse internal cut: scalar
   free-energy stress times its normal on any cut, tested on a generic jet. *)
grad=Array[vv,4]; normal=Array[nn,4];
cut=(en IdentityMatrix[4]-kap Outer[Times,grad,grad]).normal;
emit["summed smooth internal-cut traction",Simplify[cut+(cut/.Thread[normal->-normal])]];
emit["delta reaction on the rim from two signed order fluxes (external hold)",
    Coefficient[(leftSlope leftT-rightSlope rightT) dz,dz]];
(* Kirchhoff balance imposed on the material node, X=0 on a separate rim.
   Solve by constant layer flux and matched traces, not by energy. *)
matching={pp==kl(y-z),ss==ke z,rrim==contact z,pp==ss+rrim};
sol=First[Solve[matching,{pp,ss,rrim,z}]];
emit["incoming, exterior, external hold",Factor /@ ({pp,ss,-rrim}/.sol)];
emit["net local residual",FullSimplify[pp-ss-rrim/.sol]];
emit["decoupled rim and pinned rim: routed loads",{
    Factor /@ (({ss,-rrim}/.sol)/.contact->0),
    Limit[{ss,-rrim}/.sol,contact->Infinity]}];
If[res=!={0,0} || !TrueQ[FullSimplify[pp-ss-rrim/.sol]==0],Exit[1]];
