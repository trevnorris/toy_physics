(* First variation of a translated ball in material coordinates. Integrate
   across disks of constant x, then differentiate its incident-wave phase. *)
ClearAll["Global`*"];
emit[key_, v_] := Print[key, " = ", ToString[v,InputForm]];
ass = r>0 && k>0 && Element[{phase,shift,dr,dm,om,amp},Reals];
(* This is the interior minus exterior cycle-averaged Lagrangian density,
   with dr,dm the outward jumps. *)
lag[x_] := amp^2 (-dr om^2 Cos[k (x+shift)+phase]^2+
                    dm k^2 Sin[k (x+shift)+phase]^2)/4;
force = FullSimplify[Integrate[Pi (r^2-x^2) (D[lag[x],shift]/.shift->0),
    {x,-r,r},Assumptions->ass],Assumptions->ass];
emit["force per depth from displacement of the coefficient ball",force];
geom=FullSimplify[force/(Pi amp^2 (dr om^2+dm k^2) Sin[2 phase]/2),
    Assumptions->ass];
emit["geometry factor alone",geom];
emit["reflected and centered fringe plus matched-impedance controls",
  FullSimplify[{force+(force/.phase->-phase),force/.phase->0,
      force/.{dr->-rho eps,dm->mu eps,om->k Sqrt[mu/rho]}},
    Assumptions->ass && rho>0 && mu>0]];
(* Verify the same load's far-field limit on a flare by Gaussian quadrature
   through NIntegrate, with the shape measured from the mouth to either end. *)
radius[s_]:=1+(3s^2-2s^3)/2;
f[s_]=geom/.{k->1/2,r->radius[s]};
zeroth=NIntegrate[Evaluate[f[s]],{s,0,1},WorkingPrecision->35];
cosine=NIntegrate[Evaluate[f[s] Cos[Pi s]],{s,0,1},WorkingPrecision->35];
emit["midplane, end, end/midplane and cosine/common generalized load",
    N[{f[0],f[1],f[1]/f[0],cosine/zeroth},15]];
emit["long wave exit/midplane",Limit[(geom/.r->3/2)/(geom/.r->1),k->0]];
emit["constant-section cosine projection",Integrate[(geom/.r->1) Cos[Pi s],{s,0,1}]];

(* Independent finite-contrast test. Sum internal reflections by matching the
   reflected wave at a single interface, its reverse and the two propagation
   phases of a round trip. Amplitudes are referenced to the two slab faces. *)
single=(1-zi)/(1+zi);
refl=single (1-Exp[2 I d])/(1-single^2 Exp[2 I d]);
trans=(1-single^2) Exp[I d]/(1-single^2 Exp[2 I d]);
net=FullSimplify[-Im[ComplexExpand[refl Conjugate[trans]]]/2,
    Assumptions->zi>0 && Element[d,Reals]];
unitarity=FullSimplify[ComplexExpand[Abs[refl]^2+Abs[trans]^2-1],
    Assumptions->zi>0 && Element[d,Reals]];
emit["exact slab power residual",unitarity];
(* A=k=omega=rho_out=mu_out=1. Peak signed sin(2 phase)=1.
   At zero contrast the thickness phase is independent of exponent. *)
linear=FullSimplify[(D[net,zi] (-(m+1)/2)+D[net,d] ((1-m) d/2))/.zi->1,
    Assumptions->Element[{d,m},Reals]];
emit["planar linear coefficient at arbitrary exterior thickness phase",linear];
slab[e_,ex_] := net/.{zi->(1-e)^((1+ex)/2),d->Pi/2 (1-e)^((ex-1)/2)};
Do[emit["finite slab force per depletion for m="<>ToString[ex],
    N[Table[slab[e,ex]/e,{e,{1/50,1/100,1/1000}}],16]],{ex,{0,1,-1}}];
If[!TrueQ[unitarity==0] || !VectorQ[{zeroth,cosine},NumericQ],Exit[1]];
