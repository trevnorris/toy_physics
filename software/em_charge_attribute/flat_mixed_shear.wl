(* Spectral response on the full free-faced interval. Width H is hh.
   Work in one tangential polarization orthogonal to in-brane wave number k. *)
ClearAll["Global`*"];
say[tag_, val_] := Print[tag, " = ", ToString[val,InputForm]];
f0[z_] := 1/Sqrt[hh];
fo[j_,z_] := Sqrt[2/hh] Sin[(2 j+1) Pi z/hh];
fe[j_,z_] := Sqrt[2/hh] Cos[2 j Pi z/hh];
diffAt[f_] := f[b]-f[-b];
sumAt[f_] := f[b]+f[-b];
residues = Simplify[{diffAt[f0]^2,sumAt[f0]^2}/eta,hh>0];
bcs = FullSimplify[{
   D[fo[j,z],z]/.z->hh/2, D[fo[j,z],z]/.z->-hh/2,
   D[fe[j,z],z]/.z->hh/2, D[fe[j,z],z]/.z->-hh/2},
   Element[j,Integers] && hh>0];
firstMass=Pi Sqrt[mu/eta]/hh;
say["derivatives at both faces, difference/sum residues, first inverse range",
 {bcs,residues,firstMass}];
(* General mode normalization: optical inertia uses rho, static stiffness eta. *)
dispResidual = Simplify[(-mu D[fo[j,z],{z,2}]+eta k^2 fo[j,z])-
   (eta k^2+mu ((2 j+1)Pi/hh)^2) fo[j,z]];
say["mixed shear overlap with fundamental",D[f0[z],z]];
(* At b=H/4 the nonzero differences have squared weight 4/H. The nonzero
   sums have 8/H on modes numbered 4,8,... plus 4/H on mode zero. *)
odd = 4 Sum[1/(k^2+(2 j+1)^2 Pi^2), {j,0,Infinity},
             Assumptions->k>0];
even = 4/k^2+8 Sum[1/(k^2+(4 j Pi)^2), {j,1,Infinity},
                    Assumptions->k>0];
odds = FullSimplify[odd,k>0];
evens = FullSimplify[even,k>0];
say["summed responses for H=mu=eta=1 and b=H/4", {odds,evens}];
kernelSamples = N[{odds,evens}/.k->1,16];
(* Sum the massive real-space Green functions directly, divided by the
   response of the even pair's massless mode. H=1 and mu=eta. *)
distanceSamples=Table[N[Sum[Exp[-(2 j+1) Pi r],{j,0,Infinity}],18],{r,{1,2,4}}];
say["ratios at one, two and four slab widths", distanceSamples];

(* Additional longitudinal-normal control. A harmonic gradient has zero curl
   and divergence. Its virtual work relates the face tractions to an imposed
   unit normal body force at the midplane, for any static compatible solution. *)
xx=-Sinh[k z]; ww=Cosh[k z];
strains=Simplify[{D[xx,z]+k ww,D[ww,z]+k xx}];
faceJump = jump /. First[Solve[(ww/.z->0)+jump (ww/.z->hh/2)==0,jump]];
carried=1+faceJump;
say["zero-work strains, face traction jump and mixed share",{strains,faceJump,carried}];
(* Independent Dirichlet expansion for the induced shear stress from the
   normal load. The coefficient below uses integral cos times cos(0). *)
integrated = Total[Table[4 (-1)^j/((2j+1)Pi) / (1+((2j+1)Pi)^2),
                         {j,0,799}]//N];
normalSamples=N[{carried,faceJump}/.{hh->1,k->1},16];
say["spectral mixed share vs virtual-work mixed share at kH=1",
    {integrated,First[normalSamples]}];
(* Exact zero-frequency compliance and residue for the variable-profile
   control used in the other engine on -1<w<1. The static stress obeys
   stress'= -load at k=0; a unit dipole pair gives unit stress in between. *)
limits={Integrate[1/(1+z^2),{z,-1/2,1/2}],4/Integrate[2+z^2,{z,-1,1}]};
say["variable-modulus limits",limits];
samples=Flatten[{kernelSamples,distanceSamples,normalSamples,N[limits,18]}];
Print["comparison ",StringRiffle[ToString[#,InputForm]& /@ N[samples,MachinePrecision]," "]];
If[!TrueQ[bcs=={0,0,0,0} && residues=={0,4/(hh eta)} &&
  dispResidual==0 && strains=={0,0} && Abs[integrated-First[normalSamples]]<10^-9 &&
  Limit[k^2 odds,k->0]==0 && Limit[k^2 evens,k->0]==4],Exit[1]];
