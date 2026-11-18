(* Cartesian differentiation of the wave amplitude at a point on S^3.
   Only its derivatives through order two enter the stress and energy. *)
ClearAll["Global`*"];
vars = {x1, x2, x3, x4};
origin = Thread[vars -> 0];
normal = {z, Sqrt[1-z^2], 0, 0};
dxr = normal.vars;
dr = dxr + (vars.vars-dxr^2)/(2 radius);
cosine = (radius z+x1)/(radius+dr);
nBase = amp^2;
amplitude0 = amp (1+l1 dr/2+(l2/4+l1^2/8) dr^2);
amplitude1 = amplitude0 (b0+b1 dr+b2 dr^2/2) cosine/2;
velocityPotential0 = speed dr;
velocityPotential1 = (radius tangent + radial dr) cosine;
value[f_] := Factor[f/.origin];
gradient[f_] := value[Grad[f,vars]];
hessianColumn[f_] := value[Grad[D[f,x1],vars]].normal;
laplacian[f_] := value[Sum[D[f,{vars[[i]],2}],{i,4}]];

v0 = gradient[velocityPotential0]; v1 = gradient[velocityPotential1];
a0 = value[amplitude0]; a1 = value[amplitude1];
g0 = gradient[amplitude0]; g1 = gradient[amplitude1];
(* Quantum chemical potential -2 k Laplacian(amplitude)/amplitude.
   Pressure varies as n delta h with h the enthalpy per constituent. *)
dq = -2 k (laplacian[amplitude1]/a0-laplacian[amplitude0] a1/a0^2);
dp = nBase (-mass v0.v1-dq);
(* Quantum tensor 2 k [grad(a) grad(a) - a Hess(a)]. *)
stressX = 2 k (First[g1] (g0.normal) + First[g0] (g1.normal)
              - a1 hessianColumn[amplitude0]-a0 hessianColumn[amplitude1]);
carriageX = mass (2 a0 a1 First[v0] (v0.normal)
              + nBase (First[v1] (v0.normal)+ First[v0] (v1.normal)));
mean[expr_] := Factor[2/Pi Integrate[Expand[expr] Sqrt[1-z^2],{z,-1,1}]];
area = 2 Pi^2 radius^3;
channels = Factor /@ (area (mean /@ {carriageX, dp z, stressX}));
Print["outward carriage, pressure, density stress = ",ToString[channels,InputForm]];
check[name_, lhs_, rhs_] := Module[{res=Factor[lhs-rhs]},
  Print[name," | ",ToString[lhs,InputForm]," | ",ToString[rhs,InputForm],
        " | residual ",ToString[res,InputForm]];
  If[res =!= 0,Exit[1]]];
rate = area mass nBase speed;
total = Factor[Total[channels]];
reference = rate (radial+3 tangent+speed b0)/4+area k nBase (l1 b1-l2 b0)/4;
check["total flux on S^3",total,reference];
check["coefficient of density dipole curvature",Coefficient[total,b2],0];
check["coefficient of base density log curvature",Coefficient[total,l2],-area k nBase b0/4];
