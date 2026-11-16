(* Curved sleeve. Variation of the bare energy in principal normal
   coordinates and local power accounting for pure material relaxation. *)
ClearAll["Global`*"];
well[c_] := k c^2 (1-c)^2/(2 ell^2);
profile = 1/(1+Exp[d/ell]);
(* Offset along an outward unit normal in three principal directions. *)
offsetH = 2 ka/(1+d ka)+km/(1+d km);
mu = well'[profile]-k (D[profile,{d,2}]+offsetH D[profile,d]);
loss = n mob mu;
profileResidual = Simplify[loss-n mob k offsetH profile (1-profile)/ell];
(* Differentiate the order energy with respect to a shift of its contour;
   the reversible order force is mu grad(chi). *)
speed = vel /. First[Solve[(vel-u) D[profile,d] == mob mu,vel]];
speedResidual = Simplify[speed-u+mob k offsetH];
workResidual = Simplify[mu D[profile,d] (speed-u)-mob mu^2];
Print["[loss, relative speed, work] residuals = ",
      ToString[{profileResidual,speedResidual,workResidual},InputForm]];
extraLoss = n lam profile (1-profile)(1-2 profile)^2;
Print["extra loss / (n lam) at center and d=ell log(3) = ",
      ToString[Simplify[{extraLoss /. d->0,extraLoss /. d->ell Log[3]}/(n lam)],InputForm]];
file = FileNameJoin[{"software","stage1_solver","frozen","m1c",
   "834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8","wp1_background_10x8.json"}];
branch = Import[file,"RawJSON"]["config"]["branch"];
If[branch["geometry_profile"] =!= "cubic_smoothstep", Exit[1]];
{lo,hi,r0,r1} = Rationalize[Lookup[branch,{"w_min","w_max","r_mouth","r_exit"}],10^-12];
radius = With[{x=(w-lo)/(hi-lo)},r0+(r1-r0)(3x^2-2x^3)];
(* Divergence of level-set normal to r-radius in four-dimensional space. *)
normal = {1,-D[radius,w]}/Sqrt[1+D[radius,w]^2];
midH = Simplify[(D[r^2 First[normal],r]/r^2+D[Last[normal],w]) /.
                  r->radius /. w->(lo+hi)/2];
(* Evaluate mu directly at the q=1/2 surface using the midpoint curvature. *)
midLoss = Simplify[ell (well'[profile]-k (D[profile,{d,2}]+midH D[profile,d]))/k /. d->0];
Print["archived midpoint H and loss ell/(n mob k) = ", ToString[{midH,midLoss},InputForm]];
If[{profileResidual,speedResidual,workResidual} =!= {0,0,0} || !NumberQ[N[midLoss]],Exit[1]];
