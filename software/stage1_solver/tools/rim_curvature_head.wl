(* Axisymmetric two-sphere collar in four space dimensions. *)
ClearAll["Global`*"];
check[label_, expr_] := Module[{res = FullSimplify[expr, Assumptions -> a > 0 && length > 0]},
  Print[label, " = ", ToString[res, InputForm]];
  If[res =!= 0, Exit[1]]];

(* A conical frustum ends at the moving rim. Its other end is fixed;
   subtract the interior area variation to extract the rim term. *)
rad[s_] := a0 + s (a - a0);
height[s_] := b0 + s (b - b0);
v = Sqrt[(a-a0)^2 + (b-b0)^2];
frustum = 4 Pi sig v (a0^2 + a0 a + a^2)/3;
(* Variation vector s*(da,db). For the straight meridian d_t tangent = 0. *)
tangent = {a-a0, b-b0}/v;
inside = 8 Pi sig Integrate[s rad[s], {s,0,1}] (v {1,0} - (a-a0) tangent);
edge = (D[frustum, #]& /@ {a,b}) - inside;
check["frustum radial endpoint force coefficient", edge[[1]] - 4 Pi sig a^2 tangent[[1]]];
check["frustum axial endpoint force coefficient", edge[[2]] - 4 Pi sig a^2 tangent[[2]]];

(* Spherical metric area variation at common displacement. *)
sphereArea = 4 Pi (a + eps da)^2;
collarWork = tau D[sphereArea,eps]/.eps->0;
check["collar energy variation per area", collarWork/(4 Pi a^2) - 2 tau da/a];

(* Derive mechanical force from the divergence of the order stress tensor:
   one normal eigenvalue S_N and isotropic tangential eigenvalue S_T. *)
sn = kap q'[z]^2/2 - f[q[z]];
st = -kap q'[z]^2/2 - f[q[z]];
geom = ar'[z]/ar[z];
forceN = -D[sn,z] - geom (sn-st);
mu = f'[q[z]] - kap (q''[z] + geom q'[z]);
check["tensor-divergence force minus mu q prime", forceN-mu q'[z]];
velocity = flow/(ar[z] den[z]);
ode = First[Solve[velocity q'[z] == -mob[z] mu, q''[z]]];
work = Together[ar[z] velocity (forceN /. ode)];
Print["mechanical power from the order force per tube length = ", ToString[Factor[work],InputForm]];
check["power plus order dissipation", work + ar[z] (velocity q'[z])^2/mob[z]];
check["integrand of head loss", (forceN /. ode)/den[z] + flow q'[z]^2/(ar[z] den[z]^2 mob[z])];
check["static stress jump integrand with curvature", (forceN /. ode) /. flow->0];
