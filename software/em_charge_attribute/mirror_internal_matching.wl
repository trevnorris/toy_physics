(* Independent variations in real quadratures and in amplitude and phase. *)
ClearAll["Global`*"];
check[name_, expr_, ref_] := Module[{v = FullSimplify[expr]},
  Print[name, " = ", ToString[v, InputForm]];
  If[!TrueQ[FullSimplify[v == ref]], Exit[1]]
];

(* A local phase rotation or a phase rotation after conjugation. Use a
   common lab coordinate w and inward normal s. cphase is +1 or -1. *)
Do[
  ph = const + cphase theta[s w];
  dens = amp[s w]^2;
  inward = s dens D[ph, w] /. w -> s z;
  check[{"side", s, "phase sign", cphase, "inward coefficient"}, inward,
        cphase amp[z]^2 theta'[z]],
  {s, {1, -1}}, {cphase, {1, -1}}];

(* Normal derivative jump from stationarity of a real surface energy.
   h and hN are its first and second density derivatives at the trace. *)
xy = {x, y};
surface = g[x^2+y^2, chi];
normalJump = (D[surface,#] & /@ xy)/(2 cap);
check["material current jump in units hbar/m", {-y, x}.normalJump, 0];
rules = {Derivative[1,0][g][x^2+y^2,chi] -> h,
         Derivative[2,0][g][x^2+y^2,chi] -> hN,
         Derivative[1,1][g][x^2+y^2,chi] -> hChi};
complexJump = cap ({1,I}.normalJump);
direct = ((D[complexJump,x] - I D[complexJump,y])/2) /. rules;
conjugate = ((D[complexJump,x] + I D[complexJump,y])/2) /. rules;
orderTerm = D[complexJump,chi] /. rules;
check["direct coefficient", direct, h+(x^2+y^2) hN];
check["conjugate coefficient", conjugate, hN (x+I y)^2];
check["order coefficient", orderTerm, hChi (x+I y)];
check["derivative along a common phase orbit",
      direct I (x+I y)-conjugate I (x-I y), I h (x+I y)];

(* Use the quadrature energy norm of the free cavity and its two ends.
   Relative phase changes the two temporal quadratures at fixed amplitude. *)
k = Pi/(2 len);
shape[z_] := Sin[k z];
ends = Table[{shape[s w] /. w -> 0, D[shape[s w],w] /. w -> s len}, {s, {1,-1}}];
check["both cavity boundary rows", ends, {{0,0},{0,0}}];
modalMass = inertia Integrate[shape[z]^2, {z,0,len}];
modalSpring = inertia speed^2 Integrate[shape'[z]^2, {z,0,len}];
modeEnergy = (modalMass (amplitude speed k Sin[az])^2 + modalSpring (amplitude Cos[az])^2)/2;
check["standing-mode energy", modeEnergy, Pi^2 inertia speed^2 amplitude^2/(16 len)];

(* Probe the carrier time derivative with both real components. The other
   energy terms enter separate volume, surface and mode integrals. *)
psi = x + I y;
dt = xt + I yt + I q pot psi/hb;
timeL = I hb (Conjugate[psi] dt - Conjugate[dt] psi)/2;
source = -D[ComplexExpand[timeL], pot];
rimEnergy = km hh^2/2 - drive side hh;
check["static matter source", source, q (x^2+y^2)];
check["contact, rim, wave direct monopole",
      D[{surface, rimEnergy, modeEnergy}, pot], {0,0,0}];
check["mechanical rim drive at H=0", -D[rimEnergy,hh] /. hh -> 0, drive side];
(* If F' is the excess density and G'=w F', integrate mirrored lab
   profiles with fixed lab probe 1+grad w. Use positive volume on each side. *)
primitive[s_, w_] := q (s ff[s w] + grad gg[s w]);
upper = primitive[1,right]-primitive[1,left];
lower = primitive[-1,-left]-primitive[-1,-right];
check["uniform excess on both mirrors", {upper,lower} /. grad -> 0,
      {1,1} q (ff[right]-ff[left])];
check["depth-gradient half difference per gradient", D[(upper-lower)/2, grad],
      q (gg[right]-gg[left])];
