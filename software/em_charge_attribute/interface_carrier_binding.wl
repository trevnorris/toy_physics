(* Independent fixed-flux reduction of the stationary density/phase energy. *)
ClearAll["Global`*"];
check[label_, expr_, expected_] := Module[{v = Factor[expr]},
  Print[label, " = ", ToString[v, InputForm]];
  If[!TrueQ[FullSimplify[v == expected]], Exit[1]]
];

bulk = hb^2 dn^2/(8 m n) + hb^2 n phaseSlope^2/(2 m) + eos n^5/4;
phaseRule = First[Solve[D[bulk,phaseSlope]==hb j,phaseSlope]];
fixedFlux = (bulk-hb j phaseSlope) /. phaseRule;
coeff = D[fixedFlux,{n,2}] /. {n->n0,dn->0};
grad = D[fixedFlux,{dn,2}] /. n->n0;
check["local density stiffness at fixed flux", coeff, 5 eos n0^3 - m j^2/n0^3];
check["gradient coefficient", grad,hb^2/(4 m n0)];

(* Integrate the weak stationary equation against a constant. The outer
   normal gradients decay. Their jump at the contact is fixed by variation
   of the real surface energy. B is the integral of the density departure. *)
integratedVariation = k stored + grad (slopeRight - slopeLeft);
surfaceVariation = grad (slopeLeft-slopeRight) + h;
perArea = stored /. First[Solve[integratedVariation+surfaceVariation==0,stored]];
check["integrated number departure per area",perArea,-h/k];
check["fixed-intake patch stiffness",
      coeff /. j->rate/area,5 eos n0^3-m rate^2/(area^2 n0^3)];

(* The local normal and velocity both reverse on the reflected branch. *)
check["reflected outward flux difference",(-v)(-normal)-v normal,0];
pop[sign_] := area (perArea /. h->common[radius,energy]+sign odd[radius,energy]);
needed = First[Solve[(pop[1]-pop[-1])/2==target,odd[radius,energy]]];
check["required odd surface slope for a fixed half-difference",
      odd[radius,energy]/.needed,-target k/area];
check["half-difference of equal-contact mirror charges",
      q (pop[1]-pop[-1])/2 /. odd[radius,energy]->0,0];

(* For any positive density in a transmitting stream dz/v = n dz/j.
   This equality also applies to the subtracted uniform reference. *)
check["charge beyond residence contribution per depth",
      q area (density[z]-n0) - q (area j) (1/(j/density[z])-1/(j/n0)),0];
