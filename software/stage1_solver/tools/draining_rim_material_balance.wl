(* Construct a mass-current variation in Cartesian brane coordinates. The
   radial profile g equals f/r inside the band and vanishes outside it. *)
coords = {x,y,z};
rad = Sqrt[coords.coords];
potential = b g[rad] Cross[{1,0,0},coords];
current = Curl[potential,coords];
div = Simplify[Div[current,coords]];
normalCheck = Simplify[current.coords - 2 b x g[rad]];
rimProfile = {g[a]->0, g'[a]->1/a};
equator = FullSimplify[current /. {x->0,y->a,z->0},Assumptions->a>0] /. rimProfile;
profile[r_] := (r-a)(1-((r-a)/ell)^2)^4/r;
profileCheck = FullSimplify[{profile[a],profile'[a],profile[a-ell],profile'[a-ell],
                  profile[a+ell],profile'[a+ell]},Assumptions->a>ell>0];
Print["Cartesian divergence and radial identity residual: ",{div,normalCheck}];
Print["Equatorial mass-current increment: ",equator];
Print["Radial profile and slope at rim, inner edge, outer edge: ",InputForm[profileCheck]];
If[div =!= 0 || normalCheck =!= 0 || equator =!= {b,0,0} ||
   profileCheck =!= {0,1/a,0,0,0,0},Exit[1]];
