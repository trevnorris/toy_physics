(* Static shear contact in the spherical sections of the radius-one guide.
   Y is the direction cosine. A covector is written u_r Y dr + b dY + u_w Y dw. *)
ClearAll["Global`*"];
f[r_, z_] = h[z] r (r^2 - a[z]^2)/(2 a[z]^2 (1 + a'[z]^2));
ur = D[f[r,z],r]; b = f[r,z]; uw = D[f[r,z],z];
normal = FullSimplify[(ur-a'[z] uw /. r->a[z])-h[z],Assumptions->a[z]>0];
curl = Simplify[{D[ur,z]-D[uw,r], D[b,r]-ur, D[b,z]-uw}];
(* Pullback of h e_x on the sidewall: h a' Y dz + h a dY. *)
sideCurl = Expand[D[h[z] a[z],z]-h[z] a'[z]];
(* Homogeneous ordered patch, relative to a fixed reference, u_x=strain*w. *)
patchEnergy = modulus (strain-backTilt)^2/2;
stress = D[patchEnergy,strain];
print[key_, val_] := Print[key," = ",InputForm[val]];
print["error in prescribed normal displacement",normal];
print["three curl coefficients of normal-only extension",curl];
print["coefficient of dz wedge dY on fully bonded guide",sideCurl];
print["uniform patch traction with opposing tilt",stress];
If[normal=!=0 || curl=!={0,0,0} || Simplify[sideCurl-a[z] h'[z]]=!=0 || Simplify[(stress /. backTilt->0)-modulus strain]=!=0, Exit[1]];
Exit[0];
