(* Independent lab-frame readout of the raw stationary trace. Blocks are
   Re(field), Im(field), inward current, in the checkpoint's x,y,z order. *)
faceFile = "_scratch/s11c/mouth_rim_push_sympy01/resting_face.f64";
cells = 32;
spacing = 4/cells;
If[FileByteCount[faceFile] =!= 3*cells^3*8, Exit[1]];
rawTrace = BinaryReadList[faceFile, "Real64", ByteOrdering -> -1];
{re, im, inward} = ArrayReshape[rawTrace, {3, cells, cells, cells}];
density = re^2 + im^2;
coord = N[Table[(i - (cells+1)/2)*spacing, {i, cells}]];
weight = Table[Max[1-x^2-y^2-z^2, 0.]^3, {x,coord},{y,coord},{z,coord}];
mask = Map[Boole[# > 0] &, weight, {3}];
integral[grid_] := spacing^3 Total[Flatten[grid]];
dx[grid_] := (RotateLeft[grid,{1,0,0}]-RotateRight[grid,{1,0,0}])/(2 spacing);
numberRate = integral[inward];
vx = (re*dx[im] - im*dx[re])/density;
readLoad[offset_] := Module[{shift, profile, region, trial, logSlope},
  shift = {offset,0,0};
  profile = RotateRight[weight,shift];
  region = RotateRight[mask,shift];
  trial = density*RotateRight[inward/density,shift];
  logSlope = -2 profile*(1-1/Sqrt[density]);
  (* Periodic summation by parts places the mouth mask inside the
     derivative. No differentiated normal slope from the Python code is used. *)
  numberRate*integral[region*trial*vx]/integral[trial] +
    integral[logSlope*dx[region*density]]/4
];
loads = readLoad /@ {-1,0,1};
Print["inward_number_rate ", ToString[numberRate,InputForm]];
Print["frozen_brane_loads ", ExportString[loads,"RawJSON","Compact"->True]];
If[!TrueQ[Abs[numberRate-1/4]<10^-13 && Abs[loads[[2]]]<10^-13 &&
       Abs[loads[[1]]]>0 && Abs[(loads[[1]]+loads[[3]])/loads[[1]]]<10^-10], Exit[1]];

(* Solve for the carrier at complex frequency from zero initial increment.
   A general pair of coupled carrier modes exhibits the bounded resolvent;
   for any fixed finite matrix (s I - A)^(-1) = I/s + O(s^-2). *)
matrix = {{a11,a12},{a21,a22}};
source = {b1,b2};
readout = {c1,c2};
transfer = k + readout.LinearSolve[s IdentityMatrix[2]-matrix,source];
motion = accel/s^3 + jerk/s^4;
brane = Together[transfer*motion];
initialResidual = Limit[Together[s*(force/s + brane)],s->Infinity];
quadraticJet = Limit[Together[s^3*brane],s->Infinity];
cubicJet = Limit[Together[s^4*(brane-k*accel/s^3)],s->Infinity];
Print["step residual from the high-frequency limit = ", InputForm[initialResidual]];
Print["load second and third time derivatives = ",InputForm[{quadraticJet,cubicJet}]];

(* Construct the kernel of the endpoint measurements without choosing a bump
   in advance. The quartic coefficient sets the free normalization. *)
poly[z_] := p0 + p1*z + p2*z^2 + p3*z^3 + amp*z^4;
match = First@Solve[{poly[0]==0,poly'[0]==0,poly[1]==0,poly'[1]==0},{p0,p1,p2,p3}];
profile = poly[w/len] /. match;
deltaAccel = Simplify[(D[(t0+t1*w)*D[profile,w],w] /. w->0)/mu0];
Print["endpoint-compatible profile = ",InputForm[Factor[profile]]];
Print["endpoint acceleration difference = ", InputForm[deltaAccel]];
If[!TrueQ[initialResidual===force && Simplify[quadraticJet-k*accel]===0 &&
    Simplify[cubicJet-k*jerk-accel*readout.source]===0 &&
    Simplify[deltaAccel-2*amp*t0/(len^2*mu0)]===0], Exit[1]];
