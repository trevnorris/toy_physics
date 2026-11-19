(* Held triangle, rows and columns ordered by site. H20 = -t Exp[-I F].
   Positive F is a physical twist of the marked closing bond. The measured
   moments sum the full set of lines with multiplicity, without arm phases. *)
ClearAll["Global`*"];
emit[label_, value_] := Print[label, " = ", ToString[value, InputForm]];
checks = {};
check[value_] := AppendTo[checks, TrueQ[value]];
eye = IdentityMatrix[2];
pauli = Table[PauliMatrix[k], {k, 3}];
pin[p_] := Sum[p[[k]] pauli[[k]], {k, 3}];
pins = Table[Symbol["p" <> ToString[i] <> ToString[k]], {i, 3}, {k, 3}];

(* A term in Tr[H^n] is a cyclic word of n visited sites. Consecutive
   equal sites contribute a Pauli matrix. Hops are scalars. Pair a word
   winding w>0 with its reverse; its contribution to the HALF difference
   between F and -F is 2 Im[amplitude] Sin[w F]. Zero winding drops out. *)
walkTerm[word_] := Module[{steps, wind, sites, mats, n = Length[word], amp},
  steps = Transpose[{word, RotateLeft[word]}];
  wind = Count[steps, {3, 1}] - Count[steps, {1, 3}];
  If[wind <= 0, Return[0]];
  sites = Cases[steps, {i_, k_} /; i == k :> i];
  mats = (pin[pins[[#]]] &) /@ sites;
  amp = (-t)^(n-Length[sites]) (-j)^Length[sites] Tr[Fold[Dot, eye, mats]];
  2 Sin[wind f] ComplexExpand[Im[amp]]
];
odd = Table[Factor[Total[walkTerm /@ Tuples[Range[3], n]]], {n, 6}];
chi = Det[pins];
emit["Tr H^n half differences n=1..5; sixth divided by oriented triple product",
  {Most[odd], Factor[Last[odd]/chi]}];
check[Most[odd] === ConstantArray[0, 5] &&
  Expand[Last[odd] - 24 j^3 t^3 chi Sin[f]] === 0];

ham[p_, exchange_, hop_, z_] := ArrayFlatten[{
  {-exchange pin[p[[1]]], -hop eye, -hop Conjugate[z] eye},
  {-hop eye, -exchange pin[p[[2]]], -hop eye},
  {-hop z eye, -hop eye, -exchange pin[p[[3]]]}
}];
moments[m_] := Table[Tr[MatrixPower[m, n]], {n, 6}];
positive = IdentityMatrix[3];
mirrored = DiagonalMatrix[{1, -1, 1}];
planar = {{1, 0, 0}, {0, 1, 0}, {3/5, 4/5, 0}};
zForward = (3-4 I)/5;
exchange = 2; hop = 1; alpha = 3/2; baseSquared = 9;
h[p_, z_] := ham[p, exchange, hop, z];
fullDifference[p_] := moments[h[p, zForward]] - moments[h[p, Conjugate[zForward]]];
examples = fullDifference /@ {positive, mirrored, planar};
emit["finite H full differences n=1..6 for positive, mirror, coplanar pins", examples];
check[examples === {{0,0,0,0,0,1536/5}, {0,0,0,0,0,-1536/5}, ConstantArray[0,6]}];
secularDifference = Expand[CharacteristicPolynomial[h[positive,zForward],x] -
  CharacteristicPolynomial[h[positive,Conjugate[zForward]],x]];
emit["finite difference of det(x I-H)", secularDifference];
check[secularDifference === -Last[First[examples]]/6];

(* Four independent real displacements per site. Each bond compares its
   endpoint vectors after an orthogonal rotation. Kinetic energy is |q'|^2/2.
   Positive potential: (baseSquared-2 alpha hop) I - alpha exchange R(pin)
   is positive if baseSquared > alpha (2 hop + exchange Max[Norm[p_i]]).
   Add alpha hop ||q_i - R(z_ij) q_j||^2/2 for each bond i,j. *)
realBlock[m_] := ArrayFlatten[{{Re[m], -Im[m]}, {Im[m], Re[m]}}];
site[i_] := KroneckerProduct[{UnitVector[3, i]}, IdentityMatrix[4]];
bond[i_, k_, z_] := Module[{difference = site[i] - realBlock[z eye].site[k]},
  Transpose[difference].difference
];
stiffness[p_, z_] := Sum[Transpose[site[i]].(
    (baseSquared-2 alpha hop) IdentityMatrix[4] -
     alpha exchange realBlock[pin[p[[i]]]]).site[i], {i, 3}] +
  alpha hop (bond[1,2,1] + bond[2,3,1] + bond[3,1,z]);
realBySite[m_] := ArrayFlatten[Map[realBlock, Partition[m, {2,2}], {2}]];
kPlus = stiffness[positive,zForward];
kMinus = stiffness[positive,Conjugate[zForward]];
blockResidual = kPlus - baseSquared IdentityMatrix[12] - alpha realBySite[h[positive,zForward]];
doubledSecular = Expand[CharacteristicPolynomial[kPlus,x] -
  CharacteristicPolynomial[baseSquared IdentityMatrix[6]+alpha h[positive,zForward],x]^2];
lowerBound = baseSquared-alpha (2 hop + exchange Max[Norm /@ positive]);
emit["real reciprocal stiffness: symmetric, real, lower bound, block and squared secular residuals",
  {SymmetricMatrixQ[kPlus], FreeQ[kPlus,_Complex], lowerBound, Norm[Flatten[blockResidual]],doubledSecular}];
check[SymmetricMatrixQ[kPlus] && FreeQ[kPlus,_Complex] && lowerBound>0 &&
  Norm[Flatten[blockResidual]]===0 && doubledSecular===0];
dDelta = moments[kPlus]-moments[kMinus];
lines[m_] := Sort[Eigenvalues[N[m,50]]];
bothLines = lines /@ {kPlus,kMinus};
observed = Total[First[bothLines]^6] - Total[Last[bothLines]^6];
emit["full squared-frequency moment differences n=1..6, all twelve real modes",dDelta];
emit["six paired squared frequencies at +F,-F",N[bothLines[[All,1;;12;;2]],10]];
emit["full sixth difference from measured squared frequencies",N[observed,12]];
check[dDelta === 2 alpha^6 First[examples] && Abs[observed-Last[dDelta]]<10^-35];

(* Independent calibrations expressed in squared-frequency units. Compensate
   any diagonal shift when cutting bonds, or take the within-sector split.
   Which polarization is lower fixes the sign of exchange; dimer line
   splittings measure |t|. Zero-pin loop at a marked untwisted setting
   fixes the triangular hopping product sign (a single dimer cannot). *)
patch = baseSquared eye - alpha exchange pin[{0,0,1}];
barePair = baseSquared eye + alpha {{0,-hop},{-hop,0}};
refThird = 2 alpha^3 Tr[MatrixPower[h[ConstantArray[0,{3,3}],1],3]];
emit["isolated marked pin aligned/opposed squared frequencies; unpolarized dimer pair; zero-pin untwisted real third central moment",
  {Diagonal[patch],Sort[Eigenvalues[barePair]],refThird}];
check[Diagonal[patch]==={6,12} && Sort[Eigenvalues[barePair]]==={15/2,21/2} &&
  refThird===-24 (alpha hop)^3];

(* Scalar controls. A single real displacement only samples the cosine on
   the closing bond in this projection. A single COMPLEX scalar can carry
   a genuine rotation on the bond. With arbitrary real site detunings its
   characteristic polynomial is invariant under bond conjugation. The
   same holds for two uncoupled identical copies. *)
scalar[z_] := {{v1,-t,-t/z},{-t,v2,-t},{-t z,-t,v3}};
oneReal[c_] := {{v1,-t,-t c},{-t,v2,-t},{-t c,-t,v3}};
poly = Factor[CharacteristicPolynomial[scalar[z], x]];
controlResiduals = {oneReal[Cos[f]]===oneReal[Cos[-f]],
  Factor[poly-(poly /. z->1/z)],Factor[poly^2-(poly^2 /. z->1/z)]};
emit["one real scalar even; complex scalar and two identical copies' secular differences",controlResiduals];
check[controlResiduals==={True,0,0}];
emit["checks", checks];
If[!TrueQ[And@@checks], Exit[1]];
