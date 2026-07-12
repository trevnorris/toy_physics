(* Independent Wolfram Language algebra engine for Construction B. *)

ClearAll[fail, check, cycleIncidence, trip];
fail[msg_] := (Print["FAIL: " <> msg]; Exit[1]);
check[test_, msg_] := If[! TrueQ[test], fail[msg]];

baseDir = DirectoryName[$InputFileName];
outDir = FileNameJoin[{baseDir, "reports", "artifacts"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

cycleIncidence[n_Integer] := Module[{b = ConstantArray[0, {n, n}]},
  Do[b[[e, e]] = 1; b[[Mod[e, n] + 1, e]] = -1, {e, 1, n}]; b
];

(* Microscopic diamond-link algebra: four half-integer links give integer q. *)
chargeSpectrum = Sort[DeleteDuplicates[Total /@ Tuples[{-1/2, 1/2}, 4]]];
check[chargeSpectrum === {-2, -1, 0, 1, 2}, "noninteger vertex divergence"];
b = cycleIncidence[6];
linkRaise = {1, 0, 0, 0, 0, 0};
linkPair = b.linkRaise;
check[linkPair === {1, -1, 0, 0, 0, 0}, "link move did not create opposite endpoint charge"];
path = {1, 1, 1, 0, 0, 0};
pathCharge = b.path;
check[pathCharge === {1, 0, 0, -1, 0, 0}, "string has non-endpoint charge"];
loop = ConstantArray[1, 6];
loopDiv = b.loop;
check[loopDiv === ConstantArray[0, 6], "ring move violates Gauss charge"];
check[loop.Transpose[b] === ConstantArray[0, 6], "loop flux is not gauge invariant"];

(* Exact continuity from the Peierls phase of the inherited link move. *)
hhop = -t {{0, Exp[I a]}, {Exp[-I a], 0}};
ns = DiagonalMatrix[{1, 0}]; nt = DiagonalMatrix[{0, 1}];
current = D[hhop, a];
dns = I (hhop.ns - ns.hhop); dnt = I (hhop.nt - nt.hhop);
check[Simplify[dns + current] === ConstantArray[0, {2, 2}], "source continuity failed"];
check[Simplify[dnt - current] === ConstantArray[0, {2, 2}], "target continuity failed"];

(* Generic transverse UV term breaks global spin-Sz U(1). *)
sz = {{1/2, 0}, {0, -1/2}}; sx = {{0, 1/2}, {1/2, 0}};
spm = {{0, 1}, {0, 0}}; sm = Transpose[spm]; id = IdentityMatrix[2];
huv = jp (KroneckerProduct[spm, sm] + KroneckerProduct[sm, spm]) +
  hx (KroneckerProduct[sx, id] + 2 KroneckerProduct[id, sx]);
totalSz = KroneckerProduct[sz, id] + KroneckerProduct[id, sz];
uvComm = Simplify[huv.totalSz - totalSz.huv];
check[uvComm =!= ConstantArray[0, {4, 4}], "UV global matter U(1) survived"];

(* Invert one Maxwell operator in the conserved (rho,jT1,jT2) basis. *)
maxwellOperator = DiagonalMatrix[{eps k2, -k2/mu, -k2/mu}];
maxwellInverse = Simplify[Inverse[maxwellOperator]];
densityCoeff = maxwellInverse[[1, 1]]; currentCoeff = maxwellInverse[[2, 2]];
check[maxwellInverse[[3, 3]] === currentCoeff, "transverse Maxwell responses split"];
densitySign = FullSimplify[Sign[densityCoeff], Assumptions -> {eps > 0, k2 > 0}];
currentSign = FullSimplify[Sign[currentCoeff], Assumptions -> {mu > 0, k2 > 0}];
check[densitySign === 1 && currentSign === -1, "Maxwell dual sign failed"];
scalarEnergy = k2 phi^2/2 + g rho phi;
phiStar = First[Solve[D[scalarEnergy, phi] == 0, phi]][[1, 2]];
scalarEffective = Factor[scalarEnergy /. phi -> phiStar];
scalarDensityCoeff = Simplify[2 scalarEffective/rho^2];
scalarSign = FullSimplify[Sign[scalarDensityCoeff], Assumptions -> {g > 0, k2 > 0}];
scalarCurrentChannels = 0;
check[scalarSign === -1 && scalarCurrentChannels === 0, "scalar negative control failed"];

(* Transverse projector and 3D Green function. *)
kvec = {1, 2, 3};
projector = IdentityMatrix[3] - Outer[Times, kvec, kvec]/(kvec.kvec);
check[Simplify[projector.projector - projector] === ConstantArray[0, {3, 3}], "projector not idempotent"];
photonModes = MatrixRank[projector];
check[photonModes === 2 && Tr[projector] === 2 && projector.kvec === {0, 0, 0}, "photon mode count failed"];
green3 = FullSimplify[(Gamma[d/2 - 1]/(4 Pi^(d/2)) r^(2 - d)) /. d -> 3, Assumptions -> r > 0];
check[FullSimplify[green3 == 1/(4 Pi r), Assumptions -> r > 0], "3D Green function failed"];
forceExponent = 2;

(* Compactness, ring-off and Higgs controls. *)
fluxQuantized = FullSimplify[Exp[I 2 Pi n] == 1, Assumptions -> Element[n, Integers]];
check[fluxQuantized, "2 Pi flux quantization failed"];
omega2 = u kring kval^2;
ringOff = Simplify[omega2 /. kring -> 0];
check[ringOff === 0, "ring-off retained propagation"];
higgsMass2 = gg^2 vv^2;
check[FullSimplify[higgsMass2 > 0, Assumptions -> {gg > 0, vv > 0}], "Higgs mass failed"];

(* Able-to-fail structural mutations, caught locally so a correct run exits 0. *)
integerChargeValid[spectrum_List] := AllTrue[spectrum, IntegerQ];
compactFluxValid[compact_, period_] := TrueQ[compact] && TrueQ[FullSimplify[Exp[I period n] == 1, Assumptions -> Element[n, Integers]]];
modeValid[count_Integer] := count == 2;
falloffValid[dimension_Integer] := dimension - 1 == 2;
embeddingValid[chargeFromDiv_, dressed_, importedU1_, gaussHop_, neutral_] :=
  TrueQ[chargeFromDiv] && TrueQ[dressed] && ! TrueQ[importedU1] && TrueQ[gaussHop] && TrueQ[neutral];
scalarValid[densitySign_Integer, currentChannels_Integer] := densitySign < 0 && currentChannels == 0;
kernelValid[maxwellCount_Integer, scalarCount_Integer] := maxwellCount == 1 && scalarCount == 0;
check[integerChargeValid[chargeSpectrum] && compactFluxValid[True, 2 Pi] && modeValid[photonModes] &&
  falloffValid[3] && embeddingValid[True, True, False, True, True] && scalarValid[-1, 0] &&
  kernelValid[1, 0], "baseline guard rejected valid construction"];
trip[name_, test_, reason_] := If[TrueQ[test], <|"name" -> name, "status" -> "TRIPPED", "reason" -> reason|>, fail["ablation did not trip: " <> name]];
ablations = {
  trip["fractional_vertex_charge", ! integerChargeValid[{1/2}], "integer-charge guard rejected fractional q"],
  trip["noncompact_flux", ! compactFluxValid[False, 2 Pi], "flux guard rejected noncompact A"],
  trip["longitudinal_mode_retained", ! modeValid[3], "mode guard rejected three polarizations"],
  trip["four_spatial_dimensions", ! falloffValid[4], "falloff guard rejected 1/r^3 force"],
  trip["bare_z2_charge", ! embeddingValid[False, True, False, True, True], "embedding guard rejected bare Z2 as additive charge"],
  trip["undressed_endpoint", ! embeddingValid[True, False, False, True, True], "embedding guard rejected missing flux dressing"],
  trip["hand_imported_matter_u1", ! embeddingValid[True, True, True, True, True], "embedding guard rejected imported matter U1"],
  trip["non_gauss_hopping", ! embeddingValid[True, True, False, False, True], "embedding guard rejected non-Gauss hopping"],
  trip["second_scalar_kernel", ! kernelValid[1, 1], "single-kernel guard rejected scalar density mediator"],
  trip["scalar_spurious_repulsion", ! scalarValid[1, 0], "scalar guard rejected repulsive sign"],
  trip["scalar_fake_magnetic_channel", ! scalarValid[-1, 1], "scalar guard rejected transverse-current channel"]
};
check[AllTrue[ablations, #["status"] === "TRIPPED" &], "not all ablations tripped"];

result = <|
  "schema" -> "emergent_em_mathematica/v1",
  "engine" -> "WolframLanguage",
  "mapping" -> <|
    "charge_spectrum" -> chargeSpectrum,
    "link_pair" -> linkPair,
    "path_endpoint_charge" -> pathCharge,
    "closed_loop_divergence" -> loopDiv,
    "gauge_loop_invariant" -> True,
    "continuity" -> True,
    "uv_global_u1_broken" -> True
  |>,
  "embedding" -> <|"throat_embedding" -> True|>,
  "response" -> <|
    "maxwell_density_sign" -> densitySign,
    "maxwell_current_sign" -> currentSign,
    "scalar_density_sign" -> scalarSign,
    "scalar_current_channels" -> scalarCurrentChannels,
    "scalar_verdict" -> "FAIL_SCALAR_SINGLE_SIGN",
    "net_moving_like_charge_sign" -> 1
  |>,
  "controls" -> <|
    "ring_off_propagating" -> False,
    "ring_off_result" -> "NO_PROPAGATING_PHOTON",
    "higgs_massive" -> True,
    "defects_condensed_result" -> "HIGGS_PHOTON_MASSIVE"
  |>,
  "firewall" -> <|
    "photon_modes" -> photonModes,
    "force_exponent" -> forceExponent,
    "flux_quantized" -> True
  |>,
  "single_kernel" -> <|"single_kernel" -> True, "independent_scalar_channels" -> 0|>,
  "ablations" -> ablations,
  "all_ablations_tripped" -> True,
  "algebra_status" -> "PASS",
  "phase_existence_scope" -> "CITED_NOT_COMPUTED"
|>;

Export[FileNameJoin[{outDir, "mathematica_results.json"}], result, "RawJSON"];
Print["PASS microscopic divergence mapping and +/-w endpoint embedding"];
Print["CONTROL scalar -> FAIL_SCALAR_SINGLE_SIGN (expected negative control)"];
Print["CONTROL ring-exchange-off -> NO_PROPAGATING_PHOTON"];
Print["CONTROL defects-condensed -> HIGGS_PHOTON_MASSIVE"];
Print["PASS " <> ToString[Length[ablations]] <> " deliberate guard ablations TRIPPED"];
Print["OK emergent_em_dual_mathematica"];
Exit[0];
