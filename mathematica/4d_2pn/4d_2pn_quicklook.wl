(***
  Quick runner for 4d_gravity_2pn_master_harness.wl
  ------------------------------------------------
  Purpose:
    Load the main 2PN harness and print a compact digest of the key symbolic
    outputs that are likely to be useful while writing / checking the paper.
***)

ClearAll["Global`*"];

scriptDir = If[StringLength[$InputFileName] > 0, DirectoryName[$InputFileName], "."];
mainFileCandidates = DeleteDuplicates @ {
  FileNameJoin[{scriptDir, "4d_gravity_2pn_master_harness.wl"}],
  FileNameJoin[{scriptDir, "4d_gravity_2pn_master_harness_v6.wl"}],
  "4d_gravity_2pn_master_harness.wl",
  "4d_gravity_2pn_master_harness_v6.wl"
};
mainFile = SelectFirst[mainFileCandidates, FileExistsQ, Missing["NotFound"]];

If[MissingQ[mainFile],
  Print["FAIL: Could not locate a 4d_gravity_2pn_master_harness*.wl file"];
  Print["Tried: ", mainFileCandidates];
  Abort[];
];

Get[mainFile];

Print["\n--- 2PN QUICKLOOK ---"];
Print["q = ", TwoPNResults["qSolution"]];
Print["n = ", TwoPNResults["nSolution"]];
Print["kappa_PV = ", TwoPNResults["kappaPV"]];
Print["beta_1PN = ", TwoPNResults["beta1PN"]];
Print["Bernoulli index N(Phi) = ", TwoPNResults["BernoulliIndex"]];
Print["Self-sector candidate (dimensionless, U form) = ", TwoPNResults["Lself2PNScaledU"]];
Print["2PN two-body static coefficient = ", TwoPNResults["StaticCoeff2PNTwoBody"]];
Print["Isotropic pure-static fit for lambda_rho = ", TwoPNResults["lambdaRhoIsoFit"]];
Print["Direct self+static test-mass candidate (dynamic) = ", TwoPNResults["LtestMassCandidate2PNScaledDyn"]];
Print["Test-mass candidate from two-body reduction = ", TwoPNResults["LtestMassCandidate2PNScaledFromTwoBody"]];
Print["Residual after that static fit = ", TwoPNResults["ResidualCandidateToIsoStaticMatched"]];
Print["Finite-size / universal 2PN gate ratio = ", TwoPNResults["FiniteSizeToUniversal2PNRatio"]];
