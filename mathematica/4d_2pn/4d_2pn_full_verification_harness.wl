(**
  4D -> full derivation verification harness
  ------------------------------------------
  Purpose:
    Execute the verified Mathematica derivation scripts listed in wl_notes.txt
    in order, collect their PASS/FAIL output, and report one suite-level result.

  Notes:
    - 4d_2pn_quicklook.wl is explicitly excluded.
    - Each script is executed in a fresh Global` workspace.
    - Abort[] is treated as a suite failure for that script.
**)

Begin["FourDTwoPNFullVerificationHarness`Private`"];

ClearAll[
  suiteSection, suitePass, suiteFail, suiteSkip, suiteInfo,
  manifestFiles, excludedFiles, filePassCount, fileFailCount, fileSkipCount,
  totalDetectedPassLines, totalDetectedFailLines,
  scriptDirectory, readManifest, passLineQ, failLineQ, runScript, runSuite
];

filePassCount = 0;
fileFailCount = 0;
fileSkipCount = 0;
totalDetectedPassLines = 0;
totalDetectedFailLines = 0;
excludedFiles = {"4d_2pn_quicklook.wl"};

suiteSection[name_String] := Print["\n=== ", name, " ==="];
suitePass[name_String] := (filePassCount++; Print["PASS: ", name]);
suiteFail[name_String, res_: Missing["NotAvailable"]] := (
  fileFailCount++;
  If[MissingQ[res],
    Print["FAIL: ", name],
    Print["FAIL: ", name, "\n  residual -> ", res]
  ]
);
suiteSkip[name_String, msg_String : ""] := (
  fileSkipCount++;
  If[msg === "",
    Print["SKIP: ", name],
    Print["SKIP: ", name, "\n  ", msg]
  ]
);
suiteInfo[msg_String] := Print["INFO: ", msg];

scriptDirectory[] := Module[{input = $InputFileName},
  If[StringQ[input] && StringLength[input] > 0, DirectoryName[input], Directory[]]
];

readManifest[manifestPath_String] := Module[{lines},
  If[!FileExistsQ[manifestPath], Return[$Failed]];
  lines = Quiet @ Import[manifestPath, "Lines"];
  If[!ListQ[lines], Return[$Failed]];
  DeleteDuplicates @ Select[
    StringTrim /@ lines,
    # =!= "" && !StringStartsQ[#, "#"] && !MemberQ[excludedFiles, #] &
  ]
];

passLineQ[line_String] := Module[{trim = StringTrim[line]},
  StringContainsQ[trim, "PASS:"] ||
  StringContainsQ[trim, ": PASS"] ||
  (StringContainsQ[trim, " PASS"] && !StringContainsQ[trim, "PASSED"])
];

failLineQ[line_String] := Module[{trim = StringTrim[line]},
  StringContainsQ[trim, "FAIL:"] ||
  StringContainsQ[trim, ": FAIL"] ||
  (StringContainsQ[trim, " FAIL"] && !StringContainsQ[trim, "FAILED"])
];

runScript[file_String, suiteDir_String] := Module[
  {
    absPath, tempPath, stream, previousDir, outputText = "", outputLines,
    passLines, failLines, aborted = False, status, previewLines
  },
  suiteSection["RUN / " <> file];
  absPath = FileNameJoin[{suiteDir, file}];

  If[MemberQ[excludedFiles, file],
    suiteSkip[file, "excluded from the suite"];
    Return[<|"File" -> file, "Status" -> "skipped", "Aborted" -> False, "PassLineCount" -> 0, "FailLineCount" -> 0|>]
  ];

  If[!FileExistsQ[absPath],
    suiteFail[file, "missing file"];
    Return[<|"File" -> file, "Status" -> "failed", "Aborted" -> False, "PassLineCount" -> 0, "FailLineCount" -> 0|>]
  ];

  tempPath = CreateTemporary[];
  stream = OpenWrite[tempPath, CharacterEncoding -> "UTF8", PageWidth -> Infinity];
  previousDir = Directory[];

  SetDirectory[suiteDir];
  CheckAbort[
    Block[
      {
        $Context = "Global`",
        $ContextPath = {"System`", "Global`"},
        $Assumptions = True,
        $Output = {stream},
        $Messages = {stream}
      },
      Quiet @ ClearAll["Global`*"];
      Get[absPath];
    ],
    aborted = True
  ];
  SetDirectory[previousDir];
  Close[stream];

  outputText = Quiet @ Check[Import[tempPath, "Text"], ""];
  Quiet @ DeleteFile[tempPath];

  outputLines = DeleteCases[StringTrim /@ StringSplit[ToString[outputText], RegularExpression["\\r\\n|\\n|\\r"]], ""];
  passLines = Select[outputLines, passLineQ];
  failLines = Select[outputLines, failLineQ];

  totalDetectedPassLines += Length[passLines];
  totalDetectedFailLines += Length[failLines];

  suiteInfo[
    "Detected " <> ToString[Length[passLines]] <> " PASS lines and " <>
    ToString[Length[failLines]] <> " FAIL lines."
  ];

  status = Which[
    aborted, "failed",
    Length[failLines] > 0, "failed",
    Length[passLines] == 0, "failed",
    True, "passed"
  ];

  If[status === "passed",
    suitePass[file],
    previewLines = Which[
      aborted, {"script aborted before completion"},
      Length[failLines] > 0, Take[failLines, Min[5, Length[failLines]]],
      True, Take[outputLines, -Min[5, Length[outputLines]]]
    ];
    suiteFail[file, previewLines]
  ];

  <|
    "File" -> file,
    "Status" -> status,
    "Aborted" -> aborted,
    "PassLineCount" -> Length[passLines],
    "FailLineCount" -> Length[failLines],
    "FailPreview" -> If[Length[failLines] > 0, Take[failLines, Min[5, Length[failLines]]], {}]
  |>
];

runSuite[] := Module[{suiteDir, manifestPath, results},
  suiteDir = scriptDirectory[];
  manifestPath = FileNameJoin[{suiteDir, "wl_notes.txt"}];

  suiteSection["SETUP"];
  suiteInfo["Harness directory: " <> suiteDir];
  suiteInfo["Manifest: " <> manifestPath];
  suiteInfo["Excluded files: " <> StringRiffle[excludedFiles, ", "]];

  manifestFiles = readManifest[manifestPath];
  If[manifestFiles === $Failed,
    suiteFail["Unable to read wl_notes.txt", manifestPath];
    Return[<|"Manifest" -> {}, "Results" -> {}, "FilesPassed" -> filePassCount, "FilesFailed" -> fileFailCount, "FilesSkipped" -> fileSkipCount|>]
  ];

  suiteInfo["Loaded " <> ToString[Length[manifestFiles]] <> " script entries from wl_notes.txt."];
  results = runScript[#, suiteDir] & /@ manifestFiles;

  suiteSection["SUMMARY"];
  Do[
    Print[
      result["File"], " -> ", ToUpperCase[result["Status"]],
      "  (PASS lines: ", result["PassLineCount"],
      ", FAIL lines: ", result["FailLineCount"],
      If[TrueQ[result["Aborted"]], ", aborted", ""], ")"
    ],
    {result, results}
  ];

  Print["Files passed: ", filePassCount];
  Print["Files failed: ", fileFailCount];
  Print["Files skipped: ", fileSkipCount];
  Print["Detected PASS lines: ", totalDetectedPassLines];
  Print["Detected FAIL lines: ", totalDetectedFailLines];

  If[fileFailCount == 0,
    Print["\nALL SUITE CHECKS PASSED."],
    Print["\nSOME SUITE CHECKS FAILED."]
  ];

  <|
    "Manifest" -> manifestFiles,
    "Results" -> results,
    "FilesPassed" -> filePassCount,
    "FilesFailed" -> fileFailCount,
    "FilesSkipped" -> fileSkipCount,
    "DetectedPassLines" -> totalDetectedPassLines,
    "DetectedFailLines" -> totalDetectedFailLines
  |>
];

Global`VerificationSuiteResults = runSuite[];

End[];
