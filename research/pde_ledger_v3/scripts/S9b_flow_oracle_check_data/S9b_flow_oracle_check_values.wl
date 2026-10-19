(* Value comparison of the final baseline and FORM-control transcripts.
   No physical expression is transcribed here. Every operand and domain is
   read from a printed residual. Run only through the pooled guard. *)
$HistoryLength=0;
SetSystemOptions["ParallelOptions"->{"ParallelThreadNumber"->1,"MKLThreadNumber"->1}];
SetSystemOptions["InequalitySolvingOptions"->"ThreadOr"->False];
emit[tag_,x_]:=Print[tag,": ",ToString[x,InputForm,PageWidth->Infinity]];
base="/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/";
paths=<|"baseline"->base<>"S9b_flow_oracle_check.out",
 "visser_form"->base<>"S9b_flow_oracle_check_visser_form.out"|>;
emit["VALUE_COMPARISON_PROVENANCE",Append[
 Map[<|"path"->#,"sha256"->IntegerString[FileHash[#,"SHA256"],16,64]|>&,paths],
 "script_sha256"->IntegerString[FileHash[$InputFileName,"SHA256"],16,64]]];
lines=ReadList[#,String]&/@paths;
parseLine[s_]:=ToExpression[StringDrop[s,First[StringPosition[s,": "]][[2]]],InputForm];
records[ss_,tag_]:=Map[Function[i,<|"line"->i,"value"->parseLine[ss[[i]]],
 "domain"->parseLine[ss[[i-1]]],"zero_test"->parseLine[ss[[i+1]]]|>],
 Flatten[Position[StringStartsQ[#,tag<>": "]&/@ss,True]]];
checks={};
compare[variant_,tag_,indices_]:=Module[{before,after,b,a,dom,delta,test},
 before=records[lines["baseline"],tag];after=records[lines[variant],tag];
 Do[b=before[[i]];a=after[[i]];dom=b["domain"];
  emit["WORKING_ON",{variant,tag,i}];
  emit["FORM_VALUE_OPERANDS",<|"baseline"->b,"control"->a|>];
  emit["RESIDUAL_ASSUMPTIONS",dom];
  ClearSystemCache[];
  delta=FullSimplify[a["value"]-b["value"],Assumptions->dom["CAS"],TimeConstraint->Infinity];
  emit["FORM_VALUE_DELTA",delta];
  test=FullSimplify[delta==0,Assumptions->dom["CAS"],TimeConstraint->Infinity];
  emit["RESIDUAL_ZERO_TEST",<|"residual_tag"->"FORM_VALUE_DELTA","domain"->dom["id"],"value"->test|>];
  ClearSystemCache[];
  AppendTo[checks,{variant,tag,i,"same_printed_domain",SameQ[dom,a["domain"]]}],{i,indices}]];
(* Fixed selection, independent of residual values: all Visser balance cases.
   The Hamilton-Lisle FORM control is dropped by user instruction. *)
compare["visser_form","EXCHANGE_BALANCE_OPERATOR_RESIDUAL",
 Range[Length[records[lines["baseline"],"EXCHANGE_BALANCE_OPERATOR_RESIDUAL"]]]];
emit["VALUE_COMPARISON_GUARDS",checks];
If[!And@@checks[[All,-1]],Exit[2]];
emit["VALUE_COMPARISON_COMPLETION",Length[checks]];
Exit[0];
