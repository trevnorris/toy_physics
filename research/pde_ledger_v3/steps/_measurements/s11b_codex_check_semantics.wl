ClearAll[x, nested, extracted, pressurePowerResidual];
pressurePowerResidual = ComplexExpand[-Re[x]/2 + Re[x]/2, {x}];
nested = <|"OUTER" -> <|"TEST_OBJECT" -> False|>|>;
extracted = Flatten[Cases[Values[nested],
  Rule["TEST_OBJECT", value_] :> value, Infinity]];
Print["pressurePowerResidual = ", pressurePowerResidual];
Print["Cases extraction = ", extracted];
Print["And@@extraction = ", And @@ extracted];
