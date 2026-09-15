# S11c-d numerical action validation

The exact one-case operator assembly is committed at d868f949: local derivative
orders 0–3 and 80 intact nonlocal integral operators. This continues step 3 of
the variable-profile matching plan. It does not replace its later boundary,
continuum-expansion or profile-frequency-pole requirements.

1. The user approved the proposed numerical instance for the 30 inherited free
   gradient-energy coefficients. Use the separate development input. Preserve all
   existing parameter values, independent profile formulas, unit frame and
   symbolic coefficient dependence. Verify the selected extension against the
   accepted end operators and normalization inputs before use.
2. Bind the actual assembled coefficient matrices and integral operands to the
   complete input. Resolve profile derivatives and end limits from the supplied
   formulas. Reject remaining unbound physical parameters. Retain the separate
   eta/sigma symbolic grades and record the evaluated homotopy origin.
3. Evaluate local derivative actions on explicit smooth localized test fields
   in each of the five field directions. Compare with direct substitutions into
   the reduced source rows, using the inherited field/output units. Include
   multiple positions, widths and phases; a single scalar field sample does not
   establish the full operator action.
4. Evaluate every nonlocal contribution, retaining output/input/middle momentum
   roles and all measures. Keep zero-jet step profiles in coordinate space or
   use their already computed subtraction and delta/principal-value prescription.
   Test the assembly against direct reduced-row actions, preserve both operands,
   and use deliberate coefficient/measure mutations to check sensitivity.
   Shared finite quadrature agreement establishes implementation agreement;
   it does not by itself establish the infinite-domain or Abel weak limit.
5. Record grid/quadrature, spatial and momentum tails, regulator treatment,
   source-domain denominators and changes under refinement. Resolve the required
   physical limits before using the operator in two-ended boundary matching.
   Report unresolved domains explicitly. Do not substitute a fingerprint for an
   integral, or a finite-contrast result for the retained continuum expansion.

Use the accepted caches; no upstream regeneration is required for a numerical
adapter. Save operands before emission, preserve the report contract suffix,
publish validated transcripts through DataLad/git-annex and commit each step.
Use one heavy CAS job and silent local completion/error wake-ups.
