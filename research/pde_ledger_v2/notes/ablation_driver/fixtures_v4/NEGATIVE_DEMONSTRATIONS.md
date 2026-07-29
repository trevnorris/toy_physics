# Able-to-fail execution record

The A1--A9 discriminators were executed before the public driver existed. For
each item, an observation obtained from a disposable probe was accepted, its
one-property-deficient twin was rejected, and restoring only that property made
the same twin acceptable. All probes and their working trees were removed.

The same three-step execution was repeated after each repair round. The
exact-replacement, invocation-copy, fixture-owned entry-snapshot, invalid-list,
full result-audit, earned-progress, direct-execution, closed-environment, caller
tree, and live legacy-join discriminators were exercised the same way.
Disposable supervision probes completed with audited row progress visible and
restored their mirror after both normal completion and contained interruption.

| Item | Deficiency shown caught | Baseline | Deficient twin | Property restored |
|---|---|---|---|---|
| A1 | residual artifact or reset state reaches the next row | accepted | rejected | accepted |
| A2 | zero and multiple application failures collapse or execute | accepted | rejected | accepted |
| A3 | resume re-runs a banked row or differs from uninterrupted output | accepted | rejected | accepted |
| A4 | target or entry artifact is not restored on a required path | accepted | rejected | accepted |
| A5 | captures are shared, overwritten, or assigned to the wrong row | accepted | rejected | accepted |
| A6 | validation accepts without losslessly parsing the supplied list | accepted | rejected | accepted |
| A7 | a mapped retrofit field differs from the committed legacy result | accepted | rejected | accepted |
| A8 | committed evidence makes a claim unsupported by committed files | accepted | rejected | accepted |
| A9 | a row outcome overclaims what its own fields establish | accepted | rejected | accepted |

The final repair probes were also executed in disposable directories:

- Suite-local bytecode made the freeze gate reject before the oracle import.
- A generic validator accepted a runtime-generated transformed list; a canned
  digest response was rejected. A validator with no rejection branch was then
  rejected when it reported success for the invalid-list trial.
- The exact mutant observed by the producer was accepted; an extra-byte mutant
  carrying the reported exact digest was rejected. The same C-4/C-5 audit was
  invoked from every acceptance item; removing that invocation from any one item
  made the audit-call coverage probe fail.
- The caller-tree snapshot rejected content changes under `.git` as well as
  caller content, path-set, and type changes elsewhere.
- Identity-correct forged rows with inapplicable captures did not count as
  progress. Rows advanced the watchdog only after their complete result fields
  and referenced capture bytes passed the shared audit.
- Direct argv with the exact configured environment was accepted. A
  `shlex.join` plus shell-execution variant and a variant that leaked only a
  selected ambient whitelist were rejected.
- A successful A7 run whose recovery state was removed passed without a repair
  invocation; an uncatchably terminated A4 run still required and passed repair.
- Changing a mapped candidate field while recomputing its projection digest was
  rejected by the live C-9 join against the pinned legacy tables.
- The named inactivity, fresh-mirror, and post-wait A4 predicates were shown to
  be entailed by their own harness and deleted. The sweep removed six more
  construction-only or already-entailed assertions. A no-progress supervision
  probe then failed through the actual watchdog path.
- A documentation assertion confirmed the stated self-verifying decision
  function boundary, execution-environment substitution boundary,
  adoptable-oracle separation, and accepted gaps.

The execution transcript is recorded in the fixture-authoring session report;
no expected driver output is retained here.
