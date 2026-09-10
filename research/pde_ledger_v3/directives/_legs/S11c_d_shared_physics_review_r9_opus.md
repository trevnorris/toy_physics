# Fresh Claude agent (Opus) review — S11c-d SHARED PHYSICS spec v9 (round 9, DOCUMENT branch)

VERDICT: NOT-SOUND (1 must-fix, CONVERGENT with Grok F1; 1 nit). Faithful summary; full JSONL = out-of-repo task file.
Mechanically verified against the real rows: both closed rows integrate d3y (Y1-3) AND d3k (k_out/k_in/Middle),
all (2pi)-bearing, 0 DiracDelta; transfer + JetHat{1,2,3} hats at BOTH transfer and middle-leg args -> the §1c
universal-principle enumeration is exhaustive + delta-free (F2 met). Order/strong-edge re-verified symbolically.

MUST-FIX (= Grok F1, convergent): §3a J current "evaluated on the imported closed operator" + §3b pole resolvent
"of the retained first-shape-order imported operator" are the residual F1 SEAM -- 𝓛 is bound to the reduced rows
(§2/§7) but its current/flux/pole machinery can still ingest UNREDUCED 3-D content, a mixed 3-D/1-D assembly that
also breaks §3c d2(Q_par)-stripping. FIX (Opus exact global form): qualify §3a/§3b as the §1c-REDUCED operator,
or add one sentence to §2/§6: "Wherever §3 says the imported closed operator it denotes its §1c-reduced form; no
downstream current, resolvent, mode, or amplitude is built from unreduced 3-D content." Change only the Fourier
representative, keep the S11b current FORMULA + closed/nonlocal bulk content.
NIT 1: §2 -- reduced kernel is "canonical" off-diagonal; a nonzero reduced-operator-block vs reduced-kernel
residual should be a surfaced finding, not silently overridden.
All other axes (order, strong-edge, c2 honesty, N6, two channels, M2/M3, N11, completeness): sound, no regression.
