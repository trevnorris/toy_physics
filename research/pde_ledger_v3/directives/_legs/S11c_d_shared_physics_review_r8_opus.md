# Fresh Claude agent (Opus) review — S11c-d SHARED PHYSICS spec v8 (round 8, DOCUMENT branch)

VERDICT: SOUND (2 optional clarity nits). Faithful summary; full JSONL transcript is the out-of-repo task file.
Independently measured both closed rows carry the same 3-D Fourier content (operator 88 transfer + 300 jet-hat
+ 1020 middle + 218 integrals; kernel 550/408/1452/526; both 0 deltas). Confirmed all 4 focus questions YES:
the control now reduces BOTH operator AND kernel (all hat families + middle-leg args + in-plane integrals) with
one per-engine convention so 𝓛 cannot be assembled from mixed conventions; sec2 uses both reduced rows; sec7
joins both; deferral leak-safe. Order bookkeeping + strong-edge re-verified symbolically; no regression.
NIT 1 (Opus): sec1c names in-plane POSITION integrals but the rows also carry middle-momentum CONVOLUTION
integrals (d3k, (2pi)^3) -- Opus judged "already covered by the deferral". NIT 2: sec2 "y" notation overload.

ADJUDICATION NOTE (orchestrator): the Grok round-8 leg escalated Opus NIT 1 to a MUST-FIX (F2) -- correctly:
the CONTROL list is the builder checklist; naming only POSITION integrals lets d3k momentum measures stay 3-D.
Grok also found F1 (𝓛 not BOUND to the reduced rows; sec2 says "full IMPORTED operator"). Both G4-verified.
=> round-8 GATE verdict NOT-SOUND on F1+F2; both legs endorse the deferral APPROACH.
