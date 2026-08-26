# Decision-list review — S11c-a PY shifted-trace fix directive

You are reviewing an orchestrator-written build directive BEFORE the builder runs. Your job is to find
defects in the DIRECTIVE and its embedded spec change — not to build anything. A directive defect costs a
whole build round plus its reviews, so probe hard; a review that finds nothing is weak evidence.

## Read first (source of truth), then the artifact
1. Spec (the source of truth): `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`
   — especially §1 (the drain `v_bulk_normal_0`), §2a (profiles, the `(ε,η,σ_W)` multigrade + truncation),
   §2d (background state `V_s⁰=J_s⁰=0`), §3b (`δp` background zero), §3c (the trace-linearisation law).
   Form your OWN view of what the §3c trace of a traced bulk field should contain, first.
2. The raw finding evidence (literal residual + two prior derivations):
   `~/.s11_build/S11c_a_shifted_trace_finding_measurement.md`. Treat the adjudication there as a CLAIM to
   check, not a given.
3. The artifact under review: `research/pde_ledger_v3/directives/S11c_a_py_shifted_trace_fix_directive.md`.

## What to check — answer each explicitly
A. Is the directive's "correct §3c reading" actually correct against the spec? Derive it yourself (a short
   sympy script is worth more than prose — expand `f(x,h_s) = f0 + ε·δf` with `h_s = sW_bg/2 + ε·δh`,
   `W_bg=W0(1+η w1)`, to first order in ε and η, and read the terms). State which terms survive for a field
   whose background is ZERO vs one whose background is NONZERO, and confirm/deny the directive.
B. DEFECT-SCOPE BOUNDARY (both directions):
   - Under-reach: are there traced-bulk objects that consume the bulk velocity or pressure that the
     directive's field list (RELATIVE_FLUX, TRACTION, EVOLUTION_MASS_BALANCE, KINEMATIC_BALANCE, pressure
     trace) MISSES? Name any T-c…T-i object that traces a bulk field and is not covered.
   - Over-reach: the directive keeps term-2 for "a field carrying a genuinely nonzero background profile"
     (it names the conormal probe field). Is that correct — does the spec give any traced field a nonzero
     background whose term-2 must be retained? If the directive would strip term-2 from a field that needs
     it (deleting a correct term), flag it.
C. ⛔ THE KILLER DEFECT: could the stated ACCEPTANCE pass WITH either defect still present? Walk the
   acceptance items against a PY output that still (i) carries a `d_w_v_bulk_0`/`d_w_delta_p_0` term, or
   (ii) still freezes term-1 at the flat face. If any acceptance item is satisfiable with the defect intact,
   say which and how. Is the acceptance value-free (rule 5) and able-to-fail?
D. LEAK: does the directive or the proposed §3c note state an expected output value, a sign/coefficient of
   the corrected result, or otherwise let a builder iterate toward a hidden target? (Supplying the physical
   premise "the background is zero" is allowed; stating what the corrected trace equals is not.)
E. The proposed §3c spec note (PART A): is it a STRUCTURAL premise consistent with §2d/§2a/§3b/§1, or does
   it assert a computed result or contradict anything already in the spec? Would inserting it change what a
   correct engine computes for any OTHER object (a regression risk)?

## Method
Physics filter: report a finding only if it changes what the fixed engine would compute or what the
directive would let pass. If you write a derivation script, save it and its literal stdout to named
absolute paths and report them — a prose derivation is discarded. Do not edit any file.
Report per-item A–E with a verdict and any defect, then an overall PASS / defect-list.
