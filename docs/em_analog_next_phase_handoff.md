# HANDOFF — EM-analog research phase: define the `±w`-body foundations, then compute

**Read this first after `/compact`.** It resumes the research phase agreed 2026-07-12. Companion READ-FIRST: **`docs/model_definition_audit.md`** (the honest map of what's DEFINED vs IMPORTED vs OPEN; §F = the to-define list). Do not re-read the whole `conceptual_foundation.md` — it carries vision+history that repeatedly re-confused us; the audit supersedes it for "what do we actually have."

## Where we are (one paragraph)
The model is a single 4+1D compressible superfluid; our 3D space is a brane (sheet) in it; particles are throats (punctures). **Gravity (brane drain) and light (brane MacCullagh shear) are genuinely built** (calibrated GR-match; 2 transverse photons on a postulated stiffness). **EM is a *characterized-departure analog*, not exact Maxwell** — and a just-completed, dual-engine, unanimously-verified Dirac gate proved **exact/first-class `U(1)` Maxwell is NOT native** (`NATIVE_P_NO_EMERGENT_GAUSS`; the continuum little-arrows sector is second-class). The whole EM reconsideration was triggered by a concrete anomaly: the `force_visualizer` sim shows two parallel currents (like swirls) attracting, and an investigation found **that sign was imported (`pathA_39` asserted `j=sV`), not derived.** A two-reviewer fact-check (Codex+Grok) confirmed the audit is honest.

## The decision (banked)
**Target = the characterized-departure Maxwell analog.** NOT exact Maxwell (non-native by computation). **Method: define each imported/under-defined piece the model's OWN way, compute, and accept whatever departure from Maxwell falls out** — do not import Maxwell structure, do not chase an exact match. (This is the program's core philosophy: analog not derivation; characterized departures are first-class; a clean exact match would be *suspicious*.)

## The research phase (what to actually do) — see audit §F for the full list
The Maxwell analog rests on **3 imported pieces downstream of 2 definable foundations + 1 structural fact:**
- **Imported (become computed once the foundations are defined):** I1 the current law `j=sV`; I2 the electric sign "like-repel" (rests on assumed lock `G₀>0`; no-lock limit attracts); I3 the Maxwell exchange machinery + dual sign (two independently-signed kernels, electric on `h` / magnetism on `u_T`).
- **⭐ U1 (keystone) — DEFINE the `±w` throat-body as a dynamical object:** structure (puncture extended into `w`), effective inertia (mass from displaced medium), force-response.
- **U2 — DEFINE its medium-coupling** (roll vs slip: does a moving puncture drag medium or slip through it?), from the puncture's actual drain/vortex-like structure.
- **U3 — structural fact:** charge is natively `Z₂` (two signs); additive-integer charge/Gauss/conservation are computed non-native. Not upgradable natively.

**Order:** U1 → U2 → compute the **static** two-body force (→ native electric sign, replaces I2) and the **moving** two-body force (→ native magnetism sign + whether the current is really `∝V`, replaces I1/I3). **U1 is the immediate next step.**

## Why this is the right frame (don't relapse)
- We kept **importing standard-EM** whenever a piece was under-defined (the current law, the signs, even "charge is conserved"). The point of this phase is to STOP importing and define the model's own objects, then read off the physics.
- Don't re-open "can we get exact/first-class Maxwell" — that's computed non-native (native-`P` gate, verified). Don't re-raise "is the swirl a circulation vs a puncture-body" — the ontology already commits (moving 4D throat-body, NOT bulk vorticity=gravitomagnetism); the gap is the *dynamics* of that committed object, which is U1.
- "Photon count / exact 2-polarization gauge boson" is a **piggyback** — in this model light (the photon) is brane shear; the EM sector supplies *forces*, not its own photon. Don't let a reviewer kill an EM idea on photon-count grounds.

## Process / tooling (banked lessons)
- **The RIG lesson (critical):** a build CAN fake a computed verdict (hardcoded result, faked controls, tautological "search") and still pass `ENGINE_AGREE` + controls + arbiter-rerun — all theater. **Any computed verdict MUST be verified by a FRESH agent that imports the engine and independently RECOMPUTES the decisive quantity** (not just reads code). A negative verdict matching our prior needs the hardest check.
- **Gauntlet for any computation:** directive (requirement+acceptance, Codex designs the route) → Codex design-review to `DIRECTIVE_CLEAN` → build dual-engine (SymPy + independent Mathematica) with anti-fake guards (couplings must genuinely enter the brackets; controls share the code path; able to return any outcome) → fresh-agent adversarial verification → arbiter re-run.
- **Tooling:** `codex exec -c model_reasoning_effort=xhigh "…" < /dev/null` (stdin hang without `< /dev/null`); for Mathematica-running builds add `--sandbox danger-full-access`; **Grok works reading-only** (`grok --prompt-file … --output-format plain`) but its script-execution gets killed — use it for doc/code *review*, not compute; GLM via `opencode run --agent plan -m cloudflare-workers-ai/@cf/zai-org/glm-5.2`. Panels = Codex+Grok(+GLM). Background codex/grok as `run_in_background:true` (NEVER wrap codex in `&` or shell `timeout`).
- **Key files:** `docs/model_definition_audit.md` (the map), `docs/em_bulk_lattice_route.md` (the killed bulk-lattice route + panel), `software/em_charge_attribute/` (the native-`P` gate code + the vortex-sign directive `directive_native_vortex_sign.md`, which was correctly *not built* because it revealed U1/U2 are undefined).

## Immediate next step
Start **U1**: a conceptual + computed definition of the `±w` throat-body as a moving object in the medium — its structure, effective inertia, and force-response — in a way that makes sense to the model (a puncture/drain-defect in a superfluid), NOT by importing a point-charge or a Maxwell current. Talk it through conceptually with the user first (what IS the body, physically), then scope the computation and run the gauntlet.
