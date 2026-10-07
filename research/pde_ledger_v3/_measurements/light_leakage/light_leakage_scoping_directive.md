# Light-leakage scoping: fiber-loss audit and bracket inventory (directive for Codex)

**Author:** Claude (orchestrator), 2026-10-06. `AGENTS.md` applies. Folded once after one Codex + Grok pass.
Paths are relative to `research/pde_ledger_v3/` unless stated.

## Context

Light in v3 is a transverse wave of the brane (`steps/S9_light_requires_shear.md:11`). The bulk carries no
shear (`:64`), which is postulated (`:337`).

S11c asked how much light converts into other motion where the brane changes. It obtained "no accepted
nonuniform light-loss number", and in its accounting "both reflected and transmitted light count as survival"
(`steps/S11c_PARTIAL_CLOSEOUT.md:3`). Conversion routes exist in the equations, but "their size remains
unknown", and the benchmark's residual is "not a physical bound" (`:9`).

A later step will **bracket** the loss by estimating each route under worst-case and best-case parameters and
comparing with observational bounds. This task builds the inventory that bracket needs. **It computes nothing
and gives no verdicts.** Use the closeout's survival rule (`:3`) as the accounting rule throughout.

## Part 1. Fiber-loss mechanism audit

**Scope.** Cover:
- the loss mechanisms of guided light named in the loss chapters of two standard references you cite (an
  optical-fiber text and a waveguide-theory text);
- these, if the references miss them: absorption; Rayleigh, Raman and Brillouin scattering; macrobending and
  microbending; interface roughness and diameter fluctuation; tapers; leaky modes and tunneling; mode and
  polarization coupling.

Stop there.

**For each mechanism, fill independent columns:**
1. **Fiber mechanism.** The mechanism and the physical condition it needs, with a source. Mark any source you
   could not open **UNVERIFIED**.
2. **Model mapping.** One of:
   - **ABSENT**: the condition does not exist in the model as the records state it;
   - **PRESENT ANALOG**: say how it differs from the fiber case, in the records' terms, with record lines;
   - **NOT ADDRESSED** by the records.
3. **Record status of the analog,** as the records mark it: ESTABLISHED, CONDITIONAL, OPEN or UNRESOLVED.
4. **Observational constraint.** The observation, the quantity it bounds, and its source. Do not say whether the
   model must suppress the mechanism.
5. **Destination of the energy.** One of: the same transverse branch (survives, per the closeout's rule);
   another brane branch; the bulk; or an exchange of frequency.

## Part 2. Bracket inventory

1. **Routes.** Give each item a status: **equation-level route**, **candidate beyond scope**, or **open
   requirement / non-establishment**.
   - Seed the equation-level routes from:
     - `steps/S11c_b_variable_coefficient_operator.md:54–57` (off-diagonal coupling kernel);
     - `steps/S11c_c1_curved_bulk_closure.md:35–44` (bulk impedance operator, radiation regimes, dissipation
       audit);
     - `steps/S11c_c2_self_energy_fold.md:63–70` (closed coupling kernel).
   - Also include:
     - the kinematic bulk channel, which is "kinematic only" (`steps/S11_stray_longitudinal.md:97–101`);
     - A's frequency-dependent leak and the gradient-driven thickness channel
       (`steps/S11bB_interface_assembly.md:196–202`).
   - Put `steps/S11_stray_longitudinal.md:74–76`, `steps/S11c_d_profile_conditioned_scattering.md:64` and
     `R-S8-02` under the third status.
   - For each item, give the record lines, and quote what the records say about the drain `v₀` and any other
     freeze the item depends on. Starting points: `steps/S11bB_interface_assembly.md:46–48` and `:150–154`,
     and `SUBSTRATE_REQUIREMENTS.md:207–212`. Name a freeze as a freeze, and do not resolve the rest-versus-driven
     tension.
2. **Oracle candidates, never bounds.** For each route, list the prior-art results that a future
   independent v3 calculation could be checked against (`CLAUDE.md` M3). For each, give:
   - the exact equation, as the source states it;
   - its assumptions and its domain of validity;
   - the correspondence between its parameters and the records' objects.

   Never call one a v3 bound. Starting sources, all unverified:
   - Marcuse, *Bell Syst. Tech. J.* (1969), on surface imperfections;
   - Marcuse, *Bell Syst. Tech. J.* (1970), on tapers. Give exact titles and DOIs for both papers;
   - Molz & Beamish, *JASA* 99:1894 (1996);
   - Demma, Cawley & Lowe, *JASA* 113:1880 (2003);
   - Kubrusly, von der Weid & Dixon, *NDT&E* 108 (2019);
   - Peyton et al., *NDT&E* 158 (2026);
   - Gu & Fuller, *JASA* 90:2020 (1991).
3. **Route-by-regime matrix.** One row per route and regime, with these columns:
   - source and destination branch;
   - the local quantity (per encounter or per length), as the records define it, or "gap";
   - the **exposure** variable (path length, encounter count or column), and its record status;
   - frequency and polarization dependence, as recorded;
   - domain of validity;
   - parameter correlations the records note;
   - survival convention;
   - the applicable yardstick from item 4, or "gap".

   Give parameter statuses as `SUBSTRATE_REQUIREMENTS.md:543–570` gives them: input assignments, not physical
   values.
4. **Yardsticks.** For each, give:
   - the bounded quantity;
   - the bound as published, with its source;
   - the destination it constrains;
   - the exposure it integrates over.

   Include these record-named yardsticks:
   - `v₀` "bounded by cosmology" (`steps/S11bB_interface_assembly.md:62–63`);
   - the Debye loss peak, bounded by measurements of anomalous mechanical dissipation, where a numerical bound
     needs a gravity-sector coupling that is "not built" (`:65–68`);
   - the thickness channel, "bounded by BENCH-TOP OPTICS" (`:200–202`).

   Also include:
   - stellar energy budgets. State what Friedland & Giannotti, *PRL* 100:031602 (2008), actually bounds;
   - cosmological transparency;
   - the observations that exclude a gradual per-photon energy loss (tired light).

## Rules

- No calculations, no new physics, and no verdicts. Do not say what the bracket will find, or whether the model
  must suppress anything.
- Every claim about the model cites a record line. Every prior-art claim cites a source. Mark unverified
  sources.
- Quote disagreements between records; do not resolve them.

## Output

- Write `_scratch/light_leakage/light_leakage_scoping.md` (repository root path). Do not modify any other file,
  and do not commit.
- STOP report, at most half a page:
  - the number of mechanisms per model-mapping value (Part 1);
  - the number of routes per status (Part 2);
  - the unverified sources;
  - anything you could not do.
