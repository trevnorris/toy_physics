# Stage V2-29 - Superfluid Material Closure Gap

## Purpose

This note records a material-closure gap exposed by the PDE audit discussion.
It should be carried into the paper draft because it affects whether the final
PDE can be claimed as complete.

The short version:

```text
The notes contain partial superfluid EOS and sound-speed formulas, but the
audited branch machinery does not yet solve a self-consistent superfluid
material sector that determines rho, c_s(rho), and any density-dependent
effective light speed on the same frozen moving-throat branch.
```

This is not a small bookkeeping issue.  If the observable speed of light, the
tail transport speed, or the Maxwell localization speed depends on the
superfluid density, then density is a load-bearing branch variable rather than
a hidden constant.

---

## 1. What is already present

The project has not completely ignored the superfluid equation of state.
Several notes carry a frozen stiff-polytropic closure of the schematic form

```text
P(rho) = K_EOS rho^n
h(rho) = n K_EOS rho^(n-1)/(n-1)
c_s^2(rho) = (1/m) dP/drho.
```

For the commonly carried `n = 5` branch this gives

```text
c_s^2(rho) = (5 K_EOS/m) rho^4.
```

The notes also contain background-density formulas in reduced Thomas-Fermi or
stationary-branch settings, and some 2PN-side material uses normalized
sound-speed functions such as `C_s(u)=c_s^2(u)/c_s0^2`.

So the right statement is not:

```text
There is no EOS anywhere in the project.
```

The right statement is:

```text
The EOS appears as a reduced or frozen constitutive ingredient, but it is not
yet closed into the audited moving-throat branch exporter.
```

---

## 2. What is missing

The current audit and simulation packets mostly treat

```text
rho0, c_s, c, K_EOS, Z(w), W(w)
```

as symbolic constants, fixture inputs, normalization choices, or frozen branch
data.  That is acceptable for checking conditional algebra, but it is not a
completed physical PDE.

A completed branch realization needs equations that determine at least:

1. The superfluid density field:

   ```text
   rho = rho(X,t)
   ```

   including its stationary throat profile and perturbations.

2. The superfluid phase/current sector:

   ```text
   j_rho, v_s, continuity, intake, leakage, and outgoing/export terms.
   ```

3. The EOS or internal energy:

   ```text
   U(rho), P(rho), h(rho), chemical potential.
   ```

4. The local sound speed:

   ```text
   c_s^2(rho) = (1/m) dP/drho
   ```

   or its corrected equivalent if the final EOS is not the frozen
   stiff-polytropic one.

5. The effective light-speed relation, if the model assumes one:

   ```text
   c_eff = c_eff(rho, localization data, Maxwell/mixed branch data).
   ```

6. The way these material fields feed the audit readouts:

   ```text
   D0, A, C, P0, N0, N2, N4.
   ```

7. The way they feed the tail gate:

   ```text
   Theta_tail (c/c_s)^3 = 1
   ```

   or the branch-derived replacement for that gate.

Without those equations, `c_s` and `c` are parameters placed on the branch, not
outputs of the branch.

---

## 3. Why this can block the final PDE

This gap can prevent completion of the one-PDE program.

If density controls the propagation speeds, then a branch that matches the
target coefficients must also satisfy material constraints.  It is not enough
to find values of

```text
D0, C, P0, N2, N4
```

that pass the reduced audit.  The same branch must explain why the required
values of `c`, `c_s`, and `rho` are obtained and why they remain stable enough
to match observed physics.

In particular:

- If `c_eff` depends on `rho`, then density variations would generically change
  the effective light cone unless the branch has a stabilizing or screening
  mechanism.
- If `c_s` depends on `rho`, then the tail-transport gate is also a material
  equation, not just a scalar normalization condition.
- If `rho0` is left free, then source strength, port normalization, throat
  impedance, and outgoing transfer can hide untracked tuning.
- If the throat intake/output changes `rho`, then the open-system flux ledger
  must feed back into the coefficients rather than being appended after the
  coefficient extraction.

This gives a plausible reason why the reduced algebra can be clean while the
current simulations miss the target: the simulations are not solving the
material sector that would select or rule out the needed branch.

---

## 4. Referee-safe claim language

The paper should avoid saying:

```text
The present PDE derives the speed of light and sound speed from the superfluid
density.
```

unless the material sector is actually solved.

A safer statement is:

```text
The current audit treats c, c_s, and rho0 as frozen branch data or reduced
constitutive inputs.  Existing notes contain candidate EOS relations, including
a stiff-polytropic sound-speed law, but the moving-throat exporter has not yet
closed the superfluid material sector that would derive those quantities on the
same branch that produces the audit coefficients.
```

The strongest positive statement presently supported is:

```text
The audit identifies the additional material-closure equations that a completed
PDE must supply before the reduced target algebra can be promoted to a full
physical derivation.
```

---

## 5. Recommended next work package

The next derivation package should be a dedicated superfluid material closure.
It should produce a frozen branch packet only after the following are fixed
before target comparison:

```text
EOS choice or derivation
rho0(X) stationary throat profile
c_s(rho0) and perturbative c_s response
c_eff(rho0) or proof that c is density-insensitive in the relevant sector
intake/output continuity law
source and port normalization from rho and current data
tail-gate transport factor from branch data
effect of material response on D0, A, C, P0, N0, N2, N4
```

The important no-refit rule still applies: this material packet must be frozen
before comparing to the GR-like target surface.  Otherwise the material sector
would become a new tuning layer rather than a physical branch derivation.

---

## 6. Practical interpretation of the simulation miss

The current simulation miss should not be interpreted as proving the whole
program false.  It should also not be ignored.

The most accurate interpretation is:

```text
The reduced/manufactured branch families audited so far do not contain the
needed target branch.  One major missing ingredient is a self-consistent
superfluid material sector that determines density, sound speed, effective
light-speed behavior, and open-system flux feedback together with the throat
response coefficients.
```

This is exactly the kind of gap that can make a final PDE substantially harder
than the reduced algebra suggested.
