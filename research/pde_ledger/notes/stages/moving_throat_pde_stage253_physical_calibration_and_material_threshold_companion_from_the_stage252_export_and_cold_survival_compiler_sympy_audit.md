# Moving-Throat PDE — Stage 253: Physical Calibration and Material-Threshold Companion from the Stage-252 Microscopic Export and Cold-Survival Compiler

## Status

**Exact within**

1. the Stage-252 channel-resolved microscopic export / cold-survival compiler,
2. the exact safe-edge lattice event-equivalent rate
   \[
   \gamma_{\rm lat,safe}^{\rm eq}
   =
   f_{\rm lat}(s_c)\,\mu_\eta\,\frac{s_0^2-s_c^2}{s_c},
   \]
3. the explicit physical calibration dictionary
   \[
   t^{\rm phys}=t_*\,t,
   \qquad
   r^{\rm phys}=\frac{\lambda_{\rm phys}}{\lambda_{\rm ref}}\,r,
   \qquad
   E^{\rm phys}=E_*\,E,
   \]
4. a one-scalar transport projection factor
   \(\Upsilon_{\rm lat}>0\) relating the Stage-252 event-equivalent lattice export rate to a coarse lattice-turnover rate,
5. the Session-V condensed-matter identifications
   \[
   \gamma_{\rm lat,turn}^{\rm phys}=\zeta_{\rm ep}\,\lambda_{\rm ep}\,\omega_D,
   \qquad
   V_{\rm lattice}(r)=\tfrac12 k_{\rm eff}r_{\rm phys}^2,
   \qquad
   T_1T=\mathcal K_{\rm corr},
   \]
6. and, for the benchmark slice, the same reduced turning-point and crossing-time data already carried by Sessions II–V.

This stage is a **physical-calibration / materials companion**. It is not part of the core PDE theorem ladder. Its job is to translate the exact reduced Stage-252 outputs into experimentally recognizable threshold equations while keeping the unit map and coarse-graining assumptions explicit.

---

## Purpose

Stage 252 removed the ad hoc vacuum-vs-lattice heat split and replaced it by the exact microscopic channel compiler
\[
E_{\rm vac}(T),\qquad E_{\rm lat}(T),\qquad f_{\rm vac}(s),\qquad f_{\rm lat}(s),
\]
together with the exact safe-edge event-equivalent rates
\[
\gamma_{\rm vac,safe}^{\rm eq},
\qquad
\gamma_{\rm lat,safe}^{\rm eq}.
\]
But Stage 252 stopped deliberately one step before condensed-matter language. It did **not** tell us how to express those microscopic rates and energies as

- an electron-phonon turnover threshold \(\lambda_{\rm ep}\omega_D\),
- a force-matched interstitial stiffness \(k_{\rm eff}\),
- or a Korringa-limited thermal spin-survival ceiling \(T_{\max}\).

That is exactly what Stage 253 does.

The main results are:

1. the lattice-turnover threshold is an exact algebraic image of the Stage-252 safe-edge rate once one introduces a single explicit coarse transport factor \(\Upsilon_{\rm lat}\),
2. the Session-V formulas are recovered as one particular calibration slice of that exact map,
3. the geometric trigger \(\chi_\lambda\) and the force-matched stiffness \(k_{\rm eff}\) separate cleanly,
4. the Korringa ceiling collapses to one exact reduced-time prefactor,
5. and the whole condensed-matter side reduces to a short list of material-screening inequalities.

So Stage 253 does not make the corridor stronger. It makes the Stage-252 microscopic compiler **physically calibratable and screenable**.

Script-backed status:
- `scripts/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.py`
  checks the symbolic calibration identities, the legacy Session-V recovery
  slice, the harmonic trigger/stiffness compiler, the Korringa ceiling, and the
  exact material-screening ratios, then evaluates the declared Stage-252 /
  Session-V benchmark readback.
- `mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl`
  mirrors the same split in a second CAS and keeps the explicit decimal values
  confined to the benchmark-only specialization layer.
- All decimal values appearing in the benchmark subsections below are either
  declared benchmark inputs or benchmark-derived outputs. The theorem path above
  them remains symbolic.

---

## Provenance

This stage is the direct continuation Stage 252 asked for.

- Stage 250 turned the relaxed same-charge branch into a crossing-vs-collapse problem.
- Stage 251 replaced the phenomenological \(\gamma_{\rm tot}\dot V\) envelope by the microscopic cubic-plus-quintic export kernel.
- Stage 252 derived the exact vacuum/lattice heat partition and the safe-edge exported-energy / event-rate formulas.
- Session V then translated the older phenomenological lattice thresholds into the condensed-matter variables \(\lambda_{\rm ep}\omega_D\), \(k_{\rm eff}\), and \(T_{\max}\).

Stage 253 is the missing bridge between the exact Stage-252 microscopic compiler and that Session-V materials language.

---

## 0. Why this stage is needed

Before this step the stack could already say:

- the minimum exported energy needed for cold survival,
- the exact vacuum-vs-lattice partition fraction on a chosen event,
- and the microscopic event-equivalent lattice export rate.

But it still could **not** say, without one more calibration move,

- what lattice-turnover threshold a real host must satisfy,
- what interstitial stiffness corresponds to the reduced barrier force,
- or how the spin-survival condition translates into a temperature ceiling.

There was also one important conceptual mismatch left open. Session V’s phenomenological lattice damping scalar
\[
\gamma_{\rm lattice}^{\rm red}\approx 4.79562976
\]
is **not** numerically the same object as the exact Stage-252 event-equivalent lattice rate
\[
\gamma_{\rm lat,safe}^{\rm eq}.
\]
So a companion stage has to do two things at once:

1. preserve the Session-V parameter map as a valid legacy slice,
2. and make explicit the extra calibration/projection factor needed to connect it to the exact microscopic Stage-252 object.

That is the role of \(\Upsilon_{\rm lat}\) below.

---

## 1. Physical calibration dictionary

Let the reduced-to-physical bridge be
\[
\boxed{
 t^{\rm phys}=t_*\,t,
 \qquad
 r^{\rm phys}=\frac{\lambda_{\rm phys}}{\lambda_{\rm ref}}\,r,
 \qquad
 E^{\rm phys}=E_*\,E.
}
\]
Here

- \(t_*\) is the physical time unit,
- \(\lambda_{\rm phys}\) is the physical localization width assigned to the chosen reduced branch,
- \(\lambda_{\rm ref}\) is the reduced width used to attach that physical scale,
- \(E_*\) is the physical energy scale for the reduced barrier potential.

On the Session-V turning-point slice, the natural choice is
\[
\lambda_{\rm ref}=\lambda_{\rm th}(r_{\rm turn}).
\]

From Stage 252, the exact safe-edge lattice event-equivalent rate is
\[
\boxed{
\gamma_{\rm lat,safe}^{\rm eq}
=
 f_{\rm lat}(s_c)\,\mu_\eta\,\frac{s_0^2-s_c^2}{s_c}.
}
\]
Introduce one explicit coarse transport / projection factor
\[
\boxed{\Upsilon_{\rm lat}>0}
\]
by defining the coarse physical lattice-turnover rate
\[
\boxed{
\gamma_{\rm lat,turn}^{\rm phys}
:=
\frac{\gamma_{\rm lat,safe}^{\rm eq}}{\Upsilon_{\rm lat}\,t_*}.
}
\]
This is the physical rate that Stage 253 identifies with the standard condensed-matter turnover lane
\[
\boxed{
\gamma_{\rm lat,turn}^{\rm phys}
=
\zeta_{\rm ep}\,\lambda_{\rm ep}\,\omega_D.
}
\]

The point of \(\Upsilon_{\rm lat}\) is precise:

- if the microscopic event-equivalent export rate is already the correct coarse lattice turnover, then \(\Upsilon_{\rm lat}=1\),
- if the older Session-IV/Session-V envelope rate is a coarser observable, then \(\Upsilon_{\rm lat}\neq 1\) captures that down-projection exactly.

This is a calibration factor, not a new dynamical degree of freedom.

---

## 2. Microscopic lattice-turnover compiler

With the dictionary above, the Stage-252 safe-edge rate implies the exact condensed-matter threshold
\[
\boxed{
(\lambda_{\rm ep}\omega_D)_{\min}^{(253)}
=
\frac{\gamma_{\rm lat,safe}^{\rm eq}}{\Upsilon_{\rm lat}\,\zeta_{\rm ep}\,t_*}
=
\frac{f_{\rm lat}(s_c)\,\mu_\eta\,(s_0^2-s_c^2)}{\Upsilon_{\rm lat}\,\zeta_{\rm ep}\,s_c\,t_*}.
}
\]
Using
\[
 t_{\rm cross}^{\rm phys}=\frac{t_*}{s_c},
\]
the same condition can be written as the exact event-product threshold
\[
\boxed{
(\lambda_{\rm ep}\omega_D\,t_{\rm cross}^{\rm phys})_{\min}^{(253)}
=
\frac{\gamma_{\rm lat,safe}^{\rm eq}}{\Upsilon_{\rm lat}\,\zeta_{\rm ep}\,s_c}
=
\frac{f_{\rm lat}(s_c)\,\mu_\eta\,(s_0^2-s_c^2)}{\Upsilon_{\rm lat}\,\zeta_{\rm ep}\,s_c^2}.
}
\]

So the lattice side is no longer a vague materials statement. It is one exact inequality with one explicit calibration factor.

### 2.1 Legacy Session-V slice

Session V used the phenomenological reduced lattice rate
\[
\gamma_{\rm lattice}^{\rm red}\approx 4.79562976.
\]
The exact Stage-253 map reproduces that older threshold if one chooses the calibration slice
\[
\boxed{
\Upsilon_{\rm lat}^{\rm(sess)}
:=
\frac{\gamma_{\rm lat,safe}^{\rm eq}}{\gamma_{\rm lattice}^{\rm red}}.
}
\]
Then
\[
\boxed{
(\lambda_{\rm ep}\omega_D)_{\min}^{\rm(sess)}
=
\frac{4.79562976}{\zeta_{\rm ep} t_*},
}
\]
\[
\boxed{
(\lambda_{\rm ep}\omega_D\,t_{\rm cross}^{\rm phys})_{\min}^{\rm(sess)}
=
\frac{8.73618521}{\zeta_{\rm ep}},
}
\]
which is exactly the Session-V turnover product.

So the Session-V parameter map is not discarded. It is recovered as a definite calibration slice of the exact Stage-252/253 bridge.

### 2.2 Benchmark specialization

Use the Stage-252 benchmark slice
\[
 f_{\rm lat}(s_c)=\frac34,
 \qquad
 \mu_\eta=1,
 \qquad
 s_c\approx 0.5489386551,
 \qquad
 s_0\approx 6.94311167.
\]
Then
\[
\boxed{
\gamma_{\rm lat,safe}^{\rm eq}\approx 65.45193926.
}
\]
On the raw microscopic slice \(\Upsilon_{\rm lat}=1\), this implies the much stronger direct product requirement
\[
\boxed{
(\lambda_{\rm ep}\omega_D\,t_{\rm cross}^{\rm phys})_{\min}^{\rm(micro)}
\approx
\frac{187.23361317}{\zeta_{\rm ep}}.
}
\]
The legacy Session-V slice is recovered by the benchmark factor
\[
\boxed{
\Upsilon_{\rm lat}^{\rm(sess)}
\approx
\frac{65.45193926}{4.79562976}
\approx
13.64824695.
}
\]
This is the cleanest way to state the relation between the exact microscopic Stage-252 rate and the older Session-V envelope rate.

---

## 3. Harmonic interstitial trigger and force-matched stiffness compiler

Take the turning point \(r_{\rm turn}\) on the reduced branch and attach the physical radius by
\[
\boxed{
 r_{\rm turn}^{\rm phys}
 =
 \frac{r_{\rm turn}}{\lambda_{\rm ref}}\,\lambda_{\rm phys}.
}
\]
For a harmonic interstitial trap
\[
\boxed{
V_{\rm lattice}(r)=\frac12 k_{\rm eff} r_{\rm phys}^2,
\qquad
\partial_r\ln V_{\rm lattice}=\frac{2}{r_{\rm phys}}.
}
\]
Therefore the geometric trigger ratio is exactly
\[
\boxed{
\chi_{\lambda,\rm lattice}(r_{\rm turn})
=
\frac{2\lambda_{\rm phys}}{r_{\rm turn}^{\rm phys}}
=
\frac{2\lambda_{\rm ref}}{r_{\rm turn}}.
}
\]
So the \(\chi_\lambda\) criterion fixes only a geometry ratio. It does **not** fix a stiffness by itself.

### 3.1 Force-matched stiffness

Let the reduced barrier potential near the turning point be \(V_{\rm red}(r)\). Under the physical energy map
\[
V_{\rm phys}(r_{\rm phys})=E_* V_{\rm red}(r),
\qquad
r=\frac{\lambda_{\rm ref}}{\lambda_{\rm phys}}r_{\rm phys},
\]
the physical barrier force at the turning point is
\[
\boxed{
F_{\rm barrier}^{\rm phys}(r_{\rm turn})
=
\frac{E_*\lambda_{\rm ref}}{\lambda_{\rm phys}}
\Bigl|V_{\rm red}'(r_{\rm turn})\Bigr|.
}
\]
Force matching against the harmonic trap
\[
F_{\rm lattice}^{\rm phys}(r_{\rm turn})=k_{\rm eff}\,r_{\rm turn}^{\rm phys}
\]
gives the exact stiffness compiler
\[
\boxed{
k_{\rm eff,req}
=
\frac{E_*\lambda_{\rm ref}^2}{\lambda_{\rm phys}^2}
\frac{|V_{\rm red}'(r_{\rm turn})|}{r_{\rm turn}}.
}
\]
Define the reduced force-matching coefficient
\[
\boxed{
\mathcal K_{\rm turn}
:=
\frac{\lambda_{\rm ref}^2}{r_{\rm turn}}|V_{\rm red}'(r_{\rm turn})|.
}
\]
Then
\[
\boxed{
k_{\rm eff,req}=\mathcal K_{\rm turn}\,\frac{E_*}{\lambda_{\rm phys}^2}.
}
\]
If one rewrites the localization width as half an interstitial scale,
\[
\lambda_{\rm phys}=\frac{a_{\rm int}}{2},
\]
then
\[
\boxed{
k_{\rm eff,req}=4\mathcal K_{\rm turn}\,\frac{E_*}{a_{\rm int}^2}.
}
\]

So the stiffness formula is an exact companion of the reduced barrier force **only after** one adds the force-matching assumption. It is not a pure consequence of the geometric \(\chi_\lambda\) trigger.

### 3.2 Session-V benchmark slice

On the turning-point benchmark carried by the session,
\[
 r_{\rm turn}\approx 0.39096144,
 \qquad
 \lambda_{\rm ref}=\lambda_{\rm th}(r_{\rm turn})\approx 0.42826825.
\]
Therefore
\[
\boxed{
 r_{\rm turn}^{\rm phys}
 \approx
 0.9128891530\,\lambda_{\rm phys},
}
\]
\[
\boxed{
 \chi_{\lambda,\rm lattice}(r_{\rm turn})
 \approx
 2.19084649.
}
\]
If one inserts the Session-V force-matching benchmark
\[
\boxed{\mathcal K_{\rm turn}\approx 2.73855812,}
\]
then the stiffness formulas become
\[
\boxed{
 k_{\rm eff,req}
 =
 2.73855812\,
 \frac{E_*[\mathrm{eV}]}{\lambda_{\rm phys}^2[\mathrm{\AA}^2]}
 \quad \mathrm{eV/\AA^2},
}
\]
\[
\boxed{
 k_{\rm eff,req}
 =
 10.95423247\,
 \frac{E_*[\mathrm{eV}]}{a_{\rm int}^2[\mathrm{\AA}^2]}
 \quad \mathrm{eV/\AA^2},
}
\]
which is exactly the Session-V reported stiffness map.

---

## 4. Korringa-limited thermal spin-survival compiler

Use the standard Korringa relation
\[
\boxed{T_1 T = \mathcal K_{\rm corr}.}
\]
Impose the spin-survival condition on the same cold branch used in Sessions II–V,
\[
\boxed{T_1\ge t_{\rm cross}^{\rm phys}.}
\]
Then the exact ceiling temperature is
\[
\boxed{
T_{\max}
=
\frac{\mathcal K_{\rm corr}}{t_{\rm cross}^{\rm phys}}
=
\frac{s_c\,\mathcal K_{\rm corr}}{t_*}.
}
\]
Using the reduced crossing-time benchmark
\[
 s_c\approx 0.5489386551,
\]
one gets
\[
\boxed{
T_{\max}
\approx
0.548938655\,\frac{\mathcal K_{\rm corr}}{t_*}.
}
\]
So the Korringa side is the cleanest of the three condensed-matter companions: once \(t_*\) and the material-specific \(\mathcal K_{\rm corr}\) are fixed, the ceiling temperature follows immediately.

---

## 5. Material-screening companion invariants

Stage 253 packages the Session-V condensed-matter map into four exact screening ratios.

### 5.1 Lattice-turnover ratio
\[
\boxed{
\Pi_{\rm ep}
:=
\frac{\Upsilon_{\rm lat}\,\zeta_{\rm ep}\,\lambda_{\rm ep}\,\omega_D\,t_*}{\gamma_{\rm lat,safe}^{\rm eq}}
=
\frac{\lambda_{\rm ep}\omega_D}{(\lambda_{\rm ep}\omega_D)_{\min}^{(253)}}.
}
\]
The cold-survival / lattice-turnover condition is
\[
\boxed{\Pi_{\rm ep}\ge 1.}
\]

### 5.2 Geometric trigger ratio
\[
\boxed{
\Pi_{\chi}
:=
\chi_{\lambda,\rm lattice}(r_{\rm turn})
=
\frac{2\lambda_{\rm ref}}{r_{\rm turn}}.
}
\]
The formal geometry-side trigger is
\[
\boxed{\Pi_{\chi}\ge 1.}
\]

### 5.3 Force-matched stiffness ratio
\[
\boxed{
\Pi_k
:=
\frac{k_{\rm eff}\,\lambda_{\rm phys}^2}{\mathcal K_{\rm turn}E_*}
=
\frac{k_{\rm eff}}{k_{\rm eff,req}}.
}
\]
The stiffness requirement is
\[
\boxed{\Pi_k\ge 1.}
\]

### 5.4 Thermal spin-survival ratio
\[
\boxed{
\Pi_T(T)
:=
\frac{\mathcal K_{\rm corr}}{T\,t_{\rm cross}^{\rm phys}}
=
\frac{T_{\max}}{T}.
}
\]
The spin-survival requirement is
\[
\boxed{\Pi_T(T)\ge 1.}
\]

So a candidate host survives the Stage-253 companion screen only if
\[
\boxed{
\Pi_{\rm ep}\ge 1,
\qquad
\Pi_{\chi}\ge 1,
\qquad
\Pi_k\ge 1,
\qquad
\Pi_T(T)\ge 1.
}
\]

That is the exact condensed-matter triage package carried by this stage.

---

## 6. What this stage achieves physically

Stage 253 changes the barrier thread in four concrete ways.

### 6.1 It preserves the Session-V map while making the Stage-252 microscopic object explicit

The lattice-turnover threshold is no longer forced to pretend that the Stage-252 event-equivalent export rate and the older Session-IV/Session-V envelope scalar are the same numerical object. Their relation is now carried by one explicit factor \(\Upsilon_{\rm lat}\).

### 6.2 It separates geometric trigger from stiffness

The \(\chi_\lambda\) condition only says whether the turning-point geometry sits inside the formal steep-trap regime. The explicit stiffness formula needs one more assumption: force matching to the reduced barrier force. Stage 253 makes that separation exact.

### 6.3 It turns thermal spin survival into one clean ceiling formula

Once \(t_*\) and \(\mathcal K_{\rm corr}\) are fixed, the aligned-spin survival condition becomes the exact ceiling
\(T_{\max}=\mathcal K_{\rm corr}/t_{\rm cross}^{\rm phys}\).

### 6.4 It packages the materials problem as a short inequality stack

Instead of vague questions about whether a host is “good,” the condensed-matter side is now four explicit screening ratios: turnover, geometric trigger, stiffness, and Korringa survival.

---

## 7. What remains open

This stage is still a companion, not a materials theorem.

The following remain open:

1. the physical unit map \((t_*,\lambda_{\rm phys},E_*)\),
2. the transport projection factor \(\Upsilon_{\rm lat}\) unless one simply adopts the legacy Session-V slice,
3. the material-dependent constants \(\zeta_{\rm ep}\) and \(\mathcal K_{\rm corr}\),
4. and actual candidate-host screening against real databases.

So the correct reading is:

- the barrier/stability chain now produces exact reduced condensed-matter threshold equations,
- but the project still has to calibrate units and compare those equations against real hosts before making any material claim.

That is why this stage is correctly a **companion** to the core PDE ladder rather than a theorem stage of the ladder itself.

---

## 8. Immediate next step

The next honest move is now sharply defined.

1. Fix or constrain \(t_*\), \(\lambda_{\rm phys}\), and \(E_*\).
2. Decide whether the materials map is being carried on the raw microscopic slice \((\Upsilon_{\rm lat}=1)\) or on the legacy Session-V envelope slice \((\Upsilon_{\rm lat}=\Upsilon_{\rm lat}^{\rm(sess)})\).
3. Insert candidate-host values of \(\lambda_{\rm ep}\omega_D\), effective interstitial stiffness, and \(\mathcal K_{\rm corr}\).
4. Apply the exact screening ratios \((\Pi_{\rm ep},\Pi_\chi,\Pi_k,\Pi_T)\).

That is the smallest falsifier left on the condensed-matter side before returning to the actual moving-throat branch realization problem.
