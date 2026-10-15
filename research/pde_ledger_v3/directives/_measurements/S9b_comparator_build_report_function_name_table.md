# Function-name table

Native aliases are grouped in a single injective row. Meanings describe mathematical constructors; physical names and source citations are listed below. Repoint controls cover every row and alias.

| SymPy constructor | Wolfram aliases | Function denoted |
|---|---|---|
| Add | Plus | sum |
| Mul | Times | product |
| Pow | Power, densityPower | power, including the source-defined zero exponent convention |
| Equality | Equal | equality |
| Unequality | Unequal | inequality |
| StrictGreaterThan | Greater | strict greater-than |
| StrictLessThan | Less | strict less-than |
| GreaterThan | GreaterEqual | greater-than or equal |
| LessThan | LessEqual | less-than or equal |
| Sequence | List | ordered tuple |
| Integral | Integrate, Inactive/Integrate | held integral with scoped limits |
| Limit | Inactive/Limit | held limit |
| Contains | Element | set membership |
| ForAll | ForAll, Inactive/ForAll | universal quantifier over a scoped real binder |
| Exists | Exists, Inactive/Exists | existential quantifier over a scoped real binder |
| And | And | logical conjunction |
| Or | Or | logical disjunction |
| Not | Not | logical negation |
| Piecewise | Piecewise | ordered conditional branches |
| Rational | Rational | exact rational number |
| Min | Min | minimum |
| Max | Max | maximum |
| Abs | Abs | absolute value |
| sin | Sin | sine |
| cos | Cos | cosine |
| tan | Tan | tangent |
| sec | Sec | secant |
| acos | ArcCos | principal inverse cosine |
| log | Log | principal logarithm |
| exp | Exp | exponential |
| sqrt | Sqrt | principal square root |
| re | Re | real part |
| im | Im | imaginary part |

## Physical-name bindings

The following unchanged table is transcribed from the comparator source. Spec line citations refer to `S9b_SHARED_PHYSICS.md`, except entries explicitly citing the repair directive. Engine lines refer to `scripts/S9b_light_bending_sympy_audit.py` and `mathematica/S9b_light_bending_mathematica_audit.wl`. The two velocity integration constants are deliberately engine-qualified and absent from this table.

| SymPy name | Wolfram name | Spec object / lines | SymPy line | Wolfram line |
|---|---|---|---:|---:|
| `s9b_r` | `r` | radial coordinate; 141,175 | 668 | 45 |
| `b` | `b` | flyby impact parameter; 177 | 668 | 45 |
| `Z_E` | `ZE` | emitter endpoint; 178-180 | 668 | 45 |
| `Z_R` | `ZR` | reflector endpoint; 178-180 | 668 | 45 |
| `r_E` | `rE` | emitter radius; 180 | 669 | 69 |
| `r_R` | `rR` | reflector radius; 180 | 669 | 69 |
| `s9b_b_min` | `bFar` | far-zone lower bound; 172-176 | 669 | 45 |
| `GM` | `GM` | independent orbital mass parameter; 202-204 | 670 | 46 |
| `c_0` | `c0` | asymptotic optical speed; 116 | 668 | 45 |
| `K` | `K` | bulk pressure coefficient; 219-220 | 670 | 46 |
| `n` | `n` | bulk pressure exponent; 219-220 | 670 | 46 |
| `m` | `m` | bulk particle mass; 219-220 | 671 | 46 |
| `s` | `s` | optical power exponent; 246 | 670 | 46 |
| `rho_0` | `rho0` | asymptotic bulk number density; 223 | 671 | 46 |
| `Phi` | `Phi` | outward sphere flux; repair directive item 10 | 1041 | 46 |
| `s9b_gamma` | `gamma` | reference-family gamma; 208-217 | 849 | 68 |
| `rho_bulk` | `rho` | local bulk number density; 219,248 | 672 | 92 |
| `s9b_fractional_bulk_change` | `fLocal` | local fractional bulk change; 222 | 1019 | 94 |
| `s9b_delta` | `deltaProfile` | live optical speed shift; 124,141 | 397 | 50 |
| `s9b_V` | `velocityProfile` | live radial material velocity; 60,141,175 | 398 | 50 |
| `s9b_xi` | `xiProfile` | live graph displacement; 83-87,141 | 399 | 51 |
| `rho_br_profile` | `rhoBr` | live brane mass density; 57,69,80 | 400 | 63 |
| `s9b_f` | `fProfile` | live fractional bulk density change; 222 | 401 | 53 |
| `mu_perp` | `muPerp` | live transverse stiffness; 52-57 | 402 | 63 |
| `j_n` | `jn` | live normal mass exchange; 69,80 | 673 | 67 |
| `s9b_mass_V` | `vForward` | forward no-loss velocity; repair item 10 | 1042 | 651 |
| `s9b_velocity_squared` | `yProfile` | velocity-square reduction; 236-240 | 447 | 504 |
