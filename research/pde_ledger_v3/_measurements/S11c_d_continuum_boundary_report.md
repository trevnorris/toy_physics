# Continuum boundary coefficient construction

Interior coefficient matrices are accepted at7d1082df. The adapter consumes
accepted reference/end modes and polarized current operands, retaining every
required complete subspace. It computes eta/sigma mode, boundary and current
coefficients in reference-anchored coordinates. Open bases are flux-normalized
at zero only; current variation remains explicit.

The source preflight joined all three saved end pencils and current operands,
18 candidates per end,24 current sources and14 accepted input packets. The
independent two-sided matrix-polynomial test agrees exactly in all four grades.
Its multiplication-order mutation responds in the mixed coefficient; the first
attempt tested first order, where the scalar reference momentum commutes, and
was correctly rejected as insensitive. Original test logs remain preserved.
Both recursive inverse products agree to floating precision. These are
instrument tests, not a new physical scattering result.

The production attempt has900seconds, one numerical thread and2GiB. It saves
each end pencil and mode cluster before later guards. No numerical integration
or upstream physics is repeated. Full boundary/current results remain pending.
The continuum solve, physical controls, flux bookkeeping, poles and final
cases/exports follow under exploratory acceptance.
