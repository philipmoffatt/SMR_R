# SMR modular JAX translation

## Boundary

This packet starts at CPU-06 commit `90c546f`. It adds only `src/smr_jax`, JAX
tests, package metadata, and this packet record. The Perl/R/GRASS production
path and `reference/mica-creek-jan2015` remain immutable.

## Compatibility decisions

- The engine preserves float64 and the operation order in `R/smr_buildup.pl`.
- The frozen bundle did not retain the slope raster produced by GRASS
  `r.slope.aspect`. The operator authorized a Horn 3x3 JAX slope-in-degrees
  emulation with explicit NULL-neighbour behavior. It is provisional pending
  the lateral-routing parity gate; it must not be presented as a demonstrated
  general GRASS replacement.
- The first runnable integration uses completed `q.srad` fields from the
  two-day trace. This replay provider is not general radiation.
- The receiver-weighted eight-neighbour gather is tested by named orientation,
  rather than substituted with donor-weighted D-infinity routing.
- The runnable engine's active mask is the finite intersection of all required
  legacy operands, yielding the frozen reference's 68,121-cell support. It
  does not replace NULLs with zero.

## Verification sequence

The build makes the ordered engine executable. Acceptance remains ordered as:
surface/snow, vertical soil, routing/reservoir, two then 30-day CPU parity,
general radiation, and GPU parity/benchmark. No parity claim is made by the
build-only tests.
