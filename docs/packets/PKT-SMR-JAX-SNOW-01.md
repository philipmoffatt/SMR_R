---
title: Small residual snow released as next-day melt
status: active
owner: SMR_R
last_reviewed: 2026-11-12
---

# PKT-SMR-JAX-SNOW-01

Implemented in the isolated packet worktree, pending feature integration.
Base: `400765e47d270354cbadd77463181f130f06b4f8`, the current
`feature/smr-jax-plant-water` tip. The original VERIFY worktree remains intact.
Its committed trace reader/reporting tools from `35312bd` were merged only
into this isolated packet branch to support the investigation; that tooling
merge does not close VERIFY-01.

## Operator decision and behavior

The operator selected **1e-3 cm SWE** (0.01 mm water equivalent) as the cutoff
and explicitly asked to add that residual to melt **the next day**. This is a
numerical convention chosen in conversation, not a threshold fitted to the
Mica Creek trace or held-out observations.

Select it explicitly:

```powershell
python -m smr_jax.cli --melt-small-snow-next-day
python -m smr_jax.verify_surface --melt-small-snow-next-day
```

Or use `daily_step(..., melt_small_snow_next_day=True)` / the same option on
`run_days`. Bind the Python boolean statically when compiling with JAX. The
default remains legacy behavior for parity work.

In the selected mode:

1. A day's newly formed small SWE remains in the existing snow state overnight.
   No water is discarded, and no extra pending/canopy/litter store is added.
2. At the next day's start, active cells with `0 < SWE <= 1e-3 cm` release
   their entire SWE once. Their SWE, previous SWE, ice and liquid components
   are cleared before aerodynamic resistance, albedo and condensation use
   the snow state. A new snowfall can establish a new pack that same day;
   it does not absorb or postpone yesterday's scheduled release.
3. The released amount is added to that day's `snowmelt_cm`. It is excluded
   from the intervening snow calculation, so it cannot be melted twice or
   refrozen into the new pack by that calculation. An initially supplied tiny
   pack is likewise treated as an overnight residue at the first call.
4. Apply the existing road-runoff fraction to that water, then distribute the
   retained residue evenly over the four soil steps. This is necessary even
   on a day with four nonpositive temperatures: the legacy temperature-based
   melt weights would otherwise sum to zero and lose the released water.
   Ordinary energy-derived melt keeps its existing temperature weights.
5. Zero, negative, inactive and above-cutoff SWE do not trigger forced release.
   Negative legacy snow states are not silently converted to positive water.

The thermal state is retained because the legacy snow fields also represent
underlying surface heat when SWE is zero. This is a forced tiny-water melt
convention, including on cold days; it does not enforce a new latent-heat cost
or constitute an energy-conserving phase-change solver. No new soil,
interception, litter, drought or growth equations are introduced.

## Diagnostics and independent checks

Three additive daily fields identify the transfer:

- `residual_snowmelt_cm`: water removed from the prior snow state and included
  in today's total `snowmelt_cm`.
- `residual_snowmelt_to_soil_cm`: that source after the existing road fraction.
- `residual_snowmelt_road_runoff_cm`: its road-runoff share, already included
  in the existing total runoff; do not add it again.

These fields are zero in legacy mode. The CLI summary records
`snow_residue_policy` and `residual_snowmelt_m3`; existing summary fields retain
their meaning. The residue volume converts cell-cm with `0.01 * cell_area_m2`.

Test-owned arithmetic checks the snow ledger
`end_SWE - start_SWE - snowfall - condensation + total_melt` and the residue
route `released - to_soil - road_share`. Synthetic cold-column tests also
check the actual soil-storage gain, while multi-day tests show that a release
occurs exactly once. These are checks of the snow/residue change, not a claim
that the inherited whole-domain hydrology is now independently verified.

## Validation and remaining verification

On Windows CPU, JAX `0.12.0.1`, float64, using the existing read-only plant
interpreter and the isolated temporary tifffile dependency:

- Baseline before implementation: all 31 inherited JAX tests passed.
- Added the delayed-release test first; it failed on the then-missing API.
- All **41 tests passed** after implementation in 30.174 seconds, including
  ten new residue tests: below/at/above threshold, next-day timing, component
  clearing, new snowfall, resistance branch, negative/inactive values, legacy
  selection, cold-day delivery, road partition and multi-day JIT release.
- Captured all 43 two-day legacy final-state/history arrays before modifying
  flux code. Every existing array exactly equals its post-change legacy-mode
  counterpart, including NaN locations. No frozen fixture was rewritten.
- The selected two-day CLI ran successfully on 68,121 active cells. It
  released no residue on January 1 and residue in 2,884 cells on January 2,
  totaling `1.4408474413585282e-12 m3`. Both independently calculated snow and
  residue-routing residuals were exactly zero on the active cells on both
  days. Synthetic tests also exercise appreciable residues up to 1e-3 cm.
- CLI selection/help and `git diff --check` passed.

**VERIFY-01 remains open.** The frozen trace uses strict `SWE > 0` and its own
near-zero residues. Relative to that unchanged trace, selected-mode January 2
SWE RMSE/max error are `0.01353717036 / 0.09641091006 cm`; sensible energy
RMSE/max error are `2367.74737873 / 16902.05526579` in inherited energy units.
Some discrepancies grow because this convention intentionally removes tiny
packs that the frozen run treated as snow. These are reported departures,
not acceptance thresholds or evidence of field validity.

The replay shortwave field also embeds the old run's albedo. It remains a
diagnostic replay provider, not proof of self-consistent radiation under a
changed snow convention. Closing the gate still requires explicit comparison
support and acceptance assertions and a reference that applies the same
convention. Soil, routing/reservoir, 30-day, general-radiation and GPU gates
remain separate. Neither a zero reporter exit code nor these unit tests is
declared full-model parity.

## EOS scope and handoff

Live toolkit create/preflight passed with no preflight warnings before edits.
The imported VERIFY tooling and the new implementation are covered by the
adjacent YAML, identical to the controller copy under
`worktrees/packet-controls/smr-snow-20261112/`. Final report-only scope has ten
in-scope paths and no violations. The default/feature branch tips and all
reference worktrees remain unchanged by this packet. No merges into a feature
or default branch, pushes, worktree removal or global board edits were made.
The packet is implementation-complete and awaits review/integration; its YAML
stays active until that reconciliation. No follow-up is inferred as complete.
