---
title: SMR external vegetation water interface
status: active
owner: SMR_R
last_reviewed: 2026-11-12
---

# Compatibility findings and decisions

This is an opt-in **bounded two-layer coupling path**, implemented in
`src/smr_jax/plant.py`. The existing `engine.daily_step`, `run_days`, flux
functions, frozen fixtures and tests are unchanged. This update does not choose
project ownership, replace either model, port the stand-in soil, or establish
field validity. The external path is selected by calling `external_water_step`;
there is no global switch or silent change to legacy execution.

## Inspected implementation and provenance

The owning Git repository is `C:/projects/git/SMR_R`, origin
`git@github.com:philipmoffatt/SMR_R.git`, default branch `master`, primary
checkout on `vib-gpu` at `4f6bbea`. Both JAX reference worktrees point to
`4b5733ec1a9da2a86572975f3b1e1564cdbd8616`:

- `PKT-SMR-JAX-01`, branch
  `packet/PKT-SMR-JAX-01/translate-the-modular-jax-engine-from-th`, clean.
- `PKT-SMR-JAX-VERIFY-01`, branch
  `packet/PKT-SMR-JAX-VERIFY-01/verify-smr-jax-initialization-and-surfac`, with
  modified `src/smr_jax/io.py` and untracked `src/smr_jax/verify_surface.py` and
  `tests/test_jax_trace_contract.py`. Those are surface-verification work, not
  completed plant coupling. They were read only and not copied into this branch.

The plant directory has no `.git` entry and `git rev-parse` cannot establish
an owning repository or revision. Its source contents, rather than a fictional
commit, are the reference. SHA256 snapshots of the two interface files:

| File under `C:/projects/git/plant-game-of-life` | SHA256 |
|---|---|
| `contracts/smr_plant_contract.py` | `478d9ecb62552a59e860e2ce342a6ea3a0292bcb5d226682df961b811d3d28a2` |
| `src/rangepft/smr_coupling.py` | `bbbd4b907f5fbf043005cda9129a6c814f2ad09b3cf4c17e359c1f03fa69a212` |

Also inspected: plant `soilbucket.py`, `water.py`, `state.py`, `run.py`,
`docs/smr_read.md`, `docs/algorithm_manifest.md`, `docs/literature_audit.md`;
SMR `contracts.py`, `state.py`, `engine.py`, `initialization.py`, `legacy.py`,
all relevant flux modules, tests, and conductivity setup in `R/SMR_data_setup.R`.

| Topic | Actual implementation before this update | Compatibility gap / decision |
|---|---|---|
| Soil state | `TwoLayerState`: `[2,y,x]` water in cm, separately carried total | Validate total equals layers; never silently repair a mismatch. |
| Soil properties | A/B thickness and saturation; A WP plus profile WP | B WP storage = source profile WP minus A WP; validate it. No invented profile extension. |
| Plant contract | `[L,C]` and daily fluxes; marked proposed | Duck-typed mapping, no dependency or claim that the proposal is approved. L must be 2 for this backend. |
| ET | Internal land-use coefficient ET and persistent canopy storage | External path owns interception/ET exclusively and runs neither internal function. |
| Uptake | Plant step advances growth before soil; adapter subtracts all requests | Return actual and unmet amounts; checked entry refuses unsupplied requests by default. Plant needs a transactional update. |
| Vertical flow | Implicit A/B redistribution, bottom percolation only | Record actual amounts moved at every redistribution, including upward movement. These are algebraic transfers, not a new physical drainage equation. |
| Routing | Receiver-weighted 8-neighbour gather; only summed outflow returned | Measure each layer's outgoing and arriving water directly. Validate each donor's delivery fraction does not exceed 1. |
| Plant flow estimates | A drainage estimated from B net gain; lateral split by initial storage | Replace estimates with actual fluxes from this interface. |
| Plant budget | Lateral input computed as remainder; absent balance keys skipped | Independent edge routing and strict required diagnostics; no residual-derived input. |
| Snow | Legacy energy model requires replay shortwave; prototype uses degree-day snow | Supply an explicit water-conserving snow handoff. No new snow equation or guessed radiation conversion. |
| Temperature | No SMR soil-temperature model | `soil_to_plant` requires an explicit caller-supplied temperature proxy. |
| Arbitrary layers | Plant stand-in supports L layers; SMR state and equations hardwire A/B | Preserve two-layer engine. New layer count requires a separate selectable backend and scientific decisions below. |

The November 11 proposals in plant `smr_read.md` to port the stand-in and give
SMR ownership of soil remain proposals for this task. The manifest identifies
the Campbell exponent, infiltration, lateral multiplier and restricted bottom
boundary as working/provisional rules. The literature audit is evidence of
uncertainty, not acceptance or a reason to change these equations here.

## Units, shapes and soil bounds

- Storage and every daily water flux are **cm of water**, float64. A storage
  amount is not volumetric moisture: `theta = water_cm / depth_cm`.
- Core soil and demand fields are `[2,y,x]`; cell fields are `[y,x]`.
  `vertical_down_cm` and `vertical_up_cm` are `[1,y,x]`, the sole A/B interface.
- Layer thickness is positive and saturation is positive;
  `0 <= WP <= FC <= saturation <= depth` is checked in cm. Initial water may
  be below WP but must be between zero and saturation. Extraction uses only
  `max(current_water - WP, 0)` and cannot deepen an existing WP deficit.
- Profile FC and saturation must equal the sum of their layer amounts. B WP
  cannot be negative or exceed B FC. Invalid/zero-depth restricted layers are
  rejected, not extended or filled with invented properties.
- The inspected R setup converts SSURGO `ksat.r` from **micrometres/second to
  cm/day**, multiplying by 8.64 (lines 336–339). `Ksat_mpore`, `Kfc` and
  `Ksubsurface` are derived from that same scale (lines 375–390). The JAX loader
  reads them without another conversion. Supply these conductivity fields in
  **cm/day**, not the stand-in's `ksat_cm_hr`. This describes the inspected
  setup, not a guarantee about externally produced rasters.
- The external path preserves the legacy conductivity function and lateral
  expression, including the six-hour `.25` factor and both `/100` factors.
  `slope_percent` **contains degrees**. Do not convert it to slope percent or
  to a tangent. Grid length in the lateral expression is `cellsize_m`, equal
  to the legacy 30 on the reference grid. These inherited empirical factors
  need scientific review before any new physical interpretation.
- No AWC fragment correction, profile extrapolation, TWI water input, new
  pedotransfer properties, or productivity tuning is performed.

`demand_from_plant(plant, static)` maps the three water fields of a duck-typed
`PlantToSoil` to the grid. Contract cells are **packed active cells**, row-major
`np.flatnonzero(static.active_mask.ravel())`, not every grid slot. This equals
the prototype's flattening for its all-active synthetic grids. A masked-grid
plant adapter must use this packing for drivers and plant state too.
`ground_light_frac` is validated in `[0,1]` by the mapper but has no further
effect on SMR; its effect must already be in requested evaporation. There is
no second light multiplier, PET calculation or PET subtraction.

`soil_to_plant(soil, static, temperature_proxy_c, previous_diagnostics)` returns
kwargs for `SoilToPlant(**view)`. `depth_cm` means thickness, not cumulative
depth. `drain_cm` is previous-day downward A/B transfer and bottom drainage;
`lateral_cm` is previous-day actual outflow by layer. Supply explicit zero
prior flux arrays on day one. Missing previous-day flows raise. Upward
redistribution is available separately as `vertical_up_cm` and cannot be
represented in the original plant contract's downward-only field.

## Daily ordering and realized extraction

`SurfaceWater` explicitly provides daily `precip_cm`, `rain_cm`,
`snow_input_cm`, `melt_cm`, `swe_start_cm`, `swe_end_cm`, `snow_vapor_cm`.
The last is signed, positive for condensation and negative for sublimation.
All other surface amounts are nonnegative. Snow conservation and
`precip = rain + snowfall` are checked independently. With no snow, explicitly
set all five snow fields to zero; diagnostics cannot simply be omitted.
`melt_cm` excludes the liquid rain already in `rain_cm`.

1. One upstream provider partitions precipitation and releases snow water.
   The existing prototype's degree-day provider can supply this handoff; this
   API does not endorse its melt parameters. SMR's legacy `snowmelt_cm` needs
   a reviewed adapter because its expression and throughfall bookkeeping
   cannot safely be equated to this field by name alone.
2. External interception evaporates immediately, capped at rain plus snow
   release. Return requested, actual and unmet interception. There is no
   persistent canopy store. Never discard a nonzero legacy canopy store when
   changing paths; an explicit store-conserving transition is separate work.
3. Divide remaining liquid and both layer ET requests evenly among four
   six-hour substeps. These inputs are **daily totals**, not rates. A later
   wetting event cannot retrospectively supply an earlier substep's uptake;
   unused substep requests are reported unmet, not borrowed from another day.
4. In each substep, add liquid to A, move water toward the legacy A/B allocation
   target, and emit profile overflow as surface runoff/return flow. Keep both
   layers within zero/saturation. No conductivity-limited surface infiltration
   rule is introduced. `infiltration_cm` is liquid entering this soil operator
   before its saturation-excess runoff, not a measured field infiltration rate.
5. Compute bounded layer outflow and route it with literal receiver weights.
   Reconcile the A/B allocation target again; record upward/downward transfer
   and any return flow caused by arriving water.
6. Realize transpiration and soil evaporation from the available water in
   **each** layer. If combined demand exceeds that availability, scale the two
   competing requests proportionally. Never borrow another layer's water for
   an unsupplied layer request. Return actual and unmet daily amounts plus
   each substep's extraction and availability.
7. For `bottom_boundary="legacy_drainage"`, apply the existing excess-profile
   percolation limit `min(max(total-FC,0), B, Ksubsurface/4)`. For `"closed"`,
   bottom outflow is exactly zero. Then perform bounded legacy redistribution
   once more and record its transfer. No restricted-conductivity fraction,
   Campbell drainage or new hydraulic gradient is chosen.

This ordering intentionally follows the prototype's rain/routing/ET/bottom
sequence, with bounds and measured transfers. It does not claim exact legacy
parity when bounds engage or external ET differs. A regression test shows
agreement with the literal legacy soil sequence when those differences are
inactive. The original legacy engine remains available for all parity work.
The external entry takes a soil state, not a full `ModelState`: it applies no
Mica Creek lapse, road-loss, reservoir recurrence or internal surface energy
calculation. Those remain explicit integration responsibilities.

`external_water_step(...)` validates inputs and budgets, and raises on any
unmet request by default. Use `allow_unmet=True` only to obtain negotiation or
diagnostic results. Its return is `(candidate_soil, diagnostics)`; this option
does **not** authorize committing a previously calculated plant state.
`require_realized_uptake(diagnostics)` supplies the same fail-closed gate for a
caller composing a larger transaction. An unmet interception or evaporation
request also invalidates the plant's assumed water/energy split.

`external_day_kernel` is the pure JIT/scan-compatible core. Validate inputs
before compilation, declare `bottom_boundary` static, and verify every day's
result after executing a scan. A conservative budget alone does not approve
plant growth: also enforce the realized-uptake gate. The core deliberately
contains no host callbacks or hidden state.

## Diagnostics and independent water evidence

| Diagnostic | Meaning |
|---|---|
| `vertical_down_cm`, `vertical_up_cm` | Sum of actual movements at the A/B interface, separately in each direction. Do not cancel them for solute use. |
| `lateral_out_cm`, `lateral_in_cm` | Actual per-layer outgoing water and receiver gather. Inflow is never a budget remainder. |
| `lateral_export_cm` | Per-donor outflow times its unreceived routing fraction. |
| `bottom_drainage_cm` | Actual bottom export. |
| `runoff_cm` | Actual surface excess/return flow, including routing overflow. |
| `interception_evap_cm`, `transp_cm`, `soil_evap_cm` | Actual evaporation/extraction; transpiration and soil evaporation are layer fields. |
| `requested_*_cm`, `unmet_*_cm` | Daily requested and unsupplied amounts for all three vegetation processes. |
| `substep_water_cm`, `substep_available_cm`, `substep_transp_cm`, `substep_soil_evap_cm` | Four substep states, available extraction and realized components. |
| `layer_residual_cm`, `cell_residual_cm`, `domain_residual_m3` | Independent layer, cell and domain ledgers; domain uses actual routing exports. |
| `snow_residual_cm`, `surface_partition_residual_cm` | Independently checked snow storage and precipitation partition. |

Receiver weights can normalize at receivers while duplicating a donor's
water. For each donor, this interface explicitly sums all active receiver
claims using the transpose of the neighbour stencil. A claim above one is
rejected (roundoff allowance `1e-12`); weights are never silently normalized.
The unclaimed share is an **open routing export**. It includes outlet/edge
loss and unassigned interior fractions; the legacy receiver map cannot
distinguish those physically. It is not an inferred natural watershed outlet.
Inactive cells neither send nor receive. Different geometric grid areas are
not supported: the grid has one uniform `cell_area_m2`.

For layer i, with signed A-to-B transfer counted as an outgoing flux for A
and incoming for B, the tested ledger is:

```text
residual_i = final_water_i - initial_water_i - surface_input_i
             + vertical_out_i - vertical_in_i
             + lateral_out_i - lateral_in_i
             + bottom_export_i + surface_runoff_i
             + realized_transp_i + realized_soil_evap_i
```

The cell ledger independently includes change in total soil plus SWE,
precipitation, signed snow vapor, interception evaporation, extraction,
runoff, bottom drainage and measured net lateral outflow. The domain ledger
compares the domain's soil/SWE storage against precipitation/vapor and all
external sinks, using **routing exports**, not the sum of cell residuals.
All cm amounts convert to m3 with `0.01 * cell_area_m2`.
Bottom drainage, runoff and routing exports leave the modeled soil/snow
domain. If later supplied to a reservoir, they are reservoir inputs; do not
count them a second time as watershed exports in a larger ledger.

`verify_water_balance` recomputes those ledgers with NumPy, rebuilds routing
with explicit donor/receiver slices, checks substep bounds and requested vs
realized amounts, and rejects missing/nonfinite diagnostics. Default numerical
tolerance is `1e-10 cm`; domain tolerance is that depth times cell area times
active cell count. No observation or fitted threshold supplies this tolerance.
The synthetic tests also calculate their own layer balances independently of
the reported residuals and inject faulty/missing diagnostics to prove failure.

## Arbitrary layers: supported boundary and remaining science

The smallest current compatibility implementation is the two-layer interface
above. It exposes real thicknesses and explicitly rejects other layer counts;
it does not reshape ten plant layers into two or extend a shallow restricted
profile. Tests vary the A/B boundary (10/90, 30/70, 50/50 cm) while conserving
profile depth, WP/FC/saturation and initial stocks. Proportional extraction in
a closed column preserves profile totals; mixing and root access can still
depend on where that boundary lies. This is **not** evidence of convergence
with arbitrary L, nor a test of the plant's physiological layer sensitivity.

A future selectable backend would retain this water/demand and ledger boundary
but need its own `[L,y,x]` state/properties, explicit adjacent-interface fluxes,
and tests over different **layer counts** at conserved profile properties.
It must preserve `TwoLayerState` and the existing `engine` as the legacy
choice. Before implementing that backend, resolve vertical flux equations,
below-FC leakage, restricted bottom conditions, storm/infiltration timing,
lateral delivery across unequal profiles, restriction geometry, and upward
flow. The stand-in's Campbell/Darcy approximations, `bottom_k_frac=.01`, hourly
step and lateral multiplier are candidate choices, not defaults accepted here.
No dependent implementation of those unresolved choices was made.

## Required plant-side follow-ups (no plant files edited)

1. Replace the hand-composed SMR substeps in `rangepft/smr_coupling.py` with the
   checked interface and replace inferred drainage/lateral fields with actual
   diagnostics. Preserve the previous-day timing of plant leaching unless
   explicitly changing that scientific choice. Upward redistribution and
   lateral solute arrival require a separate solute design; these water
   diagnostics alone do not validate nitrate transport.
2. Split proposal from commitment. The current `plant_step` already computes
   growth using its proposed uptake. If SMR supplies less, reject that candidate
   plant state and recompute/finalize water-dependent plant processes using
   realized per-layer uptake. Aggregated water demands do not retain per-type
   allocation, so SMR cannot safely repair plant growth itself. Rerunning the
   same request without a supply-aware plant step is insufficient. Similarly,
   reconcile interception/soil-evap assumptions if those are capped.
3. Use actual thickness, water, WP and substep availability to diagnose partial
   rooting into a mixed layer. This interface does not impose a rooting model
   or allow the root-accessible fraction to be inferred from whole-layer water.
4. Keep `water_avail_rootzone` labeled as a profile-scaled index, not extractable
   cm. Do not equate it with `max(water-WP,0)` or with realized uptake.
5. Diagnose the existing small growth floor at zero water and reserve-to-shoot
   transfer separately from new carbon assimilation. SMR changes neither.
6. Keep one shared reference ET driver; do not pass plant requests as a second
   PET estimate. Supply explicit surface/snow bookkeeping, temperature proxy,
   and masked-cell packing. Review any use of legacy snow/radiation before
   bridging it into this coupling.

RUE, drought responses, rooting parameters, physiology, RAP production and
standing-biomass definitions are unchanged. No held-out observations were
used for tuning, equation choice or acceptance. Existing splits remain intact.

## Run the synthetic example and checks

With project dependencies available:

```powershell
python examples/plant_water.py
python -m unittest discover -s tests -p 'test_jax_*.py' -v
```

The example constructs synthetic inputs explicitly and shows actual uptake,
unmet requests, water residuals and next-day contract kwargs. It imports no
plant repository and makes no field-validation claim. Execution environment,
results and EOS scope evidence are recorded in
`docs/packets/PKT-SMR-JAX-PLANT-01.md`.
