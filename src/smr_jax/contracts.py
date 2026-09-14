"""Explicit static-field and daily-forcing contracts for the frozen two-layer port."""

from dataclasses import dataclass
import jax

from .state import GridSpec


@dataclass(frozen=True)
class StaticInputs:
    """Dense source rasters; NaN represents a legacy GRASS NULL in float fields."""

    grid: GridSpec
    active_mask: jax.Array
    elevation_m: jax.Array
    slope_percent: jax.Array
    landuse: jax.Array
    watershed: jax.Array
    flowunits: jax.Array
    routing_weights: jax.Array  # [N, NE, E, SE, S, SW, W, NW, y, x]
    sat_amt_cm: jax.Array
    fieldcap_amt_cm: jax.Array
    fieldcap_a_cm: jax.Array
    fieldcap_b_cm: jax.Array
    wiltpt_amt_cm: jax.Array
    soil_depth_a_cm: jax.Array
    soil_depth_b_cm: jax.Array
    sat_mc_a: jax.Array
    sat_mc_b: jax.Array
    fieldcap_mc_a: jax.Array
    fieldcap_mc_b: jax.Array
    wiltpt_mc_a: jax.Array
    ksat_matrix_a: jax.Array
    ksat_matrix_b: jax.Array
    ksat_mpore_a: jax.Array
    ksat_mpore_b: jax.Array
    kfc_a: jax.Array
    kfc_b: jax.Array
    ksubsurface: jax.Array
    max_canopy_cm: jax.Array
    canopy_cover: jax.Array
    root_zone_cm: jax.Array
    et_reduction_mc: jax.Array
    tmax_rain_c: jax.Array
    tmin_snow_c: jax.Array
    basin_ids: jax.Array
    basin_area_cells: jax.Array
    basin_reservoir_initial: jax.Array
    basin_reservoir_coeff: jax.Array


@dataclass(frozen=True)
class DailyForcing:
    """One legacy weather row. Temperatures are four 6-hour station values."""

    year: jax.Array
    month: jax.Array
    day: jax.Array
    doy: jax.Array
    tavg_c: jax.Array
    tdew_c: jax.Array
    precip_cm: jax.Array
    pet_cm: jax.Array
    temperatures_c: jax.Array
    cloud: jax.Array
    roughness_urban: jax.Array
    roughness_crop: jax.Array
    et_coeff_by_landuse: jax.Array  # landuse 1..6
    replay_shortwave: jax.Array


jax.tree_util.register_dataclass(StaticInputs, data_fields=[f for f in StaticInputs.__dataclass_fields__ if f != "grid"], meta_fields=["grid"])
jax.tree_util.register_dataclass(DailyForcing, data_fields=list(DailyForcing.__dataclass_fields__), meta_fields=[])
