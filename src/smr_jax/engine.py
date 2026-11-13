"""Ordered daily SMR JAX engine: a daily scan enclosing four soil substeps."""

from dataclasses import replace
import jax
import jax.numpy as jnp

from .fluxes.basin import basin_recurrence
from .fluxes.snow import snow_energy
from .fluxes.soil import _reallocate, allocate_input, evapotranspiration, lateral_route, percolation_and_runoff
from .fluxes.surface import canopy, precipitation_and_temperature, road_runoff
from .state import ModelState


def daily_step(state: ModelState, forcing, static, *, melt_small_snow_next_day=False):
    """Legacy ordering, optionally releasing <= 1e-3 cm SWE on the next day.

    The opt-in residue release is spread evenly over four soil steps, including
    cold days. Ordinary energy-derived melt keeps the legacy temperature
    weights. The mode is a static Python option under JIT.
    """
    surface = precipitation_and_temperature(static, forcing)
    canopy_storage, soil_pet, canopy_flux = canopy(state, static, surface["rain_cm"], surface["pet_cm"])
    snow, snow_flux = snow_energy(state, static, forcing, surface["tavg_c"], surface["tdew_c"], surface["snow_cm"], canopy_flux["throughfall_cm"],
                                  melt_small_snow_next_day=melt_small_snow_next_day)
    melt, throughfall, road = road_runoff(static, snow_flux["snowmelt_cm"], canopy_flux["throughfall_cm"])
    temp = forcing.temperatures_c[:, None, None] - (-3.20 / 1000.0) * (784.0 - static.elevation_m)[None, :, :]
    temp_sum = jnp.sum(jnp.maximum(0.0, temp), axis=0)
    zeros = jnp.zeros_like(state.soil.total_water_cm)
    residue_to_soil, residue_to_road = zeros, zeros
    if melt_small_snow_next_day:
        # Apply precisely the same road-loss fraction as the combined sources.
        # Routing residue by the ordinary positive-temperature weights would
        # silently lose the released water when all four temperatures <= 0.
        surface_input = snow_flux["snowmelt_cm"] + canopy_flux["throughfall_cm"]
        safe_input = jnp.where(surface_input == 0.0, 1.0, surface_input)
        retained_fraction = jnp.where(surface_input > 0.0, (surface_input - road) / safe_input, 1.0)
        residue_to_soil = snow_flux["residual_snowmelt_cm"] * retained_fraction
        residue_to_road = snow_flux["residual_snowmelt_cm"] - residue_to_soil

    def soil_substep(carry, temperature):
        soil, et_daily, perc_daily, runoff_daily, lateral_daily, saturation = carry
        denominator = jnp.where(temp_sum == 0.0, 1.0, temp_sum)
        water_input = melt * jnp.maximum(temperature, 0.0) / denominator + throughfall * .25
        if melt_small_snow_next_day:
            water_input = ((melt - residue_to_soil) * jnp.maximum(temperature, 0.0) / denominator
                           + throughfall * .25 + residue_to_soil * .25)
        soil = allocate_input(soil, static, water_input)
        soil, lateral = lateral_route(soil, static)
        a, b = soil.layer_water_cm
        from .state import TwoLayerState
        a, b = _reallocate(soil.total_water_cm, a, static)
        soil = TwoLayerState(soil.total_water_cm, jnp.stack((a, b)), soil.active_mask)
        soil, et = evapotranspiration(soil, static, forcing, soil_pet)
        soil, perc, runoff, saturation = percolation_and_runoff(soil, static)
        return (soil, et_daily + et, perc_daily + perc, runoff_daily + runoff, lateral_daily + lateral, saturation), None

    init = (state.soil, zeros, zeros, zeros, zeros, zeros)
    (soil, et, perc, runoff, lateral, saturation), _ = jax.lax.scan(soil_substep, init, temp)
    runoff = runoff + road
    basin, discharge = basin_recurrence(state.basin, static, runoff, perc)
    next_state = replace(state, soil=soil, canopy_storage_cm=canopy_storage, snow=snow, basin=basin,
                         annual_runoff_cm=state.annual_runoff_cm + runoff,
                         annual_saturation_count=state.annual_saturation_count + (saturation >= 100.0))
    diagnostics = {**surface, **canopy_flux, **snow_flux, "road_runoff_cm": road,
                   "actual_et_cm": et, "percolation_cm": perc, "runoff_cm": runoff,
                   "lateral_flow_cm": lateral, "saturation_percent": saturation,
                   "discharge_cms": discharge, "storage_amt_cm": soil.total_water_cm,
                   "storage_amt_a_cm": soil.layer_water_cm[0], "storage_amt_b_cm": soil.layer_water_cm[1],
                   "residual_snowmelt_to_soil_cm": residue_to_soil,
                   "residual_snowmelt_road_runoff_cm": residue_to_road}
    return next_state, diagnostics


def run_days(initial_state, forcings, static, *, melt_small_snow_next_day=False):
    """Daily scan; bind the optional snow-residue mode statically under JIT."""
    return jax.lax.scan(lambda carry, day: daily_step(carry, day, static,
                        melt_small_snow_next_day=melt_small_snow_next_day), initial_state, forcings)
