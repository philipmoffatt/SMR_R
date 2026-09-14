"""Ordered daily SMR JAX engine: a daily scan enclosing four soil substeps."""

from dataclasses import replace
import jax
import jax.numpy as jnp

from .fluxes.basin import basin_recurrence
from .fluxes.snow import snow_energy
from .fluxes.soil import _reallocate, allocate_input, evapotranspiration, lateral_route, percolation_and_runoff
from .fluxes.surface import canopy, precipitation_and_temperature, road_runoff
from .state import ModelState


def daily_step(state: ModelState, forcing, static):
    """One exact legacy ordering pass. Radiation is injected by its provider."""
    surface = precipitation_and_temperature(static, forcing)
    canopy_storage, soil_pet, canopy_flux = canopy(state, static, surface["rain_cm"], surface["pet_cm"])
    snow, snow_flux = snow_energy(state, static, forcing, surface["tavg_c"], surface["tdew_c"], surface["snow_cm"], canopy_flux["throughfall_cm"])
    melt, throughfall, road = road_runoff(static, snow_flux["snowmelt_cm"], canopy_flux["throughfall_cm"])
    temp = forcing.temperatures_c[:, None, None] - (-3.20 / 1000.0) * (784.0 - static.elevation_m)[None, :, :]
    temp_sum = jnp.sum(jnp.maximum(0.0, temp), axis=0)
    zeros = jnp.zeros_like(state.soil.total_water_cm)

    def soil_substep(carry, temperature):
        soil, et_daily, perc_daily, runoff_daily, lateral_daily, saturation = carry
        denominator = jnp.where(temp_sum == 0.0, 1.0, temp_sum)
        water_input = melt * jnp.maximum(temperature, 0.0) / denominator + throughfall * .25
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
                   "storage_amt_a_cm": soil.layer_water_cm[0], "storage_amt_b_cm": soil.layer_water_cm[1]}
    return next_state, diagnostics


def run_days(initial_state, forcings, static):
    """Execute arbitrary daily forcing stacks through one outer ``lax.scan``."""
    return jax.lax.scan(lambda carry, day: daily_step(carry, day, static), initial_state, forcings)
