"""Precipitation lapse, partitioning, canopy interception, and road runoff."""

import jax.numpy as jnp


def precipitation_and_temperature(static, forcing):
    """Legacy daily lapse/partition expressions, returning arrays in cm and C."""
    elevation = static.elevation_m
    precip = (elevation * .09195 + 15.54) / 130.0 * forcing.precip_cm
    tavg = forcing.tavg_c + (-3.20 / 1000.0) * (elevation - 784.0)
    # The legacy dewpoint expression algebraically cancels elevation. Preserve it.
    tdew = forcing.tdew_c + (-3.20 / 1000.0) * (elevation - 784.0 - elevation)
    rain = jnp.where(tavg > static.tmax_rain_c, precip,
                     jnp.where(tavg < static.tmin_snow_c, 0.0,
                               (tavg - static.tmin_snow_c) /
                               (static.tmax_rain_c - static.tmin_snow_c) * precip))
    pet_data = 1.8281 * forcing.pet_cm - .2486
    pet = (pet_data - forcing.pet_cm) * (elevation - 784.0) / (784.0 - 900.0) + pet_data
    return {"precip_cm": precip, "rain_cm": rain, "snow_cm": precip - rain,
            "tavg_c": tavg, "tdew_c": tdew, "pet_cm": pet}


def canopy(state, static, rain_cm, pet_cm):
    """Sequential canopy updates; evaporation consumes PET before soil ET."""
    storage = state.canopy_storage_cm + rain_cm
    throughfall = jnp.where(storage > static.max_canopy_cm, storage - static.max_canopy_cm, 0.0)
    storage = storage - throughfall
    canopy_evap = jnp.minimum(storage, pet_cm)  # legacy branch has equal arms
    storage = storage - canopy_evap
    return storage, jnp.maximum(0.0, pet_cm - canopy_evap), {"throughfall_cm": throughfall, "canopy_evap_cm": canopy_evap}


def road_runoff(static, snowmelt_cm, throughfall_cm):
    """Preserve the road-only loss and proportional source adjustment."""
    road = jnp.where(static.landuse == 2.0, jnp.maximum(0.0, snowmelt_cm + throughfall_cm - 1.5) / 30.0, 0.0)
    total = snowmelt_cm + throughfall_cm
    safe_total = jnp.where(total == 0.0, 1.0, total)
    adjusted_melt = jnp.where(total > 0.0, snowmelt_cm * (total - road) / safe_total, snowmelt_cm)
    adjusted_throughfall = jnp.where(total > 0.0, throughfall_cm * (total - road) / safe_total, throughfall_cm)
    return adjusted_melt, adjusted_throughfall, road
