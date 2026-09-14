"""Canopy-to-ground snow energy, refreeze, and melt in legacy update order."""

from dataclasses import replace
import jax.numpy as jnp

from ..state import SnowState


def snow_energy(state, static, forcing, tavg_c, tdew_c, snow_cm, throughfall_cm):
    """Run one daily SAM update; shortwave comes only from the replay provider."""
    snow = state.snow
    rh = jnp.where(snow.swe_cm > 0.0, forcing.roughness_urban / (1.0 - static.canopy_cover / 3.0),
                   forcing.roughness_crop / (1.0 - static.canopy_cover / 3.0))
    age = jnp.where((snow_cm > 0.0) & (throughfall_cm == 0.0), 1.0, snow.age_days + 1.0)
    albedo = jnp.where(snow.previous_swe_cm + snow_cm > 0.0,
                       jnp.minimum(.95, .7383 * age ** (-.1908)), .2)
    q_srad = forcing.replay_shortwave
    cloud = forcing.cloud
    q_lw = 5.67e-8 * 3600.0 * 24.0 / 1000.0 * (
        jnp.maximum(.72, ((9.2e-6 * (tavg_c + 273.15) ** 2) * (1.0 - .84 * cloud) + .84 * cloud)
        * (1 - static.canopy_cover) + .92 * static.canopy_cover) * (tavg_c + 273.15) ** 4
        - .98 * (snow.surface_temp_c + 273.15) ** 4)
    vap_sat = jnp.exp((16.78 * snow.surface_temp_c - 116.9) / (snow.surface_temp_c + 237.3)) / (273.15 + snow.surface_temp_c) / .4615
    vap_snow = jnp.where(snow.surface_temp_c < 0.0, vap_sat * jnp.exp(335.0 * snow.surface_temp_c / (.4615 * (273.15 + snow.surface_temp_c) ** 2)), vap_sat)
    vap_air = jnp.exp((16.78 * tdew_c - 116.9) / (tdew_c + 237.3)) / (tdew_c + 273.15) / .4615
    q_latent = (2500.0 + 335.0) * (vap_air - vap_snow) / (rh / (3600.0 * 24.0))
    condens = jnp.where(snow.previous_swe_cm + snow_cm > 0.0, q_latent / ((2500.0 + 335.0) * 1000.0) * 100.0, 0.0)
    q_sensible = 1.0 * 1.29 * (tavg_c - snow.surface_temp_c) / (rh / (3600.0 * 24.0))
    q_rain = 4.2 * 1000.0 * (throughfall_cm / 100.0) * (jnp.maximum(0.0, tavg_c) - snow.surface_temp_c) + 2.1 * 1000.0 * (snow_cm / 100.0) * (jnp.minimum(0.0, tavg_c) - snow.surface_temp_c)
    q_total = q_srad + q_lw + q_sensible + q_latent + q_rain
    energy_term = 100.0 * (snow.total_energy + q_total) / (335.0 * 1000.0)
    refreeze = jnp.minimum(throughfall_cm + snow.liquid_cm, jnp.maximum(0.0, -energy_term)) - jnp.minimum(snow_cm + snow.ice_cm + condens, jnp.maximum(0.0, energy_term))
    total_energy = snow.total_energy + q_total + refreeze / 100.0 * 335.0 * 1000.0
    ice = jnp.where(total_energy > .000001, 0.0, snow.ice_cm + snow_cm + refreeze + condens)
    liquid = jnp.where((total_energy > .000001) | (total_energy < -.000001), 0.0,
                       jnp.minimum(.03 / .97 * ice, snow.liquid_cm + throughfall_cm - refreeze))
    swe = liquid + ice
    snowmelt = (snow.previous_swe_cm - swe) + snow_cm + condens
    factor = jnp.minimum(swe / 100.0, .02) / .05 + jnp.maximum(.02 - swe / 100.0, 0.0) / .3
    surface_energy = jnp.where(total_energy == 0.0, 0.0,
                               jnp.where(swe > 0.0, jnp.minimum(0.0, snow.surface_energy + (q_total + refreeze / 100.0 * 335.0 * 1000.0) * factor), snow.surface_energy + (q_total + refreeze / 100.0 * 335.0 * 1000.0) * factor))
    denominator = 24.0 * (2.1 * 1000.0 * jnp.minimum(swe / 100.0, .02) + 1000.0 * 2.1 * jnp.maximum(0.0, .02 - swe / 100.0))
    temp = surface_energy / denominator
    next_snow = replace(snow, swe_cm=swe, previous_swe_cm=swe, ice_cm=ice, liquid_cm=liquid,
                        age_days=age, total_energy=total_energy, surface_energy=surface_energy, surface_temp_c=temp)
    return next_snow, {"snowmelt_cm": snowmelt, "swe_cm": swe, "condens_cm": condens,
                       "q_srad": q_srad, "q_lw": q_lw, "q_latent": q_latent,
                       "q_sensible": q_sensible, "q_total": q_total, "albedo": albedo}
