"""Canopy-to-ground snow energy, refreeze, and melt in legacy update order."""

from dataclasses import replace
import jax.numpy as jnp

from ..state import SnowState
from ..radiation import shortwave_from_components


SNOW_RESIDUE_CM = 1e-3  # Operator-selected convention, not a fitted tolerance.


def snow_energy(state, static, forcing, tavg_c, tdew_c, snow_cm, throughfall_cm, *,
                melt_small_snow_next_day=False, exponential=jnp.exp):
    """Daily SAM update with an optional, explicit next-day residue convention.

    The default preserves legacy behavior. In the selected mode, positive
    start-of-day SWE <= 1e-3 cm is released as melt before snow-presence
    decisions. Yesterday's water remains stored overnight; today's newly
    formed residue is left in the returned state until tomorrow. No canopy or
    litter store is involved. Thermal state is retained because it also
    represents the underlying surface, including when snow is absent.

    This forced small-water release is an operator-authorized approximation,
    not an energy-conserving phase-change calculation. The mode is static
    when used under JIT. Shortwave still comes from the replay provider.
    """
    snow = state.snow
    residual_melt = jnp.zeros_like(snow.swe_cm)
    if melt_small_snow_next_day:
        release = static.active_mask & (snow.swe_cm > 0.0) & (snow.swe_cm <= SNOW_RESIDUE_CM)
        residual_melt = jnp.where(release, snow.swe_cm, 0.0)
        clear = lambda value: jnp.where(release, 0.0, value)
        snow = replace(snow, swe_cm=clear(snow.swe_cm), previous_swe_cm=clear(snow.previous_swe_cm),
                       ice_cm=clear(snow.ice_cm), liquid_cm=clear(snow.liquid_cm))
    rh = jnp.where(snow.swe_cm > 0.0, forcing.roughness_urban / (1.0 - static.canopy_cover / 3.0),
                   forcing.roughness_crop / (1.0 - static.canopy_cover / 3.0))
    age = jnp.where((snow_cm > 0.0) & (throughfall_cm == 0.0), 1.0, snow.age_days + 1.0)
    albedo = jnp.where(snow.previous_swe_cm + snow_cm > 0.0,
                       jnp.minimum(.95, .7383 * age ** (-.1908)), .2)
    components = getattr(forcing, "radiation_components", None)
    if components is None:
        if forcing.replay_shortwave is None:
            raise ValueError("provide replay shortwave or prepared radiation components")
        q_srad = forcing.replay_shortwave
    else:
        if forcing.replay_shortwave is not None:
            raise ValueError("choose one radiation provider; replay and components cannot be mixed")
        q_srad = shortwave_from_components(components, albedo, static.canopy_cover)
    cloud = forcing.cloud
    q_lw = 5.67e-8 * 3600.0 * 24.0 / 1000.0 * (
        jnp.maximum(.72, (((9.2 * (10.0 ** -6.0)) * (tavg_c + 273.15) ** 2) * (1.0 - .84 * cloud) + .84 * cloud)
        * (1 - static.canopy_cover) + .92 * static.canopy_cover) * (tavg_c + 273.15) ** 4
        - .98 * (snow.surface_temp_c + 273.15) ** 4)
    vap_sat = exponential((16.78 * snow.surface_temp_c - 116.9) / (snow.surface_temp_c + 237.3)) / (273.15 + snow.surface_temp_c) / .4615
    vap_snow = jnp.where(snow.surface_temp_c < 0.0, vap_sat * exponential(335.0 * snow.surface_temp_c / (.4615 * (273.15 + snow.surface_temp_c) ** 2)), vap_sat)
    vap_air = exponential((16.78 * tdew_c - 116.9) / (tdew_c + 237.3)) / (tdew_c + 273.15) / .4615
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
    denominator = (2.1 * 1000.0 * jnp.minimum(swe / 100.0, .02) + 1000.0 * 2.1 * jnp.maximum(0.0, .02 - swe / 100.0))
    # The source divides sequentially; regrouping changes wet-period branches.
    temp = surface_energy / 24.0 / denominator
    next_snow = replace(snow, swe_cm=swe, previous_swe_cm=swe, ice_cm=ice, liquid_cm=liquid,
                        age_days=age, total_energy=total_energy, surface_energy=surface_energy, surface_temp_c=temp)
    return next_snow, {"snowmelt_cm": snowmelt + residual_melt, "swe_cm": swe, "condens_cm": condens,
                       "residual_snowmelt_cm": residual_melt,
                       "q_srad": q_srad, "q_lw": q_lw, "q_latent": q_latent,
                       "q_sensible": q_sensible, "q_total": q_total, "q_rain_ground": q_rain, "albedo": albedo}
