"""Legacy two-layer allocation, lateral outflow, ET, percolation, and runoff."""

import jax.numpy as jnp

from .routing import legacy_receive_gather
from ..state import TwoLayerState


def _reallocate(total, layer_a, static):
    """The repeated Perl A/B redistribution expressions, preserving total separately."""
    tmp_a = jnp.where(total < static.fieldcap_amt_cm, jnp.minimum(layer_a, static.fieldcap_a_cm), layer_a)
    layer_b = jnp.where(total < static.fieldcap_amt_cm, total - tmp_a,
                        jnp.where(total - static.fieldcap_a_cm < static.sat_mc_b * static.soil_depth_b_cm,
                                  total - static.fieldcap_a_cm, static.sat_mc_b * static.soil_depth_b_cm))
    layer_a = jnp.where(total < static.fieldcap_amt_cm, tmp_a,
                        jnp.where(total - layer_b < static.sat_mc_a * static.soil_depth_a_cm,
                                  total - layer_b, static.sat_mc_a * static.soil_depth_a_cm))
    return layer_a, layer_b


def allocate_input(soil, static, input_cm):
    total = soil.total_water_cm + input_cm
    tmp_a = jnp.where(total < static.fieldcap_amt_cm, jnp.minimum(input_cm + soil.layer_water_cm[0], static.fieldcap_a_cm), soil.layer_water_cm[0])
    a, b = _reallocate(total, tmp_a, static)
    return TwoLayerState(total, jnp.stack((a, b)), soil.active_mask)


def _effective_conductivity(storage, fieldcap, sat_mc, depth, matrix, mpore, kfc):
    safe_depth = jnp.where(depth == 0.0, 1.0, depth)
    value = jnp.where(storage < fieldcap, matrix * jnp.exp((-13.0 / sat_mc) * (sat_mc - storage / safe_depth)),
                      jnp.where(storage >= sat_mc * depth, mpore,
                                (mpore - kfc) * (storage / safe_depth - fieldcap / safe_depth) /
                                (sat_mc - fieldcap / safe_depth) + kfc))
    return value * .25


def lateral_route(soil, static):
    a, b = soil.layer_water_cm
    eff_a = _effective_conductivity(a, static.fieldcap_a_cm, static.sat_mc_a, static.soil_depth_a_cm, static.ksat_matrix_a, static.ksat_mpore_a, static.kfc_a)
    eff_b = _effective_conductivity(b, static.fieldcap_b_cm, static.sat_mc_b, static.soil_depth_b_cm, static.ksat_matrix_b, static.ksat_mpore_b, static.kfc_b)
    flow_a = jnp.minimum(eff_a * static.slope_percent / 100.0 * static.soil_depth_a_cm / 30.0 / 100.0, a)
    flow_b = jnp.minimum(eff_b * static.slope_percent / 100.0 * static.soil_depth_b_cm / 30.0 / 100.0, b)
    flow_a, flow_b = jnp.where(static.flowunits == 0.0, 0.0, flow_a), jnp.where(static.flowunits == 0.0, 0.0, flow_b)
    flow = flow_a + flow_b
    a = a - flow_a + legacy_receive_gather(flow_a, static.routing_weights, soil.active_mask)
    b = b - flow_b + legacy_receive_gather(flow_b, static.routing_weights, soil.active_mask)
    total = soil.total_water_cm - flow + legacy_receive_gather(flow, static.routing_weights, soil.active_mask)
    return TwoLayerState(total, jnp.stack((a, b)), soil.active_mask), flow


def evapotranspiration(soil, static, forcing, pet_cm):
    a = soil.layer_water_cm[0]
    root_storage = jnp.where(static.landuse == 3.0, a, soil.total_water_cm)
    root = static.root_zone_cm
    safe_root = jnp.where(root == 0.0, 1.0, root)
    landuse_index = static.landuse.astype(jnp.int32)
    known_landuse = (static.landuse == landuse_index) & (landuse_index >= 1) & (landuse_index <= 6)
    # The nested Perl expression falls back to cc_water outside classes 1..6.
    coefficient_index = jnp.where(known_landuse, jnp.clip(landuse_index, 1, 6) - 1, 0)
    et_coeff = forcing.et_coeff_by_landuse[coefficient_index]
    threshold = root * static.et_reduction_mc
    high = jnp.minimum(jnp.maximum(0.0, root_storage - static.wiltpt_amt_cm), pet_cm * et_coeff * .25)
    reduction_denominator = static.et_reduction_mc - static.wiltpt_amt_cm / safe_root
    safe_reduction_denominator = jnp.where(reduction_denominator == 0.0, 1.0, reduction_denominator)
    reduction = (root_storage / safe_root - static.wiltpt_amt_cm / safe_root) / safe_reduction_denominator
    low = jnp.minimum(jnp.maximum(0.0, root_storage - static.wiltpt_amt_cm), jnp.maximum(0.0, pet_cm * reduction * et_coeff * .25))
    # GRASS selects the high-storage arm even when root_zone is zero. Protect
    # the unselected ratios without suppressing that legacy ET flux.
    et = jnp.where(root_storage > threshold, high, jnp.where(root_storage > static.wiltpt_amt_cm, low, 0.0))
    total = soil.total_water_cm - et
    replacement_a = jnp.where(total < static.fieldcap_amt_cm, static.fieldcap_a_cm * total / jnp.where(static.fieldcap_amt_cm == 0.0, 1.0, static.fieldcap_amt_cm), static.fieldcap_a_cm)
    a = jnp.where(static.landuse == 3.0, a - et, jnp.where(known_landuse, replacement_a, 0.0))
    return TwoLayerState(total, jnp.stack((a, total - a)), soil.active_mask), et


def percolation_and_runoff(soil, static):
    a, b = soil.layer_water_cm
    perc = jnp.minimum(jnp.minimum(jnp.maximum(soil.total_water_cm - static.fieldcap_amt_cm, 0.0), b), static.ksubsurface * .25)
    total = soil.total_water_cm - perc
    b = b - perc
    a, b = _reallocate(total, a, static)
    runoff = jnp.where(total > static.sat_amt_cm, total - static.sat_amt_cm, 0.0)
    total = total - runoff
    return TwoLayerState(total, jnp.stack((a, b)), soil.active_mask), perc, runoff, total / static.sat_amt_cm * 100.0
