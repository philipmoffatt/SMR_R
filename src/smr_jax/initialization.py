"""Reference-compatible initial state construction for the two-layer target."""

import jax.numpy as jnp

from .contracts import StaticInputs
from .legacy import TwoLayerLegacyInputs
from .state import BasinState, ModelState, SnowState, TwoLayerState


def initialize_two_layer_reference_state(inputs: TwoLayerLegacyInputs) -> TwoLayerState:
    """Mirror the legacy initial total and A-layer expressions without fluxes.

    The legacy engine retains total storage separately from A/B storage during
    later sequential updates. Both representations therefore remain explicit.
    B is inferred only for the initial two-layer test state.
    """

    sat_amt = jnp.asarray(inputs.sat_amt_cm, dtype=jnp.float64)
    total = sat_amt * 0.4
    wilt = jnp.asarray(inputs.wiltpt_mc_a, dtype=jnp.float64)
    fieldcap = jnp.asarray(inputs.fieldcap_mc_a, dtype=jnp.float64)
    moisture_range = (fieldcap-wilt).astype(jnp.float32).astype(jnp.float64)
    layer_a = (wilt+moisture_range*.2)*jnp.asarray(inputs.soil_depth_a_cm, dtype=jnp.float64)
    layer_b = total - layer_a
    mask = jnp.isfinite(sat_amt)
    return TwoLayerState(
        total_water_cm=total,
        layer_water_cm=jnp.stack((layer_a, layer_b)),
        active_mask=mask,
    )


def initialize_legacy_state(static: StaticInputs) -> ModelState:
    """All mutable maps initialized in the order and values in ``smr_buildup.pl``."""
    total = static.sat_amt_cm * .4
    # The inner difference has two imported FCELL operands. GRASS rounds it
    # before the decimal 0.2 promotes the enclosing expression to DCELL.
    moisture_range = (static.fieldcap_mc_a - static.wiltpt_mc_a).astype(jnp.float32).astype(jnp.float64)
    a = (static.wiltpt_mc_a + moisture_range * .2) * static.soil_depth_a_cm
    soil = TwoLayerState(total, jnp.stack((a, total - a)), static.active_mask)
    zero = jnp.zeros_like(total)
    snow = SnowState(zero, zero, zero, zero, jnp.ones_like(total),
                     24.0 * -.5 * (2.1 * 1000.0 * zero / 100.0 + 2.1 * 1000.0 * .4),
                     24.0 * -2.0 * (2.1 * 1000.0 * jnp.minimum(zero / 100.0, .02) + 1000.0 * 2.1 * jnp.maximum(0.0, .02 - zero / 100.0)),
                     jnp.full_like(total, -2.0))
    return ModelState(soil, zero, snow, BasinState(static.basin_reservoir_initial, jnp.zeros_like(static.basin_reservoir_initial)), zero, zero)
