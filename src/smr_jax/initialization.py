"""Reference-compatible initial state construction for the two-layer target."""

import jax.numpy as jnp

from .legacy import TwoLayerLegacyInputs
from .state import TwoLayerState


def initialize_two_layer_reference_state(inputs: TwoLayerLegacyInputs) -> TwoLayerState:
    """Mirror the legacy initial total and A-layer expressions without fluxes.

    The legacy engine retains total storage separately from A/B storage during
    later sequential updates. Both representations therefore remain explicit.
    B is inferred only for the initial two-layer test state.
    """

    sat_amt = jnp.asarray(inputs.sat_amt_cm, dtype=jnp.float64)
    total = sat_amt * 0.4
    layer_a = jnp.asarray(
        (inputs.wiltpt_mc_a + (inputs.fieldcap_mc_a - inputs.wiltpt_mc_a) * 0.2)
        * inputs.soil_depth_a_cm,
        dtype=jnp.float64,
    )
    layer_b = total - layer_a
    mask = jnp.isfinite(sat_amt)
    return TwoLayerState(
        total_water_cm=total,
        layer_water_cm=jnp.stack((layer_a, layer_b)),
        active_mask=mask,
    )
