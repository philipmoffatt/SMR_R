"""Literal GRASS receiving-cell weighted eight-neighbour routing gather."""

import jax.numpy as jnp


_OFFSETS = ((-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1))


def _neighbor(field, row_offset: int, col_offset: int):
    """Return ``field[row + offset]`` with GRASS NULL borders represented by zero."""
    h, w = field.shape
    padded = jnp.pad(field, 1, constant_values=jnp.nan)
    shifted = padded[1 + row_offset:1 + row_offset + h, 1 + col_offset:1 + col_offset + w]
    return jnp.where(jnp.isfinite(shifted), shifted, 0.0)


def legacy_receive_gather(lateral_flow, routing_weights, active_mask):
    """Perl's exact receiving-cell expression, not donor-weighted D-infinity.

    ``routing_weights[k, row, col]`` belongs to the receiver.  For example,
    the receiver's ``north`` weight multiplies lateral flow from ``[-1, 0]``.
    Border/NULL neighbour flow is zero before multiplication, exactly as each
    ``if(isnull(...), 0.0, ...)`` term in ``smr_buildup.pl``.
    """
    # GRASS MASK excludes donors as well as receivers. A finite backing-array
    # value outside support must not leak into an active neighboring cell.
    masked_flow = jnp.where(active_mask, lateral_flow, jnp.nan)
    gathered = sum(_neighbor(masked_flow, dr, dc) * routing_weights[k]
                   for k, (dr, dc) in enumerate(_OFFSETS))
    return jnp.where(active_mask, gathered, 0.0)
