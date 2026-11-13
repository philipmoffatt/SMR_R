"""Legacy reporting-coupled watershed reservoir and baseflow recurrence."""

import jax.numpy as jnp
from ..state import BasinState


def basin_recurrence(basin, static, runoff_cm, percolation_cm):
    """Preserve Perl's asymmetric units and update order, per listed watershed."""
    def one(basin_id, volume, previous_baseflow, coefficient):
        mask = (static.watershed == basin_id) & static.active_mask
        runoff_sum, perc_sum = jnp.sum(jnp.where(mask, runoff_cm, 0.0)), jnp.sum(jnp.where(mask, percolation_cm, 0.0))
        runoff_cms = runoff_sum / 100.0 * static.grid.cell_area_m2 / (24.0 * 3600.0)
        perc_cms = perc_sum / 100.0 / (24.0 * 3600.0)
        volume = volume + perc_cms - previous_baseflow
        baseflow = coefficient * volume
        return volume, baseflow, baseflow + runoff_cms
    volume, baseflow, discharge = jnp.vectorize(one, signature="(),(),(),()->(),(),()")(
        static.basin_ids, basin.reservoir_volume, basin.baseflow_cms, static.basin_reservoir_coeff)
    return BasinState(volume, baseflow), discharge
