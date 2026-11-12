"""Synthetic two-layer interface example. Run: python examples/plant_water.py."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import jax.numpy as jnp

from smr_jax.contracts import StaticInputs
from smr_jax.plant import (SurfaceWater, VegetationDemand, external_water_step,
                           require_realized_uptake, soil_to_plant)
from smr_jax.state import GridSpec, TwoLayerState


def synthetic_case(shape=(1, 1), depths=(30.0, 70.0), moisture=.20):
    """Explicit synthetic properties only; these are not defaults for real soil."""
    f = lambda v: jnp.full(shape, v, dtype=jnp.float64)
    a, b = depths
    static = StaticInputs(
        grid=GridSpec(*shape, 30.0, 0.0, 0.0), active_mask=jnp.ones(shape, dtype=bool),
        elevation_m=f(784.0), slope_percent=f(0.0), landuse=f(5.0), watershed=f(1.0),
        flowunits=f(1.0), routing_weights=jnp.zeros((8,) + shape, dtype=jnp.float64),
        sat_amt_cm=f(.45 * (a + b)), fieldcap_amt_cm=f(.28 * (a + b)),
        fieldcap_a_cm=f(.28 * a), fieldcap_b_cm=f(.28 * b), wiltpt_amt_cm=f(.10 * (a + b)),
        soil_depth_a_cm=f(a), soil_depth_b_cm=f(b), sat_mc_a=f(.45), sat_mc_b=f(.45),
        fieldcap_mc_a=f(.28), fieldcap_mc_b=f(.28), wiltpt_mc_a=f(.10),
        ksat_matrix_a=f(50.0), ksat_matrix_b=f(20.0), ksat_mpore_a=f(200.0), ksat_mpore_b=f(80.0),
        kfc_a=f(5.0), kfc_b=f(2.0), ksubsurface=f(.5), max_canopy_cm=f(.2), canopy_cover=f(.5),
        root_zone_cm=f(a + b), et_reduction_mc=f(.224), tmax_rain_c=f(1.5), tmin_snow_c=f(-.5),
        basin_ids=jnp.asarray([1.0]), basin_area_cells=jnp.asarray([float(shape[0] * shape[1])]),
        basin_reservoir_initial=jnp.asarray([0.0]), basin_reservoir_coeff=jnp.asarray([.01]))
    water = jnp.stack((f(moisture * a), f(moisture * b)))
    soil = TwoLayerState(water.sum(axis=0), water, static.active_mask)
    demand = VegetationDemand(jnp.stack((f(.08), f(.12))), jnp.stack((f(.04), f(0.0))), f(.02))
    surface = SurfaceWater(f(.4), f(.4), f(0.0), f(0.0), f(0.0), f(0.0), f(0.0))
    return soil, static, demand, surface


def main():
    soil, static, demand, surface = synthetic_case()
    # In production this demand comes from a candidate plant step. The checked
    # call refuses unmet demand by default; plant state can be committed only
    # after this transaction succeeds (or after explicit plant reconciliation).
    after, d = external_water_step(soil, static, demand, surface)
    require_realized_uptake(d)
    # Tomorrow's SoilToPlant kwargs use today's actual flows and an explicit
    # air-temperature proxy. No plant repository import is required.
    tomorrow = soil_to_plant(after, static, jnp.full(static.grid.shape, 15.0), d)
    print("backend: bounded legacy two-layer coupling (synthetic, not field validation)")
    for kind in ("transp", "soil_evap", "interception"):
        actual = d["interception_evap_cm"] if kind == "interception" else d[f"{kind}_cm"]
        print(f"{kind}: requested={float(d[f'requested_{kind}_cm'].sum()):.6f}, "
              f"actual={float(actual.sum()):.6f}, unmet={float(d[f'unmet_{kind}_cm'].sum()):.6f} cm")
    print(f"max layer residual: {float(jnp.max(jnp.abs(d['layer_residual_cm']))):.3e} cm")
    print(f"domain residual: {float(d['domain_residual_m3']):.3e} m3")
    print(f"tomorrow's packed water [layers, active cells]: {tomorrow['water_cm'].tolist()}")


if __name__ == "__main__":
    main()
