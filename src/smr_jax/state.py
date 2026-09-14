"""JAX state contracts; no hydrologic flux equations live in this module."""

from dataclasses import dataclass

import jax
import jax.numpy as jnp

# Reference parity begins in float64. This must be set before arrays are made.
jax.config.update("jax_enable_x64", True)


@dataclass(frozen=True)
class GridSpec:
    """Static raster geometry shared by every dynamic state field."""

    nrows: int
    ncols: int
    cellsize_m: float
    xllcorner: float
    yllcorner: float

    @property
    def shape(self) -> tuple[int, int]:
        return (self.nrows, self.ncols)

    @property
    def cell_area_m2(self) -> float:
        return self.cellsize_m * self.cellsize_m


@dataclass(frozen=True)
class TwoLayerState:
    """Dynamic soil state with a leading layer axis and explicit legacy total."""

    total_water_cm: jax.Array
    layer_water_cm: jax.Array
    active_mask: jax.Array

    def __post_init__(self) -> None:
        if self.layer_water_cm.ndim != 3 or self.layer_water_cm.shape[0] != 2:
            raise ValueError("TwoLayerState requires layer_water_cm with shape [2, y, x]")
        if self.total_water_cm.shape != self.layer_water_cm.shape[1:]:
            raise ValueError("total_water_cm must match the raster shape of layer_water_cm")
        if self.active_mask.shape != self.total_water_cm.shape:
            raise ValueError("active_mask must match total_water_cm")

    @property
    def n_layers(self) -> int:
        return int(self.layer_water_cm.shape[0])


jax.tree_util.register_dataclass(
    TwoLayerState,
    data_fields=["total_water_cm", "layer_water_cm", "active_mask"],
    meta_fields=[],
)
