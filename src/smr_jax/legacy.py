"""Read the frozen legacy ASCII bundle without invoking GRASS."""

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import jax.numpy as jnp

from .contracts import StaticInputs
from .state import GridSpec


@dataclass(frozen=True)
class TwoLayerLegacyInputs:
    """Minimum static fields required for reference-compatible initialization."""

    grid: GridSpec
    sat_amt_cm: np.ndarray
    wiltpt_mc_a: np.ndarray
    fieldcap_mc_a: np.ndarray
    soil_depth_a_cm: np.ndarray
    soil_depth_b_cm: np.ndarray


def _read_esri_ascii(archive: ZipFile, member: str) -> tuple[GridSpec, np.ndarray]:
    with archive.open(member) as stream:
        lines = stream.read().decode("utf-8").splitlines()
    if len(lines) < 7:
        raise ValueError(f"{member} is not a complete ESRI ASCII raster")

    header = {}
    for line in lines[:6]:
        key, value = line.split(maxsplit=1)
        header[key.lower()] = value
    grid = GridSpec(
        nrows=int(header["nrows"]),
        ncols=int(header["ncols"]),
        cellsize_m=float(header["cellsize"]),
        xllcorner=float(header["xllcorner"]),
        yllcorner=float(header["yllcorner"]),
    )
    values = np.loadtxt(StringIO("\n".join(lines[6:])), dtype=np.float64)
    nodata = float(header.get("nodata_value", "nan"))
    if np.isfinite(nodata):
        values[values == nodata] = np.nan
    if values.shape != grid.shape:
        raise ValueError(f"{member} has shape {values.shape}, expected {grid.shape}")
    return grid, values


def load_two_layer_legacy_inputs(bundle_path: str | Path) -> TwoLayerLegacyInputs:
    """Load the frozen Mica Creek fields needed before flux translation begins."""

    with ZipFile(bundle_path) as archive:
        fields = {}
        grid = None
        for name in (
            "sat_amt",
            "wiltpt_mc_A",
            "fieldcap_mc_A",
            "soil_depth_A",
            "soil_depth_B",
        ):
            loaded_grid, fields[name] = _read_esri_ascii(archive, f"{name}.asc")
            if grid is None:
                grid = loaded_grid
            elif loaded_grid != grid:
                raise ValueError(f"Grid metadata differs for {name}.asc")

    assert grid is not None
    return TwoLayerLegacyInputs(
        grid=grid,
        sat_amt_cm=fields["sat_amt"],
        wiltpt_mc_a=fields["wiltpt_mc_A"],
        fieldcap_mc_a=fields["fieldcap_mc_A"],
        soil_depth_a_cm=fields["soil_depth_A"],
        soil_depth_b_cm=fields["soil_depth_B"],
    )


_RASTERS = (
    "el", "ETreduction_mc", "fieldcap_amt", "fieldcap_amt_A", "fieldcap_amt_B",
    "fieldcap_mc_A", "fieldcap_mc_B", "flowunits", "Kfc_A", "Kfc_B",
    "Ksat_matrix_A", "Ksat_matrix_B", "Ksat_mpore_A", "Ksat_mpore_B",
    "Ksubsurface", "landuse", "n", "ne", "nw", "s", "sat_amt", "sat_mc_A",
    "sat_mc_B", "se", "soil_depth", "soil_depth_A", "soil_depth_B", "sw", "w", "e",
    "watershed", "wiltpt_amt", "wiltpt_mc_A", "wiltpt_mc_B",
)


def _landcover_properties(landuse: np.ndarray, soil_depth: np.ndarray, soil_depth_a: np.ndarray,
                          fieldcap_amt: np.ndarray, fieldcap_a: np.ndarray,
                          wiltpt_amt: np.ndarray) -> tuple[np.ndarray, ...]:
    """The nested Perl land-use expressions, evaluated once as static inputs."""
    lu = landuse
    def by_lu(v6, v5, v4, v3, v2, v1, other=0.0):
        return np.select([lu == 6, lu == 5, lu == 4, lu == 3, lu == 2, lu == 1],
                         [v6, v5, v4, v3, v2, v1], default=other)
    max_canopy = by_lu(.147, .15, .15, .31, .1, 0.)
    canopy_cover = by_lu(.50, .60, .60, .80, 0., 0.)
    root_zone = by_lu(soil_depth, soil_depth_a, soil_depth, soil_depth, 0., 0.)
    et_reduction = by_lu(fieldcap_amt * .8 / soil_depth, fieldcap_a * .8 / soil_depth_a,
                         fieldcap_amt * .8 / soil_depth, fieldcap_amt * .7 / soil_depth,
                         0., 0.)
    tmax_rain = by_lu(1.8, -.5, 1.8, 3.1, -.5, -.5)
    tmin_snow = by_lu(-1., -3., -1., .9, -3., -3.)
    return max_canopy, canopy_cover, root_zone, et_reduction, tmax_rain, tmin_snow


def _reservoirs(archive: ZipFile) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows = []
    with archive.open("wshed_res_properties.ini") as stream:
        for line in stream.read().decode("utf-8").splitlines():
            if line.strip():
                rows.append(tuple(float(x) for x in line.split()))
    values = np.asarray(rows, dtype=np.float64)
    return values[:, 0], values[:, 1], values[:, 2], values[:, 3]


def _grass_default_slope_degrees(elevation_m: np.ndarray, cellsize_m: float) -> np.ndarray:
    """Horn 3x3 slope in degrees, matching GRASS ``r.slope.aspect`` defaults.

    The Perl engine requests the default degree output, then uses that numeric
    value in its historic ``slope / 100`` expression. A derivative needs a
    complete 3x3 finite neighbourhood; NULL centres and borders remain NULL
    rather than being filled, extrapolated, or silently treated as flat.

    This is an authorized emulation, not evidence of final parity. The routing
    verification gate compares its downstream effect before this is promoted
    from the replay implementation.
    """
    padded = np.pad(elevation_m, 1, constant_values=np.nan)
    z1, z2, z3 = padded[:-2, :-2], padded[:-2, 1:-1], padded[:-2, 2:]
    z4, z6 = padded[1:-1, :-2], padded[1:-1, 2:]
    z7, z8, z9 = padded[2:, :-2], padded[2:, 1:-1], padded[2:, 2:]
    valid = np.isfinite(np.stack((z1, z2, z3, z4, z6, z7, z8, z9))).all(axis=0)
    dz_dx = ((z3 + 2 * z6 + z9) - (z1 + 2 * z4 + z7)) / (8 * cellsize_m)
    dz_dy = ((z7 + 2 * z8 + z9) - (z1 + 2 * z2 + z3)) / (8 * cellsize_m)
    result = np.degrees(np.arctan(np.sqrt(dz_dx * dz_dx + dz_dy * dz_dy)))
    return np.where(valid, result, np.nan)


def load_legacy_bundle(bundle_path: str | Path) -> StaticInputs:
    """Load all frozen, model-ready fields without GRASS; never infer missing data."""
    with ZipFile(bundle_path) as archive:
        fields: dict[str, np.ndarray] = {}
        grid = None
        for name in _RASTERS:
            current, fields[name] = _read_esri_ascii(archive, f"{name}.asc")
            if grid is None:
                grid = current
            elif current != grid:
                raise ValueError(f"Grid metadata differs for {name}.asc")
        basin_ids, basin_area, basin_initial, basin_coeff = _reservoirs(archive)
    assert grid is not None
    slope = _grass_default_slope_degrees(fields["el"], grid.cellsize_m)
    # Integer watershed NULLs have already become NaN. Never use their original
    # uint sentinel as a basin ID.
    watershed = fields["watershed"]
    props = _landcover_properties(fields["landuse"], fields["soil_depth"], fields["soil_depth_A"],
                                  fields["fieldcap_amt"], fields["fieldcap_amt_A"], fields["wiltpt_amt"])
    ordered_weights = np.stack([fields[n] for n in ("n", "ne", "e", "se", "s", "sw", "w", "nw")])
    # GRASS r.mask and map algebra propagate NULL. The runnable engine uses the
    # intersection of every operand required before the soil loop, rather than
    # treating an input NULL as numerical zero. The CPU-06 foundation retains
    # its deliberately wider sat_amt-only initialization contract.
    active = np.isfinite(slope) & np.isfinite(watershed) & np.isfinite(ordered_weights).all(axis=0)
    for values in fields.values():
        active &= np.isfinite(values)
    f = lambda name: jnp.asarray(fields[name], dtype=jnp.float64)
    return StaticInputs(
        grid, jnp.asarray(active), f("el"), jnp.asarray(slope, dtype=jnp.float64),
        f("landuse"), jnp.asarray(watershed, dtype=jnp.float64), f("flowunits"), jnp.asarray(ordered_weights, dtype=jnp.float64),
        f("sat_amt"), f("fieldcap_amt"), f("fieldcap_amt_A"), f("fieldcap_amt_B"), f("wiltpt_amt"),
        f("soil_depth_A"), f("soil_depth_B"), f("sat_mc_A"), f("sat_mc_B"), f("fieldcap_mc_A"), f("fieldcap_mc_B"),
        f("wiltpt_mc_A"), f("Ksat_matrix_A"), f("Ksat_matrix_B"), f("Ksat_mpore_A"), f("Ksat_mpore_B"),
        f("Kfc_A"), f("Kfc_B"), f("Ksubsurface"), *(jnp.asarray(x, dtype=jnp.float64) for x in props),
        jnp.asarray(basin_ids), jnp.asarray(basin_area), jnp.asarray(basin_initial), jnp.asarray(basin_coeff),
    )
