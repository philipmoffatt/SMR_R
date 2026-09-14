"""Read the frozen legacy ASCII bundle without invoking GRASS."""

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np

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
