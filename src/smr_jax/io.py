"""Host-side frozen-bundle, weather, and replay-radiation readers."""

from dataclasses import fields
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import jax
import jax.numpy as jnp
import numpy as np

from .contracts import DailyForcing
from .state import GridSpec


class FrozenTraceProvider:
    """Read completed legacy trace maps and align their exported extent to the model grid."""

    def __init__(self, trace_zip: str | Path, target_grid: GridSpec):
        self.archive_path = Path(trace_zip)
        self.target_grid = target_grid

    def raster(self, date: str, map_name: str) -> np.ndarray:
        try:
            import tifffile
        except ImportError as exc:  # pragma: no cover - environment diagnosis
            raise RuntimeError("Replay radiation requires the project dependency tifffile") from exc
        member = f"trace/{date}--{map_name.replace('.', '_')}.tif"
        with ZipFile(self.archive_path) as archive:
            try:
                # The frozen archive was built on Windows. zip members can
                # retain backslashes even though archive listings normalize
                # them, so resolve by a normalized name rather than assuming
                # a host-path separator.
                actual = next(name for name in archive.namelist()
                              if name.replace("\\", "/") == member)
                with tifffile.TiffFile(BytesIO(archive.read(actual))) as image:
                    page = image.pages[0]
                    values = page.asarray().astype(np.float64)
                    scale = page.tags["ModelPixelScaleTag"].value
                    tiepoint = page.tags["ModelTiepointTag"].value
                if not np.isclose(scale[0], self.target_grid.cellsize_m) or not np.isclose(scale[1], self.target_grid.cellsize_m):
                    raise ValueError(f"Replay resolution {scale[:2]} differs from frozen grid cell size {self.target_grid.cellsize_m}")
                col = round((tiepoint[3] - self.target_grid.xllcorner) / self.target_grid.cellsize_m)
                grid_top = self.target_grid.yllcorner + self.target_grid.nrows * self.target_grid.cellsize_m
                row = round((grid_top - tiepoint[4]) / self.target_grid.cellsize_m)
                if row < 0 or col < 0 or row + values.shape[0] > self.target_grid.nrows or col + values.shape[1] > self.target_grid.ncols:
                    raise ValueError(f"Replay field {actual} does not fit frozen grid at row={row}, col={col}")
                aligned = np.full(self.target_grid.shape, np.nan, dtype=np.float64)
                aligned[row:row + values.shape[0], col:col + values.shape[1]] = values
                return aligned
            except KeyError as exc:
                raise KeyError(f"No frozen trace field for {map_name} on {date}") from exc


class ReplayRadiationProvider(FrozenTraceProvider):
    """Read completed legacy ``q.srad`` maps; this is deliberately not general radiation."""

    def shortwave(self, date: str) -> np.ndarray:
        return self.raster(date, "q.srad")


def _weather_rows(path: str | Path):
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.strip():
            yield line.split()


def load_daily_forcings(weather_path: str | Path, radiation: ReplayRadiationProvider, limit: int | None = None) -> DailyForcing:
    """Convert the space-delimited frozen GridMET forcing into stacked JAX input."""
    days = []
    for number, row in enumerate(_weather_rows(weather_path)):
        if limit is not None and number >= limit:
            break
        date = row[0]
        scalar = lambda index: jnp.asarray(float(row[index]), dtype=jnp.float64)
        coeff = jnp.asarray([float(row[i]) for i in (17, 18, 19, 20, 21, 22)], dtype=jnp.float64)
        days.append(DailyForcing(scalar(1), scalar(2), scalar(3), scalar(4), scalar(7), scalar(8), scalar(9), scalar(10),
                                 jnp.asarray([float(row[i]) for i in (11, 12, 13, 14)], dtype=jnp.float64),
                                 scalar(16), scalar(24), scalar(28), coeff,
                                 jnp.asarray(radiation.shortwave(date), dtype=jnp.float64)))
    if not days:
        raise ValueError("No forcing rows were loaded")
    return jax.tree_util.tree_map(lambda *parts: jnp.stack(parts), *days)
