"""Prepared r.sun terrain/weather components with evolving JAX snow albedo.

This decomposition uses GRASS during preparation only. Its float32 reflected
basis introduces bounded FCELL quantization; it is not a native r.sun port or
a completed shortwave replay. Units entering the kernel are Wh/m2/day; the
legacy absorbed shortwave is kJ/m2/day.
"""

from hashlib import sha256
import json
from pathlib import Path

import jax.numpy as jnp
import numpy as np

from .contracts import RadiationComponents
from .state import GridSpec


def _fcell(value):
    return value.astype(jnp.float32).astype(jnp.float64)


def shortwave_from_components(components, albedo, canopy_cover):
    """Apply simulated albedo and preserve source FCELL input/output boundaries.

    r.sun reads its ground albedo as FCELL. Its beam/diffuse/reflected outputs
    and their raster-only additions are also FCELL. The outside absorption
    factor uses the original DCELL albedo. Do not clip or silently fill NULLs.
    """
    reflected = _fcell(_fcell(albedo) * components.reflected_unit)
    incident = _fcell(_fcell(components.beam + components.diffuse) + reflected)
    return (1.-canopy_cover)*(1.-albedo)*incident*3600./1000.


def radiation_error_bound(components, albedo, canopy_cover):
    """Conservative forward-error bound from FCELL unit roundoff, not data fitting.

    Eight float32 epsilons cover rounding of the reflected basis, scaling and
    both incident additions in the two evaluation paths. Float64 algebra gets
    the existing 512-epsilon budget. This bounds radiation, not model states.
    """
    scale = np.abs((1.-np.asarray(canopy_cover))*(1.-np.asarray(albedo)))*3.6
    incident = sum(np.abs(np.asarray(v)) for v in
                   (components.beam, components.diffuse, components.reflected_unit))
    return 8.*np.finfo(np.float32).eps*scale*incident + 512.*np.finfo(np.float64).eps*np.maximum(1., scale*incident)


class PreparedRadiationProvider:
    """Validate an explicit terrain/forcing cache; never fall back to replay."""

    def __init__(self, directory, grid, bundle_path, *, required_mask=None):
        directory = Path(directory)
        metadata = json.loads((directory/"manifest.json").read_text(encoding="utf-8"))
        if (metadata.get("schema_version") != 2 or metadata.get("provider") != "grass-rsun-albedo-decomposition"
            or metadata.get("albedo_support") != "finite-landuse-and-elevation"
            or metadata.get("units") != "Wh/m2/day" or metadata.get("albedo_basis") != 1.
            or metadata.get("r_sun_step_hours") != .5 or metadata.get("terrain_shadowing") is not True):
            raise ValueError("unsupported radiation cache conventions")
        prepared_grid = GridSpec(**metadata["grid"])
        # GRASS serializes projected origins to fewer decimal places than the
        # ASCII bundle. Use the existing TIFF reader's 1e-8-cell alignment floor.
        aligned_origin = np.allclose(
            ((prepared_grid.xllcorner-grid.xllcorner)/grid.cellsize_m,
             (prepared_grid.yllcorner-grid.yllcorner)/grid.cellsize_m),
            (0., 0.), rtol=0., atol=1e-8)
        if (prepared_grid.shape != grid.shape or prepared_grid.cellsize_m != grid.cellsize_m or not aligned_origin):
            raise ValueError("radiation cache grid differs from model grid")
        if sha256(Path(bundle_path).read_bytes()).hexdigest() != metadata["bundle_sha256"]:
            raise ValueError("radiation cache terrain/bundle provenance differs")
        artifact = directory/"components.npz"
        if sha256(artifact.read_bytes()).hexdigest() != metadata["components_sha256"]:
            raise ValueError("radiation component checksum differs")
        cases = metadata["cases"]
        dates = [c["date"] for c in cases]
        if not cases or len(set(dates)) != len(dates):
            raise ValueError("radiation cache needs unique nonempty dates")
        with np.load(artifact, allow_pickle=False) as archive:
            if set(archive.files) != {"beam", "diffuse", "reflected_unit"}:
                raise ValueError("radiation cache component inventory differs")
            self.arrays = {key: archive[key] for key in archive.files}
        support = None
        for value in self.arrays.values():
            if value.shape != (len(cases),)+grid.shape or value.dtype != np.float32:
                raise ValueError("radiation components need complete [date,y,x] FCELL arrays")
            finite = np.isfinite(value)
            if not finite.any(axis=(1, 2)).all() or np.any(np.isinf(value)) or np.any(value[finite] < 0.):
                raise ValueError("radiation components require finite nonnegative supported values")
            if support is not None and not np.array_equal(support, finite):
                raise ValueError("radiation component supports differ")
            support = finite
        if required_mask is not None:
            required = np.asarray(required_mask)
            if required.dtype != np.bool_ or required.shape != grid.shape or not required.any():
                raise ValueError("required radiation support must be a nonempty boolean model mask")
            if np.any(required[None] & ~support):
                raise ValueError("radiation cache is missing required model support")
        self.metadata, self.cases = metadata, {c["date"]: (i, c) for i, c in enumerate(cases)}

    def components(self, date, *, cloud, linke, doy):
        if date not in self.cases:
            raise ValueError(f"no prepared radiation for {date}; prepare the requested forcing dates")
        index, case = self.cases[date]
        if any(case[key] != value for key, value in (("cloud", cloud), ("linke", linke), ("doy", doy))):
            raise ValueError(f"{date}: prepared radiation weather differs; regenerate components")
        return RadiationComponents(*(jnp.asarray(self.arrays[key][index], dtype=jnp.float64)
                                     for key in ("beam", "diffuse", "reflected_unit")))
