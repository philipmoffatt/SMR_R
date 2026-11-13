"""Independent r.sun comparisons under changing date, sky and albedo.

Run in a NEW disposable GRASS project. Full capture stays outside Git; the
small deterministic representative-cell arrays may be versioned as fixtures.
"""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import grass.script as gs
import grass.script.array as ga

from prepare_jax_radiation import prepare_terrain, sky_coefficients, solar_components, digest


def capture(bundle, source, output):
    grid = prepare_terrain(bundle, output)
    output = Path(output)
    with ZipFile(bundle) as archive:
        path = output/"inputs/landuse.asc"
        path.write_bytes(archive.read("landuse.asc"))
    gs.run_command("r.in.gdal", flags="o", input=str(path), output="landuse", quiet=True)
    # Use the original source expression as the independent canopy oracle.
    import re
    expressions = re.findall(r'expression=(canopy_cover = [^"\r\n]+)', Path(source).read_text())
    if len(expressions) != 1:
        raise ValueError("missing/ambiguous source canopy expression")
    gs.mapcalc(expressions[0], overwrite=True, quiet=True)
    cover = np.asarray(ga.array("canopy_cover", null=np.nan))
    samples = None
    arrays = {name: [] for name in ("beam", "diffuse", "reflected_unit", "albedo", "canopy", "expected_q", "expected_reflected")}
    cases = []
    for doy in (1, 80, 172, 355):
        for cloud in (0., .5, 1.):
            linke = 3.+2.*cloud
            sky_coefficients(cloud)
            basis = solar_components(doy, linke)
            if samples is None:
                support = np.isfinite(cover) & np.isfinite(basis).all(axis=0)
                indices = np.flatnonzero(support)
                samples = indices[np.linspace(0, len(indices)-1, min(128, len(indices)), dtype=int)]
            for albedo in (0., .2, .7383, .95, 1.):
                gs.mapcalc(f"albedo = {albedo!r}", overwrite=True, quiet=True)
                beam, diffuse, reflected = solar_components(doy, linke, albedo_map="albedo")
                gs.mapcalc("q = (1.0-canopy_cover)*(1.0-albedo)*(rad_beam+rad_diffuse+rad_reflected)*3600.0/1000.0",
                           overwrite=True, quiet=True)
                actual = np.asarray(ga.array("q", null=np.nan))
                # Beam/diffuse independence must be demonstrated, not assumed.
                np.testing.assert_array_equal(beam, basis[0])
                np.testing.assert_array_equal(diffuse, basis[1])
                values = (*basis, np.full(cover.shape, albedo), cover, actual, reflected)
                for name, value in zip(arrays, values):
                    arrays[name].append(value.ravel()[samples])
                cases.append({"doy": doy, "cloud": cloud, "linke": linke, "albedo": albedo})
            print(f"captured doy={doy} cloud={cloud}", flush=True)
    artifact = output/"radiation.npz"
    np.savez_compressed(artifact, **{k: np.asarray(v, dtype=np.float64) for k, v in arrays.items()}, indices=samples)
    metadata = {"schema_version": 1, "captured_at": datetime.now(timezone.utc).isoformat(),
                "grid": grid, "cases": cases, "representatives": len(samples), "units": "Wh/m2/day; q in kJ/m2/day",
                "grass": gs.parse_command("g.version", flags="g"), "bundle_sha256": digest(bundle),
                "source_sha256": digest(source), "capture_script_sha256": digest(__file__),
                "preparation_script_sha256": digest(Path(__file__).with_name("prepare_jax_radiation.py")),
                "artifact_sha256": digest(artifact), "artifact_bytes": artifact.stat().st_size,
                "note": "Independent full r.sun at each albedo; basis contains only unit-albedo components."}
    (output/"radiation.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    capture(args.bundle, args.source, args.output)
