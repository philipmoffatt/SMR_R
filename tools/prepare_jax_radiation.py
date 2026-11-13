"""Prepare albedo-independent r.sun components in an EMPTY disposable GRASS mapset.

Run with grass --tmp-project <CRS> --exec python this_script ... . Output must
be a new directory outside Git. The cache contains beam, diffuse and reflected
irradiation at unit albedo, never completed q.srad or evolving snow states.
"""

import argparse
from datetime import date
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import grass.script as gs
import grass.script.array as ga


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def prepare_terrain(bundle, output):
    if Path(output).exists():
        raise ValueError("output directory must be NEW; refusing overwrite")
    if gs.list_strings(type="raster", mapset="."):
        raise ValueError("preparation requires an empty disposable GRASS mapset")
    output = Path(output)
    (output/"inputs").mkdir(parents=True)
    with ZipFile(bundle) as archive:
        for name in ("el", "watershed", "landuse"):
            path = output/"inputs"/(name+".asc")
            path.write_bytes(archive.read(name+".asc"))
            gs.run_command("r.in.gdal", flags="o", input=str(path), output=name, quiet=True)
    gs.run_command("g.region", raster="watershed")
    gs.run_command("r.mask", raster="watershed", quiet=True)
    gs.run_command("r.slope.aspect", elevation="el", slope="slope", aspect="aspect", quiet=True)
    # Albedo NULLs remove terrain from r.sun shadow casting as well as output.
    gs.mapcalc("unit_albedo = if(isnull(landuse)||isnull(el),null(),1.0)", overwrite=True, quiet=True)
    region = gs.region()
    if region["nsres"] != region["ewres"]:
        raise ValueError("SMR radiation requires square cells")
    return {"nrows": region["rows"], "ncols": region["cols"], "cellsize_m": region["nsres"],
            "xllcorner": region["w"], "yllcorner": region["s"]}


def sky_coefficients(cloud):
    if not np.isfinite(cloud) or not 0. <= cloud <= 1.:
        raise ValueError("cloud must be finite in [0,1]")
    value = repr(float(cloud))
    # Preserve the literal legacy real-sky equations, including coefdh > 1.
    gs.mapcalc(f"t_diff = if({value}==1.0,0.1,0.75*(1-{value})*(1-exp(-0.6*{value}/((0.75-0.4)*(1-{value})))))*watershed",
               overwrite=True, quiet=True)
    gs.mapcalc("coefdh = t_diff/0.10", overwrite=True, quiet=True)
    gs.mapcalc(f"coefbh = max(0.0,(1-{value})-t_diff/0.75)", overwrite=True, quiet=True)


def solar_components(doy, linke, *, albedo_map=None):
    if not 1 <= doy <= 365 or not np.isfinite(linke) or linke <= 0:
        raise ValueError("r.sun requires day 1..365 and positive finite Linke turbidity")
    options = {"albedo": albedo_map} if albedo_map else {"albedo": "unit_albedo"}
    gs.run_command("r.sun", elevation="el", slope="slope", aspect="aspect",
                   day=doy, linke_value=linke, coeff_bh="coefbh", coeff_dh="coefdh",
                   beam_rad="rad_beam", diff_rad="rad_diffuse", refl_rad="rad_reflected",
                   nprocs=1, overwrite=True, quiet=True, **options)
    return tuple(np.asarray(ga.array(name, null=np.nan)).astype(np.float32)
                 for name in ("rad_beam", "rad_diffuse", "rad_reflected"))


def prepare(bundle, weather, output, limit=None):
    grid = prepare_terrain(bundle, output)
    output = Path(output)
    cases, data = [], []
    rows = [line.split() for line in Path(weather).read_text().splitlines() if line.strip()]
    if limit is not None:
        rows = rows[:limit]
    for row in rows:
        when = date.fromisoformat(row[0])
        doy, cloud, linke = int(row[4]), float(row[16]), float(row[15])
        if doy != when.timetuple().tm_yday or row[0] in {c["date"] for c in cases}:
            raise ValueError("inconsistent day-of-year or duplicate weather date")
        sky_coefficients(cloud)
        data.append(solar_components(doy, linke))
        cases.append({"date": row[0], "doy": doy, "cloud": cloud, "linke": linke})
        print(f"prepared {row[0]}", flush=True)
    if not cases:
        raise ValueError("weather contains no preparation dates")
    arrays = np.asarray(data, dtype=np.float32)
    artifact = output/"components.npz"
    np.savez_compressed(artifact, beam=arrays[:, 0], diffuse=arrays[:, 1], reflected_unit=arrays[:, 2])
    metadata = {"schema_version": 2, "albedo_support": "finite-landuse-and-elevation", "provider": "grass-rsun-albedo-decomposition", "grid": grid,
                "units": "Wh/m2/day", "cases": cases, "grass": gs.parse_command("g.version", flags="g"),
                "projection": gs.parse_command("g.proj", flags="g"),
                "bundle_sha256": digest(bundle), "weather_sha256": digest(weather),
                "preparation_script_sha256": digest(__file__), "components_sha256": digest(artifact),
                "albedo_basis": 1., "r_sun_step_hours": .5, "terrain_shadowing": True,
                "note": "Only terrain/date/weather components; model albedo is evaluated inside JAX."}
    (output/"manifest.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--weather", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    prepare(args.bundle, args.weather, args.output, args.limit)


if __name__ == "__main__":
    main()
