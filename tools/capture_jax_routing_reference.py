"""Run inside a NEW temporary GRASS EPSG:32611 location to capture VERIFY-03.

Example: grass --tmp-location EPSG:32611 --exec python this_script.py
Outputs are refused if already present. No existing SMR runtime is touched.
"""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np
import grass.script as gs
import grass.script.array as ga


def capture(root, output):
    binary, manifest = output.with_suffix(".npz"), output.with_suffix(".json")
    if binary.exists() or manifest.exists():
        raise ValueError("refusing to overwrite existing reference capture")
    reference = root / "reference/mica-creek-jan2015"
    source_path = root / "R/smr_buildup.pl"
    source = source_path.read_text(encoding="utf-8")
    artifacts = json.loads((reference / "manifest.json").read_text(encoding="utf-8-sig"))["artifacts"]
    for item in artifacts:
        path = reference / item["path"]
        if path.stat().st_size != item["bytes"] or sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"corrupt frozen artifact: {path.name}")
    aliases = dict(n="north", ne="northeast", e="east", se="southeast", s="south", sw="southwest", w="west", nw="northwest")
    names = ["el", "watershed", "landuse", "wiltpt_mc_A", "wiltpt_mc_B", "soil_depth_A", "soil_depth_B",
             "soil_depth", "fieldcap_amt", "fieldcap_amt_A", "fieldcap_amt_B", "sat_mc_A", "sat_mc_B",
             "Ksat_matrix_A", "Ksat_matrix_B", "Ksat_mpore_A", "Ksat_mpore_B", "Kfc_A", "Kfc_B", "flowunits", *aliases]
    arrays, metadata = {}, {}
    read = lambda name: np.asarray(ga.array(name, null=np.nan), dtype=np.float64).copy()
    with TemporaryDirectory(prefix="smr-verify03-inputs-") as directory:
        directory = Path(directory)
        with ZipFile(reference / "model_inputs.zip") as bundle:
            for name in names:
                path = directory / f"{name}.asc"
                path.write_bytes(bundle.read(f"{name}.asc"))
                gs.run_command("r.in.gdal", flags="o", input=str(path), output=aliases.get(name, name), quiet=True)
            metadata["reservoir_ini"] = bundle.read("wshed_res_properties.ini").decode("utf-8")
        gs.run_command("g.region", raster="watershed")
        gs.mapcalc("watershed = if(watershed==1.0,1,1)", overwrite=True, quiet=True)
        gs.run_command("r.mask", raster="watershed", quiet=True)
        gs.run_command("r.slope.aspect", elevation="el", slope="slope", aspect="aspect", quiet=True)
        arrays["slope"] = read("slope")
        properties = source.split("# Land cover maps", 1)[1].split("# Initial conditions", 1)[0]
        for expression in re.findall(r'expression=([^"\r\n]+)', properties):
            target = expression.split(" = ", 1)[0]
            if target in ("root_zone", "wiltpt_amt", "ETreduction_mc"):
                gs.mapcalc(expression, overwrite=True, quiet=True)
                arrays[target] = read(target)
        metadata["raster_types"] = {name: gs.raster_info(name)["datatype"] for name in
                                    ("el", "slope", "root_zone", "wiltpt_amt", "ETreduction_mc")}
        block = source.split("# Lateral flow (Darcy + Bresler)", 1)[1].split("# Re-distribute between A and B", 1)[0]
        expressions = re.findall(r'expression=([^"\r\n]+)', block)
        if len(expressions) != 10:
            raise ValueError("changed source lateral block")
        with ZipFile(reference / "two_day_internal_trace.zip") as trace:
            for date in ("2015-01-01", "2015-01-02"):
                for name in ("storage_amt", "storage_amt_A", "storage_amt_B", "saturation"):
                    path = directory / f"{name}.tif"
                    path.write_bytes(trace.read(f"trace/{date}--{name}.tif"))
                    gs.run_command("r.in.gdal", flags="o", input=str(path), output=name, overwrite=True, quiet=True)
                support = np.isfinite(read("saturation"))
                if "support" not in arrays:
                    arrays["support"] = support
                    lu = read("landuse")
                    points = []
                    for code in np.unique(lu[support]):
                        flat = np.flatnonzero(support & (lu == code))
                        points.extend(flat[np.unique(np.linspace(0, len(flat)-1, min(64, len(flat)), dtype=int))])
                    arrays["sample_flat"] = np.asarray(sorted(set(points)), dtype=np.int64)
                elif not np.array_equal(arrays["support"], support):
                    raise ValueError("frozen support changed between days")
                for expression in expressions:
                    gs.mapcalc(expression.replace("$temp_time_step", "6").replace("$gridsize", "30"), overwrite=True, quiet=True)
                for name in ("effK_A", "effK_B", "lateral_flow_A", "lateral_flow_B", "lateral_flow",
                             "storage_amt", "storage_amt_A", "storage_amt_B"):
                    values = read(name)
                    if not np.isfinite(values[support]).all():
                        raise ValueError(f"nonfinite captured lateral field: {date} {name}")
                    arrays[f"{date}_{name}"] = values.ravel()[arrays["sample_flat"]]
    binary.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(binary, **arrays)
    metadata.update(schema_version=1, captured_at=datetime.now(timezone.utc).isoformat(),
                    grass=gs.parse_command("g.version", flags="g"), projection="EPSG:32611",
                    source_sha256=sha256(source_path.read_bytes()).hexdigest(),
                    capture_script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                    frozen_artifacts={item["path"]: item["sha256"] for item in artifacts},
                    binary_sha256=sha256(binary.read_bytes()).hexdigest(), binary_bytes=binary.stat().st_size,
                    required_cells=int(arrays["support"].sum()), representative_cells=len(arrays["sample_flat"]),
                    note="Fresh GRASS slope, static properties and isolated lateral stages from frozen daily starting states; not integrated daily truth.")
    manifest.write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"bytes": binary.stat().st_size, "cells": metadata["required_cells"],
                      "representatives": metadata["representative_cells"], "raster_types": metadata["raster_types"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    capture(args.root.resolve(), args.output.resolve())
