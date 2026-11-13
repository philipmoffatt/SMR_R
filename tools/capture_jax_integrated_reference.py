"""Capture a full thirty-day reference inside a NEW temporary GRASS project.

Only a temporary Perl copy receives additional state exports. Its original
scientific expressions remain unchanged. Output must be a new directory outside
Git; the large lossless trace is rebuildable from this script and frozen inputs.
"""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
from zipfile import ZipFile, ZIP_DEFLATED

import grass.script as gs


EXTRA_MAPS = ("canopy_storage_amt", "ice.content", "liquid.water", "snow.age", "u.total", "u.surface",
              "tsnow_surf", "swe.yesterday", "road_runoff", "q.rain.ground", "runoff_2015", "Psat_2015")


def capture(root, output, perl):
    if output.exists():
        raise ValueError("capture requires a NEW output directory; refusing to overwrite")
    reference = root/"reference/mica-creek-jan2015"
    source_path = root/"R/smr_buildup.pl"
    source_bytes = source_path.read_bytes()
    source = source_bytes.decode("utf-8")
    inventory = json.loads((reference/"manifest.json").read_text(encoding="utf-8-sig"))["artifacts"]
    for item in inventory:
        path = reference/item["path"]
        if path.stat().st_size != item["bytes"] or sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"corrupt frozen artifact: {path.name}")
    output.mkdir(parents=True)
    inputs, trace = output/"inputs", output/"trace"
    inputs.mkdir()
    trace.mkdir()
    aliases = dict(n="north", ne="northeast", e="east", se="southeast", s="south", sw="southwest", w="west", nw="northwest")
    names = ["watershed", "el", "strms_30m", "landuse", "Ksat_matrix_A", "Ksat_matrix_B", "Ksat_mpore_A", "Ksat_mpore_B",
             "Kfc_A", "Kfc_B", "sat_mc_A", "sat_mc_B", "sat_amt", "wiltpt_mc_A", "wiltpt_mc_B", "wiltpt_amt",
             "fieldcap_mc_A", "fieldcap_mc_B", "fieldcap_amt_A", "fieldcap_amt_B", "fieldcap_amt", "soil_depth_A",
             "soil_depth_B", "soil_depth", "ETreduction_mc", "Ksubsurface", "flowunits", *aliases]
    with ZipFile(reference/"model_inputs.zip") as bundle:
        for name in names:
            path = inputs/f"{name}.asc"
            path.write_bytes(bundle.read(f"{name}.asc"))
            gs.run_command("r.in.gdal", flags="o", input=str(path), output=aliases.get(name, name), quiet=True)
        (inputs/"wshed_res_properties.ini").write_bytes(bundle.read("wshed_res_properties.ini"))
    # Add outputs only; preserve the complete sequence of scientific commands.
    match = re.search(r"my @trace_maps = qw\((.*?)\);", source, re.S)
    if match is None:
        raise ValueError("missing original trace inventory")
    trace_names = match[1].split()+list(EXTRA_MAPS)
    instrumented = source[:match.start(1)]+"\n    "+" ".join(trace_names)+"\n"+source[match.end(1):]
    initial_export = '''
if ($trace_enabled) {
    foreach my $verify_map (qw(storage_amt storage_amt_A)) {
        system('r.out.gdal', "input=$verify_map", "output=$trace_dir/initial--$verify_map.tif",
               'format=GTiff', 'createopt=COMPRESS=DEFLATE', '--overwrite', '--quiet');
    }
}
'''
    instrumented = instrumented.replace("# Weather loop", initial_export+"\n# Weather loop", 1)
    # Four first-day soil checkpoints expose integration ordering errors without
    # requiring changes to the production engine or its original frozen source.
    helper = '''
sub verify_checkpoint {
    my ($label, $step, $date) = @_;
    return unless $date eq '2015-01-01';
    foreach my $verify_map (qw(storage_amt storage_amt_A storage_amt_B)) {
        system('r.out.gdal', "input=$verify_map", "output=$trace_dir/checkpoint-$step-$label--$verify_map.tif",
               'format=GTiff', 'createopt=COMPRESS=DEFLATE', '--overwrite', '--quiet');
    }
}
'''
    instrumented = instrumented.replace("# Weather loop", helper+"\n# Weather loop", 1)
    loops = instrumented.split("foreach my $hrly_tmp (@hourly_tmp_array) {")
    if len(loops) != 3:
        raise ValueError("changed source substep loops")
    instrumented = loops[0]+"foreach my $hrly_tmp (@hourly_tmp_array) {"+loops[1]+"my $verify_step = 0;\n    foreach my $hrly_tmp (@hourly_tmp_array) {"+loops[2]
    for marker, label in (("# Lateral flow (Darcy + Bresler)", "allocation"),
                           ("# Re-distribute between A and B", "lateral"),
                           ("# Percolation", "et"), ("# Saturation-excess runoff", "percolation")):
        instrumented = instrumented.replace(marker, f'verify_checkpoint("{label}", $verify_step, $date);\n        '+marker, 1)
    instrumented = instrumented.replace('    print `r.mapcalc "expression=runoff_daily_flow = runoff_daily_flow+road_runoff"',
        '    print `r.mapcalc "expression=runoff_daily_flow = runoff_daily_flow+road_runoff"', 1)
    instrumented = instrumented.replace('        print `r.mapcalc "expression=runoff_$year = runoff_flow+runoff_$year" --overwrite`;',
        '        print `r.mapcalc "expression=runoff_$year = runoff_flow+runoff_$year" --overwrite`;\n        $verify_step++;', 1)
    script = output/"capture-engine.pl"
    script.write_text(instrumented, encoding="utf-8", newline="\n")
    environment = dict(os.environ, SMR_TRACE_DIR=str(trace), SMR_OUTPUT_MODE="full", SMR_RSUN_NPROCS="4")
    print(f"Starting 30-day frozen Perl capture in {output}", flush=True)
    with (output/"capture.log").open("w", encoding="utf-8") as log:
        subprocess.run([str(perl), str(script), str(reference/"gridmet_weather_2015-01.csv"), str(output/"daily.csv")],
                       cwd=inputs, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
    archive_path = output/"daily_trace.zip"
    with ZipFile(archive_path, "x", compression=ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(trace.glob("*.tif")):
            archive.write(path, f"trace/{path.name}")
        archive.write(output/"daily.csv", "daily.csv")
    metadata = dict(schema_version=1, captured_at=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha256(source_bytes).hexdigest(), instrumented_source_sha256=sha256(script.read_bytes()).hexdigest(),
        capture_script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), grass=gs.parse_command("g.version", flags="g"),
        frozen_artifacts={item["path"]: item["sha256"] for item in inventory}, days=30, daily_maps=trace_names,
        artifacts={path.name: dict(bytes=path.stat().st_size, sha256=sha256(path.read_bytes()).hexdigest())
                   for path in (archive_path, output/"daily.csv", script)},
        note="Additional exports only, temporary GRASS project; original frozen source and artifacts unchanged. Validate anchors before using this extension.")
    (output/"manifest.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"trace_bytes": archive_path.stat().st_size, "daily_maps": len(trace_names), "days": 30}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--perl", type=Path, default=Path("C:/Strawberry/perl/bin/perl.exe"))
    args = parser.parse_args()
    capture(args.root.resolve(), args.output.resolve(), args.perl.resolve())
