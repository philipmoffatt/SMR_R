"""Export a bounded exact-commit GPU execution snapshot; never overwrite outputs."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tarfile


PATHS = ("src/smr_jax", "tests", "examples", "reference/mica-creek-jan2015",
         "R/smr_buildup.pl", "tools/capture_jax_integrated_reference.py", "pyproject.toml")


def package(repository, revision, output):
    repository, output = Path(repository).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError("snapshot output must be a NEW directory")
    revision = subprocess.check_output(["git", "-C", str(repository), "rev-parse", "--verify", revision+"^{commit}"], text=True).strip()
    subprocess.run(["git", "-C", str(repository), "cat-file", "-e", revision+":src/smr_jax/verify_gpu.py"], check=True)
    output.mkdir(parents=True)
    archive = output/"source.tar"
    subprocess.run(["git", "-C", str(repository), "archive", "--format=tar", "--output", str(archive), revision, "--", *PATHS], check=True)
    hashes = {}
    with tarfile.open(archive) as bundle:
        for member in bundle:
            relative = Path(member.name)
            if relative.is_absolute() or ".." in relative.parts or not (member.isfile() or member.isdir()):
                raise ValueError("unexpected source snapshot member")
            if member.isfile():
                hashes[member.name] = sha256(bundle.extractfile(member).read()).hexdigest()
    metadata = {"schema_version": 1, "source_revision": revision, "archive_sha256": sha256(archive.read_bytes()).hexdigest(),
                "files": hashes, "note": "Execution replica only; canonical edits remain in the Windows packet worktree."}
    (output/"snapshot.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    return metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--revision", default="HEAD")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = package(args.repository, args.revision, args.output)
    print(json.dumps({k: v for k, v in result.items() if k != "files"}, indent=2))
