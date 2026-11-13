"""Build/load an explicit native exp option for the frozen Windows CPU gate.

Ordinary model runs keep JAX math. This optional typed FFI executes C exp on XLA
CPU buffers, without Python callbacks, GRASS, reference states or global XLA
flags. One immutable library binding is allowed per process. The library and
its build record are generated at an explicitly supplied location outside Git.
"""

import argparse
from copy import deepcopy
import ctypes
from functools import partial
from hashlib import sha256
import json
from pathlib import Path
import platform
import subprocess
import sys

import jax
import jax.numpy as jnp


_SOURCE = Path(__file__).with_name("legacy_exp.cc")
_binding = None


def build_library(compiler, output):
    """Build with an operator-supplied MinGW compiler; refuse existing outputs."""
    output = Path(output).resolve()
    record = output.with_suffix(output.suffix+".json")
    if output.exists() or record.exists():
        raise ValueError("native math build requires new library and record paths")
    if sys.platform != "win32" or platform.machine().lower() not in ("amd64", "x86_64"):
        raise ValueError("native legacy math is qualified only for Windows x64")
    machine = subprocess.check_output([str(compiler), "-dumpmachine"], text=True).strip()
    if machine != "x86_64-w64-mingw32":
        raise ValueError("native legacy math requires an x86_64-w64-mingw32 compiler")
    command = [str(compiler), "-std=c++17", "-shared", "-static", "-O2", "-fno-fast-math",
               "-ffp-contract=off", "-fno-builtin-exp",
               "-I", jax.ffi.include_dir(), str(_SOURCE), "-o", str(output)]
    output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(command, check=True)
    metadata = {"schema_version": 1, "profile": "windows-mingw64-scalar-exp",
                "jax_version": jax.__version__, "source_sha256": sha256(_SOURCE.read_bytes()).hexdigest(),
                "library_sha256": sha256(output.read_bytes()).hexdigest(),
                "compiler_version": subprocess.check_output([str(compiler), "--version"], text=True).splitlines()[0],
                "command": command}
    with record.open("x", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")
    return metadata


def load_library(path):
    """Verify build provenance and register one immutable CPU-only FFI target."""
    global _binding
    path = Path(path).resolve()
    metadata = json.loads(path.with_suffix(path.suffix+".json").read_text(encoding="utf-8"))
    if (metadata.get("schema_version") != 1 or metadata.get("profile") != "windows-mingw64-scalar-exp"
        or metadata.get("jax_version") != jax.__version__
        or metadata.get("source_sha256") != sha256(_SOURCE.read_bytes()).hexdigest()
        or metadata.get("library_sha256") != sha256(path.read_bytes()).hexdigest()):
        raise ValueError("native math library/build provenance differs; rebuild with this JAX and source")
    if _binding is not None:
        if _binding[0] != path or _binding[1] != metadata:
            raise ValueError("native math library is already bound; use a fresh process to change it")
        return deepcopy(metadata)
    library = ctypes.CDLL(str(path))
    target = "smr_legacy_exp_"+metadata["library_sha256"]
    jax.ffi.register_ffi_target(target, jax.ffi.pycapsule(library.SmrLegacyExp), platform="cpu")
    _binding = (path, metadata, library, target)  # Retain the DLL for compiled executable lifetimes.
    return deepcopy(metadata)


@jax.custom_jvp
def native_exp(value):
    """Elementwise native exp, available only after explicit library loading."""
    if value.dtype != jnp.float64:
        raise ValueError("native legacy exp requires float64")
    if _binding is None:
        raise ValueError("load a native legacy math library before selecting native_exp")
    return jax.ffi.ffi_call(_binding[3], jax.ShapeDtypeStruct(value.shape, value.dtype),
                            vmap_method="expand_dims")(value)


@native_exp.defjvp
def _native_exp_jvp(primals, tangents):
    result = native_exp(primals[0])
    return result, result*tangents[0]


def compile_native_cpu(function, *, static_argnames=()):
    """Bind the explicit exponential option and the existing CPU compilation profile."""
    from .fluxes.legacy_math import compile_legacy_cpu
    return compile_legacy_cpu(partial(function, exponential=native_exp), static_argnames=static_argnames)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    print(json.dumps(build_library(args.compiler, args.output), indent=2))


if __name__ == "__main__":
    main()
