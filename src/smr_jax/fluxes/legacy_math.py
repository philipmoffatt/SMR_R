"""Explicit CPU compilation contract for legacy map-by-map float64 arithmetic.

GRASS map algebra materializes successive double-precision operations. XLA's
algebraic simplification can replace division by multiplication by a rounded
reciprocal; fusion permits further reassociation across those boundaries.
Both passes must be disabled together for this legacy parity profile.
This changes compilation, not equations, states, or reference data.
"""

from functools import wraps

import jax


LEGACY_CPU_COMPILER_OPTIONS = {
    "xla_cpu_enable_fast_math": False,
    "xla_disable_hlo_passes": "algsimp,fusion",
}


def compile_legacy_cpu(function, *, static_argnames=()):
    """Compile an explicit host-entry CPU parity profile; do not nest under JIT.

    Keep this profile distinct from optimized/GPU performance runs. The
    installed JAX must support per-compilation ``compiler_options``. No global
    XLA flags, NumPy callbacks, or expected-value substitutions are used.
    """
    compiled = jax.jit(function, compiler_options=dict(LEGACY_CPU_COMPILER_OPTIONS),
                       static_argnames=static_argnames)

    @wraps(function)
    def run(*args, **kwargs):
        if not jax.config.jax_enable_x64:
            raise ValueError("legacy CPU parity requires jax_enable_x64=True")
        for leaf in jax.tree_util.tree_leaves((args, kwargs)):
            if isinstance(leaf, jax.Array) and any(d.platform != "cpu" for d in leaf.devices()):
                raise ValueError("legacy CPU parity requires CPU inputs; GPU parity is a separate gate")
        with jax.default_device(jax.devices("cpu")[0]):
            return compiled(*args, **kwargs)

    return run
