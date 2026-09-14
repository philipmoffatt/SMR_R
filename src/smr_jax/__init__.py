"""Two-layer, behavior-preserving foundation for the future SMR JAX engine."""

from .engine import run_days
from .initialization import initialize_legacy_state, initialize_two_layer_reference_state
from .legacy import load_legacy_bundle, load_two_layer_legacy_inputs
from .state import GridSpec, ModelState, TwoLayerState

__all__ = [
    "GridSpec",
    "ModelState",
    "TwoLayerState",
    "initialize_legacy_state",
    "initialize_two_layer_reference_state",
    "load_legacy_bundle",
    "load_two_layer_legacy_inputs",
    "run_days",
]
