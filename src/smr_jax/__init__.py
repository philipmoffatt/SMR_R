"""Two-layer, behavior-preserving foundation for the future SMR JAX engine."""

from .initialization import initialize_two_layer_reference_state
from .legacy import load_two_layer_legacy_inputs
from .state import GridSpec, TwoLayerState

__all__ = [
    "GridSpec",
    "TwoLayerState",
    "initialize_two_layer_reference_state",
    "load_two_layer_legacy_inputs",
]
