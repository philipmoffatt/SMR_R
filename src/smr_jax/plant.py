"""Opt-in, bounded two-layer vegetation coupling; the legacy engine is unchanged.

Water and daily fluxes are cm, rasters [y, x], layer fields [2, y, x].
Call ``external_water_step`` at the host boundary. ``external_day_kernel`` is
the JIT/scan-compatible numerical core for already validated inputs.
"""

from dataclasses import dataclass, fields

import jax
import jax.numpy as jnp
import numpy as np

from .fluxes.routing import _OFFSETS, _neighbor, legacy_receive_gather
from .fluxes.soil import _effective_conductivity, _reallocate
from .state import TwoLayerState


@dataclass(frozen=True)
class VegetationDemand:
    """Daily requests, NOT realized uptake. No SMR internal canopy or ET runs."""

    transp_cm: jax.Array
    soil_evap_cm: jax.Array
    interception_cm: jax.Array


@dataclass(frozen=True)
class SurfaceWater:
    """Explicit daily surface/snow handoff from ONE upstream provider.

    ``snow_vapor_cm`` is signed: condensation positive, sublimation negative.
    ``melt_cm`` is snow-store release, excluding rain already in ``rain_cm``.
    With no snow all snow fields must explicitly be zero. No canopy store is
    carried here: intercepted liquid evaporates on the same day.
    """

    precip_cm: jax.Array
    rain_cm: jax.Array
    snow_input_cm: jax.Array
    melt_cm: jax.Array
    swe_start_cm: jax.Array
    swe_end_cm: jax.Array
    snow_vapor_cm: jax.Array


for _cls in (VegetationDemand, SurfaceWater):
    jax.tree_util.register_dataclass(
        _cls, data_fields=[f.name for f in fields(_cls)], meta_fields=[])


def _routing_delivery(static):
    """Fraction of EACH DONOR actually delivered by the receiver-weight map."""
    mask = static.active_mask
    weights = jnp.where(mask[None], static.routing_weights, 0.0)
    delivered = sum(_neighbor(weights[k], -dr, -dc)
                    for k, (dr, dc) in enumerate(_OFFSETS))
    return jnp.where(mask, delivered, 0.0)


def validate_external_inputs(soil, static, demand, surface, bottom_boundary):
    """Host-only validation; no inferred soil properties or profile extension."""
    if bottom_boundary not in ("legacy_drainage", "closed"):
        raise ValueError("bottom_boundary must be 'legacy_drainage' or 'closed'")
    shape = static.grid.shape
    mask = np.asarray(static.active_mask)
    if mask.shape != shape or mask.dtype != np.bool_ or not mask.any():
        raise ValueError("active_mask must be a nonempty boolean grid")
    if not np.isfinite(static.grid.cellsize_m) or static.grid.cellsize_m <= 0:
        raise ValueError("cellsize_m must be positive and finite")
    if not np.array_equal(np.asarray(soil.active_mask), mask):
        raise ValueError("soil and static active masks must agree")

    def array(name, value, layer=False, nonnegative=True):
        value = np.asarray(value)
        expected = (2,) + shape if layer else shape
        if value.shape != expected or value.dtype != np.float64:
            raise ValueError(f"{name} must be float64 with shape {expected}")
        active = value[:, mask] if layer else value[mask]
        if not np.isfinite(active).all() or (nonnegative and (active < 0).any()):
            raise ValueError(f"{name} must be finite" + (" and nonnegative" if nonnegative else ""))
        return value

    names = ("soil_depth_a_cm", "soil_depth_b_cm", "fieldcap_a_cm", "fieldcap_b_cm",
             "sat_mc_a", "sat_mc_b", "wiltpt_mc_a", "wiltpt_amt_cm",
             "fieldcap_amt_cm", "sat_amt_cm", "slope_percent", "flowunits",
             "ksat_matrix_a", "ksat_matrix_b", "ksat_mpore_a", "ksat_mpore_b",
             "kfc_a", "kfc_b", "ksubsurface")
    for name in names:
        array(name, getattr(static, name))
    depth, wilt, fc, sat = _layer_properties(static, np)
    if np.any(depth[:, mask] <= 0):
        raise ValueError("layer thicknesses must be positive; no extension through restrictions")
    if np.any((wilt[:, mask] < 0) | (wilt[:, mask] > fc[:, mask]) |
              (fc[:, mask] > sat[:, mask]) | (sat[:, mask] <= 0) | (sat[:, mask] > depth[:, mask])):
        raise ValueError("need 0 <= WP <= FC <= saturation <= thickness, with positive saturation")
    for total, layers in ((static.fieldcap_amt_cm, fc), (static.sat_amt_cm, sat)):
        if not np.allclose(np.asarray(total)[mask], layers.sum(axis=0)[mask], rtol=0, atol=1e-10):
            raise ValueError("profile and layer bounds disagree")
    water = array("layer_water_cm", soil.layer_water_cm, layer=True)
    total = array("total_water_cm", soil.total_water_cm)
    if np.any(water[:, mask] > sat[:, mask]):
        raise ValueError("initial layer water exceeds saturation")
    if not np.allclose(total[mask], water.sum(axis=0)[mask], rtol=0, atol=1e-10):
        raise ValueError("legacy total/layer mismatch must be resolved explicitly before coupling")
    for name in ("transp_cm", "soil_evap_cm"):
        array(name, getattr(demand, name), layer=True)
    array("interception_cm", demand.interception_cm)
    for f in fields(surface):
        array(f.name, getattr(surface, f.name), nonnegative=f.name != "snow_vapor_cm")
    weights = np.asarray(static.routing_weights)
    if weights.shape != (8,) + shape or weights.dtype != np.float64:
        raise ValueError("routing_weights must be float64 [8, y, x]")
    if not np.isfinite(weights[:, mask]).all() or np.any(weights[:, mask] < 0):
        raise ValueError("active routing weights must be finite and nonnegative")
    if np.any(np.asarray(_routing_delivery(static))[mask] > 1.0 + 1e-12):
        raise ValueError("receiver weights duplicate donor water; fix routing upstream")


def _layer_properties(static, xp=jnp):
    depth = xp.stack((xp.asarray(static.soil_depth_a_cm), xp.asarray(static.soil_depth_b_cm)))
    wilt_a = xp.asarray(static.wiltpt_mc_a) * depth[0]
    # This is a source total minus its known component, not an invented B soil.
    wilt = xp.stack((wilt_a, xp.asarray(static.wiltpt_amt_cm) - wilt_a))
    fc = xp.stack((xp.asarray(static.fieldcap_a_cm), xp.asarray(static.fieldcap_b_cm)))
    sat = xp.stack((xp.asarray(static.sat_mc_a), xp.asarray(static.sat_mc_b))) * depth
    return depth, wilt, fc, sat


def _redistribute(water, static, sat):
    """Book the legacy allocation target as an explicit, bounded A/B transfer.

    Positive transfer is A -> B, negative is B -> A. Surface excess leaves A;
    excess arriving in B must first move upward and is recorded as such.
    This is legacy algebraic redistribution, NOT Darcy/Campbell drainage.
    """
    total = water.sum(axis=0)
    runoff = jnp.maximum(total - sat.sum(axis=0), 0.0)
    retained = total - runoff
    _, target_b = _reallocate(retained, water[0], static)
    target_b = jnp.clip(target_b, jnp.maximum(retained - sat[0], 0.0),
                        jnp.minimum(retained, sat[1]))
    transfer = target_b - water[1]
    return jnp.stack((retained - target_b, target_b)), transfer, runoff


def _lateral(water, static):
    flows = []
    for layer, suffix in enumerate(("a", "b")):
        eff = _effective_conductivity(
            water[layer], getattr(static, f"fieldcap_{suffix}_cm"),
            getattr(static, f"sat_mc_{suffix}"), getattr(static, f"soil_depth_{suffix}_cm"),
            getattr(static, f"ksat_matrix_{suffix}"), getattr(static, f"ksat_mpore_{suffix}"),
            getattr(static, f"kfc_{suffix}"))
        # Preserve legacy degrees/100 and the additional /100 factor literally.
        capacity = eff * static.slope_percent / 100.0 * getattr(static, f"soil_depth_{suffix}_cm") / static.grid.cellsize_m / 100.0
        flows.append(jnp.where(static.active_mask & (static.flowunits != 0),
                               jnp.minimum(jnp.maximum(capacity, 0.0), water[layer]), 0.0))
    out = jnp.stack(flows)
    weights = jnp.where(static.active_mask[None], static.routing_weights, 0.0)
    incoming = jnp.stack([legacy_receive_gather(out[i], weights, static.active_mask) for i in range(2)])
    export = out * jnp.maximum(1.0 - _routing_delivery(static), 0.0)[None]
    return water - out + incoming, out, incoming, export


def external_day_kernel(soil, static, demand, surface, *, bottom_boundary="legacy_drainage"):
    """Pure core: four equal six-hour steps, no internal canopy/ET/snow/basin.

    Validate before using under jit/scan; check ledgers and realized uptake
    before committing any candidate plant state. ``bottom_boundary`` is static.
    """
    mask = static.active_mask
    clean = lambda value: jnp.where(mask, value, 0.0)
    start = clean(soil.layer_water_cm)
    _, wilt, _, sat = _layer_properties(static)
    liquid = clean(surface.rain_cm + surface.melt_cm)
    interception = jnp.minimum(clean(demand.interception_cm), liquid)
    net = liquid - interception
    transp_request = clean(demand.transp_cm)
    evap_request = clean(demand.soil_evap_cm)

    def substep(water, _):
        water = water.at[0].add(net * .25)
        water, transfer1, runoff1 = _redistribute(water, static, sat)
        water, out, incoming, export = _lateral(clean(water), static)
        water, transfer2, runoff2 = _redistribute(water, static, sat)
        want_t, want_e = transp_request * .25, evap_request * .25
        want = want_t + want_e
        available = jnp.maximum(water - wilt, 0.0)
        scale = jnp.minimum(1.0, available / jnp.where(want > 0.0, want, 1.0))
        actual_t, actual_e = want_t * scale, want_e * scale
        water = water - actual_t - actual_e
        if bottom_boundary == "closed":
            bottom = jnp.zeros_like(net)
        else:
            bottom = jnp.minimum(jnp.minimum(jnp.maximum(water.sum(axis=0) - static.fieldcap_amt_cm, 0.0),
                                            water[1]), static.ksubsurface * .25)
        water = water.at[1].add(-bottom)
        water, transfer3, runoff3 = _redistribute(water, static, sat)
        transfers = jnp.stack((transfer1, transfer2, transfer3))
        water = clean(water)
        flux = {"vertical_down_cm": clean(jnp.maximum(transfers, 0.0).sum(axis=0))[None],
                "vertical_up_cm": clean(jnp.maximum(-transfers, 0.0).sum(axis=0))[None],
                "lateral_out_cm": out, "lateral_in_cm": incoming, "lateral_export_cm": export,
                "bottom_drainage_cm": clean(bottom), "runoff_cm": clean(runoff1 + runoff2 + runoff3),
                "transp_cm": clean(actual_t), "soil_evap_cm": clean(actual_e),
                "available_cm": clean(available)}
        return water, (water, flux)

    water, (history, hourly) = jax.lax.scan(substep, start, None, length=4)
    d = {name: values.sum(axis=0) for name, values in hourly.items() if name != "available_cm"}
    d.update({f.name: clean(getattr(surface, f.name)) for f in fields(surface)})
    d.update(interception_evap_cm=interception, requested_interception_cm=clean(demand.interception_cm),
             unmet_interception_cm=clean(demand.interception_cm) - interception,
             requested_transp_cm=transp_request, requested_soil_evap_cm=evap_request,
             unmet_transp_cm=jnp.maximum(transp_request - d["transp_cm"], 0.0),
             unmet_soil_evap_cm=jnp.maximum(evap_request - d["soil_evap_cm"], 0.0),
             infiltration_cm=net, substep_water_cm=history,
             substep_transp_cm=hourly["transp_cm"], substep_soil_evap_cm=hourly["soil_evap_cm"],
             substep_available_cm=hourly["available_cm"])
    signed = d["vertical_down_cm"][0] - d["vertical_up_cm"][0]
    d["drain_cm"] = jnp.stack((d["vertical_down_cm"][0], d["bottom_drainage_cm"]))
    layer_input = jnp.stack((net, jnp.zeros_like(net)))
    vertical_net_out = jnp.stack((signed, -signed))
    surface_and_bottom = jnp.stack((d["runoff_cm"], d["bottom_drainage_cm"]))
    d["layer_residual_cm"] = clean(water - start - layer_input + vertical_net_out +
                                  d["lateral_out_cm"] - d["lateral_in_cm"] +
                                  surface_and_bottom + d["transp_cm"] + d["soil_evap_cm"])
    d["snow_residual_cm"] = clean(surface.swe_end_cm - surface.swe_start_cm -
                                  surface.snow_input_cm - surface.snow_vapor_cm + surface.melt_cm)
    d["surface_partition_residual_cm"] = clean(surface.precip_cm - surface.rain_cm - surface.snow_input_cm)
    d["cell_residual_cm"] = clean(water.sum(axis=0) - start.sum(axis=0) +
        surface.swe_end_cm - surface.swe_start_cm - surface.precip_cm - surface.snow_vapor_cm +
        interception + (d["transp_cm"] + d["soil_evap_cm"]).sum(axis=0) + d["runoff_cm"] +
        d["bottom_drainage_cm"] + (d["lateral_out_cm"] - d["lateral_in_cm"]).sum(axis=0))
    # Domain ledger uses measured exports, independently of the cell residuals.
    d["domain_residual_m3"] = (jnp.sum(water - start) +
        jnp.sum(clean(surface.swe_end_cm - surface.swe_start_cm - surface.precip_cm - surface.snow_vapor_cm)) +
        jnp.sum(interception) + jnp.sum(d["transp_cm"] + d["soil_evap_cm"]) +
        jnp.sum(d["runoff_cm"] + d["bottom_drainage_cm"]) +
        jnp.sum(d["lateral_export_cm"])) * static.grid.cell_area_m2 * .01
    return TwoLayerState(water.sum(axis=0), water, mask), d


_LAYER_FLUXES = ("lateral_out_cm", "lateral_in_cm", "lateral_export_cm", "transp_cm", "soil_evap_cm",
                 "requested_transp_cm", "requested_soil_evap_cm", "unmet_transp_cm", "unmet_soil_evap_cm", "drain_cm")
_CELL_FLUXES = tuple(f.name for f in fields(SurfaceWater)) + (
    "bottom_drainage_cm", "runoff_cm", "interception_evap_cm", "requested_interception_cm",
    "unmet_interception_cm", "infiltration_cm")


def verify_water_balance(before, after, d, static, *, atol_cm=1e-10):
    """Fail closed on missing/nonfinite diagnostics; recompute ledgers with NumPy.

    Returns independent residuals. Tolerance is numerical, not fitted to data.
    Does not approve plant growth: call ``require_realized_uptake`` separately.
    """
    if not np.isfinite(atol_cm) or atol_cm <= 0:
        raise ValueError("atol_cm must be positive and finite")
    mask = np.asarray(static.active_mask)
    shape = mask.shape
    required = {name: (2,) + shape for name in _LAYER_FLUXES}
    required.update({name: shape for name in _CELL_FLUXES})
    required.update({name: (1,) + shape for name in ("vertical_down_cm", "vertical_up_cm")})
    required.update(layer_residual_cm=(2,) + shape, cell_residual_cm=shape,
                    snow_residual_cm=shape, surface_partition_residual_cm=shape, domain_residual_m3=(),
                    substep_water_cm=(4, 2) + shape, substep_transp_cm=(4, 2) + shape,
                    substep_soil_evap_cm=(4, 2) + shape, substep_available_cm=(4, 2) + shape)
    for name, expected in required.items():
        if name not in d:
            raise ValueError(f"missing required water diagnostic: {name}")
        value = np.asarray(d[name])
        if value.shape != expected or not np.isfinite(value).all():
            raise ValueError(f"invalid water diagnostic: {name}; expected finite {expected}")
    v = {key: np.asarray(value) for key, value in d.items()}
    w0 = np.where(mask, np.asarray(before.layer_water_cm), 0.0)
    w1 = np.where(mask, np.asarray(after.layer_water_cm), 0.0)
    _, _, _, sat = _layer_properties(static, np)
    for state in (before, after):
        if not np.array_equal(np.asarray(state.active_mask), mask):
            raise ValueError("state mask changed during coupling")
        if not np.isfinite(np.asarray(state.layer_water_cm)[:, mask]).all():
            raise FloatingPointError("nonfinite layer state")
        if not np.allclose(np.asarray(state.total_water_cm)[mask],
                           np.asarray(state.layer_water_cm).sum(axis=0)[mask], rtol=0, atol=atol_cm):
            raise FloatingPointError("total/layer storage mismatch")
    history = v["substep_water_cm"][:, :, mask]
    if np.any(history < -atol_cm) or np.any(history > sat[:, mask][None] + atol_cm):
        raise FloatingPointError("substep storage bound violated")
    if np.any(w1[:, mask] < -atol_cm) or np.any(w1[:, mask] > sat[:, mask] + atol_cm):
        raise FloatingPointError("final storage bound violated")
    for name in _LAYER_FLUXES + tuple(n for n in _CELL_FLUXES if n != "snow_vapor_cm") + ("vertical_down_cm", "vertical_up_cm"):
        if np.any(v[name] < -atol_cm):
            raise FloatingPointError(f"negative water flux: {name}")
    extracted = v["substep_transp_cm"] + v["substep_soil_evap_cm"]
    if np.any(extracted < -atol_cm) or np.any(extracted > v["substep_available_cm"] + atol_cm):
        raise FloatingPointError("substep extraction exceeds available water")
    for kind in ("transp", "soil_evap"):
        if not np.allclose(v[f"substep_{kind}_cm"].sum(axis=0), v[f"{kind}_cm"], rtol=0, atol=atol_cm):
            raise FloatingPointError(f"substep/daily extraction mismatch: {kind}")
    if not np.allclose(v["infiltration_cm"], v["rain_cm"] + v["melt_cm"] - v["interception_evap_cm"], rtol=0, atol=atol_cm):
        raise FloatingPointError("liquid/interception/infiltration mismatch")
    if not np.allclose(v["drain_cm"][0], v["vertical_down_cm"][0], rtol=0, atol=atol_cm) or not np.allclose(
            v["drain_cm"][1], v["bottom_drainage_cm"], rtol=0, atol=atol_cm):
        raise FloatingPointError("plant drainage view disagrees with actual vertical fluxes")
    # Rebuild routing from explicit donor -> receiver edges, independently of
    # the JAX gather and of EVERY water-budget residual.
    expected_in = np.zeros_like(v["lateral_in_cm"])
    delivered = np.zeros(shape)
    weights = np.asarray(static.routing_weights)
    rows, cols = shape
    for k, (dr, dc) in enumerate(_OFFSETS):
        r0, r1 = max(0, -dr), min(rows, rows - dr)
        c0, c1 = max(0, -dc), min(cols, cols - dc)
        recv = np.s_[r0:r1, c0:c1]
        donor = np.s_[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
        coeff = np.where(mask[recv] & mask[donor], weights[k][recv], 0.0)
        expected_in[(slice(None),) + recv] += v["lateral_out_cm"][(slice(None),) + donor] * coeff
        delivered[donor] += coeff
    expected_export = v["lateral_out_cm"] * (1.0 - delivered)
    if not np.allclose(expected_in, v["lateral_in_cm"], rtol=0, atol=atol_cm) or not np.allclose(
            expected_export, v["lateral_export_cm"], rtol=0, atol=atol_cm):
        raise FloatingPointError("routing diagnostics disagree with donor/receiver edges")
    signed = v["vertical_down_cm"][0] - v["vertical_up_cm"][0]
    layer = w1 - w0
    layer[0] += signed + v["runoff_cm"] - v["infiltration_cm"]
    layer[1] += v["bottom_drainage_cm"] - signed
    layer += v["lateral_out_cm"] - v["lateral_in_cm"] + v["transp_cm"] + v["soil_evap_cm"]
    snow = v["swe_end_cm"] - v["swe_start_cm"] - v["snow_input_cm"] - v["snow_vapor_cm"] + v["melt_cm"]
    partition = v["precip_cm"] - v["rain_cm"] - v["snow_input_cm"]
    cell = (w1.sum(axis=0) + v["swe_end_cm"]) - (w0.sum(axis=0) + v["swe_start_cm"])
    cell += (-v["precip_cm"] - v["snow_vapor_cm"] + v["interception_evap_cm"] +
             v["transp_cm"].sum(axis=0) + v["soil_evap_cm"].sum(axis=0) +
             v["runoff_cm"] + v["bottom_drainage_cm"] +
             v["lateral_out_cm"].sum(axis=0) - v["lateral_in_cm"].sum(axis=0))
    # Independent whole-domain storage vs external source/sinks, not sum(cell).
    domain = (w1.sum() + v["swe_end_cm"].sum() - w0.sum() - v["swe_start_cm"].sum() -
              v["precip_cm"].sum() - v["snow_vapor_cm"].sum() + v["interception_evap_cm"].sum() +
              v["transp_cm"].sum() + v["soil_evap_cm"].sum() + v["runoff_cm"].sum() +
              v["bottom_drainage_cm"].sum() + v["lateral_export_cm"].sum()) * static.grid.cell_area_m2 * .01
    residuals = {"layer_residual_cm": layer, "cell_residual_cm": cell, "snow_residual_cm": snow,
                 "surface_partition_residual_cm": partition, "domain_residual_m3": domain}
    for name, value in residuals.items():
        tolerance = atol_cm * static.grid.cell_area_m2 * .01 * mask.sum() if name.endswith("m3") else atol_cm
        if np.any(np.abs(value) > tolerance) or not np.allclose(v[name], value, rtol=0, atol=tolerance):
            raise FloatingPointError(f"water balance failed: {name} = {np.max(np.abs(value)):.3e}")
    for kind in ("transp", "soil_evap", "interception"):
        actual = v["interception_evap_cm"] if kind == "interception" else v[f"{kind}_cm"]
        if np.any(actual < -atol_cm) or np.any(v[f"unmet_{kind}_cm"] < -atol_cm) or not np.allclose(
                v[f"requested_{kind}_cm"], actual + v[f"unmet_{kind}_cm"], rtol=0, atol=atol_cm):
            raise FloatingPointError(f"requested/realized/unmet mismatch: {kind}")
    if not np.allclose(w1[:, mask], v["substep_water_cm"][-1][:, mask], rtol=0, atol=atol_cm):
        raise FloatingPointError("final storage differs from final substep")
    return residuals


def require_realized_uptake(diagnostics, *, atol_cm=1e-10):
    """Refuse committing candidate plant growth based on unsupplied water."""
    if not np.isfinite(atol_cm) or atol_cm <= 0:
        raise ValueError("atol_cm must be positive and finite")
    for kind in ("transp", "soil_evap", "interception"):
        key = f"unmet_{kind}_cm"
        if key not in diagnostics:
            raise ValueError(f"missing required uptake diagnostic: {key}")
        value = np.asarray(diagnostics[key])
        if not np.isfinite(value).all() or np.any(value > atol_cm):
            raise ValueError(f"unmet {kind}: reconcile plant water use and growth before committing state")


def external_water_step(soil, static, demand, surface, *, bottom_boundary="legacy_drainage", allow_unmet=False):
    """Checked daily entry point. Unmet demand raises unless explicitly requested.

    ``allow_unmet=True`` is for demand negotiation/diagnosis, not approval of a
    candidate plant state. It returns all requested, actual and unmet amounts.
    """
    validate_external_inputs(soil, static, demand, surface, bottom_boundary)
    after, diagnostics = external_day_kernel(soil, static, demand, surface, bottom_boundary=bottom_boundary)
    verify_water_balance(soil, after, diagnostics, static)
    if not allow_unmet:
        require_realized_uptake(diagnostics)
    return after, diagnostics


def demand_from_plant(plant, static):
    """Duck-type PlantToSoil, packed active cells in row-major order; no dependency."""
    mask = np.asarray(static.active_mask)
    count = int(mask.sum())

    def unpack(name, layer):
        value = np.asarray(getattr(plant, name), dtype=np.float64)
        if value.shape != ((2, count) if layer else (count,)):
            raise ValueError(f"PlantToSoil.{name} must use two layers and {count} packed active cells")
        grid = np.zeros((2,) + mask.shape if layer else mask.shape, dtype=np.float64)
        if layer:
            grid[:, mask] = value
        else:
            grid[mask] = value
        return jnp.asarray(grid)

    light = np.asarray(plant.ground_light_frac)
    if light.shape != (count,) or not np.isfinite(light).all() or np.any((light < 0) | (light > 1)):
        raise ValueError("ground_light_frac must be finite [active cells] within [0, 1]")
    return VegetationDemand(unpack("transp_cm", True), unpack("soil_evap_cm", True), unpack("interception_cm", False))


def soil_to_plant(soil, static, temperature_proxy_c, previous_diagnostics):
    """Return SoilToPlant kwargs, packed active cells, with previous-day ACTUAL flows.

    Temperature is an explicitly caller-supplied proxy; SMR has no soil-temperature
    model. Pass explicit zero prior fluxes on day one, never omit diagnostics.
    Also inspect vertical_up_cm; drain_cm alone cannot describe upward mixing.
    """
    mask = np.asarray(static.active_mask)
    depth, wilt, fc, sat = _layer_properties(static, np)
    for name in ("drain_cm", "lateral_out_cm", "vertical_up_cm"):
        if name not in previous_diagnostics:
            raise ValueError(f"missing previous-day actual flux: {name}")
        expected = ((1,) if name == "vertical_up_cm" else (2,)) + mask.shape
        value = np.asarray(previous_diagnostics[name])
        if value.shape != expected or not np.isfinite(value[:, mask]).all() or np.any(value[:, mask] < 0):
            raise ValueError(f"invalid previous-day flux: {name}")
    temperature = np.asarray(temperature_proxy_c)
    if temperature.shape != mask.shape or not np.isfinite(temperature[mask]).all():
        raise ValueError("temperature_proxy_c must be an explicit finite grid on active cells")
    return dict(water_cm=jnp.asarray(soil.layer_water_cm[:, mask]), wilt_cm=jnp.asarray(wilt[:, mask]),
                fieldcap_cm=jnp.asarray(fc[:, mask]), sat_cm=jnp.asarray(sat[:, mask]), depth_cm=jnp.asarray(depth[:, mask]),
                tsoil_c=jnp.asarray(temperature[mask]),
                drain_cm=jnp.asarray(previous_diagnostics["drain_cm"][:, mask]),
                lateral_cm=jnp.asarray(previous_diagnostics["lateral_out_cm"][:, mask]))
