"""Independent radiation evidence, evolving snow feedback and cache boundaries."""

from dataclasses import asdict, replace
from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import jax
import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
sys.path.insert(0, str(ROOT/"examples"))

from plant_water import synthetic_case
from smr_jax.contracts import DailyForcing, RadiationComponents
from smr_jax.engine import run_days
from smr_jax.fluxes.legacy_math import compile_legacy_cpu
from smr_jax.initialization import initialize_legacy_state
from smr_jax.io import load_daily_forcings
from smr_jax.radiation import PreparedRadiationProvider, shortwave_from_components, radiation_error_bound
from smr_jax.state import GridSpec


class RadiationKernelTest(unittest.TestCase):
    def test_units_absorption_reflection_and_canopy_extremes(self):
        components = RadiationComponents(jnp.asarray([[100., 100., 100., 0.]]),
            jnp.asarray([[50., 50., 50., 0.]]), jnp.asarray([[10., 10., 10., 0.]]))
        albedo = jnp.asarray([[0., 1., .5, .2]])
        canopy = jnp.asarray([[0., 0., 1., 0.]])
        actual = compile_legacy_cpu(shortwave_from_components)(components, albedo, canopy)
        np.testing.assert_array_equal(actual, [[540., 0., 0., 0.]])
        half = compile_legacy_cpu(shortwave_from_components)(components, jnp.full((1, 4), .5), jnp.zeros((1, 4)))
        np.testing.assert_array_equal(half, [[279., 279., 279., 0.]])

    def test_quantization_bound_covers_independently_rounded_reflection(self):
        # Emulate a double-precision physical coefficient known before either
        # FCELL output, rather than using the already rounded unit basis as truth.
        rng = np.random.default_rng(36)
        beam = rng.uniform(0., 7000., (9, 13)).astype(np.float32).astype(np.float64)
        diffuse = rng.uniform(0., 3000., (9, 13)).astype(np.float32).astype(np.float64)
        coefficient = rng.uniform(0., 1200., (9, 13))
        unit = coefficient.astype(np.float32).astype(np.float64)
        albedo = rng.uniform(0., 1., (9, 13))
        canopy = rng.uniform(0., 1., (9, 13))
        reflected = (albedo.astype(np.float32).astype(np.float64)*coefficient).astype(np.float32)
        incident = (beam.astype(np.float32)+diffuse.astype(np.float32)+reflected).astype(np.float64)
        expected = (1.-canopy)*(1.-albedo)*incident*3600./1000.
        components = RadiationComponents(*map(jnp.asarray, (beam, diffuse, unit)))
        actual = np.asarray(compile_legacy_cpu(shortwave_from_components)(components, jnp.asarray(albedo), jnp.asarray(canopy)))
        self.assertTrue(np.all(np.abs(actual-expected) <= radiation_error_bound(components, albedo, canopy)))

    def test_snow_age_changes_shortwave_inside_daily_scan_without_replay(self):
        _, static, _, _ = synthetic_case()
        zero = jnp.asarray(0., dtype=jnp.float64)
        fields = lambda value: jnp.full((1, 1), value, dtype=jnp.float64)
        components = RadiationComponents(fields(100.), fields(20.), fields(5.))
        forcing = DailyForcing(zero+2015., zero+1., zero+1., zero+1., zero-5., zero-5., zero, zero,
            jnp.full((4,), -5.), zero, zero+100., zero+10., jnp.ones((6,)), None, components)
        days = jax.tree_util.tree_map(lambda value: jnp.stack((value, value)), forcing)
        initial = initialize_legacy_state(static)
        snow = replace(initial.snow, swe_cm=fields(2.), previous_swe_cm=fields(2.), ice_cm=fields(2.),
            age_days=fields(1.), surface_temp_c=fields(-5.))
        initial = replace(initial, snow=snow)
        _, history = compile_legacy_cpu(run_days)(initial, days, static)
        self.assertIsNone(days.replay_shortwave)
        self.assertGreater(float(history["albedo"][0, 0, 0]), float(history["albedo"][1, 0, 0]))
        self.assertLess(float(history["q_srad"][0, 0, 0]), float(history["q_srad"][1, 0, 0]))
        bare = replace(initial, snow=replace(snow, swe_cm=fields(0.), previous_swe_cm=fields(0.), ice_cm=fields(0.)))
        _, bare_history = compile_legacy_cpu(run_days)(bare, days, static)
        np.testing.assert_array_equal(bare_history["albedo"], np.full((2, 1, 1), .2))
        self.assertGreater(float(bare_history["q_srad"][0, 0, 0]), float(history["q_srad"][0, 0, 0]))

    def test_mixed_and_missing_radiation_sources_fail(self):
        from smr_jax.fluxes.snow import snow_energy
        from types import SimpleNamespace
        _, static, _, _ = synthetic_case()
        state = initialize_legacy_state(static)
        field = jnp.zeros((1, 1), dtype=jnp.float64)
        component = RadiationComponents(field, field, field)
        for replay, components in ((None, None), (field, component)):
            forcing = SimpleNamespace(replay_shortwave=replay, radiation_components=components,
                roughness_urban=jnp.asarray(100.), roughness_crop=jnp.asarray(10.))
            with self.assertRaisesRegex(ValueError, "radiation|provide|choose"):
                snow_energy(state, static, forcing, field, field, field, field)


class PreparedCacheTest(unittest.TestCase):
    def make_cache(self, folder):
        folder = Path(folder)
        grid = GridSpec(2, 2, 30., 100., 200.)
        bundle = folder/"bundle.zip"
        bundle.write_bytes(b"synthetic terrain provenance")
        arrays = {key: np.full((1, 2, 2), value, dtype=np.float32)
                  for key, value in (("beam", 100.), ("diffuse", 50.), ("reflected_unit", 10.))}
        metadata = {"schema_version": 2, "albedo_support": "finite-landuse-and-elevation",
            "provider": "grass-rsun-albedo-decomposition", "grid": asdict(grid),
            "units": "Wh/m2/day", "albedo_basis": 1., "r_sun_step_hours": .5, "terrain_shadowing": True,
            "cases": [{"date": "2015-01-01", "doy": 1, "cloud": .1, "linke": 3.}],
            "bundle_sha256": sha256(bundle.read_bytes()).hexdigest()}
        self.save_cache(folder, arrays, metadata)
        return grid, bundle, arrays, metadata

    def save_cache(self, folder, arrays, metadata):
        np.savez_compressed(Path(folder)/"components.npz", **arrays)
        metadata["components_sha256"] = sha256((Path(folder)/"components.npz").read_bytes()).hexdigest()
        (Path(folder)/"manifest.json").write_text(json.dumps(metadata), encoding="utf-8")

    def test_explicit_forcing_binding_and_no_replay_fallback(self):
        with TemporaryDirectory() as folder:
            grid, bundle, _, _ = self.make_cache(folder)
            cache = PreparedRadiationProvider(folder, grid, bundle, required_mask=np.ones((2, 2), dtype=bool))
            values = cache.components("2015-01-01", cloud=.1, linke=3., doy=1)
            self.assertEqual(values.beam.dtype, jnp.float64)
            for kwargs in ({"cloud": .2, "linke": 3., "doy": 1}, {"cloud": .1, "linke": 4., "doy": 1},
                           {"cloud": .1, "linke": 3., "doy": 2}):
                with self.assertRaisesRegex(ValueError, "weather differs"):
                    cache.components("2015-01-01", **kwargs)
            with self.assertRaisesRegex(ValueError, "no prepared"):
                cache.components("2015-01-02", cloud=.1, linke=3., doy=2)

    def test_changed_grid_bundle_and_corrupt_archive_fail(self):
        with TemporaryDirectory() as folder:
            grid, bundle, _, _ = self.make_cache(folder)
            with self.assertRaisesRegex(ValueError, "grid differs"):
                PreparedRadiationProvider(folder, replace(grid, xllcorner=115.), bundle)
            original = bundle.read_bytes()
            bundle.write_bytes(b"changed terrain")
            with self.assertRaisesRegex(ValueError, "provenance differs"):
                PreparedRadiationProvider(folder, grid, bundle)
            bundle.write_bytes(original)
            with (Path(folder)/"components.npz").open("ab") as stream:
                stream.write(b"corrupt")
            with self.assertRaisesRegex(ValueError, "checksum differs"):
                PreparedRadiationProvider(folder, grid, bundle)

    def test_missing_negative_nonfinite_wrong_type_and_wrong_shape_fail(self):
        for corruption in ("missing", "negative", "inf", "float64", "shape", "support", "required", "old-schema"):
            with self.subTest(corruption=corruption), TemporaryDirectory() as folder:
                grid, bundle, arrays, metadata = self.make_cache(folder)
                if corruption == "missing": arrays.pop("diffuse")
                elif corruption == "negative": arrays["beam"][0, 0, 0] = -1.
                elif corruption == "inf": arrays["beam"][0, 0, 0] = np.inf
                elif corruption == "float64": arrays["beam"] = arrays["beam"].astype(np.float64)
                elif corruption == "shape": arrays["beam"] = arrays["beam"][:, :1]
                elif corruption == "support": arrays["beam"][0, 0, 0] = np.nan
                elif corruption == "required":
                    for value in arrays.values(): value[0, 0, 0] = np.nan
                else: metadata["schema_version"] = 1
                self.save_cache(folder, arrays, metadata)
                with self.assertRaises(ValueError):
                    PreparedRadiationProvider(folder, grid, bundle, required_mask=np.ones((2, 2), dtype=bool))

    def test_weather_loader_supplies_components_and_explicitly_absent_replay(self):
        with TemporaryDirectory() as folder:
            grid, bundle, _, metadata = self.make_cache(folder)
            row = (ROOT/"reference/mica-creek-jan2015/gridmet_weather_2015-01.csv").read_text().splitlines()[0].split()
            metadata["cases"][0].update(cloud=float(row[16]), linke=float(row[15]))
            (Path(folder)/"manifest.json").write_text(json.dumps(metadata), encoding="utf-8")
            weather = Path(folder)/"weather.txt"
            weather.write_text(" ".join(row), encoding="utf-8")
            forcing = load_daily_forcings(weather, PreparedRadiationProvider(folder, grid, bundle))
            self.assertIsNone(forcing.replay_shortwave)
            self.assertEqual(forcing.radiation_components.beam.shape, (1, 2, 2))
            with self.assertRaisesRegex(ValueError, "requested 2 forcing days but found 1"):
                load_daily_forcings(weather, PreparedRadiationProvider(folder, grid, bundle), 2)


class GrassRadiationEvidenceTest(unittest.TestCase):
    def test_independent_grass_date_cloud_and_albedo_cases(self):
        from smr_jax.verify_radiation import component_report
        report = component_report(ROOT/"tests/fixtures")
        self.assertEqual(report["cases"], 60)
        self.assertEqual(report["representatives_per_case"], 128)
        self.assertLessEqual(report["max_fraction_of_bound"], 1.)

    def test_default_jit_meets_the_independent_source_precision_bound(self):
        with np.load(ROOT/"tests/fixtures/jax_verify05_radiation.npz", allow_pickle=False) as archive:
            components = RadiationComponents(*(jnp.asarray(archive[k]) for k in ("beam", "diffuse", "reflected_unit")))
            albedo, canopy, expected = (archive[k] for k in ("albedo", "canopy", "expected_q"))
        actual = np.asarray(jax.jit(shortwave_from_components)(components, jnp.asarray(albedo), jnp.asarray(canopy)))
        self.assertTrue(np.all(np.abs(actual-expected) <= radiation_error_bound(components, albedo, canopy)))

    def test_bad_or_incomplete_native_evidence_cannot_pass(self):
        from smr_jax.verify_radiation import component_report
        for corrupted in ("value", "cases"):
            with self.subTest(corrupted=corrupted), TemporaryDirectory() as directory:
                path = Path(directory)/"jax_verify05_radiation.npz"
                metadata = json.loads((ROOT/"tests/fixtures/jax_verify05_radiation.json").read_text())
                with np.load(ROOT/"tests/fixtures/jax_verify05_radiation.npz", allow_pickle=False) as archive:
                    arrays = {key: archive[key] for key in archive.files}
                if corrupted == "value": arrays["expected_q"][0, 0] += 1.
                else: metadata["cases"].pop()
                np.savez_compressed(path, **arrays)
                metadata.update(artifact_sha256=sha256(path.read_bytes()).hexdigest(), artifact_bytes=path.stat().st_size)
                path.with_suffix(".json").write_text(json.dumps(metadata), encoding="utf-8")
                with self.assertRaises((AssertionError, ValueError)):
                    component_report(directory)


if __name__ == "__main__":
    unittest.main()
