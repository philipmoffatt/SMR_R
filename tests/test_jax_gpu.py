"""GPU gate failure detection and synchronized harness checks, runnable on CPU."""

import contextlib
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import jax
import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

from smr_jax import verify_gpu as gate


class ComparisonBoundaryTest(unittest.TestCase):
    def setUp(self):
        self.mask = np.asarray([[True, False]])
        self.expected = {"water_cm": np.asarray([[1., np.nan]], dtype=np.float64)}

    def test_same_outputs_pass_and_changed_water_does_not(self):
        self.assertTrue(gate.compare_outputs(self.expected, self.expected, self.mask)["passed"])
        changed = {"water_cm": np.asarray([[1.+1e-8, np.nan]], dtype=np.float64)}
        result = gate.compare_outputs(self.expected, changed, self.mask)
        self.assertFalse(result["passed"])
        self.assertEqual(result["failed_fields"], 1)
        self.assertEqual(result["failures"][0]["outside_budget"], 1)

    def test_masks_shapes_types_and_required_finiteness_cannot_silently_pass(self):
        bad = ({}, {"water_cm": np.ones((1, 1))}, {"water_cm": np.asarray([[1., np.nan]], dtype=np.float32)},
               {"water_cm": np.asarray([[np.nan, np.nan]])}, {"water_cm": np.asarray([[1., 0.]])})
        for actual in bad:
            with self.subTest(actual=actual), self.assertRaises(ValueError):
                gate.compare_outputs(self.expected, actual, self.mask)

    def test_counters_and_booleans_are_exact(self):
        for key, before, after in (("snow_age_days", np.asarray([[2.]]), np.asarray([[np.nextafter(2., 3.)]])),
                                  ("active_mask", np.asarray([[True]]), np.asarray([[False]]))):
            result = gate.compare_outputs({key: before}, {key: after}, np.asarray([[True]]))
            self.assertFalse(result["passed"])

    def test_unknown_diagnostic_has_no_implicit_budget(self):
        with self.assertRaisesRegex(ValueError, "no declared"):
            gate.compare_outputs({"new_quantity": np.ones((1, 1))}, {"new_quantity": np.ones((1, 1))}, np.asarray([[True]]))

    def test_failed_parity_never_produces_a_validated_speedup(self):
        cpu_value = {"water_cm": np.asarray([[1.]])}
        gpu_value = {"water_cm": np.asarray([[1.1]])}
        with patch.object(gate, "execute", side_effect=[(cpu_value, {}, None), (gpu_value, {}, None)]), \
             patch.object(gate, "warmed_timing") as benchmark:
            report, _, _ = gate.verify_case(None, (), None, None, np.asarray([[True]]), repeats=5)
        benchmark.assert_not_called()
        self.assertFalse(report["comparison"]["passed"])
        self.assertNotIn("validated_compute_speedup", report["comparison"])


class ExecutionHarnessTest(unittest.TestCase):
    def test_dataclass_arguments_lower_and_warm_without_changing_first_timings(self):
        cpu = jax.devices("cpu")[0]
        with jax.default_device(cpu):
            state, forcing, static = gate.synthetic_model()
        function = lambda a, b, c: gate.run_days(a, b, c)
        with patch.object(gate.jax, "block_until_ready", wraps=jax.block_until_ready) as synchronize:
            output, metrics, resident = gate.execute(function, gate._host_tree((state, forcing, static)), cpu)
            first_compile = metrics["compilation_s"]
            gate.warmed_timing(resident, cpu, metrics, 3)
            self.assertGreaterEqual(synchronize.call_count, 5)
        self.assertEqual(metrics["compilation_s"], first_compile)
        self.assertEqual(len(metrics["warmed_execution_s"]), 3)
        self.assertGreater(metrics["output_bytes"], 0)
        self.assertTrue(gate.compare_outputs(output, output, np.asarray(static.active_mask))["passed"])

    def test_coupling_fixtures_pass_independent_water_ledgers(self):
        cpu = jax.devices("cpu")[0]
        with jax.default_device(cpu):
            _, _, static = gate.synthetic_model()
            cases = list(gate.coupling_cases(static))
            run = jax.jit(gate.external_day_kernel, static_argnames=("bottom_boundary",))
            for args, boundary in cases:
                after, diagnostics = run(*args, bottom_boundary=boundary)
                residuals = gate.verify_water_balance(args[0], after, diagnostics, args[1])
                self.assertLess(float(np.max(np.abs(residuals["layer_residual_cm"]))), 1e-10)
        self.assertEqual(len(cases), 10)


class ProvenanceAndCliTest(unittest.TestCase):
    def test_snapshot_checksums_and_path_boundary_are_enforced(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root/"module.py"
            source.write_text("value = 1\n", encoding="utf-8")
            metadata = {"schema_version": 1, "source_revision": "a"*40,
                        "files": {"module.py": sha256(source.read_bytes()).hexdigest()}}
            witness = root/"snapshot.json"
            witness.write_text(json.dumps(metadata), encoding="utf-8")
            self.assertEqual(gate.verify_snapshot(root, witness), "a"*40)
            source.write_text("value = 2\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "checksum differs"):
                gate.verify_snapshot(root, witness)
            metadata["files"] = {"../outside.py": "0"*64}
            witness.write_text(json.dumps(metadata), encoding="utf-8")
            with self.assertRaises(ValueError):
                gate.verify_snapshot(root, witness)

    def test_missing_cuda_fails_instead_of_falling_back(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)/"report.json"
            with patch.object(gate, "require_devices", side_effect=RuntimeError("CUDA GPU required")):
                with contextlib.redirect_stdout(StringIO()):
                    code = gate.main(["--radiation-dir", directory, "--report-path", str(output)])
            self.assertEqual(code, 1)
            self.assertFalse(output.exists())

    def test_existing_report_is_preserved_before_running(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)/"report.json"
            output.write_text("preserve", encoding="utf-8")
            with patch.object(gate, "verification_report") as verify, contextlib.redirect_stderr(StringIO()):
                with self.assertRaises(SystemExit):
                    gate.main(["--radiation-dir", directory, "--report-path", str(output)])
            verify.assert_not_called()
            self.assertEqual(output.read_text(), "preserve")


if __name__ == "__main__":
    unittest.main()
